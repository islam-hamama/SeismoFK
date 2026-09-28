"""
pmcc.py — Progressive Multi-Channel Correlation detector for SeismoFK v1.2.1.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.

This module implements the *core* of the PMCC method (Cansi, 1995) that is the
de-facto standard for IMS infrasound arrays: in each of a set of narrow
frequency bands, the inter-sensor time delays are measured by cross-correlation
and kept only where they are mutually *consistent* around closed sensor
triplets.  Consistent time–frequency points ("pixels") are converted to a
slowness vector (back-azimuth + trace velocity) by least squares.

:func:`pmcc` returns the pixels; :func:`pmcc_families` then aggregates pixels
that are adjacent in time and frequency and agree in back-azimuth and trace
velocity into *families* — one family per coherent arrival, the detection
product used for bulletins (Cansi, 1995; Cansi & Le Pichon, 2008).

References:
    Cansi, Y. (1995). An automatic seismic event processing for detection and
    location: The P.M.C.C. method. Geophys. Res. Lett., 22(9), 1021-1024.
    Cansi, Y., & Le Pichon, A. (2008). Infrasound event detection using the
    progressive multi-channel correlation algorithm. In Handbook of Signal
    Processing in Acoustics, Springer, 1425-1435.
"""

import itertools

import numpy as np
import pandas as pd
from scipy.signal import correlate, correlation_lags

from fk_analysis import station_xy_km


def log_bands(fmin, fmax, n_bands):
    """Return *n_bands* (f_low, f_high) tuples spaced logarithmically."""
    edges = np.geomspace(fmin, fmax, n_bands + 1)
    return [(float(edges[i]), float(edges[i + 1])) for i in range(n_bands)]


def circular_mean_deg(angles, weights=None):
    """Mean direction in degrees, preserving wraparound at north."""
    angles = np.asarray(angles, dtype=float)
    if angles.size == 0:
        return np.nan
    weights = np.ones_like(angles) if weights is None else np.asarray(weights, float)
    radians = np.deg2rad(angles)
    x = np.sum(weights * np.cos(radians))
    y = np.sum(weights * np.sin(radians))
    if np.hypot(x, y) <= 1e-12:
        return np.nan
    return float(np.degrees(np.arctan2(y, x)) % 360.0)


def circular_median_deg(angles):
    """
    Median direction in degrees that respects wraparound at north: angles
    are unwrapped around their circular mean before taking the median, so
    359° and 1° give 0°, not 180°.
    """
    angles = np.asarray(angles, dtype=float)
    centre = circular_mean_deg(angles)
    if np.isnan(centre):
        return np.nan
    offsets = (angles - centre + 180.0) % 360.0 - 180.0
    return float((centre + np.median(offsets)) % 360.0)


def _aligned_stream(st, stime=None, etime=None):
    """Sample all sensors on the same time grid over their common overlap."""
    if len(st) < 3:
        raise ValueError("PMCC needs at least 3 sensors for triplet consistency.")
    start = max(tr.stats.starttime for tr in st)
    end = min(tr.stats.endtime for tr in st)
    if stime is not None:
        start = max(start, stime)
    if etime is not None:
        end = min(end, etime)
    if end <= start:
        raise ValueError("PMCC sensors have no common time interval in the selected window.")

    fs = min(float(tr.stats.sampling_rate) for tr in st)
    npts = int(np.floor((end - start) * fs)) + 1
    if npts < 8:
        raise ValueError("PMCC common time interval is too short.")
    aligned = st.copy()
    for tr in aligned:
        if tr.stats.sampling_rate > fs * 1.01:
            tr.filter("lowpass", freq=0.45 * fs, corners=4, zerophase=True)
        if (tr.stats.starttime == start and
                np.isclose(tr.stats.sampling_rate, fs) and
                tr.stats.npts >= npts):
            tr.data = tr.data[:npts].copy()
        else:
            tr.interpolate(sampling_rate=fs, starttime=start, npts=npts,
                           method="linear")
    return aligned


def _pair_delays(seg, fs, max_lag_samp):
    """
    Cross-correlation time delays and correlation coefficients for every
    sensor pair within one window.

    Convention: b_ij = lag(i, j) / fs = (arrival at i) − (arrival at j)
                     = s · (r_i − r_j)        [seconds]
    so that the closure b_ij + b_jk + b_ki = 0 for a perfect plane wave.

    Returns
    -------
    pairs   : list of (i, j) with i < j
    delays  : (npairs,) seconds, b_ij as defined above
    corrs   : (npairs,) normalised peak correlation coefficient in [−1, 1]
    overlap : (npairs,) fraction of the window that overlaps at the peak lag,
              (n − |lag|) / n.  ``corrs / overlap`` corrects the energy
              normalisation for the shorter overlap and is the coherence
              estimate used for the delay uncertainty.
    """
    N = seg.shape[0]
    # Demean and pre-compute energies for normalisation.
    x = seg - seg.mean(axis=1, keepdims=True)
    energy = np.sqrt(np.sum(x ** 2, axis=1))
    lags_full = correlation_lags(seg.shape[1], seg.shape[1], mode='full')
    keep = np.abs(lags_full) <= max_lag_samp

    n_samp = seg.shape[1]
    pairs, delays, corrs, overlap = [], [], [], []
    for i, j in itertools.combinations(range(N), 2):
        cc = correlate(x[i], x[j], mode='full')[keep]
        lags = lags_full[keep]
        # A negative peak is an anti-correlation and gives a half-cycle timing
        # error on narrow-band signals; only positive peaks represent arrivals.
        k = int(np.argmax(cc))
        lag = float(lags[k])
        if 0 < k < len(cc) - 1:
            curvature = cc[k - 1] - 2 * cc[k] + cc[k + 1]
            if curvature < 0:
                offset = 0.5 * (cc[k - 1] - cc[k + 1]) / curvature
                lag += float(np.clip(offset, -1.0, 1.0))
        denom = energy[i] * energy[j]
        cmax = cc[k] / denom if denom > 0 else 0.0
        pairs.append((i, j))
        delays.append(lag / fs)
        corrs.append(cmax)
        overlap.append(max(n_samp - abs(lag), 1.0) / n_samp)
    return pairs, np.asarray(delays), np.asarray(corrs), np.asarray(overlap)


def _consistency(pairs, delays):
    """
    RMS of the triplet closure residuals (seconds).  Lower = more plane-wave-like.
    Uses the antisymmetry b_ji = −b_ij so all triplets are reachable from the
    i<j pairs.
    """
    lut = {p: d for p, d in zip(pairs, delays)}

    def b(a, c):
        return lut[(a, c)] if (a, c) in lut else -lut[(c, a)]

    idx = sorted({k for p in pairs for k in p})
    res = []
    for i, j, k in itertools.combinations(idx, 3):
        res.append(b(i, j) + b(j, k) + b(k, i))
    if not res:
        return 0.0
    return float(np.sqrt(np.mean(np.square(res))))


def _slowness_lstsq(positions, pairs, delays):
    """
    Solve b_ij = s · (r_i − r_j) for the 2-D slowness vector s (s/km) by least
    squares, then convert to back-azimuth (deg) and trace velocity (m/s).
    """
    G = np.array([positions[i] - positions[j] for (i, j) in pairs])  # (npairs,2) km
    s, *_ = np.linalg.lstsq(G, delays, rcond=None)                   # s/km
    sx, sy = float(s[0]), float(s[1])
    slow = float(np.hypot(sx, sy))
    baz = float(np.degrees(np.arctan2(-sx, -sy)) % 360.0)
    app_vel = float(1e3 / slow) if slow > 0 else np.inf
    return baz, app_vel, slow, sx, sy


def _delay_crlb(corrs, window_sec, f_lo, f_hi):
    """
    Cramér–Rao lower bound on a cross-correlation delay variance (s²) for a
    flat spectrum over [f_lo, f_hi] with coherence taken as the peak
    correlation ρ (Carter 1987):

        var(τ) = 3 (1 − ρ²) / (8 π² T ρ² (f_hi³ − f_lo³))

    Returns the mean over sensor pairs.
    """
    rho = np.clip(np.asarray(corrs, float), 0.05, 0.999)
    var = 3.0 * (1.0 - rho ** 2) / (
        8.0 * np.pi ** 2 * window_sec * rho ** 2 * (f_hi ** 3 - f_lo ** 3))
    return float(np.mean(var))


def _slowness_uncertainty(positions, pairs, delays, s_vec, pair_crlb,
                          timing_err=0.0):
    """
    Covariance of the slowness vector and the propagated back-azimuth and
    trace-velocity standard deviations for one pixel.

    Error model (Szuberla & Olson, 2004, arrival-time form): every sensor's
    arrival time has an independent error σₐ, so pair delays are correlated
    and Cov(s) = σₐ² (XcᵀXc)⁻¹ with Xc the centred sensor positions (km).
    (The naive pair form σ²(GᵀG)⁻¹ is too small by N/2.)  σₐ² is the
    largest of: the fit residual estimate RSS / (N (N − 3)) (needs N ≥ 4);
    half the pair-delay CRLB; and ``timing_err``² (clock / position error).

    Returns dict: delay_err (σₐ, s), s_cov (2×2, s²/km²), baz_std (deg),
    vel_std (m/s).
    """
    n = positions.shape[0]
    G = np.array([positions[i] - positions[j] for (i, j) in pairs])
    rss = float(np.sum((G @ s_vec - delays) ** 2))
    var_a = max(rss / (n * (n - 3)) if n > 3 else 0.0,
                pair_crlb / 2.0, timing_err ** 2)
    xc = positions - positions.mean(axis=0)
    cov = var_a * np.linalg.inv(xc.T @ xc)
    baz_std, vel_std = _baz_vel_std(s_vec, cov)
    return dict(delay_err=float(np.sqrt(var_a)), s_cov=cov,
                baz_std=baz_std, vel_std=vel_std)


def _baz_vel_std(s_vec, cov):
    """First-order propagation of a slowness covariance (s²/km²) to the
    back-azimuth (deg) and trace-velocity (m/s) standard deviations."""
    sx, sy = s_vec
    s2 = sx * sx + sy * sy
    if s2 <= 0:
        return np.nan, np.nan
    j_baz = np.array([sy, -sx]) / s2                       # rad per (s/km)
    j_vel = -1e3 * np.array([sx, sy]) / s2 ** 1.5          # (m/s) per (s/km)
    return (float(np.degrees(np.sqrt(j_baz @ cov @ j_baz))),
            float(np.sqrt(j_vel @ cov @ j_vel)))


def _combine_slowness(fam, pixel_overlap):
    """
    Inverse-variance (generalised least squares) combination of the pixel
    slowness vectors of one family.

    Pixels overlap in time (hop < window), so they are not independent: the
    formal covariance (Σ Wₖ)⁻¹ is inflated by ``pixel_overlap`` (window /
    hop), and further by the Birge ratio χ²/dof when the pixels scatter more
    than their errors allow.  χ² and its degrees of freedom are both reduced
    by ``pixel_overlap`` for the consistency p-value; a small p-value means
    the pixels are not one plane wave (e.g. two overlapping sources).

    Returns (s_vec, cov, chi2_p).
    """
    from scipy.stats import chi2 as chi2_dist

    s = fam[['sx', 'sy']].to_numpy(float)
    covs = np.stack([[fam['s_cxx'], fam['s_cxy']],
                     [fam['s_cxy'], fam['s_cyy']]]).transpose(2, 0, 1)
    weights = np.linalg.inv(covs)
    info = weights.sum(axis=0)
    cov_formal = np.linalg.inv(info)
    s_hat = cov_formal @ np.einsum('kij,kj->i', weights, s)
    resid = s - s_hat
    chi2 = float(np.einsum('ki,kij,kj->', resid, weights, resid))
    overlap = max(float(pixel_overlap), 1.0)
    dof = max(2.0 * (len(s) / overlap - 1.0), 1.0)
    chi2_eff = chi2 / overlap
    cov = cov_formal * overlap * max(1.0, chi2_eff / dof)
    return s_hat, cov, float(chi2_dist.sf(chi2_eff, dof))


def _direction_change(fam, pixel_overlap, n_sigma=3.0, min_half=3):
    """
    Largest significant back-azimuth difference (deg) between two halves of
    a family, split by time (early / late pixels) and by frequency (low /
    high bands); 0 when no split differs by more than ``n_sigma`` combined
    standard errors.

    Two merged sources show up as a systematic change of direction within
    the family.  Scatter that is merely larger than the random error model
    (wavefront curvature, sensor-position errors) inflates both halves'
    errors through the Birge factor and does not trigger it — unlike a
    plain χ² test, which rejects most real arrivals.
    """
    worst = 0.0
    n = len(fam)
    if n < 2 * min_half:
        return worst
    for key in ('time', 'f_center'):
        ordered = fam.sort_values(key)
        halves = (ordered.iloc[:n // 2], ordered.iloc[n // 2:])
        if (halves[0][key].iloc[-1] == halves[1][key].iloc[0]
                and key == 'f_center'):
            continue                     # all pixels in one band: no split
        estimates = []
        for half in halves:
            s_hat, cov, _ = _combine_slowness(half, pixel_overlap)
            baz = np.degrees(np.arctan2(-s_hat[0], -s_hat[1])) % 360.0
            estimates.append((baz, _baz_vel_std(s_hat, cov)[0]))
        (b1, e1), (b2, e2) = estimates
        diff = abs((b1 - b2 + 180.0) % 360.0 - 180.0)
        if diff > n_sigma * np.hypot(e1, e2):
            worst = max(worst, diff)
    return float(worst)


def pmcc(st, bands, window_sec, step_sec=None, consistency_max=0.05,
         corr_min=0.5, min_app_vel=200.0, max_app_vel=600.0,
         stime=None, etime=None, verbose=False, diagnostics=None,
         timing_err=0.0):
    """
    Run the PMCC-style detector.

    Parameters
    ----------
    st              : Stream with `tr.stats.coordinates` attached
    bands           : list of (f_low, f_high) tuples (see ``log_bands``)
    window_sec      : analysis window length per band (s)
    step_sec        : window hop (s); default = window_sec / 2
    consistency_max : max RMS closure residual to accept a pixel (s)
    corr_min        : min mean pairwise correlation to accept a pixel
    min/max_app_vel : trace-velocity gate (m/s) to reject incoherent solutions
    diagnostics     : optional dict populated with candidate/rejection counts,
                      plus ``windows``: a DataFrame with one row per candidate
                      window and band (time, f_center, mean_corr, consistency,
                      accepted) — including rejected windows, which show the
                      noise level the accepted pixels stand out from.

    timing_err      : per-sensor timing / position error expressed in
                      seconds, added to the uncertainty floor (default 0)

    Returns
    -------
    DataFrame: time, f_min, f_max, f_center, baz, app_vel, slowness,
               consistency, mean_corr, n_pairs, and per-pixel uncertainty:
               baz_std (deg), vel_std (m/s), delay_err (s), sx, sy (s/km)
               and the slowness covariance s_cxx, s_cxy, s_cyy (s²/km²)
               — see :func:`_slowness_uncertainty`.
    """
    aligned = _aligned_stream(st, stime=stime, etime=etime)
    positions, _, _ = station_xy_km(aligned)
    if np.linalg.matrix_rank(positions - positions.mean(axis=0)) < 2:
        raise ValueError("PMCC needs non-collinear sensor coordinates to estimate direction and velocity.")
    fs = float(aligned[0].stats.sampling_rate)
    if any(f_hi >= fs / 2 for _, f_hi in bands):
        raise ValueError(f"PMCC maximum frequency must be below Nyquist ({fs / 2:g} Hz).")
    rows = []
    counts = dict(candidates=0, low_correlation=0, inconsistent=0,
                  outside_velocity=0, accepted=0)
    windows = [] if diagnostics is not None else None

    for (f_lo, f_hi) in bands:
        a = aligned.copy()
        a.detrend('demean')
        a.filter('bandpass', freqmin=f_lo, freqmax=f_hi, corners=4,
                 zerophase=True)

        data = np.array([tr.data.astype(float) for tr in a])
        npts = data.shape[1]
        win_samp = int(window_sec * fs)
        if win_samp < 8 or win_samp > npts:
            if verbose:
                print(f"[WARN] band {f_lo:.2f}-{f_hi:.2f} Hz: window unusable")
            continue
        step = step_sec if step_sec is not None else window_sec / 2.0
        step_samp = max(1, int(step * fs))

        # Largest physically meaningful lag: array aperture / slowest velocity.
        aperture_km = np.max(np.linalg.norm(
            positions[:, None, :] - positions[None, :, :], axis=-1))
        max_lag_samp = int(np.ceil(aperture_km / (min_app_vel / 1e3) * fs)) + 2

        t0 = aligned[0].stats.starttime
        f_center = float(np.sqrt(f_lo * f_hi))

        for w0 in range(0, npts - win_samp + 1, step_samp):
            counts['candidates'] += 1
            seg = data[:, w0:w0 + win_samp]
            t_center = t0 + (w0 + win_samp / 2.0) / fs
            pairs, delays, corrs, overlap = _pair_delays(seg, fs,
                                                         max_lag_samp)
            mean_corr = float(np.mean(corrs))
            cons = np.nan
            accepted = False
            try:
                if mean_corr < corr_min:
                    counts['low_correlation'] += 1
                    continue
                cons = _consistency(pairs, delays)
                if cons > consistency_max:
                    counts['inconsistent'] += 1
                    continue
                baz, app_vel, slow, sx, sy = _slowness_lstsq(positions, pairs,
                                                             delays)
                if not (min_app_vel <= app_vel <= max_app_vel):
                    counts['outside_velocity'] += 1
                    continue
                accepted = True
                unc = _slowness_uncertainty(
                    positions, pairs, delays, np.array([sx, sy]),
                    _delay_crlb(corrs / overlap, window_sec, f_lo, f_hi),
                    timing_err)
            finally:
                if windows is not None:
                    windows.append((t_center.datetime, f_center, mean_corr,
                                    cons, accepted))
            counts['accepted'] += 1
            rows.append({
                'time':        t_center.datetime,
                'f_min':       f_lo,
                'f_max':       f_hi,
                'f_center':    f_center,
                'baz':         baz,
                'app_vel':     app_vel,
                'slowness':    slow,
                'consistency': cons,
                'mean_corr':   mean_corr,
                'n_pairs':     len(pairs),
                'baz_std':     unc['baz_std'],
                'vel_std':     unc['vel_std'],
                'delay_err':   unc['delay_err'],
                'sx':          sx,
                'sy':          sy,
                's_cxx':       float(unc['s_cov'][0, 0]),
                's_cxy':       float(unc['s_cov'][0, 1]),
                's_cyy':       float(unc['s_cov'][1, 1]),
            })

        if verbose:
            print(f"[INFO] band {f_lo:.2f}-{f_hi:.2f} Hz → "
                  f"{sum(1 for r in rows if r['f_min'] == f_lo)} pixels")

    cols = ['time', 'f_min', 'f_max', 'f_center', 'baz', 'app_vel',
            'slowness', 'consistency', 'mean_corr', 'n_pairs',
            'baz_std', 'vel_std', 'delay_err', 'sx', 'sy',
            's_cxx', 's_cxy', 's_cyy']
    if diagnostics is not None:
        diagnostics.update(counts)
        diagnostics['windows'] = pd.DataFrame(
            windows, columns=['time', 'f_center', 'mean_corr',
                              'consistency', 'accepted'])
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows)[cols].sort_values('time').reset_index(drop=True)


FAMILY_COLUMNS = ['family', 't_start', 't_end', 'duration_s', 'n_pixels',
                  'n_bands', 'f_min', 'f_max', 'f_center',
                  'baz', 'baz_ci95', 'baz_spread',
                  'app_vel', 'vel_ci95', 'vel_spread',
                  'chi2_p', 'mixed', 'mixed_dbaz',
                  'ell_major', 'ell_minor', 'ell_azimuth',
                  'mean_corr', 'max_corr', 'consistency']

#: χ²(2) quantile for a 95 % confidence ellipse.
_CHI2_95_2DOF = 5.991464547107979


def pmcc_families(pixels, time_tol, baz_tol=10.0, vel_tol=0.15, band_gap=1,
                  min_pixels=5, pixel_overlap=1.0):
    """
    Aggregate PMCC pixels into families (one family per coherent arrival).

    Two pixels are linked when all of the following hold, and a family is a
    connected group of linked pixels (single linkage):

    * |Δt| <= ``time_tol`` seconds (use the window length, so overlapping and
      consecutive windows connect);
    * their frequency bands are at most ``band_gap`` bands apart;
    * |Δ back-azimuth| <= ``baz_tol`` degrees (wraps at north);
    * |Δ velocity| <= ``vel_tol`` × the smaller trace velocity (relative).

    Groups smaller than ``min_pixels`` are isolated pixels, not families.

    Returns ``(pixels_out, families)``: a copy of *pixels* with a ``family``
    column (0 = not in a family) and a DataFrame with one row per family,
    numbered 1, 2, … in order of start time:

    * time span, pixel and band counts, frequency range;
    * ``baz`` / ``app_vel`` — inverse-variance estimate from the pixel
      slowness covariances (:func:`_combine_slowness`), with ``baz_ci95`` /
      ``vel_ci95`` half-widths of the 95 % confidence intervals and a 95 %
      slowness ellipse (``ell_major``, ``ell_minor`` in s/km, ``ell_azimuth``
      of the major axis in degrees from north);
    * ``baz_spread`` / ``vel_spread`` — scatter of the pixels themselves;
    * ``chi2_p`` — consistency of the pixels with one plane wave within
      their random errors (a diagnostic: real arrivals often score low
      because of wavefront curvature and position errors);
    * ``mixed`` — True when the direction changes significantly within the
      family (early vs late or low vs high band pixels differ by > 3 σ),
      i.e. probably two sources merged; ``mixed_dbaz`` is that change (deg).

    ``pixel_overlap`` is window / hop of the PMCC run (overlapping pixels are
    not independent).  Pixel tables without covariance columns (older
    archives) fall back to correlation-weighted means with NaN intervals.
    """
    out = pixels.copy()
    out['family'] = 0
    if pixels is None or pixels.empty:
        return out, pd.DataFrame(columns=FAMILY_COLUMNS)

    order = np.argsort(pd.to_datetime(pixels['time']).to_numpy())
    px = pixels.iloc[order].reset_index()
    t = (pd.to_datetime(px['time']) - pd.Timestamp(0)).dt.total_seconds() \
        .to_numpy()
    band_edges = np.sort(px['f_min'].unique())
    band = np.searchsorted(band_edges, px['f_min'].to_numpy())
    baz = px['baz'].to_numpy(float)
    vel = px['app_vel'].to_numpy(float)

    parent = np.arange(len(px))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(px)):
        j = i + 1
        while j < len(px) and t[j] - t[i] <= time_tol:
            dbaz = abs((baz[j] - baz[i] + 180.0) % 360.0 - 180.0)
            if (abs(band[j] - band[i]) <= band_gap and dbaz <= baz_tol
                    and abs(vel[j] - vel[i]) <= vel_tol * min(vel[i], vel[j])):
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
            j += 1

    roots = np.array([find(i) for i in range(len(px))])
    groups = [np.flatnonzero(roots == r) for r in np.unique(roots)]
    groups = [g for g in groups if len(g) >= min_pixels]
    groups.sort(key=lambda g: t[g].min())

    has_cov = {'sx', 'sy', 's_cxx', 's_cxy', 's_cyy'} <= set(px.columns)
    rows = []
    for number, g in enumerate(groups, start=1):
        fam = px.iloc[g]
        w = fam['mean_corr'].to_numpy(float)
        rad = np.deg2rad(fam['baz'].to_numpy(float))
        resultant = np.hypot(np.sum(w * np.sin(rad)),
                             np.sum(w * np.cos(rad))) / np.sum(w)
        times = pd.to_datetime(fam['time'])
        row = {
            'family': number,
            't_start': times.min(), 't_end': times.max(),
            'duration_s': float((times.max() - times.min()).total_seconds()),
            'n_pixels': len(g),
            'n_bands': int(fam['f_min'].nunique()),
            'f_min': float(fam['f_min'].min()),
            'f_max': float(fam['f_max'].max()),
            'f_center': float(np.exp(np.average(np.log(fam['f_center']),
                                                weights=w))),
            'baz': circular_mean_deg(fam['baz'], weights=w),
            'baz_ci95': np.nan,
            'baz_spread': float(np.degrees(np.sqrt(
                -2.0 * np.log(min(max(resultant, 1e-12), 1.0))))),
            'app_vel': float(fam['app_vel'].median()),
            'vel_ci95': np.nan,
            'vel_spread': float(fam['app_vel'].std(ddof=0)),
            'chi2_p': np.nan, 'mixed': False, 'mixed_dbaz': 0.0,
            'ell_major': np.nan, 'ell_minor': np.nan, 'ell_azimuth': np.nan,
            'mean_corr': float(w.mean()),
            'max_corr': float(w.max()),
            'consistency': float(fam['consistency'].mean()),
        }
        if has_cov and len(g) >= 2:
            s_hat, cov, chi2_p = _combine_slowness(fam, pixel_overlap)
            sx, sy = s_hat
            baz_std, vel_std = _baz_vel_std(s_hat, cov)
            evals, evecs = np.linalg.eigh(cov)
            major = evecs[:, 1]                         # (east, north)
            row.update(
                baz=float(np.degrees(np.arctan2(-sx, -sy)) % 360.0),
                app_vel=float(1e3 / np.hypot(sx, sy)),
                baz_ci95=1.96 * baz_std, vel_ci95=1.96 * vel_std,
                chi2_p=chi2_p,
                ell_major=float(np.sqrt(evals[1] * _CHI2_95_2DOF)),
                ell_minor=float(np.sqrt(evals[0] * _CHI2_95_2DOF)),
                ell_azimuth=float(np.degrees(np.arctan2(major[0], major[1]))
                                  % 180.0))
            row['mixed_dbaz'] = _direction_change(fam, pixel_overlap)
            row['mixed'] = row['mixed_dbaz'] > 0
        rows.append(row)
        out.loc[fam['index'].to_numpy(), 'family'] = number
    return out, pd.DataFrame(rows, columns=FAMILY_COLUMNS)
