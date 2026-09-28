"""
fk_analysis.py — Core FK / beamforming routines for SeismoFK (v1.2.1).

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.

This is the enhanced v1.2.1 analysis module.  It is a *self-contained* copy of
the v1.2.0 routines (load_inventory, get_coordinates_safe, compute_beam,
fk_array) plus the following new capabilities for the v1.2.1 release:

    • High-resolution adaptive beamforming  (Capon / MVDR and MUSIC) in
      addition to the classical Bartlett (conventional FK) estimator.
    • Cross-spectral density matrix (CSDM) estimation via Welch snapshots.
    • A windowed scanner (fk_scan) that produces one detection row per
      sliding window — the workhorse for the long-term batch CLI.
    • Bootstrap uncertainty on back-azimuth / apparent velocity, plus a
      slowness-space confidence ellipse.
    • Array Response Function (ARF) for the current geometry.

Acknowledgements:
    Conventional array processing methodology adapted in part from
    Jelle Assink (jelle.assink@knmi.nl)
    https://github.com/roseseismo/roses2021/blob/main/unit08/array_processing.py
"""

import os

import numpy as np
import pandas as pd

from obspy import read_inventory, UTCDateTime, geodetics
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import array_processing

try:
    from scipy.stats import circmean, circstd
except Exception:  # pragma: no cover - scipy is a hard dependency, but stay safe
    circmean = circstd = None

from scipy.signal import detrend as _detrend

__version__ = "1.2.1"

# Earth radius used for the local flat-Earth station projection (metres / km).
_R_EARTH_M = 6_371_000.0
_R_EARTH_KM = 6_371.0


# ─────────────────────────────────────────────────────────────────────────────
#  Inventory loader
# ─────────────────────────────────────────────────────────────────────────────

def load_inventory(inv_path, verbose=True):
    """
    Accept an already-loaded Inventory object or a path to a StationXML /
    RESP / dataless SEED file.  Returns an obspy Inventory.
    """
    from obspy.core.inventory import Inventory
    if isinstance(inv_path, Inventory):
        return inv_path

    if not os.path.exists(inv_path):
        raise FileNotFoundError(f"Inventory file not found: {inv_path}")

    filename = os.path.basename(inv_path).lower()
    ext      = os.path.splitext(filename)[1]

    if ext in ('.xml', '.stationxml'):
        try:
            inv = read_inventory(inv_path, format='STATIONXML')
            if verbose:
                print(f"[INFO] Loaded StationXML: {inv_path}")
            return inv
        except Exception as e:
            if verbose:
                print(f"[WARN] StationXML failed ({e}), trying auto-detect ...")

    if 'resp' in filename or ext == '.resp':
        try:
            return read_inventory(inv_path, format='RESP')
        except Exception as e:
            if verbose:
                print(f"[WARN] RESP failed ({e})")

    if ext in ('.seed', '.dataless') or 'dataless' in filename:
        try:
            return read_inventory(inv_path, format='SEED')
        except Exception as e:
            if verbose:
                print(f"[WARN] SEED failed ({e})")

    try:
        inv = read_inventory(inv_path)
        if verbose:
            print(f"[INFO] Loaded inventory (auto-detected): {inv_path}")
        return inv
    except Exception as e:
        raise ValueError(
            f"Could not read inventory '{inv_path}'.\n"
            f"Supported formats: StationXML, RESP, dataless SEED.\n"
            f"Error: {e}"
        )


# ─────────────────────────────────────────────────────────────────────────────
#  Coordinate lookup
# ─────────────────────────────────────────────────────────────────────────────

def get_coordinates_safe(inv, trace_id, datetime=None, verbose=True):
    """
    Return a coordinate dict for *trace_id*, trying several location-code
    fallbacks and ultimately station-level coordinates.
    """
    try:
        return inv.get_coordinates(trace_id, datetime=datetime)
    except Exception:
        pass

    net, sta, loc, cha = trace_id.split('.')
    for loc_try in ['', '--', '00', '10']:
        try:
            coords = inv.get_coordinates(f"{net}.{sta}.{loc_try}.{cha}",
                                         datetime=datetime)
            if verbose:
                print(f"[INFO] Matched {trace_id} with location '{loc_try}'")
            return coords
        except Exception:
            pass

    for network in inv:
        if network.code != net:
            continue
        for station in network:
            if station.code != sta:
                continue
            for channel in station:
                if channel.code == cha:
                    return {
                        'latitude':  channel.latitude  or station.latitude,
                        'longitude': channel.longitude or station.longitude,
                        'elevation': channel.elevation or station.elevation,
                    }
            if verbose:
                print(f"[WARN] Using station-level coords for {trace_id}")
            return {
                'latitude':  station.latitude,
                'longitude': station.longitude,
                'elevation': station.elevation,
            }

    raise LookupError(
        f"No coordinates found for {trace_id}.\n"
        f"Available: {[f'{n.code}.{s.code}' for n in inv for s in n]}"
    )


def attach_coordinates(st, inv, verbose=True):
    """
    Attach `tr.stats.coordinates` (lat/lon/elev) to every trace in *st* using
    *inv* (path or Inventory).  Returns the same stream for chaining.
    """
    inv = load_inventory(inv, verbose=verbose)
    for tr in st:
        coords = get_coordinates_safe(inv, tr.id,
                                      datetime=tr.stats.starttime,
                                      verbose=verbose)
        tr.stats.coordinates = AttribDict({
            'latitude':  coords['latitude'],
            'elevation': coords['elevation'],
            'longitude': coords['longitude'],
        })
    return st


def epoch_choices(inv, st):
    """
    Pick one StationXML channel epoch per trace in *st*.

    Some inventories (e.g. IMS metadata) contain a long nominal epoch that
    overlaps the dated calibration history, so several epochs cover the same
    time.  ObsPy's get_response / get_coordinates / remove_response silently
    use the first match; here the epoch with the latest start date — the most
    recent calibration — wins.

    Returns {trace_id: (network, station, channel, n_matching)} for every trace
    with at least one matching epoch.
    """
    choices = {}
    for tr in st:
        if tr.id in choices:
            continue
        net, sta, loc, cha = tr.id.split('.')
        t = tr.stats.starttime
        matches = [(n, s, c)
                   for n in inv.select(network=net, station=sta, location=loc,
                                       channel=cha, time=t)
                   for s in n for c in s]
        if matches:
            n, s, c = max(matches,
                          key=lambda m: m[2].start_date or UTCDateTime(0))
            choices[tr.id] = (n, s, c, len(matches))
    return choices


def resolve_epochs(inv, st):
    """
    Return an Inventory with exactly one channel epoch per trace of *st*
    (see :func:`epoch_choices`), so every later ObsPy lookup is unambiguous.
    Traces without an exact match keep their whole station from *inv*, so
    the location-code and station-level fallbacks of get_coordinates_safe
    still apply.
    """
    import copy
    from obspy.core.inventory import Inventory

    choices = epoch_choices(inv, st)
    networks = {}

    def _add(network, station, channels):
        net_copy = networks.get(network.code)
        if net_copy is None:
            net_copy = copy.copy(network)
            net_copy.stations = []
            networks[network.code] = net_copy
        sta_copy = copy.copy(station)
        sta_copy.channels = list(channels)
        net_copy.stations.append(sta_copy)

    for n, s, c, _ in choices.values():
        _add(n, s, [c])
    for tr in st:
        if tr.id not in choices:
            net, sta = tr.id.split('.')[:2]
            for n in inv.select(network=net, station=sta):
                for s in n:
                    _add(n, s, s.channels)
    return Inventory(networks=list(networks.values()),
                     source=inv.source, sender=inv.sender)


def scalar_calibrate(st, inv):
    """
    Divide every trace of *st* (in place) by its StationXML overall
    sensitivity and return the resulting unit label (e.g. 'Pa').

    All-or-nothing: if any trace lacks a usable sensitivity, or the traces
    would end up in different units, *st* is left in counts and 'counts' is
    returned, so calibrated and raw channels are never mixed.
    """
    factors, units = [], set()
    for tr in st:
        try:
            value = float(inv.get_response(tr.id, tr.stats.starttime)
                          .instrument_sensitivity.value)
        except Exception:
            return 'counts'
        if not (np.isfinite(value) and value > 0):
            return 'counts'
        factors.append(value)
        units.add(calibrated_units(inv, tr, 'DEF'))
    if len(units) != 1:
        return 'counts'
    for tr, value in zip(st, factors):
        tr.data = tr.data.astype(float) / value
    return units.pop()


def require_response_stages(inv, st):
    """
    Raise ValueError naming the traces whose StationXML response has only an
    overall sensitivity and no stages (common in IMS metadata).  ObsPy's
    remove_response fails on these with an unhelpful IndexError; callers use
    this to explain the fallback to scalar sensitivity division.
    """
    missing = []
    for tr in st:
        try:
            if not inv.get_response(tr.id, tr.stats.starttime).response_stages:
                missing.append(tr.id)
        except Exception:
            missing.append(tr.id)
    if missing:
        raise ValueError("StationXML has an overall sensitivity but no "
                         "response stages for " + ", ".join(missing))


# StationXML input-unit spellings → display label.
_UNIT_LABELS = {
    'PA': 'Pa', 'PASCAL': 'Pa', 'PASCALS': 'Pa',
    'M': 'm', 'M/S': 'm/s', 'M/S**2': 'm/s²', 'M/S2': 'm/s²', 'M/S/S': 'm/s²',
}
_MOTION_OUTPUT_UNITS = {'DISP': 'm', 'VEL': 'm/s', 'ACC': 'm/s²'}


def calibrated_units(inv, tr, output='VEL'):
    """
    Physical units of *tr* after calibration with *inv*.

    ObsPy converts ground-motion responses (m, m/s, m/s²) to the requested
    *output*, but applies a pressure (Pa) response as-is, so infrasound data
    stays in Pa even with output='VEL'.  Pass output='DEF' for scalar
    sensitivity division, which always leaves the sensor's input units.
    Returns 'counts' if the inventory has no response for *tr*.
    """
    try:
        response = inv.get_response(tr.id, datetime=tr.stats.starttime)
        raw = (response.instrument_sensitivity.input_units or '').strip()
    except Exception:
        return 'counts'
    if not raw:
        return 'uncalibrated units'
    label = _UNIT_LABELS.get(raw.upper(), raw)
    if label in ('m', 'm/s', 'm/s²') and output in _MOTION_OUTPUT_UNITS:
        return _MOTION_OUTPUT_UNITS[output]
    return label


def station_xy_km(st):
    """
    Local flat-Earth projection of station coordinates relative to the array
    centroid.  Returns (positions, ref_lat, ref_lon) where *positions* is an
    (N, 2) array of [x_east_km, y_north_km].
    """
    lats = np.array([tr.stats.coordinates.latitude  for tr in st], float)
    lons = np.array([tr.stats.coordinates.longitude for tr in st], float)
    ref_lat = float(np.mean(lats))
    ref_lon = float(np.mean(lons))
    x = np.deg2rad(lons - ref_lon) * _R_EARTH_KM * np.cos(np.deg2rad(ref_lat))
    y = np.deg2rad(lats - ref_lat) * _R_EARTH_KM
    return np.column_stack([x, y]), ref_lat, ref_lon


# ─────────────────────────────────────────────────────────────────────────────
#  Beamforming (time-domain delay-and-sum)
# ─────────────────────────────────────────────────────────────────────────────

def compute_beam(st, baz_deg, app_vel_ms, stime, etime, verbose=True):
    """
    Time-domain delay-and-sum beam for the given back-azimuth and apparent
    velocity.  Returns (beam_array, matplotlib_times).
    """
    a = st.copy()
    a.trim(starttime=stime, endtime=etime)

    if len(a) == 0 or any(len(tr.data) == 0 for tr in a):
        raise ValueError(
            "compute_beam: stream is empty after trimming. "
            "Check that start_time is within the data window."
        )

    a.detrend('demean')

    baz_rad = np.deg2rad(baz_deg)
    sx = -np.sin(baz_rad) / app_vel_ms
    sy = -np.cos(baz_rad) / app_vel_ms

    ref_lat = np.mean([tr.stats.coordinates.latitude  for tr in a])
    ref_lon = np.mean([tr.stats.coordinates.longitude for tr in a])

    beam = np.zeros(len(a[0].data))
    dt   = a[0].stats.delta

    for tr in a:
        dlat = np.deg2rad(tr.stats.coordinates.latitude  - ref_lat) * _R_EARTH_M
        dlon = (np.deg2rad(tr.stats.coordinates.longitude - ref_lon) * _R_EARTH_M *
                np.cos(np.deg2rad(ref_lat)))
        delay_samp = int(round((sx * dlon + sy * dlat) / dt))
        shifted = np.roll(tr.data, -delay_samp)
        if delay_samp > 0:
            shifted[-delay_samp:] = 0.0
        elif delay_samp < 0:
            shifted[:-delay_samp] = 0.0
        beam += shifted

    beam /= len(a)
    # Remove any residual linear trend (and mean) introduced by the delay-and-sum
    # stacking so the returned beam is centred on zero.
    beam = _detrend(beam, type='linear')
    return beam, a[0].times('matplotlib')


# ─────────────────────────────────────────────────────────────────────────────
#  Conventional FK (ObsPy array_processing) — unchanged v1.2.0 behaviour
# ─────────────────────────────────────────────────────────────────────────────

def fk_array(st, inv, freq_min, freq_max, win_length, overlap,
             start_t, duration, arr_name,
             source_lat, source_lon,
             max_slowness_skm=4.0, slowness_step_skm=0.16, prewhiten=0,
             abs_semb_thresh=0.35, save_csv=True,
             verbose=True):
    """
    Run conventional frequency-wavenumber array processing on *st*.

    New in v1.2.1: the slowness grid (``max_slowness_skm``, ``slowness_step_skm``)
    and ``prewhiten`` are now parameters instead of hard-coded constants, and
    CSV writing is optional (``save_csv``).  Defaults reproduce v1.2.0 output
    exactly (±4 s/km grid, step 0.16 s/km, no prewhitening).
    ``overlap`` is the legacy argument name for ObsPy's ``win_frac``: the
    fractional window step, so 0.1 means a 10% step and 90% overlap.

    ``source_lat``/``source_lon`` may be None (no presumed source), in which
    case ``expected_bazi`` is NaN.  ``med_baz`` is a circular median, so
    detections either side of north do not average to south.

    Returns the same result dict as v1.2.0 (see keys in the return statement).
    """
    from pmcc import circular_median_deg   # pmcc imports this module

    inv = load_inventory(inv, verbose=verbose)
    attach_coordinates(st, inv, verbose=verbose)

    # Expected back-azimuth from a presumed source.
    if source_lat is None or source_lon is None:
        bazz = np.nan
    else:
        bazz = geodetics.gps2dist_azimuth(
            source_lat, source_lon,
            st[0].stats.coordinates.latitude,
            st[0].stats.coordinates.longitude,
            a=6378137.0, f=0.0033528106647474805
        )[2]
        if verbose:
            print(f'[INFO] Expected Back-Azimuth >>> {bazz:.2f}°')

    stime = UTCDateTime(start_t)
    etime = stime + duration

    a = st.copy()
    a.trim(starttime=stime, endtime=etime)

    if len(a) == 0 or any(len(tr.data) == 0 for tr in a):
        raise ValueError(
            f"Stream is empty after trimming to [{stime}, {etime}].  "
            f"Verify that start_time '{start_t}' falls within the data window "
            f"[{st[0].stats.starttime}, {st[0].stats.endtime}]."
        )

    a.detrend('demean')

    rates = [tr.stats.sampling_rate for tr in a]
    if len(set(rates)) > 1:
        target_rate = min(rates)
        if verbose:
            print(f"[INFO] Mixed sampling rates {sorted(set(rates))} Hz — "
                  f"resampling all to {target_rate} Hz")
        for tr in a:
            if tr.stats.sampling_rate != target_rate:
                tr.resample(target_rate)

    a.filter('bandpass', freqmin=freq_min, freqmax=freq_max,
             corners=4, zerophase=True)

    actual_etime = min(tr.stats.endtime for tr in a)
    if actual_etime < etime:
        etime = actual_etime

    slm = max_slowness_skm
    out = array_processing(
        a,
        win_len=win_length,
        win_frac=overlap,
        frqlow=freq_min,
        frqhigh=freq_max,
        prewhiten=prewhiten,
        sll_x=-slm, slm_x=slm,
        sll_y=-slm, slm_y=slm,
        sl_s=slowness_step_skm,
        semb_thres=-1e9,
        vel_thres=-1e9,
        timestamp='julsec',
        stime=stime,
        etime=etime,
    )

    if out is None or len(out) == 0:
        raise ValueError(
            "array_processing returned no results. "
            "Check window length vs. duration and frequency range."
        )

    n_instr   = len(a)
    semblance = out[:, 1]
    fk_power  = out[:, 2]
    bazi      = out[:, 3] % 360
    slowness  = out[:, 4]           # s/km
    app_vel   = 1e3 / slowness      # m/s

    fisher = (n_instr - 1) * semblance / (1.0 - semblance + 1e-12)

    mask = semblance >= abs_semb_thresh
    if mask.sum() == 0:
        mask = semblance >= np.percentile(semblance, 75)
        if verbose:
            print(f'[WARN] No windows above semblance={abs_semb_thresh}; '
                  f'falling back to top-25 % percentile')
    med_baz = circular_median_deg(bazi[mask])
    med_vel = float(np.median(app_vel[mask]))
    if verbose:
        print(f'[INFO] Median beam → baz={med_baz:.1f}°, app_vel={med_vel:.0f} m/s '
              f'({mask.sum()}/{len(semblance)} windows used)')

    beam_wave, beam_times = compute_beam(
        a, med_baz, med_vel, stime, actual_etime, verbose=verbose
    )

    positions, ref_lat, ref_lon = station_xy_km(a)
    station_x_km = positions[:, 0]
    station_y_km = positions[:, 1]
    sta_ids = [tr.stats.station for tr in a]

    if save_csv:
        csv_file = f"Output_{arr_name}.csv"
        pd.DataFrame({
            'time':      pd.to_datetime(out[:, 0], unit='s'),
            'semblance': semblance,
            'fk_power':  fk_power,
            'fisher':    fisher,
            'bazi':      bazi,
            'app_vel':   app_vel,
            'slowness':  slowness,
        }).to_csv(csv_file, index=False)
        if verbose:
            print(f'[INFO] CSV saved → {csv_file}')

    return {
        'time':           pd.to_datetime(out[:, 0], unit='s'),
        'fisher':         fisher,
        'semblance':      semblance,
        'bazi':           bazi,
        'app_vel':        app_vel,
        'slowness':       slowness,
        'expected_bazi':  bazz,
        'med_baz':        med_baz,
        'med_vel':        med_vel,
        'beam_waveform':  beam_wave,
        'beam_times':     beam_times,
        'waveform_times': a[0].times('matplotlib'),
        'waveform_data':  a[0].data,
        'station_name':   a[0].stats.station,
        'n_stations':     n_instr,
        'station_x_km':   station_x_km,
        'station_y_km':   station_y_km,
        'station_ids':    sta_ids,
        'array_lat':      ref_lat,
        'array_lon':      ref_lon,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  NEW (v1.2.1): high-resolution adaptive beamforming
# ─────────────────────────────────────────────────────────────────────────────

def cross_spectral_matrix(st, fmin, fmax, n_segments=16, overlap=0.5,
                          seg_len=None):
    """
    Estimate the per-frequency cross-spectral density matrix (CSDM) using
    overlapping Welch snapshots.

    Returns
    -------
    freqs : (F,) ndarray         frequency bins inside [fmin, fmax]
    R     : (F, N, N) ndarray    complex CSDM, one Hermitian matrix per bin
    n_snap: int                  number of snapshots averaged

    The number of snapshots controls the rank of R; Capon/MUSIC need
    ``n_snap >= N`` for a well-conditioned matrix (diagonal loading in
    ``fk_beamform`` relaxes this).
    """
    data = np.array([tr.data.astype(float) for tr in st])
    N, npts = data.shape
    fs = float(st[0].stats.sampling_rate)

    if seg_len is None:
        seg_len = max(32, npts // n_segments)
    seg_len = min(seg_len, npts)
    hop = max(1, int(seg_len * (1.0 - overlap)))
    window = np.hanning(seg_len)

    freqs_all = np.fft.rfftfreq(seg_len, d=1.0 / fs)
    band = np.where((freqs_all >= fmin) & (freqs_all <= fmax))[0]
    if band.size == 0:  # fall back to the single nearest bin
        band = np.array([np.argmin(np.abs(freqs_all - 0.5 * (fmin + fmax)))])

    R = np.zeros((band.size, N, N), dtype=complex)
    n_snap = 0
    for s0 in range(0, npts - seg_len + 1, hop):
        seg = data[:, s0:s0 + seg_len] * window
        X = np.fft.rfft(seg, axis=1)[:, band]   # (N, F)
        for fi in range(band.size):
            v = X[:, fi]
            R[fi] += np.outer(v, np.conj(v))
        n_snap += 1

    if n_snap > 0:
        R /= n_snap
    return freqs_all[band], R, n_snap


def _slowness_grid(smax, ds):
    """Symmetric slowness axis in s/km."""
    return np.arange(-smax, smax + ds / 2.0, ds)


def _steering_grid(positions, sx_axis, sy_axis, freq):
    """
    Steering matrix for every (sx, sy) grid node at a single frequency.

    Returns A of shape (G, N) where G = len(sy_axis) * len(sx_axis), ordered
    so that A.reshape(len(sy_axis), len(sx_axis), N) indexes [iy, ix].
    """
    SX, SY = np.meshgrid(sx_axis, sy_axis)          # (ny, nx)
    sx = SX.ravel()
    sy = SY.ravel()
    # A signal of slowness s arrives at station r delayed by τ = s·r, i.e. with
    # spectral phase exp(-j2πf τ).  The matched steering vector therefore uses
    # the conjugate sign so that the Bartlett/Capon/MUSIC map peaks at the true
    # slowness (and back-azimuth, matching compute_beam / PMCC conventions).
    phase = 2.0 * np.pi * freq * (np.outer(sx, positions[:, 0]) +
                                  np.outer(sy, positions[:, 1]))
    return np.exp(-1j * phase)                      # (G, N)


def fk_beamform(st, fmin, fmax, method='capon', smax=4.0, ds=0.1,
                n_segments=16, n_sources=1, diag_load=0.01, verbose=False):
    """
    Frequency-domain slowness power map by Bartlett, Capon (MVDR) or MUSIC.

    Parameters
    ----------
    method     : 'bartlett' | 'capon' | 'music'
    smax, ds   : slowness grid half-extent and step (s/km)
    n_sources  : number of signal eigenvectors for MUSIC (default 1)
    diag_load  : diagonal-loading fraction (× mean eigenvalue) for numerical
                 stability of Capon / MUSIC when snapshots are few.

    Returns
    -------
    sx_axis, sy_axis : (nx,), (ny,) slowness axes (s/km)
    P                : (ny, nx) power map (linear), frequency-averaged
    """
    method = method.lower()
    if method not in ('bartlett', 'capon', 'music'):
        raise ValueError(f"Unknown beamforming method '{method}'")

    positions, _, _ = station_xy_km(st)
    N = positions.shape[0]
    if method == 'music' and not (1 <= n_sources < N):
        raise ValueError("MUSIC needs 1 <= n_sources < number of sensors")
    freqs, R, n_snap = cross_spectral_matrix(st, fmin, fmax,
                                             n_segments=n_segments)
    if freqs.size == 0:
        raise ValueError("No frequency bins in the requested band.")

    sx_axis = _slowness_grid(smax, ds)
    sy_axis = _slowness_grid(smax, ds)
    G = sy_axis.size * sx_axis.size
    P = np.zeros(G)

    eye = np.eye(N)
    for fi, f in enumerate(freqs):
        Rf = R[fi]
        load = diag_load * np.real(np.trace(Rf)) / N
        Rf = Rf + load * eye
        A = _steering_grid(positions, sx_axis, sy_axis, f)   # (G, N)

        if method == 'bartlett':
            # P = aᴴ R a   (real, positive)
            P += np.real(np.einsum('gi,ij,gj->g', np.conj(A), Rf, A))
        elif method == 'capon':
            Rinv = np.linalg.pinv(Rf)
            denom = np.real(np.einsum('gi,ij,gj->g', np.conj(A), Rinv, A))
            P += 1.0 / np.maximum(denom, 1e-30)
        else:  # music
            w, V = np.linalg.eigh(Rf)          # ascending eigenvalues
            idx = np.argsort(w)[::-1]
            noise = V[:, idx[n_sources:]]       # (N, N-n_sources)
            En = noise @ noise.conj().T         # noise-subspace projector
            denom = np.real(np.einsum('gi,ij,gj->g', np.conj(A), En, A))
            P += 1.0 / np.maximum(denom, 1e-30)

    P /= freqs.size
    P = P.reshape(sy_axis.size, sx_axis.size)
    if verbose:
        print(f"[INFO] {method} map: {N} sensors, {freqs.size} freq bins, "
              f"{n_snap} snapshots")
    return sx_axis, sy_axis, P


def estimate_slowness_peak(sx_axis, sy_axis, P):
    """
    Locate the maximum of a slowness power map and convert it to physical
    quantities.

    Returns dict: baz (deg), app_vel (m/s), slowness (s/km), sx, sy, power,
    coherence (peak / mean power, a sharpness proxy).
    """
    iy, ix = np.unravel_index(np.argmax(P), P.shape)
    sx = float(sx_axis[ix])
    sy = float(sy_axis[iy])
    slow = float(np.hypot(sx, sy))                       # s/km
    app_vel = float(1e3 / slow) if slow > 0 else np.inf  # m/s
    baz = float(np.degrees(np.arctan2(-sx, -sy)) % 360.0)
    pmax = float(P[iy, ix])
    pmean = float(np.mean(P))
    return {
        'baz': baz, 'app_vel': app_vel, 'slowness': slow,
        'sx': sx, 'sy': sy, 'power': pmax,
        'coherence': pmax / pmean if pmean > 0 else np.nan,
    }


def fk_scan(st, fmin, fmax, win_length, overlap=0.5, method='capon',
            smax=4.0, ds=0.1, stime=None, etime=None,
            n_segments=8, n_sources=1, verbose=False):
    """
    Slide a window across the stream and run a high-resolution beamformer in
    each window.  Produces one detection row per window — the core routine
    behind the long-term batch CLI.

    Returns a DataFrame: time, baz, app_vel, slowness, power, coherence, method.
    *time* is the window-centre UTC timestamp.
    """
    a = st.copy()
    if stime is not None or etime is not None:
        a.trim(starttime=stime, endtime=etime)
    if len(a) == 0 or any(len(tr.data) == 0 for tr in a):
        return pd.DataFrame(columns=['time', 'baz', 'app_vel', 'slowness',
                                     'power', 'coherence', 'method'])
    a.detrend('demean')
    a.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)

    t0 = max(tr.stats.starttime for tr in a)
    t1 = min(tr.stats.endtime for tr in a)
    step = max(win_length * (1.0 - overlap), 1.0 / a[0].stats.sampling_rate)

    rows = []
    first_error = None
    t = t0
    while t + win_length <= t1 + 1e-6:
        win = a.slice(t, t + win_length)
        if len(win) == len(a) and all(len(tr.data) > 8 for tr in win):
            try:
                sx_axis, sy_axis, P = fk_beamform(
                    win, fmin, fmax, method=method, smax=smax, ds=ds,
                    n_segments=n_segments, n_sources=n_sources)
                pk = estimate_slowness_peak(sx_axis, sy_axis, P)
                pk['time'] = (t + win_length / 2.0).datetime
                pk['method'] = method
                rows.append(pk)
            except Exception as e:
                if first_error is None:
                    first_error = e
                if verbose:
                    print(f"[WARN] window @ {t}: {e}")
        t += step

    if not rows:
        if first_error is not None:
            raise RuntimeError(f"Every FK window failed: {first_error}") from first_error
        return pd.DataFrame(columns=['time', 'baz', 'app_vel', 'slowness',
                                     'power', 'coherence', 'method'])
    df = pd.DataFrame(rows)
    return df[['time', 'baz', 'app_vel', 'slowness', 'power',
               'coherence', 'method']]


# ─────────────────────────────────────────────────────────────────────────────
#  NEW (v1.2.1): bootstrap uncertainty
# ─────────────────────────────────────────────────────────────────────────────

def bootstrap_beam(bazi, app_vel, slowness, mask=None, n_boot=1000,
                   ci=95.0, seed=0, block=1, slowness_step=0.0):
    """
    Bootstrap the back-azimuth / apparent-velocity estimate over the selected
    detection windows and derive a slowness-space confidence ellipse.

    Parameters
    ----------
    bazi, app_vel, slowness : per-window arrays (deg, m/s, s/km)
    mask                    : boolean selection (e.g. semblance >= thresh);
                              None → use all windows
    n_boot                  : number of bootstrap resamples
    ci                      : central confidence interval percentage
    block                   : moving-block length in windows.  Overlapping
                              FK windows are strongly correlated (a 0.1 step
                              means 90 % overlap), so resampling single
                              windows as if independent gives intervals that
                              are far too narrow.  Pass ≈ 1 / step fraction
                              (window / hop) to resample runs of consecutive
                              selected windows instead.
    slowness_step           : FK grid step (s/km).  Every window picks the
                              nearest grid node, so windows share the same
                              quantisation error and the bootstrap cannot see
                              it; each interval is widened to at least the
                              half-diagonal of a grid cell (ds/√2) converted to
                              back-azimuth and velocity.  ``resolution_limited``
                              reports when this floor, not the scatter, sets
                              the interval.

    Returns
    -------
    dict with baz_mean/std/ci, vel_mean/std/ci, n_used, and ``ellipse``
    (centre sx/sy, semi-axes a/b in s/km, orientation deg) at the requested CI.
    """
    bazi = np.asarray(bazi, float)
    app_vel = np.asarray(app_vel, float)
    slowness = np.asarray(slowness, float)
    if mask is None:
        mask = np.ones(bazi.shape, bool)
    baz_sel = bazi[mask]
    vel_sel = app_vel[mask]
    slow_sel = slowness[mask]
    n = baz_sel.size
    if n < 2:
        raise ValueError("bootstrap_beam needs at least 2 selected windows.")

    rng = np.random.default_rng(seed)
    block = int(min(max(1, round(block)), n))
    starts_max = n - block + 1
    n_blocks = int(np.ceil(n / block))
    lo = (100.0 - ci) / 2.0
    hi = 100.0 - lo

    def _cmean(deg):
        if circmean is not None:
            return float(np.degrees(circmean(np.deg2rad(deg))))
        # fallback: vector mean
        return float(np.degrees(np.arctan2(np.mean(np.sin(np.deg2rad(deg))),
                                           np.mean(np.cos(np.deg2rad(deg))))) % 360)

    boot_baz = np.empty(n_boot)
    boot_vel = np.empty(n_boot)
    for b in range(n_boot):
        starts = rng.integers(0, starts_max, n_blocks)
        idx = (starts[:, None] + np.arange(block)).ravel()[:n]
        boot_baz[b] = _cmean(baz_sel[idx])
        boot_vel[b] = np.median(vel_sel[idx])

    # Circular spread of the bootstrap back-azimuths.
    if circstd is not None:
        baz_std = float(np.degrees(circstd(np.deg2rad(boot_baz))))
    else:
        baz_std = float(np.std(boot_baz))

    # Slowness-vector scatter of the selected windows → covariance ellipse.
    sx = -np.sin(np.deg2rad(baz_sel)) * slow_sel
    sy = -np.cos(np.deg2rad(baz_sel)) * slow_sel
    cov = np.cov(np.vstack([sx, sy]))
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    evals = np.clip(evals[order], 0, None)
    evecs = evecs[:, order]
    # Chi-square scale for a 2-D confidence region.
    from scipy.stats import chi2
    scale = np.sqrt(chi2.ppf(ci / 100.0, df=2))
    ellipse = {
        'cx': float(np.mean(sx)), 'cy': float(np.mean(sy)),
        'a':  float(scale * np.sqrt(evals[0])),
        'b':  float(scale * np.sqrt(evals[1])),
        'theta_deg': float(np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0]))),
    }

    baz_mean = _cmean(boot_baz)
    # Percentiles must be taken on angles unwrapped around the circular mean;
    # ordinary percentiles are wrong when samples straddle north (0/360 deg).
    baz_offsets = (boot_baz - baz_mean + 180.0) % 360.0 - 180.0
    baz_ci = (baz_mean + np.percentile(baz_offsets, [lo, hi])) % 360.0

    # Grid-resolution floor (see ``slowness_step``).
    resolution_limited = False
    if slowness_step > 0:
        s_med = float(np.median(slow_sel))
        v_med = float(np.median(vel_sel))
        half = slowness_step / np.sqrt(2.0)                    # s/km
        baz_floor = float(np.degrees(half / s_med)) if s_med > 0 else 180.0
        vel_floor = v_med ** 2 * half / 1e3                    # m/s
        lo_off = (baz_ci[0] - baz_mean + 180.0) % 360.0 - 180.0
        hi_off = (baz_ci[1] - baz_mean + 180.0) % 360.0 - 180.0
        if lo_off > -baz_floor or hi_off < baz_floor:
            resolution_limited = True
            baz_ci = (baz_mean + np.array([min(lo_off, -baz_floor),
                                           max(hi_off, baz_floor)])) % 360.0
        vel_ci = [float(np.percentile(boot_vel, lo)),
                  float(np.percentile(boot_vel, hi))]
        vel_mid = float(np.mean(boot_vel))
        if vel_ci[0] > vel_mid - vel_floor or vel_ci[1] < vel_mid + vel_floor:
            resolution_limited = True
            vel_ci = [min(vel_ci[0], vel_mid - vel_floor),
                      max(vel_ci[1], vel_mid + vel_floor)]
    else:
        vel_ci = [float(np.percentile(boot_vel, lo)),
                  float(np.percentile(boot_vel, hi))]

    return {
        'baz_mean': baz_mean,
        'baz_std':  baz_std,
        'baz_ci':   [float(baz_ci[0]), float(baz_ci[1])],
        'vel_mean': float(np.mean(boot_vel)),
        'vel_std':  float(np.std(boot_vel)),
        'vel_ci':   [float(vel_ci[0]), float(vel_ci[1])],
        'n_used':   int(n),
        'block':    block,
        'resolution_limited': resolution_limited,
        'ci':       float(ci),
        'ellipse':  ellipse,
    }


# ─────────────────────────────────────────────────────────────────────────────
#  NEW (v1.2.1): array response function
# ─────────────────────────────────────────────────────────────────────────────

def array_response(st, fmin, fmax, smax=4.0, ds=0.05, n_freq=15):
    """
    Theoretical array response function (ARF) for the current geometry,
    frequency-averaged over [fmin, fmax].  Peaks away from the origin reveal
    spatial aliasing; main-lobe width sets the slowness resolution.

    Returns sx_axis, sy_axis (s/km) and the normalised ARF (peak = 1 at origin).
    """
    positions, _, _ = station_xy_km(st)
    N = positions.shape[0]
    sx_axis = _slowness_grid(smax, ds)
    sy_axis = _slowness_grid(smax, ds)
    freqs = np.linspace(fmin, fmax, max(1, n_freq))

    arf = np.zeros(sy_axis.size * sx_axis.size)
    for f in freqs:
        A = _steering_grid(positions, sx_axis, sy_axis, f)   # (G, N)
        arf += np.abs(A.sum(axis=1) / N) ** 2
    arf /= freqs.size
    arf = arf.reshape(sy_axis.size, sx_axis.size)
    arf /= arf.max() if arf.max() > 0 else 1.0
    return sx_axis, sy_axis, arf
