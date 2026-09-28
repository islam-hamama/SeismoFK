"""
noise_levels.py — time-domain RMS noise levels for long array records.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.

For every sensor the band-passed pressure is cut into consecutive windows
(seconds to an hour) on an absolute UTC grid, and each window's RMS is
expressed as a sound pressure level,

    L = 20 · log10(p_rms / p_ref),   p_ref = 20 µPa  (dB re 20 µPa)

when the data are calibrated in Pa (otherwise dB re 1 count, clearly
labelled).  Summary statistics follow acoustic practice:

* ``L90`` / ``L50`` / ``L10`` — levels exceeded 90 / 50 / 10 % of the time
  (L90 = background floor, L10 = windy periods and events);
* ``Leq`` — equivalent continuous level, an average of *power* (p²), never
  of dB values;
* offset of each sensor from the median of the *other* sensors, per window,
  summarised by its median and MAD (dB).  Leaving the sensor out matters: a
  faulty sensor would otherwise pull the array median toward itself and
  hide part of its own error.  A persistent offset of ≥ ``warn_db`` (3 dB) is a
  warning and ≥ ``fault_db`` (6 dB, a factor of 2 in pressure) a fault —
  gain / metadata errors, wind-noise-reducer faults or dead sensors.

The module is Qt-free and shared by the GUI noise window and
``seismofk_cli --method noise``.
"""

import numpy as np
import pandas as pd

P_REF_PA = 20e-6

LEVEL_COLUMNS = ['time', 'trace_id', 'station', 'rms', 'peak', 'level_db',
                 'coverage']
STAT_COLUMNS = ['trace_id', 'station', 'n_windows', 'L90', 'L50', 'L10',
                'Leq', 'min_db', 'max_db', 'std_db', 'offset_db',
                'offset_mad_db', 'flag']


def reference(units):
    """(reference value, label) for levels in *units*."""
    if units == 'Pa':
        return P_REF_PA, 'dB re 20 µPa'
    return 1.0, f'dB re 1 {units or "count"}'


def _gap_mask(data, fs, min_gap_sec=1.0):
    """True where samples belong to a zero-filled gap: a run of exact zeros
    lasting at least ``min_gap_sec`` (how SeismoFK fills merge gaps)."""
    zero = np.asarray(data) == 0
    mask = np.zeros(zero.size, dtype=bool)
    if not zero.any():
        return mask
    edges = np.flatnonzero(np.diff(np.concatenate(([0], zero.view(np.int8),
                                                   [0]))))
    min_len = max(int(min_gap_sec * fs), 2)
    for start, stop in zip(edges[::2], edges[1::2]):
        if stop - start >= min_len:
            mask[start:stop] = True
    return mask


def _zero_phase_filter(data, fs, fmin, fmax):
    """4-pole Butterworth, forward-backward, with odd (reflected) padding of
    a few periods of the lowest frequency instead of a taper, so the first
    and last windows of a record keep their full energy."""
    from scipy.signal import butter, sosfiltfilt

    nyquist = fs / 2.0
    fmax = min(fmax, 0.95 * nyquist) if fmax else None
    if fmin and fmax:
        sos = butter(4, [fmin, fmax], btype='bandpass', fs=fs, output='sos')
    elif fmin:
        sos = butter(4, fmin, btype='highpass', fs=fs, output='sos')
    elif fmax:
        sos = butter(4, fmax, btype='lowpass', fs=fs, output='sos')
    else:
        return data
    lowest = fmin or fmax
    padlen = min(data.size - 1, int(10.0 * fs / lowest))
    return sosfiltfilt(sos, data, padtype='odd', padlen=padlen)


def noise_levels(st, window_sec, fmin=None, fmax=None, t0=None, t1=None,
                 units='Pa', min_coverage=0.9):
    """
    RMS level of every sensor in consecutive windows of ``window_sec``.

    Windows lie on an absolute grid (multiples of ``window_sec`` since the
    epoch) and only complete windows inside [t0, t1] are returned.  The data
    are band-passed to [fmin, fmax] Hz (both None = broadband, demeaned
    only).  Samples in zero-filled gaps are excluded; windows with less than
    ``min_coverage`` real data get NaN levels.

    Returns a DataFrame with one row per window and trace: ``time`` (window
    start, UTC), ``trace_id``, ``station``, ``rms`` and ``peak`` (in
    *units*), ``level_db`` (see :func:`reference`) and ``coverage``.
    """
    if window_sec <= 0:
        raise ValueError("window_sec must be positive")
    p_ref, _ = reference(units)
    rows = []
    for tr in st:
        fs = float(tr.stats.sampling_rate)
        gaps = _gap_mask(tr.data, fs)
        work = tr.copy()
        work.data = work.data.astype(float)
        if gaps.any():
            real = work.data[~gaps]
            work.data[gaps] = real.mean() if real.size else 0.0
        work.detrend('demean')
        data = _zero_phase_filter(work.data, fs, fmin, fmax)
        # The last sample covers one more sample interval: a 600 s record at
        # 20 Hz ends at 599.95 s but holds ten complete 60 s windows.
        data_end = tr.stats.endtime + tr.stats.delta
        start = max(tr.stats.starttime, t0) if t0 is not None else tr.stats.starttime
        end = min(data_end, t1) if t1 is not None else data_end
        first = np.ceil(start.timestamp / window_sec) * window_sec
        n_win = int(np.floor((end.timestamp - first) / window_sec + 1e-9))
        per_win = int(round(window_sec * fs))
        for k in range(n_win):
            w0 = first + k * window_sec
            i0 = int(round((w0 - tr.stats.starttime.timestamp) * fs))
            i1 = i0 + per_win
            if i0 < 0 or i1 > data.size:
                continue
            good = ~gaps[i0:i1]
            coverage = float(good.mean())
            if coverage >= min_coverage:
                seg = data[i0:i1][good]
                rms = float(np.sqrt(np.mean(seg ** 2)))
                peak = float(np.max(np.abs(seg)))
                level = 20.0 * np.log10(rms / p_ref) if rms > 0 else np.nan
            else:
                rms = peak = level = np.nan
            rows.append((pd.Timestamp(w0, unit='s'), tr.id,
                         tr.stats.station, rms, peak, level, coverage))
    return pd.DataFrame(rows, columns=LEVEL_COLUMNS)


def noise_statistics(levels, warn_db=3.0, fault_db=6.0, units='Pa'):
    """
    Per-sensor statistics of a :func:`noise_levels` table.

    ``offset_db`` is the median, over windows, of the sensor level minus the
    median level of the other sensors in the same window; ``offset_mad_db`` is its median
    absolute deviation (a steady offset suggests calibration, a varying one
    wind exposure).  ``flag`` is 'fault' if |offset| ≥ ``fault_db``,
    'warning' if ≥ ``warn_db``, otherwise 'ok'.  Leq averages power.
    """
    p_ref, _ = reference(units)
    valid = levels.dropna(subset=['level_db'])
    if valid.empty:
        return pd.DataFrame(columns=STAT_COLUMNS)
    grid = valid.pivot_table(index='time', columns='trace_id',
                             values='level_db')
    rows = []
    for trace_id, lv in valid.groupby('trace_id'):
        db = lv['level_db'].to_numpy()
        rms = lv['rms'].to_numpy()
        others = grid.drop(columns=trace_id).median(axis=1)
        offsets = (grid[trace_id] - others).dropna().to_numpy() \
            if grid.shape[1] > 1 else np.zeros(1)
        if offsets.size == 0:
            offsets = np.zeros(1)
        offset = float(np.median(offsets))
        flag = ('fault' if abs(offset) >= fault_db else
                'warning' if abs(offset) >= warn_db else 'ok')
        rows.append({
            'trace_id': trace_id, 'station': lv['station'].iloc[0],
            'n_windows': int(db.size),
            'L90': float(np.percentile(db, 10)),
            'L50': float(np.percentile(db, 50)),
            'L10': float(np.percentile(db, 90)),
            'Leq': float(10.0 * np.log10(np.mean(rms ** 2) / p_ref ** 2)),
            'min_db': float(db.min()), 'max_db': float(db.max()),
            'std_db': float(db.std()),
            'offset_db': offset,
            'offset_mad_db': float(np.median(np.abs(offsets - offset))),
            'flag': flag,
        })
    return pd.DataFrame(rows, columns=STAT_COLUMNS)


def hourly_levels(levels):
    """Median level per hour of day (UTC) and trace — the diurnal (wind)
    cycle.  Returns a DataFrame indexed 0–23 with one column per trace."""
    valid = levels.dropna(subset=['level_db'])
    if valid.empty:
        return pd.DataFrame(index=range(24))
    hours = pd.to_datetime(valid['time']).dt.hour
    return (valid.assign(hour=hours)
            .pivot_table(index='hour', columns='trace_id', values='level_db',
                         aggfunc='median')
            .reindex(range(24)))
