"""
fk_analysis.py — Core FK / beamforming routines for the SeismoFK GUI.

Copyright (c) 2024-2025 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.

Acknowledgements:
    Array processing methodology adapted in part from:
    Jelle Assink (jelle.assink@knmi.nl)
    https://github.com/roseseismo/roses2021/blob/main/unit08/array_processing.py

"""

import numpy as np
import pandas as pd
import os

from obspy import read_inventory, UTCDateTime, geodetics
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import array_processing


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

def get_coordinates_safe(inv, trace_id, verbose=True):
    """
    Return a coordinate dict for *trace_id*, trying several location-code
    fallbacks and ultimately station-level coordinates.
    """
    try:
        return inv.get_coordinates(trace_id)
    except Exception:
        pass

    net, sta, loc, cha = trace_id.split('.')
    for loc_try in ['', '--', '00', '10']:
        try:
            coords = inv.get_coordinates(f"{net}.{sta}.{loc_try}.{cha}")
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


# ─────────────────────────────────────────────────────────────────────────────
#  Beamforming
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
        dlat = np.deg2rad(tr.stats.coordinates.latitude  - ref_lat) * 6371000
        dlon = (np.deg2rad(tr.stats.coordinates.longitude - ref_lon) * 6371000 *
                np.cos(np.deg2rad(ref_lat)))
        delay_samp = int(round((sx * dlon + sy * dlat) / dt))
        shifted = np.roll(tr.data, -delay_samp)
        if delay_samp > 0:
            shifted[-delay_samp:] = 0.0
        elif delay_samp < 0:
            shifted[:-delay_samp] = 0.0
        beam += shifted

    beam /= len(a)
    return beam, a[0].times('matplotlib')


# ─────────────────────────────────────────────────────────────────────────────
#  Main FK routine
# ─────────────────────────────────────────────────────────────────────────────

def fk_array(st, inv, freq_min, freq_max, win_length, overlap,
             start_t, duration, arr_name,
             source_lat, source_lon,
             verbose=True):
    """
    Run frequency-wavenumber array processing on *st*.

    Parameters
    ----------
    st            : obspy Stream (sensitivity-corrected, e.g. in Pa)
    inv           : inventory path or obspy Inventory object
    freq_min/max  : bandpass corner frequencies (Hz)
    win_length    : FK window length (s)
    overlap       : fractional window overlap (0 < overlap < 1)
    start_t       : analysis start time (ISO string)
    duration      : analysis duration (s)
    arr_name      : label for CSV filename
    source_lat/lon: expected source location for reference back-azimuth
    verbose       : print progress messages

    Returns
    -------
    dict with keys: time, fisher, semblance, bazi, app_vel, slowness,
                    expected_bazi, beam_waveform, beam_times,
                    waveform_times, waveform_data, station_name,
                    med_baz, med_vel, n_stations,
                    station_x_km, station_y_km, station_ids
    """
    # ── 1. Load inventory and attach coordinates ──────────────────────────
    inv = load_inventory(inv, verbose=verbose)

    for tr in st:
        coords = get_coordinates_safe(inv, tr.id, verbose=verbose)
        tr.stats.coordinates = AttribDict({
            'latitude':  coords['latitude'],
            'elevation': coords['elevation'],
            'longitude': coords['longitude'],
        })

    # ── 2. Expected back-azimuth ──────────────────────────────────────────
    bazz = geodetics.gps2dist_azimuth(
        source_lat, source_lon,
        st[0].stats.coordinates.latitude,
        st[0].stats.coordinates.longitude,
        a=6378137.0, f=0.0033528106647474805
    )[2]
    if verbose:
        print(f'[INFO] Expected Back-Azimuth >>> {bazz:.2f}°')

    # ── 3. Trim, detrend, filter ──────────────────────────────────────────
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

    # Resample to common rate if traces differ (e.g. HL2 at 50 Hz vs others at 100 Hz)
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

    # ── 4. array_processing ───────────────────────────────────────────────
    # Clamp etime to actual trace endtime (resampling can shift it by ±1 sample)
    #Obspy Array Function
    actual_etime = min(tr.stats.endtime for tr in a)
    if actual_etime < etime:
        etime = actual_etime

    out = array_processing(
        a,
        win_len=win_length,
        win_frac=overlap,
        frqlow=freq_min,
        frqhigh=freq_max,
        prewhiten=0,
        sll_x=-1e3 / 250.0,  slm_x=1e3 / 250.0,
        sll_y=-1e3 / 250.0,  slm_y=1e3 / 250.0,
        sl_s=(1e3 / 250.0) / 25,
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

    # ── 5. Unpack results ─────────────────────────────────────────────────
    n_instr   = len(a)
    semblance = out[:, 1]
    fk_power  = out[:, 2]
    bazi      = out[:, 3] % 360
    slowness  = out[:, 4]           # s/km
    app_vel   = 1e3 / slowness      # m/s

    # Fisher ratio — guard against divide-by-zero
    fisher = (n_instr - 1) * semblance / (1.0 - semblance + 1e-12)

    # ── 6. Median beam direction (absolute semblance threshold) ──────────
    #The threshold value of the best beam = 0.35....
    ABS_SEMB_THRESH = 0.35
    mask = semblance >= ABS_SEMB_THRESH
    if mask.sum() == 0:          # fallback: nothing clears the bar → use top-25 %
        mask = semblance >= np.percentile(semblance, 75)
        if verbose:
            print('[WARN] No windows above semblance=0.35; falling back to top-25 % percentile')
    med_baz = float(np.median(bazi[mask]))
    med_vel = float(np.median(app_vel[mask]))
    if verbose:
        print(f'[INFO] Median beam → baz={med_baz:.1f}°, app_vel={med_vel:.0f} m/s '
              f'({mask.sum()}/{len(semblance)} windows used)')

    # ── 7. Delay-and-sum beam ─────────────────────────────────────────────
    beam_wave, beam_times = compute_beam(
        a, med_baz, med_vel, stime, actual_etime, verbose=verbose
    )

    # ── 8. Station geometry (offsets from centroid in km) ─────────────────
    lats    = np.array([tr.stats.coordinates.latitude  for tr in a])
    lons    = np.array([tr.stats.coordinates.longitude for tr in a])
    sta_ids = [tr.stats.station for tr in a]

    ref_lat = float(np.mean(lats))
    ref_lon = float(np.mean(lons))

    station_x_km = (np.deg2rad(lons - ref_lon) * 6371.0 *
                    np.cos(np.deg2rad(ref_lat)))
    station_y_km = np.deg2rad(lats - ref_lat) * 6371.0

    # ── 9. Save CSV ───────────────────────────────────────────────────────
    csv_file = f"Output_{arr_name}.csv"
    outfile  = pd.DataFrame({
        'time':      pd.to_datetime(out[:, 0], unit='s'),
        'semblance': semblance,
        'fk_power':  fk_power,
        'fisher':    fisher,
        'bazi':      bazi,
        'app_vel':   app_vel,
        'slowness':  slowness,
    })
    outfile.to_csv(csv_file, index=False)
    if verbose:
        print(f'[INFO] CSV saved → {csv_file}')

    # ── 10. Return ────────────────────────────────────────────────────────
    return {
        'time':           outfile['time'],
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
