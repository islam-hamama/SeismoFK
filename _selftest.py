"""Synthetic plane-wave self-test for the v1.2.1 array-analysis modules."""
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from obspy.core.util import AttribDict

import fk_analysis as fk
import pmcc as pmcc_mod


def check_direction(baz, vel, label, baz_tol=5.0, vel_tol=40.0):
    """Fail the self-test if a method misses the known synthetic wave."""
    baz_error = abs((baz - true_baz + 180.0) % 360.0 - 180.0)
    assert baz_error <= baz_tol, f"{label}: back-azimuth error {baz_error:.1f}°"
    assert abs(vel - true_vel * 1000) <= vel_tol, (
        f"{label}: apparent velocity {vel:.1f} m/s")

np.random.seed(0)

# ── Build a 6-element array (coords in deg around a centre) ──────────────────
lat0, lon0 = 30.0, 31.0
fs = 50.0
true_baz = 120.0          # degrees
true_vel = 0.34           # km/s  (=340 m/s)
fmin, fmax = 1.0, 8.0

# random small aperture array (~1 km)
rng = np.random.default_rng(1)
dx_km = rng.uniform(-0.6, 0.6, 6)
dy_km = rng.uniform(-0.6, 0.6, 6)
R = 6371.0
dlon = dx_km / (R * np.cos(np.deg2rad(lat0))) * 180 / np.pi
dlat = dy_km / R * 180 / np.pi
lats = lat0 + dlat
lons = lon0 + dlon

# slowness vector (s/km) for the true plane wave
sx = -np.sin(np.deg2rad(true_baz)) / true_vel
sy = -np.cos(np.deg2rad(true_baz)) / true_vel

npts = int(120 * fs)
t = np.arange(npts) / fs
# broadband wavelet (filtered noise burst) centred at 60 s
base = np.zeros(npts)
burst = (t > 40) & (t < 80)
base[burst] = np.random.randn(burst.sum())

start = UTCDateTime(2025, 5, 1)
st = Stream()
for i in range(6):
    tau = sx * dx_km[i] + sy * dy_km[i]          # arrival delay (s) = s·r
    shift = int(round(tau * fs))
    data = np.roll(base, shift) + 0.05 * np.random.randn(npts)
    tr = Trace(data=data)
    tr.stats.sampling_rate = fs
    tr.stats.starttime = start
    tr.stats.network = "XX"
    tr.stats.station = f"S{i+1}"
    tr.stats.channel = "BDF"
    tr.stats.coordinates = AttribDict(
        {"latitude": lats[i], "longitude": lons[i], "elevation": 0.0})
    st += tr

print(f"True baz={true_baz}°, vel={true_vel*1000:.0f} m/s\n")

# ── 1. High-resolution beamformers ───────────────────────────────────────────
for method in ("bartlett", "capon", "music"):
    win = st.copy().trim(start + 45, start + 75)
    win.detrend("demean")
    win.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
    sxa, sya, P = fk.fk_beamform(win, fmin, fmax, method=method, smax=4, ds=0.05)
    pk = fk.estimate_slowness_peak(sxa, sya, P)
    print(f"{method:9s}: baz={pk['baz']:6.1f}°  vel={pk['app_vel']:6.0f} m/s  "
          f"coh={pk['coherence']:.1f}")
    check_direction(pk['baz'], pk['app_vel'], method)

# ── 2. Array response function ───────────────────────────────────────────────
sxa, sya, arf = fk.array_response(st, fmin, fmax, smax=4, ds=0.1)
print(f"\nARF: peak={arf.max():.2f} at origin, shape={arf.shape}")
assert np.isclose(arf.max(), 1.0) and np.isclose(arf[40, 40], 1.0)

# ── 3. fk_scan windowed ──────────────────────────────────────────────────────
df = fk.fk_scan(st, fmin, fmax, win_length=20, overlap=0.5, method="capon",
                stime=start + 40, etime=start + 85)
print(f"\nfk_scan: {len(df)} windows; median baz={df['baz'].median():.1f}°, "
      f"vel={df['app_vel'].median():.0f} m/s")
assert len(df) >= 2, "fk_scan produced too few windows"
check_direction(df['baz'].median(), df['app_vel'].median(), "fk_scan")

# ── 4. bootstrap ─────────────────────────────────────────────────────────────
if len(df) >= 2:
    boot = fk.bootstrap_beam(df["baz"].to_numpy(), df["app_vel"].to_numpy(),
                             df["slowness"].to_numpy(), n_boot=500,
                             slowness_step=0.1)          # fk_scan default ds
    print(f"bootstrap: baz={boot['baz_mean']:.1f}° "
          f"[{boot['baz_ci'][0]:.2f}, {boot['baz_ci'][1]:.2f}] (95% CI)  "
          f"vel={boot['vel_mean']:.0f} m/s "
          f"[{boot['vel_ci'][0]:.1f}, {boot['vel_ci'][1]:.1f}]  "
          f"ellipse a={boot['ellipse']['a']:.3f} b={boot['ellipse']['b']:.3f} s/km")
    check_direction(boot['baz_mean'], boot['vel_mean'], "bootstrap")
    b_lo, b_hi = boot['baz_ci']
    assert (true_baz - b_lo) % 360 <= (b_hi - b_lo) % 360, \
        "bootstrap 95% back-azimuth interval excludes the true value"
    assert boot['vel_ci'][0] <= true_vel * 1000 <= boot['vel_ci'][1], \
        "bootstrap 95% velocity interval excludes the true value"
    assert boot['n_used'] == len(df)

# ── 5. PMCC ──────────────────────────────────────────────────────────────────
bands = pmcc_mod.log_bands(fmin, fmax, 4)
pdf = pmcc_mod.pmcc(st, bands, window_sec=20, step_sec=10,
                    consistency_max=0.1, corr_min=0.4,
                    min_app_vel=200, max_app_vel=600,
                    stime=start + 40, etime=start + 85)
print(f"\nPMCC: {len(pdf)} pixels")
assert len(pdf) > 0, "PMCC produced no detection pixels"
if len(pdf):
    print(f"      median baz={pdf['baz'].median():.1f}°, "
          f"vel={pdf['app_vel'].median():.0f} m/s, "
          f"consistency={pdf['consistency'].median():.4f}s")
    check_direction(pdf['baz'].median(), pdf['app_vel'].median(), "PMCC")

# ── 6. PMCC families: the single synthetic arrival forms one family ─────────
_, fam = pmcc_mod.pmcc_families(pdf, time_tol=20, min_pixels=3,
                                pixel_overlap=2.0)   # window 20 s / hop 10 s
print(f"\nPMCC families: {len(fam)}")
assert len(fam) == 1, f"expected one PMCC family, got {len(fam)}"
print(f"      F1 baz={fam.baz[0]:.1f}° ± {fam.baz_ci95[0]:.2f}° (95%), "
      f"vel={fam.app_vel[0]:.0f} ± {fam.vel_ci95[0]:.1f} m/s, "
      f"{fam.n_pixels[0]} pixels, mixed={fam.mixed[0]}")
check_direction(fam.baz[0], fam.app_vel[0], "PMCC family")
assert not fam.mixed[0], "a single plane wave was flagged as mixed sources"

print("\nSELF-TEST OK")
