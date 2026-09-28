"""Focused release checks for archive integrity and CLI input safety."""

import csv
import io
import re
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np
import pandas as pd
from obspy import Stream, Trace, UTCDateTime
from obspy.core.util import AttribDict
from obspy.core.inventory import Channel, Inventory, Network, Station
from obspy.core.inventory.response import Response
from scipy.signal import butter, sosfiltfilt

import results_io
import seismofk_cli
import fk_analysis
import pmcc
from data_readiness import assess_stream
from plot_spectrogram_window import compute_spectrogram


HERE = Path(__file__).resolve().parent


def _write_array(root, offsets_km, input_units=None, gap=None, seconds=40.0):
    """Write a synthetic plane-wave array to *root*; return (mseed, xml) paths.

    *offsets_km* are (east, north) sensor offsets.  With *input_units* each
    channel gets a flat response with that input unit; *gap* is a (t0, t1)
    second range removed from every trace.
    """
    fs = 20.0
    start = UTCDateTime(2025, 1, 1)
    base = np.random.default_rng(12).normal(size=int(seconds * fs))
    stream, stations = Stream(), []
    for i, (east, north) in enumerate(offsets_km):
        code = f"S{i + 1}"
        delay = -(east * 0.8 + north * 2.8)
        tr = Trace(data=np.roll(base, round(delay * fs)))
        tr.stats.update(dict(network="XX", station=code, channel="BDF",
                             sampling_rate=fs, starttime=start))
        if gap is None:
            stream += tr
        else:
            stream += tr.slice(start, start + gap[0])
            stream += tr.slice(start + gap[1], None)
        latitude, longitude = 30 + north / 111.2, 31 + east / 96.3
        response = None
        if input_units:
            response = Response.from_paz(zeros=[], poles=[], stage_gain=100.0,
                                         input_units=input_units,
                                         output_units="COUNTS")
        channel = Channel(code="BDF", location_code="", latitude=latitude,
                          longitude=longitude, elevation=0, depth=0,
                          sample_rate=fs, response=response)
        stations.append(Station(code=code, latitude=latitude,
                                longitude=longitude, elevation=0,
                                creation_date=start, channels=[channel]))
    waveform, inventory = root / "wave.mseed", root / "inventory.xml"
    stream.write(str(waveform), format="MSEED")
    Inventory(networks=[Network(code="XX", stations=stations)],
              source="release test").write(str(inventory), format="STATIONXML")
    return waveform, inventory


SQUARE_KM = ((0, 0), (0.1, 0), (0, 0.1), (0.1, 0.1))


def _plane_wave(seed, baz, vel_ms, snr, burst=(90, 150), n_sta=6,
                aperture_km=1.2, fs=20.0, seconds=240, band=(0.5, 4.0)):
    """In-memory array recording of one band-limited plane wave with
    fractional (Fourier-shifted) delays and independent sensor noise at
    amplitude signal-to-noise ratio *snr*.  Same *seed* = same geometry."""
    rng = np.random.default_rng(seed)
    dx = rng.uniform(-aperture_km / 2, aperture_km / 2, n_sta)
    dy = rng.uniform(-aperture_km / 2, aperture_km / 2, n_sta)
    sx = -np.sin(np.radians(baz)) / (vel_ms / 1e3)
    sy = -np.cos(np.radians(baz)) / (vel_ms / 1e3)
    n = int(seconds * fs)
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    source = np.zeros(n)
    i0, i1 = int(burst[0] * fs), int(burst[1] * fs)
    source[i0:i1] = rng.standard_normal(i1 - i0) * np.hanning(i1 - i0)
    source = sosfiltfilt(sos, source)
    freqs = np.fft.rfftfreq(n, 1 / fs)
    spectrum = np.fft.rfft(source)
    level = source[i0:i1].std()
    stream = Stream()
    for k in range(n_sta):
        tau = sx * dx[k] + sy * dy[k]
        noise = sosfiltfilt(sos, rng.standard_normal(n))
        tr = Trace(data=np.fft.irfft(spectrum * np.exp(-2j * np.pi * freqs * tau), n)
                   + noise * level / snr / noise.std())
        tr.stats.update(dict(network="XX", station=f"S{k}", channel="BDF",
                             sampling_rate=fs, starttime=UTCDateTime(2025, 1, 1)))
        tr.stats.coordinates = AttribDict(
            latitude=30 + dy[k] / 111.195,
            longitude=31 + dx[k] / (111.195 * np.cos(np.radians(30))),
            elevation=0.0)
        stream += tr
    return stream


def _families(stream, stime, etime, window=20.0, hop=5.0):
    pixels = pmcc.pmcc(stream, pmcc.log_bands(0.5, 4, 4), window_sec=window,
                       step_sec=hop, stime=stime, etime=etime, corr_min=0.3,
                       consistency_max=0.2)
    return pmcc.pmcc_families(pixels, time_tol=window,
                              pixel_overlap=window / hop)


class ReleaseChecks(unittest.TestCase):
    def test_readiness_rejects_collinear_inventory(self):
        start = UTCDateTime(2025, 1, 1)
        stream = Stream()
        stations = []
        for index in range(3):
            station_code = f"S{index + 1}"
            trace = Trace(data=np.arange(400, dtype=float))
            trace.stats.network = "XX"
            trace.stats.station = station_code
            trace.stats.channel = "BDF"
            trace.stats.sampling_rate = 20
            trace.stats.starttime = start
            stream += trace
            longitude = 31 + index * 0.001
            channel = Channel(code="BDF", location_code="", latitude=30,
                              longitude=longitude, elevation=0, depth=0,
                              sample_rate=20)
            stations.append(Station(code=station_code, latitude=30,
                                    longitude=longitude, elevation=0,
                                    creation_date=start, channels=[channel]))
        inventory = Inventory(networks=[Network(code="XX", stations=stations)],
                              source="readiness test")
        findings = assess_stream(stream, inventory, fmin=0.5, fmax=6,
                                 start=start, duration=15)
        self.assertIn("Collinear array geometry",
                      {item.title for item in findings if item.severity == "error"})

    def test_readiness_identifies_bad_band_and_missing_metadata(self):
        start = UTCDateTime(2025, 1, 1)
        traces = Stream()
        for station in ("S1", "S2", "S3"):
            tr = Trace(data=np.ones(400))
            tr.stats.network = "XX"
            tr.stats.station = station
            tr.stats.channel = "BDF"
            tr.stats.sampling_rate = 20
            tr.stats.starttime = start
            traces += tr
        findings = assess_stream(traces, None, fmin=0.5, fmax=12,
                                 start=start, duration=15)
        titles = {item.title for item in findings if item.severity == "error"}
        self.assertIn("Band exceeds Nyquist", titles)
        self.assertIn("No StationXML selected", titles)

    def test_spectrogram_counts_fallback_resolves_short_signal(self):
        fs = 20
        trace = Trace(data=100 * np.sin(2 * np.pi * 2 * np.arange(400) / fs))
        trace.stats.sampling_rate = fs
        trace.stats.starttime = UTCDateTime(2025, 1, 1)
        freqs, times, power = compute_spectrogram(
            trace, None, filter_freqmin=0.5, filter_freqmax=6,
            freq_max=6, nperseg=100, noverlap=75, units="counts")
        self.assertGreaterEqual(len(times), 10)
        self.assertAlmostEqual(freqs[np.argmax(power.mean(axis=1))], 2, delta=0.21)
        _, _, smoothed = compute_spectrogram(
            trace, None, filter_freqmin=0.5, filter_freqmax=6,
            freq_max=6, nperseg=100, noverlap=75, units="counts",
            smooth_bins=1.0)
        self.assertEqual(smoothed.shape, power.shape)
        self.assertFalse(np.array_equal(smoothed, power))
        self.assertAlmostEqual(
            freqs[np.argmax(smoothed.mean(axis=1))], 2, delta=0.21)

    def test_inventory_coordinates_use_trace_epoch(self):
        early = Channel(code="BDF", location_code="", latitude=30,
                        longitude=31, elevation=0, depth=0, sample_rate=20,
                        start_date=UTCDateTime(2020, 1, 1),
                        end_date=UTCDateTime(2021, 1, 1))
        late = Channel(code="BDF", location_code="", latitude=31,
                       longitude=32, elevation=0, depth=0, sample_rate=20,
                       start_date=UTCDateTime(2025, 1, 1))
        station = Station(code="S1", latitude=30, longitude=31, elevation=0,
                          creation_date=UTCDateTime(2020, 1, 1),
                          channels=[early, late])
        inventory = Inventory(networks=[Network(code="XX", stations=[station])],
                              source="release test")
        trace = Trace(data=np.ones(200))
        trace.stats.network = "XX"
        trace.stats.station = "S1"
        trace.stats.channel = "BDF"
        trace.stats.starttime = UTCDateTime(2025, 5, 1)
        trace.stats.sampling_rate = 20
        fk_analysis.attach_coordinates(Stream([trace]), inventory, verbose=False)
        self.assertEqual(trace.stats.coordinates.latitude, 31)
        self.assertEqual(trace.stats.coordinates.longitude, 32)

    def test_pmcc_direction_and_velocity_with_staggered_mixed_rate_traces(self):
        start = UTCDateTime(2025, 1, 1)
        baz, velocity = 120.0, 340.0
        rng = np.random.default_rng(4)
        base = sosfiltfilt(
            butter(4, [0.5, 7], btype="bandpass", fs=100, output="sos"),
            rng.normal(size=6000))
        base_time = np.arange(len(base)) / 100
        stream = Stream()
        for east, north, fs, offset in ((0, 0, 20, 0), (0.2, 0, 25, 0.17),
                                        (0, 0.2, 20, 0.11), (0.2, 0.2, 25, 0.26)):
            delay = (-np.sin(np.deg2rad(baz)) * east -
                     np.cos(np.deg2rad(baz)) * north) * 1000 / velocity
            samples = offset + np.arange(int((40 - offset) * fs)) / fs
            data = np.interp(samples - delay, base_time, base)
            trace = Trace(data=data)
            trace.stats.starttime = start + offset
            trace.stats.sampling_rate = fs
            trace.stats.coordinates = AttribDict(
                latitude=30 + north / 111.2,
                longitude=31 + east / 96.3, elevation=0)
            stream += trace
        pixels = pmcc.pmcc(stream, pmcc.log_bands(0.5, 7, 4),
                           window_sec=12, step_sec=6, stime=start,
                           etime=start + 38)
        self.assertGreater(len(pixels), 10)
        direction = pmcc.circular_mean_deg(pixels["baz"], pixels["mean_corr"])
        self.assertLess(abs((direction - baz + 180) % 360 - 180), 5)
        self.assertLess(abs(pixels["app_vel"].median() - velocity), 25)

    def test_pmcc_direction_summary_wraps_at_north(self):
        direction = pmcc.circular_mean_deg([358, 359, 1, 2])
        self.assertLess(abs((direction + 180) % 360 - 180), 2)

    def test_pmcc_rejects_collinear_array_geometry(self):
        stream = Stream()
        for index in range(3):
            trace = Trace(data=np.ones(200))
            trace.stats.starttime = UTCDateTime(2025, 1, 1)
            trace.stats.sampling_rate = 20
            trace.stats.coordinates = AttribDict(
                latitude=30, longitude=31 + index * 0.001, elevation=0)
            stream += trace
        with self.assertRaisesRegex(ValueError, "non-collinear"):
            pmcc.pmcc(stream, pmcc.log_bands(0.5, 4, 2), window_sec=5)

    def test_bootstrap_azimuth_interval_wraps_at_north(self):
        baz = np.array([358.0, 359.0, 1.0, 2.0])
        result = fk_analysis.bootstrap_beam(
            baz, np.full(4, 340.0), np.full(4, 1000 / 340), n_boot=200)
        for endpoint in result["baz_ci"]:
            self.assertLess(abs((endpoint + 180) % 360 - 180), 5)

    def test_cli_archives_synthetic_pmcc_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fs = 20.0
            start = UTCDateTime(2025, 1, 1)
            rng = np.random.default_rng(12)
            base = rng.normal(size=800)
            stream = Stream()
            stations = []
            for i, (east, north) in enumerate(((0, 0), (0.1, 0),
                                                (0, 0.1), (0.1, 0.1))):
                station_code = f"S{i + 1}"
                delay = -(east * 0.8 + north * 2.8)
                tr = Trace(data=np.roll(base, round(delay * fs)))
                tr.stats.network = "XX"
                tr.stats.station = station_code
                tr.stats.channel = "BDF"
                tr.stats.sampling_rate = fs
                tr.stats.starttime = start
                stream += tr
                latitude = 30 + north / 111.2
                longitude = 31 + east / 96.3
                channel = Channel(code="BDF", location_code="",
                                  latitude=latitude, longitude=longitude,
                                  elevation=0, depth=0, sample_rate=fs)
                stations.append(Station(code=station_code, latitude=latitude,
                                        longitude=longitude, elevation=0,
                                        creation_date=start, channels=[channel]))
            waveform = root / "wave.mseed"
            inventory = root / "inventory.xml"
            stream.write(str(waveform), format="MSEED")
            Inventory(networks=[Network(code="XX", stations=stations)],
                      source="release test").write(str(inventory),
                                                    format="STATIONXML")
            out_dir = root / "results"
            args = ["--mseed", str(waveform), "--inventory", str(inventory),
                    "--out-dir", str(out_dir), "--start", str(start),
                    "--end", str(start + 30), "--method", "pmcc",
                    "--fmin", "0.5", "--fmax", "4", "--pmcc-bands", "2",
                    "--win-length", "10", "--corr-min", "0.1",
                    "--consistency", "0.2", "--no-response", "--csv"]
            seismofk_cli.main(args)
            with open(out_dir / "index.csv", newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 1)
            self.assertGreater(int(rows[0]["n_detections"]), 0)
            self.assertTrue((out_dir / rows[0]["detections_path"]).exists())

    def test_index_header_extends_without_losing_rows(self):
        with tempfile.TemporaryDirectory() as out_dir:
            results_io.append_index(out_dir, {"run_id": "first", "method": "capon"})
            results_io.append_index(out_dir, {
                "run_id": "second", "method": "capon", "baz_std": 2.5,
            })
            results_io.append_index(out_dir, {"run_id": "third", "method": "pmcc"})
            with open(f"{out_dir}/index.csv", newline="") as fh:
                reader = csv.DictReader(fh)
                rows = list(reader)
                self.assertIn("baz_std", reader.fieldnames)
                self.assertEqual([r["run_id"] for r in rows],
                                 ["first", "second", "third"])
                self.assertEqual(rows[1]["baz_std"], "2.5")
                self.assertEqual(rows[2]["baz_std"], "")

    def test_repeated_invocations_get_distinct_ids(self):
        start, end = UTCDateTime(2025, 1, 1), UTCDateTime(2025, 1, 2)
        a = seismofk_cli.make_run_tag("I31", start, end, "capon")
        b = seismofk_cli.make_run_tag("I31", start, end, "capon")
        self.assertNotEqual(a, b)
        self.assertTrue(a.startswith("I31_20250101T000000_20250102T000000_capon_"))

    def test_invalid_cli_parameters_fail_before_reading_files(self):
        common = ["--mseed", "missing.mseed", "--inventory", "missing.xml",
                  "--start", "2025-01-01", "--end", "2025-01-02"]
        for invalid in (["--chunk", "0"], ["--overlap", "1"],
                        ["--fmin", "4", "--fmax", "1"],
                        ["--array-name", "../outside"]):
            with self.subTest(invalid=invalid), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    seismofk_cli.main(common + invalid)
                self.assertEqual(error.exception.code, 2)

    def test_calibrated_units_follow_stationxml_input_units(self):
        trace = Trace(data=np.ones(10))
        trace.stats.update(dict(network="XX", station="S1", channel="BDF",
                                starttime=UTCDateTime(2025, 1, 1)))
        with tempfile.TemporaryDirectory() as tmp:
            for input_units, output, expected in (
                    ("PA", "VEL", "Pa"),          # ObsPy applies Pa as-is
                    ("M/S", "VEL", "m/s"),
                    ("M/S**2", "VEL", "m/s"),     # converted to requested output
                    ("M/S", "DEF", "m/s")):
                _, xml = _write_array(Path(tmp), [(0, 0)], input_units)
                inv = fk_analysis.load_inventory(str(xml), verbose=False)
                with self.subTest(units=input_units, output=output):
                    self.assertEqual(
                        fk_analysis.calibrated_units(inv, trace, output), expected)
            _, xml = _write_array(Path(tmp), [(0, 0)])
            inv = fk_analysis.load_inventory(str(xml), verbose=False)
            self.assertEqual(fk_analysis.calibrated_units(inv, trace), "counts")

    def test_cli_records_units_calibration_and_gaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM, "PA",
                                               gap=(24.0, 28.0))
            out_dir = root / "results"
            start = UTCDateTime(2025, 1, 1)
            with redirect_stdout(io.StringIO()):
                seismofk_cli.main([
                    "--mseed", str(waveform), "--inventory", str(inventory),
                    "--out-dir", str(out_dir), "--start", str(start),
                    "--end", str(start + 30), "--method", "pmcc",
                    "--fmin", "0.5", "--fmax", "4", "--pmcc-bands", "2",
                    "--win-length", "10", "--corr-min", "0.1",
                    "--consistency", "0.2", "--csv"])
            with open(out_dir / "index.csv", newline="") as fh:
                row = next(csv.DictReader(fh))
            self.assertEqual(row["units"], "Pa")
            self.assertEqual(row["calibration"], "full_response")
            self.assertEqual(row["n_channels"], "4")
            self.assertAlmostEqual(float(row["gap_fraction"]), 4 / 30, places=2)

    def test_cli_skips_chunks_above_max_gap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM, gap=(5.0, 25.0))
            out_dir = root / "results"
            start = UTCDateTime(2025, 1, 1)
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                seismofk_cli.main([
                    "--mseed", str(waveform), "--inventory", str(inventory),
                    "--out-dir", str(out_dir), "--start", str(start),
                    "--end", str(start + 30), "--method", "pmcc",
                    "--win-length", "10", "--no-response", "--csv"])
            self.assertIn("skipped", stdout.getvalue())
            with open(out_dir / "index.csv", newline="") as fh:
                row = next(csv.DictReader(fh))
            self.assertEqual(row["status"], "gap_skipped")
            self.assertEqual(row["n_detections"], "0")

    def test_cli_preflight_blocks_collinear_array(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(
                root, ((0, 0), (0.1, 0), (0.2, 0)))
            out_dir = root / "results"
            common = ["--mseed", str(waveform), "--inventory", str(inventory),
                      "--out-dir", str(out_dir), "--start", "2025-01-01",
                      "--end", "2025-01-01T00:00:30", "--method", "pmcc",
                      "--win-length", "10", "--no-response"]
            for extra in ([], ["--check-only"]):
                with self.subTest(extra=extra), redirect_stdout(io.StringIO()), \
                        redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as error:
                        seismofk_cli.main(common + extra)
                    self.assertEqual(error.exception.code, 2)
            self.assertFalse(out_dir.exists())

    def test_cli_check_only_passes_good_array(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM)
            with redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    seismofk_cli.main([
                        "--mseed", str(waveform), "--inventory", str(inventory),
                        "--out-dir", str(root / "results"),
                        "--start", "2025-01-01", "--end", "2025-01-01T00:00:30",
                        "--check-only"])
            self.assertEqual(error.exception.code, 0)
            self.assertFalse((root / "results").exists())

    def test_overlapping_epochs_use_latest_calibration(self):
        # IMS-style metadata: a nominal 2017–2050 epoch overlapping a newer
        # calibrated epoch.  ObsPy alone returns the first (nominal) match.
        def epoch(start, sensitivity, latitude):
            response = Response.from_paz(zeros=[], poles=[],
                                         stage_gain=sensitivity,
                                         input_units="PA", output_units="COUNTS")
            return Channel(code="BDF", location_code="", latitude=latitude,
                           longitude=31, elevation=0, depth=0, sample_rate=20,
                           start_date=UTCDateTime(start),
                           end_date=UTCDateTime(2050, 1, 1), response=response)
        station = Station(code="S1", latitude=30, longitude=31, elevation=0,
                          creation_date=UTCDateTime(2017, 1, 1),
                          channels=[epoch(2017, 10000.0, 30.0),
                                    epoch(2022, 7874.02, 30.001)])
        inventory = Inventory(networks=[Network(code="XX", stations=[station])],
                              source="release test")
        trace = Trace(data=np.ones(400))
        trace.stats.update(dict(network="XX", station="S1", channel="BDF",
                                sampling_rate=20,
                                starttime=UTCDateTime(2026, 6, 1)))
        stream = Stream([trace])
        resolved = fk_analysis.resolve_epochs(inventory, stream)
        start = trace.stats.starttime
        self.assertAlmostEqual(resolved.get_response(trace.id, start)
                               .instrument_sensitivity.value, 7874.02)
        self.assertAlmostEqual(
            resolved.get_coordinates(trace.id, start)["latitude"], 30.001)
        findings = assess_stream(stream, inventory, fmin=0.5, fmax=4,
                                 start=start, duration=10)
        self.assertIn("Overlapping StationXML epochs",
                      {item.title for item in findings})

    def test_sensitivity_only_response_is_explained(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, xml = _write_array(Path(tmp), SQUARE_KM[:1], "PA")
            inv = fk_analysis.load_inventory(str(xml), verbose=False)
            trace = Trace(data=np.ones(10))
            trace.stats.update(dict(network="XX", station="S1", channel="BDF",
                                    starttime=UTCDateTime(2025, 1, 1)))
            inv[0][0][0].response.response_stages = []
            with self.assertRaisesRegex(ValueError, "no response stages"):
                fk_analysis.require_response_stages(inv, Stream([trace]))

    def _run_cli(self, argv):
        with redirect_stdout(io.StringIO()):
            seismofk_cli.main(argv)

    def _index_rows(self, out_dir):
        with open(Path(out_dir) / "index.csv", newline="") as fh:
            return list(csv.DictReader(fh))

    def test_chunk_boundaries_do_not_lose_windows(self):
        # The same interval split into 10 s chunks must give the windows of a
        # single 30 s chunk: the window straddling each boundary is kept once.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM)
            common = ["--mseed", str(waveform), "--inventory", str(inventory),
                      "--start", "2025-01-01T00:00:00",
                      "--end", "2025-01-01T00:00:30", "--method", "pmcc",
                      "--fmin", "0.5", "--fmax", "4", "--pmcc-bands", "2",
                      "--win-length", "4", "--step", "2", "--corr-min", "0.1",
                      "--consistency", "0.2", "--no-response", "--csv"]
            times = {}
            for chunk in ("30", "10"):
                out_dir = root / f"chunk{chunk}"
                self._run_cli(common + ["--chunk", chunk,
                                        "--out-dir", str(out_dir)])
                pixels = []
                for row in self._index_rows(out_dir):
                    if row["status"] == "ok":
                        with open(out_dir / row["detections_path"],
                                  newline="") as fh:
                            pixels += [(d["time"], d["f_center"])
                                       for d in csv.DictReader(fh)]
                times[chunk] = sorted(pixels)
            self.assertEqual(times["10"], times["30"])

    def test_cli_saves_chunk_and_summary_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM, "PA")
            for method in ("pmcc", "bartlett", "fk"):
                out_dir = root / method
                self._run_cli([
                    "--mseed", str(waveform), "--inventory", str(inventory),
                    "--out-dir", str(out_dir), "--start", "2025-01-01",
                    "--end", "2025-01-01T00:00:30", "--chunk", "15",
                    "--method", method, "--fmin", "0.5", "--fmax", "4",
                    "--pmcc-bands", "2", "--win-length", "10",
                    "--corr-min", "0.1", "--consistency", "0.2", "--csv",
                    "--save-figures", "--figure-dpi", "40"])
                rows = self._index_rows(out_dir)
                with self.subTest(method=method):
                    self.assertEqual(len(rows), 2)
                    for row in rows:
                        figure = out_dir / row["figure_path"]
                        self.assertTrue(figure.is_file())
                        self.assertGreater(figure.stat().st_size, 1000)
                    self.assertEqual(
                        len(list(out_dir.glob("*_summary.png"))), 1)

    def test_cli_fk_matches_gui_pipeline(self):
        # --method fk runs fk_array like the GUI; archived times are window
        # centres and the median direction is the circular median.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM)
            out_dir = root / "results"
            self._run_cli([
                "--mseed", str(waveform), "--inventory", str(inventory),
                "--out-dir", str(out_dir), "--start", "2025-01-01",
                "--end", "2025-01-01T00:00:30", "--method", "fk",
                "--fmin", "0.5", "--fmax", "4", "--win-length", "5",
                "--no-response", "--csv", "--semb-thresh", "0.1"])
            row = self._index_rows(out_dir)[0]
            self.assertEqual(row["status"], "ok")
            with open(out_dir / row["detections_path"], newline="") as fh:
                first = next(csv.DictReader(fh))
            self.assertEqual(first["time"][11:19], "00:00:02")  # 0 s + 5 s/2
            self.assertIn("semblance", first)

    def test_min_coherence_gates_fk_windows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM)
            out_dir = root / "results"
            self._run_cli([
                "--mseed", str(waveform), "--inventory", str(inventory),
                "--out-dir", str(out_dir), "--start", "2025-01-01",
                "--end", "2025-01-01T00:00:30", "--method", "bartlett",
                "--win-length", "10", "--no-response", "--csv",
                "--min-coherence", "1e9"])
            self.assertEqual(self._index_rows(out_dir)[0]["status"],
                             "no_detections")

    def test_circular_median_wraps_at_north(self):
        self.assertAlmostEqual(pmcc.circular_median_deg([359, 1, 2]), 1.0)
        self.assertAlmostEqual(pmcc.circular_median_deg([350, 355, 10]), 355.0)

    def test_pmcc_families_group_arrivals_and_drop_isolated_pixels(self):
        start = pd.Timestamp("2025-01-01")
        rows = []

        def arrival(t0, baz, vel, n=8):
            for k in range(n):
                band = k % 3
                rows.append(dict(time=start + pd.Timedelta(seconds=t0 + 5 * (k // 3)),
                                 f_min=[0.5, 1.0, 2.0][band],
                                 f_max=[1.0, 2.0, 4.0][band],
                                 f_center=[0.7, 1.4, 2.8][band],
                                 baz=(baz + (k % 2) * 1.5) % 360, app_vel=vel,
                                 slowness=1e3 / vel, consistency=0.01,
                                 mean_corr=0.9, n_pairs=6))

        arrival(100, 120.0, 340.0)            # arrival 1
        arrival(400, 359.0, 330.0)            # arrival 2, wraps at north
        arrival(400, 250.0, 340.0, n=2)       # too small: isolated pixels
        arrival(100, 120.0, 520.0, n=3)       # same time/azimuth, other velocity
        pixels = pd.DataFrame(rows)
        out, families = pmcc.pmcc_families(pixels, time_tol=10, min_pixels=5)
        self.assertEqual(len(families), 2)
        first, second = families.itertuples()
        self.assertAlmostEqual(first.baz, 120.75, delta=0.5)
        self.assertEqual(first.n_pixels, 8)
        self.assertLess(abs((second.baz + 180) % 360 - 180), 1.5)   # ~0°
        self.assertEqual(int((out["family"] == 0).sum()), 5)

    def test_pmcc_uncertainty_covers_truth(self):
        # Coverage check: 95 % intervals must contain the true back-azimuth
        # and velocity in about 95 % of noisy synthetic arrivals, for pixels
        # and for families, without being grossly over-conservative.
        t0 = UTCDateTime(2025, 1, 1)
        rng = np.random.default_rng(2026)
        pixel_z, family_hits, family_z = [], [], []
        for trial in range(16):
            baz, vel = rng.uniform(0, 360), rng.uniform(300, 380)
            pixels, families = _families(
                _plane_wave(trial, baz, vel, snr=2.0), t0 + 70, t0 + 170)
            dbaz = (pixels["baz"] - baz + 180) % 360 - 180
            pixel_z += list(np.abs(dbaz / pixels["baz_std"]))
            fam = families.loc[families["n_pixels"].idxmax()]
            fdbaz = abs((fam.baz - baz + 180) % 360 - 180)
            family_hits.append(fdbaz <= fam.baz_ci95
                               and abs(fam.app_vel - vel) <= fam.vel_ci95)
            family_z.append(fdbaz / (fam.baz_ci95 / 1.96))
            self.assertFalse(fam.mixed)
        self.assertGreaterEqual(np.mean(np.array(pixel_z) < 1.96), 0.90)
        self.assertGreaterEqual(np.mean(family_hits), 0.85)
        # Not over-conservative: median |z| near 0.67 for calibrated errors.
        self.assertGreater(np.median(pixel_z), 0.35)
        self.assertGreater(np.median(family_z), 0.30)

    def test_pmcc_family_flags_merged_sources(self):
        # Two back-to-back arrivals 8° apart merge into one family and must
        # be flagged; the same pair from one direction must not be.
        t0 = UTCDateTime(2025, 1, 1)
        for separation, expected in ((0.0, False), (8.0, True)):
            first = _plane_wave(3, 100.0, 340.0, 5.0, burst=(90, 125))
            second = _plane_wave(3, 100.0 + separation, 340.0, 5.0,
                                 burst=(120, 155))
            for a, b in zip(first, second):
                a.data = a.data + b.data
            _, families = _families(first, t0 + 70, t0 + 175)
            fam = families.loc[families["n_pixels"].idxmax()]
            with self.subTest(separation=separation):
                self.assertEqual(bool(fam.mixed), expected)
                if expected:
                    self.assertGreater(fam.mixed_dbaz, 3.0)

    def test_cli_archives_pmcc_families(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM)
            out_dir = root / "results"
            self._run_cli([
                "--mseed", str(waveform), "--inventory", str(inventory),
                "--out-dir", str(out_dir), "--start", "2025-01-01",
                "--end", "2025-01-01T00:00:30", "--method", "pmcc",
                "--fmin", "0.5", "--fmax", "4", "--pmcc-bands", "2",
                "--win-length", "10", "--corr-min", "0.1",
                "--consistency", "0.2", "--no-response", "--csv",
                "--family-min-pixels", "3"])
            row = self._index_rows(out_dir)[0]
            self.assertGreaterEqual(int(row["n_families"]), 1)
            with open(out_dir / row["families_path"], newline="") as fh:
                family = next(csv.DictReader(fh))
            self.assertEqual(family["run_id"], row["run_id"])
            self.assertGreaterEqual(int(family["n_pixels"]), 3)

    def _sensor(self, data, station="S1", fs=20.0):
        tr = Trace(data=np.asarray(data, float))
        tr.stats.update(dict(network="XX", station=station, channel="BDF",
                             sampling_rate=fs,
                             starttime=UTCDateTime(2025, 1, 1)))
        return tr

    def test_noise_level_of_known_sine_in_db_re_20upa(self):
        import noise_levels as nl
        fs, seconds = 20.0, 600
        t = np.arange(int(seconds * fs)) / fs
        levels = nl.noise_levels(Stream([self._sensor(np.sin(2 * np.pi * t))]),
                                 60, 0.5, 4.0, units="Pa")
        # 1 Pa amplitude → RMS 1/√2 Pa → 20·log10(0.7071 / 20e-6) = 90.97 dB
        self.assertEqual(len(levels), 10)             # absolute 1-min grid
        np.testing.assert_allclose(levels["level_db"], 90.97, atol=0.05)

    def test_noise_leq_averages_power_and_gaps_are_excluded(self):
        import noise_levels as nl
        fs = 20.0
        n = int(60 * fs)
        rng = np.random.default_rng(1)
        quiet = rng.standard_normal(n) * 0.002         # ~40 dB
        loud = rng.standard_normal(n) * 0.02           # ~60 dB
        gap = np.concatenate([rng.standard_normal(n // 4) * 0.002,
                              np.zeros(n - n // 4)])   # 75 % zero-filled
        stream = Stream([self._sensor(np.concatenate([quiet, loud, gap]))])
        levels = nl.noise_levels(stream, 60, units="Pa")
        self.assertTrue(np.isnan(levels["level_db"].iloc[2]))
        self.assertAlmostEqual(levels["coverage"].iloc[2], 0.25, places=2)
        stats = nl.noise_statistics(levels, units="Pa").iloc[0]
        db = levels["level_db"].iloc[:2].to_numpy()
        expected_leq = 10 * np.log10(np.mean(10 ** (db / 10)))
        self.assertAlmostEqual(stats.Leq, expected_leq, places=2)
        self.assertGreater(stats.Leq, np.mean(db) + 2.5)  # not a dB mean

    def test_noise_flags_sensor_offsets(self):
        import noise_levels as nl
        rng = np.random.default_rng(5)
        common = rng.standard_normal(int(1800 * 20.0)) * 0.01
        gains = {"S1": 1.0, "S2": 1.0, "S3": 1.0, "S4": 1.0,
                 "S5": 1.5, "S6": 2.0}              # +3.5 dB, +6.0 dB
        stream = Stream([
            self._sensor(common * g + rng.standard_normal(common.size) * 1e-4,
                         station=name) for name, g in gains.items()])
        stats = nl.noise_statistics(nl.noise_levels(stream, 60, units="Pa"),
                                    units="Pa").set_index("station")
        self.assertEqual(stats.loc["S1", "flag"], "ok")
        self.assertEqual(stats.loc["S5", "flag"], "warning")
        self.assertEqual(stats.loc["S6", "flag"], "fault")
        self.assertAlmostEqual(stats.loc["S6", "offset_db"], 6.02, delta=0.1)

    def test_cli_noise_levels_do_not_depend_on_chunking(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            waveform, inventory = _write_array(root, SQUARE_KM, "PA",
                                               seconds=60.0)
            tables = {}
            for chunk in ("60", "20"):
                out_dir = root / f"chunk{chunk}"
                self._run_cli([
                    "--mseed", str(waveform), "--inventory", str(inventory),
                    "--out-dir", str(out_dir), "--start", "2025-01-01",
                    "--end", "2025-01-01T00:01:00", "--method", "noise",
                    "--noise-window", "5", "--chunk", chunk, "--csv",
                    "--fmin", "0.5", "--fmax", "4"])
                rows = self._index_rows(out_dir)
                self.assertTrue(all(r["units"] == "Pa" for r in rows))
                frames = [pd.read_csv(out_dir / r["detections_path"])
                          for r in rows if r["status"] == "ok"]
                tables[chunk] = pd.concat(frames).sort_values(
                    ["time", "trace_id"]).reset_index(drop=True)
                self.assertEqual(len(list(out_dir.glob("*_noise_stats.csv"))), 1)
            self.assertEqual(len(tables["60"]), len(tables["20"]))
            np.testing.assert_allclose(tables["20"]["level_db"],
                                       tables["60"]["level_db"], atol=0.3)

    def test_block_bootstrap_widens_interval_for_correlated_windows(self):
        # Heavily overlapping windows give smoothly varying estimates; an
        # i.i.d. bootstrap treats them as independent and is far too narrow.
        rng = np.random.default_rng(3)
        slow_drift = np.convolve(rng.standard_normal(400), np.ones(10) / 10,
                                 mode="same")
        baz = 120.0 + 3.0 * slow_drift
        vel = np.full(baz.size, 340.0)
        slowness = np.full(baz.size, 1e3 / 340.0)
        iid = fk_analysis.bootstrap_beam(baz, vel, slowness, n_boot=400)
        blocked = fk_analysis.bootstrap_beam(baz, vel, slowness, n_boot=400,
                                             block=10)
        width = lambda b: (b["baz_ci"][1] - b["baz_ci"][0]) % 360
        self.assertEqual(blocked["block"], 10)
        self.assertGreater(width(blocked), 2.0 * width(iid))

    def test_fk_figure_reports_percentile_interval_not_std(self):
        from matplotlib.figure import Figure
        import result_plots

        n = 60
        rng = np.random.default_rng(4)
        r = {
            "time": pd.date_range("2025-01-01", periods=n, freq="2s"),
            "semblance": np.full(n, 0.8), "fisher": np.full(n, 20.0),
            "bazi": 120 + 0.05 * rng.standard_normal(n),
            "app_vel": 340 + 0.2 * rng.standard_normal(n),
            "slowness": np.full(n, 1e3 / 340.0), "expected_bazi": 120.0,
            "beam_waveform": np.zeros(10),
            "beam_times": pd.date_range("2025-01-01", periods=10, freq="12s"),
            "n_stations": 4, "station_x_km": np.array([0, 0.1, 0, 0.1]),
            "station_y_km": np.array([0, 0, 0.1, 0.1]),
            "station_ids": ["S1", "S2", "S3", "S4"], "beam_units": "Pa"}
        fig = Figure(figsize=(15, 11))
        result_plots.fk_results_figure(
            fig, r, {"bootstrap": True, "overlap": 0.1}, event_name="T",
            fmin=0.5, fmax=4, semb_thresh=0.3)
        title = fig._suptitle.get_text()
        self.assertIn("95% CI", title)
        self.assertIn("blocks of 10", title)
        self.assertRegex(title, r"baz 1\d\d\.\d° \[\d")
        self.assertNotIn("± 0", title)
        self.assertNotIn("±0", title)

    def test_version_strings_agree(self):
        pyproject = (HERE / "pyproject.toml").read_text()
        citation = (HERE / "CITATION.cff").read_text()
        version = fk_analysis.__version__
        self.assertIn(f'version = "{version}"', pyproject)
        self.assertRegex(citation, rf'(?m)^version: "{re.escape(version)}"$')


if __name__ == "__main__":
    unittest.main()
