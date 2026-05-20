#!/usr/bin/env python3
"""
convert_ims.py — Convert all_IMS_sts.xml → one XML file per IMS station in XML_IM/

Copyright (c) 2024-2025 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.

Output format matches XML/I57US.xml reference.

Rules:
  - Take ALL BDF channels per sensor element
  - Group by sensor name (e.g. I57H1, I57H2 …) → each becomes a <Station>
  - Normalise sensor description:
      already starts with array code  → keep as-is
      otherwise                       → "{array_code} {FIRST_TWO_WORDS_UPPER}"
  - <Site><Name> = description + " " + element_suffix  (sensor_name[3:])
  - Station lat/lon/elev  = first channel's coordinates
  - Station startDate/endDate = min/max of all channel dates
  - Every channel: locationCode="" , add <Azimuth>/<Dip> if absent, reorder
    <InstrumentSensitivity> children as: Value, Frequency, InputUnits, OutputUnits
"""

import argparse
import xml.etree.ElementTree as ET
import os
import re
from datetime import datetime

# Portable defaults — resolved relative to this script so the tool works
# regardless of where the repo is checked out. Override with CLI args.
_SCRIPT_DIR      = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SRC_FILE = os.path.join(_SCRIPT_DIR, 'all_IMS_sts.xml')
DEFAULT_OUT_DIR  = os.path.join(_SCRIPT_DIR, 'XML_IM')
NS       = 'http://www.fdsn.org/xml/station/1'
TAG      = lambda t: f'{{{NS}}}{t}'   # helper: TAG('Station') → '{...}Station'

# ── helpers ──────────────────────────────────────────────────────────────────

def norm_desc(array_code, raw):
    """Return normalised sensor description."""
    if not raw:
        return array_code
    if raw.upper().startswith(array_code.upper()):
        return raw
    words = re.findall(r'[A-Za-z]+', raw)[:2]
    return f"{array_code} {' '.join(w.upper() for w in words)}" if words else array_code


def txt(el: ET.Element, tag: str) -> str:
    """Find child text, return '' if absent."""
    child = el.find(TAG(tag))
    return (child.text or '') if child is not None else ''


def sensitivity_parts(ch: ET.Element):
    """Return (value, freq, in_units, out_units) from <InstrumentSensitivity>."""
    resp = ch.find(TAG('Response'))
    if resp is None:
        return '', '', '', ''
    sens = resp.find(TAG('InstrumentSensitivity'))
    if sens is None:
        return '', '', '', ''
    val   = txt(sens, 'Value')
    freq  = txt(sens, 'Frequency')
    in_u_el  = sens.find(TAG('InputUnits'))
    out_u_el = sens.find(TAG('OutputUnits'))
    in_u  = txt(in_u_el,  'Name') if in_u_el  is not None else ''
    out_u = txt(out_u_el, 'Name') if out_u_el is not None else ''
    return val, freq, in_u, out_u


# ── per-station conversion ────────────────────────────────────────────────────

def write_station_xml(station_el: ET.Element, out_dir: str):
    array_code  = station_el.get('code')
    array_start = station_el.get('startDate')
    array_end   = station_el.get('endDate')

    # Collect all BDF channels
    all_ch = station_el.findall(TAG('Channel'))
    bdf    = [c for c in all_ch if c.get('code') == 'BDF']

    if not bdf:
        print(f'  SKIP {array_code}: no BDF channels')
        return

    # Group by sensor name, preserving first-seen order
    order   = []
    groups  = {}
    for ch in bdf:
        name = ch.get('name') or ''
        if name not in groups:
            order.append(name)
            groups[name] = []
        groups[name].append(ch)

    # Normalised description per sensor (use first Sensor/Description found)
    descs = {}
    for name, chs in groups.items():
        raw = None
        for ch in chs:
            s = ch.find(TAG('Sensor'))
            if s is not None:
                d = s.find(TAG('Description'))
                if d is not None and d.text:
                    raw = d.text
                    break
        descs[name] = norm_desc(array_code, raw)

    # ── build XML text ────────────────────────────────────────────────────────
    L = []
    def w(*args):
        L.append(''.join(str(a) for a in args))

    w('<?xml version="1.0" encoding="UTF-8"?>')
    w(f'<FDSNStationXML xmlns="{NS}" schemaVersion="1">')
    w('    <Source>CTBTO</Source>')
    w('    <Created>2026-04-14T00:00:00</Created>')
    w(f'    <Network code="IM" startDate="{array_start}" endDate="{array_end}">')
    w(f'        <TotalNumberStations>{len(order)}</TotalNumberStations>')

    for sensor_name in order:
        chs   = groups[sensor_name]
        desc  = descs[sensor_name]
        elem  = sensor_name[3:]          # e.g. "H1" from "I57H1"
        site  = f'{desc} {elem}'

        # Station date span
        st_start = min(c.get('startDate', '') for c in chs)
        st_end   = max(c.get('endDate',   '') for c in chs)

        # Station coordinates from first channel
        fc   = chs[0]
        lat  = txt(fc, 'Latitude')
        lon  = txt(fc, 'Longitude')
        elev = txt(fc, 'Elevation')

        w(f'        <Station code="{sensor_name}" startDate="{st_start}" endDate="{st_end}">')
        w(f'            <Latitude>{lat}</Latitude>')
        w(f'            <Longitude>{lon}</Longitude>')
        w(f'            <Elevation>{elev}</Elevation>')
        w( '            <Site>')
        w(f'                <Name>{site}</Name>')
        w( '            </Site>')

        for ch in chs:
            ch_start = ch.get('startDate', '')
            ch_end   = ch.get('endDate',   '')
            ch_lat   = txt(ch, 'Latitude')
            ch_lon   = txt(ch, 'Longitude')
            ch_elev  = txt(ch, 'Elevation')
            ch_depth = txt(ch, 'Depth')   or '0'
            ch_az    = txt(ch, 'Azimuth') or '0'
            ch_dip   = txt(ch, 'Dip')     or '0'
            ch_sr    = txt(ch, 'SampleRate')
            val, freq, in_u, out_u = sensitivity_parts(ch)

            w(f'            <Channel code="BDF" locationCode="" startDate="{ch_start}" endDate="{ch_end}">')
            w( '                <Sensor>')
            w(f'                    <Description>{desc}</Description>')
            w( '                </Sensor>')
            w(f'                <Latitude>{ch_lat}</Latitude>')
            w(f'                <Longitude>{ch_lon}</Longitude>')
            w(f'                <Elevation>{ch_elev}</Elevation>')
            w(f'                <Depth>{ch_depth}</Depth>')
            w(f'                <Azimuth>{ch_az}</Azimuth>')
            w(f'                <Dip>{ch_dip}</Dip>')
            w(f'                <SampleRate>{ch_sr}</SampleRate>')
            w( '                <Response>')
            w( '                    <InstrumentSensitivity>')
            if val:   w(f'                        <Value>{val}</Value>')
            if freq:  w(f'                        <Frequency>{freq}</Frequency>')
            if in_u:
                w( '                        <InputUnits>')
                w(f'                            <Name>{in_u}</Name>')
                w( '                        </InputUnits>')
            if out_u:
                w( '                        <OutputUnits>')
                w(f'                            <Name>{out_u}</Name>')
                w( '                        </OutputUnits>')
            w( '                    </InstrumentSensitivity>')
            w( '                </Response>')
            w( '            </Channel>')

        w('        </Station>')

    w('    </Network>')
    w('</FDSNStationXML>')

    out_path = os.path.join(out_dir, f'{array_code}.xml')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    print(f'  {array_code}.xml  ({len(order)} elements, {len(bdf)} BDF channels)')


# ── main ─────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Convert all_IMS_sts.xml → one XML file per IMS station.',
    )
    parser.add_argument(
        'src_file', nargs='?', default=DEFAULT_SRC_FILE,
        help='Source IMS station XML file '
             '(default: all_IMS_sts.xml next to this script).',
    )
    parser.add_argument(
        '-o', '--out', dest='out_dir', default=DEFAULT_OUT_DIR,
        help='Output directory for per-station XML files '
             '(default: XML_IM/ next to this script).',
    )
    return parser.parse_args(argv)


def main(src_file=DEFAULT_SRC_FILE, out_dir=DEFAULT_OUT_DIR):
    if not os.path.exists(src_file):
        raise SystemExit(
            f"Source file not found: {src_file}\n"
            f"Pass a path as the first argument, or place 'all_IMS_sts.xml' "
            f"next to this script."
        )
    os.makedirs(out_dir, exist_ok=True)

    print(f'Parsing {src_file} …')
    tree = ET.parse(src_file)
    root = tree.getroot()

    network  = root.find(TAG('Network'))
    stations = network.findall(TAG('Station'))
    print(f'Found {len(stations)} stations. Writing to {out_dir}/ …\n')

    for st in stations:
        write_station_xml(st, out_dir)

    written = len(os.listdir(out_dir))
    print(f'\nDone — {written} files in {out_dir}/')


if __name__ == '__main__':
    args = parse_args()
    main(args.src_file, args.out_dir)
