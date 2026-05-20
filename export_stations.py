#!/usr/bin/env python3
"""
export_stations.py — Export individual IMS stations from all_IMS_sts.xml to separate XML files.

Copyright (c) 2024-2025 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.

Each file follows the same structure as I51GB.xml/I57US.xml:
- Network code is derived from the station code prefix (e.g., I51GB -> IM, I57US -> IM)
- Each unique channel name becomes a separate Station with its own Site and coordinates
"""

import xml.etree.ElementTree as ET
import os
import copy

# Portable defaults — resolved relative to this script so the tool works
# regardless of the current working directory. See also convert_ims.py.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE  = os.path.join(_SCRIPT_DIR, "all_IMS_sts.xml")
OUTPUT_DIR  = os.path.join(_SCRIPT_DIR, "XML_IM")

NS = {
    'fdsn': 'http://www.fdsn.org/xml/station/1',
    'dtk': 'http://www.fdsn.org/xml/dtk/1'
}

def strip_ns(elem):
    """Recursively strip XML namespace prefixes from an element's tag and attributes."""
    if '}' in elem.tag:
        elem.tag = elem.tag.split('}', 1)[-1]
    for k in list(elem.attrib):
        if '}' in k:
            elem.set(k.split('}', 1)[-1], elem.attrib.pop(k))
    for child in elem:
        strip_ns(child)


def indent_xml(elem, level=0):
    """Add indentation to XML elements."""
    indent = "\n" + level * "    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = indent + "    "
        if not elem.tail or not elem.tail.strip():
            elem.tail = indent
        for i, child in enumerate(elem):
            indent_xml(child, level + 1)
            if not child.tail or not child.tail.strip():
                child.tail = indent + "    " if i < len(elem) - 1 else indent
        if not elem[-1].tail or not elem[-1].tail.strip():
            elem[-1].tail = indent
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = indent


def extract_network_code(station_code):
    """Extract network code from station code (e.g., I51GB -> IM, I57US -> IM)."""
    # Pattern: IXXYY -> IM where XX is the number
    if station_code.startswith('I') and len(station_code) >= 3:
        number_part = station_code[1:3]
        return f"I{number_part}"
    return "IM"


def main():
    # Fail early with a clear message if the source file is missing
    if not os.path.exists(INPUT_FILE):
        raise SystemExit(
            f"Source file not found: {INPUT_FILE}\n"
            f"Place 'all_IMS_sts.xml' next to this script, or edit INPUT_FILE."
        )

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Parse the big XML file
    tree = ET.parse(INPUT_FILE)
    root = tree.getroot()

    # Register namespaces
    ET.register_namespace('', 'http://www.fdsn.org/xml/station/1')
    ET.register_namespace('dtk', 'http://www.fdsn.org/xml/dtk/1')

    # Find all Network elements
    networks = root.findall('fdsn:Network', NS)

    total_stations = 0
    total_files = 0

    for network in networks:
        stations = network.findall('fdsn:Station', NS)

        for station in stations:
            station_code = station.get('code')
            network_code_out = "IM"  # IMS network code

            # Group channels by their 'name' attribute (each becomes a station)
            channels = station.findall('fdsn:Channel', NS)
            channel_groups = {}
            for ch in channels:
                ch_name = ch.get('name')
                if ch_name not in channel_groups:
                    channel_groups[ch_name] = []
                channel_groups[ch_name].append(ch)

            # Create one XML file per station code (containing all channel-group stations)
            new_root = ET.Element('FDSNStationXML')
            new_root.set('xmlns', 'http://www.fdsn.org/xml/station/1')
            new_root.set('xmlns:dtk', 'http://www.fdsn.org/xml/dtk/1')
            new_root.set('schemaVersion', '1')

            # Add Source and Created elements (required by ObsPy)
            ET.SubElement(new_root, 'Source').text = "CTBTO"
            ET.SubElement(new_root, 'Created').text = "2024-01-01T00:00:00.000000Z"

            new_network = ET.SubElement(new_root, 'Network')
            new_network.set('code', network_code_out)

            total_stations += len(channel_groups)

            # Create a Station for each channel group
            for ch_name, ch_list in sorted(channel_groups.items()):
                # Use first channel for coordinates and station info
                first_ch = ch_list[0]
                ch_lat = first_ch.find('fdsn:Latitude', NS)
                ch_lon = first_ch.find('fdsn:Longitude', NS)
                ch_elev = first_ch.find('fdsn:Elevation', NS)
                ch_start = first_ch.get('startDate', '')
                ch_end = first_ch.get('endDate', '')

                # Build station element
                new_station = ET.SubElement(new_network, 'Station')
                new_station.set('code', ch_name)
                new_station.set('startDate', ch_start)
                new_station.set('endDate', ch_end)

                # Add Site element
                site = ET.SubElement(new_station, 'Site')
                name_elem = ET.SubElement(site, 'Name')
                name_elem.text = f"{station_code} {ch_name}"

                # Copy Latitude, Longitude, Elevation
                if ch_lat is not None:
                    lat_elem = ET.SubElement(new_station, 'Latitude')
                    lat_elem.text = ch_lat.text
                if ch_lon is not None:
                    lon_elem = ET.SubElement(new_station, 'Longitude')
                    lon_elem.text = ch_lon.text
                if ch_elev is not None:
                    elev_elem = ET.SubElement(new_station, 'Elevation')
                    elev_elem.text = ch_elev.text

                # Copy all channels (strip namespace from tags and attributes)
                for ch in ch_list:
                    new_ch = copy.deepcopy(ch)
                    strip_ns(new_ch)
                    new_station.append(new_ch)

            # Write to file
            output_file = os.path.join(OUTPUT_DIR, f"{station_code}.xml")

            # Add indentation
            indent_xml(new_root)
            if new_root.tail:
                new_root.tail = "\n"

            new_tree = ET.ElementTree(new_root)
            new_tree.write(output_file, encoding="utf-8", xml_declaration=True)

            total_files += 1
            print(f"Exported: {station_code} -> {len(channel_groups)} stations ({network_code_out})")

    print(f"\nTotal files exported: {total_files}")
    print(f"Total stations across all files: {total_stations}")
    print(f"Output directory: {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()
