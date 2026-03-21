#!/usr/bin/env python3
"""
Migrate old pyMDMix .config solvent files to new JSON format.

Reads INI-style .config files and converts them to the new JSON format
expected by pymdmix.core.solvent.
"""
from __future__ import annotations

import configparser
import json
import shutil
from pathlib import Path


# Probe type descriptions
PROBE_TYPE_DESC = {
    'Wat': 'Water oxygen',
    'Don': 'Hydrogen bond donor',
    'Acc': 'Hydrogen bond acceptor',
    'Hyd': 'Hydrophobic',
    'Pos': 'Positive charge',
    'Neg': 'Negative charge',
}


def parse_probe_spec(spec: str) -> tuple[str, list[str]]:
    """
    Parse probe specification like 'WAT@O' or 'ETA@O1,O2'.
    
    Returns (residue_name, [atom_names])
    """
    residue, atoms = spec.split('@')
    atom_list = [a.strip() for a in atoms.split(',')]
    return residue, atom_list


def convert_config_to_json(config_path: Path) -> dict:
    """Convert a .config file to JSON-compatible dictionary."""
    
    config = configparser.ConfigParser()
    config.read(config_path)
    
    # GENERAL section
    name = config.get('GENERAL', 'name').strip()
    info = config.get('GENERAL', 'info', fallback='').replace('%%', '%').strip()
    off_file = config.get('GENERAL', 'objectfile', fallback=None)
    box_unit = config.get('GENERAL', 'boxunit', fallback=None)
    water_model = config.get('GENERAL', 'watermodel', fallback='TIP3P')
    
    # PROBES section
    probes = []
    probe_types = {}
    
    # Get types if available
    if config.has_section('TYPES'):
        for probe_name, type_str in config.items('TYPES'):
            probe_types[probe_name.upper()] = type_str
    
    # Parse probes
    if config.has_section('PROBES'):
        for probe_name, spec in config.items('PROBES'):
            if not spec or spec.startswith('#'):
                continue
            
            probe_name = probe_name.upper()
            residue, atoms = parse_probe_spec(spec)
            
            # Build description from type
            probe_type = probe_types.get(probe_name, '')
            type_parts = [t.strip() for t in probe_type.split(',')]
            descriptions = [PROBE_TYPE_DESC.get(t, t) for t in type_parts if t]
            description = ', '.join(descriptions) if descriptions else ''
            
            probes.append({
                'name': probe_name,
                'residue': residue,
                'atoms': atoms,
                'description': description,
            })
    
    # Build residues list (simplified - could parse OFF for accurate counts)
    residue_names = set()
    for probe in probes:
        residue_names.add(probe['residue'])
    
    residues = [{'name': r, 'count': 0, 'description': ''} for r in sorted(residue_names)]
    
    return {
        'name': name,
        'full_name': name,
        'description': info,
        'off_file': off_file,
        'box_unit': box_unit,
        'water_model': water_model,
        'residues': residues,
        'probes': probes,
    }


def migrate_solvents(
    old_dir: Path,
    new_dir: Path,
    copy_off_files: bool = True,
) -> list[str]:
    """
    Migrate all .config files to JSON format.
    
    Parameters
    ----------
    old_dir : Path
        Directory containing old .config and .off files
    new_dir : Path
        Target directory for JSON files and .off copies
    copy_off_files : bool
        Whether to copy .off files to new location
        
    Returns
    -------
    list[str]
        Names of migrated solvents
    """
    new_dir.mkdir(parents=True, exist_ok=True)
    migrated = []
    
    for config_file in old_dir.glob('*.config'):
        print(f"Processing: {config_file.name}")
        
        try:
            data = convert_config_to_json(config_file)
            name = data['name']
            
            # Write JSON
            json_path = new_dir / f"{name.lower()}.json"
            with open(json_path, 'w') as f:
                json.dump(data, f, indent=2)
            print(f"  → Created: {json_path.name}")
            
            # Copy OFF file if requested
            if copy_off_files and data.get('off_file'):
                off_src = old_dir / data['off_file']
                if off_src.exists():
                    off_dst = new_dir / data['off_file']
                    if not off_dst.exists():
                        shutil.copy2(off_src, off_dst)
                        print(f"  → Copied: {data['off_file']}")
                else:
                    print(f"  ⚠ OFF file not found: {off_src}")
            
            migrated.append(name)
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    return migrated


def main():
    """Run migration."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Migrate pyMDMix solvents to JSON')
    parser.add_argument('--old-dir', type=Path, 
                        default=Path(__file__).parent.parent / 'pyMDMix/data/solventlib',
                        help='Source directory with .config files')
    parser.add_argument('--new-dir', type=Path,
                        default=Path(__file__).parent.parent / 'pymdmix/data/solvents',
                        help='Target directory for JSON files')
    parser.add_argument('--no-copy-off', action='store_true',
                        help='Do not copy .off files')
    
    args = parser.parse_args()
    
    print(f"Migrating solvents from: {args.old_dir}")
    print(f"                    to: {args.new_dir}")
    print()
    
    migrated = migrate_solvents(
        args.old_dir,
        args.new_dir,
        copy_off_files=not args.no_copy_off,
    )
    
    print()
    print(f"✓ Migrated {len(migrated)} solvents: {', '.join(migrated)}")


if __name__ == '__main__':
    main()
