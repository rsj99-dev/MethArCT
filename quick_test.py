#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Minimal installation check for MethArCT.

Verifies that the package imports, the configuration resolves a database
directory and the expected reference/model files are present. This does not
require DIAMOND, CheckM2 or the optional prediction models.
"""

import os
import sys

try:
    import metharct
    from metharct import MethArCTAnalyzer, Config
except Exception as exc:  # pragma: no cover - import error is the failure
    print(f"[FAIL] Cannot import metharct: {exc}")
    sys.exit(1)


def main() -> int:
    print(f"MethArCT version: {metharct.__version__}")

    config = Config()
    db_dir = config.get('databases.base_dir')
    print(f"Database directory: {db_dir}")

    failures = []
    if not os.path.isdir(db_dir):
        failures.append(f"database directory not found: {db_dir}")

    expected_db_files = [
        os.path.join('diamond', 'SR.dmnd'),
        os.path.join('cultivation', 'path.aa'),
        os.path.join('kangshengsu', 'Methanobacterium.faa'),
    ]
    for rel_path in expected_db_files:
        if not os.path.exists(os.path.join(db_dir, rel_path)):
            failures.append(f"missing reference data: {rel_path}")

    cultivation_dirs = [
        'amino_acid', 'vitamin', 'nucleotide', 'atp'
    ]
    for pathway_type in cultivation_dirs:
        try:
            path = config.get_cultivation_pathway_dir(pathway_type)
        except ValueError as exc:
            failures.append(str(exc))
            continue
        if not os.path.isdir(path):
            failures.append(f"missing cultivation directory: {path}")

    # Public API smoke test (lazy import)
    assert MethArCTAnalyzer is not None

    if failures:
        for failure in failures:
            print(f"[FAIL] {failure}")
        return 1

    print("[OK] Package importable")
    print("[OK] Database directory resolved")
    print("[OK] Reference databases and cultivation directories found")
    return 0


if __name__ == '__main__':
    sys.exit(main())
