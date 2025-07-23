# Copyright (C) 2025 ISIS Rutherford Appleton Laboratory UKRI
# SPDX - License - Identifier: GPL-3.0-or-later
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import mi

PROCESS_FILE = Path("process.json")
OUT_DIR = Path("/tmp/test_output")


def main():
    print(f"Mantid Imaging {mi.version()}")

    process = load_process(PROCESS_FILE)

    dataset = mi.load_dataset(Path(process['dataset_dir']))
    mi.show_dataset(dataset)

    for operation in process['operation_history']:
        print(operation['name'], operation['kwargs'])
        if "Recon" not in operation['name']:
            mi.run_operation(dataset, operation['name'], operation['kwargs'])
        else:
            settings = operation['kwargs'] | {'pixel_size': process['pixel_size']}
            recon = mi.run_recon(dataset.sample, settings)
            dataset.add_recon(recon)

    # mi.show_stack(dataset.recons[0])
    mi.save_stack(dataset.recons[0], OUT_DIR)


def load_process(file_name: Path) -> dict[str, Any]:
    return json.load(file_name.open())


if __name__ == '__main__':
    main()
