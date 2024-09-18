# Copyright (C) 2024 ISIS Rutherford Appleton Laboratory UKRI
# SPDX - License - Identifier: GPL-3.0-or-later
from __future__ import annotations
import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy

from mantidimaging.versions import package_version
from mantidimaging.core.io.loader import loader
from mantidimaging.core.io.saver import image_save
from mantidimaging.core.io.filenames import FilenameGroup
from mantidimaging.core.data import ImageStack
from mantidimaging.core.utility.execution_timer import ExecutionTimer
from matplotlib import pyplot as plt


def main() -> None:
    print(f"Using: Mantid Imaging {package_version}")
    args = get_args()
    print(f"{args=}")

    shutters = get_shutters(args.sample)

    with ExecutionTimer(msg="load_stack Sample"):
        sample = load_stack(args.sample)
    with ExecutionTimer(msg="quick_overlap Sample"):
        sample_cor_mi = ImageStack(quick_overlap(sample.data, shutters))
    with ExecutionTimer(msg="save_stack Sample"):
        name_prefix, out_format = get_save_params(args, sample.filenames)
        save_path = args.sample / args.output
        save_stack(save_path, name_prefix, out_format, sample_cor_mi)

    sample_spec = get_spectrum(sample)
    sample_cor_mi_spec = get_spectrum(sample_cor_mi)

    if args.open:
        shutters_open = get_shutters(args.open)
        with ExecutionTimer(msg="load_stack Open Beam"):
            open_beam = load_stack(args.open)
        with ExecutionTimer(msg="quick_overlap Open Beam"):
            open_beam_cor_mi = ImageStack(quick_overlap(open_beam.data, shutters_open))
        with ExecutionTimer(msg="save_stack Open Beam"):
            name_prefix, out_format = get_save_params(args, open_beam.filenames)
            save_path = args.open / args.output
            save_stack(save_path, name_prefix, out_format, open_beam_cor_mi)

        open_beam_spec = get_spectrum(open_beam)
        open_beam_cor_mi_spec = get_spectrum(open_beam_cor_mi)

    if args.open:
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True)
    else:
        fig, ax1 = plt.subplots(nrows=1, sharex=True)

    ax1.set_title("Sample")
    ax1.plot(sample_spec, label="uncorrected")
    ax1.plot(sample_cor_mi_spec, label="corrected (MI)")
    #ax1.plot(sample_cor_spec, label="corrected")
    ax1.legend()

    if args.open:
        ax2.set_title("Open")
        ax2.plot(open_beam_spec, label="uncorrected")
        ax2.plot(open_beam_cor_mi_spec, label="corrected (MI)")

        ax3.set_title("Normalised")
        norm_spectrum = normalised_spectrum(sample_spec, open_beam_spec)
        ax3.plot(norm_spectrum, label="uncorrected")

        norm_spectrum_cor = normalised_spectrum(sample_cor_mi_spec, open_beam_cor_mi_spec)
        norm_spectrum_cor_min, norm_spectrum_cor_ymax = numpy.nanmin(norm_spectrum_cor), numpy.nanmax(norm_spectrum_cor)
        ax3.plot(norm_spectrum_cor, label="corrected (MI)")
        ax3.set_ylim(norm_spectrum_cor_min - 0.1, norm_spectrum_cor_ymax + 0.5)

        ax2.legend()
        ax3.legend()
    plt.show()


def get_save_params(args: ArgType, orig_names: list[str]) -> tuple[str, str]:
    orig_name = Path(orig_names[0]).name
    name_prefix = orig_name.rpartition("_")[0]
    out_format = orig_name.rpartition(".")[2]
    return name_prefix, out_format


def load_stack(path: Path) -> ImageStack:
    fname = sorted(path.glob("*.fits"))[0]
    filename_group = FilenameGroup.from_file(Path(fname))
    filename_group.find_all_files()
    image_stack = loader.load(filename_group)
    print(f"Loaded: {image_stack.data.shape}")
    return image_stack


def save_stack(path: Path, name_prefix: str, out_format: str, stack: ImageStack) -> None:
    print(f"{path=}")
    image_save(stack, output_dir=str(path), name_prefix=name_prefix, out_format=out_format, zfill_len=5)


def quick_overlap(raw_data: numpy.ndarray, shutters: list[ShutterInfo]) -> numpy.ndarray:
    print(f"{raw_data.max()=}")
    prob_occupied = numpy.zeros_like(raw_data, dtype=numpy.float32)

    for shutter in shutters:
        ss, se = shutter.start_index, shutter.end_index
        prob_occupied[ss + 1:se] = numpy.cumsum(raw_data[ss:se - 1], axis=0) / shutter.count

    #print(f"{prob_occupied.max()=}")
    #print(f"{prob_occupied[0]=}")
    corrected = raw_data / (1 - prob_occupied)
    return corrected


@dataclass
class ShutterInfo:
    number: int
    count: int
    start_time: float = 0
    end_time: float = 0
    start_index: int = 0
    end_index: int = 0


def get_shutters(data_dir: Path) -> list[ShutterInfo]:
    shuter_count_file = sorted(data_dir.glob("*ShutterCount.txt"))[0]
    shutter_times_file = sorted(data_dir.glob("*ShutterTimes.txt"))[0]
    spectra_file = sorted(data_dir.glob("*Spectra.txt"))[0]

    shuter_count = numpy.loadtxt(shuter_count_file, dtype=int)
    shutter_times = numpy.loadtxt(shutter_times_file)
    spectra = numpy.loadtxt(spectra_file)
    # print(shuter_count.shape, shutter_times.shape, spectra.shape)

    shutters = []
    prev_time = 0.0
    for number, count in shuter_count:
        if count == 0:
            break

        this_shutter = ShutterInfo(number, count)
        delay = shutter_times[number, 1]
        duration = shutter_times[number, 2]

        this_shutter.start_time = prev_time + delay
        this_shutter.end_time = this_shutter.start_time + duration
        prev_time = this_shutter.end_time
        this_shutter.start_index = int(numpy.searchsorted(spectra[:, 0], this_shutter.start_time))
        this_shutter.end_index = int(numpy.searchsorted(spectra[:, 0], this_shutter.end_time))
        shutters.append(this_shutter)
        print(this_shutter)

    return shutters


def get_spectrum(stack: ImageStack) -> numpy.ndarray:
    stack_spectrum = stack.data[:, 5:-5, 5:-5].mean(axis=(1, 2))
    return stack_spectrum


def normalised_spectrum(sample_spec: numpy.ndarray, open_spec: numpy.ndarray) -> numpy.ndarray:
    return numpy.divide(sample_spec, open_spec, out=numpy.zeros_like(sample_spec), where=open_spec != 0)


class ArgType(argparse.Namespace):
    sample: Path
    open: Path | None
    output: str


def get_args() -> ArgType:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", type=Path, help="Sample directory", required=True)
    parser.add_argument("--open", type=Path, help="Open beam directory")
    parser.add_argument("--output", type=str, help="Output directory", default="corrected_mi")
    return parser.parse_args(namespace=ArgType())


console_handler = logging.StreamHandler(sys.stdout)
perf_logger = logging.getLogger('perf')
perf_logger.setLevel(1)
perf_logger.addHandler(console_handler)

if __name__ == "__main__":
    main()
