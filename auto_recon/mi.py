# Copyright (C) 2025 ISIS Rutherford Appleton Laboratory UKRI
# SPDX - License - Identifier: GPL-3.0-or-later
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from mantidimaging import __version__
from mantidimaging.core.data import ImageStack
from mantidimaging.core.data.dataset import Dataset
from mantidimaging.core.io.loader import loader
from mantidimaging.core.io.loader.loader import create_loading_parameters_for_file_path, ImageParameters
from mantidimaging.core.operations.divide import DivideFilter
from mantidimaging.core.operations.loader import load_filter_packages
from mantidimaging.core.reconstruct import get_reconstructor_for
from mantidimaging.core.rotation import CorTiltDataModel
from mantidimaging.core.utility.data_containers import ReconstructionParameters, ScalarCoR, Degrees, FILE_TYPES
from mantidimaging.core.io.saver import image_save

FILTERS = {f.__name__: f for f in load_filter_packages()}


def version() -> str:
    return __version__


def load_dataset(file_path: Path) -> Dataset:
    parameters = create_loading_parameters_for_file_path(file_path)

    # Based on MainWindowModel.do_load_dataset()
    def load(im_param: ImageParameters) -> ImageStack:
        return loader.load_stack_from_image_params(im_param, None, dtype=parameters.dtype)

    sample = load(parameters.image_stacks[FILE_TYPES.SAMPLE])
    sample.set_geometry()
    ds = Dataset(sample=sample)
    sample._is_sinograms = parameters.sinograms
    sample.pixel_size = parameters.pixel_size

    for file_type in [
            FILE_TYPES.FLAT_BEFORE,
            FILE_TYPES.FLAT_AFTER,
            FILE_TYPES.DARK_BEFORE,
            FILE_TYPES.DARK_AFTER,
            FILE_TYPES.PROJ_180,
    ]:
        if im_param := parameters.image_stacks.get(file_type):
            image_stack = load(im_param)
            ds.set_stack(file_type, image_stack)

    return ds


def show_dataset(ds: Dataset):
    print(f"Dataset: {ds}")
    for stack in ds.all:
        print(f"  name: {stack.name} id: {stack.id}: {stack.data.shape}")


def run_operation(dataset: Dataset, op_name: str, params: dict[str, Any]):
    op_class = FILTERS[op_name]
    op_func = op_class.filter_func
    apply_to_dataset = True

    match op_name:
        case 'FlatFieldFilter':
            params = setup_flat_field(dataset, params)
            apply_to_dataset = False

    op_func(dataset.sample, **params)
    if apply_to_dataset:
        for stack in [dataset.flat_before, dataset.flat_after, dataset.dark_before, dataset.dark_after]:
            if stack:
                op_func(stack, **params)


def setup_flat_field(dataset: Dataset, params: dict[str, Any]) -> dict[str, Any]:
    params = dict(params)
    if dataset.flat_before:
        params['flat_before'] = dataset.flat_before
    if dataset.dark_before:
        params['dark_before'] = dataset.dark_before
    if dataset.flat_after:
        params['flat_after'] = dataset.flat_after
    if dataset.dark_after:
        params['dark_after'] = dataset.dark_after
    return params


default_settings = {'algorithm': 'FBP_CUDA', 'filter_name': 'ram-lak', 'cor': 1, 'tilt': 0, 'max_projection_angle': 360}


def run_recon(image_stack, settings=None):
    if settings is None:
        settings = {}
    settings = default_settings | settings

    do_clip = settings.pop('clip', False)

    reconstructor = get_reconstructor_for(settings['algorithm'])

    settings['cor'] = ScalarCoR(settings['cor'])
    settings['tilt'] = Degrees(settings['tilt'])

    params = ReconstructionParameters(**settings)

    cor_tilt = CorTiltDataModel()
    cor_tilt.set_precalculated(params.cor, params.tilt)

    cor_list = cor_tilt.get_all_cors_from_regression(image_stack.height)

    recon = reconstructor.full(image_stack, cor_list, params, progress=None)
    recon = DivideFilter.filter_func(recon, value=params.pixel_size, unit="micron", progress=None)

    if do_clip:
        np.clip(recon.data, a_min=0, a_max=None, out=recon.data)
    return recon


def save_stack(image_stack: ImageStack, out_dir: Path):
    image_save(image_stack, out_dir)


def show_stack(image_stack: ImageStack):
    from pyqtgraph.Qt import QtWidgets
    import pyqtgraph as pg
    pg.setConfigOptions(imageAxisOrder="row-major")
    pg.mkQApp("image stack")
    win = QtWidgets.QMainWindow()
    iv = pg.image(image_stack.data)
    win.setCentralWidget(iv)
    win.show()
    pg.exec()
