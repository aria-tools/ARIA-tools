"""Shared workflow helpers for extract/timeseries modernization."""

from __future__ import annotations

import contextlib
import glob
import json
import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import ARIAtools
import ARIAtools.extractProduct
import ARIAtools.product
import ARIAtools.util.s3
import ARIAtools.util.dem
import ARIAtools.util.mask
import ARIAtools.util.runlog
import ARIAtools.util.vrt
import osgeo.gdal
from ARIAtools.constants import ARIA_LAYERS


def create_runlog(args: Any, *, routine_name: str):
    """Create and populate the standard ARIA run log for a command."""

    runlog = ARIAtools.util.runlog.RunLog(args.workdir)
    runlog.update("aria_version", ARIAtools.__version__)
    runlog.update("aria_routine", routine_name)
    runlog.update("args", args)
    return runlog


def build_standard_product_info(args: Any, *, runlog: Any):
    """Build the shared ARIA product model from command arguments."""

    return ARIAtools.product.Product(
        args.imgfile,
        bbox=args.bbox,
        projection=args.projection,
        workdir=args.workdir,
        num_threads=args.num_threads,
        url_version=args.version,
        nc_version=args.nc_version,
        verbose=args.verbose,
        tropo_models=args.tropo_models,
        layers=args.layers,
        croptounion=args.croptounion,
        runlog=runlog,
        demfile=args.demfile,
        mask=args.mask,
    )


def merge_product_bounding_boxes(
    args: Any,
    *,
    product_info: Any,
    runlog: Any,
):
    """Merge product bounding boxes for extract/timeseries setup."""

    return ARIAtools.extractProduct.merged_productbbox(
        product_info.products[0],
        product_info.products[1],
        os.path.join(args.workdir, "productBoundingBox"),
        product_info.bbox_file,
        args.croptounion,
        num_threads=args.num_threads,
        minimumOverlap=args.minimumOverlap,
        verbose=args.verbose,
        runlog=runlog,
    )


def collect_mask_source_products(
    product_layer_dicts: Iterable[dict[str, list[str]]],
    *,
    is_nisar_file: bool,
) -> list[str]:
    """Collect unique mask-source products while preserving first-seen order."""

    layer_name = "coherence" if is_nisar_file else "amplitude"
    seen: set[str] = set()
    products: list[str] = []

    for product_layers in product_layer_dicts:
        for item in product_layers.get(layer_name, []):
            if item in seen:
                continue
            seen.add(item)
            products.append(item)

    return products


def prepare_mask(
    args: Any,
    *,
    product_info: Any,
    is_nisar_file: bool,
    bbox_file: str,
    prods_TOTbbox: Any,
    proj: Any,
    arrres: Any,
    runlog: Any,
):
    """Prepare the optional workflow mask from shared command arguments."""

    if product_info.mask is None:
        return None

    return ARIAtools.util.mask.prep_mask(
        product_dict=collect_mask_source_products(
            product_info.products[1],
            is_nisar_file=is_nisar_file,
        ),
        maskfilename=product_info.mask,
        bbox_file=bbox_file,
        prods_TOTbbox=prods_TOTbbox,
        proj=proj,
        amp_thresh=args.amp_thresh,
        arrres=arrres,
        workdir=args.workdir,
        outputFormat=args.outputFormat,
        num_threads=args.num_threads,
        multilooking=args.multilooking,
        rankedResampling=args.rankedResampling,
        runlog=runlog,
    )


def prepare_dem(
    args: Any,
    *,
    product_info: Any,
    bbox_file: str,
    prods_TOTbbox: Any,
    prods_TOTbbox_metadatalyr: Any,
    proj: Any,
    arrres: Any,
    runlog: Any,
):
    """Prepare the shared workflow DEM inputs for extract/timeseries."""

    return ARIAtools.util.dem.prep_dem(
        demfilename=product_info.demfile,
        bbox_file=bbox_file,
        prods_TOTbbox=prods_TOTbbox,
        prods_TOTbbox_metadatalyr=prods_TOTbbox_metadatalyr,
        proj=proj,
        arrres=arrres,
        workdir=args.workdir,
        outputFormat=args.outputFormat,
        num_threads=args.num_threads,
        multilooking=args.multilooking,
        rankedResampling=args.rankedResampling,
        runlog=runlog,
    )


def normalize_runtime_num_threads(num_threads: Any) -> Any:
    """Normalize CLI-style worker settings for runtime helpers."""

    if isinstance(num_threads, str) and num_threads.lower() == "all":
        return "ALL_CPUS"
    return num_threads


@dataclass
class WorkflowContext:
    """Shared context for extract/timeseries workflow orchestration."""

    runlog: Any
    product_info: Any
    prods_TOTbbox: Any
    prods_TOTbbox_metadatalyr: Any
    arrres: Any
    proj: Any
    update_mode: Any
    is_nisar_file: bool
    demfile: Any
    demfile_expanded: Any
    lat: Any
    lon: Any
    maskfilename: Any


def prepare_workflow_context(
    args: Any,
    *,
    routine_name: str,
    build_dem: bool,
) -> WorkflowContext:
    """Prepare the shared setup context for extract/timeseries workflows."""

    runlog = create_runlog(args, routine_name=routine_name)
    product_info = build_standard_product_info(args, runlog=runlog)
    args.num_threads = normalize_runtime_num_threads(args.num_threads)
    (
        product_info.products[0],
        product_info.products[1],
        product_info.bbox_file,
        prods_TOTbbox,
        prods_TOTbbox_metadatalyr,
        arrres,
        proj,
        update_mode,
        is_nisar_file,
    ) = merge_product_bounding_boxes(
        args,
        product_info=product_info,
        runlog=runlog,
    )
    maskfilename = prepare_mask(
        args,
        product_info=product_info,
        is_nisar_file=is_nisar_file,
        bbox_file=product_info.bbox_file,
        prods_TOTbbox=prods_TOTbbox,
        proj=proj,
        arrres=arrres,
        runlog=runlog,
    )
    if build_dem:
        demfile, demfile_expanded, lat, lon = prepare_dem(
            args,
            product_info=product_info,
            bbox_file=product_info.bbox_file,
            prods_TOTbbox=prods_TOTbbox,
            prods_TOTbbox_metadatalyr=prods_TOTbbox_metadatalyr,
            proj=proj,
            arrres=arrres,
            runlog=runlog,
        )
    else:
        demfile, demfile_expanded, lat, lon = None, None, None, None

    return WorkflowContext(
        runlog=runlog,
        product_info=product_info,
        prods_TOTbbox=prods_TOTbbox,
        prods_TOTbbox_metadatalyr=prods_TOTbbox_metadatalyr,
        arrres=arrres,
        proj=proj,
        update_mode=update_mode,
        is_nisar_file=is_nisar_file,
        demfile=demfile,
        demfile_expanded=demfile_expanded,
        lat=lat,
        lon=lon,
        maskfilename=maskfilename,
    )


def run_extract_workflow(args: Any, *, logger: Any, valid_layers: Any) -> None:
    """Run the importable extract workflow."""

    args.tropo_total = False
    model_names = []
    runlog = create_runlog(args, routine_name="ariaExtract.py")
    product_info = build_standard_product_info(args, runlog=runlog)
    args.layers, args.tropo_total, model_names = ARIAtools.util.vrt.layerCheck(
        product_info.products[1],
        args.layers,
        args.nc_version,
        args.tropo_models,
        extract_or_ts="extract",
    )
    args.num_threads = normalize_runtime_num_threads(args.num_threads)
    (
        product_info.products[0],
        product_info.products[1],
        product_info.bbox_file,
        prods_TOTbbox,
        prods_TOTbbox_metadatalyr,
        arrres,
        proj,
        update_mode,
        is_nisar_file,
    ) = merge_product_bounding_boxes(
        args,
        product_info=product_info,
        runlog=runlog,
    )
    maskfilename = prepare_mask(
        args,
        product_info=product_info,
        is_nisar_file=is_nisar_file,
        bbox_file=product_info.bbox_file,
        prods_TOTbbox=prods_TOTbbox,
        proj=proj,
        arrres=arrres,
        runlog=runlog,
    )
    if product_info.demfile is not None:
        demfile, demfile_expanded, lat, lon = prepare_dem(
            args,
            product_info=product_info,
            bbox_file=product_info.bbox_file,
            prods_TOTbbox=prods_TOTbbox,
            prods_TOTbbox_metadatalyr=prods_TOTbbox_metadatalyr,
            proj=proj,
            arrres=arrres,
            runlog=runlog,
        )
    else:
        demfile, demfile_expanded, lat, lon = None, None, None, None

    logger.info(
        "Thread count specified for gdal multiprocessing = %s",
        args.num_threads,
    )
    if update_mode == "crop_only":
        args.layers = ARIAtools.extractProduct.track_existing_outputs(
            args.workdir,
            args.layers,
            valid_layers,
            [],
        )

    logger.info("Extracting products")
    ARIAtools.extractProduct.export_products(
        full_product_dict=product_info.products[1],
        bbox_file=product_info.bbox_file,
        prods_TOTbbox=prods_TOTbbox,
        proj=proj,
        layers=args.layers,
        iono_filter=args.iono_filter,
        is_nisar_file=is_nisar_file,
        arrres=arrres,
        rankedResampling=args.rankedResampling,
        demfile=demfile,
        demfile_expanded=demfile_expanded,
        lat=lat,
        lon=lon,
        maskfile=maskfilename,
        outDir=args.workdir,
        outputFormat=args.outputFormat,
        verbose=args.verbose,
        num_threads=args.num_threads,
        multilooking=args.multilooking,
        tropo_total=args.tropo_total,
        model_names=model_names,
        runlog=runlog,
        multiproc_method="processes",
    )


def run_timeseries_workflow(
    args: Any,
    *,
    logger: Any,
    generate_stack_func: Any,
    stack_defaults: list[str],
    stack_outputs: dict[str, str],
) -> None:
    """Run the importable timeseries workflow."""

    extract_bperp_layer = False
    with contextlib.suppress(TypeError):
        extract_bperp_layer = "bPerpendicular" in args.layers
        if extract_bperp_layer:
            layers = [layer.strip() for layer in args.layers.split(",")]
            layers.pop(layers.index("bPerpendicular"))
            args.layers = ",".join(layers)

    if args.layers.lower() == "standard":
        logger.debug("Using standard layers: %s", ARIAtools.constants.ARIA_STANDARD_LAYERS)
        args.layers = ",".join(ARIAtools.constants.ARIA_STANDARD_LAYERS)

    if "tropo" not in args.layers:
        args.tropo_models = None

    context = prepare_workflow_context(
        args,
        routine_name="ariaTSsetup.py",
        build_dem=True,
    )

    logger.info(
        "Thread count specified for gdal multiprocessing = %s",
        args.num_threads,
    )
    export_kwargs = {
        "proj": context.proj,
        "bbox_file": context.product_info.bbox_file,
        "prods_TOTbbox": context.prods_TOTbbox,
        "demfile": context.demfile,
        "demfile_expanded": context.demfile_expanded,
        "iono_filter": args.iono_filter,
        "is_nisar_file": context.is_nisar_file,
        "arrres": context.arrres,
        "lat": context.lat,
        "lon": context.lon,
        "maskfile": context.maskfilename,
        "outDir": args.workdir,
        "outputFormat": args.outputFormat,
        "verbose": args.verbose,
        "num_threads": args.num_threads,
        "multilooking": args.multilooking,
    }

    intf_layers = ARIAtools.constants.ARIA_STANDARD_INTF_LAYERS
    logger.info("Extracting %s for each interferogram pair", intf_layers)
    ref_arr_record = ARIAtools.extractProduct.export_products(
        context.product_info.products[1],
        tropo_total=False,
        layers=intf_layers,
        rankedResampling=args.rankedResampling,
        multiproc_method="processes",
        runlog=context.runlog,
        **export_kwargs,
    )

    osgeo.gdal.VSICurlClearCache()
    first_pair_dict = context.product_info.products[1][0]
    geom_layers = ARIAtools.constants.ARIA_STANDARD_GEOM_LAYERS
    logger.info(
        "Extracting single %s files valid over common interferometric grid",
        geom_layers,
    )
    prod_arr_record = ARIAtools.extractProduct.export_products(
        [first_pair_dict],
        tropo_total=False,
        layers=geom_layers,
        multiproc_method="single",
        runlog=context.runlog,
        **export_kwargs,
    )
    ARIAtools.util.vrt.dim_check(ref_arr_record, prod_arr_record)

    if extract_bperp_layer:
        logger.info(
            "Extracting perpendicular baseline grids for each interferogram pair"
        )
        prod_arr_record = ARIAtools.extractProduct.export_products(
            context.product_info.products[1],
            tropo_total=False,
            layers=["bPerpendicular"],
            multiproc_method="processes",
            runlog=context.runlog,
            **export_kwargs,
        )
        ARIAtools.util.vrt.dim_check(ref_arr_record, prod_arr_record)
    else:
        bperp_dict = ARIAtools.extractProduct.extract_bperp_dict(
            context.product_info.products[1],
            num_threads=args.num_threads,
        )
        bperp_outdir = os.path.join(args.workdir, "bPerpendicular")
        with contextlib.suppress(FileExistsError):
            os.mkdir(bperp_outdir)
        with open(os.path.join(bperp_outdir, "bperp.json"), "w") as ofp:
            json.dump(bperp_dict, ofp)

    layers, args.tropo_total, model_names = ARIAtools.util.vrt.layerCheck(
        context.product_info.products[1],
        args.layers,
        args.nc_version,
        args.tropo_models,
        extract_or_ts="tssetup",
    )
    if context.update_mode == "crop_only":
        ignore_names = [
            "unwrappedPhase",
            "connectedComponents",
            "incidenceAngle",
            "azimuthAngle",
            "coherence",
            "bPerpendicular",
        ]
        layers = ARIAtools.extractProduct.track_existing_outputs(
            args.workdir,
            layers,
            ARIA_LAYERS,
            ignore_names,
        )

    if layers or args.tropo_total is True:
        if layers:
            logger.info(
                "Extracting optional, user-specified layers %s for each interferogram pair",
                layers,
            )
        if args.tropo_total is True:
            logger.info(
                "Extracting, %s for each applicable interferogram pair",
                "troposphereTotal",
            )
        prod_arr_record = ARIAtools.extractProduct.export_products(
            context.product_info.products[1],
            tropo_total=args.tropo_total,
            model_names=model_names,
            layers=layers,
            multiproc_method="processes",
            runlog=context.runlog,
            **export_kwargs,
        )
        ARIAtools.util.vrt.dim_check(ref_arr_record, prod_arr_record)

    ARIAtools.util.s3.fixup_vrt_s3_paths(args.workdir)
    ref_dlist = generate_stack_func(
        context.product_info,
        "unwrappedPhase",
        "unwrapStack",
        workdir=args.workdir,
        is_nisar_file=context.is_nisar_file,
    )

    layers += stack_defaults
    layers.remove("unwrappedPhase")
    layers = sorted(set(layers))
    available_layers: list[str] = []
    for layer in layers:
        if os.path.exists(os.path.join(args.workdir, layer)):
            available_layers.append(layer)
    layers = available_layers
    if args.tropo_total is False and "troposphereTotal" in layers:
        layers.remove("troposphereTotal")

    stack_kwargs = {
        "workdir": args.workdir,
        "ref_dlist": ref_dlist,
        "is_nisar_file": context.is_nisar_file,
    }
    for layer in layers:
        if layer not in stack_outputs:
            logger.warning(
                "Selected layer %s not supported in tsSetupAvailable layers are: %s",
                layer,
                list(stack_outputs.keys()),
            )
            continue

        if "tropo" in layer and not context.is_nisar_file:
            model_dirs = glob.glob(args.workdir + f"/{layer}/*", recursive=True)
            model_dirs = [os.path.basename(i) for i in model_dirs]
            for sublyr in model_dirs:
                generate_stack_func(
                    context.product_info,
                    sublyr,
                    stack_outputs[sublyr],
                    ref_tropokey=layer,
                    **stack_kwargs,
                )
            continue

        generate_stack_func(
            context.product_info,
            layer,
            stack_outputs[layer],
            **stack_kwargs,
        )
