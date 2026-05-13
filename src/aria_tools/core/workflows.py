"""Shared workflow helpers for extract/timeseries modernization."""

from __future__ import annotations

import contextlib
import copy
import datetime
import glob
import json
import logging
import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import ARIAtools
import ARIAtools.extractProduct
import ARIAtools.product
import ARIAtools.util.dem
import ARIAtools.util.mask
import ARIAtools.util.misc
import ARIAtools.util.runlog
import ARIAtools.util.s3
import ARIAtools.util.vrt
import osgeo.gdal
from ARIAtools.constants import (
    ARIA_EXTERNAL_CORRECTIONS,
    ARIA_LAYERS,
    ARIA_TROPO_MODELS,
)

LOGGER = logging.getLogger(__name__)


def create_runlog(args: Any, *, routine_name: str) -> Any:
    """Create and populate the standard ARIA run log for a command."""

    runlog = ARIAtools.util.runlog.RunLog(args.workdir)
    runlog.update("aria_version", ARIAtools.__version__)
    runlog.update("aria_routine", routine_name)
    runlog.update("args", args)
    return runlog


def build_standard_product_info(args: Any, *, runlog: Any) -> Any:
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
) -> tuple[Any, ...]:
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
) -> Any | None:
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
) -> tuple[Any, Any, Any, Any]:
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


def normalize_runtime_num_threads(num_threads: int | str) -> int | str:
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
        logger.debug(
            "Using standard layers: %s", ARIAtools.constants.ARIA_STANDARD_LAYERS
        )
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
                "Extracting optional, user-specified layers %s for each "
                "interferogram pair",
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
    ref_dlist = generate_stack(
        context.product_info,
        "unwrappedPhase",
        "unwrapStack",
        workdir=args.workdir,
        is_nisar_file=context.is_nisar_file,
    )

    stack_layers = resolve_stack_layers_for_generation(
        layers,
        workdir=args.workdir,
        stack_defaults=stack_defaults,
        tropo_total=args.tropo_total,
    )
    for stack_layer, output_name, ref_tropokey in iter_stack_generation_requests(
        stack_layers,
        workdir=args.workdir,
        is_nisar_file=context.is_nisar_file,
        stack_outputs=stack_outputs,
    ):
        generate_stack(
            context.product_info,
            stack_layer,
            output_name,
            workdir=args.workdir,
            ref_dlist=ref_dlist,
            ref_tropokey=ref_tropokey,
            is_nisar_file=context.is_nisar_file,
        )


def extract_bperp_dict_ts(domain_name, aria_prod):
    """Extract mean perpendicular baseline values from products."""

    os.environ["GDAL_PAM_ENABLED"] = "NO"
    meta = {}
    for item in aria_prod:
        pair_name = item[-21:-4]
        stat = 0
        if domain_name == "unwrappedPhase":
            b_perp = item.split("/")
            b_perp[-2] = "bPerpendicular"
            b_perp = "/".join(b_perp)
            if os.path.exists(b_perp):
                data_set = None
                try:
                    data_set = osgeo.gdal.Open(b_perp, osgeo.gdal.GA_ReadOnly)
                    if data_set is not None:
                        band = data_set.GetRasterBand(1)
                        try:
                            stat = band.GetStatistics(True, True)[2]
                        except Exception:
                            stat = band.GetStatistics(False, True)[2]
                        band = None
                finally:
                    data_set = None

        meta[pair_name] = stat

    return meta


def extract_utc_time(aria_dates, aztime_list):
    """Extract UTC time metadata for stack records."""

    utc_dict = {}
    utc_time = None
    for index, pair_name in enumerate(aria_dates):
        if ([aztime_list[0]] * len(aztime_list) != aztime_list) or utc_time is None:
            mid_time = []
            for time_str in aztime_list[index]:
                mid_time.append(
                    datetime.datetime.strptime(time_str, "%Y-%m-%dT%H:%M:%S.%f")
                )
            min_mid_time = min(mid_time)
            max_mid_time = max(mid_time)
            time_delta = (max_mid_time - min_mid_time) / 2
            utc_time = (min_mid_time + time_delta).time()

        utc_dict[pair_name] = utc_time.strftime("%H:%M:%S.%f")
    return utc_dict


def generate_stack(
    aria_prod,
    stack_layer,
    output_file_name,
    workdir="./",
    ref_tropokey=None,
    ref_dlist=None,
    is_nisar_file=False,
):
    """Generate a time-series VRT stack for the requested layer."""

    os.environ["GDAL_PAM_ENABLED"] = "YES"
    stack_dir = os.path.join(workdir, "stack")
    if not os.path.exists(stack_dir):
        os.makedirs(stack_dir)

    domain_name = copy.deepcopy(stack_layer)
    data_type = "Int16" if domain_name == "connectedComponents" else "Float32"

    if domain_name in ARIA_TROPO_MODELS:
        stack_layer = f"{ref_tropokey}/{stack_layer}"
        stack_dir = os.path.join(stack_dir, ref_tropokey)
        if not os.path.exists(stack_dir):
            os.makedirs(stack_dir)

    if (
        domain_name in ARIA_EXTERNAL_CORRECTIONS or domain_name in ARIA_TROPO_MODELS
    ) and not is_nisar_file:
        stack_layer = f"{stack_layer}/dates"

    aria_dates = sorted([prod["pair_name"][0] for prod in aria_prod.products[0]])
    if not is_nisar_file and (
        domain_name in ARIA_EXTERNAL_CORRECTIONS or domain_name in ARIA_TROPO_MODELS
    ):
        aria_indiv_dates = []
        rejected_dates = []
        for aria_date in aria_dates:
            dates = aria_date.split("_")
            dt1_fname = os.path.join(workdir, stack_layer, dates[0] + ".vrt")
            if os.path.exists(dt1_fname):
                aria_indiv_dates += [dates[0]]
            else:
                rejected_dates += [dates[0]]

            dt2_fname = os.path.join(workdir, stack_layer, dates[1] + ".vrt")
            if os.path.exists(dt2_fname):
                aria_indiv_dates += [dates[1]]
            else:
                rejected_dates += [dates[1]]

        aria_dates = sorted(list(set(aria_indiv_dates)))
        rejected_dates = sorted(list(set(rejected_dates)))
        if rejected_dates:
            LOGGER.warning(
                "The following %d date(s) lack %s layers: %s",
                len(rejected_dates),
                domain_name,
                ", ".join(rejected_dates),
            )

    dlist = sorted(
        [
            os.path.join(workdir, stack_layer, aria_date + ".vrt")
            for aria_date in aria_dates
        ]
    )

    prog_bar = ARIAtools.util.misc.ProgressBar(
        maxValue=len(dlist),
        prefix=f"Exporting {output_file_name}: ",
    )

    b_perp = []
    new_dlist = [os.path.basename(item).split(".vrt")[0] for item in dlist]
    if is_nisar_file or (
        domain_name not in ARIA_EXTERNAL_CORRECTIONS
        and domain_name not in ARIA_TROPO_MODELS
    ):
        aztime_list = [
            product["azimuthZeroDopplerMidTime"] for product in aria_prod.products[0]
        ]
        b_perp_json_file = os.path.join(workdir, "bPerpendicular", "bperp.json")
        if os.path.isfile(b_perp_json_file):
            with open(b_perp_json_file) as ifp:
                b_perp = json.loads(ifp.read())
        else:
            b_perp = extract_bperp_dict_ts(domain_name, dlist)

        if ref_dlist and new_dlist != ref_dlist:
            subset_ind = []
            for index, stack_name in enumerate(new_dlist):
                if stack_name in ref_dlist:
                    subset_ind.append(index)
            new_dlist = [new_dlist[index] for index in subset_ind]
            dlist = [dlist[index] for index in subset_ind]
    else:
        aztime_list = len(aria_dates) * [
            aria_prod.products[0][0]["azimuthZeroDopplerMidTime"]
        ]

    utc_time = extract_utc_time(aria_dates, aztime_list)
    width, height, geo_trans, projection, no_data = ARIAtools.util.vrt.get_basic_attrs(
        dlist[0]
    )
    ymin, ymax, xmin, xmax = [0, height, 0, width]
    xsize = xmax - xmin
    ysize = ymax - ymin

    wvl = aria_prod.products[0][0]["wavelength"][0]
    start_range = aria_prod.products[0][0]["slantRangeStart"][0]
    end_range = aria_prod.products[0][0]["slantRangeEnd"][0]
    range_spacing = aria_prod.products[0][0]["slantRangeSpacing"][0]
    if is_nisar_file:
        orbit_direction = str.split(os.path.basename(aria_prod.files[0]), "_")[6]
        platform = "NISAR"
    else:
        orbit_direction = str.split(os.path.basename(aria_prod.files[0]), "-")[2]
        platform = "Sen"

    with open(os.path.join(stack_dir, output_file_name + ".vrt"), "w") as fid:
        fid.write(
            f'<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">\n'
            f"        <SRS>{projection}</SRS>\n"
            f"        <GeoTransform>"
            f"{geo_trans[0]},{geo_trans[1]},{geo_trans[2]},"
            f"{geo_trans[3]},{geo_trans[4]},{geo_trans[5]}"
            f"</GeoTransform>\n\n"
        )

        for index, data in enumerate(dlist, start=1):
            dates = data.split("/")[-1][:-4]
            prog_bar.update(index, suffix=dates)
            try:
                acq = utc_time[dates]
            except BaseException:
                continue

            if orbit_direction == "D":
                orb_dir = "DESCENDING"
            elif orbit_direction == "A":
                orb_dir = "ASCENDING"
            else:
                orb_dir = "UNKNOWN"

            path = os.path.relpath(os.path.abspath(data), start=stack_dir)
            outstr = f"""  <VRTRasterBand dataType="{data_type}" band="{index}">
        <NoDataValue>{no_data}</NoDataValue>
        <SimpleSource>
            <SourceFilename relativeToVRT="1">{path}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{width}" RasterYSize="{height}"
                DataType="{data_type}"/>
            <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}"
                ySize="{ysize}"/>
            <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
        </SimpleSource>
        <Metadata domain='{domain_name}'>
            <MDI key="Dates">{dates}</MDI>
            <MDI key="Wavelength (m)">{wvl}</MDI>
            <MDI key="UTCTime (HH:MM:SS.ss)">{acq}</MDI>
            <MDI key="startRange">{start_range}</MDI>
            <MDI key="endRange">{end_range}</MDI>
            <MDI key="slantRangeSpacing">{range_spacing}</MDI>
            <MDI key="orbitDirection">{orb_dir}</MDI>
            <MDI key="PLATFORM">{platform}</MDI>"""
            fid.write(outstr)
            if b_perp:
                fid.write(
                    f"""
            <MDI key="perpendicularBaseline">{b_perp[dates]}</MDI>"""
                )
            fid.write(
                """
        </Metadata>
    </VRTRasterBand>\n"""
            )
        fid.write("</VRTDataset>\n")
        prog_bar.close()

    return new_dlist


def resolve_stack_layers_for_generation(
    layers: list[str],
    *,
    workdir: str,
    stack_defaults: list[str],
    tropo_total: bool,
) -> list[str]:
    """Resolve the set of stack layers that should be generated."""

    resolved_layers = sorted(set(layers + stack_defaults))
    if "unwrappedPhase" in resolved_layers:
        resolved_layers.remove("unwrappedPhase")

    available_layers = []
    for layer in resolved_layers:
        if os.path.exists(os.path.join(workdir, layer)):
            available_layers.append(layer)

    if not tropo_total and "troposphereTotal" in available_layers:
        available_layers.remove("troposphereTotal")

    return available_layers


def iter_stack_generation_requests(
    layers: list[str],
    *,
    workdir: str,
    is_nisar_file: bool,
    stack_outputs: dict[str, str],
) -> list[tuple[str, str, str | None]]:
    """Expand stack layers into concrete stack-generation requests."""

    requests = []
    for layer in layers:
        if "tropo" in layer and not is_nisar_file:
            model_dirs = glob.glob(workdir + f"/{layer}/*", recursive=True)
            model_dirs = [os.path.basename(path) for path in model_dirs]
            for sublayer in model_dirs:
                if sublayer in stack_outputs:
                    requests.append((sublayer, stack_outputs[sublayer], layer))
            continue

        if layer not in stack_outputs:
            continue

        requests.append((layer, stack_outputs[layer], None))

    return requests
