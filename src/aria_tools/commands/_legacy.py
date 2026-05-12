"""Adapters from modern subcommands to legacy script entry points."""

from __future__ import annotations

import importlib.util
import runpy
import sys
import textwrap
from collections.abc import Sequence
from types import ModuleType

from aria_tools.cli.common import legacy_script_path

LEGACY_HELP = {
    "ariaDownload.py": (
        "download",
        "Download or list Sentinel-1/NISAR GUNW products from ASF.",
        "-o/--output, -t/--track, -b/--bbox, -w/--workdir, -s/--start, "
        "-e/--end, -u/--user, -p/--pass, --mission, -l/--daysless, "
        "-m/--daysmore, -nt/--num_threads, -i/--ifg, -d/--direction, "
        "--version, -v/--verbose, --log-level",
    ),
    "ariaExtract.py": (
        "extract",
        "Extract product layers from ARIA GUNW products.",
        "-f/--file, -w/--workdir, -l/--layers, -tm/--tropo_models, "
        "-d/--demfile, -p/--projection, -b/--bbox, -m/--mask, "
        "-at/--amp_thresh, -nt/--num_threads, -of/--outputFormat, "
        "-croptounion/--croptounion, -ml/--multilooking, "
        "-rr/--rankedResampling, -mo/--minimumOverlap, -if/--iono_filter, "
        "--version, --nc_version, -verbose/--verbose, --log-level",
    ),
    "ariaTSsetup.py": (
        "timeseries",
        "Prepare ARIA products for time-series processing.",
        "-f/--file, -w/--workdir, -l/--layers, -tm/--tropo_models, "
        "-d/--demfile, -p/--projection, -b/--bbox, -m/--mask, "
        "-at/--amp_thresh, -nt/--num_threads, -of/--outputFormat, "
        "-croptounion/--croptounion, -ml/--multilooking, "
        "-rr/--rankedResampling, -mo/--minimumOverlap, -if/--iono_filter, "
        "--version, --nc_version, -verbose/--verbose, --log-level",
    ),
    "ariaPlot.py": (
        "plot",
        "Generate quality-control and baseline plots.",
        "-f/--file, -w/--workdir, -b/--bbox, -m/--mask, -at/--amp_thresh, "
        "-nt/--num_threads, -of/--outputFormat, -croptounion/--croptounion, "
        "-plottracks/--plottracks, -plotbperp/--plotbperp, "
        "-plotbperpcoh/--plotbperpcoh, -plotcoh/--plotcoh, "
        "-makeavgoh/--makeavgoh, -plotall/--plotall, -mo/--minimumOverlap, "
        "--figwidth, --version, --nc_version, -v/--verbose, --log-level",
    ),
    "ariaOrderASF.py": (
        "order",
        "Discover frames, generate pairs, order HyP3 jobs, or check status.",
        "--getframes, --getpairs, --orderpairs, --statusjobs, -b/--bbox, "
        "-t/--track, -d/--direction, --frame, --network, --num-neighbors, "
        "--seasonal-window, -s/--start, -e/--end, --pairs-file, --job-name, "
        "--dry-run, --limit, --status-name, --job-ids, -w/--workdir, "
        "-v/--verbose, --log-level",
    ),
    "ariaMisclosure.py": (
        "misclosure",
        "Compute and analyze phase triplet misclosure.",
        "-f/--file, -w/--workdir, --startdate, --enddate, --exclude-pairs, "
        "--plot-pairs, --mintime, --maxtime, --print-triplets, "
        "--plot-triplets, -refX, -refY, -refLon, -refLat, --queryX, "
        "--queryY, --queryLon, --queryLat, -v/--verbose, --pctmin, "
        "--pctmax, --plot-time-intervals, --log-level",
    ),
    "ariaAOIassist.py": (
        "aoi",
        "Use ASF CSV metadata to assist AOI creation.",
        "-f/--file, -w/--workdir, -t/--tracks, -l/--lat_bounds, "
        "-s/--start_date, -e/--end_date, -x/--exclude_dates, --plot_raw, "
        "--flag_partial_coverage, --remove_incomplete_dates, "
        "--approximate_AOI, -v/--verbose",
    ),
    "ariaKml2box.py": (
        "kml2box",
        "Convert KML/KMZ polygons to GeoJSON bounding boxes.",
        "-w/--workdir, -f/--file, -o/--outfile, --log-level",
    ),
}

SCRIPT_REPLACEMENTS = {
    "ariaDownload.py": "aria-tools download",
    "ariaExtract.py": "aria-tools extract",
    "ariaTSsetup.py": "aria-tools timeseries",
    "ariaPlot.py": "aria-tools plot",
    "ariaOrderASF.py": "aria-tools order",
    "ariaMisclosure.py": "aria-tools misclosure",
    "ariaAOIassist.py": "aria-tools aoi",
    "ariaKml2box.py": "aria-tools kml2box",
    "aria*.py": "aria-tools commands",
}


def _wants_help(argv: Sequence[str] | None) -> bool:
    return argv is not None and any(arg in {"-h", "--help"} for arg in argv)


def _print_static_help(script_name: str) -> None:
    command, summary, flags = LEGACY_HELP[script_name]
    print(f"usage: aria-tools {command} [options]\n")
    print(summary)
    print("\nOptions:")
    print(
        textwrap.fill(
            flags,
            width=78,
            subsequent_indent="  ",
            break_on_hyphens=False,
        )
    )
    print(
        "\nInstall the command dependencies to see full option descriptions "
        "and examples."
    )


def _is_missing_optional_dependency(exc: ModuleNotFoundError) -> bool:
    return exc.name is not None and not exc.name.startswith("aria_tools")


def _load_legacy_module(script_name: str) -> ModuleType:
    script_path = legacy_script_path(script_name)
    module_name = f"aria_tools._delegated_{script_path.stem}"
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import legacy command script: {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _legacy_parser(module: ModuleType):
    if hasattr(module, "createParser"):
        return module.createParser()
    if hasattr(module, "create_parser"):
        return module.create_parser()
    return None


def _modernize_help(help_text: str) -> str:
    for old, new in SCRIPT_REPLACEMENTS.items():
        help_text = help_text.replace(old, new)
    return help_text


def _print_parser_help(script_name: str) -> bool:
    command = LEGACY_HELP.get(script_name, (script_name,))[0]
    module = _load_legacy_module(script_name)
    parser = _legacy_parser(module)
    if parser is None:
        return False

    parser.prog = f"aria-tools {command}"
    print(_modernize_help(parser.format_help()), end="")
    return True


def run_legacy_script(script_name: str, argv: Sequence[str] | None = None) -> None:
    """Execute a legacy ``tools/bin`` script as ``__main__``.

    This is a temporary compatibility bridge. It keeps the command behavior
    identical while domain logic is migrated into importable modules.
    """

    if _wants_help(argv) and script_name in LEGACY_HELP:
        try:
            if _print_parser_help(script_name):
                return
        except ModuleNotFoundError as exc:
            if _is_missing_optional_dependency(exc):
                _print_static_help(script_name)
                return
            raise

    script_path = legacy_script_path(script_name)
    command = LEGACY_HELP.get(script_name, (script_name,))[0]
    old_argv = sys.argv[:]
    sys.argv = [
        f"aria-tools {command}",
        *(list(argv) if argv is not None else sys.argv[1:]),
    ]
    try:
        runpy.run_path(str(script_path), run_name="__main__")
    except ModuleNotFoundError as exc:
        if (
            _wants_help(argv)
            and script_name in LEGACY_HELP
            and (_is_missing_optional_dependency(exc))
        ):
            _print_static_help(script_name)
            return
        raise
    finally:
        sys.argv = old_argv
