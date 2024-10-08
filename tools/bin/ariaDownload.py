#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga, Emre Havazli
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse
import getpass
import logging
import math
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from typing import List, Tuple

import ARIAtools.util.log
import asf_search
import shapely
from ARIAtools.util.shp import open_shp
from ARIAtools.util.url import url_versions
from tqdm import tqdm

logging.getLogger("asf_search").setLevel("ERROR")
LOGGER = logging.getLogger("ariaDownload.py")


def create_parser() -> argparse.ArgumentParser:
    """Create and return the argument parser."""
    parser = argparse.ArgumentParser(
        description='Download ARIA products using asf_search',
        epilog="Examples of use:\n"
        '\t ariaDownload.py --track 004 --output count\n'
        '\t ariaDownload.py --bbox "36.75 37.225 -76.655 -75.928"\n'
        '\t ariaDownload.py -t 004,077 --start 20190101 -o count\n'
        '\t ariaDownload.py --mission NISAR -o count',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-o",
        "--output",
        default="Download",
        type=str.title,
        choices=("Download", "Count", "Url"),
        help='Output type. Default="Download".',
    )
    parser.add_argument(
        "-t",
        "--track",
        default=None,
        type=str,
        help='Track to download; single number or comma separated.',
    )
    parser.add_argument(
        "-b",
        "--bbox",
        default=None,
        type=str,
        help=(
            'Lat/Lon Bounding SNWE, or GDAL-readable file containing '
            'POLYGON geometry.'
        ),
    )
    parser.add_argument(
        "-w",
        "--workdir",
        dest="wd",
        default="./products",
        type=str,
        help=(
            'Directory to deposit outputs. Default is "products" '
            'in local directory.'
        ),
    )
    parser.add_argument(
        "-s",
        "--start",
        default="20100101",
        type=str,
        help='Start date as YYYYMMDD.'
    )
    parser.add_argument(
        "-e",
        "--end",
        default="21000101",
        type=str,
        help='End date as YYYYMMDD.'
    )
    parser.add_argument(
        "-u",
        "--user",
        default=None,
        type=str,
        help='NASA Earthdata URS user login.'
    )
    parser.add_argument(
        "-p",
        "--pass",
        dest="passw",
        default=None,
        type=str,
        help='NASA Earthdata URS user password.',
    )
    parser.add_argument(
        "--mission",
        default="S1",
        type=str.upper,
        choices=("S1", "NISAR"),
        help='Sentinel-1 (S1) or NISAR. Default is S1',
    )
    parser.add_argument(
        "-l",
        "--daysless",
        dest="dayslt",
        default=math.inf,
        type=int,
        help='Take pairs with temporal baseline less than this value.',
    )
    parser.add_argument(
        "-m",
        "--daysmore",
        dest="daysgt",
        default=0,
        type=int,
        help='Take pairs with temporal baseline greater than this value.',
    )
    parser.add_argument(
        "-nt",
        "--num_threads",
        default="1",
        type=str,
        help='Number of threads for download. Default="1". Can specify "All".',
    )
    parser.add_argument(
        "-i",
        "--ifg",
        default=None,
        type=str,
        help=(
            'Retrieve one interferogram by its start/end date, '
            'as YYYYMMDD_YYYYMMDD.'
        ),
    )
    parser.add_argument(
        "-d",
        "--direction",
        dest="flightdir",
        default=None,
        type=str,
        help='Flight direction, options: ascending, a, descending, d',
    )
    parser.add_argument(
        "--version",
        default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods.',
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help='Print products to be downloaded to stdout',
    )
    parser.add_argument(
        "--log-level",
        default="warning",
        help="Logger log level"
    )

    return parser


def make_bbox(inp_bbox: str) -> shapely.geometry.Polygon:
    """Make a WKT from SNWE or a shapefile."""
    if inp_bbox is None:
        return None

    if os.path.exists(os.path.abspath(inp_bbox)):
        ring = open_shp(inp_bbox, 0, 0).exterior
        return shapely.geometry.Polygon(ring)

    try:
        S, N, W, E = map(float, inp_bbox.split())
        W, E = [x - 360 if x > 180 else x for x in (W, E)]
        N, S = [x - 90 if x > 90 else x for x in (N, S)]
        return shapely.geometry.Polygon([(W, N), (W, S), (E, S), (E, N)])
    except ValueError:
        raise ValueError(
            "Invalid bbox format. Use 'S N W E' or a valid shapefile path."
        )


def get_url_ifg(
    scenes: List[asf_search.ASFProduct],
) -> Tuple[List[str], List[str], bool]:
    """Get url, ifg of fetched ASF scene."""
    urls, ifgs = [], []
    for scene in scenes:
        s = scene.geojson()["properties"]
        urls.append(s["url"])

        if s["fileID"].startswith("NISAR_"):
            f = s["fileID"].split("_")
            pairname = f"{f[11][:8]}_{f[13][:8]}"
        else:
            f = s["fileID"].split("-")
            pairname = f[6]
        ifgs.append(pairname)

    is_nisar_file = any("/NISAR_" in url for url in urls)
    return urls, ifgs, is_nisar_file


def fmt_dst(args: argparse.Namespace) -> str:
    """Format the save name."""
    fn_track = f"track{args.track}".replace(",", "-") if args.track else ""

    if args.bbox:
        bbox = make_bbox(args.bbox)
        WSEN = [
            math.floor(bbox.bounds[i]) if i < 2 else math.ceil(bbox.bounds[i])
            for i in range(4)
        ]
        fn_bbox = f"_bbox{WSEN[0]}W{WSEN[1]}S{WSEN[2]}E{WSEN[3]}N"
    else:
        fn_bbox = ""

    base_name = f"{fn_track}{fn_bbox}_0.txt".lstrip("_")
    dst = os.path.join(args.wd, base_name)

    count = 1
    while os.path.exists(dst):
        base_name = f"{os.path.splitext(base_name)[0][:-1]}{count}.txt"
        dst = os.path.join(args.wd, base_name)
        count += 1

    return dst


class Downloader:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.args.output = self.args.output.title()
        self.args.wd = os.path.abspath(self.args.wd)
        os.makedirs(self.args.wd, exist_ok=True)
        LOGGER.setLevel(logging.DEBUG if self.args.verbose else logging.INFO)

    def __call__(self):
        scenes = self.query_asf()
        urls, ifgs, is_nisar_file = get_url_ifg(scenes)

        # Subset everything by version
        urls = url_versions(urls, self.args.version, self.args.wd)
        scenes = [scene for scene, url in zip(scenes, urls) if url in urls]
        ifgs = [ifg for ifg, url in zip(ifgs, urls) if url in urls]

        # Filter scenes based on date and elapsed time criteria
        scenes, urls, ifgs = self.filter_scenes(
            scenes,
            urls,
            ifgs,
            is_nisar_file
        )

        if self.args.output == "Count":
            LOGGER.info("Found -- %d -- products", len(scenes))
        elif self.args.output == "Url":
            self.write_urls(urls)
        elif self.args.output == "Download":
            self.download_scenes(scenes)

        if self.args.verbose:
            for scene in scenes:
                LOGGER.info(scene.geojson()["properties"]["sceneName"])

    def query_asf(self) -> List[asf_search.ASFProduct]:
        """Query ASF for scenes."""
        bbox = make_bbox(self.args.bbox)
        bbox_wkt = bbox.wkt if bbox else None

        flight_direction = None
        if self.args.flightdir:
            flight_direction = (
                "ascending"
                if self.args.flightdir.lower().startswith("a")
                else "descending"
            )

        tracks = (
            [int(track) for track in self.args.track.split(",")]
            if self.args.track
            else None
        )

        start = self.args.start - timedelta(days=1)
        end = self.args.end + timedelta(days=1)

        if self.args.mission.upper() == "S1":
            return asf_search.geo_search(
                collections=["C2859376221-ASF", "C1261881077-ASF"],
                dataset=asf_search.DATASET.ARIA_S1_GUNW,
                processingLevel=asf_search.constants.GUNW_STD,
                relativeOrbit=tracks,
                flightDirection=flight_direction,
                intersectsWith=bbox_wkt,
                start=start,
                end=end,
            )
        elif self.args.mission.upper() == "NISAR":
            session = asf_search.ASFSession()
            session.auth_with_token(getpass.getpass("EDL Token:"))
            LOGGER.info("Token accepted.")

            search_opts = asf_search.ASFSearchOptions(
                shortName="NISAR_L2_GUNW_BETA_V1",
                intersectsWith=bbox_wkt,
                start=start,
                end=end,
                session=session,
            )
            scenes = asf_search.search(opts=search_opts, maxResults=250)

            LOGGER.info("Found %d NISAR GUNW Betas.", len(scenes))
            return scenes

    def filter_scenes(self, scenes, urls, ifgs, is_nisar_file):
        filtered_scenes, filtered_urls, filtered_ifgs = [], [], []
        for scene, url, ifg in zip(scenes, urls, ifgs):
            eni, sti = self.parse_dates(ifg, is_nisar_file)

            if self.args.ifg:
                if self.match_single_ifg(sti, eni):
                    filtered_scenes.append(scene)
                    filtered_urls.append(url)
                    filtered_ifgs.append(ifg)
            elif self.match_date_criteria(sti, eni):
                filtered_scenes.append(scene)
                filtered_urls.append(url)
                filtered_ifgs.append(ifg)

        return filtered_scenes, filtered_urls, filtered_ifgs

    def parse_dates(self, ifg, is_nisar_file):
        if is_nisar_file:
            sti, eni = [datetime.strptime(d, "%Y%m%d") for d in ifg.split("_")]
        else:
            eni, sti = [datetime.strptime(d, "%Y%m%d") for d in ifg.split("_")]
        return eni, sti

    def match_single_ifg(self, sti, eni):
        dates = [
            datetime.strptime(i, "%Y%m%d").date()
            for i in self.args.ifg.split("_")
        ]
        st1, en1 = sorted(dates)
        return st1 == sti.date() and en1 == eni.date()

    def match_date_criteria(self, sti, eni):
        sten_chk = sti >= self.args.start and eni <= self.args.end
        elap = (eni - sti).days
        elap_chk = self.args.daysgt <= elap <= self.args.dayslt
        return sten_chk and elap_chk

    def write_urls(self, urls):
        dst = fmt_dst(self.args)
        with open(dst, "w") as fh:
            for url in urls:
                print(url, file=fh)
        LOGGER.info("Wrote -- %d -- product urls to: %s", len(urls), dst)

    def download_scenes(self, scenes):
        scenes = asf_search.ASFSearchResults(scenes)
        nt = int(self.args.num_threads)
        LOGGER.info("Downloading %d products...", len(scenes))

        session = asf_search.ASFSession()
        if self.args.user:
            session.auth_with_creds(self.args.user, self.args.passw)

        def download_file(url):
            local_filename = url.split("/")[-1]
            filepath = os.path.join(self.args.wd, local_filename)

            # Skip download if file already exists
            if os.path.exists(filepath):
                LOGGER.info("Product already in directory: %s", filepath)
                pbar.update(1)
                return filepath

            response = session.get(url, stream=True)
            response.raise_for_status()

            with open(filepath, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            return filepath

        # Create a progress bar
        pbar = tqdm(total=len(scenes), unit="file", desc="Downloading")
        try:
            urls = [scene.properties["url"] for scene in scenes]
            with ThreadPoolExecutor(max_workers=nt) as executor:
                future_to_url = {
                    executor.submit(download_file, url): url for url in urls
                }
                for future in as_completed(future_to_url):
                    url = future_to_url[future]
                    try:
                        filepath = future.result()
                        LOGGER.debug(f"Downloaded: {filepath}")
                        pbar.update(1)
                    except Exception as exc:
                        LOGGER.error(f"{url} generated an exception: {exc}")
        finally:
            pbar.close()

        LOGGER.info(
            "Download complete. Wrote -- %d -- products to: %s",
            len(scenes),
            self.args.wd,
        )


def main():
    parser = create_parser()
    args = parser.parse_args()

    log_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_level, int):
        raise ValueError(f"Invalid log level: {args.log_level}")
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    args.start = datetime.strptime(args.start, "%Y%m%d")
    args.end = datetime.strptime(args.end, "%Y%m%d")

    if not args.track and not args.bbox and args.mission.upper() != "NISAR":
        raise ValueError("Must specify either a bbox or track")

    Downloader(args)()


if __name__ == "__main__":
    main()
