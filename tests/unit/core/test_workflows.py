"""Unit tests for shared extract/timeseries workflow helpers."""

from __future__ import annotations

from types import SimpleNamespace

from aria_tools.core import workflows


def build_args() -> SimpleNamespace:
    return SimpleNamespace(
        imgfile="products/*.nc",
        bbox="1 2 3 4",
        projection="4326",
        workdir="/tmp/job",
        num_threads="2",
        version="all",
        nc_version="1b",
        verbose=False,
        tropo_models=None,
        layers="coherence",
        croptounion=False,
        demfile="Download",
        mask="mask.tif",
        minimumOverlap=0.0081,
        amp_thresh=None,
        outputFormat="VRT",
        multilooking=None,
        rankedResampling=False,
    )


def test_create_runlog_populates_standard_fields(monkeypatch) -> None:
    calls: list[tuple[str, object]] = []

    class FakeRunLog:
        def __init__(self, workdir: str) -> None:
            self.workdir = workdir

        def update(self, key: str, value: object) -> None:
            calls.append((key, value))

    args = build_args()
    monkeypatch.setattr(workflows.ARIAtools.util.runlog, "RunLog", FakeRunLog)
    monkeypatch.setattr(workflows.ARIAtools, "__version__", "2.0-test")

    runlog = workflows.create_runlog(args, routine_name="ariaExtract.py")

    assert runlog.workdir == "/tmp/job"
    assert calls == [
        ("aria_version", "2.0-test"),
        ("aria_routine", "ariaExtract.py"),
        ("args", args),
    ]


def test_build_standard_product_info_forwards_arguments(monkeypatch) -> None:
    captured: dict[str, object] = {}
    args = build_args()

    def fake_product(*positional, **kwargs):
        captured["positional"] = positional
        captured["kwargs"] = kwargs
        return "product-info"

    monkeypatch.setattr(workflows.ARIAtools.product, "Product", fake_product)

    result = workflows.build_standard_product_info(args, runlog="runlog")

    assert result == "product-info"
    assert captured["positional"] == ("products/*.nc",)
    assert captured["kwargs"] == {
        "bbox": "1 2 3 4",
        "projection": "4326",
        "workdir": "/tmp/job",
        "num_threads": "2",
        "url_version": "all",
        "nc_version": "1b",
        "verbose": False,
        "tropo_models": None,
        "layers": "coherence",
        "croptounion": False,
        "runlog": "runlog",
        "demfile": "Download",
        "mask": "mask.tif",
    }


def test_merge_product_bounding_boxes_forwards_expected_values(
    monkeypatch,
) -> None:
    captured: dict[str, object] = {}
    args = build_args()
    product_info = SimpleNamespace(
        products=[["meta"], ["layers"]],
        bbox_file="bbox.json",
    )

    def fake_merge(*positional, **kwargs):
        captured["positional"] = positional
        captured["kwargs"] = kwargs
        return ("merged",)

    monkeypatch.setattr(
        workflows.ARIAtools.extractProduct,
        "merged_productbbox",
        fake_merge,
    )

    result = workflows.merge_product_bounding_boxes(
        args,
        product_info=product_info,
        runlog="runlog",
    )

    assert result == ("merged",)
    assert captured["positional"] == (
        ["meta"],
        ["layers"],
        "/tmp/job/productBoundingBox",
        "bbox.json",
        False,
    )
    assert captured["kwargs"] == {
        "num_threads": "2",
        "minimumOverlap": 0.0081,
        "verbose": False,
        "runlog": "runlog",
    }


def test_collect_mask_source_products_preserves_order_for_s1() -> None:
    result = workflows.collect_mask_source_products(
        [
            {"amplitude": ["a", "b", "a"]},
            {"amplitude": ["b", "c"]},
        ],
        is_nisar_file=False,
    )

    assert result == ["a", "b", "c"]


def test_collect_mask_source_products_uses_coherence_for_nisar() -> None:
    result = workflows.collect_mask_source_products(
        [
            {"coherence": ["c1", "c2"], "amplitude": ["ignore"]},
            {"coherence": ["c2", "c3"]},
        ],
        is_nisar_file=True,
    )

    assert result == ["c1", "c2", "c3"]


def test_prepare_mask_returns_none_when_mask_is_not_requested() -> None:
    args = build_args()
    product_info = SimpleNamespace(mask=None, products=[None, []])

    assert (
        workflows.prepare_mask(
            args,
            product_info=product_info,
            is_nisar_file=False,
            bbox_file="bbox.json",
            prods_TOTbbox="bbox",
            proj="proj",
            arrres="arrres",
            runlog="runlog",
        )
        is None
    )


def test_prepare_mask_builds_expected_kwargs(monkeypatch) -> None:
    captured: dict[str, object] = {}
    args = build_args()
    product_info = SimpleNamespace(
        mask="mask.tif",
        products=[None, [{"amplitude": ["a", "b", "a"]}]],
    )

    monkeypatch.setattr(
        workflows.ARIAtools.util.mask,
        "prep_mask",
        lambda **kwargs: captured.setdefault("kwargs", kwargs) or "mask-out",
    )

    workflows.prepare_mask(
        args,
        product_info=product_info,
        is_nisar_file=False,
        bbox_file="bbox.json",
        prods_TOTbbox="bbox",
        proj="proj",
        arrres="arrres",
        runlog="runlog",
    )

    assert captured["kwargs"] == {
        "product_dict": ["a", "b"],
        "maskfilename": "mask.tif",
        "bbox_file": "bbox.json",
        "prods_TOTbbox": "bbox",
        "proj": "proj",
        "amp_thresh": None,
        "arrres": "arrres",
        "workdir": "/tmp/job",
        "outputFormat": "VRT",
        "num_threads": "2",
        "multilooking": None,
        "rankedResampling": False,
        "runlog": "runlog",
    }


def test_prepare_dem_builds_expected_kwargs(monkeypatch) -> None:
    captured: dict[str, object] = {}
    args = build_args()
    product_info = SimpleNamespace(demfile="Download")

    monkeypatch.setattr(
        workflows.ARIAtools.util.dem,
        "prep_dem",
        lambda **kwargs: captured.setdefault("kwargs", kwargs) or ("dem",),
    )

    workflows.prepare_dem(
        args,
        product_info=product_info,
        bbox_file="bbox.json",
        prods_TOTbbox="bbox",
        prods_TOTbbox_metadatalyr="metadata",
        proj="proj",
        arrres="arrres",
        runlog="runlog",
    )

    assert captured["kwargs"] == {
        "demfilename": "Download",
        "bbox_file": "bbox.json",
        "prods_TOTbbox": "bbox",
        "prods_TOTbbox_metadatalyr": "metadata",
        "proj": "proj",
        "arrres": "arrres",
        "workdir": "/tmp/job",
        "outputFormat": "VRT",
        "num_threads": "2",
        "multilooking": None,
        "rankedResampling": False,
        "runlog": "runlog",
    }
