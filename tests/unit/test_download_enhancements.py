"""Offline regression tests for downloader validation and resume behavior."""

from __future__ import annotations

import hashlib
import importlib.util
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tools" / "bin" / "ariaDownload.py"


class DummyTqdm:
    def __init__(self, *args, **kwargs) -> None:
        self.args = args
        self.kwargs = kwargs

    def update(self, *_args, **_kwargs) -> None:
        return None

    def close(self) -> None:
        return None

    def set_postfix_str(self, *_args, **_kwargs) -> None:
        return None


class DummyResponse:
    def __init__(
        self,
        *,
        headers=None,
        status_code: int = 200,
        exc: Exception | None = None,
        body_chunks=None,
    ):
        self.headers = headers or {}
        self.status_code = status_code
        self.exc = exc
        self.body_chunks = list(body_chunks or [])

    def raise_for_status(self) -> None:
        if self.exc is not None:
            raise self.exc
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size=8192):
        del chunk_size
        yield from self.body_chunks

    def close(self) -> None:
        return None


class DummySession:
    def __init__(self, *, head=None, get=None):
        self._head = list(head or [])
        self._get = list(get or [])

    def head(self, *_args, **_kwargs):
        item = self._head.pop(0)
        if isinstance(item, Exception):
            raise item
        return item

    def get(self, *_args, **_kwargs):
        item = self._get.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


class FakeFuture:
    def __init__(self, result):
        self._result = result

    def result(self):
        return self._result


class FakeExecutor:
    def __init__(self, *, max_workers: int) -> None:
        self.max_workers = max_workers
        self.submitted: list[FakeFuture] = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None

    def submit(self, worker, item):
        future = FakeFuture(worker(item))
        self.submitted.append(future)
        return future


def _install_download_stubs(monkeypatch: pytest.MonkeyPatch) -> None:
    aria_pkg = types.ModuleType("ARIAtools")
    util_pkg = types.ModuleType("ARIAtools.util")
    log_mod = types.ModuleType("ARIAtools.util.log")
    s3_mod = types.ModuleType("ARIAtools.util.s3")
    shp_mod = types.ModuleType("ARIAtools.util.shp")
    url_mod = types.ModuleType("ARIAtools.util.url")
    asf_search_mod = types.ModuleType("asf_search")
    shapely_mod = types.ModuleType("shapely")
    tqdm_mod = types.ModuleType("tqdm")
    osgeo_pkg = types.ModuleType("osgeo")
    gdal_mod = types.ModuleType("osgeo.gdal")

    shp_mod.open_shp = lambda *_args, **_kwargs: None
    url_mod.url_versions = lambda urls, *_args, **_kwargs: urls
    tqdm_mod.tqdm = DummyTqdm
    gdal_mod.GA_ReadOnly = 0
    gdal_mod.Open = lambda *_args, **_kwargs: object()

    util_pkg.log = log_mod
    util_pkg.s3 = s3_mod
    util_pkg.shp = shp_mod
    util_pkg.url = url_mod
    aria_pkg.util = util_pkg

    osgeo_pkg.gdal = gdal_mod
    shapely_mod.geometry = types.SimpleNamespace(Polygon=lambda *_args, **_kwargs: None)
    shapely_mod.wkt = types.SimpleNamespace(loads=lambda *_args, **_kwargs: None)
    asf_search_mod.ASFSession = type("ASFSession", (), {})
    asf_search_mod.ASFSearchResults = list

    stub_modules = {
        "ARIAtools": aria_pkg,
        "ARIAtools.util": util_pkg,
        "ARIAtools.util.log": log_mod,
        "ARIAtools.util.s3": s3_mod,
        "ARIAtools.util.shp": shp_mod,
        "ARIAtools.util.url": url_mod,
        "asf_search": asf_search_mod,
        "shapely": shapely_mod,
        "tqdm": tqdm_mod,
        "osgeo": osgeo_pkg,
        "osgeo.gdal": gdal_mod,
    }

    for name, module in stub_modules.items():
        monkeypatch.setitem(sys.modules, name, module)


def load_download_module(monkeypatch: pytest.MonkeyPatch, *, module_name: str):
    _install_download_stubs(monkeypatch)
    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_remote_file_info_parser_understands_content_range(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_probe")
    response = DummyResponse(
        headers={"Content-Range": "bytes 0-0/10", "Accept-Ranges": "bytes"},
        status_code=206,
    )

    expected_size, supports_range = module._remote_file_info_from_response(response)

    assert expected_size == 10
    assert supports_range is True


def test_validate_resume_info_resumes_partial_file_when_remote_size_is_known(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_partial")
    filepath = tmp_path / "partial.nc"
    filepath.write_bytes(b"abc")
    monkeypatch.setattr(module, "_probe_remote_file_info", lambda *_args: (10, True))

    is_complete, resume_from, expected_size = module.validate_and_get_resume_info(
        str(filepath),
        "https://example.test/file.nc",
        DummySession(),
    )

    assert is_complete is False
    assert resume_from == 3
    assert expected_size == 10


def test_validate_resume_info_probes_remote_size_for_new_file(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_new_file")
    filepath = tmp_path / "fresh.nc"

    monkeypatch.setattr(module, "_probe_remote_file_info", lambda *_args: (25, True))

    is_complete, resume_from, expected_size = module.validate_and_get_resume_info(
        str(filepath),
        "https://example.test/file.nc",
        DummySession(),
    )

    assert is_complete is False
    assert resume_from == 0
    assert expected_size == 25


def test_validate_resume_info_accepts_checksum_verified_file_without_remote_size(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_checksum")
    payload = b"valid file payload"
    filepath = tmp_path / "complete.nc"
    filepath.write_bytes(payload)
    checksum = hashlib.md5(payload).hexdigest()

    monkeypatch.setattr(module, "_probe_remote_file_info", lambda *_args: (0, False))
    monkeypatch.setattr(module.osgeo.gdal, "Open", lambda *_args, **_kwargs: object())

    is_complete, resume_from, expected_size = module.validate_and_get_resume_info(
        str(filepath),
        "https://example.test/file.nc",
        DummySession(),
        expected_checksum=checksum,
        checksum_type="md5",
    )

    assert is_complete is True
    assert resume_from == filepath.stat().st_size
    assert expected_size == filepath.stat().st_size


def test_validate_resume_info_does_not_trust_gdal_without_remote_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_untrusted")
    filepath = tmp_path / "unknown.nc"
    filepath.write_bytes(b"maybe valid")

    monkeypatch.setattr(module, "_probe_remote_file_info", lambda *_args: (0, False))
    monkeypatch.setattr(module.osgeo.gdal, "Open", lambda *_args, **_kwargs: object())

    is_complete, resume_from, expected_size = module.validate_and_get_resume_info(
        str(filepath),
        "https://example.test/file.nc",
        DummySession(),
    )

    assert is_complete is False
    assert resume_from == 0
    assert expected_size == 0


def test_validate_resume_info_rejects_checksum_mismatch_even_if_sizes_match(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_mismatch")
    payload = b"not the expected payload"
    filepath = tmp_path / "mismatch.nc"
    filepath.write_bytes(payload)

    monkeypatch.setattr(
        module,
        "_probe_remote_file_info",
        lambda *_args: (filepath.stat().st_size, True),
    )
    monkeypatch.setattr(module.osgeo.gdal, "Open", lambda *_args, **_kwargs: object())

    is_complete, resume_from, expected_size = module.validate_and_get_resume_info(
        str(filepath),
        "https://example.test/file.nc",
        DummySession(),
        expected_checksum="0" * 32,
        checksum_type="md5",
    )

    assert is_complete is False
    assert resume_from == 0
    assert expected_size == filepath.stat().st_size


def test_validate_resume_info_uses_cached_local_validation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_cached")
    payload = b"validated earlier"
    filepath = tmp_path / "cached.nc"
    filepath.write_bytes(payload)

    monkeypatch.setattr(module, "_probe_remote_file_info", lambda *_args: (0, False))
    monkeypatch.setattr(module.osgeo.gdal, "Open", lambda *_args, **_kwargs: object())
    module._persist_validation_metadata(str(filepath), checksum_type="md5")

    is_complete, resume_from, expected_size = module.validate_and_get_resume_info(
        str(filepath),
        "https://example.test/file.nc",
        DummySession(),
    )

    assert is_complete is True
    assert resume_from == filepath.stat().st_size
    assert expected_size == filepath.stat().st_size


def test_download_file_resumable_uses_response_size_for_fresh_progress_bar(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_response_bar")

    created_bars = []

    class RecordingProgressBar:
        def __init__(self, filename, total_size, initial_size=0, position=1) -> None:
            self.filename = filename
            self.total_size = total_size
            self.initial_size = initial_size
            self.position = position
            self.updates: list[int] = []
            created_bars.append(self)

        def update(self, amount) -> None:
            self.updates.append(amount)

        def close(self) -> None:
            return None

        def set_postfix(self, _text) -> None:
            return None

    filepath = tmp_path / "fresh-download.nc"
    response = DummyResponse(
        headers={"Content-Length": "6"},
        body_chunks=[b"abc", b"def"],
    )
    session = DummySession(get=[response])

    monkeypatch.setattr(
        module,
        "validate_and_get_resume_info",
        lambda *_args, **_kwargs: (False, 0, 0),
    )
    monkeypatch.setattr(module, "DownloadProgressBar", RecordingProgressBar)
    monkeypatch.setattr(module.osgeo.gdal, "Open", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        module, "_persist_validation_metadata", lambda *_args, **_kwargs: None
    )

    success = module.download_file_resumable(
        str(filepath),
        "https://example.test/file.nc",
        session,
        show_progress=True,
    )

    assert success is True
    assert len(created_bars) == 1
    assert created_bars[0].total_size == 6
    assert created_bars[0].updates == [3, 3]


def test_download_scenes_counts_overall_progress_once_per_scene_and_reports_successes(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    module = load_download_module(monkeypatch, module_name="aria_download_progress")

    progress_bars = []

    class RecordingTqdm(DummyTqdm):
        def __init__(self, *args, **kwargs) -> None:
            super().__init__(*args, **kwargs)
            self.total = kwargs.get("total")
            self.position = kwargs.get("position", 0)
            self.updates: list[int] = []
            progress_bars.append(self)

        def update(self, amount=1, *_args, **_kwargs) -> None:
            self.updates.append(amount)

    monkeypatch.setattr(module.tqdm, "tqdm", RecordingTqdm)
    monkeypatch.setattr(
        module.ARIAtools.util.s3,
        "is_on_aws",
        lambda: False,
        raising=False,
    )
    monkeypatch.setattr(module, "_get_scene_checksum", lambda _scene: (None, None))
    monkeypatch.setattr(
        module,
        "validate_and_get_resume_info",
        lambda *_args, **_kwargs: (False, 0, 0),
    )
    monkeypatch.setattr(
        module,
        "download_file_resumable",
        lambda filepath, *_args, **_kwargs: not filepath.endswith("second.nc"),
    )

    class FakeASFSession:
        def auth_with_creds(self, *_args, **_kwargs) -> None:
            return None

    monkeypatch.setattr(module.asf_search, "ASFSession", FakeASFSession)
    monkeypatch.setattr(module.asf_search, "ASFSearchResults", list)

    executor = FakeExecutor(max_workers=2)
    monkeypatch.setattr(
        module.concurrent.futures,
        "ThreadPoolExecutor",
        lambda max_workers: executor,
    )
    monkeypatch.setattr(
        module.concurrent.futures,
        "as_completed",
        lambda futures: iter(futures),
    )

    log_messages: list[str] = []
    warn_messages: list[str] = []
    error_messages: list[str] = []
    monkeypatch.setattr(
        module.LOGGER,
        "info",
        lambda msg, *args: log_messages.append(msg % args if args else msg),
    )
    monkeypatch.setattr(
        module.LOGGER,
        "warning",
        lambda msg, *args: warn_messages.append(msg % args if args else msg),
    )
    monkeypatch.setattr(
        module.LOGGER,
        "error",
        lambda msg, *args: error_messages.append(msg % args if args else msg),
    )
    monkeypatch.setattr(module.LOGGER, "debug", lambda *_args, **_kwargs: None)

    args = types.SimpleNamespace(
        output="Download",
        wd=str(tmp_path),
        verbose=False,
        user=None,
        passw=None,
        num_threads="2",
    )
    downloader = module.Downloader(args)

    class FakeScene:
        def __init__(self, url: str) -> None:
            self.properties = {"url": url}

    scenes = [
        FakeScene("https://example.test/first.nc"),
        FakeScene("https://example.test/second.nc"),
    ]

    downloader.download_scenes(scenes)

    overall_bar = next(bar for bar in progress_bars if bar.position == 0)
    assert overall_bar.total == 2
    assert overall_bar.updates == [1, 1]
    assert any("Wrote -- 1/2 -- products" in message for message in log_messages)
    assert any("Failed to download" in message for message in error_messages)
    assert any("1 failed products" in message for message in warn_messages)
