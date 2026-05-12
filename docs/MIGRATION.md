# ARIA-tools v2 Migration Guide

This guide summarizes the practical migration from the legacy script-first
interface to the v2 `aria-tools <command>` interface.

The key principle is simple: the primary command surface changed, but the core
workflow arguments and outputs were intentionally kept familiar.

## What Changed

- Use `aria-tools <command>` as the primary interface.
- Keep legacy script entry points installed during the compatibility period.
- Use modern `pyproject.toml` packaging with optional extras.
- Use native Python process execution instead of GNU `parallel`.
- Run tests against the installed package rather than a source-path shim.

## Command Mapping

| Legacy script | v2 command |
| --- | --- |
| `ariaDownload.py` | `aria-tools download` |
| `ariaExtract.py` | `aria-tools extract` |
| `ariaTSsetup.py` | `aria-tools timeseries` |
| `ariaPlot.py` | `aria-tools plot` |
| `ariaOrderASF.py` | `aria-tools order` |
| `ariaMisclosure.py` | `aria-tools misclosure` |
| `ariaAOIassist.py` | `aria-tools aoi` |
| `ariaKml2box.py` | `aria-tools kml2box` |

For most workflows, the practical migration is just swapping the entry point.

## Installation Changes

v2 uses `pyproject.toml` as the authoritative pip install surface.

```bash
python -m pip install -e .
```

Optional functionality is installed with extras:

```bash
python -m pip install -e ".[aws]"
python -m pip install -e ".[order]"
python -m pip install -e ".[plot]"
python -m pip install -e ".[dev]"
```

Current extras:

| Extra | Purpose |
| --- | --- |
| `aws` | AWS/S3 direct-access support |
| `order` | HyP3 ordering dependencies |
| `plot` | Plotting dependencies |
| `dev` | Contributor tooling such as `pytest`, `ruff`, and `pre-commit` |

For repo work and CI, `environment.yml` remains intentionally broader than the
base pip install surface.

## Execution Model Changes

The biggest runtime change is the removal of GNU `parallel`.

- `extract` and `timeseries` now rely on native Python worker execution.
- The current supported backend is the Python process-pool path.
- Users still control concurrency with the familiar `--num_threads` option.

What did not change:

- core extract/timeseries arguments
- output structure expectations
- the legacy `ARIAtools` Python package import surface

## Virtual Access Notes

The `.txt` URL-list workflow remains supported:

```bash
aria-tools download --output url
aria-tools extract -f urls.txt -w workdir
```

The local metadata cache sidecar, `aria_meta_cache.json`, is still part of that
workflow. The first remote run is expected to be slower because metadata has to
be discovered before it can be reused.

Validation status:

- text-file URL inputs are still part of the supported interface
- metadata-cache helper behavior is covered by unit tests
- full remote parity still requires network access and appropriate credentials

## Testing and Contributor Workflow

The modernized test surface is install-first:

```bash
python -m pytest tests/unit -q
python -m pytest tests/integration -q
```

Opt-in markers are used for slower or environment-dependent checks:

- `offline`
- `slow`
- `network_required`
- `credentialed`

Examples:

```bash
python -m pytest tests/regression -q --run-slow --run-network --run-credentialed
python -m ruff check src/aria_tools tests/unit tests/integration
python -m ruff format --check src/aria_tools tests/unit tests/integration
```

## Backward Compatibility

The compatibility period is still active.

- Legacy scripts remain installed.
- `import ARIAtools` remains supported.
- Existing notebooks and internal scripts do not have to migrate all at once.

The recommended direction for new usage and updated documentation is:

```bash
aria-tools <command> [options]
```

## Known Boundaries

The v2 branch is modernized, but not fully re-architected.

- `product.py` and `extractProduct.py` remain large legacy cores.
- Secondary commands such as `download`, `plot`, and `order` still lean on
  compatibility wrappers.
- Full networked parity for credentialed workflows still needs environment
  access beyond the default local quick-test lane.

## Troubleshooting

### `aria-tools: command not found`

Install the package in the active environment:

```bash
python -m pip install -e .
```

### `ModuleNotFoundError: No module named 'ARIAtools'`

The repo is not installed in the current environment:

```bash
python -m pip install -e .
```

### First virtual-access run is slower than expected

That is normal while `aria_meta_cache.json` is being populated.

### A command help path works but a live workflow still fails

That usually means the interface is wired correctly but the workflow still
depends on data, credentials, or optional dependencies not present in the
current environment.

## Migration Checklist

- [ ] Install from `pyproject.toml` instead of relying on legacy setup flows
- [ ] Prefer `aria-tools <command>` in new docs, scripts, and examples
- [ ] Install only the extras your workflow needs
- [ ] Re-run your key extract/timeseries workflows with real data
- [ ] Keep legacy wrappers only where you still need compatibility
