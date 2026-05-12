# ARIA-tools v2 Architecture

This note describes the current architecture of the v2 branch as it exists
during the modernization transition. It is intentionally practical rather than
aspirational.

## Architecture Summary

The repo now has a modern command surface, but it is still a hybrid codebase.

- `src/aria_tools/` contains the modern CLI, routing, workflow seams, and
  execution helpers.
- `tools/ARIAtools/` still contains large legacy domain modules that the modern
  layer delegates into.
- `tools/bin/` legacy scripts remain installed for compatibility while the new
  command surface becomes the default.

The target is a testable, mergeable vertical slice, not a total rewrite of
every historical module.

## Current Package Layout

```text
src/aria_tools/
  __main__.py              python -m aria_tools entry point
  cli/
    app.py                 top-level aria-tools router
    common.py              shared CLI helpers
    types.py               light CLI typing support
  commands/
    *.py                   modern command wrappers
    _legacy.py             compatibility helpers for legacy-backed commands
    legacy/
      extract.py           extract parser + validation + delegation
      timeseries.py        timeseries parser + delegation
  core/
    workflows.py           shared extract/timeseries orchestration
  execution/
    workers.py             native Python worker helpers
  config/
    logging.py             modern logging setup

tools/ARIAtools/
  product.py               product loading and metadata-heavy legacy core
  extractProduct.py        extraction/export orchestration legacy core
  util/                    legacy shared utilities

tools/bin/
  aria*.py                 compatibility-era script entry points
```

## Command Flow

### Modern command path

```text
aria-tools extract ...
  -> aria_tools.cli.app
  -> aria_tools.commands.extract
  -> aria_tools.commands.legacy.extract
  -> aria_tools.core.workflows.run_extract_workflow()
  -> ARIAtools.extractProduct / ARIAtools.product / ARIAtools.util...
```

The same overall pattern applies to `timeseries`.

### Compatibility path

```text
ariaExtract.py ...
  -> legacy script entry point
  -> legacy parser/main
  -> aria_tools.core.workflows.run_extract_workflow()
```

That means the modern and legacy entrypoints now share more of the same
underlying orchestration surface than they did before the refactor.

## Responsibility Boundaries

### `aria_tools.cli`

Handles:

- top-level parser construction
- shared flags such as `--log-level` and `--workdir`
- command dispatch
- top-level expected-error rendering

### `aria_tools.commands`

Handles:

- command-specific routing
- compatibility shims for legacy-backed commands
- lightweight validation and argv forwarding

These modules should stay thin.

### `aria_tools.core.workflows`

Handles:

- extract/timeseries workflow setup
- runlog creation
- product loading
- bbox merge logic
- DEM and mask preparation
- shared stack-generation orchestration

This is the main importable seam introduced by the modernization work.

### `aria_tools.execution.workers`

Handles:

- runtime worker-count normalization
- process-pool execution
- simple progress callback integration

The supported parallelization story is now the native Python process backend.

### `ARIAtools.*`

Still handles a lot of domain-heavy behavior, especially:

- product discovery and metadata extraction
- layer export orchestration
- GDAL/VRT manipulation
- correction handling
- detailed legacy workflow semantics

These modules are still real dependencies of the modern path.

## Execution Model

The runtime no longer depends on GNU `parallel`.

Current supported story:

- serial mode for intentionally single-worker paths
- Python process-pool execution for multi-worker export paths

That execution logic is exercised through:

- `src/aria_tools/execution/workers.py`
- `tools/ARIAtools/extractProduct.py`
- unit coverage for backend dispatch and worker behavior

## Testing Model

The modern test structure is:

```text
tests/unit/
tests/integration/
tests/regression/
```

The quick lanes are install-first and offline by default. Marker-gated paths
exist for slower or environment-dependent checks:

- `offline`
- `slow`
- `network_required`
- `credentialed`

This keeps the fast modernization surface easy to run while preserving room for
heavier workflow validation.

## What Is Fully Modernized vs Transitional

### Strongly modernized

- top-level `aria-tools` command router
- modern command wrappers
- extract/timeseries shared workflow layer
- native Python execution helpers
- install-first unit/integration test lanes
- lint-oriented CI structure

### Still transitional

- large portions of `product.py`
- large portions of `extractProduct.py`
- most secondary commands beyond the core extract/timeseries slice
- full live parity coverage for networked and credentialed workflows

## Design Intent

The modernization strategy is deliberately incremental:

- prefer vertical slices over broad scaffolding
- keep compatibility while reducing risk
- extract reusable seams only where they buy testing or clarity
- avoid blocking release readiness on perfect decomposition

## Near-Term Follow-On Work

Good candidates after the current v2 milestone:

- deeper extraction from `product.py`
- deeper extraction from `extractProduct.py`
- stronger live validation for `download`, `order`, and virtual-access paths
- broader type coverage on extracted modules
- gradual reduction of legacy wrapper reliance for secondary commands
