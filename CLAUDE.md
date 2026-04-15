# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Environment Setup

This project runs inside a conda environment called `SAT`. On the cluster, load the module and activate before running anything:

```bash
module load miniforge/25.9.1
eval "$(conda shell.bash hook)"
conda activate SAT
```

After activation, always use `$CONDA_PREFIX/bin/python` instead of bare `python` — the system `python` on PATH may point to the base miniforge Python (3.12) rather than the SAT environment Python (3.10).

- Conda environment: `SAT`
- Environment location: `/mnt/labs/home/jnomburg/.conda/envs/SAT`
- Python version: 3.10.5

---

## Common Commands

**Run the tool:**
```bash
$CONDA_PREFIX/bin/python sat/sat.py <subcommand> [args...]
# Example:
$CONDA_PREFIX/bin/python sat/sat.py aln_filter -a alignment.m8 -o filtered.m8
```

**Run all tests (from project root):**
```bash
cd /mnt/labs/data/nomburg/shared/code/SAT2
$CONDA_PREFIX/bin/python -m pytest
# Or via poetry (used in CI):
poetry run pytest
```

**Run a single test file:**
```bash
$CONDA_PREFIX/bin/python -m pytest tests/test_aln_filter.py
```

**Skip slow ete3 tests (require network/database):**
```bash
$CONDA_PREFIX/bin/python -m pytest -m "not ete3"
```

**Lint and format:**
```bash
flake8 sat/
black sat/
```

---

## Architecture

SAT is a CLI toolkit for structural biology analysis. The entry point is `sat/sat.py`, which uses `argparse` subparsers to dispatch to individual script modules.

### Dispatch pattern

`sat/sat.py` has two sections:
1. **Top half** — one `subparsers.add_parser(...)` block per subcommand, in alphabetical order.
2. **Bottom half** — one `call_<subcommand>(args)` function per subcommand, in alphabetical order, each doing a lazy import and calling `<subcommand>_main(args)`. Imports use `from scripts.<name>` (no leading dot).

### Script files (`sat/scripts/`)

Each subcommand lives in its own file: `sat/scripts/<subcommand_name>.py`. Key conventions:
- Imports from local utils use **relative** imports: `from .utils.misc import talk_to_me, make_output_dir`
- The entry function is named `<subcommand_name>_main(args)` and accepts the argparse `args` namespace.
- `talk_to_me(msg)` is used for all progress output.
- `make_output_dir(path)` is called before writing any output file.
- Files end with an `if __name__ == "__main__":` guard that raises `ValueError`.

### Subcommand naming

Subcommand names use snake_case with a category prefix:
- `aln_` — alignment-focused (Foldseek/DALI output)
- `struc_` — structure-focused (PDB manipulation)
- `seq_` — sequence-focused (FASTA manipulation)
- `tab_` — tabular file manipulation
- `plot_` — plotting

The script filename, main function, and argparse subcommand string must all match (e.g. `struc_rebase` → `struc_rebase.py` → `struc_rebase_main` → `call_struc_rebase`).

### Utilities (`sat/scripts/utils/`)

Shared logic lives here:
- `misc.py` — `talk_to_me`, `make_output_dir`, `read_fasta_to_memory`, `read_tsv` (generator), `arg_str2bool`
- `alignments.py` — alignment parsing helpers
- `structure.py` — PDB/structure manipulation
- `clusters.py` — clustering logic
- `dali.py` — DALI alignment parsing
- `ete3_taxonomy.py` — taxonomy lookups (requires network; tests marked `ete3`)
- `uniprot.py` — UniProt API queries
- `Foldseek_Dataset.py` — Foldseek result handling

### Tests

Tests live in `tests/` and import directly from `sat.scripts.<subcommand>`. They typically instantiate a plain `Args` class (not argparse) and call `<subcommand>_main(args)` directly. Test data lives in `tests/test_data/`.

---

## Adding a New Subcommand

Three files to modify, one to create — all insertions in **alphabetical order**:

1. **Create** `sat/scripts/<name>.py` — follow the script structure above.
2. **Modify** `sat/sat.py` — add the parser block (top half) and `call_<name>` function (bottom half).
3. **Modify** `README.md` — add a one-liner in the overview list and a detailed section (both alphabetical within their section).

The detailed README section must include the RICH-CODEX screenshot directive:
```markdown
<!-- RICH-CODEX hide_command: true -->
![`poetry run .github/tmp/sat_codex.py <name> -h`](.github/img/<name>.png)
```

### Argparse conventions

- Each argument uses both `-x` (short) and `--long_name` (long) flags.
- `help` strings use triple-quoted multiline strings.
- Required args: `required=True`, no `default`.
- Optional args: `required=False` with a `default`. String defaults meaning "not provided" use `default=""`.
- Boolean flags use `action="store_true"` (no `type`, `required`, or `default`).
- Block ends with `.set_defaults(func=call_<name>)`.
