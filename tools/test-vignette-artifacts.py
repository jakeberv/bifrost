#!/usr/bin/env python3
"""Run vignette artifact integration checks."""

from __future__ import annotations

from pathlib import Path

from vignette_artifact_checks.common import run
from vignette_artifact_checks.contracts import run_repository_contract_checks
from vignette_artifact_checks.notebooks import run_notebook_checks
from vignette_artifact_checks.planning import run_planning_checks


def main() -> None:
    source = Path(__file__).resolve().parents[1]
    all_slugs = run(
        source, "Rscript", "tools/vignette_artifacts.R", "slugs", "--sep", "space"
    ).stdout.split()
    run_notebook_checks(source, all_slugs)
    run_planning_checks(source, all_slugs)
    run_repository_contract_checks(source, all_slugs)
    print("Vignette artifact integration checks passed.")


if __name__ == "__main__":
    main()
