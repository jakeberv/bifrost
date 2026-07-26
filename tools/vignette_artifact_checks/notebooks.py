"""Notebook generation, normalization, freshness, and Colab policy checks."""

from __future__ import annotations

import json
import re
import shutil
import tempfile
from pathlib import Path

try:
    from colab_dependencies import (
        BASE_R_PACKAGES,
        COMMON_COLAB_PACKAGES,
        description_hard_dependencies,
        referenced_r_packages,
    )
except ModuleNotFoundError:
    from tools.colab_dependencies import (
        BASE_R_PACKAGES,
        COMMON_COLAB_PACKAGES,
        description_hard_dependencies,
        referenced_r_packages,
    )
from .common import run

__all__ = ["run_notebook_checks"]


def bootstrapped_packages(setup: str) -> set[str]:
    match = re.search(r"colab_packages\s*<-\s*c\((.*?)\)", setup, re.DOTALL)
    if not match:
        return set()
    return set(re.findall(r'["\']([A-Za-z][A-Za-z0-9.]*)["\']', match.group(1)))


def missing_notebook_dependencies(
    repo: Path,
    notebook: dict,
    hard_dependencies: set[str],
) -> set[str]:
    setup = notebook["cells"][1]["source"]
    executable_code = "\n".join(
        cell["source"]
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
    referenced = referenced_r_packages(repo, executable_code)
    available = (
        BASE_R_PACKAGES
        | hard_dependencies
        | bootstrapped_packages(setup)
        | {"bifrost"}
    )
    return referenced - available


def run_clean_notebook_setup_probe(
    repo: Path,
    notebook: dict,
    *,
    setup_token: str,
    assertions: str,
    label: str,
) -> None:
    setup_cells = [
        cell["source"]
        for cell in notebook["cells"]
        if cell["cell_type"] == "code" and setup_token in cell["source"]
    ]
    if len(setup_cells) != 1:
        raise AssertionError(
            f"{label} expected one setup cell containing {setup_token!r}, "
            f"found {len(setup_cells)}"
        )
    with tempfile.TemporaryDirectory(prefix="bifrost-colab-clean-session-") as temp:
        probe = Path(temp) / "probe.R"
        probe.write_text(setup_cells[0] + "\n" + assertions + "\n")
        run(repo, "Rscript", "--vanilla", str(probe))


def run_notebook_reporting_probe(
    repo: Path,
    notebook: dict,
    *,
    fixture: str,
    reporting_tokens: tuple[str, ...],
    expected_output: tuple[str, ...],
    label: str,
) -> None:
    setup_cells = [
        cell["source"]
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
        and "render_preview_table <- function" in cell["source"]
    ]
    if len(setup_cells) != 1:
        raise AssertionError(f"{label} expected one reporting helper cell")

    reporting_cells = []
    for token in reporting_tokens:
        matches = [
            cell["source"]
            for cell in notebook["cells"]
            if cell["cell_type"] == "code" and token in cell["source"]
        ]
        if len(matches) != 1:
            raise AssertionError(
                f"{label} expected one reporting cell containing {token!r}"
            )
        reporting_cells.append(matches[0])

    assignments = []
    for index, source in enumerate(reporting_cells, start=1):
        assignments.append(f"reported_{index} <- local({{\n{source}\n}})")
    reported_names = ", ".join(
        f"reported_{index}" for index in range(1, len(reporting_cells) + 1)
    )
    expected_vector = ", ".join(json.dumps(value) for value in expected_output)
    assertions = (
        f"rendered <- paste(capture.output(print(list({reported_names}))), "
        'collapse = "\\n")\n'
        f"expected <- c({expected_vector})\n"
        "stopifnot(all(vapply(expected, grepl, logical(1L), "
        "x = rendered, fixed = TRUE)))\n"
    )

    with tempfile.TemporaryDirectory(prefix="bifrost-colab-reporting-") as temp:
        probe = Path(temp) / "probe.R"
        probe.write_text(
            setup_cells[0]
            + "\n"
            + fixture
            + "\n"
            + "\n".join(assignments)
            + "\n"
            + assertions
        )
        run(repo, "Rscript", "--vanilla", str(probe))


def run_notebook_checks(source: Path, all_slugs: list[str]) -> None:
    required_simulation_slugs = {
        "simulation-study-part-1",
        "simulation-study-part-2",
    }
    with tempfile.TemporaryDirectory(prefix="bifrost-colab-policy-") as temp:
        repo = Path(temp) / "repo"
        (repo / "tools").mkdir(parents=True)
        (repo / "vignettes").mkdir()
        shutil.copy2(source / "DESCRIPTION", repo / "DESCRIPTION")
        shutil.copy2(source / "_pkgdown.yml", repo / "_pkgdown.yml")
        shutil.copy2(
            source / "tools/build-colab-notebook.py",
            repo / "tools/build-colab-notebook.py",
        )
        shutil.copy2(
            source / "tools/colab_dependencies.py",
            repo / "tools/colab_dependencies.py",
        )
        (repo / "vignettes/colab-policy-probe.Rmd").write_text(
            "---\n"
            'title: "Colab policy probe"\n'
            'colab-packages: "manualPackage, manualPackageTwo"\n'
            "---\n\n"
            "```{r setup, include=FALSE}\n"
            "knitr::opts_chunk$set(eval = FALSE)\n"
            "evaluated_setup_value <- 41L\n"
            "evaluated_setup_helper <- function() {\n"
            "  if (knitr::is_html_output()) 1L else 1L\n"
            "}\n"
            "```\n\n"
            "```{r visible, eval=FALSE}\n"
            "visible_probe <- TRUE\n"
            "library(autoLibrary)\n"
            'requireNamespace("autoRequired", quietly = TRUE)\n'
            'loadNamespace("autoLoaded")\n'
            "autoNamespace::run()\n"
            "ape::Ntip(NULL)\n"
            "stats::lm(visible_probe ~ 1)\n"
            "```\n\n"
            "```{r dependent}\n"
            "evaluated_setup_result <- evaluated_setup_value + evaluated_setup_helper()\n"
            "stopifnot(identical(evaluated_setup_result, 42L))\n"
            "```\n\n"
            "Inline result: `r evaluated_setup_value + 1L`.\n\n"
            "```{r visible-asis, results='asis'}\n"
            "visible_asis_probe <- data.frame(metric = 'BA', score = 0.75)\n"
            "visible_asis_probe\n"
            "```\n\n"
            "```{r hidden, include=FALSE, eval=FALSE}\n"
            "hidden_probe <- TRUE\n"
            "hiddenPackage::run()\n"
            "```\n"
        )
        run(
            repo,
            "python3",
            "tools/build-colab-notebook.py",
            "--slug",
            "colab-policy-probe",
        )
        policy_notebook = json.loads(
            (repo / "vignettes/colab/colab-policy-probe.ipynb").read_text()
        )
        policy_code = "\n".join(
            cell["source"]
            for cell in policy_notebook["cells"]
            if cell["cell_type"] == "code"
        )
        if "visible_probe <- TRUE" not in policy_code:
            raise AssertionError("visible eval=FALSE chunks must be executable in Colab")
        if "evaluated_setup_value <- 41L" not in policy_code:
            raise AssertionError(
                "evaluated include=FALSE setup chunks must remain executable in Colab"
            )
        if "hidden_probe <- TRUE" in policy_code:
            raise AssertionError("hidden unevaluated chunks must be omitted from Colab")
        if "visible_asis_probe <-" not in policy_code:
            raise AssertionError(
                "visible results='asis' chunks must remain executable in Colab"
            )
        policy_markdown = "\n".join(
            cell["source"]
            for cell in policy_notebook["cells"]
            if cell["cell_type"] == "markdown"
        )
        if re.search(r"`r\s+.+?`", policy_markdown):
            raise AssertionError(
                "inline R expressions must not remain literal notebook Markdown"
            )
        if "evaluated_setup_value + 1L" not in policy_code:
            raise AssertionError(
                "inline R expressions must become executable notebook code"
            )
        execution_cells = [
            cell["source"]
            for cell in policy_notebook["cells"]
            if cell["cell_type"] == "code"
            and (
                "evaluated_setup_value <- 41L" in cell["source"]
                or "evaluated_setup_result <-" in cell["source"]
            )
        ]
        execution_probe = Path(temp) / "colab-clean-session-probe.R"
        execution_probe.write_text("\n".join(execution_cells))
        run(repo, "Rscript", "--vanilla", str(execution_probe))
        policy_setup_packages = bootstrapped_packages(
            policy_notebook["cells"][1]["source"]
        )
        expected_policy_packages = {
            "remotes",
            "knitr",
            "autoLibrary",
            "autoLoaded",
            "autoNamespace",
            "autoRequired",
            "manualPackage",
            "manualPackageTwo",
        }
        if policy_setup_packages != expected_policy_packages:
            raise AssertionError(
                "automatic dependency probe expected "
                f"{sorted(expected_policy_packages)}, got "
                f"{sorted(policy_setup_packages)}"
            )

    theoretical_rmd = (
        source / "vignettes/theoretical-background-vignette.Rmd"
    ).read_text()
    theoretical_notebook = (
        source / "vignettes/colab/theoretical-background-vignette.ipynb"
    ).read_text()
    if "set.seed(0.1)" in theoretical_rmd or "set.seed(0.1)" in theoretical_notebook:
        raise AssertionError("theoretical vignette must use an explicit integer seed")

    hard_dependencies = description_hard_dependencies(source / "DESCRIPTION")
    required_colab_packages = {
        "avian-skeleton-part-3": {"univariateML", "evd"},
        "avian-skeleton-part-5": {"phylolm"},
    }
    dependency_probe: dict | None = None
    for notebook_path in sorted((source / "vignettes/colab").glob("*.ipynb")):
        notebook = json.loads(notebook_path.read_text())
        if notebook_path.stem == "quick-start-vignette":
            dependency_probe = json.loads(json.dumps(notebook))
        runtime_note = notebook["cells"][0]["source"]
        for phrase in (
            "Recommended Colab runtime",
            "v5e-1 TPU",
            "host CPUs",
            "runtime with the most CPUs",
        ):
            if phrase not in runtime_note:
                raise AssertionError(
                    f"{notebook_path.name} runtime note is missing {phrase!r}"
                )
        fenced_r_examples = [
            cell
            for cell in notebook["cells"]
            if cell["cell_type"] == "markdown" and "```r\n" in cell["source"]
        ]
        if fenced_r_examples:
            raise AssertionError(
                f"{notebook_path.name} contains R examples that are not executable cells"
            )
        hidden_maintenance_cells = [
            cell
            for cell in notebook["cells"]
            if cell["cell_type"] == "code"
            and "rate_map_save_arc_figure(" in cell["source"]
        ]
        if hidden_maintenance_cells:
            raise AssertionError(
                f"{notebook_path.name} contains hidden vignette maintenance code"
            )
        setup = notebook["cells"][1]["source"]
        missing_required = required_colab_packages.get(
            notebook_path.stem, set()
        ) - bootstrapped_packages(setup)
        if missing_required:
            raise AssertionError(
                f"{notebook_path.name} Colab setup is missing declared runtime "
                f"packages: {', '.join(sorted(missing_required))}"
            )
        if "parallel::detectCores(logical = TRUE)" not in setup:
            raise AssertionError(
                f"{notebook_path.name} setup must report detected logical CPUs"
            )
        if "git clone --depth 1 " not in setup:
            raise AssertionError(
                f"{notebook_path.name} setup must use a shallow Git clone"
            )
        if "dependencies = NA" not in setup or "dependencies = TRUE" in setup:
            raise AssertionError(
                f"{notebook_path.name} setup must install only hard package dependencies"
            )
        if '"knitr"' not in setup:
            raise AssertionError(
                f"{notebook_path.name} setup must install the shared knitr dependency"
            )

        body_code = "\n".join(
            cell["source"]
            for cell in notebook["cells"][2:]
            if cell["cell_type"] == "code"
        )
        body_markdown = "\n".join(
            cell["source"]
            for cell in notebook["cells"][2:]
            if cell["cell_type"] == "markdown"
        )
        required_code = {
            "avian-skeleton-part-2": (
                "library(bifrost)",
                "lineage_decay_widget <- list(",
            ),
            "simulation-study-part-1": ("pkg_file <- function(...)",),
            "simulation-study-part-2": ("pkg_file <- function(...)",),
        }
        for required in required_code.get(notebook_path.stem, ()):
            if required not in body_code:
                raise AssertionError(
                    f"{notebook_path.name} is missing required setup code: {required}"
                )
        simulation_notebook_code = {
            "simulation-study-part-1": (
                "render_preview_table(\n  fixed_null_display",
                "render_preview_table(\n  fixed_recovery_display",
                '"Fuzzy recall"',
                '"Fuzzy specificity"',
                '"Fuzzy F1"',
                '"Fuzzy balanced accuracy"',
            ),
            "simulation-study-part-2": (
                "render_preview_table(\n  gic_preview_table",
                "render_preview_table(\n  bic_preview_table",
                "render_preview_table(\n  preview_recommendations",
                '"Prop. BA"',
                '"Int.-rate BA"',
                '"Score"',
            ),
        }
        for required in simulation_notebook_code.get(notebook_path.stem, ()):
            if required not in body_code:
                raise AssertionError(
                    f"{notebook_path.name} omits executable reporting code: {required}"
                )
        if notebook_path.stem in required_simulation_slugs:
            if re.search(r"`r\s+.+?`", body_markdown):
                raise AssertionError(
                    f"{notebook_path.name} contains literal inline R Markdown"
                )
        if notebook_path.stem == "simulation-study-part-2":
            if body_code.count(
                "scenario_weights = c(proportional = 0.50, correlation = 0.50)"
            ) < 4:
                raise AssertionError(
                    "Part 2 notebook must execute explicit equal scenario weights"
                )
            required_hard_stops = (
                "!gic_preview_tuned$used_all_settings",
                "!bic_preview_tuned$used_all_settings",
                "stopifnot(!gic_tuned$used_all_settings)",
                "stopifnot(!bic_tuned$used_all_settings)",
            )
            for required in required_hard_stops:
                if required not in body_code:
                    raise AssertionError(
                        "Part 2 notebook is missing mandatory selector hard-stop: "
                        + required
                    )
        reporting_probes = {
            "simulation-study-part-1": {
                "fixture": (
                    "fixed_null_display <- data.frame(\n"
                    "  IC = 'GIC', `Mean FP` = 0.0123, `Any FP` = 0.04,\n"
                    "  `Mean shifts` = 0.02, check.names = FALSE\n"
                    ")\n"
                    "fixed_recovery_display <- data.frame(\n"
                    "  IC = 'GIC', Scenario = 'Proportional',\n"
                    "  `Fuzzy recall` = 0.731, `Fuzzy specificity` = 0.887,\n"
                    "  `Fuzzy F1` = 0.809, `Fuzzy balanced accuracy` = 0.809,\n"
                    "  `Mean shifts` = 4.1, check.names = FALSE\n"
                    ")"
                ),
                "reporting_tokens": (
                    "render_preview_table(\n  fixed_null_display",
                    "render_preview_table(\n  fixed_recovery_display",
                ),
                "expected_output": (
                    "Mean FP", "Fuzzy rec.", "Fuzzy spec.", "Fuzzy F1",
                    "Fuzzy BA", "0.731", "0.887", "0.809",
                ),
            },
            "simulation-study-part-2": {
                "fixture": (
                    "gic_preview_table <- data.frame(\n"
                    "  Threshold = 10, `Min clade` = 10, `Null FP` = 0.01,\n"
                    "  `Prop. BA` = 0.811, `Int.-rate BA` = 0.722, Score = 0.765,\n"
                    "  check.names = FALSE\n"
                    ")\n"
                    "bic_preview_table <- gic_preview_table\n"
                    "preview_recommendations <- data.frame(\n"
                    "  IC = 'GIC', Threshold = 10, `Min clade` = 10,\n"
                    "  `Null FP` = 0.01, `Prop. BA` = 0.811,\n"
                    "  `Int.-rate BA` = 0.722, Score = 0.765, check.names = FALSE\n"
                    ")"
                ),
                "reporting_tokens": (
                    "render_preview_table(\n  gic_preview_table",
                    "render_preview_table(\n  bic_preview_table",
                    "render_preview_table(\n  preview_recommendations",
                ),
                "expected_output": (
                    "Prop. BA", "Int.-rate BA", "Score",
                    "0.811", "0.722", "0.765",
                ),
            },
        }
        if notebook_path.stem in reporting_probes:
            probe = reporting_probes[notebook_path.stem]
            run_notebook_reporting_probe(
                source,
                notebook,
                fixture=probe["fixture"],
                reporting_tokens=probe["reporting_tokens"],
                expected_output=probe["expected_output"],
                label=notebook_path.name,
            )
        clean_setup_probes = {
            "avian-skeleton-part-2": (
                "lineage_decay_widget <- list(",
                'stopifnot("package:bifrost" %in% search())\n'
                "stopifnot(identical(lineage_decay_widget$defaults$half_life, 5))",
            ),
            "simulation-study-part-1": (
                "pkg_file <- function(...)",
                "stopifnot(file.exists(pkg_file(\n"
                '  "extdata", "avian-skeleton", "passerine_bodyplan_tree.tre"\n'
                ")))",
            ),
            "simulation-study-part-2": (
                "pkg_file <- function(...)",
                "stopifnot(file.exists(pkg_file(\n"
                '  "extdata", "avian-skeleton", "passerine_bodyplan_tree.tre"\n'
                ")))",
            ),
        }
        if notebook_path.stem in clean_setup_probes:
            setup_token, assertions = clean_setup_probes[notebook_path.stem]
            run_clean_notebook_setup_probe(
                source,
                notebook,
                setup_token=setup_token,
                assertions=assertions,
                label=notebook_path.name,
            )
        expected_optional = referenced_r_packages(source, body_code) - (
            BASE_R_PACKAGES
            | hard_dependencies
            | set(COMMON_COLAB_PACKAGES)
            | {"bifrost"}
        )
        expected_optional |= required_colab_packages.get(notebook_path.stem, set())
        actual_optional = bootstrapped_packages(setup) - set(COMMON_COLAB_PACKAGES)
        if actual_optional != expected_optional:
            raise AssertionError(
                f"{notebook_path.name} expected automatic optional dependencies "
                f"{sorted(expected_optional)}, got {sorted(actual_optional)}"
            )
        missing = missing_notebook_dependencies(source, notebook, hard_dependencies)
        if missing:
            raise AssertionError(
                f"{notebook_path.name} executable cells use packages unavailable in "
                f"Colab setup: {', '.join(sorted(missing))}. Automatic dependency "
                "detection did not add them to the setup cell."
            )

    if dependency_probe is None:
        raise AssertionError("quick-start notebook is required for dependency audit probe")
    dependency_probe["cells"][1]["source"] += "htmltools::tagList()\n"
    dependency_probe["cells"].append(
        {
            "cell_type": "code",
            "source": (
                "plotly::plot_ly()\n"
                "library(RColorBrewer)\n"
                'requireNamespace("geomorph", quietly = TRUE)\n'
            ),
        }
    )
    probe_missing = missing_notebook_dependencies(
        source, dependency_probe, hard_dependencies
    )
    expected_probe_missing = {"RColorBrewer", "geomorph", "htmltools", "plotly"}
    if probe_missing != expected_probe_missing:
        raise AssertionError(
            "dependency audit probe expected "
            f"{sorted(expected_probe_missing)}, got {sorted(probe_missing)}"
        )
