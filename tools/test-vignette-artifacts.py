#!/usr/bin/env python3
"""Integration checks for per-vignette artifact planning."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from colab_dependencies import (
    BASE_R_PACKAGES,
    COMMON_COLAB_PACKAGES,
    description_hard_dependencies,
    referenced_r_packages,
)


def run(
    repo: Path,
    *args: str,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        args,
        cwd=repo,
        check=False,
        text=True,
        capture_output=True,
        env=env,
    )
    if check and result.returncode != 0:
        raise AssertionError(
            f"Command failed ({result.returncode}): {' '.join(args)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def pkgdown_template_value(config_text: str, key: str) -> str | None:
    in_template = False
    key_pattern = re.compile(
        rf"^  {re.escape(key)}:\s*(?P<value>[^#]*?)(?:\s+#.*)?$"
    )
    for line in config_text.splitlines():
        stripped = line.strip()
        if not in_template:
            if re.fullmatch(r"template:\s*(?:#.*)?", line):
                in_template = True
            continue
        if not stripped or stripped.startswith("#"):
            continue
        if not line[0].isspace():
            break
        match = key_pattern.fullmatch(line)
        if match:
            return match.group("value").strip().strip("\"'")
    return None


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


def commit(repo: Path, message: str) -> str:
    run(repo, "git", "add", "-A")
    run(repo, "git", "commit", "-qm", message)
    return run(repo, "git", "rev-parse", "HEAD").stdout.strip()


def plan(
    repo: Path,
    target: str,
    artifact_dir: Path,
    base: str,
    *,
    cache_dir: Path | None = None,
    ignore_missing: bool = True,
) -> list[str]:
    args = [
        "Rscript",
        "tools/vignette_artifacts.R",
        "plan",
        "--target",
        target,
        "--artifact-dir",
        str(artifact_dir),
        "--base",
        base,
        "--sep",
        "space",
    ]
    if cache_dir is not None:
        args.extend(["--cache-dir", str(cache_dir)])
    if ignore_missing:
        args.append("--ignore-missing")
    return run(repo, *args).stdout.split()


def assert_plan(actual: list[str], expected: list[str], label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected}, got {actual}")


def make_fixture(source: Path, destination: Path) -> None:
    destination.mkdir()
    for filename in ["DESCRIPTION", "_pkgdown.yml"]:
        shutil.copy2(source / filename, destination / filename)
    for dirname in ["R", "inst", "vignettes"]:
        shutil.copytree(source / dirname, destination / dirname)

    tools = destination / "tools"
    tools.mkdir()
    for filename in [
        "build-colab-notebook.py",
        "colab_dependencies.py",
        "render-vignette-pdf.R",
        "test-vignette-artifacts.py",
        "validate-empirical-artifacts.py",
        "vignette_artifacts.R",
    ]:
        shutil.copy2(source / "tools" / filename, tools / filename)

    workflows = destination / ".github/workflows"
    workflows.mkdir(parents=True)
    shutil.copy2(
        source / ".github/workflows/pkgdown.yml",
        workflows / "pkgdown.yml",
    )

    run(destination, "git", "init", "-q")
    run(destination, "git", "config", "user.email", "test@example.com")
    run(destination, "git", "config", "user.name", "Artifact Test")
    commit(destination, "baseline")


def main() -> None:
    source = Path(__file__).resolve().parents[1]
    all_slugs = run(
        source, "Rscript", "tools/vignette_artifacts.R", "slugs", "--sep", "space"
    ).stdout.split()
    required_simulation_slugs = {
        "simulation-study-part-1",
        "simulation-study-part-2",
    }
    missing_simulation_slugs = required_simulation_slugs - set(all_slugs)
    if missing_simulation_slugs:
        raise AssertionError(
            "missing split simulation vignettes: "
            f"{sorted(missing_simulation_slugs)}"
        )
    if "simulation-study-vignette" in all_slugs:
        raise AssertionError("superseded simulation-study-vignette slug remains")

    with tempfile.TemporaryDirectory(prefix="bifrost-r-profile-") as temp:
        profile = Path(temp) / "Rprofile"
        profile.write_text(
            "invisible(utils::capture.output(trace(\n"
            "  base::loadNamespace,\n"
            "  tracer = quote(if (identical(package, 'yaml')) {\n"
            "    stop('yaml hidden for artifact test')\n"
            "  }),\n"
            "  print = FALSE\n"
            ")))\n"
            "requireNamespace <- function(package, ...) {\n"
            "  if (identical(package, 'yaml')) return(FALSE)\n"
            "  base::requireNamespace(package, ...)\n"
            "}\n"
        )
        env = os.environ.copy()
        env["R_PROFILE_USER"] = str(profile)
        fallback_slugs = run(
            source,
            "Rscript",
            "tools/vignette_artifacts.R",
            "slugs",
            "--sep",
            "space",
            env=env,
        ).stdout.split()
        assert_plan(fallback_slugs, all_slugs, "slug discovery without yaml")

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

    with tempfile.TemporaryDirectory(prefix="bifrost-artifact-tests-") as temp:
        repo = Path(temp) / "repo"
        make_fixture(source, repo)

        fixture_manifest_path = repo / "inst/extdata/empirical-artifacts.json"
        fixture_manifest = json.loads(fixture_manifest_path.read_text())
        simulation_record = next(
            artifact
            for artifact in fixture_manifest["artifacts"]
            if artifact["path"]
            == "inst/extdata/simulation-study-cache/passerine_preview_tables.rds"
        )
        simulation_record["source_location"] = (
            "Schema-2 cache constructed from validated empirical grids."
        )
        simulation_record["transformation"]["method"] = (
            "Construct the schema-2 simulation vignette cache."
        )
        fixture_manifest_path.write_text(
            json.dumps(fixture_manifest, indent=2, ensure_ascii=False) + "\n"
        )
        stale_schema = run(
            repo,
            "python3",
            "tools/validate-empirical-artifacts.py",
            check=False,
        )
        if stale_schema.returncode == 0:
            raise AssertionError(
                "empirical artifact validator must reject stale schema-2 cache metadata"
            )
        shutil.copy2(
            source / "inst/extdata/empirical-artifacts.json",
            fixture_manifest_path,
        )

        base = run(repo, "git", "rev-parse", "HEAD").stdout.strip()
        fake_bin = Path(temp) / "fake-bin"
        fake_bin.mkdir()
        fake_git = fake_bin / "git"
        fake_git.write_text(
            "#!/bin/sh\n"
            "printf '%s\\n' 'vignettes/quick-start-vignette.Rmd'\n"
            "printf '%s\\n' 'R/warning-from-stderr.R' >&2\n"
        )
        fake_git.chmod(0o755)
        fake_env = os.environ.copy()
        fake_env["PATH"] = f"{fake_bin}{os.pathsep}{fake_env['PATH']}"
        fake_profile = Path(temp) / "fake-git.Rprofile"
        fake_profile.write_text(
            f"Sys.setenv(PATH = {json.dumps(fake_env['PATH'])})\n"
        )
        fake_env["R_PROFILE_USER"] = str(fake_profile)
        stderr_probe = run(
            repo,
            "Rscript",
            "tools/vignette_artifacts.R",
            "plan",
            "--target",
            "colab",
            "--artifact-dir",
            "vignettes/colab",
            "--base",
            "probe-base",
            "--ignore-missing",
            "--sep",
            "space",
            env=fake_env,
        ).stdout.split()
        assert_plan(stderr_probe, ["quick-start-vignette"], "Git stderr isolation")

        pdf_dir = Path(temp) / "pdfs"
        cache_dir = Path(temp) / "cache"
        pdf_dir.mkdir()
        for slug in all_slugs:
            (pdf_dir / f"{slug}.pdf").write_bytes(b"%PDF-artifact-test")
        run(
            repo,
            "Rscript",
            "tools/vignette_artifacts.R",
            "update-cache",
            "--pdf-dir",
            str(pdf_dir),
            "--cache-dir",
            str(cache_dir),
        )
        deleted_asset = repo / "vignettes/jaw-shape/IC_decay.png"
        deleted_asset.unlink()
        assert_plan(
            plan(repo, "pdf", pdf_dir, "", cache_dir=cache_dir),
            ["jaw-shape-vignette"],
            "deleted dependency invalidates manifest",
        )
        shutil.copy2(source / "vignettes/jaw-shape/IC_decay.png", deleted_asset)

        pkgdown_workflow = repo / ".github/workflows/pkgdown.yml"
        with pkgdown_workflow.open("a") as handle:
            handle.write("\n# Artifact cache workflow-change probe.\n")
        assert_plan(
            plan(repo, "pdf", pdf_dir, "", cache_dir=cache_dir),
            all_slugs,
            "production workflow invalidates all cached PDFs",
        )
        shutil.copy2(source / ".github/workflows/pkgdown.yml", pkgdown_workflow)

        with (repo / "vignettes/quick-start-vignette.Rmd").open("a") as handle:
            handle.write("\nArtifact planner source-change probe.\n")
        commit(repo, "change one vignette")
        assert_plan(
            plan(repo, "colab", repo / "vignettes/colab", base),
            ["quick-start-vignette"],
            "single vignette source",
        )

        base = run(repo, "git", "rev-parse", "HEAD").stdout.strip()
        with (repo / "vignettes/jaw-shape/branch_rates.png").open("ab") as handle:
            handle.write(b"artifact-test")
        commit(repo, "change shared vignette asset")
        assert_plan(
            plan(repo, "colab", repo / "vignettes/colab", base),
            ["jaw-shape-vignette", "rate-map-jaw-shape-vignette"],
            "shared vignette asset",
        )

        base = run(repo, "git", "rev-parse", "HEAD").stdout.strip()
        (repo / "vignettes/jaw-shape/IC_decay.png").unlink()
        commit(repo, "delete referenced vignette asset")
        assert_plan(
            plan(repo, "colab", repo / "vignettes/colab", base),
            ["jaw-shape-vignette"],
            "deleted vignette asset",
        )

        base = run(repo, "git", "rev-parse", "HEAD").stdout.strip()
        with (repo / "vignettes/colab/quick-start-vignette.ipynb").open("a") as handle:
            handle.write("\n")
        commit(repo, "change generated notebook")
        assert_plan(
            plan(repo, "colab", repo / "vignettes/colab", base),
            ["quick-start-vignette"],
            "stale generated notebook",
        )
        assert_plan(
            plan(repo, "pdf", repo / "pdf-artifacts", base),
            [],
            "notebook change does not select PDF",
        )

        base = run(repo, "git", "rev-parse", "HEAD").stdout.strip()
        package_file = sorted((repo / "R").glob("*.R"))[0]
        with package_file.open("a") as handle:
            handle.write("\n# Artifact planner shared-change probe.\n")
        commit(repo, "change shared package source")
        assert_plan(
            plan(repo, "colab", repo / "vignettes/colab", base),
            all_slugs,
            "shared package input",
        )

        invalid = run(
            repo,
            "Rscript",
            "tools/vignette_artifacts.R",
            "plan",
            "--target",
            "colab",
            "--artifact-dir",
            str(repo / "vignettes/colab"),
            "--base",
            "does-not-exist",
            "--ignore-missing",
            check=False,
        )
        if invalid.returncode == 0:
            raise AssertionError("invalid Git base must fail artifact planning")

        base = run(repo, "git", "rev-parse", "HEAD").stdout.strip()
        (repo / "vignettes/colab/quick-start-vignette.ipynb").unlink()
        assert_plan(
            plan(
                repo,
                "colab",
                repo / "vignettes/colab",
                base,
                ignore_missing=False,
            ),
            ["quick-start-vignette"],
            "missing generated notebook",
        )

    workflow = (source / ".github/workflows/pkgdown.yml").read_text()
    for pattern in ["vignettes/**/*.rds", "vignettes/**/*.RDS"]:
        if pattern not in workflow:
            raise AssertionError(f"PDF cache key is missing {pattern}")
    if "'.github/workflows/pkgdown.yml'" not in workflow:
        raise AssertionError("PDF cache key must include its production workflow")
    if "'tools/colab_dependencies.py'" not in workflow:
        raise AssertionError("PDF cache key must include Colab dependency detection")
    if "  pull_request:\n" not in workflow:
        raise AssertionError("pkgdown workflow must build pull requests")
    upload_gate = (
        "      - name: Upload site artifact for Pages\n"
        "        if: github.event_name != 'pull_request'\n"
    )
    if upload_gate not in workflow:
        raise AssertionError(
            "pkgdown Pages artifact upload must be disabled for pull requests"
        )

    deploy_gate = (
        "  deploy:\n"
        "    if: github.event_name != 'pull_request'\n"
    )
    if deploy_gate not in workflow:
        raise AssertionError(
            "pkgdown deploy job must be disabled for pull requests"
        )

    check_workflow = (source / ".github/workflows/R-CMD-check.yaml").read_text()
    if '      NOT_CRAN: "false"' not in check_workflow:
        raise AssertionError("R CMD checks must exercise the CRAN-style cheap-test path")

    coverage_workflow = (source / ".github/workflows/test-coverage.yaml").read_text()
    coverage_gate = (
        "          coverage <- covr::percent_coverage(cov)\n"
        "          uncovered <- covr::zero_coverage(cov)\n"
        "          if (coverage != 100 || nrow(uncovered) != 0L) {\n"
        "            stop("
    )
    if coverage_gate not in coverage_workflow:
        raise AssertionError("coverage workflow must enforce 100% coverage with no uncovered rows")
    coverage_stop = (
        '            stop("Coverage gate failed: require exactly 100% coverage '
        'and zero uncovered rows.")'
    )
    coverage_stop_position = coverage_workflow.find(coverage_stop)
    coverage_report_positions = {
        "coverage summary": coverage_workflow.find("          print(cov)"),
        "Cobertura report": coverage_workflow.find(
            "          covr::to_cobertura(cov)"
        ),
    }
    late_or_missing_reports = [
        label
        for label, position in coverage_report_positions.items()
        if position == -1 or position > coverage_stop_position
    ]
    if coverage_stop_position == -1 or late_or_missing_reports:
        raise AssertionError(
            "coverage workflow must print coverage and write Cobertura before "
            "the failure gate stops: "
            + ", ".join(late_or_missing_reports)
        )
    codecov_step = (
        "      - uses: codecov/codecov-action@v6\n"
        "        if: always()\n"
    )
    if codecov_step not in coverage_workflow:
        raise AssertionError("Codecov upload must run after a failed coverage gate")

    renderer = (source / "tools/render-vignette-pdf.R").read_text()
    if "rmarkdown::resolve_output_format(" not in renderer:
        raise AssertionError("PDF renderer must resolve each vignette's YAML format")

    pkgdown_config = (source / "_pkgdown.yml").read_text()
    pkgdown_math_parser_cases = (
        ("template:\n  math-rendering: katex\n", "katex"),
        ("template:\n  math-rendering: 'katex' # renderer\n", "katex"),
        ("template:\n  # math-rendering: katex\n", None),
        (
            "template:\n"
            "  includes:\n"
            "    in_header: |\n"
            "      math-rendering: katex\n"
            "navbar:\n"
            "  math-rendering: katex\n",
            None,
        ),
    )
    for config_text, expected in pkgdown_math_parser_cases:
        actual = pkgdown_template_value(config_text, "math-rendering")
        if actual != expected:
            raise AssertionError(
                "pkgdown template parser returned "
                f"{actual!r}; expected {expected!r}"
            )
    math_renderer = pkgdown_template_value(pkgdown_config, "math-rendering")
    if math_renderer != "katex":
        raise AssertionError(
            "pkgdown config must use KaTeX so equations render across reference "
            "pages and articles"
        )
    if "bifrost.goatcounter.com/count" not in pkgdown_config:
        raise AssertionError("pkgdown config must inline the GoatCounter header include")
    if (source / "pkgdown/extra-head.html").exists():
        raise AssertionError("unused pkgdown header include must not remain tracked")
    if "const encodedSlug = encodeURIComponent(slug);" not in pkgdown_config:
        raise AssertionError("pkgdown artifact links must URL-encode vignette slugs")
    if "if (sourcePath.includes('/articles/')) return;" not in pkgdown_config:
        raise AssertionError(
            "pkgdown artifact links must skip website-only articles"
        )
    if "pdf.href = './' + encodedSlug + '.pdf';" not in pkgdown_config:
        raise AssertionError("pkgdown PDF link must use the encoded vignette slug")
    if "encodedSlug + '.ipynb';" not in pkgdown_config:
        raise AssertionError("pkgdown Colab link must use the encoded vignette slug")
    if "actions.setAttribute('role', 'group');" not in pkgdown_config:
        raise AssertionError("pkgdown artifact action label must describe a group")

    artifact_tool = (source / "tools/vignette_artifacts.R").read_text()
    if '"MISSING"' not in artifact_tool:
        raise AssertionError("artifact hashes must mark missing dependency paths")
    if '"tools/colab_dependencies.py"' not in artifact_tool:
        raise AssertionError("artifact hashes must include Colab dependency detection")

    part2_source = (source / "vignettes/avian-skeleton-part-2.Rmd").read_text()
    part2_widget = part2_source[
        part2_source.index("<figure id=\"lineage-decay-widget-part2\"") :
        part2_source.index("</figure>")
    ]
    indented_block_tags = re.findall(
        r"(?mi)^[ \t]+</?(?:div|figure|script|style|details|summary|svg|p|br)\b",
        part2_widget,
    )
    if "~~~{=html}" in part2_source or indented_block_tags:
        raise AssertionError(
            "Part 2 must emit its HTML widget directly and omit leading "
            "indentation so Pandoc can match the widget's block-level Div tags"
        )

    manifest_validator = source / "tools/validate-empirical-artifacts.py"
    if not manifest_validator.exists():
        raise AssertionError("empirical artifact checksum validator is missing")
    run(source, "python3", str(manifest_validator))

    pr_workflow = (source / ".github/workflows/vignette-artifacts.yml").read_text()
    if "\npermissions:\n  contents: write\n" in pr_workflow:
        raise AssertionError("PR artifact workflow must not grant write access globally")
    if "  update-colab:\n" not in pr_workflow:
        raise AssertionError("PR artifact workflow must isolate Colab updates in a job")
    if "    permissions:\n      contents: write\n" not in pr_workflow:
        raise AssertionError("Colab update job must declare its write permission locally")
    if "ref: ${{ github.event.pull_request.head.sha }}" not in pr_workflow:
        raise AssertionError("PR artifact checks must pin checkout to the event SHA")
    if "ref: ${{ github.event.pull_request.head.ref }}" in pr_workflow:
        raise AssertionError("PR artifact checks must not checkout a mutable branch ref")
    if "      - tools/colab_dependencies.py" not in pr_workflow:
        raise AssertionError("PR artifact workflow must watch Colab dependency detection")
    if "      - tools/validate-empirical-artifacts.py" not in pr_workflow:
        raise AssertionError("PR artifact workflow must watch the artifact validator")
    if "      - tools/avian-skeleton/**" not in pr_workflow:
        raise AssertionError("PR artifact workflow must watch artifact generators")
    generate_step_name = "      - name: Generate changed Colab notebooks"
    audit_step_name = "      - name: Test vignette artifacts"
    generate_step = pr_workflow.find(generate_step_name)
    audit_step = pr_workflow.find(audit_step_name)
    missing_steps = [
        name.strip()
        for name, position in (
            (generate_step_name, generate_step),
            (audit_step_name, audit_step),
        )
        if position == -1
    ]
    if missing_steps:
        raise AssertionError(
            "PR artifact workflow is missing required steps: "
            + ", ".join(missing_steps)
        )
    if audit_step < generate_step:
        raise AssertionError(
            "PR artifact workflow must audit dependencies after notebook generation"
        )

    required_pdf_step = "      - name: Smoke-render all manuscript vignette PDFs"
    if required_pdf_step not in pr_workflow:
        raise AssertionError(
            "PR artifact workflow must smoke-render all seven manuscript vignettes"
        )
    required_pdf_slugs = {
        "avian-skeleton-part-1",
        "avian-skeleton-part-2",
        "avian-skeleton-part-3",
        "avian-skeleton-part-4",
        "avian-skeleton-part-5",
        "simulation-study-part-1",
        "simulation-study-part-2",
    }
    missing_pdf_slugs = sorted(
        slug for slug in required_pdf_slugs if slug not in pr_workflow
    )
    if missing_pdf_slugs:
        raise AssertionError(
            "PR artifact workflow is missing required PDF renders: "
            + ", ".join(missing_pdf_slugs)
        )
    for slug in required_pdf_slugs:
        rmd_text = (source / "vignettes" / f"{slug}.Rmd").read_text()
        if re.search(r'fig\.cap\s*=\s*"\*\*Figure', rmd_text):
            raise AssertionError(
                f"{slug} duplicates the renderer's figure number in fig.cap"
            )

    print("Vignette artifact integration checks passed.")


if __name__ == "__main__":
    main()
