"""Temporary Git fixtures and vignette artifact planner checks."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path

from .common import run

__all__ = ["assert_plan", "commit", "make_fixture", "plan", "run_planning_checks"]


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
    for dirname in ["R", "data-remote", "inst", "vignettes"]:
        shutil.copytree(source / dirname, destination / dirname)

    tools = destination / "tools"
    tools.mkdir()
    for filename in [
        "build-colab-notebook.py",
        "colab_dependencies.py",
        "render-vignette-pdf.R",
        "test-vignette-artifacts.py",
        "validate-empirical-artifacts.py",
        "validate-pkgdown-config.R",
        "vignette_artifacts.R",
    ]:
        shutil.copy2(source / "tools" / filename, tools / filename)
    shutil.copytree(
        source / "tools/vignette_artifact_checks",
        tools / "vignette_artifact_checks",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )

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


def run_planning_checks(source: Path, all_slugs: list[str]) -> None:
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

    with tempfile.TemporaryDirectory(prefix="bifrost-artifact-tests-") as temp:
        repo = Path(temp) / "repo"
        make_fixture(source, repo)

        fixture_manifest_path = repo / "data-remote/empirical-artifacts.json"
        fixture_manifest = json.loads(fixture_manifest_path.read_text())
        # Supplementary files are integrity-checked without extending the
        # package's fixed downloader identifier contract.
        run(repo, "python3", "tools/validate-empirical-artifacts.py")
        for fault in ("checksum", "size", "downloader", "missing", "type"):
            broken = json.loads(json.dumps(fixture_manifest))
            records = broken["supplementary_artifacts"]
            if fault == "checksum":
                records[0]["sha256"] = "0" * 64
            elif fault == "size":
                records[0]["size_bytes"] += 1
            elif fault == "downloader":
                records[0]["artifact_id"] = "not-a-downloader-entry"
            elif fault == "missing":
                broken["supplementary_artifacts"] = []
            else:
                broken["supplementary_artifacts"] = {}
            fixture_manifest_path.write_text(json.dumps(broken) + "\n")
            rejected = run(repo, "python3", "tools/validate-empirical-artifacts.py", check=False)
            if rejected.returncode == 0:
                raise AssertionError(f"supplementary artifact {fault} corruption was accepted")
        fixture_manifest_path.write_text(json.dumps(fixture_manifest) + "\n")
        simulation_record = next(
            artifact
            for artifact in fixture_manifest["artifacts"]
            if artifact["path"]
            == "data-remote/simulation-study-cache/passerine_preview_tables.rds"
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
            source / "data-remote/empirical-artifacts.json",
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
