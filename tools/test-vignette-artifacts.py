#!/usr/bin/env python3
"""Integration checks for per-vignette artifact planning."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


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
        "render-vignette-pdf.R",
        "test-vignette-artifacts.py",
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

    theoretical_rmd = (
        source / "vignettes/theoretical-background-vignette.Rmd"
    ).read_text()
    theoretical_notebook = (
        source / "vignettes/colab/theoretical-background-vignette.ipynb"
    ).read_text()
    if "set.seed(0.1)" in theoretical_rmd or "set.seed(0.1)" in theoretical_notebook:
        raise AssertionError("theoretical vignette must use an explicit integer seed")

    for notebook_path in sorted((source / "vignettes/colab").glob("*.ipynb")):
        notebook = json.loads(notebook_path.read_text())
        setup = notebook["cells"][1]["source"]
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

        expects_phylolm = (
            notebook_path.stem == "pca-model-selection-and-bifrost-vignette"
        )
        if ('"phylolm"' in setup) != expects_phylolm:
            raise AssertionError(
                f"{notebook_path.name} setup has the wrong vignette-specific dependencies"
            )

    with tempfile.TemporaryDirectory(prefix="bifrost-artifact-tests-") as temp:
        repo = Path(temp) / "repo"
        make_fixture(source, repo)

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

    renderer = (source / "tools/render-vignette-pdf.R").read_text()
    if "rmarkdown::resolve_output_format(" not in renderer:
        raise AssertionError("PDF renderer must resolve each vignette's YAML format")

    pkgdown_config = (source / "_pkgdown.yml").read_text()
    if "bifrost.goatcounter.com/count" not in pkgdown_config:
        raise AssertionError("pkgdown config must inline the GoatCounter header include")
    if (source / "pkgdown/extra-head.html").exists():
        raise AssertionError("unused pkgdown header include must not remain tracked")
    if "const encodedSlug = encodeURIComponent(slug);" not in pkgdown_config:
        raise AssertionError("pkgdown artifact links must URL-encode vignette slugs")
    if "pdf.href = './' + encodedSlug + '.pdf';" not in pkgdown_config:
        raise AssertionError("pkgdown PDF link must use the encoded vignette slug")
    if "encodedSlug + '.ipynb';" not in pkgdown_config:
        raise AssertionError("pkgdown Colab link must use the encoded vignette slug")
    if "actions.setAttribute('role', 'group');" not in pkgdown_config:
        raise AssertionError("pkgdown artifact action label must describe a group")

    artifact_tool = (source / "tools/vignette_artifacts.R").read_text()
    if '"MISSING"' not in artifact_tool:
        raise AssertionError("artifact hashes must mark missing dependency paths")

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

    print("Vignette artifact integration checks passed.")


if __name__ == "__main__":
    main()
