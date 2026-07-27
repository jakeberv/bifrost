"""Repository, workflow, package, coverage, and provenance contracts."""

from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path

from .common import run

__all__ = ["run_repository_contract_checks"]


def run_repository_contract_checks(source: Path, all_slugs: list[str]) -> None:
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

    workflow = (source / ".github/workflows/pkgdown.yml").read_text()
    for pattern in ["vignettes/**/*.rds", "vignettes/**/*.RDS"]:
        if pattern not in workflow:
            raise AssertionError(f"PDF cache key is missing {pattern}")
    if "'.github/workflows/pkgdown.yml'" not in workflow:
        raise AssertionError("PDF cache key must include its production workflow")
    if "'tools/colab_dependencies.py'" not in workflow:
        raise AssertionError("PDF cache key must include Colab dependency detection")
    for cache_input in (
        "'pkgdown/**'",
        "'tools/vignette_artifact_checks/**'",
        "'tools/validate-pkgdown-config.R'",
    ):
        if cache_input not in workflow:
            raise AssertionError(f"PDF cache key is missing {cache_input}")
    if "  pull_request:\n" not in workflow:
        raise AssertionError("pkgdown workflow must build pull requests")
    browser_smoke_path = "tools/pkgdown-browser-smoke/**"
    if f"'{browser_smoke_path}'" not in workflow:
        raise AssertionError(
            "PDF cache key must include the pkgdown browser smoke suite"
        )
    workflow_preamble = workflow[: workflow.index("\njobs:\n")]
    if "\nconcurrency:\n" in workflow_preamble:
        raise AssertionError(
            "pkgdown workflow must not make PR builds contend for Pages concurrency"
        )
    build_schedule = (
        "  build:\n"
        "    if: >-\n"
        "      github.event_name != 'pull_request' ||\n"
        "      github.actor != 'github-actions[bot]' ||\n"
        "      github.event.pull_request.user.login == 'github-actions[bot]'\n"
        "    concurrency:\n"
        "      group: pkgdown-build-${{ github.event.pull_request.number || "
        "github.ref }}\n"
        "      cancel-in-progress: ${{ github.event_name == 'pull_request' }}\n"
    )
    if build_schedule not in workflow:
        raise AssertionError(
            "pkgdown build job must skip bot notebook follow-ups on human-authored "
            "PRs, allow bot-authored maintenance PRs, and cancel only stale PR builds"
        )
    upload_gate = (
        "      - name: Upload site artifact for Pages\n"
        "        if: github.event_name != 'pull_request'\n"
    )
    if upload_gate not in workflow:
        raise AssertionError(
            "pkgdown Pages artifact upload must be disabled for pull requests"
        )
    pkgdown_dependencies = (
        "      - name: Setup Node\n"
        "        uses: actions/setup-node@v6\n"
        "        with:\n"
        "          node-version: 24\n"
        "          cache: npm\n"
        "          cache-dependency-path: tools/pkgdown-browser-smoke/package-lock.json\n"
        "\n"
        "      - name: Install R dependencies (incl. pkgdown)\n"
        "        uses: r-lib/actions/setup-r-dependencies@v2\n"
        "        with:\n"
        "          extra-packages: any::pkgdown, local::.\n"
    )
    if pkgdown_dependencies not in workflow:
        raise AssertionError(
            "pkgdown workflow must set up the cached Node 24 browser-test "
            "environment before retaining any::pkgdown"
        )
    pkgdown_browser_smoke_steps = (
        "      - name: Build pkgdown site into docs/\n"
        "        shell: Rscript {0}\n"
        "        run: |\n"
        "          pkgdown::build_site_github_pages(\n"
        "            new_process = FALSE,\n"
        "            install = FALSE\n"
        "          )\n"
        "\n"
        "      - name: Install pkgdown browser smoke dependencies\n"
        "        run: npm ci --prefix tools/pkgdown-browser-smoke\n"
        "\n"
        "      - name: Install Chromium for pkgdown browser smoke tests\n"
        "        run: npx --prefix tools/pkgdown-browser-smoke playwright install --with-deps chromium\n"
        "\n"
        "      - name: Test built pkgdown site in Chromium\n"
        "        run: npm test --prefix tools/pkgdown-browser-smoke\n"
    )
    if pkgdown_browser_smoke_steps not in workflow:
        raise AssertionError(
            "pkgdown workflow must run the Chromium browser smoke suite after "
            "building docs/"
        )

    deploy_schedule = (
        "  deploy:\n"
        "    if: github.event_name != 'pull_request'\n"
        "    concurrency:\n"
        "      group: pages\n"
        "      cancel-in-progress: false\n"
    )
    if deploy_schedule not in workflow:
        raise AssertionError(
            "pkgdown deploy job must be disabled for pull requests and serialized"
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

    pkgdown_config_path = source / "_pkgdown.yml"
    pkgdown_config = pkgdown_config_path.read_text()
    pkgdown_validator = source / "tools/validate-pkgdown-config.R"
    pkgdown_validator_cases = (
        ("block template", "template:\n  math-rendering: katex\n", True),
        ("inline template", "template: {math-rendering: katex}\n", True),
        ("missing template", "navbar:\n  structure: [left, right]\n", False),
        ("plural template key", "templates:\n  math-rendering: katex\n", False),
        (
            "extended template key",
            "template-extra:\n  math-rendering: katex\n",
            False,
        ),
        ("missing math-rendering", "template:\n  bootstrap: 5\n", False),
        ("non-katex renderer", "template:\n  math-rendering: mathjax\n", False),
        (
            "inline sequence renderer",
            "template:\n  math-rendering: [katex]\n",
            False,
        ),
        (
            "block sequence renderer",
            "template:\n  math-rendering:\n    - katex\n",
            False,
        ),
        ("malformed YAML", "template: [\n", False),
    )
    with tempfile.TemporaryDirectory(prefix="bifrost-pkgdown-config-") as temp:
        for label, config_text, should_succeed in pkgdown_validator_cases:
            config_path = Path(temp) / f"{label.replace(' ', '-')}.yml"
            config_path.write_text(config_text)
            result = run(
                source,
                "Rscript",
                "--vanilla",
                str(pkgdown_validator),
                str(config_path),
                check=False,
            )
            if (result.returncode == 0) != should_succeed:
                raise AssertionError(
                    f"pkgdown validator case {label!r} returned "
                    f"{result.returncode}; stdout:\n{result.stdout}\n"
                    f"stderr:\n{result.stderr}"
                )
    run(
        source,
        "Rscript",
        "--vanilla",
        str(pkgdown_validator),
        str(pkgdown_config_path),
    )
    pkgdown_assets = {
        "stylesheet": source / "pkgdown/extra.css",
        "script": source / "pkgdown/extra.js",
    }
    missing_pkgdown_assets = [
        label for label, path in pkgdown_assets.items() if not path.is_file()
    ]
    if missing_pkgdown_assets:
        raise AssertionError(
            "missing extracted pkgdown assets: " + ", ".join(missing_pkgdown_assets)
        )
    pkgdown_script = pkgdown_assets["script"].read_text()
    artifact_slug_match = re.search(
        r"const ARTICLE_ARTIFACT_SLUGS = new Set\(\[\s*"
        r'(?P<slugs>(?:"[^"]+"\s*,?\s*)*)\]\);',
        pkgdown_script,
    )
    if artifact_slug_match is None:
        raise AssertionError(
            "pkgdown artifact links must define ARTICLE_ARTIFACT_SLUGS"
        )
    artifact_slugs = set(
        re.findall(r'"([^\"]+)"', artifact_slug_match.group("slugs"))
    )
    notebook_slugs = {
        path.stem for path in (source / "vignettes/colab").glob("*.ipynb")
    }
    if artifact_slugs != notebook_slugs:
        raise AssertionError(
            "pkgdown artifact slug set must match committed Colab notebooks; "
            f"expected {sorted(notebook_slugs)}, got {sorted(artifact_slugs)}"
        )
    if "function getArticleArtifactSlug(pathname) {" not in pkgdown_script:
        raise AssertionError(
            "pkgdown artifact links must resolve eligibility from the pathname"
        )
    if (
        "const slug = getArticleArtifactSlug(window.location.pathname);"
        not in pkgdown_script
    ):
        raise AssertionError(
            "pkgdown artifact links must use pathname-based eligibility"
        )
    for selector in ("small.dont-index code", ".name code"):
        if selector in pkgdown_script:
            raise AssertionError(
                "pkgdown artifact links must not read internal source-label markup: "
                + selector
            )
    goatcounter_header = (
        "  includes:\n"
        "    in_header: |\n"
        "      <script data-goatcounter=\"https://bifrost.goatcounter.com/count\"\n"
        "              async src=\"https://gc.zgo.at/count.js\"></script>\n"
    )
    if goatcounter_header not in pkgdown_config:
        raise AssertionError(
            "pkgdown config must inline the GoatCounter header include"
        )
    if "in_header: pkgdown/extra-head.html" in pkgdown_config:
        raise AssertionError(
            "pkgdown config must not reference pkgdown/extra-head.html"
        )
    if (source / "pkgdown/extra-head.html").exists():
        raise AssertionError("pkgdown/extra-head.html must be removed")
    mermaid_import = (
        "import('https://cdn.jsdelivr.net/npm/mermaid@10.9.1/dist/"
        "mermaid.esm.min.mjs?v=1091')"
    )
    if mermaid_import not in pkgdown_script:
        raise AssertionError("pkgdown script must dynamically import pinned Mermaid 10.9.1")
    if "window.mermaid = mermaid;" not in pkgdown_script:
        raise AssertionError("pkgdown script must expose Mermaid on window")
    if "const encodedSlug = encodeURIComponent(slug);" not in pkgdown_script:
        raise AssertionError("pkgdown artifact links must URL-encode vignette slugs")
    if "pdf.href = './' + encodedSlug + '.pdf';" not in pkgdown_script:
        raise AssertionError("pkgdown PDF link must use the encoded vignette slug")
    if "encodedSlug + '.ipynb';" not in pkgdown_script:
        raise AssertionError("pkgdown Colab link must use the encoded vignette slug")
    if "actions.setAttribute('role', 'group');" not in pkgdown_script:
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
    manifest_path = source / "data-remote/empirical-artifacts.json"
    manifest = json.loads(manifest_path.read_text())
    stale_location_sources = [
        artifact["path"]
        for artifact in manifest["artifacts"]
        if re.search(
            r"(?i)\bpackage-local\b",
            artifact.get("transformation", {}).get("method", ""),
        )
    ]
    simulation_producer = (
        source / "data-raw/run_simulation_study_vignette_grids.R"
    ).read_text()
    if re.search(r"(?i)\bpackaged passerine\b", simulation_producer):
        stale_location_sources.append(
            "data-raw/run_simulation_study_vignette_grids.R"
        )
    if stale_location_sources:
        raise AssertionError(
            "authoritative sources retain stale package-local passerine wording: "
            + ", ".join(stale_location_sources)
        )
    expected_ids = {
        "jaw-tree", "jaw-landmarks", "passerine-tree", "passerine-traits",
        "passerine-search", "passerine-sensitivity", "passerine-posthoc",
        "simulation-preview-tables",
    }

    source_vignette_paths = [
        source / "vignettes/jaw-shape-vignette.Rmd",
        source / "vignettes/rate-map-jaw-shape-vignette.Rmd",
        source / "vignettes/rate-map-jaw-shape-part-2-comparisons.Rmd",
        source / "vignettes/avian-skeleton-part-1.Rmd",
        source / "vignettes/avian-skeleton-part-2.Rmd",
        source / "vignettes/avian-skeleton-part-3.Rmd",
        source / "vignettes/avian-skeleton-part-4.Rmd",
        source / "vignettes/avian-skeleton-part-5.Rmd",
        source / "vignettes/simulation-study-part-1.Rmd",
        source / "vignettes/simulation-study-part-2.Rmd",
    ]
    source_vignettes = {
        path.name: path.read_text() for path in source_vignette_paths
    }
    forbidden_source_patterns = {
        "system.file()": r"\bsystem\.file\s*\(",
        "pkg_file helper": r"\bpkg_file\s*<-\s*function\b",
        "inst/extdata path": r"\binst/extdata\b",
        "stale packaged-data claim": (
            r"(?i)\b(?:bundled|packaged|package-local)\s+"
            r"(?:empirical\s+)?(?:data|dataset|artifact|file|copy|inventory)\b"
        ),
    }
    for filename, text in source_vignettes.items():
        forbidden = [
            label
            for label, pattern in forbidden_source_patterns.items()
            if re.search(pattern, text)
        ]
        if forbidden:
            raise AssertionError(
                f"{filename} retains package-local empirical-data loading: "
                + ", ".join(forbidden)
            )
        if "bifrost_example_file(" not in text:
            raise AssertionError(
                f"{filename} must resolve empirical data with bifrost_example_file()"
            )
        if re.search(r"refresh\s*=\s*TRUE", text) is None:
            raise AssertionError(
                f"{filename} must show the explicit refresh = TRUE update check"
            )

    resolved_ids = set()
    for text in source_vignettes.values():
        resolved_ids.update(
            re.findall(r'bifrost_example_file\(\s*"([^"]+)"', text)
        )
    missing_source_ids = expected_ids - resolved_ids
    if missing_source_ids:
        raise AssertionError(
            "source vignettes do not resolve every empirical artifact identifier: "
            + ", ".join(sorted(missing_source_ids))
        )
    if (source / "inst/extdata").exists():
        raise AssertionError("ordinary empirical data must live outside inst/extdata")
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
    if f"      - {browser_smoke_path}" not in pr_workflow:
        raise AssertionError(
            "PR artifact workflow must watch the pkgdown browser smoke suite"
        )
    for watched_path in (
        "      - pkgdown/**",
        "      - tools/vignette_artifact_checks/**",
        "      - tools/validate-pkgdown-config.R",
    ):
        if watched_path not in pr_workflow:
            raise AssertionError(
                f"PR artifact workflow must watch {watched_path.strip()[2:]}"
            )
    if "      - tools/validate-empirical-artifacts.py" not in pr_workflow:
        raise AssertionError("PR artifact workflow must watch the artifact validator")
    if "      - tools/avian-skeleton/**" not in pr_workflow:
        raise AssertionError("PR artifact workflow must watch artifact generators")
    artifact_schedule = (
        "  vignette-artifacts:\n"
        "    if: >-\n"
        "      github.actor != 'github-actions[bot]' ||\n"
        "      github.event.pull_request.user.login == 'github-actions[bot]'\n"
        "    concurrency:\n"
        "      group: vignette-artifacts-${{ github.event.pull_request.number }}\n"
        "      cancel-in-progress: true\n"
    )
    if artifact_schedule not in pr_workflow:
        raise AssertionError(
            "PR artifact job must skip bot notebook follow-ups on human-authored "
            "PRs, allow bot-authored maintenance PRs, and cancel stale runs per PR"
        )
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
