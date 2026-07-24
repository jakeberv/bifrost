"""Repository, workflow, package, coverage, and provenance contracts."""

from __future__ import annotations

import re
from pathlib import Path

from .common import run

__all__ = ["run_repository_contract_checks"]


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
