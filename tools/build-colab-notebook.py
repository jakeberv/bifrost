#!/usr/bin/env python3
"""Convert bifrost R Markdown vignettes into lightweight Colab notebooks."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.parse
from pathlib import Path

from colab_dependencies import (
    BASE_R_PACKAGES,
    COMMON_COLAB_PACKAGES,
    description_hard_dependencies,
    referenced_r_packages,
)


REPO = "https://github.com/jakeberv/bifrost"
COLAB_BASE = "https://colab.research.google.com/github/jakeberv/bifrost/blob/main/vignettes/colab"
COLAB_BADGE = "https://colab.research.google.com/assets/colab-badge.svg"
RAW_BASE = "https://raw.githubusercontent.com/jakeberv/bifrost/main/vignettes"
PKGDOWN_ARTICLES = "https://jakeberv.com/bifrost/articles"
RUNTIME_NOTE = (
    "> **Recommended Colab runtime:** For compute-intensive examples, choose "
    "**Runtime > Change runtime type > v5e-1 TPU** when available. This notebook "
    "uses the runtime's host CPUs, not the TPU itself. Otherwise, choose the "
    "available runtime with the most CPUs."
)
INLINE_R = re.compile(r"`r\s+(.+?)`")


def find_repo_root(start: Path) -> Path:
    path = start.resolve()
    for candidate in [path, *path.parents]:
        if (candidate / "DESCRIPTION").exists() and (candidate / "_pkgdown.yml").exists():
            return candidate
    raise SystemExit(f"Could not find repository root from {start}")


def discover_slugs(repo_root: Path) -> list[str]:
    rmds = sorted((repo_root / "vignettes").glob("*.Rmd"))
    return [path.stem for path in rmds]


def declared_colab_packages(value: str) -> set[str]:
    packages = {package.strip() for package in value.split(",") if package.strip()}
    invalid = sorted(
        package
        for package in packages
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9.]*", package) is None
    )
    if invalid:
        raise SystemExit(
            "Invalid package name in colab-packages: " + ", ".join(invalid)
        )
    return packages


def colab_packages(
    repo_root: Path,
    cells: list[dict],
    declared: set[str] | None = None,
) -> tuple[str, ...]:
    code = "\n".join(
        cell["source"] for cell in cells if cell["cell_type"] == "code"
    )
    referenced = referenced_r_packages(repo_root, code)
    available = (
        BASE_R_PACKAGES
        | description_hard_dependencies(repo_root / "DESCRIPTION")
        | set(COMMON_COLAB_PACKAGES)
        | {"bifrost"}
    )
    optional = sorted((referenced | (declared or set())) - available)
    return COMMON_COLAB_PACKAGES + tuple(optional)


def markdown_cell(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source}


def code_cell(source: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": source}


def inline_r_code(markdown: str) -> str:
    arguments: list[str] = []
    cursor = 0
    for match in INLINE_R.finditer(markdown):
        if match.start() > cursor:
            arguments.append(json.dumps(markdown[cursor : match.start()]))
        arguments.append(f"as.character({match.group(1).strip()})")
        cursor = match.end()
    if cursor < len(markdown):
        arguments.append(json.dumps(markdown[cursor:]))
    return "cat(paste0(\n  " + ",\n  ".join(arguments) + "\n))\n"


def append_markdown_cells(cells: list[dict], source: str) -> None:
    markdown: list[str] = []

    def flush_markdown() -> None:
        if markdown and "".join(markdown).strip():
            cells.append(markdown_cell("".join(markdown)))
        markdown.clear()

    for line in source.splitlines(keepends=True):
        if INLINE_R.search(line):
            flush_markdown()
            cells.append(code_cell(inline_r_code(line)))
        else:
            markdown.append(line)
    flush_markdown()


def colab_url(slug: str) -> str:
    return f"{COLAB_BASE}/{urllib.parse.quote(slug)}.ipynb"


def colab_badge_markdown(slug: str) -> str:
    return f"[![Open In Colab]({COLAB_BADGE})]({colab_url(slug)})"


def parse_front_matter(lines: list[str]) -> tuple[dict[str, str], list[str]]:
    if not lines or lines[0].strip() != "---":
        return {}, lines
    end = None
    for idx in range(1, len(lines)):
        if lines[idx].strip() == "---":
            end = idx
            break
    if end is None:
        return {}, lines
    yaml_lines = lines[1:end]
    meta: dict[str, str] = {}
    for line in yaml_lines:
        match = re.match(r"^([A-Za-z0-9_.-]+):\s*(.*)$", line)
        if match:
            key, value = match.groups()
            value = value.strip().strip("\"'")
            if value:
                meta[key] = value
    return meta, lines[end + 1 :]


def parse_chunks(lines: list[str]) -> list[tuple[str, str, str]]:
    blocks: list[tuple[str, str, str]] = []
    markdown: list[str] = []
    idx = 0
    chunk_start = re.compile(r"^```\{r(.*)\}\s*$")
    while idx < len(lines):
        match = chunk_start.match(lines[idx])
        if not match:
            markdown.append(lines[idx])
            idx += 1
            continue
        if markdown:
            blocks.append(("markdown", "", "".join(markdown).strip()))
            markdown = []
        header = match.group(1).strip()
        idx += 1
        code: list[str] = []
        while idx < len(lines) and not re.match(r"^```\s*$", lines[idx]):
            code.append(lines[idx])
            idx += 1
        if idx < len(lines):
            idx += 1
        blocks.append(("r", header, "".join(code).rstrip()))
    if markdown:
        blocks.append(("markdown", "", "".join(markdown).strip()))
    return blocks


def header_has(header: str, option: str, value: str | None = None) -> bool:
    if value is None:
        return re.search(rf"\b{re.escape(option)}\b", header, re.IGNORECASE) is not None
    return re.search(
        rf"\b{re.escape(option)}\s*=\s*['\"]?{re.escape(value)}['\"]?",
        header,
        re.IGNORECASE,
    ) is not None


def is_eval_false(header: str) -> bool:
    return header_has(header, "eval", "FALSE")


def is_include_false(header: str) -> bool:
    return header_has(header, "include", "FALSE")


def is_hidden_unevaluated(header: str) -> bool:
    return is_include_false(header) and is_eval_false(header)


def is_html_or_pkgdown(header: str, code: str) -> bool:
    lowered = code.lower()
    return any(
        marker in lowered
        for marker in [
            "<style",
            "<script",
            "knitr::is_html_output",
            "sys.getenv(\"in_pkgdown\")",
            "sys.getenv('in_pkgdown')",
            "htmltools::",
        ]
    )


def image_paths_from_code(code: str) -> list[str]:
    paths: list[str] = []
    for pattern in [
        r"(?:knitr::)?include_graphics\(\s*['\"]([^'\"]+)['\"]",
        r"(?:knitr::)?image_uri\(\s*['\"]([^'\"]+)['\"]",
    ]:
        paths.extend(re.findall(pattern, code))
    paths.extend(
        f"rate-map-jaw-shape/{match}"
        for match in re.findall(
            r"(?:knitr::)?include_graphics\(\s*rate_map_artifact\(\s*['\"]([^'\"]+)['\"]",
            code,
        )
    )
    return paths


def raw_image_markdown(paths: list[str]) -> str:
    lines = []
    for path in paths:
        normalized = path.lstrip("./")
        quoted = urllib.parse.quote(normalized)
        label = os.path.basename(normalized)
        lines.append(f"![{label}]({RAW_BASE}/{quoted})")
    return "\n\n".join(lines)


def rewrite_markdown_links(markdown: str) -> str:
    pattern = re.compile(
        r"(\]\()((?!https?://|mailto:|#)[^\s)]+\.html(?:#[^\s)]*)?)(\))",
        re.IGNORECASE,
    )
    return pattern.sub(
        lambda match: f"{match.group(1)}{PKGDOWN_ARTICLES}/{match.group(2)}{match.group(3)}",
        markdown,
    )


def setup_source(packages: tuple[str, ...]) -> str:
    package_vector = ", ".join(json.dumps(package) for package in packages)
    return f"""message("Detected logical CPUs: ", parallel::detectCores(logical = TRUE))
if (!dir.exists("/content/bifrost")) {{
  system("git clone --depth 1 https://github.com/jakeberv/bifrost.git /content/bifrost")
}}
colab_packages <- c({package_vector})
missing_packages <- colab_packages[
  !vapply(colab_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {{
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}}
remotes::install_local("/content/bifrost", dependencies = NA, upgrade = "never")
setwd("/content/bifrost/vignettes")
"""


def convert(slug: str, repo_root: Path) -> dict:
    rmd_path = repo_root / "vignettes" / f"{slug}.Rmd"
    if not rmd_path.exists():
        raise SystemExit(f"Missing vignette source: {rmd_path}")
    lines = rmd_path.read_text(encoding="utf-8").splitlines(keepends=True)
    meta, body = parse_front_matter(lines)
    title = meta.get("title", slug)

    cells: list[dict] = [
        markdown_cell(
            f"# {title}\n\n"
            f"{colab_badge_markdown(slug)}\n\n"
            f"Converted from [`vignettes/{slug}.Rmd`]({REPO}/blob/main/vignettes/{slug}.Rmd).\n\n"
            f"{RUNTIME_NOTE}"
        ),
    ]

    for block_type, header, content in parse_chunks(body):
        if block_type == "markdown":
            if content:
                append_markdown_cells(
                    cells,
                    rewrite_markdown_links(content) + "\n",
                )
            continue

        if is_hidden_unevaluated(header):
            continue

        # Evaluated include=FALSE chunks often define setup objects used by
        # later visible cells. Preserve their code even when they also contain
        # output-format helpers; only genuinely unevaluated maintenance chunks
        # are omitted above.
        if not is_include_false(header):
            image_paths = image_paths_from_code(content)
            if image_paths:
                cells.append(markdown_cell(raw_image_markdown(image_paths) + "\n"))
                continue

            if is_html_or_pkgdown(header, content):
                cells.append(markdown_cell("_This pkgdown-only HTML chunk was omitted from the Colab notebook._\n"))
                continue

        if content.strip():
            cells.append(code_cell(content.strip() + "\n"))

    declared = declared_colab_packages(meta.get("colab-packages", ""))
    cells.insert(
        1,
        code_cell(setup_source(colab_packages(repo_root, cells, declared))),
    )

    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "R",
                "language": "R",
                "name": "ir",
            },
            "language_info": {"name": "R"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--slug", action="append", dest="slugs", help="Vignette slug to convert")
    parser.add_argument("--all", action="store_true", help="Convert all top-level vignette Rmd files")
    parser.add_argument("--output-dir", default="vignettes/colab")
    args = parser.parse_args(argv)

    repo_root = find_repo_root(Path.cwd())
    slugs = args.slugs or ([] if not args.all else discover_slugs(repo_root))
    if not slugs:
        parser.error("Provide --slug at least once, or use --all")

    output_dir = repo_root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    for slug in slugs:
        notebook = convert(slug, repo_root)
        out = output_dir / f"{slug}.ipynb"
        out.write_text(json.dumps(notebook, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        try:
            display_path = out.relative_to(repo_root)
        except ValueError:
            display_path = out
        print(f"Wrote {display_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
