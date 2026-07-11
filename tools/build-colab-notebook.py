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


REPO = "https://github.com/jakeberv/bifrost"
COLAB_BASE = "https://colab.research.google.com/github/jakeberv/bifrost/blob/main/vignettes/colab"
COLAB_BADGE = "https://colab.research.google.com/assets/colab-badge.svg"
RAW_BASE = "https://raw.githubusercontent.com/jakeberv/bifrost/main/vignettes"
PKGDOWN_ARTICLES = "https://jakeberv.com/bifrost/articles"


def find_repo_root(start: Path) -> Path:
    path = start.resolve()
    for candidate in [path, *path.parents]:
        if (candidate / "DESCRIPTION").exists() and (candidate / "_pkgdown.yml").exists():
            return candidate
    raise SystemExit(f"Could not find repository root from {start}")


def discover_slugs(repo_root: Path) -> list[str]:
    rmds = sorted((repo_root / "vignettes").glob("*.Rmd"))
    return [path.stem for path in rmds]


def markdown_cell(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source}


def code_cell(source: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": source}


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


def is_eval_true(header: str) -> bool:
    return header_has(header, "eval", "TRUE")


def sets_default_eval_false(code: str) -> bool:
    return (
        "opts_chunk$set" in code
        and re.search(r"\beval\s*=\s*FALSE\b", code, re.IGNORECASE) is not None
    )


def is_html_or_pkgdown(header: str, code: str) -> bool:
    lowered = code.lower()
    if re.search(r"\bresults\s*=\s*['\"]asis['\"]", header, re.IGNORECASE):
        return True
    return any(
        marker in lowered
        for marker in [
            "<style",
            "<script",
            "knitr::asis_output",
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


def setup_source() -> str:
    return """if (!dir.exists("/content/bifrost")) {
  system("git clone --depth 1 https://github.com/jakeberv/bifrost.git /content/bifrost")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_local("/content/bifrost", dependencies = TRUE, upgrade = "never")
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
            f"Converted from [`vignettes/{slug}.Rmd`]({REPO}/blob/main/vignettes/{slug}.Rmd)."
        ),
        code_cell(setup_source()),
    ]

    default_eval_false = False
    for block_type, header, content in parse_chunks(body):
        if block_type == "markdown":
            if content:
                cells.append(markdown_cell(rewrite_markdown_links(content) + "\n"))
            continue

        image_paths = image_paths_from_code(content)
        effective_eval_false = is_eval_false(header) or (default_eval_false and not is_eval_true(header))

        if image_paths and not effective_eval_false:
            cells.append(markdown_cell(raw_image_markdown(image_paths) + "\n"))
            if sets_default_eval_false(content):
                default_eval_false = True
            continue

        if effective_eval_false:
            if content.strip():
                cells.append(markdown_cell("```r\n" + content.strip() + "\n```\n"))
            if sets_default_eval_false(content):
                default_eval_false = True
            continue

        if is_html_or_pkgdown(header, content):
            cells.append(markdown_cell("_This pkgdown-only HTML chunk was omitted from the Colab notebook._\n"))
            if sets_default_eval_false(content):
                default_eval_false = True
            continue

        if content.strip():
            cells.append(code_cell(content.strip() + "\n"))

        if sets_default_eval_false(content):
            default_eval_false = True

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
