"""Detect R package dependencies used by executable Colab cells."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path


BASE_R_PACKAGES = frozenset(
    {
        "base",
        "compiler",
        "datasets",
        "graphics",
        "grDevices",
        "grid",
        "methods",
        "parallel",
        "splines",
        "stats",
        "stats4",
        "tcltk",
        "tools",
        "utils",
    }
)
COMMON_COLAB_PACKAGES = ("remotes", "knitr")

R_PACKAGE_REFERENCE_SCRIPT = r"""
code <- paste(readLines(file("stdin"), warn = FALSE), collapse = "\n")
expressions <- parse(text = code)
packages <- character()

walk <- function(node) {
  if (is.call(node)) {
    operator <- node[[1L]]
    if (is.symbol(operator)) {
      operator_name <- as.character(operator)
      if (operator_name %in% c("::", ":::") && length(node) >= 3L) {
        packages <<- c(packages, as.character(node[[2L]]))
      }
      if (operator_name %in% c("library", "require") &&
          length(node) >= 2L &&
          (is.symbol(node[[2L]]) || is.character(node[[2L]]))) {
        packages <<- c(packages, as.character(node[[2L]]))
      }
      if (operator_name %in% c("requireNamespace", "loadNamespace") &&
          length(node) >= 2L && is.character(node[[2L]])) {
        packages <<- c(packages, as.character(node[[2L]]))
      }
    }
    invisible(lapply(as.list(node), walk))
  } else if (is.expression(node) || is.pairlist(node)) {
    invisible(lapply(as.list(node), walk))
  }
}

walk(expressions)
cat(sort(unique(packages)), sep = "\n")
"""


def description_hard_dependencies(path: Path) -> set[str]:
    fields: dict[str, list[str]] = {}
    current = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:1].isspace() and current:
            fields[current].append(line.strip())
            continue
        if ":" not in line:
            current = ""
            continue
        current, value = line.split(":", 1)
        fields[current] = [value.strip()]

    packages: set[str] = set()
    for field in ("Depends", "Imports", "LinkingTo"):
        for entry in ",".join(fields.get(field, [])).split(","):
            match = re.match(r"\s*([A-Za-z][A-Za-z0-9.]*)", entry)
            if match and match.group(1) != "R":
                packages.add(match.group(1))
    return packages


def referenced_r_packages(repo_root: Path, code: str) -> set[str]:
    try:
        result = subprocess.run(
            ["Rscript", "--vanilla", "-e", R_PACKAGE_REFERENCE_SCRIPT],
            cwd=repo_root,
            input=code,
            text=True,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError as error:
        raise SystemExit("Rscript is required to detect Colab dependencies") from error
    if result.returncode != 0:
        raise SystemExit(f"Could not parse vignette R code:\n{result.stderr}")
    return set(result.stdout.split())
