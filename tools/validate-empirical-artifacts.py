#!/usr/bin/env python3
"""Validate or refresh checksums in the empirical-artifact manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path, PurePosixPath


MANIFEST = Path("inst/extdata/empirical-artifacts.json")
SHA256 = re.compile(r"^[0-9a-f]{64}$")


def find_repo_root(start: Path) -> Path:
    for candidate in [start.resolve(), *start.resolve().parents]:
        if (candidate / "DESCRIPTION").exists() and (candidate / MANIFEST).exists():
            return candidate
    raise SystemExit(f"Could not find package root and {MANIFEST} from {start}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_text(mapping: dict, key: str, context: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise AssertionError(f"{context} requires non-empty {key!r}")
    return value


def scoped_files(root: Path, scope: dict) -> set[str]:
    directories = scope.get("directories")
    extensions = scope.get("extensions")
    if not isinstance(directories, list) or not directories:
        raise AssertionError("manifest scope requires a non-empty directories list")
    if not isinstance(extensions, list) or not extensions:
        raise AssertionError("manifest scope requires a non-empty extensions list")
    normalized_extensions = {str(extension).lower() for extension in extensions}
    files: set[str] = set()
    for directory in directories:
        relative = PurePosixPath(require_text({"directory": directory}, "directory", "scope"))
        if relative.is_absolute() or ".." in relative.parts:
            raise AssertionError(f"unsafe scope directory: {directory}")
        full = root.joinpath(*relative.parts)
        if not full.is_dir():
            raise AssertionError(f"scope directory does not exist: {directory}")
        for path in full.rglob("*"):
            if path.is_file() and path.suffix.lower() in normalized_extensions:
                files.add(path.relative_to(root).as_posix())
    return files


def validate(root: Path, manifest: dict, update_checksums: bool) -> int:
    if manifest.get("schema_version") != 1:
        raise AssertionError("manifest schema_version must be 1")
    require_text(manifest, "description", "manifest")

    sources = manifest.get("sources")
    licenses = manifest.get("license_records")
    artifacts = manifest.get("artifacts")
    if not isinstance(sources, dict) or not sources:
        raise AssertionError("manifest requires source records")
    if not isinstance(licenses, dict) or not licenses:
        raise AssertionError("manifest requires license records")
    if not isinstance(artifacts, list) or not artifacts:
        raise AssertionError("manifest requires artifact records")

    for source_id, source in sources.items():
        for key in ("title", "doi", "url", "version"):
            require_text(source, key, f"source {source_id!r}")
    for license_id, license_record in licenses.items():
        for key in (
            "repository_license",
            "upstream_status",
            "evidence_url",
            "reuse_action",
        ):
            require_text(license_record, key, f"license {license_id!r}")

    paths: list[str] = []
    for index, artifact in enumerate(artifacts, start=1):
        context = f"artifact #{index}"
        path_text = require_text(artifact, "path", context)
        relative = PurePosixPath(path_text)
        if relative.is_absolute() or ".." in relative.parts:
            raise AssertionError(f"{context} has unsafe path: {path_text}")
        path = root.joinpath(*relative.parts)
        if not path.is_file():
            raise AssertionError(f"manifest artifact does not exist: {path_text}")
        paths.append(path_text)

        source_id = require_text(artifact, "source_id", context)
        if source_id not in sources:
            raise AssertionError(f"{context} references unknown source_id {source_id!r}")
        source_location = require_text(artifact, "source_location", context)
        license_id = require_text(artifact, "license_id", context)
        if license_id not in licenses:
            raise AssertionError(f"{context} references unknown license_id {license_id!r}")

        transformation = artifact.get("transformation")
        if not isinstance(transformation, dict):
            raise AssertionError(f"{context} requires a transformation record")
        method = require_text(transformation, "method", f"{context} transformation")
        require_text(transformation, "script", f"{context} transformation")
        if path_text == (
            "inst/extdata/simulation-study-cache/passerine_preview_tables.rds"
        ):
            if "simulation-study-vignette.Rmd" in source_location:
                raise AssertionError(
                    f"{context} references the deleted simulation vignette slug"
                )
            if (
                "schema-3" not in source_location.lower()
                or "schema-3" not in method.lower()
            ):
                raise AssertionError(
                    f"{context} must describe the simulation vignette cache as schema-3"
                )
            if "Extract eight named HTML preview tables" in method:
                raise AssertionError(
                    f"{context} retains the obsolete preview-table description"
                )

        actual = sha256(path)
        recorded = artifact.get("sha256")
        if update_checksums:
            artifact["sha256"] = actual
        elif not isinstance(recorded, str) or SHA256.fullmatch(recorded) is None:
            raise AssertionError(f"{context} requires a lowercase SHA-256 checksum")
        elif actual != recorded:
            raise AssertionError(
                f"checksum mismatch for {path_text}: expected {recorded}, got {actual}"
            )

    if len(paths) != len(set(paths)):
        raise AssertionError("manifest contains duplicate artifact paths")
    expected = scoped_files(root, manifest.get("scope", {}))
    recorded_paths = set(paths)
    missing = sorted(expected - recorded_paths)
    extra = sorted(recorded_paths - expected)
    if missing or extra:
        details = []
        if missing:
            details.append("unrecorded scoped artifacts: " + ", ".join(missing))
        if extra:
            details.append("records outside scope: " + ", ".join(extra))
        raise AssertionError("; ".join(details))
    return len(paths)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--update-checksums",
        action="store_true",
        help="rewrite recorded checksums after intentional artifact changes",
    )
    args = parser.parse_args()

    root = find_repo_root(Path.cwd())
    manifest_path = root / MANIFEST
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    count = validate(root, manifest, args.update_checksums)
    if args.update_checksums:
        manifest_path.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        validate(root, manifest, False)
        print(f"Updated and validated {count} empirical artifact checksums.")
    else:
        print(f"Validated {count} empirical artifact checksums and provenance records.")


if __name__ == "__main__":
    main()
