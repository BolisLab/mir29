#!/usr/bin/env python3
"""Read-only publication-safety and source-layout checks; standard library only.

This checks repository packaging, not the scientific models or their results.
No files are edited, installed, downloaded, rendered or executed. Exit status is
nonzero for errors. Warnings flag expected checkpoint boundaries and manual
review needs; they do not assert that redistributing the inputs is authorized.
"""

from __future__ import annotations

import argparse
import csv
from html import unescape
import io
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit


SKIP_DIRS = {".git", ".quarto", "_site", "__pycache__", "node_modules", ".venv"}
FORBIDDEN_ENDINGS = (
    ".rds", ".rdata", ".rda", ".rhistory", ".fastq", ".fastq.gz",
    ".fq", ".fq.gz", ".bam", ".cram", ".sam", ".idat", ".idat.gz",
    ".zip", ".tar", ".tar.gz", ".tgz",
)
TEXT_SUFFIXES = {
    ".md", ".rmd", ".qmd", ".r", ".py", ".sh", ".txt", ".csv", ".tsv",
    ".json", ".yml", ".yaml", ".toml", ".css", ".js", ".html", ".svg", ".cff",
}
MAX_FILE_BYTES = 50 * 1024 * 1024
SECRET_PATTERNS = {
    "GitHub access token": re.compile(r"\bgh[pousr]_[A-Za-z0-9]{30,}\b"),
    "GitHub fine-grained token": re.compile(r"\bgithub_pat_[A-Za-z0-9_]{30,}\b"),
    "AWS access key": re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b"),
    "private key block": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----"),
    "OpenAI key": re.compile(r"\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{32,}\b"),
}
PERSONAL_PATH = re.compile(r"(?:/Users/[^\s/]+/|/home/[^\s/]+/|/g100_scratch/|~/Desktop/)")
GTEX_ID = re.compile(r"\bGTEX-[A-Z0-9]{4,}(?:-[A-Z0-9]+)*\b")
EXACT_AGE_HEADERS = {"age", "age_years", "age_in_years", "chronological_age", "donor_age"}


def visible_markdown(text: str) -> str:
    """Mask fenced code so literal examples are not mistaken for live images."""
    out = []
    fence = None
    for line in text.splitlines():
        match = re.match(r"^\s*(`{3,}|~{3,})", line)
        if match:
            marker = match.group(1)
            if fence is None:
                fence = (marker[0], len(marker))
            elif marker[0] == fence[0] and len(marker) >= fence[1]:
                fence = None
            out.append("")
        else:
            out.append(line if fence is None else "")
    return "\n".join(out)


def image_targets(text: str):
    text = visible_markdown(text)
    # Inline Markdown images, with an optional title, and quoted HTML src.
    for match in re.finditer(r"!\[[^\]]*\]\(\s*(<[^>]+>|[^\s)]+)(?:\s+[^)]*)?\)", text):
        yield match.group(1).strip("<>")
    for match in re.finditer(r"<img\b[^>]*?\bsrc\s*=\s*([\"'])(.*?)\1", text, re.I | re.S):
        yield match.group(2)
    definitions = {
        key.strip().casefold(): target.strip("<>")
        for key, target in re.findall(r"^\s*\[([^\]]+)\]:\s*(<[^>]+>|\S+)", text, re.M)
    }
    for label, key in re.findall(r"!\[([^\]]*)\]\[([^\]]*)\]", text):
        normalized = (key or label).strip().casefold()
        if normalized in definitions:
            yield definitions[normalized]


def inspect_repository(root: Path) -> dict:
    errors = []
    warnings = []
    counts = {"files": 0, "bytes": 0, "text_files": 0, "local_image_references": 0,
              "tables_inspected": 0, "modules": 0}
    for i in range(1, 6):
        name = f"Figure{i}_and_S{i}"
        module = root / "analyses" / name
        required = [module / "READABLE.md", module / f"{name}.Rmd"]
        for path in required:
            if not path.is_file():
                errors.append(f"Missing required module file: {path.relative_to(root)}")
        if all(path.is_file() for path in required):
            counts["modules"] += 1

    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if any(part in SKIP_DIRS for part in relative.parts):
            continue
        if "private_audit" in relative.parts:
            errors.append(f"Private audit material inside repository: {relative}")
            continue
        if path.is_symlink():
            errors.append(f"Symlink requires review before publication: {relative}")
            continue
        if not path.is_file():
            continue
        counts["files"] += 1
        size = path.stat().st_size
        counts["bytes"] += size
        if str(path).lower().endswith(FORBIDDEN_ENDINGS):
            errors.append(f"Disallowed raw/checkpoint/archive file in Git payload: {relative}")
            continue
        if size > MAX_FILE_BYTES:
            errors.append(f"File exceeds 50 MiB source-repository budget: {relative}")
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES and path.name not in {".gitignore", "LICENSE"}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeError, OSError):
            errors.append(f"Unreadable text file: {relative}")
            continue
        counts["text_files"] += 1
        for label, pattern in SECRET_PATTERNS.items():
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"Possible {label}: {relative}:{line} (value not displayed)")
        for match in PERSONAL_PATH.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            if path.resolve() == Path(__file__).resolve() and text.splitlines()[line - 1].startswith("PERSONAL_PATH = "):
                # The detector's own regex contains literal path prefixes.
                continue
            errors.append(f"Personal/server absolute path: {relative}:{line}")
        if GTEX_ID.search(text):
            warnings.append(f"GTEx identifier present; inspect disclosure context: {relative}")
        if path.suffix.lower() in {".md", ".rmd", ".qmd"}:
            for target in image_targets(text):
                target = unescape(target)
                parsed = urlsplit(target)
                if parsed.scheme or target.startswith(("//", "#")):
                    continue
                if any(marker in target for marker in ("`r ", "{{", "${")):
                    warnings.append(f"Dynamic image target not resolved statically: {relative}")
                    continue
                counts["local_image_references"] += 1
                resolved = (path.parent / unquote(parsed.path)).resolve()
                if not resolved.is_relative_to(root):
                    errors.append(f"Image target escapes repository: {relative} -> {target}")
                elif not resolved.is_file():
                    errors.append(f"Unresolved local image: {relative} -> {target}")
        if path.suffix.lower() in {".csv", ".tsv"}:
            counts["tables_inspected"] += 1
            try:
                rows = csv.DictReader(io.StringIO(text), delimiter="\t" if path.suffix.lower() == ".tsv" else ",")
                headers = {(h or "").strip().casefold() for h in (rows.fieldnames or [])}
                if headers & EXACT_AGE_HEADERS:
                    if GTEX_ID.search(text):
                        errors.append(f"Potential individual GTEx age disclosure in table: {relative}")
                    else:
                        warnings.append(f"Age column requires source-specific sharing review: {relative}")
            except csv.Error:
                errors.append(f"Cannot parse delimited table: {relative}")

    catalogue_path = root / "data" / "input-catalogue.json"
    try:
        catalogue = json.loads(catalogue_path.read_text(encoding="utf-8"))
        modules = catalogue.get("modules", [])
        if len(modules) != 5:
            errors.append("Input catalogue must contain five modules.")
        for module in modules:
            module_path = root / module.get("module_path", "")
            notebook = module_path / module.get("notebook", "")
            if not notebook.is_file():
                errors.append(f"Catalogue notebook path does not resolve: {module.get('id', '?')}")
            if not module.get("public_release_ready"):
                warnings.append(f"{module.get('id', '?')}: checkpoint distribution is pending; fresh checkout is not self-contained.")
        warnings.append("Static checks do not certify scientific reproducibility, redistribution rights or absence of all sensitive data.")
    except (OSError, ValueError, TypeError) as exc:
        errors.append(f"Input catalogue cannot be parsed: {type(exc).__name__}")
    return {"status": "PASS" if not errors else "FAIL", "counts": counts,
            "errors": sorted(set(errors)), "warnings": sorted(set(warnings))}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--json", action="store_true", help="Print a machine-readable report to stdout.")
    args = parser.parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        parser.error("Repository root does not exist.")
    report = inspect_repository(root)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"Repository packaging: {report['status']}")
        print("; ".join(f"{key}={value}" for key, value in report["counts"].items()))
        for error in report["errors"]:
            print(f"ERROR: {error}")
        for warning in report["warnings"]:
            print(f"WARNING: {warning}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
