#!/usr/bin/env python3
"""Build a local documentation copy without running any scientific analysis.

Requires Python 3.9+. Generated files live only in docs/figures and docs/analyses.
The source notebooks, recorded results and manuscript files are never modified.
"""
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
CHAPTERS = [
    (1, "Ageing and rejuvenation", "Transcriptional convergence across physiological ageing, somatic reprogramming and early development."),
    (2, "An evolutionary perspective", "Species, seeds, miRNA loci and conservation of a regulatory CpG."),
    (3, "Transcriptional responses", "miR-29 perturbation, transcriptional-age readouts and independent ageing programmes."),
    (4, "Methylome remodelling", "Dynamics, clocks and the 123-CpG fibroblast reference."),
    (5, "The relationship in cancer", "Tumour and cancer-cell-line analyses of miR-29, molecular age and methylation divergence."),
]


def transform_body(text, module):
    """Keep recorded code/output; make source-relative links work in the site."""
    text = re.sub(r"\A---\s*\n.*?\n---\s*\n", "", text, count=1, flags=re.S)
    def rewrite(match):
        prefix, target, suffix = match.groups()
        if target.startswith(("#", "http:", "https:", "mailto:", "data:", "../")):
            return match.group(0)
        return prefix + "../analyses/" + module + "/" + target + suffix
    text = re.sub(r"(!?\[[^\]\n]*\]\()([^\s)]+)(\))", rewrite, text)
    text = re.sub(r'(src=["\'])(?!https?:|data:|/)([^"\']+)(["\'])', rewrite, text)
    # These are static notebook outputs. Disable Quarto's executable-chunk
    # interpretation, while retaining the language for syntax highlighting.
    text = re.sub(r"^```\{r[^\n]*\}\s*$", "```r", text, flags=re.M)
    text = re.sub(r"^```\{python[^\n]*\}\s*$", "```python", text, flags=re.M)
    # Fold recorded source without executing it; plain output blocks stay open.
    lines, result, in_code, fold = text.splitlines(), [], False, False
    for line in lines:
        if line.startswith("```"):
            if not in_code:
                language = line[3:].strip()
                fold = language in {"r", "R", "python", "bash", "sh"} or language.startswith("{.r")
                if fold:
                    result.extend(["<details class=\"recorded-code\">", "<summary>Inspect recorded code</summary>", ""])
                in_code = True
            else:
                result.append(line)
                if fold:
                    result.extend(["", "</details>", ""])
                in_code, fold = False, False
                continue
        result.append(line)
    text = "\n".join(result) + "\n"
    return text


def main():
    figures = ROOT / "docs" / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for number, title, description in CHAPTERS:
        module = f"Figure{number}_and_S{number}"
        source = ROOT / "analyses" / module
        readable = source / "READABLE.md"
        if not readable.is_file():
            raise SystemExit(f"Missing required chapter: {readable}")
        destination = ROOT / "docs" / "analyses" / module
        shutil.copytree(source, destination, dirs_exist_ok=True)
        # Module instructions are served beneath docs/analyses in the site.
        # Repository-level data documentation remains in GitHub, not docs/data.
        module_readme = destination / "README.md"
        module_readme.write_text(
            module_readme.read_text(encoding="utf-8").replace(
                "](../../data/",
                "](https://github.com/BolisLab/mir29/blob/main/data/",
            ),
            encoding="utf-8",
        )
        header = (
            f'---\ntitle: "{title}"\n'
            f'subtitle: "Figure {number} + Supplementary Figure S{number}"\n'
            f'description: "{description}"\n'
            'toc: true\ntoc-depth: 2\nexecute:\n  eval: false\n---\n\n'
            f'[Module files and input requirements](../analyses/{module}/README.md)'
            f' · [R Markdown source](../analyses/{module}/{module}.Rmd)\n\n'
        )
        (figures / f"figure{number}.qmd").write_text(
            header + transform_body(readable.read_text(), module), encoding="utf-8"
        )
        print(f"Prepared Figure {number} + S{number} (recorded output; no analysis executed)")


if __name__ == "__main__":
    main()
