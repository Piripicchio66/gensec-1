#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# setup_licensing.py — Download the AGPL-3.0 license text, optionally
#                       replace the project-name placeholder, and verify
#                       that all required licensing files are present.
#
# Usage:
#     python setup_licensing.py                   # uses default name
#     python setup_licensing.py MyProject         # replaces __PROJECT_NAME__
#     python setup_licensing.py GenSec --header    # also prepends license
#                                                 # headers to all .py files
#
# Works on Windows, macOS, and Linux — requires only the standard library.
# ---------------------------------------------------------------------------
from __future__ import annotations

import os
import sys
import textwrap
import urllib.request
from pathlib import Path

AGPL_URL = "https://www.gnu.org/licenses/agpl-3.0.txt"
PLACEHOLDER = "__PROJECT_NAME__"

REQUIRED_FILES = [
    "LICENSE",
    "NOTICE",
    "CLA.md",
    "CONTRIBUTING.md",
    "COMMERCIAL_LICENSE.md",
    os.path.join(".github", "pull_request_template.md"),
    os.path.join(".github", "workflows", "cla.yml"),
]

LICENSE_HEADER = textwrap.dedent("""\
    # ---------------------------------------------------------------------------
    # {project} — Copyright (c) 2026 Andrea Albero
    #
    # This file is part of {project}.
    #
    # {project} is free software: you can redistribute it and/or modify it under
    # the terms of the GNU Affero General Public License as published by the
    # Free Software Foundation, either version 3 of the License, or (at your
    # option) any later version.
    #
    # {project} is distributed in the hope that it will be useful, but WITHOUT
    # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    # FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public
    # License for more details.
    #
    # You should have received a copy of the GNU Affero General Public License
    # along with {project}.  If not, see <https://www.gnu.org/licenses/>.
    # ---------------------------------------------------------------------------
""")


def download_license(repo_root: Path) -> None:
    """Download the official AGPL-3.0 text into ``LICENSE``."""
    dest = repo_root / "LICENSE"
    print(f"==> Downloading AGPL-3.0 from {AGPL_URL} ...")
    try:
        urllib.request.urlretrieve(AGPL_URL, dest)
        lines = dest.read_text(encoding="utf-8").splitlines()
        print(f"    LICENSE written ({len(lines)} lines).")
    except Exception as exc:
        print(f"    [ERROR] Download failed: {exc}")
        print("    Download manually from:", AGPL_URL)
        print("    and save it as LICENSE in the repo root.")


def replace_placeholder(repo_root: Path, project_name: str) -> None:
    """Replace ``__PROJECT_NAME__`` in all licensing files."""
    targets = ["NOTICE", "CLA.md", "CONTRIBUTING.md", "COMMERCIAL_LICENSE.md"]
    count = 0
    for name in targets:
        path = repo_root / name
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        if PLACEHOLDER in text:
            path.write_text(
                text.replace(PLACEHOLDER, project_name), encoding="utf-8"
            )
            count += 1
            print(f"    {name}: replaced {PLACEHOLDER} -> {project_name}")
    if count == 0:
        print(f"    No {PLACEHOLDER} placeholders found (files may already "
              f"use the final project name).")


def prepend_headers(repo_root: Path, project_name: str) -> None:
    """Prepend the AGPL license header to all ``.py`` files under ``src/``."""
    header = LICENSE_HEADER.format(project=project_name)
    marker = "This file is part of"
    src_dir = repo_root / "src"
    if not src_dir.exists():
        # Fall back: look for .py files directly in repo root.
        src_dir = repo_root

    py_files = sorted(src_dir.rglob("*.py"))
    added, skipped = 0, 0
    for py in py_files:
        content = py.read_text(encoding="utf-8")
        if marker in content:
            skipped += 1
            continue
        py.write_text(header + content, encoding="utf-8")
        added += 1
        print(f"    + {py.relative_to(repo_root)}")

    print(f"    Header added to {added} files, {skipped} already had it.")


def verify(repo_root: Path) -> bool:
    """Check that all required licensing files exist."""
    print("==> Verifying licensing files ...")
    all_ok = True
    for rel in REQUIRED_FILES:
        path = repo_root / rel
        tag = "[OK]" if path.exists() else "[MISSING]"
        if not path.exists():
            all_ok = False
        print(f"    {tag} {rel}")
    return all_ok


def main() -> None:
    repo_root = Path(__file__).resolve().parent
    project_name = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else None
    add_headers = "--header" in sys.argv or "--headers" in sys.argv

    # Step 1 — download LICENSE.
    download_license(repo_root)

    # Step 2 — replace placeholders if a project name was given.
    if project_name:
        print(f"==> Replacing placeholder with project name: {project_name}")
        replace_placeholder(repo_root, project_name)

    # Step 3 — optionally add headers to .py files.
    if add_headers:
        name = project_name or "GenSec"
        print(f"==> Adding license headers to .py files (project: {name}) ...")
        prepend_headers(repo_root, name)

    # Step 4 — verify.
    all_ok = verify(repo_root)

    print()
    if all_ok:
        print("All licensing files in place.")
    else:
        print("[!] Some files are missing — copy them from the template.")

    print()
    print("Next steps:")
    print("  1. Replace [your-email@example.com] in COMMERCIAL_LICENSE.md")
    print("  2. Review NOTICE and update contributor entries if needed")
    print("  3. git add . && git commit -m 'chore(license): AGPL-3.0 dual licensing'")
    print("  4. Push to GitHub — the CLA bot activates automatically on PRs")


if __name__ == "__main__":
    main()