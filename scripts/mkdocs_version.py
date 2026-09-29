"""MkDocs hook: stamp the site with the version it was built from.

The version comes from Cargo.toml (single source; pyproject.toml is dynamic).
Docs are only deployed from the release workflow, so this always matches the
published package versions.
"""

import re
from pathlib import Path

from mkdocs.config.defaults import MkDocsConfig

ROOT = Path(__file__).resolve().parents[1]
CARGO_TOML = ROOT / "Cargo.toml"


def read_version() -> str:
    match = re.search(r'^version = "([^"]+)"', CARGO_TOML.read_text(), re.MULTILINE)
    if match is None:
        raise ValueError(f"No `version = \"...\"` found in {CARGO_TOML}")
    return match.group(1)


def on_config(config: MkDocsConfig) -> MkDocsConfig:
    version = read_version()
    release_url = f"{config.repo_url}/releases/tag/v{version}"
    config.extra["immunum_version"] = version
    config.copyright = (
        f'Documentation for <a href="{release_url}">immunum v{version}</a>'
    )
    return config
