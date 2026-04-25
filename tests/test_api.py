"""
Tests for the public GenSec API facade.

These tests cover the *contract* of ``gensec.api``: input validation,
caching behaviour, and Pydantic-model shape.  They use a monkeypatched
``_Session`` so they run without the real solver being implemented.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from gensec import api


EXAMPLES = Path(__file__).parent.parent / "examples"


# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------

def test_load_yaml_text_requires_exactly_one_source():
    with pytest.raises(ValueError):
        api._load_yaml_text()
    with pytest.raises(ValueError):
        api._load_yaml_text(yaml_text="a: 1", yaml_path="x.yaml")


def test_load_yaml_text_from_path(tmp_path: Path):
    p = tmp_path / "t.yaml"
    p.write_text("foo: bar\n", encoding="utf-8")
    assert api._load_yaml_text(yaml_path=p) == "foo: bar\n"


def test_load_yaml_text_path_missing_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        api._load_yaml_text(yaml_path=tmp_path / "nope.yaml")


# ---------------------------------------------------------------------------
# Hashing and normalisation
# ---------------------------------------------------------------------------

def test_yaml_key_is_stable_across_formatting():
    """Whitespace/comment/key-order changes must not change the hash."""
    a = "a: 1\nb: 2\n"
    b = "# comment\nb: 2\na: 1\n"
    assert api.yaml_key(a) == api.yaml_key(b)


def test_yaml_key_changes_on_value_change():
    assert api.yaml_key("a: 1") != api.yaml_key("a: 2")


# ---------------------------------------------------------------------------
# Caching
# ---------------------------------------------------------------------------

def test_session_cache_hits_on_identical_input(monkeypatch):
    """Second call with the same YAML must reuse the cached Session."""
    api.clear_cache()
    build_calls = []

    def fake_build(normalised_yaml):
        build_calls.append(normalised_yaml)
        return MagicMock(spec=api._Session)

    monkeypatch.setattr(api._Session, "build", classmethod(
        lambda cls, norm: fake_build(norm)
    ))

    yaml_text = "section:\n  B: 300\n  H: 500\n"
    api._get_session(api.yaml_key(yaml_text), api._normalise_yaml(yaml_text))
    api._get_session(api.yaml_key(yaml_text), api._normalise_yaml(yaml_text))
    assert len(build_calls) == 1


def test_clear_cache_evicts(monkeypatch):
    api.clear_cache()
    info = api._get_session.cache_info()
    assert info.currsize == 0


# ---------------------------------------------------------------------------
# Example YAML smoke tests  (skipped until _Session.build is wired)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not (EXAMPLES / "biaxial_column.yaml").exists(),
    reason="examples directory not present",
)

@pytest.mark.timeout(60)
def test_analyze_returns_valid_model(tmp_path):
    """Smoke test on a low-resolution biaxial column."""
    src = EXAMPLES / "biaxial_column.yaml"
    if not src.exists():
        pytest.skip("examples/biaxial_column.yaml missing")

    import yaml as _yaml
    data = _yaml.safe_load(src.read_text())
    data.setdefault("output", {})
    data["output"]["n_points"] = 40          # small N-M grid
    data["output"]["generate_3d_surface"] = False   # don't build hull
    data["output"]["generate_mx_my"] = False
    data["output"]["generate_moment_curvature"] = False
    data["output"]["generate_polar_ductility"] = False
    data["output"]["generate_3d_moment_curvature"] = False

    quick = tmp_path / "quick.yaml"
    quick.write_text(_yaml.safe_dump(data))

    res = api.analyze(yaml_path=quick)
    assert res.section.B_mm == 400
    assert res.section.H_mm == 600
    assert len(res.demands) >= 1
    assert res.meta.elapsed_ms >= 0
