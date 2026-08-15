#!/usr/bin/env python3
"""Generate the keyword inventory (JSON) from the full StdFace parser registry.

Population: core's ``_COMMON_KEYWORDS`` + the keyword tables of every
plugin enumerated via the solver registry (no hardcoding of classes),
plus lattice-name/model-name aliases (from the lattice/model registry).

The output (aggregated by canonical unit, stably sorted) consists of two
kinds of entries:

- ``keyword`` kind::

    {"keyword_canonical": "J0x", "keyword_lower": "j0x",
     "sources": ["common"], "kind": "keyword"}

- alias kind (lattice_alias / model_alias)::

    {"keyword_canonical": "chain", "sources": ["lattice_registry"],
     "kind": "lattice_alias", "canonical_target": "chain"}

Solvers are enumerated by scanning the registered plugins in
``stdface.plugin`` (de-duplicated by id); the four concrete class names
are not hardcoded into this script. This means the inventory picks up
future plugins automatically once they are added.

``model_alias`` scans every key of
``stdface.core.stdface_main.MODEL_ALIASES`` +
``MODEL_ALIASES_HPHI_BOOST`` (``fermionhubbard``, ``hubbardgc``,
``spingc``, ``kondolattice``, ``kondogc``, ``spingcboost``, etc. — every
alias including the GC/Boost extensions) — it does not hardcode only
the 3 canonical ``ModelType`` values (spin/hubbard/kondo) (the A1 fix).
Entries originating from ``MODEL_ALIASES_HPHI_BOOST`` are tagged with
``"model_registry_hphi_boost"`` in ``sources`` so that HPhi-only
aliases can be distinguished.

If importing the registry itself fails (e.g. ``python/`` is not on
PYTHONPATH), or if the resolved result is empty, exit non-zero without
printing an incomplete inventory (A6-iv).
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

try:
    from stdface.core.keyword_parser import _COMMON_KEYWORDS  # noqa: E402
    from stdface.lattice import get_all_lattices  # noqa: E402
    from stdface import plugin as _solver_plugin_mod  # noqa: E402
    # Model-name alias registry (A1: this is the authoritative population
    # for model_alias).
    from stdface.core.stdface_main import (  # noqa: E402
        MODEL_ALIASES as _MODEL_ALIASES,
        MODEL_ALIASES_HPHI_BOOST as _MODEL_ALIASES_HPHI_BOOST,
    )
except ImportError as e:  # pragma: no cover - environment/setup error
    sys.exit(
        "keyword_inventory: failed to import the stdface registry "
        f"(check that python/ is on PYTHONPATH): {e!r}"
    )

_CANON_EXCEPTIONS = {"2s": "2S", "2sz": "2Sz", "gamma": "Gamma",
                     "gamma_y": "Gamma_y", "u": "U", "v": "V", "j": "J",
                     "d": "D", "k": "K", "l": "L", "w": "W", "h": "h",
                     "mu": "mu", "height": "Height"}


def canon(kw: str) -> str:
    """Canonicalise a lower-cased StdFace keyword for inventory display.

    Parameters
    ----------
    kw : str
        Raw (lower-cased) keyword as it appears in a keyword dispatch table.

    Returns
    -------
    str
        The canonical display form (e.g. ``"j0'x"`` -> ``"J0'x"``).
    """
    if kw in _CANON_EXCEPTIONS:
        return _CANON_EXCEPTIONS[kw]
    if kw and kw[0] in "jv" and len(kw) > 1:      # J/V family: capitalize first letter
        return kw[0].upper() + kw[1:]
    return kw


def iter_solver_plugins():
    """Yield ``(name, keyword_table)`` for every *uniquely registered*
    solver plugin in :mod:`stdface.plugin`.

    The registry maps every alias (e.g. ``HWAVE``, ``UHFR``, ``UHFK``) to
    the same plugin *instance*; this generator de-duplicates by identity
    so each concrete plugin is visited exactly once, regardless of how
    many aliases it is registered under.

    Yields
    ------
    tuple[str, dict[str, tuple]]
        The plugin's canonical name and its ``keyword_table``.
    """
    _solver_plugin_mod._discover_plugins()
    seen: set[int] = set()
    for pl in _solver_plugin_mod._plugins.values():
        if id(pl) in seen:
            continue
        seen.add(id(pl))
        name = pl.name.value if hasattr(pl.name, "value") else str(pl.name)
        yield name, pl.keyword_table


def build_inventory() -> list[dict]:
    """Build the full keyword/alias inventory as a list of entry dicts.

    Returns
    -------
    list[dict]
        Stably sorted, de-duplicated inventory entries. See module
        docstring for the entry shapes.
    """
    agg: dict[tuple[str, str], set[str]] = defaultdict(set)   # (canonical, kind) -> sources
    targets: dict[tuple[str, str], str] = {}                  # (canonical, kind) -> canonical_target

    # --- core common keywords ---------------------------------------
    for kw in _COMMON_KEYWORDS:
        agg[(canon(kw), "keyword")].add("common")

    # --- solver plugin keyword tables (registry-driven) -------------
    for solver_name, table in iter_solver_plugins():
        for kw in table:
            agg[(canon(kw), "keyword")].add(solver_name)

    # --- lattice aliases (registry-driven) ---------------------------
    for lat in get_all_lattices():
        for alias in lat.aliases:
            key = (alias, "lattice_alias")
            agg[key].add("lattice_registry")
            targets[key] = lat.name

    # --- model aliases (registry-driven; MODEL_ALIASES + the HPhi-only ----
    #     Boost-extension aliases in MODEL_ALIASES_HPHI_BOOST) -------------
    for alias, cfg in _MODEL_ALIASES.items():
        key = (alias, "model_alias")
        agg[key].add("model_registry")
        targets[key] = str(cfg.model.value)
    for alias, cfg in _MODEL_ALIASES_HPHI_BOOST.items():
        key = (alias, "model_alias")
        agg[key].add("model_registry_hphi_boost")
        targets[key] = str(cfg.model.value)

    entries = []
    for (c, k), s in sorted(agg.items()):
        entry = {"keyword_canonical": c, "kind": k, "sources": sorted(s)}
        if k == "keyword":
            entry["keyword_lower"] = c.lower()
        else:
            entry["canonical_target"] = targets[(c, k)]
        entries.append(entry)
    return entries


def main() -> None:
    """Print the keyword/alias inventory as JSON to stdout.

    Exits non-zero (without printing a partial inventory) if the solver
    plugin registry or the lattice registry resolved to zero entries —
    that always indicates a broken registry lookup rather than a
    legitimately empty catalog, so a silent partial inventory would be
    worse than a hard failure (A6-iv).
    """
    entries = build_inventory()
    lattice_aliases = [e for e in entries if e["kind"] == "lattice_alias"]
    model_aliases = [e for e in entries if e["kind"] == "model_alias"]
    solver_names = [name for name, _table in iter_solver_plugins()]
    if not lattice_aliases or not model_aliases or not solver_names:
        sys.exit(
            "keyword_inventory: the registry returned an empty result "
            f"(lattice_alias={len(lattice_aliases)}, "
            f"model_alias={len(model_aliases)}, "
            f"solver_plugins={len(solver_names)}) — "
            "exiting without printing an incomplete inventory"
        )
    json.dump(entries, sys.stdout, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
