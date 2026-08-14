#!/usr/bin/env python3
"""StdFace パーサーレジストリ全体からキーワード目録(JSON)を生成する。

母集合: core の ``_COMMON_KEYWORDS`` + ソルバーレジストリ経由で列挙した
全プラグインのキーワードテーブル(クラスのハードコード禁止)+
格子名/模型名の alias(lattice/model registry から)。

出力(canonical 単位に集約、安定ソート)は以下の 2 種類のエントリからなる:

- keyword 種別::

    {"keyword_canonical": "J0x", "keyword_lower": "j0x",
     "sources": ["common"], "kind": "keyword"}

- alias 種別(lattice_alias / model_alias)::

    {"keyword_canonical": "chain", "sources": ["lattice_registry"],
     "kind": "lattice_alias", "canonical_target": "chain"}

ソルバーの列挙は ``stdface.plugin`` の登録済みプラグイン(id で重複排除)
を走査して行い、4 つの具象クラス名をこのスクリプトにハードコードしない。
これにより将来プラグインが追加された場合も自動的に目録へ反映される。
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from stdface.core.keyword_parser import _COMMON_KEYWORDS  # noqa: E402
from stdface.core.stdface_vals import ModelType  # noqa: E402
from stdface.lattice import get_all_lattices  # noqa: E402
from stdface import plugin as _solver_plugin_mod  # noqa: E402

_CANON_EXCEPTIONS = {"2s": "2S", "2sz": "2Sz", "gamma": "Gamma",
                     "gamma_y": "Gamma_y", "u": "U", "v": "V", "j": "J",
                     "d": "D", "k": "K", "l": "L", "w": "W", "h": "h",
                     "mu": "mu"}


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
    if kw and kw[0] in "jv" and len(kw) > 1:      # J/V 族は先頭大文字
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

    # --- model aliases -------------------------------------------------
    for model in ModelType:
        key = (str(model.value), "model_alias")
        agg[key].add("model_registry")
        targets[key] = str(model.value)

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
    """Print the keyword/alias inventory as JSON to stdout."""
    entries = build_inventory()
    json.dump(entries, sys.stdout, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
