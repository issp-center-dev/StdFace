"""Wannier90 reader: leading-orbital projection, record-indexed distance
cutoff and layout tolerance (C: 0dd1c41, 6319bbb; upstream
``test/check_wannier_reader.cmake``)."""
from __future__ import annotations

import pytest

from stdface import generate

_TYPES = ("hr", "ur", "jr", "dr")
_GEOM_ORIGIN = "1 0 0\n0 1 0\n0 0 1\n{n}\n" + "0 0 0\n" * 2
_GEOM_OFFSET = "1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0.45 0 0\n"


def _rows(n, coordinates, values):
    """Records in RESPACK order (j fastest), as upstream ``write_wannier``."""
    rows, k = [], 0
    for r in coordinates:
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                rows.append(f"{r} {i} {j} {values[k]}")
                k += 1
    return rows


def _text(n, coordinates, rows, degeneracy=None, comment="Wannier fixture"):
    if degeneracy is None:
        degeneracy = [" ".join("1" for _ in coordinates)]
    return "\n".join([comment, str(n), str(len(coordinates)), *degeneracy, *rows]) + "\n"


def _two_orbital_files():
    return {
        "hr": ["-0.2 0", "0.15 0.1", "0.15 -0.1", "-0.3 0"],
        "ur": ["4 0", "0 0", "0 0", "8 0"],
        "jr": ["0 0", "0 0", "0 0", "0 0"],
        "dr": ["0.5 0", "0.1 0.2", "0.1 -0.2", "0.75 0"],
    }


def _run(work, monkeypatch, files, geom, extra=None):
    work.mkdir()
    (work / "zvo_geom.dat").write_text(geom)
    for t, text in files.items():
        (work / f"zvo_{t}.dat").write_bytes(text.encode())
    monkeypatch.chdir(work)
    data = {"lattice": "wannier90", "model": "HubbardGC", "W": 3, "L": 1,
            "Height": 1, "method": "CG", "doublecounting": "hartree"}
    data.update(extra or {})
    out = work / "out"
    generate(data, solver="HPhi", output_dir=out)
    return {p.name: p.read_text() for p in out.iterdir() if p.suffix == ".def"}


def _valid_files():
    return {t: _text(2, ["0 0 0"], _rows(2, ["0 0 0"], v))
            for t, v in _two_orbital_files().items()}


def test_geometry_selects_leading_orbitals(tmp_path, monkeypatch):
    """A 1-site geometry on 2-orbital files equals genuine 1-orbital files."""
    one = {t: _text(1, ["0 0 0"], _rows(1, ["0 0 0"], [v[0]]))
           for t, v in _two_orbital_files().items()}
    expected = _run(tmp_path / "one", monkeypatch, one,
                    "1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n")
    # The unused second centre stays in the geometry file, as upstream.
    projected = _run(tmp_path / "projection", monkeypatch, _valid_files(),
                     "1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n0.2 0 0\n")
    assert projected == expected
    assert "trans.def" in projected


_RS = ["-1 0 0", "0 0 0", "1 0 0"]
_CUT = {"cutoff_length_t": 0.6, "cutoff_length_U": 0.6, "cutoff_length_J": 0.6}
_NOCUT = {"cutoff_length_t": -1, "cutoff_length_U": -1, "cutoff_length_J": -1}


def _offset_values(filtered):
    if filtered:
        neg, pos = ["0 0", "0.3 0.2", "0 0", "0 0"], ["0 0", "0 0", "0.3 -0.2", "0 0"]
    else:
        neg, pos = ["0.1 0", "0.3 0.2", "0.4 -0.1", "0.2 0"], ["0.1 0", "0.4 0.1", "0.3 -0.2", "0.2 0"]
    centre = {"hr": ["-0.2 0", "0.15 0.1", "0.15 -0.1", "-0.3 0"],
              "ur": ["4 0", "0.2 0", "0.2 0", "8 0"],
              "jr": ["0 0", "0.05 0", "0.05 0", "0 0"]}
    values = {t: neg + c + pos for t, c in centre.items()}
    values["dr"] = ["0.1 0", "0.2 0.1", "0.3 -0.1", "0.2 0", "0.5 0", "0.1 0.2",
                    "0.1 -0.2", "0.75 0", "0.1 0", "0.3 0.1", "0.2 -0.1", "0.2 0"]
    return values


def _reorder(rows, order):
    return [rows[4 * b + k] for b in range(3) for k in order]


@pytest.mark.parametrize("order, degeneracy", [
    ((0, 1, 2, 3), None),                 # RESPACK: j fastest
    ((0, 2, 1, 3), ["1 1 1"]),            # native Wannier90: i fastest
    ((3, 2, 1, 0), ["1", "2 1"]),         # reversed rows, wrapped degeneracy
], ids=["respack", "native", "reversed"])
def test_distance_cutoff_uses_record_orbitals(tmp_path, monkeypatch, order, degeneracy):
    """With offset centres, the length cutoff follows the record's i, j,
    not the row's position inside its R block."""
    filtered = {t: _text(2, _RS, _rows(2, _RS, v))
                for t, v in _offset_values(True).items()}
    expected = _run(tmp_path / "filtered", monkeypatch, filtered, _GEOM_OFFSET, _NOCUT)
    files = {t: _text(2, _RS, _reorder(_rows(2, _RS, v), order), degeneracy)
             for t, v in _offset_values(False).items()}
    assert _run(tmp_path / "cut", monkeypatch, files, _GEOM_OFFSET, _CUT) == expected


@pytest.mark.parametrize("ftype", _TYPES)
def test_layout_does_not_change_output(tmp_path, monkeypatch, ftype):
    """Reversed rows, CRLF, blank lines, a long comment and no final
    newline read the same as the plain file."""
    expected = _run(tmp_path / "valid", monkeypatch, _valid_files(),
                    _GEOM_ORIGIN.format(n=2))
    files = _valid_files()
    rows = _rows(2, ["0 0 0"], _two_orbital_files()[ftype])[::-1]
    comment = "long header " + "0123456789" * 40
    files[ftype] = "\r\n".join([comment, "2", "", "1", "2", *rows]).replace(
        rows[1], rows[1] + "\r\n  \t")
    assert _run(tmp_path / "layout", monkeypatch, files,
                _GEOM_ORIGIN.format(n=2)) == expected


def test_geometry_blank_lines(tmp_path, monkeypatch):
    expected = _run(tmp_path / "plain", monkeypatch, _valid_files(),
                    _GEOM_ORIGIN.format(n=2))
    geom = "\n1 0 0\n0 1 0\r\n\n0 0 1\n2\n\n0 0 0\n0 0 0"
    assert _run(tmp_path / "blank", monkeypatch, _valid_files(), geom) == expected


@pytest.mark.parametrize("ftype", _TYPES)
def test_fewer_orbitals_than_geometry_is_an_input_error(tmp_path, monkeypatch, ftype):
    """Absent orbitals are not filled with zeros."""
    files = _valid_files()
    files[ftype] = _text(1, ["0 0 0"], _rows(1, ["0 0 0"], ["0.5 0"]))
    with pytest.raises(ValueError, match=rf"zvo_{ftype}\.dat: file has fewer orbitals"):
        _run(tmp_path / "few", monkeypatch, files, _GEOM_ORIGIN.format(n=2))


@pytest.mark.parametrize("ftype", _TYPES)
def test_truncated_body_is_an_input_error(tmp_path, monkeypatch, ftype):
    files = _valid_files()
    files[ftype] = "\n".join(files[ftype].splitlines()[:-1]) + "\n"
    with pytest.raises(ValueError, match=rf"zvo_{ftype}\.dat: unexpected end of file"):
        _run(tmp_path / "eof", monkeypatch, files, _GEOM_ORIGIN.format(n=2))
