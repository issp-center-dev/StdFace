"""Keyword parsing helpers and solver-specific parsers.

This module contains functions for parsing keyword-value pairs from
StdFace input files, including duplicate-checking storage helpers and
solver-specific keyword dispatchers for HPhi, mVMC, UHF, and H-wave.
"""
from __future__ import annotations

import cmath
import logging
import math

from .stdface_vals import StdIntList, NaN_i

logger = logging.getLogger(__name__)


_TRIM_TABLE = str.maketrans("", "", " \t:;\"\\\b\v\n\r\0")
"""Translation table for :func:`trim_space_quote`.

Removes: space, tab, colon, semicolon, double-quote, backslash,
backspace (``\\b``), vertical-tab (``\\v``), newline, CR, and null.
"""


def trim_space_quote(text: str) -> str:
    """Remove whitespace, colons, semicolons, quotes and backslashes from *text*.

    Uses a precomputed :func:`str.maketrans` table for efficient
    character deletion.

    Parameters
    ----------
    text : str
        Raw keyword or value string from an input file.

    Returns
    -------
    str
        The cleaned string with the above characters removed.
    """
    return text.translate(_TRIM_TABLE)


def _fail_duplicate(keyword: str) -> None:
    """Log a duplicate-keyword error and raise.

    Parameters
    ----------
    keyword : str
        The duplicated keyword name.

    Raises
    ------
    ValueError
        Always raised after logging the error message.
    """
    msg = f"Keyword {keyword} is duplicated."
    logger.error(msg)
    raise ValueError(msg)


def store_with_check_dup_s(keyword: str, value: str, current: str) -> str:
    """Store a string value after checking for duplicate assignment.

    If *current* has already been assigned (i.e. it is not ``None``),
    the program prints an error and exits.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file.
    current : str
        The current stored value (``None`` means unset).

    Returns
    -------
    str
        The accepted value (equal to *value*).

    Raises
    ------
    ValueError
        If *current* is not the sentinel, indicating a duplicate keyword.
    """
    if current is not None:
        _fail_duplicate(keyword)
    return value


def store_with_check_dup_sl(
    keyword: str, value: str, current: str, maxlen: int = 256
) -> str:
    """Store a string value (forced lower-case) after checking for duplicates.

    Behaves like :func:`store_with_check_dup_s` but additionally
    converts *value* to lower case and truncates it to *maxlen*
    characters.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file.
    current : str
        The current stored value (``None`` means unset).
    maxlen : int, optional
        Maximum number of characters to keep (default 256).

    Returns
    -------
    str
        The accepted value, lower-cased and truncated.

    Raises
    ------
    ValueError
        If *current* is not the sentinel, indicating a duplicate keyword.
    """
    if current is not None:
        _fail_duplicate(keyword)
    return value[:maxlen].lower()


def store_with_check_dup_i(keyword: str, value: str, current: int) -> int:
    """Store an integer value after checking for duplicate assignment.

    Scalar integer fields are unset as ``None``; integer matrix elements
    (``box`` / ``cutoff_*R`` / ``boxsub``) are unset as the sentinel
    ``NaN_i`` (``2147483647``).  Either state means "not yet specified".

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file (will be converted to int).
    current : int
        The current stored value (``None`` or ``NaN_i`` means unset).

    Returns
    -------
    int
        The parsed integer value.

    Raises
    ------
    ValueError
        If *current* is already set, indicating a duplicate keyword.
    """
    if current is not None and current != NaN_i:
        _fail_duplicate(keyword)
    # C sscanf("%d") truncates floats like "2.0" -> 2
    return int(float(value))


def store_with_check_dup_d(keyword: str, value: str, current: float) -> float:
    """Store a float value after checking for duplicate assignment.

    If *current* is **not** NaN the program prints an error and exits.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file (will be converted to float).
    current : float
        The current stored value (``NaN`` means unset).

    Returns
    -------
    float
        The parsed float value.

    Raises
    ------
    ValueError
        If *current* is not NaN, indicating a duplicate keyword.
    """
    if current is not None and not math.isnan(current):
        _fail_duplicate(keyword)
    return float(value)


def _safe_float(s: str) -> float:
    """Parse a string to float; an empty string means 0.0.

    The empty-string default supports the ``",imag"`` / ``"real,"``
    complex forms where one part is omitted.

    Parameters
    ----------
    s : str
        String to parse.

    Returns
    -------
    float
        Parsed value, or 0.0 if *s* is empty.

    Raises
    ------
    ValueError
        If *s* is non-empty and not a valid number.
    """
    s = s.strip()
    if not s:
        return 0.0
    try:
        return float(s)
    except ValueError:
        msg = f"Invalid numerical value: {s!r}"
        logger.error(msg)
        raise ValueError(msg) from None


def store_with_check_dup_c(keyword: str, value: str, current: complex) -> complex:
    """Store a complex value after checking for duplicate assignment.

    The input string *value* may be in one of the following forms:

    * ``"real,imag"`` -- both parts specified
    * ``"real"``      -- imaginary part defaults to 0
    * ``",imag"``     -- real part defaults to 0

    If the real part of *current* is **not** NaN the program prints an
    error and exits.

    Parameters
    ----------
    keyword : str
        The keyword name (used only for the error message).
    value : str
        The new value read from the input file.
    current : complex
        The current stored value (real-part ``NaN`` means unset).

    Returns
    -------
    complex
        The parsed complex value.

    Raises
    ------
    ValueError
        If *current* is already set, indicating a duplicate keyword.
    """
    if current is not None and not cmath.isnan(current):
        _fail_duplicate(keyword)

    # Split on comma, mirroring the C strtok(",") logic
    parts = value.split(",", 1) if "," in value else [value]
    real_part = _safe_float(parts[0])
    imag_part = _safe_float(parts[1]) if len(parts) > 1 else 0.0

    return complex(real_part, imag_part)


# =====================================================================
#  Common keyword dispatch table
# =====================================================================


def _j_matrix_keywords(prefix: str, scalar_field: str,
                       matrix_field: str) -> dict[str, tuple]:
    """Generate keyword entries for a 3x3 spin-coupling matrix and its isotropic scalar.

    Each J-family (J, J0, J0', J0'', J1, ..., J2'', J', J'') has an
    isotropic scalar keyword (e.g. ``"j0"``) and 9 anisotropic component
    keywords (e.g. ``"j0x"``, ``"j0xy"``, ``"j0xz"``, ...).

    Parameters
    ----------
    prefix : str
        Lowered keyword prefix (e.g. ``"j0"``, ``"j0'"``, ``"j2''"``, ``"j"``).
    scalar_field : str
        StdIntList attribute for the isotropic scalar (e.g. ``"J0All"``).
    matrix_field : str
        StdIntList attribute for the 3x3 matrix (e.g. ``"J0"``).

    Returns
    -------
    dict[str, tuple]
        10 keyword entries (1 scalar + 9 array-element).
    """
    # Component suffixes → (row, col) in the 3×3 matrix
    _COMPONENTS = [
        ("x",  (0, 0)), ("xy", (0, 1)), ("xz", (0, 2)),
        ("y",  (1, 1)), ("yx", (1, 0)), ("yz", (1, 2)),
        ("z",  (2, 2)), ("zx", (2, 0)), ("zy", (2, 1)),
    ]
    d: dict[str, tuple] = {
        prefix: (store_with_check_dup_d, scalar_field),
    }
    for suffix, idx in _COMPONENTS:
        d[prefix + suffix] = (store_with_check_dup_d, matrix_field, idx, float)
    return d


def _grid3x3_keywords(
    fmt: str, field: str, store_func: object, cast: type,
) -> dict[str, tuple]:
    """Generate 9 keyword entries for a 3x3 array indexed by (a0/a1/a2) × (w/l/h).

    Parameters
    ----------
    fmt : str
        Format string with ``{a}`` and ``{c}`` placeholders for the row name
        and column name, e.g. ``"{a}{c}"`` → ``"a0w"``, ``"a1l"``, etc.
        or ``"cutoff_j_{a}{c}"`` → ``"cutoff_j_a0w"``, etc.
    field : str
        StdIntList attribute name for the 3x3 array (e.g. ``"box"``).
    store_func : callable
        Store-with-duplicate-check function (e.g. ``store_with_check_dup_i``).
    cast : type
        Type cast for the stored value (``int`` or ``float``).

    Returns
    -------
    dict[str, tuple]
        9 keyword entries mapping ``fmt.format(a=..., c=...)`` to
        ``(store_func, field, (row, col), cast)``.
    """
    d: dict[str, tuple] = {}
    for row, aname in enumerate(("a0", "a1", "a2")):
        for col, comp in enumerate(("w", "l", "h")):
            d[fmt.format(a=aname, c=comp)] = (
                store_func, field, (row, col), cast
            )
    return d


_COMMON_KEYWORDS: dict[str, tuple] = {
    # --- scalar a --------------------------------------------------------
    "a": (store_with_check_dup_d, "a"),
    # --- box (supercell) -------------------------------------------------
    **_grid3x3_keywords("{a}{c}", "box", store_with_check_dup_i, int),
    # --- cutoff J --------------------------------------------------------
    "cutoff_j":        (store_with_check_dup_d, "cutoff_j"),
    "cutoff_jw":       (store_with_check_dup_i, "cutoff_JR", 0, int),
    "cutoff_jl":       (store_with_check_dup_i, "cutoff_JR", 1, int),
    "cutoff_jh":       (store_with_check_dup_i, "cutoff_JR", 2, int),
    **_grid3x3_keywords("cutoff_j_{a}{c}", "cutoff_JVec", store_with_check_dup_d, float),
    "cutoff_length_j": (store_with_check_dup_d, "cutoff_length_J"),
    "cutoff_length_u": (store_with_check_dup_d, "cutoff_length_U"),
    "cutoff_length_t": (store_with_check_dup_d, "cutoff_length_t"),
    # --- cutoff t --------------------------------------------------------
    "cutoff_t":        (store_with_check_dup_d, "cutoff_t"),
    "cutoff_tw":       (store_with_check_dup_i, "cutoff_tR", 0, int),
    "cutoff_tl":       (store_with_check_dup_i, "cutoff_tR", 1, int),
    "cutoff_th":       (store_with_check_dup_i, "cutoff_tR", 2, int),
    **_grid3x3_keywords("cutoff_t_{a}{c}", "cutoff_tVec", store_with_check_dup_d, float),
    # --- cutoff U --------------------------------------------------------
    "cutoff_u":        (store_with_check_dup_d, "cutoff_u"),
    "cutoff_uw":       (store_with_check_dup_i, "cutoff_UR", 0, int),
    "cutoff_ul":       (store_with_check_dup_i, "cutoff_UR", 1, int),
    "cutoff_uh":       (store_with_check_dup_i, "cutoff_UR", 2, int),
    **_grid3x3_keywords("cutoff_u_{a}{c}", "cutoff_UVec", store_with_check_dup_d, float),
    # --- lambda, alpha, D ------------------------------------------------
    "lambda":          (store_with_check_dup_d, "lambda_"),
    "lambda_u":        (store_with_check_dup_d, "lambda_U"),
    "lambda_j":        (store_with_check_dup_d, "lambda_J"),
    "alpha":           (store_with_check_dup_d, "alpha"),
    "d":               (store_with_check_dup_d, "D", (2, 2), float),
    "doublecounting":  (store_with_check_dup_sl, "double_counting_mode"),
    # --- magnetic field and related --------------------------------------
    "gamma":           (store_with_check_dup_d, "Gamma"),
    "h":               (store_with_check_dup_d, "h"),
    "gamma_y":         (store_with_check_dup_d, "Gamma_y"),
    "height":          (store_with_check_dup_i, "Height"),
    "hlength":         (store_with_check_dup_d, "length", 2, float),
    "hx":              (store_with_check_dup_d, "direct", (2, 0), float),
    "hy":              (store_with_check_dup_d, "direct", (2, 1), float),
    "hz":              (store_with_check_dup_d, "direct", (2, 2), float),
    # --- J exchange couplings --------------------------------------------
    **_j_matrix_keywords("j",     "JAll",    "J"),
    **_j_matrix_keywords("j0",    "J0All",   "J0"),
    **_j_matrix_keywords("j0'",   "J0pAll",  "J0p"),
    **_j_matrix_keywords("j0''",  "J0ppAll", "J0pp"),
    **_j_matrix_keywords("j1",    "J1All",   "J1"),
    **_j_matrix_keywords("j1'",   "J1pAll",  "J1p"),
    **_j_matrix_keywords("j1''",  "J1ppAll", "J1pp"),
    **_j_matrix_keywords("j2",    "J2All",   "J2"),
    **_j_matrix_keywords("j2'",   "J2pAll",  "J2p"),
    **_j_matrix_keywords("j2''",  "J2ppAll", "J2pp"),
    **_j_matrix_keywords("j'",    "JpAll",   "Jp"),
    **_j_matrix_keywords("j''",   "JppAll",  "Jpp"),
    # --- K, L, lattice, length -------------------------------------------
    "k":               (store_with_check_dup_d,  "K"),
    "l":               (store_with_check_dup_i,  "L"),
    "lattice":         (store_with_check_dup_sl, "lattice"),
    "llength":         (store_with_check_dup_d, "length", 1, float),
    "lx":              (store_with_check_dup_d, "direct", (1, 0), float),
    "ly":              (store_with_check_dup_d, "direct", (1, 1), float),
    "lz":              (store_with_check_dup_d, "direct", (1, 2), float),
    # --- model, mu, ncond, outputmode ------------------------------------
    "model":           (store_with_check_dup_sl, "model"),
    "mu":              (store_with_check_dup_d,  "mu"),
    "ncond":           (store_with_check_dup_i,  "ncond"),
    "nelec":           (store_with_check_dup_i,  "ncond"),  # alias
    "outputmode":      (store_with_check_dup_sl, "outputmode"),
    # --- phase -----------------------------------------------------------
    "phase0":          (store_with_check_dup_d, "phase", 0, float),
    "phase1":          (store_with_check_dup_d, "phase", 1, float),
    "phase2":          (store_with_check_dup_d, "phase", 2, float),
    # --- hopping t -------------------------------------------------------
    "t":               (store_with_check_dup_c, "t"),
    "t0":              (store_with_check_dup_c, "t0"),
    "t0'":             (store_with_check_dup_c, "t0p"),
    "t0''":            (store_with_check_dup_c, "t0pp"),
    "t1":              (store_with_check_dup_c, "t1"),
    "t1'":             (store_with_check_dup_c, "t1p"),
    "t1''":            (store_with_check_dup_c, "t1pp"),
    "t2":              (store_with_check_dup_c, "t2"),
    "t2'":             (store_with_check_dup_c, "t2p"),
    "t2''":            (store_with_check_dup_c, "t2pp"),
    "t'":              (store_with_check_dup_c, "tp"),
    "t''":             (store_with_check_dup_c, "tpp"),
    # --- U, V ------------------------------------------------------------
    "u":               (store_with_check_dup_d, "U"),
    "v":               (store_with_check_dup_d, "V"),
    "v0":              (store_with_check_dup_d, "V0"),
    "v0'":             (store_with_check_dup_d, "V0p"),
    "v0''":            (store_with_check_dup_d, "V0pp"),
    "v1":              (store_with_check_dup_d, "V1"),
    "v1'":             (store_with_check_dup_d, "V1p"),
    "v1''":            (store_with_check_dup_d, "V1pp"),
    "v2":              (store_with_check_dup_d, "V2"),
    "v2'":             (store_with_check_dup_d, "V2p"),
    "v2''":            (store_with_check_dup_d, "V2pp"),
    "v'":              (store_with_check_dup_d, "Vp"),
    "v''":             (store_with_check_dup_d, "Vpp"),
    # --- W, wlength, wx/wy/wz -------------------------------------------
    "w":               (store_with_check_dup_i, "W"),
    "wlength":         (store_with_check_dup_d, "length", 0, float),
    "wx":              (store_with_check_dup_d, "direct", (0, 0), float),
    "wy":              (store_with_check_dup_d, "direct", (0, 1), float),
    "wz":              (store_with_check_dup_d, "direct", (0, 2), float),
    # --- 2Sz -------------------------------------------------------------
    "2sz":             (store_with_check_dup_i, "Sz2"),
}
"""Common keyword dispatch table shared by all solvers."""


def parse_common_keyword(keyword: str, value: str, StdI: StdIntList) -> bool:
    """Parse a keyword common to all solvers.

    Looks up *keyword* in :data:`_COMMON_KEYWORDS` and applies the
    corresponding store operation to *StdI*.

    Parameters
    ----------
    keyword : str
        The lowered keyword string.
    value : str
        The raw value string from the input file.
    StdI : StdIntList
        The global parameter structure, modified in place.

    Returns
    -------
    bool
        True if the keyword was recognised, False otherwise.
    """
    return _apply_keyword_table(_COMMON_KEYWORDS, keyword, value, StdI)


def _apply_keyword_table(
    table: dict[str, tuple],
    keyword: str,
    value: str,
    StdI: StdIntList,
) -> bool:
    """Look up *keyword* in *table* and apply the corresponding store operation.

    Parameters
    ----------
    table : dict[str, tuple]
        Keyword dispatch table. Each entry is either:
        - ``(store_func, field_name)`` for scalar fields, or
        - ``(store_func, array_name, index, cast)`` for array elements.
    keyword : str
        The lowered keyword string.
    value : str
        The raw value string from the input file.
    StdI : StdIntList
        The global parameter structure, modified in place.

    Returns
    -------
    bool
        True if *keyword* was found in *table*, False otherwise.
    """
    entry = table.get(keyword)
    if entry is None:
        return False
    if len(entry) == 2:
        # Scalar field: (store_func, field_name)
        store_func, field_name = entry
        setattr(StdI, field_name,
                store_func(keyword, value, getattr(StdI, field_name)))
    else:
        # Array element: (store_func, array_name, index, cast)
        store_func, array_name, index, cast = entry
        arr = getattr(StdI, array_name)
        arr[index] = store_func(keyword, value, cast(arr[index]))
    return True


