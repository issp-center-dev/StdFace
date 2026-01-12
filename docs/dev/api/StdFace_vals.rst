StdFace_vals.h
==============

Purpose
-------

Defines the central configuration structure ``StdIntList`` that holds all model
and lattice parameters. This structure is the primary data container passed
throughout the StdFace processing pipeline.

**Source reference**: ``src/StdFace_vals.h``

Public Types
------------

struct StdIntList
^^^^^^^^^^^^^^^^^

Main configuration structure for lattice models and interactions. Contains
approximately 125 fields organized into categories below.

**Initialization Fields**

``int NaN_i``
   Sentinel value for unspecified integer parameters. Set in ``StdFace_ResetVals()``.

``double pi``
   Mathematical constant pi (3.14159...).

**Lattice Parameters**

``char lattice[256]``
   Name of the lattice type (e.g., "chain", "square", "honeycomb").

``double a``
   Lattice constant.

``double length[3]``
   Anisotropic lattice constants (wlength, llength, hlength).

``int W``
   Number of sites along the 1st axis.

``int L``
   Number of sites along the 2nd axis.

``int Height``
   Number of sites along the 3rd axis.

``double direct[3][3]``
   Unit direct lattice vectors. Set in ``StdFace_InitSite()``.

``int box[3][3]``
   Super-cell shape. Input parameters a0W, a0L, a0H, etc.

``int rbox[3][3]``
   Inverse of box matrix. Set in ``StdFace_InitSite()``.

``int NCell``
   Number of unit cells in the super-cell. Set in ``StdFace_InitSite()``.

``int **Cell``
   Cell positions in fractional coordinates. Allocated in ``StdFace_InitSite()``.

``int NsiteUC``
   Number of sites in the unit cell.

``double **tau``
   Cell-internal site positions in fractional coordinates.

**Model Parameters**

``char model[256]``
   Name of the model (e.g., "hubbard", "spin", "kondo").

``double mu``
   Chemical potential.

``double complex t``
   Nearest-neighbor hopping.

``double complex tp``
   2nd-nearest hopping.

``double complex t0``, ``t0p``, ``t0pp``
   Anisotropic hopping (direction 0).

``double complex t1``, ``t1p``, ``t1pp``
   Anisotropic hopping (direction 1).

``double complex t2``, ``t2p``, ``t2pp``
   Anisotropic hopping (direction 2).

``double complex tpp``
   3rd-nearest hopping.

``double U``
   On-site Coulomb potential.

``double V``, ``Vp``, ``Vpp``
   Off-site Coulomb potential (1st, 2nd, 3rd neighbor).

``double V0``, ``V0p``, ``V0pp``
   Anisotropic Coulomb potential (direction 0).

``double V1``, ``V1p``, ``V1pp``
   Anisotropic Coulomb potential (direction 1).

``double V2``, ``V2p``, ``V2pp``
   Anisotropic Coulomb potential (direction 2).

**Exchange Parameters**

``double JAll``, ``JpAll``, ``JppAll``
   Isotropic diagonal spin coupling (1st, 2nd, 3rd neighbor).

``double J0All``, ``J0pAll``, ``J0ppAll``
   Anisotropic diagonal spin coupling (direction 0).

``double J1All``, ``J1pAll``, ``J1ppAll``
   Anisotropic diagonal spin coupling (direction 1).

``double J2All``, ``J2pAll``, ``J2ppAll``
   Anisotropic diagonal spin coupling (direction 2).

``double J[3][3]``, ``Jp[3][3]``, ``Jpp[3][3]``
   Full 3x3 exchange coupling matrices.

``double J0[3][3]``, ``J0p[3][3]``, ``J0pp[3][3]``
   Anisotropic 3x3 exchange matrices (direction 0).

``double J1[3][3]``, ``J1p[3][3]``, ``J1pp[3][3]``
   Anisotropic 3x3 exchange matrices (direction 1).

``double J2[3][3]``, ``J2p[3][3]``, ``J2pp[3][3]``
   Anisotropic 3x3 exchange matrices (direction 2).

``double D[3][3]``
   Single-ion anisotropy. Only D[2][2] is used.

``double h``
   Longitudinal magnetic field.

``double Gamma``
   Transverse magnetic field (x-direction).

``double Gamma_y``
   Transverse magnetic field (y-direction).

**Boundary Phase**

``double phase[3]``
   Boundary phase in degrees (phase0, phase1, phase2).

``double complex ExpPhase[3]``
   Computed as exp(i \* pi \* phase / 180).

``int AntiPeriod[3]``
   Set to 1 if corresponding phase = 180.

**Transfer and Interaction Arrays**

``int nsite``
   Total number of sites.

``int *locspinflag``
   LocSpin flags per site. Array of size nsite.

``int ntrans``
   Number of transfer (one-body) terms.

``int **transindx``
   Site/spin indices for transfer terms. Shape [ntrans][4].

``double complex *trans``
   Coefficients for transfer terms. Array of size ntrans.

``int nintr``
   Number of interaction (two-body) terms.

``int **intrindx``
   Site/spin indices for interaction terms. Shape [nintr][8].

``double complex *intr``
   Coefficients for interaction terms. Array of size nintr.

**Specialized Interaction Arrays**

``int NCintra``, ``int **CintraIndx``, ``double *Cintra``
   Intra-site Coulomb terms.

``int NCinter``, ``int **CinterIndx``, ``double *Cinter``
   Inter-site Coulomb terms.

``int NHund``, ``int **HundIndx``, ``double *Hund``
   Hund coupling terms.

``int NEx``, ``int **ExIndx``, ``double *Ex``
   Exchange terms.

``int NPairLift``, ``int **PLIndx``, ``double *PairLift``
   Pair-lift terms.

``int NPairHopp``, ``int **PHIndx``, ``double *PairHopp``
   Pair-hopping terms.

**Calculation Conditions**

``int ncond``
   Number of electrons.

``int lGC``
   Grand canonical ensemble switch (1 = enabled).

``int S2``
   Total spin 2S of local spins.

``char outputmode[256]``
   Output verbosity mode.

``char CDataFileHead[256]``
   Header for output files.

``int Sz2``
   Total Sz.

**Wannier90 Mode Parameters**

``double cutoff_t``, ``cutoff_u``, ``cutoff_j``
   Cutoff values for hopping, Coulomb, and Hund terms.

``double cutoff_length_t``, ``cutoff_length_U``, ``cutoff_length_J``
   Cutoff lengths for R vectors.

``double lambda``, ``lambda_U``, ``lambda_J``
   Tuning parameters for U and J.

``char double_counting_mode[256]``
   Mode for double counting correction.

``double alpha``
   Chemical potential correction parameter.

**Mode-Specific Fields**

Additional fields are defined conditionally based on build mode:

- **HPhi mode** (``_HPhi``): Lanczos parameters, spectrum calculation settings,
  time evolution parameters
- **mVMC mode** (``_mVMC``): Variational parameters, sublattice settings,
  orbital indices
- **UHF mode** (``_UHF``): Mean-field iteration parameters, sublattice settings
- **HWAVE mode** (``_HWAVE``): Similar to UHF plus Wannier90 export settings

Public Macros/Constants
-----------------------

None.

Public Functions
----------------

None defined in this header. Functions that operate on ``StdIntList`` are declared
in ``StdFace_ModelUtil.h``.

Ownership and Lifetime Rules
----------------------------

- The ``StdIntList`` structure is allocated in ``StdFace_main()`` via ``malloc()``.
- Fields are initialized by ``StdFace_ResetVals()`` (defined in ``StdFace_main.c``).
- Internal arrays (``trans``, ``intr``, ``Cell``, etc.) are allocated by
  ``StdFace_MallocInteractions()`` and lattice constructors.
- All allocated memory is freed at the end of ``StdFace_main()``.

Error Handling
--------------

No error handling in this header (structure definition only).

Thread / MPI Safety
-------------------

Unspecified in the current code.

Source Reference
----------------

- Header: ``src/StdFace_vals.h``
- Initialization: ``StdFace_ResetVals()`` in ``src/StdFace_main.c``
