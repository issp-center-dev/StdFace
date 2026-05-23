StdFace_ModelUtil.h
===================

Purpose
-------

Core model utilities for lattice construction and physical interactions. This header
provides functions to:

- Construct various lattice geometries (chain, honeycomb, kagome, etc.)
- Set up hopping and interaction terms
- Handle input parameter validation and printing
- Manage site labeling and memory allocation for interaction arrays

**Source reference**: ``src/StdFace_ModelUtil.h``

Public Types
------------

None defined in this header. Uses ``struct StdIntList`` from ``StdFace_vals.h``.

Public Macros/Constants
-----------------------

None.

Public Functions
----------------

Error/Exit Functions
^^^^^^^^^^^^^^^^^^^^

.. c:function:: void StdFace_exit(int errorcode)

   Terminates the program with the given error code.

   :param errorcode: Exit status code

Interaction Generator Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void StdFace_intr(struct StdIntList *StdI, double complex intr0, int site1, int spin1, int site2, int spin2, int site3, int spin3, int site4, int spin4)

   Adds a general two-body interaction term to the interaction list.

   :param StdI: Pointer to the configuration structure
   :param intr0: Complex coefficient of the interaction
   :param site1: First site index
   :param spin1: First spin index
   :param site2: Second site index
   :param spin2: Second spin index
   :param site3: Third site index
   :param spin3: Third spin index
   :param site4: Fourth site index
   :param spin4: Fourth spin index

.. c:function:: void StdFace_Hopping(struct StdIntList *StdI, double complex trans0, int isite, int jsite, double *dR)

   Adds a hopping term between two sites.

   :param StdI: Pointer to the configuration structure
   :param trans0: Complex hopping amplitude
   :param isite: Source site index
   :param jsite: Destination site index
   :param dR: Displacement vector (3 elements)

.. c:function:: void StdFace_trans(struct StdIntList *StdI, double complex trans0, int isite, int ispin, int jsite, int jspin)

   Adds a transfer (one-body) term with explicit spin indices.

   :param StdI: Pointer to the configuration structure
   :param trans0: Complex transfer amplitude
   :param isite: Source site index
   :param ispin: Source spin index
   :param jsite: Destination site index
   :param jspin: Destination spin index

.. c:function:: void StdFace_HubbardLocal(struct StdIntList *StdI, double mu0, double h0, double Gamma0, double Gamma0_y, double U0, int isite)

   Sets up local Hubbard model terms (chemical potential, magnetic field, on-site U).

   :param StdI: Pointer to the configuration structure
   :param mu0: Chemical potential
   :param h0: Longitudinal magnetic field
   :param Gamma0: Transverse magnetic field (x)
   :param Gamma0_y: Transverse magnetic field (y)
   :param U0: On-site Coulomb interaction
   :param isite: Site index

.. c:function:: void StdFace_MagField(struct StdIntList *StdI, int S2, double h, double Gamma, double Gamma_y, int isite)

   Adds magnetic field terms for spin models.

   :param StdI: Pointer to the configuration structure
   :param S2: Total spin 2S
   :param h: Longitudinal field
   :param Gamma: Transverse field (x)
   :param Gamma_y: Transverse field (y)
   :param isite: Site index

.. c:function:: void StdFace_Coulomb(struct StdIntList *StdI, double V, int isite, int jsite)

   Adds inter-site Coulomb interaction.

   :param StdI: Pointer to the configuration structure
   :param V: Coulomb interaction strength
   :param isite: First site index
   :param jsite: Second site index

.. c:function:: void StdFace_GeneralJ(struct StdIntList *StdI, double J[3][3], int Si2, int Sj2, int isite, int jsite)

   Adds general exchange interaction with full 3x3 J matrix.

   :param StdI: Pointer to the configuration structure
   :param J: 3x3 exchange coupling matrix
   :param Si2: Spin 2S at site i
   :param Sj2: Spin 2S at site j
   :param isite: First site index
   :param jsite: Second site index

Print Utility Functions
^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void StdFace_PrintVal_d(char *valname, double *val, double val0)

   Prints a double parameter value; sets default if unspecified.

   :param valname: Name of the parameter
   :param val: Pointer to the value (may be modified)
   :param val0: Default value

.. c:function:: void StdFace_PrintVal_dd(char *valname, double *val, double val0, double val1)

   Prints a double parameter with two possible defaults.

   :param valname: Name of the parameter
   :param val: Pointer to the value
   :param val0: Primary default value
   :param val1: Secondary default value

.. c:function:: void StdFace_PrintVal_c(char *valname, double complex *val, double complex val0)

   Prints a complex parameter value.

   :param valname: Name of the parameter
   :param val: Pointer to the complex value
   :param val0: Default complex value

.. c:function:: void StdFace_PrintVal_i(char *valname, int *val, int val0)

   Prints an integer parameter value.

   :param valname: Name of the parameter
   :param val: Pointer to the integer value
   :param val0: Default value

Validation Functions
^^^^^^^^^^^^^^^^^^^^

.. c:function:: void StdFace_NotUsed_d(char *valname, double val)

   Warns if a double parameter was specified but is not used in the current model.

   :param valname: Name of the parameter
   :param val: The specified value

.. c:function:: void StdFace_NotUsed_i(char *valname, int val)

   Warns if an integer parameter was specified but is not used.

   :param valname: Name of the parameter
   :param val: The specified value

.. c:function:: void StdFace_NotUsed_c(char *valname, double complex val)

   Warns if a complex parameter was specified but is not used.

   :param valname: Name of the parameter
   :param val: The specified value

.. c:function:: void StdFace_NotUsed_J(char *valname, double JAll, double J[3][3])

   Warns if exchange parameters were specified but not used.

   :param valname: Name of the J parameter
   :param JAll: Isotropic J value
   :param J: 3x3 J matrix

.. c:function:: void StdFace_RequiredVal_i(char *valname, int val)

   Checks that a required integer parameter was specified; exits if not.

   :param valname: Name of the parameter
   :param val: The value to check

Input Processing Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void StdFace_InputSpinNN(double J[3][3], double JAll, double J0[3][3], double J0All, char *J0name)

   Processes nearest-neighbor spin coupling input.

   :param J: Output 3x3 matrix
   :param JAll: Isotropic coupling
   :param J0: Anisotropic coupling matrix
   :param J0All: Anisotropic isotropic coupling
   :param J0name: Parameter name for messages

.. c:function:: void StdFace_InputSpin(double Jp[3][3], double JpAll, char *Jpname)

   Processes spin coupling input (general).

   :param Jp: Output 3x3 matrix
   :param JpAll: Isotropic coupling
   :param Jpname: Parameter name for messages

.. c:function:: void StdFace_InputCoulombV(double V, double *V0, char *V0name)

   Processes Coulomb interaction input.

   :param V: General V value
   :param V0: Pointer to specific V value
   :param V0name: Parameter name

.. c:function:: void StdFace_InputHopp(double complex t, double complex *t0, char *t0name)

   Processes hopping parameter input.

   :param t: General hopping value
   :param t0: Pointer to specific hopping value
   :param t0name: Parameter name

Site/Lattice Utility Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void StdFace_InitSite(struct StdIntList *StdI, FILE *fp, int dim)

   Initializes site arrays and lattice geometry.

   :param StdI: Pointer to the configuration structure
   :param fp: Output file pointer for geometry information
   :param dim: Lattice dimensionality (1, 2, or 3)

.. c:function:: void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp, int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC, int *isite, int *jsite, int connect, double complex *Cphase, double *dR)

   Sets site labels and computes connection information (2D lattices).

   :param StdI: Configuration structure
   :param fp: Output file pointer
   :param iW: Cell index (W direction)
   :param iL: Cell index (L direction)
   :param diW: Displacement in W direction
   :param diL: Displacement in L direction
   :param isiteUC: Source unit cell site index
   :param jsiteUC: Target unit cell site index
   :param isite: Output source global site index
   :param jsite: Output target global site index
   :param connect: Connection type
   :param Cphase: Output boundary phase
   :param dR: Output displacement vector

.. c:function:: void StdFace_PrintGeometry(struct StdIntList *StdI)

   Outputs geometry information to file.

   :param StdI: Configuration structure

.. c:function:: void StdFace_MallocInteractions(struct StdIntList *StdI, int ntransMax, int nintrMax)

   Allocates memory for transfer and interaction arrays.

   :param StdI: Configuration structure
   :param ntransMax: Maximum number of transfer terms
   :param nintrMax: Maximum number of interaction terms

.. c:function:: void StdFace_FindSite(struct StdIntList *StdI, int iW, int iL, int iH, int diW, int diL, int diH, int isiteUC, int jsiteUC, int *isite, int *jsite, double complex *Cphase, double *dR)

   Finds site indices with full 3D support.

   :param StdI: Configuration structure
   :param iW: Cell index (W direction)
   :param iL: Cell index (L direction)
   :param iH: Cell index (H direction)
   :param diW: Displacement in W direction
   :param diL: Displacement in L direction
   :param diH: Displacement in H direction
   :param isiteUC: Source unit cell site index
   :param jsiteUC: Target unit cell site index
   :param isite: Output source global site index
   :param jsite: Output target global site index
   :param Cphase: Output boundary phase
   :param dR: Output displacement vector

.. c:function:: void StdFace_PrintXSF(struct StdIntList *StdI)

   Outputs XSF format file for visualization.

   :param StdI: Configuration structure

Lattice Constructor Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void StdFace_Chain(struct StdIntList *StdI)

   Constructs a 1D chain lattice.

   :param StdI: Configuration structure (modified in place)

.. c:function:: void StdFace_Ladder(struct StdIntList *StdI)

   Constructs a 1D ladder lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_Tetragonal(struct StdIntList *StdI)

   Constructs a 2D square/tetragonal lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_Triangular(struct StdIntList *StdI)

   Constructs a 2D triangular lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_Honeycomb(struct StdIntList *StdI)

   Constructs a 2D honeycomb lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_Kagome(struct StdIntList *StdI)

   Constructs a 2D Kagome lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_Orthorhombic(struct StdIntList *StdI)

   Constructs a 3D orthorhombic/cubic lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_FCOrtho(struct StdIntList *StdI)

   Constructs a 3D face-centered orthorhombic/cubic lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_Pyrochlore(struct StdIntList *StdI)

   Constructs a 3D pyrochlore lattice.

   :param StdI: Configuration structure

.. c:function:: void StdFace_Wannier90(struct StdIntList *StdI)

   Constructs a lattice from Wannier90 tight-binding data.

   :param StdI: Configuration structure

Conditional Functions (Mode-specific)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions are only available in specific build modes:

**HPhi mode** (``_HPhi`` defined):

- ``StdFace_Chain_Boost()``
- ``StdFace_Ladder_Boost()``
- ``StdFace_Honeycomb_Boost()``
- ``StdFace_Kagome_Boost()``

**mVMC mode** (``_mVMC`` defined):

- ``StdFace_generate_orb()``
- ``StdFace_Proj()``
- ``PrintJastrow()``

Ownership and Lifetime Rules
----------------------------

- The ``StdIntList *StdI`` parameter is never allocated or freed by these functions;
  the caller (``StdFace_main()``) is responsible for its lifetime.
- ``StdFace_MallocInteractions()`` allocates internal arrays within ``StdI``; these
  are freed by the caller after output generation.
- Memory allocation uses functions from ``setmemory.h``.

Error Handling
--------------

- ``StdFace_exit(errorcode)`` terminates the program with the given code.
- ``StdFace_RequiredVal_i()`` calls ``StdFace_exit()`` if a required parameter is missing.
- Invalid parameter combinations result in error messages printed to stdout followed
  by program termination.

Thread / MPI Safety
-------------------

No threading constructs are present.  MPI usage is limited to abort/finalize
on fatal exit in ``StdFace_exit()`` (guarded by ``#ifdef MPI``).

Source Reference
----------------

- Header: ``src/StdFace_ModelUtil.h``
- Implementation: ``src/StdFace_ModelUtil.c``
- Lattice implementations: ``src/ChainLattice.c``, ``src/SquareLattice.c``, etc.
