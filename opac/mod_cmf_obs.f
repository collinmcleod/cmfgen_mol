!
! Module containing basic data for:
!                   (i) Each species
!          	   (ii) Model atom data (populations etc)
!		  (iii) Basic atmospheric structure.
!
	MODULE MOD_CMF_OBS
	USE SET_KIND_MODULE
!
! Number of atomic species (e.g. H, C, N is 3 species).
!
	INTEGER*4, PARAMETER :: NUM_SPECIES=22
!
! Maximum number of ionization stages per species. For H, need at this number
! has to be 2 or higher (as I and II). A setting of 10 implies that we can treat
! full atoms for ION_IX.
!
	INTEGER*4, PARAMETER :: MAX_IONS_PER_SPECIES=15
	INTEGER*4, PARAMETER :: MAX_NUM_IONS=NUM_SPECIES*MAX_IONS_PER_SPECIES
!
! Maximum number of photoionization routes for each species.
!
	INTEGER*4, PARAMETER :: NPHOT_MAX=4
!
! Actual number of ions in calculation. Stored sequentially.
!
	INTEGER*4 NUM_IONS
!
	REAL(KIND=LDP) AT_MASS(NUM_SPECIES)		!Atomic mass in amu
	REAL(KIND=LDP) AT_NO(NUM_SPECIES)		!Atomic number of species
	REAL(KIND=LDP) AT_ABUND(NUM_SPECIES)		!Fractional species abundance
	REAL(KIND=LDP) ABUND_SCALE_FAC(NUM_SPECIES)	!To scale species abundance
	REAL(KIND=LDP) SOL_MASS_FRAC(NUM_SPECIES)	!Solar mass fraction
	REAL(KIND=LDP) SOL_ABUND_HSCL(NUM_SPECIES)	!Solar abundance with H=12.0
!
! Total population density of each species (#/cm^3).
!
	REAL(KIND=LDP), ALLOCATABLE :: POP_SPECIES(:,:)
!
! Conservation equation for each species.
!
	INTEGER*4 EQ_SPECIES(NUM_SPECIES)
!
! Indicate at what location in the ion arrays each species starts
! an ends. This includes the highest ionization stage (e.g. H+).
!
	INTEGER*4 SPECIES_BEG_ID(NUM_SPECIES)
	INTEGER*4 SPECIES_END_ID(NUM_SPECIES)
!
! Principal species abbreviation (e.g. CARB for carbon)
!
	CHARACTER*10 SPECIES(NUM_SPECIES)
!
! Abbreviation for species used to identify ions (e.g. Ca for calcium).
!
	CHARACTER*2  SPECIES_ABR(NUM_SPECIES)
!
	CHARACTER*5  GEN_ION_ID(MAX_IONS_PER_SPECIES)
!
! Link from ion identification to parent species.
!
	INTEGER*4 SPECIES_LNK(MAX_NUM_IONS)
!
! Identification of ion.
!
	CHARACTER*12 ION_ID(MAX_NUM_IONS)
!
! Indicates whether a species (e.g. Carbon) is included.
!
	LOGICAL SPECIES_PRES(NUM_SPECIES)
!
! 
! Data arrays for full atom. These are ordered according to variable type.
! This helps to enure that they fall on correct boudaries etc.
!
	TYPE MODEL_ATOM_DATA
!
	  REAL(KIND=LDP), POINTER :: XzV_F(:,:)		!Level populations in FULL atom
	  REAL(KIND=LDP), POINTER :: XzVLTE_F(:,:)	!LTE level populations in FULL atom
	  REAL(KIND=LDP), POINTER :: W_XzV_F(:,:)	!Level dissolution factors
	  REAL(KIND=LDP), POINTER :: DXzV_F(:)		!Ion population for full atom
	  REAL(KIND=LDP), POINTER :: AXzV_F(:,:)	!Oscillator strength (A(I,j), i<j)
	  REAL(KIND=LDP), POINTER :: EDGEXzV_F(:)	!Ionization energy to g.s. (10^15 Hz)
	  REAL(KIND=LDP), POINTER :: GXzV_F(:)		!Level statistical weights in full atom
	  REAL(KIND=LDP), POINTER :: ARAD(:)		!Inverse radiative lifetime of level
	  REAL(KIND=LDP), POINTER :: GAM2(:)		!Collisional profile parameter.
	  REAL(KIND=LDP), POINTER :: GAM4(:)		!Collisional profile parameter.
!
	  REAL(KIND=LDP), POINTER :: DXzV(:)		!Ion population for super level
	  REAL(KIND=LDP), POINTER :: XzV(:,:)		!Level population in SL atom
	  REAL(KIND=LDP), POINTER :: XzVLTE(:,:)	!LTE populations in SL atom
	  REAL(KIND=LDP), POINTER :: dlnXzVLTE_dlnT(:,:)
!
	  REAL(KIND=LDP) ZXzV			!Charge on ion (=1 for HI)
	  REAL(KIND=LDP) GIONXzV_F		!Statistical weight of ion
!
! Identifications corresponding to each photoionization route.
!
	  INTEGER*4 XzV_ION_LEV_ID(NPHOT_MAX)
	  INTEGER*4, POINTER :: F_TO_S_XzV(:)	!Link of full levels to super levels
	  INTEGER*4, POINTER :: INT_SEQ_XzV(:)
!
	  INTEGER*4 NXzV_F		!Number of levels in full atom
	  INTEGER*4 NXzV		!Number of levels in SL atom
	  INTEGER*4 EQXzV		!Equation in BA matrix for g.s. of atom
	  INTEGER*4 N_XzV_PHOT		!Number of states species can ionize to.
!
	  LOGICAL, POINTER :: OBSERVED_LEVEL(:)	!Link of full levels to super levels
	
! Indicates whether a species is present. The final ionization state is regarded
! as not present (e.g. HII_PRES is ALWAYS false, even though we treat H+ when
! we treat HI).
!
	  LOGICAL XzV_PRES		!indicates
!
! Dielectronic variables.
!
	  LOGICAL DIE_AUTO_XzV
	  LOGICAL DIE_WI_XzV
!
	  CHARACTER*30, POINTER :: XzVLEVNAME_F(:)	!Level name
!
	  CHARACTER*6  XzV_TRANS_TYPE	!Transition type (e.g. Blank, Sobolev etc)
	  CHARACTER*10 XzV_PROF_TYPE	!Type of profile (e.g. Doppler, Stark)
	  CHARACTER*11 XzV_OSCDATE
	  CHARACTER*11 NEW_XzV_OSCDATE
!
	END TYPE MODEL_ATOM_DATA
!
!
	INTEGER*4 EQNE		!Electron conservation equation
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)		!Radius in units of 10^10 cm
	REAL(KIND=LDP), ALLOCATABLE :: V(:)		!V in units of km/s
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)		!dlnV/dlnR-1
	REAL(KIND=LDP), ALLOCATABLE :: T(:)		!Temperature in units of 10^4 K
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)		!Electron density (#/cm^3)
!
	REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)	!Rosseland mean opacity.
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)	!Flux mean opacity
	REAL(KIND=LDP), ALLOCATABLE :: POP_ATOM(:)	!Total atom density (#/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: DENSITY(:)	!Mass density (gm/cm^3)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC(:)	!Volume filling factor for clumps
	REAL(KIND=LDP), ALLOCATABLE :: POPION(:)	!Ion density
!
	REAL(KIND=LDP) STARS_MASS			!In Msun
	REAL(KIND=LDP) STARS_LUM			!In Lsun
!
	TYPE (MODEL_ATOM_DATA) ATM(NUM_SPECIES*MAX_IONS_PER_SPECIES)
!
! Indicates generic ionization names.
!
	DATA GEN_ION_ID /'I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV'/
!
	END MODULE MOD_CMF_OBS
