	SUBROUTINE RD_NUC_DECAY_DATA(INCL_RAD_DECAYS,ND,LU)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
! Altered 23-Jan-2015 : Code checks for Date record -- for bookeeping.
! Altered 06-Jan-2015 : Code can read in the kinetic energy of the positrons.
!                       It assumes that the kinetic energy appears after the total decay energy.
!                       Old format file can still be read.
!
! Altered 15-Jul-2010 : Error check that NUM_ISOTOPES and NUM_DEACY_PATHS are compatible
!                         with the maximum values set in NUC_ISO_MOD inserted.
!
	INTEGER ND
	INTEGER LU
!
	INTEGER IN		!Index for NUC
	INTEGER IN_LOOP
	INTEGER IS		!Index for ISO
	INTEGER IP		!Index for PAR
	INTEGER J
	INTEGER L
	CHARACTER(LEN=132) STRING
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) ELECTRON_VOLT
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU,ELECTRON_VOLT
	LOGICAL UNRECOGNIZED_ISO
	LOGICAL KINETIC_ENERGY_AVAILABLE
	LOGICAL VERBOSE
	LOGICAL INCL_RAD_DECAYS
!
	LUER=ERROR_LU()
	DO_RAD_DECAYS=INCL_RAD_DECAYS
	IF(.NOT. DO_RAD_DECAYS)RETURN
	CALL GET_VERBOSE_INFO(VERBOSE)
!
	OPEN(UNIT=LU,FILE='NUC_DECAY_DATA',STATUS='OLD',ACTION='READ')
	STRING=' '
	DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	  READ(LU,'(A)')STRING
	END DO
	WRITE(LUER,*)'Accessing isotope data from NUC_DATA_FILE'
!
	LUER=ERROR_LU()
	NUM_DECAY_PATHS=0
	NUM_ISOTOPES=0
	IF(INDEX(STRING,'06-Jan-2015') .NE. 0)THEN
	  STRING='22-Sep-2006'
	  KINETIC_ENERGY_AVAILABLE=.TRUE.
	ELSE
	  KINETIC_ENERGY_AVAILABLE=.FALSE.
	  WRITE(6,*)' '
	  WRITE(6,*)'Warning -- NUC_DECAY_DATA has the old format'
	  WRITE(6,*)'This may affect results if usig ABS_TRANS option for gamma-ray transport'
	  WRITE(6,*)' '
	END IF
!
	IF(INDEX(STRING,'22-Sep-2006') .NE. 0)THEN
	  DO WHILE(NUM_DECAY_PATHS .EQ. 0)
	    READ(LU,'(A)')STRING
	    IF(INDEX(STRING,'!Date') .NE. 0)THEN
	       STRING=ADJUSTL(STRING); J=INDEX(STRING,'!Date')
	       WRITE(LUER,*)'Date associated with nuclear data file is',TRIM(STRING(1:J-1))
	    ELSE IF(INDEX(STRING,'Number of species') .NE. 0)THEN
	    ELSE IF(INDEX(STRING,'!Total number of isotopes') .NE. 0)THEN
	      READ(STRING,*)NUM_ISOTOPES
	      IF(NUM_ISOTOPES .GT. MAX_NUM_ISOTOPES)THEN
	        WRITE(LUER,*)'Error in RD_NUCLEAR_DECAY_DATA - MAX_NUM_ISOTOPES in NUC_IS_MOD too small'
	        WRITE(LUER,*)'Number of isotopes to be read is:',NUM_ISOTOPES
	        WRITE(LUER,*)'Maximum number of isotopes set in NUC_ISO_MOD is:',MAX_NUM_ISOTOPES
	        STOP
	      END IF
	    ELSE IF(INDEX(STRING,'!Maximum number of isotopes/species') .NE. 0)THEN
	    ELSE IF(INDEX(STRING,'!Number of reactions') .NE. 0)THEN
	      READ(STRING,*)NUM_DECAY_PATHS
	      IF(NUM_DECAY_PATHS.GT. MAX_NUM_DECAY_PATHS)THEN
	        WRITE(LUER,*)'Error in RD_NUCLEAR_DECAY_DATA - MAX_NUM_DECAY_PATHS in NUC_IS_MOD too small'
	        WRITE(LUER,*)'Number of decay paths (reactions) to be read is:',NUM_DECAY_PATHS
	        WRITE(LUER,*)'Maximum number of decay paths in NUC_ISO_MOD is:',MAX_NUM_DECAY_PATHS
	        STOP
	      END IF
	    ELSE IF(STRING(1:1) .EQ. ' ' .OR. STRING(1:1) .EQ. '!')THEN
	    ELSE
	      WRITE(LUER,*)'Error in RD_NUCLEAR_DECAY_DATA'
	      WRITE(LUER,*)'Unrecognized statment -- statement is:'
	      WRITE(LUER,*)TRIM(STRING)
	      STOP
	    END IF
	  END DO
	ELSE
	  WRITE(LUER,*)'Error in RD_NUCLEAR_DECAY_DATA'
	  WRITE(LUER,*)'Unrecognized date or date statement out of order'
	  STOP
	END IF
	IF(NUM_ISOTOPES .EQ. 0)THEN
	  WRITE(LUER,*)'Error in RD_NUCLEAR_DECAY_DATA'
	  WRITE(LUER,*)'The number of isotopes was not read in'
	  STOP
	END IF
!
	DO IS=1,NUM_ISOTOPES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LU,'(A)')STRING
	  END DO
	  J=INDEX(STRING,'  ')
	  ISO(IS)%SPECIES=STRING(1:J-1)
	  STRING=STRING(J+1:); STRING=ADJUSTL(STRING)
	  READ(STRING,*)ISO(IS)%MASS
	  ISO(IS)%BARYON_NUMBER=NINT(ISO(IS)%MASS)
	  J=INDEX(STRING,'  '); STRING=ADJUSTL(STRING(J+1:))
	  IF(STRING(1:1) .EQ. 's')THEN
	    ISO(IS)%STABLE=.TRUE.
	  ELSE IF(STRING(1:1) .EQ. 'u')THEN
	    ISO(IS)%STABLE=.FALSE.
	  ELSE
	    WRITE(LUER,*)'Error in RD_NUCLEAR_DECAY_DATA'
	    WRITE(LUER,*)'Stability identifier not recognized'
	    WRITE(LUER,*)'Identifier=',STRING(1:1)
	    STOP
	  END IF
	  IF(VERBOSE)THEN
	    WRITE(LUER,*)IS, ISO(IS)%SPECIES, ISO(IS)%MASS, ISO(IS)%BARYON_NUMBER
	  END IF
	END DO
	WRITE(LUER,*)'Successfully read in isotope data from NUC_DATA_FILE'
!
! Determine the link between ISOTOPE and SPECIES (as specified in MOD_CMFGEN).
!
	UNRECOGNIZED_ISO=.FALSE.
	DO IS=1,NUM_ISOTOPES
	  DO J=1,NUM_SPECIES
	    IF(ISO(IS)%SPECIES .EQ. SPECIES(J))THEN
	      ISO(IS)%ISPEC=J
	      EXIT
	    END IF
	  END DO
	  IF(ISO(IS)%ISPEC .EQ. 0)UNRECOGNIZED_ISO=.TRUE.
	END DO
!
	IF(UNRECOGNIZED_ISO)THEN
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)'Error in RD_NUCLEAR_DECAY_DATA: unrecognized isotopes'
	  WRITE(LUER,*)'Add the ISOTOPE to the SN_HYDRO_FILE (preferred)'
	  WRITE(LUER,*)'or delete ISOTOPE (and reactions?) from NUC_DECAY_DATA'
	  DO IS=1,NUM_ISOTOPES
	    IF(ISO(IS)%ISPEC .EQ. 0)WRITE(LUER,*)ISO(IS)%SPECIES
	  END DO
	  STOP
	END IF
!
	IN=0
	DO IN_LOOP=1,NUM_DECAY_PATHS
	  IN=IN+1
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LU,'(A)')STRING
	  END DO
!
	  STRING=ADJUSTL(STRING); L=INDEX(STRING,' ')
	  NUC(IN)%SPECIES=STRING(1:L-1)
	  STRING=STRING(L+1:)
!
	  STRING=ADJUSTL(STRING); L=INDEX(STRING,' ')
	  READ(STRING,*)NUC(IN)%MASS
	  NUC(IN)%BARYON_NUMBER=NINT(NUC(IN)%MASS)
	  STRING=STRING(L+1:)
!
	  STRING=ADJUSTL(STRING); L=INDEX(STRING,' ')
	  READ(STRING,*)NUC(IN)%HALF_LIFE
	  NUC(IN)%HALF_LIFE=NUC(IN)%HALF_LIFE*24.0_LDP*3600.0_LDP		!Convert to seconds
	  STRING=STRING(L+1:)
	  T1=2.0_LDP
	  NUC(IN)%DECAY_CONST=LOG(T1)/NUC(IN)%HALF_LIFE			!s^-1
!
	  STRING=ADJUSTL(STRING); L=INDEX(STRING,' ')
	  NUC(IN)%DAUGHTER=STRING(1:L-1)
	  STRING=STRING(L+1:)
!
	  STRING=ADJUSTL(STRING); L=INDEX(STRING,' ')
	  READ(STRING,*)NUC(IN)%DAUGHTER_MASS
	  NUC(IN)%DAUGHTER_BARYON_NUMBER=NINT(NUC(IN)%DAUGHTER_MASS)
	  STRING=STRING(L+1:)
!
	  STRING=ADJUSTL(STRING); L=INDEX(STRING,' ')
	  READ(STRING,*)NUC(IN)%ENERGY_PER_DECAY
	  NUC(IN)%ENERGY_PER_DECAY=1.0E+06_LDP*ELECTRON_VOLT()*NUC(IN)%ENERGY_PER_DECAY	!now in ergs
!
	  IF(KINETIC_ENERGY_AVAILABLE)THEN
	    STRING=STRING(L+1:)
	    STRING=ADJUSTL(STRING); L=INDEX(STRING,' ')
	    READ(STRING,*)NUC(IN)%KINETIC_PER_DECAY
	    NUC(IN)%KINETIC_PER_DECAY=1.0E+06_LDP*ELECTRON_VOLT()*NUC(IN)%KINETIC_PER_DECAY	!now in ergs
	  ELSE
	    NUC(IN)%KINETIC_PER_DECAY=0.0_LDP
	  END IF
!
	  STRING=STRING(L+1:)
	  STRING=ADJUSTL(STRING)
	  NUC(IN)%SEQUENCE=STRING(1:1)
	  IF(NUC(IN)%SEQUENCE .NE. 'F' .AND.  NUC(IN)%SEQUENCE .NE. 'S' .AND. NUC(IN)%SEQUENCE .NE. 'E')THEN
	    WRITE(LUER,*)'Unrecognized sequence string in RD_NUC_DECAY_DATA'
	    WRITE(LUER,*)'Sequence key read in is: ',NUC(IN)%SEQUENCE
	    WRITE(LUER,*)'STRING is: ',TRIM(STRING)
	    STOP
	  END IF
!
! Determine links between NUCLEAR reactant and the corresponding ISOTOPE.
! We need to do both the parent and the isotope.
!
	  DO IS=1,NUM_ISOTOPES
	    IF(NUC(IN)%SPECIES .EQ. ISO(IS)%SPECIES)THEN
	       NUC(IN)%ISPEC=ISO(IS)%ISPEC
	      IF(NUC(IN)%BARYON_NUMBER .EQ. ISO(IS)%BARYON_NUMBER)THEN
	        NUC(IN)%LNK_TO_ISO=IS
	        EXIT
	      END IF
	    END IF
	  END DO
!
	  DO IS=1,NUM_ISOTOPES
	    IF(NUC(IN)%DAUGHTER .EQ. ISO(IS)%SPECIES)THEN
	       NUC(IN)%DAUGHTER_ISPEC=ISO(IS)%ISPEC
	       IF(NUC(IN)%DAUGHTER_BARYON_NUMBER .EQ. ISO(IS)%BARYON_NUMBER)THEN
	         NUC(IN)%DAUGHTER_LNK_TO_ISO=IS
	         EXIT
	      END IF
	    END IF
	  END DO
!
	  IF(NUC(IN)%LNK_TO_ISO .EQ. 0 .OR. NUC(IN)%DAUGHTER_LNK_TO_ISO .EQ. 0)THEN
	    NUC(IN)%LNK_TO_ISO=0
	    NUC(IN)%DAUGHTER_LNK_TO_ISO=0
	    WRITE(6,*)'Warning: No match was found for the following nuclear reaction'
	    WRITE(6,'(2(A5,I6))')NUC(IN)%SPECIES,NUC(IN)%BARYON_NUMBER,NUC(IN)%DAUGHTER,NUC(IN)%DAUGHTER_BARYON_NUMBER
	    IN=IN-1
	  END IF
!
	END DO
	CLOSE(LU)
	NUM_DECAY_PATHS=IN
	WRITE(6,'(A,/)')' Successfully read in NUCLEAR data'
!
	IF(VERBOSE)THEN
	  DO IN=1,NUM_DECAY_PATHS
	    WRITE(6,'(I5,2A10,3I8)')IN,TRIM(NUC(IN)%SPECIES),TRIM(NUC(IN)%DAUGHTER),NUC(IN)%BARYON_NUMBER,
	1                  NUC(IN)%LNK_TO_ISO,NUC(IN)%DAUGHTER_LNK_TO_ISO
	  END DO
	END IF
!
! Determine link bewteen ISOTOPES (in ISO) and the PARENT group (in PAR)
! This determines the storage order in PAR.
!
	NUM_PARENTS=1
	PAR(1)%ISPEC=ISO(1)%ISPEC
	ISO(1)%LNK_TO_PAR=NUM_PARENTS
	DO IS=2,NUM_ISOTOPES
	  IF(PAR(NUM_PARENTS)%ISPEC .NE. ISO(IS)%ISPEC)THEN
	    NUM_PARENTS=NUM_PARENTS+1
	    PAR(NUM_PARENTS)%ISPEC=ISO(IS)%ISPEC
	  END IF
	  ISO(IS)%LNK_TO_PAR=NUM_PARENTS
	END DO
!
! Allocate all vectors
!
	DO IP=1,NUM_PARENTS
	  ALLOCATE (PAR(IP)%OLD_POP(ND));	   PAR(IP)%OLD_POP(:)=0.0_LDP
	  ALLOCATE (PAR(IP)%OLD_POP_DECAY(ND));    PAR(IP)%OLD_POP_DECAY(:)=0.0_LDP
	END DO
!
	DO IS=1,NUM_ISOTOPES
	  ALLOCATE (ISO(IS)%POP(ND));              ISO(IS)%POP=0.0_LDP
	  ALLOCATE (ISO(IS)%OLD_POP(ND));          ISO(IS)%OLD_POP=0.0_LDP
	  ALLOCATE (ISO(IS)%OLD_POP_DECAY(ND));    ISO(IS)%OLD_POP_DECAY=0.0_LDP
	END DO
!
	ALLOCATE (RADIOACTIVE_DECAY_ENERGY(ND));   RADIOACTIVE_DECAY_ENERGY=0.0_LDP
	ALLOCATE (KINETIC_DECAY_ENERGY(ND));       KINETIC_DECAY_ENERGY=0.0_LDP
!
	RETURN
	END
