	PROGRAM MASS_SC
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
C
C Routine to read in a set of abundances (generally solar) from a file.
C
C Abundances of individual species can be adjusted. Species not adjusted
C are held fixed at a constant mass-fraction, and a revised X/He abundace
C output.
C
	INTEGER*4, PARAMETER :: MAX_EL=92
C
	REAL*8 AT_NO(MAX_EL)			!Atomic number
	REAL*8 AT_MASS(MAX_EL)			!Atomic mass(amu)
	REAL*8 ABUND(MAX_EL)			!Fractional abundace (relative)
	REAL*8 MASS_FRAC(MAX_EL)		!Mass fractional abundace.
	CHARACTER*2 SYMB(MAX_EL)		!Element symbol
	CHARACTER*20 NAME(MAX_EL)		!Element name
C
C
	REAL*8 NEW_ABUND(MAX_EL)
	REAL*8 NEW_MASS_FRAC(MAX_EL)
	LOGICAL*1 ALTERED_ABUND(MAX_EL)		!Indicates revised abundace
C
	CHARACTER*2 SPEC			!Used for IO
	LOGICAL*1 SYMB_OK
	CHARACTER*80 STRING
	REAL*4 VAL
C
	REAL*8 MASS
	REAL*8 OLD_MASS_SUM
	REAL*8 NEW_MASS_SUM
	REAL*8 NHE
	INTEGER*4 I,N
C
C Input abundance data from approproiately formated file. Strings must
C be enclosed in quotes to allow free format read.
C
C We skip comments by finding the first record with "Hydrogen" in it.
C
	OPEN(UNIT=20,FILE='SOL_ABUND',STATUS='OLD',ACTION='READ')
	 STRING=' '
	 DO WHILE(INDEX(STRING,'Hydrogen') .EQ. 0)
	   READ(20,'(A)')STRING
	 END DO
	 BACKSPACE(20)
C
	 I=0
	 DO WHILE(1 .EQ. 1)
	   READ(20,*,END=100)AT_NO(I+1),SYMB(I+1),
	1             NAME(I+1),AT_MASS(I+1),ABUND(I+1)
	   I=I+1
	 END DO
100	 CONTINUE
	 N=I
	CLOSE(UNIT=20)
!
! Decide on format for the abundances.
!
	IF(ABS(ABUND(1)-12.0D0) .LT. 0.1)THEN
	  DO I=1,N
	    ABUND(I)=10.0D0**(ABUND(I)-12.0D0)
	  END DO
	END IF
C
C Compute the mass fractions of all species.
C
	MASS=0.0D0
	DO I=1,N
	  MASS_FRAC(I)=AT_MASS(I)*ABUND(I)
	  MASS=MASS+MASS_FRAC(I)
	END DO
	MASS_FRAC(1:N)=MASS_FRAC(1:N)/MASS
C
C Output for check.
C
	DO I=1,N
	  WRITE(6,200)AT_NO(I),SYMB(I),NAME(I),AT_MASS(I),ABUND(I),
	1              MASS_FRAC(I)
	END DO
	WRITE(6,'(A)')' '
200	FORMAT(X,F4.1,3X,A2,3X,A10,4X,F5.1,3X,1PE8.2,3X,E8.2)
C
C We now begin the section to allow individual abundances to be varied.
C Each species can be separately changed by using it chenical abbreviation.
C
	NEW_ABUND(1:N)=0.0D0
	SPEC=' '
	DO WHILE(SPEC(1:1) .NE. 'E')
	  SPEC='E'
	  CALL GEN_IN(SPEC,'Chemical symbol')
	  SYMB_OK=.FALSE.
	  DO I=1,N
	    IF(SPEC(1:2) .EQ. SYMB(I))THEN
	      VAL=ABUND(I)
	      CALL GEN_IN(VAL,'Abund for '//NAME(I))
	      NEW_ABUND(I)=VAL
	      ALTERED_ABUND(I)=.TRUE.
	      SYMB_OK=.TRUE.
	    END IF
	  END DO
	  IF(.NOT. SYMB_OK .AND. SPEC(1:1) .NE. 'E')
	1        WRITE(6,*)'Invalid chemical symbol'
	END DO
C
C Bow adjust the new abundances so that their combined mass fraction is the 
C same as C in the solar data. As a consequence all other species will have 
C the same-mass fraction.
C
C Compute the scale factors.
C
	OLD_MASS_SUM=0.0D0
	NEW_MASS_SUM=0.0D0
	DO I=1,N
	  IF(ALTERED_ABUND(I))THEN
	    OLD_MASS_SUM=OLD_MASS_SUM+MASS_FRAC(I)
	    NEW_MASS_FRAC(I)=AT_MASS(I)*NEW_ABUND(I)
	    NEW_MASS_SUM=NEW_MASS_SUM+NEW_MASS_FRAC(I)
	  END IF
	END DO
C
C Do the actual scaling.
C
	MASS=0.0D0
	DO I=1,N
	  IF(ALTERED_ABUND(I))THEN
	    NEW_MASS_FRAC(I)=NEW_MASS_FRAC(I)*OLD_MASS_SUM/NEW_MASS_SUM
	  ELSE
	    NEW_MASS_FRAC(I)=MASS_FRAC(I)              
	  END IF
	  MASS=MASS+NEW_MASS_FRAC(I)
	END DO
	WRITE(6,*)MASS
C
C Copute the relative fractional populations.
C
	DO I=1,N
	  NEW_ABUND(I)=NEW_MASS_FRAC(I)*MASS/AT_MASS(I)
	  IF(SYMB(I) .EQ. 'He')NHE=NEW_ABUND(I)
	END DO
	IF(NHE .EQ. 0.0D0)NHE=MAXVAL(NEW_ABUND)
C
C Normalize the abundances so the the He abundace is 1.0
C
	NEW_ABUND(1:N)=NEW_ABUND(1:N)/NHE
C
C With the new fractional abundace we compute the revised mass-fractions. 
C Acts as a check.
C
	MASS=0.0D0
	DO I=1,N
	  NEW_MASS_FRAC(I)=AT_MASS(I)*NEW_ABUND(I)
	  MASS=MASS+NEW_MASS_FRAC(I)
	END DO
	NEW_MASS_FRAC(1:N)=NEW_MASS_FRAC(1:N)/MASS
C
	WRITE(6,'(34X,A,A)')'  N(old)     M(old)  ',
	1                   '   N(new)     M(new)  '
	DO I=1,N
	  WRITE(6,300)AT_NO(I),SYMB(I),NAME(I),AT_MASS(I),ABUND(I),
	1              MASS_FRAC(I),NEW_ABUND(I),NEW_MASS_FRAC(I)
	END DO
300	FORMAT(X,F4.1,3X,A2,3X,A10,4X,F5.1,3X,1PE8.2,3X,E8.2,3X,E8.2,
	1       3X,E8.2)
C
	OPEN(UNIT=10,FILE='REV_ABUND',STATUS='UNKNOWN')
	  WRITE(10,'(34X,A,A)')'  N(old)     M(old)  ',
	1                   '   N(new)     M(new)  '
	  DO I=1,N
	    WRITE(10,300)AT_NO(I),SYMB(I),NAME(I),AT_MASS(I),ABUND(I),
	1              MASS_FRAC(I),NEW_ABUND(I),NEW_MASS_FRAC(I)
	  END DO
	CLOSE(UNIT=10)
C
	STOP
	END
