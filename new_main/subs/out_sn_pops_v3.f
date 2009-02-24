!
! Designed to ouptut SN data for the next model in a time
! dependent SN sequence.
!
	SUBROUTINE OUT_SN_POPS_V3(FILENAME,SN_AGE_DAYS,USE_OLD_MF_OUTPUT,ND,LU)
	USE MOD_CMFGEN
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
! Altered 11-Feb-2009 : Use isotope data, when available, to compute mass fractions.
!                       Added USE_OLD_MF_OUTPUT to allow consistency checks with older
!                         models.
! Created 21-Oct-2007
!
	REAL*8 SN_AGE_DAYS
	LOGICAL USE_OLD_MF_OUTPUT
	INTEGER ND
	INTEGER LU
	CHARACTER(LEN=*)FILENAME
!
! Local variables
!
	REAL*8 TMP_VEC(ND)
	INTEGER I,IS
	INTEGER ICOUNT
	CHARACTER*120 TMP_STR
	LOGICAL DONE_ISO
!
	ICOUNT=0
	DO I=1,NUM_SPECIES
	  IF(SUM(POP_SPECIES(I,:)) .NE. 0.0D0)ICOUNT=ICOUNT+1
	END DO
!
	OPEN(UNIT=LU,FILE=FILENAME,STATUS='UNKNOWN',ACTION='WRITE')
	WRITE(LU,'(/,A,I5)')'Number of data points:        ',ND
	WRITE(LU,'(A,I5)')  'Number of mass fractions:     ',ICOUNT
	WRITE(LU,'(A,I5)')  'Number of isotopes:           ',NUM_ISOTOPES
	WRITE(LU,'(A,F13.7,/)')'Time(days) since explosion:   ',SN_AGE_DAYS
!
	CALL OUT_SN_VEC(R,ND,'Radius grid (10^10cm)',LU)
	CALL OUT_SN_VEC(V,ND,'Velocity (km/s)',LU)
	CALL OUT_SN_VEC(SIGMA,ND,'Sigma (dlnV/dlnr-1)',LU)
	CALL OUT_SN_VEC(T,ND,'Temperature (10^4 K)',LU)
	CALL OUT_SN_VEC(DENSITY,ND,'Density (gm/cm^3)',LU)
	CALL OUT_SN_VEC(POP_ATOM,ND,'Atom density (/cm^3)',LU)
	CALL OUT_SN_VEC(ED,ND,'Electron density (/cm^3)',LU)
	CALL OUT_SN_VEC(ROSS_MEAN,ND,'Rosseland mean opacty (10^{-10} cm^{-1})',LU)
	TMP_VEC=1.0D-10*ROSS_MEAN/DENSITY
	CALL OUT_SN_VEC(TMP_VEC,ND,'Kappa (cm^2/gm)',LU)
!
! We output mass-fractions to RD_SN_DATA. As species have different atomic
! masses, we sume the individual mass fractions when available.
!
	DO I=1,NUM_SPECIES
	  TMP_STR=TRIM(SPECIES(I))//' mass fraction'
	  IF(SUM(POP_SPECIES(:,I))  .NE. 0.0D0)THEN
	    TMP_VEC=0.0D0
	    DONE_ISO=.FALSE.
	    DO IS=1,NUM_ISOTOPES
	      IF(ISO(IS)%ISPEC .EQ. I)THEN
	        TMP_VEC=TMP_VEC+1.66D-24*ISO(IS)%POP*ISO(IS)%MASS/DENSITY
	        DONE_ISO=.TRUE.
	      END IF
	    END DO
	    IF(USE_OLD_MF_OUTPUT .OR. .NOT. DONE_ISO)THEN
	      TMP_VEC=1.66D-24*POP_SPECIES(:,I)*AT_MASS(I)/DENSITY
	    END IF
	    CALL OUT_SN_VEC(TMP_VEC,ND,TMP_STR,LU)
	  ELSE
	    WRITE(LU,'(/,A)')TRIM(TMP_STR)
	    WRITE(LU,'(2X,I5,A)')ND,'*0.00000D0'
	  END IF
	END DO
!
	DO I=1,NUM_ISOTOPES
	  WRITE(TMP_STR(1:3),'(I3)')ISO(I)%BARYON_NUMBER
	  TMP_STR=TRIM(ISO(I)%SPECIES)//TMP_STR(1:3)//' mass fraction'
	  TMP_VEC=1.66D-24*ISO(I)%POP*ISO(I)%MASS/DENSITY
	  CALL OUT_SN_VEC(TMP_VEC,ND,TMP_STR,LU)
	END DO
!
	CLOSE(LU)
	RETURN
	END
!
!^L
!
	SUBROUTINE OUT_SN_VEC(X,ND,HEADER,LU)
	IMPLICIT NONE
	INTEGER ND,LU
	REAL*8 X(ND)
	CHARACTER(LEN=*) HEADER
	INTEGER I
!
	WRITE(LU,'(/,A)')TRIM(HEADER)
	WRITE(LU,'(1X,1P8E16.7)')(X(I),I=1,ND)
!
	RETURN
	END
