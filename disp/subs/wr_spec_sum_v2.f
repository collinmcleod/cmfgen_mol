!
! Designed to output additional information (XV) to aid in
! reading mass-fractions plots.
!
	SUBROUTINE WR_SPEC_SUM_V2(NUM_FRAC,MASS_FRAC,XV,ND)
	USE MOD_DISP
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Created: 24_mar-2012
!          Based on WR_SPEC_CUM
!
	REAL(10) XV(ND)
	INTEGER ND
	LOGICAL NUM_FRAC
	LOGICAL MASS_FRAC
!
	REAL(10) T1
	REAL(10) ATOMIC_MASS_UNIT
	REAL(10) VALUE(NSPEC)
	INTEGER ISPEC_STORE(NSPEC)
	INTEGER IWORK(NSPEC)
	INTEGER INDX(NSPEC)
!
	EXTERNAL ATOMIC_MASS_UNIT
	INTEGER, PARAMETER :: NSTR=10
	CHARACTER(LEN=160) STRING(NSTR)
	CHARACTER(LEN=160) HEADER_A
	CHARACTER(LEN=160) HEADER_B
!
	INTEGER, PARAMETER :: NCOL=4
	INTEGER ID(NCOL)
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	INTEGER I,J,K
	INTEGER DPTH_INDX
	INTEGER ISPEC
!
	STRING=' '
	HEADER_A=' '
	HEADER_B=' '
	DPTH_INDX=1
!
	T1=NCOL-1.0D0
	DO K=1,NCOL
	  ID(K)=NINT(1+(K-1)*(ND-1)/T1)
	END DO
	IF(XV(ND) .LT. XV(1))THEN
	  DO K=1,NCOL/2
	    J=ID(K)
	    ID(K)=ID(NCOL+1-K)
	    ID(NCOL+1-K)=J
	  END DO
	END IF
	WRITE(6,*)ID(1:NCOL)
!
	DO K=1,NCOL
!
	  DPTH_INDX=ID(K)
	  J=LEN_TRIM(HEADER_A)+1
	  IF(K .GT. 1)J=J+11
	  WRITE(HEADER_A(J:),'(7X,A,I3,A,F8.1)')'V(',DPTH_INDX,')=',V(DPTH_INDX)
	  IF(ABS(XV(DPTH_INDX)) .GT. 0.1 .AND. ABS(XV(DPTH_INDX))  .LT. 99.999D0)THEN
	    WRITE(HEADER_B(J:),'(7X,A,I3,A,F8.3)')'X(',DPTH_INDX,')=',XV(DPTH_INDX)
	  ELSE
	    WRITE(HEADER_B(J:),'(7X,A,I3,A,ES8.2)')'X(',DPTH_INDX,')=',XV(DPTH_INDX)
	  END IF
	
	  VALUE=0.0D0
	  DO ISPEC=1,NSPEC
	    IF(POPDUM(DPTH_INDX,ISPEC) .GT. 0.0D0)THEN
	      IF(NUM_FRAC)THEN
	        VALUE(ISPEC)=POPDUM(DPTH_INDX,ISPEC)/POP_ATOM(DPTH_INDX)+1.0D-100
	      ELSE IF(MASS_FRAC)THEN
	        T1=AT_MASS(ISPEC)*ATOMIC_MASS_UNIT()
	        VALUE(ISPEC)=T1*POPDUM(DPTH_INDX,ISPEC)/MASS_DENSITY(DPTH_INDX)+1.0D-100
	      ELSE
	        VALUE(ISPEC)=POPDUM(DPTH_INDX,ISPEC)+1.0D-100
	      END IF
	    END IF
	    ISPEC_STORE(ISPEC)=ISPEC
	  END DO
	  CALL INDEXX(NSPEC,VALUE,INDX,L_FALSE)
	  CALL SORTINT(NSPEC,ISPEC_STORE,INDX,IWORK)
!
	  DO I=1,NSTR
	    IF(ISPEC_STORE(I) .EQ. 0)EXIT
	    J=LEN_TRIM(STRING(I))+1
	    IF(J .NE. 1)J=J+5
	    ISPEC=ISPEC_STORE(I)
	    T1=VALUE(ISPEC)
	    IF(NUM_FRAC .OR. MASS_FRAC)THEN
	      WRITE(STRING(I)(J:),'(A6,2X,F9.6,2X,ES9.2,4X)')TRIM(SPECIES(ISPEC)),T1,LOG10(T1)
	    ELSE
	      WRITE(STRING(I)(J:),'(A6,2X,ES9.3,3X,ES9.2,5X)')TRIM(SPECIES(ISPEC)),T1,LOG10(T1)
	    END IF
	    END DO
	END DO
!
        WRITE(6,'(A)')' '
	IF(NUM_FRAC)THEN
	  WRITE(6,'(A)')BLUE_PEN//' Number fraction (NF) and Log(NF) at inner boundary'//DEF_PEN
	ELSE IF(MASS_FRAC)THEN
	  WRITE(6,'(A)')BLUE_PEN//' Mass fraction (MF) and Log(MF) at inner boundary'//DEF_PEN
	ELSE
	  WRITE(6,'(A)')BLUE_PEN//' Number denisty (ND) and Log(ND) at inner boundary'//DEF_PEN
	END IF
	WRITE(6,'(A)')' '
!
	WRITE(6,'(A)')TRIM(HEADER_A)
	WRITE(6,'(A)')TRIM(HEADER_B)
	DO I=1,NSTR
	  IF(STRING(I) .NE. ' ')THEN
	    WRITE(6,'(A)')TRIM(STRING(I))
	  END IF
	END DO
!
	RETURN
	END
