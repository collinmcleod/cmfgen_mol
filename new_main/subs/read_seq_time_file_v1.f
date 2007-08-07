!
! This subroutine simply returns the number of depth points used for the OLD_MODEL.
! This subsequently allows to allocate the necessary storage for reading in the data.
!
! Created 14-Mar-2007.
!
	SUBROUTINE GET_ND_SEQ_MODEL_FILE(ND,LU)
	IMPLICIT NONE
	INTEGER ND
	INTEGER LU
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=11) DATE
!
	OPEN(UNIT=LU,FILE='OLD_MODEL_DATA',FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
	  READ(LU)DATE
	  IF(DATE .NE. '17-Mar-2007')THEN
	    WRITE(ERROR_LU(),*)'Invalid format date in READ_SEQ_TIME_V1'
	    WRITE(ERROR_LU(),*)'Date is ',DATE
	    STOP
	  END IF
	  READ(LU)ND
	CLOSE(UNIT=LU)
!
	RETURN
	END
!
! Read in sequential unformatted model file with old model data. The format of
! this file has been designed so that number of species, levels, and depths can be
! changed in a time sequence. There are still issues if we add a species.
! 
	SUBROUTINE READ_SEQ_TIME_FILE_V1(OLD_R,OLD_V,OLD_SIGMA,OLD_POP_ATOM,OLD_DENSITY,
	1              POPS,OLD_ION_STAGE_PRES,SN_AGE,ND,NT,LU)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 14-Mar-2007  Based on READ_TIME_MODEL_V2
!
	INTEGER NT
	INTEGER ND
	INTEGER LU
!
	REAL*8 OLD_R(ND)
	REAL*8 OLD_V(ND)
	REAL*8 OLD_SIGMA(ND)
	REAL*8 OLD_POP_ATOM(ND)
	REAL*8 OLD_DENSITY(ND)
	REAL*8 POPS(NT,ND)
	REAL*8 SN_AGE
!
	LOGICAL OLD_ION_STAGE_PRES(NUM_IONS)
!
! Local variables.
!
	REAL*8, ALLOCATABLE :: OLD_XzV(:,:)
	REAL*8, ALLOCATABLE :: TMP_XzV(:,:)
	REAL*8, ALLOCATABLE :: TMP_DXzV(:)
	REAL*8 HDKT
	REAL*8 T1
!
	INTEGER ID
	INTEGER ISPEC
	INTEGER NX
	INTEGER I
	INTEGER J
	INTEGER L
	INTEGER ID_BEG,ID_END
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=11) DATE
	CHARACTER(LEN=10) SPECIES_NAME
!
	LUER=ERROR_LU()
        HDKT=4.7994145D0					!1.0D+15*H/k/1.0D+04
	POPS=0.0D0
	OLD_ION_STAGE_PRES(1:NUM_IONS)=.FALSE.
!
	OPEN(UNIT=LU,FILE='OLD_MODEL_DATA',FORM='UNFORMATTED',STATUS='UNKNOWN',ACTION='READ')
!
	  READ(LU)DATE
	  IF(DATE .NE. '17-Mar-2007')THEN
	    WRITE(LUER,*)'Invalid format date in READ_SEQ_TIME_V1'
	    WRITE(LUER,*)'Date is ',DATE
	    STOP
	  END IF
	  READ(LU)I
	  IF(I .NE. ND)THEN
	    WRITE(LUER,*)'Number of depth points doesn''t match in READ_SEQ_TIME_V1'
	    STOP
	  END IF
!
	  READ(LU)SN_AGE
	  READ(LU)OLD_R
	  READ(LU)OLD_V
	  READ(LU)OLD_SIGMA
	  READ(LU)POPS(NT,:)		!Temperature
	  READ(LU)POPS(NT-1,:)		!Electron density
	  READ(LU)OLD_POP_ATOM
	  READ(LU)OLD_DENSITY
!
! Loop over all species, and all ionization stages, in the file.
!
	  ALLOCATE (TMP_DXzV(ND))
	  DO ISPEC=1,NUM_SPECIES
	    READ(LU)ID_BEG,ID_END,SPECIES_NAME
	    IF(SPECIES_NAME .NE. SPECIES(ISPEC))THEN
	      WRITE(LUER,*)'Error in READ_SEQ_TIME_FILE_V1: Species mismatch'
	      WRITE(LUER,*)'   Species name in file=',TRIM(SPECIES_NAME)
	      WRITE(LUER,*)'Species name in program=',TRIM(SPECIES(ISPEC))
              STOP
	    END IF
!
! NX is the number of levels in the old model atom.
! We use TMP_XzV for the populations in the old model as read from the file.
! We use OLD_XzV for the old populations, with the number of levels adjusted for 
!                                         the new model.
!
	    DO ID=ID_BEG,ID_END				!Ionization stages
	      READ(LU)I,NX
	      IF(I .NE. ID)THEN
	        WRITE(LUER,*)'Error in READ_SEQ_TIME_FILE_V1: ID mismatch'
	        WRITE(LUER,*)'      ID in file=',I
	        WRITE(LUER,*)'   ID in program=',ID
                STOP
	      END IF
	      ALLOCATE(TMP_XzV(NX,ND))
	      READ(LU)TMP_XzV
	      READ(LU)TMP_DXzV
!
	      IF(ID .GE. SPECIES_BEG_ID(ISPEC) .AND. ID .LE. SPECIES_END_ID(ISPEC)-1)THEN
!
! If level is not present, set departure coefficient equal to that of highest level.
! Ignores level dissolution.
!
	        ALLOCATE (OLD_XzV(ATM(ID)%NXzV_F,ND))
	        OLD_ION_STAGE_PRES(ID)=.TRUE.
	        J=MIN(NX,ATM(ID)%NXzV_F)
	        OLD_XzV(1:J,1:ND)=TMP_XzV(1:J,1:ND)
	        DO L=1,ND
	          DO I=J+1,ATM(ID)%NXzV_F
		     T1=HDKT*(ATM(ID)%EDGEXZV_F(J)-ATM(ID)%EDGEXZV_F(I))/T(L)
	             OLD_XzV(I,L)=OLD_XzV(J,L)*ATM(ID)%GXZV_F(I)/ATM(ID)%GXZV_F(J)*EXP(T1)
	          END DO
	        END DO
!
! Store population in the POPS array. POPS has been previously zeroed.
!
	     	DO L=1,ND
	          DO I=1,ATM(ID)%NXzV_F
	            J=ATM(ID)%EQXzV-1+ATM(ID)%F_TO_S_XzV(I)
	            POPS(J,L)=POPS(J,L)+OLD_XzV(I,L)
	          END DO
	        END DO
!
	        IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
	          J=ATM(ID)%EQXzV+ATM(ID)%NXzV		!Small NXzV
	          POPS(J,:)=TMP_DXzV
	        END IF
!
	        DEALLOCATE(OLD_XzV)
	      END IF
	      DEALLOCATE(TMP_XzV)
!
	    END DO		!Loop over ionization stage.
	  END DO		!Loop over species
	  DEALLOCATE(TMP_DXzV)
!
	CLOSE(LU)
!
	RETURN
	END
