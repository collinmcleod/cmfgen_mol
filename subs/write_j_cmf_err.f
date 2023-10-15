	SUBROUTINE WRITE_J_CMF_ERR(ITERATION_NO)
	IMPLICIT NONE
C
	INTEGER I,ITERATION_NO
	INTEGER N_ERR_MAX,FG_ERR_CNT,MOM_ERR_CNT,FG_ERR_TYPE
	PARAMETER (N_ERR_MAX=1000)
	REAL(10) FG_ERR_ON_FREQ,MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	COMMON /FG_J_CMF_ERR/FG_ERR_ON_FREQ(N_ERR_MAX),
	1                    FG_ERR_TYPE(N_ERR_MAX),FG_ERR_CNT
C
	IF(MOM_ERR_CNT .NE. 0 .OR. FG_ERR_CNT .NE. 0)
	1       WRITE(47,'(//,1X,A,I3)')'Iteration No',ITERATION_NO
C
	IF(MOM_ERR_CNT .NE. 0)THEN
	  WRITE(47,'(/,A)')' Error in MOM_J_CMF for the following frequencies:'
	  WRITE(47,'(A,I5)')' Number of frequencies for which error occurs: ',MOM_ERR_CNT
	  WRITE(47,'(1X,5ES14.6)')(MOM_ERR_ON_FREQ(I),I=1,MIN(MOM_ERR_CNT,N_ERR_MAX))
	END IF
C
	IF(FG_ERR_CNT .NE. 0)THEN
	  WRITE(47,'(/,A)')' Error in FG_J_CMF for the following frequencies:'
	  WRITE(47,'(1X,4(ES14.6,A,I2,A))')
	1      (FG_ERR_ON_FREQ(I),'(',FG_ERR_TYPE(I),')',I=1,FG_ERR_CNT)
	END IF
C
	RETURN
	END
