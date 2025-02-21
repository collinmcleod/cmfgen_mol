	SUBROUTINE DETERM_REDDENING(OBS_SPEC,OBS_NU,NOBS,
	1                  FIT_SPEC,FIT_NU,NFIT,
	1                  R_MIN,R_MAX)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NOBS,NFIT
	REAL(KIND=LDP) OBS_SPEC(NOBS),OBS_NU(NOBS)
	REAL(KIND=LDP) FIT_SPEC(NFIT),FIT_NU(NFIT)
	REAL(KIND=LDP) R_MIN,R_MAX
!
	REAL(KIND=LDP) RHS(2)
	REAL(KIND=LDP) MAT(2,2)
!
	REAL(KIND=LDP), ALLOCATABLE :: NU_ST(:)
	REAL(KIND=LDP), ALLOCATABLE :: NU_END(:)
	REAL(KIND=LDP), ALLOCATABLE :: LAM_ST(:)
	REAL(KIND=LDP), ALLOCATABLE :: LAM_END(:)
	REAL(KIND=LDP), ALLOCATABLE :: OBS(:)
	REAL(KIND=LDP), ALLOCATABLE :: FIT(:)
	REAL(KIND=LDP), ALLOCATABLE :: WEIGHT(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAL_RED_LAW(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAL_RED_LAW_SAVE(:)
!
	INTEGER IOS
	INTEGER I,K,L
	INTEGER LST,LEND
	REAL(KIND=LDP) R_EXT
	REAL(KIND=LDP) RAX,RBX
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) EBMV
	REAL(KIND=LDP) LOG_DIST
	REAL(KIND=LDP) CHI_SQ
	REAL(KIND=LDP) ANG_TO_HZ
!
	REAL(KIND=LDP) CHI_SQ_SVE
	REAL(KIND=LDP) EBMV_SVE
	REAL(KIND=LDP) DIST_SVE
	REAL(KIND=LDP) R_SVE
!
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER NUM_BNDS
	INTEGER GET_INDX_DP
	EXTERNAL GET_INDX_DP
!
	CHI_SQ_SVE=1.0E+200_LDP
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0E-07_LDP      !10^8/10^15
!
	OPEN(UNIT=10,FILE='RED_BANDS',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error opening RED_BANDS -- must be in same directory as issuing PLT_SPEC command'
	    RETURN
	  END IF
	  READ(10,*)NUM_BNDS
	  ALLOCATE (LAM_ST(NUM_BNDS))
	  ALLOCATE (LAM_END(NUM_BNDS))
	  ALLOCATE (WEIGHT(NUM_BNDS))
	  DO I=1,NUM_BNDS
	    READ(10,*)LAM_ST(I),LAM_END(I),WEIGHT(I)
	  END DO
	CLOSE(UNIT=10)
!
	ALLOCATE (NU_ST(NUM_BNDS))
	ALLOCATE (NU_END(NUM_BNDS))
	ALLOCATE (OBS(NUM_BNDS))
	ALLOCATE (FIT(NUM_BNDS))
	ALLOCATE (GAL_RED_LAW(NUM_BNDS))
	ALLOCATE (GAL_RED_LAW_SAVE(NUM_BNDS))
!
	DO I=1,NUM_BNDS
	  NU_ST(I)=ANG_TO_HZ/LAM_ST(I)
	  NU_END(I)=ANG_TO_HZ/LAM_END(I)
	END DO
!
	DO I=1,NUM_BNDS
	  OBS(I)=0.0_LDP
	  LST=GET_INDX_DP(NU_ST(I),OBS_NU,NOBS)
	  LEND=GET_INDX_DP(NU_END(I),OBS_NU,NOBS)
	  DO L=MIN(LST,LEND),MAX(LST,LEND)
	    OBS(I)=OBS(I)+OBS_SPEC(L)
	  END DO
	  OBS(I)=OBS(I)/(ABS(LST-LEND)+1)
	END DO
!
	DO I=1,NUM_BNDS
	  FIT(I)=0.0_LDP
	  LST=GET_INDX_DP(NU_ST(I),FIT_NU,NFIT)
	  LEND=GET_INDX_DP(NU_END(I),FIT_NU,NFIT)
	  DO L=MIN(LST,LEND),MAX(LST,LEND)
	    FIT(I)=FIT(I)+FIT_SPEC(L)
	  END DO
	  FIT(I)=FIT(I)/(ABS(LST-LEND)+1)
	END DO
!
	WRITE(6,*)' '
	WRITE(6,'(A,4(6X,A))')'  I',' Lam(st)','Lam(end)','  OBS(I)','  MOD(I)'
	DO I=1,NUM_BNDS
	  WRITE(6,'(I3,4ES14.4)')I,LAM_ST(I),LAM_END(I),OBS(I),FIT(I)
	END DO
!
	WRITE(6,'(A)')
	WRITE(6,'(5(8X,A))')'     R','     d','E(B-V)',' Chi^2','dCHI^2'
!
! Peform the minization. We loop over R_EXT.
!
	DO L=1,NINT(1+(R_MAX-R_MIN)/0.05_LDP)
	  R_EXT=R_MIN+(L-1)*0.05_LDP
!
! Get Galactic reddenining curve.
!
	  DO I=1,NUM_BNDS
	    T1=(LAM_ST(I)+LAM_END(I))/2.0_LDP
	    T1=10000.0_LDP/T1                    !1/Lambda(um)
	   IF(T1 .LT. 1.1_LDP)THEN
	      RAX=0.574_LDP*(T1**1.61_LDP)
	      RBX=-0.527_LDP*(T1**1.61_LDP)
	    ELSE IF(T1. LT. 3.3_LDP)THEN
	      T2=T1-1.82_LDP
	      RAX=1+T2*(0.17699_LDP-T2*(0.50447_LDP+T2*(0.02427_LDP-T2*(0.72085_LDP
	1                 +T2*(0.01979_LDP-T2*(0.77530_LDP-0.32999_LDP*T2))))))
	      RBX=T2*(1.41338_LDP+T2*(2.28305_LDP+T2*(1.07233_LDP-T2*(5.38434_LDP
	1                +T2*(0.62251_LDP-T2*(5.30260_LDP-2.09002_LDP*T2))))))
	    ELSE IF(T1 .lT. 5.9_LDP)THEN
	      RAX=1.752_LDP-0.316_LDP*T1-0.104_LDP/((T1-4.67_LDP)**2+0.341_LDP)
	      RBX=-3.090_LDP+1.825_LDP*T1+1.206_LDP/((T1-4.62_LDP)**2+0.263_LDP)
	    ELSE IF(T1 .LT. 8.0_LDP)THEN
  	      T2=T1-5.9
	      RAX=1.752_LDP-0.316_LDP*T1-0.104_LDP/((T1-4.67_LDP)**2+0.341_LDP) -
	1                       T2*T2*(0.04773_LDP+0.009779_LDP*T2)
	      RBX=-3.090_LDP+1.825_LDP*T1+1.206_LDP/((T1-4.62_LDP)**2+0.263_LDP)+
	1                       T2*T2*(0.2130_LDP+0.1207_LDP*T2)
	    ELSE IF(T1 .LT. 10)THEN
	      T2=T1-8
	      RAX=-1.073_LDP-T2*(0.628_LDP-T2*(0.137_LDP-0.070_LDP*T2))
	      RBX=13.670_LDP+T2*(4.257_LDP-T2*(0.420_LDP-0.374_LDP*T2))
	    ELSE
	      T1=10
	      T2=T1-8
	      RAX=-1.073_LDP-T2*(0.628_LDP-T2*(0.137_LDP-0.070_LDP*T2))
	      RBX=13.670_LDP+T2*(4.257_LDP-T2*(0.420_LDP-0.374_LDP*T2))
	    END IF
            GAL_RED_LAW(I)=R_EXT*(RAX+RBX/R_EXT)
	  END DO
!
	  GAL_RED_LAW(1:NUM_BNDS)=0.921_LDP*GAL_RED_LAW(1:NUM_BNDS)
	  MAT=0.0_LDP
	  RHS=0.0_LDP
	  DO I=1,NUM_BNDS
	    RHS(1)=RHS(1)-WEIGHT(I)*2.0_LDP*LOG(OBS(I)/FIT(I))
	    RHS(2)=RHS(2)-WEIGHT(I)*LOG(OBS(I)/FIT(I))*GAL_RED_LAW(I)
	    MAT(1,1)=MAT(1,1)+WEIGHT(I)*4.0_LDP
	    MAT(1,2)=MAT(1,2)+WEIGHT(I)*2.0_LDP*GAL_RED_LAW(I)
	    MAT(2,1)=MAT(2,1)+WEIGHT(I)*2.0_LDP*GAL_RED_LAW(I)
	    MAT(2,2)=MAT(2,2)+WEIGHT(I)*GAL_RED_LAW(I)*GAL_RED_LAW(I)
	  END DO
!
! Determine parameters:
!                  Log d(kpc) and E(B-V)
!
	  I=2
	  CALL SIMQ(MAT,RHS,I,K)
!
! Get fit value.
!
	  LOG_DIST=RHS(1)
	  EBMV=RHS(2)
	  CHI_SQ=0.0_LDP
	  DO I=1,NUM_BNDS
	    T1=WEIGHT(I)*(LOG(OBS(I)/FIT(I))+2*LOG_DIST+GAL_RED_LAW(I)*EBMV)**2
	    CHI_SQ=CHI_SQ+T1
	  END DO
	  WRITE(6,'(5ES14.4)')R_EXT, EXP(LOG_DIST), EBMV, CHI_SQ, 1.0D+04*CHI_SQ/NUM_BNDS
	  IF(CHI_SQ .LT.  CHI_SQ_SVE)THEN
	    CHI_SQ_SVE=CHI_SQ
	    GAL_RED_LAW_SAVE=GAL_RED_LAW
	    EBMV_SVE=EBMV
	    DIST_SVE=EXP(LOG_DIST)
	    R_SVE=R_EXT
	  END IF
	END DO
!
	DO I=1,NUM_BNDS
	  T1=WEIGHT(I)*(LOG(OBS(I)/FIT(I))+2*LOG(DIST_SVE)+GAL_RED_LAW_SAVE(I)*EBMV_SVE)**2
	  WRITE(6,'(2F14.2,3X,F3.0,ES12.3)')LAM_ST(I),LAM_END(I),WEIGHT(I),T1/CHI_SQ_SVE
	END DO
!
	WRITE(6,'(A)')
	WRITE(6,'(A,F10.3)')'        The optimal E(B-V) is:',EBMV_SVE
	WRITE(6,'(A,F10.3)')'       The optimal R value is:',R_SVE
	WRITE(6,'(A,F10.3)')'The optimal distance (kpc) is:',DIST_SVE
	WRITE(6,'(A)')
!
	DEALLOCATE (LAM_ST,LAM_END,NU_ST,NU_END,OBS,FIT,WEIGHT,GAL_RED_LAW)
!
	RETURN
	END
