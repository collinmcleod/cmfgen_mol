	SUBROUTINE SET_SN_R_GRID(R,OLD_R,OLD_TAU,IB_RAT,OB_RAT,N_IB_INS,N_OB_INS,ND,NS)
	IMPLICIT NONE
!
	INTEGER NS
	INTEGER ND
!
	REAL*8 R(ND)
	REAL*8 OLD_R(NS)
	REAL*8 OLD_TAU(NS)
!
	REAL*8 ZN(2*ND)
	REAL*8 XN(ND)
	REAL*8 LOG_R(2*ND)
	REAL*8 LOG_OLD_R(NS)
	REAL*8 LOG_OLD_TAU(NS)
!
	REAL*8 IB_RAT
	REAL*8 OB_RAT
!
	INTEGER N_IB_INS
	INTEGER N_OB_INS
!
	REAL*8 dTAU
	REAL*8 dLOGR
	REAL*8 LOG_TAU_MIN
	REAL*8 LOG_R_MAX
	REAL*8 TAU_BEG,TAU_END
	REAL*8 T1,T2
	REAL*8 NEXT_R
	REAL*8 OB_RAT_LOC
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LUER=6
	INTEGER I,J,J_SAV
	INTEGER ND_TMP
!
! Estimate spacing to get required grid spacing.
!
	LOG_OLD_TAU=LOG(OLD_TAU)
	LOG_OLD_R=LOG(OLD_R)
!
	dTAU=(LOG_OLD_TAU(NS)-LOG_OLD_TAU(1))/(ND-1)
	OB_RAT_LOC=MIN(OB_RAT,EXP(dTAU))
!
	J=ND-1-N_OB_INS-N_IB_INS
	dTAU=(LOG_OLD_TAU(NS)-LOG_OLD_TAU(1)-N_OB_INS*LOG(OB_RAT_LOC))/J
	WRITE(6,*)LOG_OLD_TAU(ND),LOG_OLD_TAU(1)
	WRITE(6,*)'dTAU=',dTAU
	WRITE(6,*)'OB_RAT_LOC=',OB_RAT_LOC
!
	LOG_TAU_MIN=LOG_OLD_TAU(1)+N_OB_INS*LOG(OB_RAT_LOC)
	I=1
	DO WHILE(LOG_TAU_MIN .GT. LOG_OLD_TAU(I+1))
	  I=I+1
	END DO
	T1=(LOG_TAU_MIN-LOG_OLD_TAU(I))/(LOG_OLD_TAU(I+1)-LOG_OLD_TAU(I))
	LOG_R_MAX=T1*LOG_OLD_R(I+1)+(1.0D0-T1)*LOG_OLD_R(I)
	dLOGR=(LOG_R_MAX-LOG_OLD_R(NS))/J
	WRITE(6,*)'dLOGR=',dLOGR
	WRITE(6,*)'LOG_TAU_MIN=',LOG_TAU_MIN
	WRITE(6,*)'R_MAX=',EXP(LOG_R_MAX)
!
! Define the new radius grid. The step size in R corresponds to the smaller of
! dLOGR and dLOG_TAU. We create a "uniform" grid. The finer grid at the inner and
! outer boudaries is now only set when we set the FINAL grid.
!
	J=1; I=1
	LOG_R(1)=LOG_R_MAX
	DO WHILE(1 .EQ. 1)
	  I=I+1
	  IF(I .GT. 2*ND)THEN
	    WRITE(LUER,*)'Error in SET_SN_R_GRID --- LOG_R and TAU vectors too small'
	    WRITE(LUER,*)'I=',I,'J=',J
	    WRITE(LUER,*)'Log R(I)=',LOG_R(I)
	    WRITE(LUER,*)'Log old R(I)=',LOG_OLD_R(I)
	    WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL_V3 --- LOG_R and TAU vectors too small'
	    STOP
	  END IF
!
	  DO WHILE(LOG_OLD_R(J+1) .GT. LOG_R(I-1))
	    J=J+1
	  END DO
	  T1=(LOG_R(I-1)-LOG_OLD_R(J+1))/(LOG_OLD_R(J)-LOG_OLD_R(J+1))
	  TAU_BEG=T1*LOG_OLD_TAU(J)+(1.0D0-T1)*LOG_OLD_TAU(J+1)
!
! Compute dTAU step, and check to see if we need to use that.
!
	  J_SAV=J
	  NEXT_R=LOG_R(I-1)-dLOGR
	  DO WHILE(LOG_OLD_R(J+1) .GT. NEXT_R)
	    J=J+1
	  END DO
	  T1=(NEXT_R-LOG_OLD_R(J+1))/(LOG_OLD_R(J)-LOG_OLD_R(J+1))
	  TAU_END=T1*LOG_OLD_TAU(J)+(1.0D0-T1)*LOG_OLD_TAU(J+1)
!
	  IF(TAU_END-TAU_BEG .GT. dTAU)THEN
	    NEXT_R=LOG_R(I-1)-dLOGR*dTAU/(TAU_END-TAU_BEG)
	  END IF
	  J=J_SAV
	  LOG_R(I)=NEXT_R
!
! Check whether close enough to inner bondary.
!
	  IF(LOG_R(I)-1.5D0*(LOG_R(I-1)-LOG_R(I)) .LT. LOG_OLD_R(NS))EXIT
	END DO
	T1=LOG_R(I-1)-LOG_OLD_R(NS)
	LOG_R(I)=LOG_R(I-1)-0.5D0*T1
	I=I+1
	LOG_R(I)=LOG_OLD_R(NS)
	ND_TMP=I
!
! We now rescale the grid to have the correct number of grid points.
!
	J=ND-N_IB_INS-N_OB_INS
	DO I=1,ND_TMP; ZN(I)=I; END DO
	DO I=1,J
	  XN(I)=1.0D0+(I-1.0D0)*(ND_TMP-1.0D0)/(J-1.0D0)
	END DO
	CALL MON_INTERP(R,J,IONE,XN,J,LOG_R,ND_TMP,ZN,ND_TMP)
	LOG_R=R
	ND_TMP=J
	WRITE(6,*)EXP(R(1))
!
! Add extra points at inner boundary.
!
	I=ND_TMP
	T1=LOG_R(ND_TMP-1)-LOG_R(ND_TMP)
	IF(N_IB_INS .EQ. 1)THEN
	  LOG_R(I)=LOG_OLD_R(NS)+0.2D0*T1
	  ND_TMP=I+1
	ELSE IF(N_IB_INS .EQ. 2)THEN
	  LOG_R(I+1)=LOG_OLD_R(NS)+0.1D0*T1     !0.1D0
	  LOG_R(I)=LOG_OLD_R(NS)+0.4D0*T1      !0.4D0
	  ND_TMP=I+2
	ELSE IF(N_IB_INS .EQ. 3)THEN
	  LOG_R(I+2)=LOG_OLD_R(NS)+0.06D0*T1
	  LOG_R(I+1)=LOG_OLD_R(NS)+0.16D0*T1
	  LOG_R(I)=LOG_OLD_R(NS)+0.4D0*T1
	  ND_TMP=I+3
	ELSE
	  T2=1.0
	  DO J=1,N_IB_INS
	    T2=IB_RAT*T2+1.0D0
	  END DO
	  T1=T1/T2
	  DO J=N_IB_INS,1,-1
	    LOG_R(I+J-1)=LOG_OLD_R(NS)+T1
	    T1=T1*IB_RAT
	  END DO
	  ND_TMP=I+N_IB_INS
	END IF
	LOG_R(ND_TMP)=LOG_OLD_R(NS)
!
! Add finer grid at outer boundary.
!
! Shift grid to allow for insertion of extra ponts
!
	WRITE(6,*)EXP(LOG_R(1))
	DO I=ND_TMP,1,-1
	  LOG_R(I+N_OB_INS+1)=LOG_R(I)
	END DO
	WRITE(6,*)EXP(LOG_R(5))
!
	T1=(LOG_OLD_R(1)-LOG_R_MAX)/(N_OB_INS)
	WRITE(6,*)T1
	DO I=1,N_OB_INS-1
	   LOG_R(N_OB_INS+2-I)=LOG_R(N_OB_INS+2)+I*T1
	END DO
	LOG_R(1)=LOG_OLD_R(1)
	LOG_R(2)=LOG_R(1)-0.05D0*T1
!
	R=EXP(LOG_R)
	R(1)=OLD_R(1)
	R(ND)=OLD_R(NS)
	WRITE(LUER,*)'Computed R grid in SET_SN_R_GRID'
!
	RETURN
	END
