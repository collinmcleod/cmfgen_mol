!
! Subroutine to read a HYDRO model from SN_HYDRO_DATA and genrate
! the R grid, and the V and SIGMA vectors.
!
	SUBROUTINE SET_RV_HYDRO_MODEL_V3(R,V,SIGMA,RMAX,RCORE,RMAX_ON_RCORE,
	1                  SN_AGE_DAYS,PURE_HUBBLE_FLOW,N_IB_INS,N_OB_INS,RDINR,ND,LU)
	IMPLICIT NONE
!
! Altered 10-Jul-2012 : Modified computation of R grid, especially the grid near the boundaries.
! Altered 28-Jan-2012 : Initialize dLOGR and dTAU=0.0D when using R grdid from RDINR as used for diagnostic output.
! Altered 15-Jul-2010 : Adjusted scalinf the R grid when read in from RDINR.
! Altered 24-Jun-2009 : Changed to V3: PURE_HUBBLE_FLOW passed in call.
! Altered 10-Apr-2009 : RDINR now adopts R = Rold + V(r).dt
! Altered 31-Jan-2008 : Minor bug fix.
! Altered 29-Dec-2008 : Based on SET_RV_HYDRO_MODEL_V2
!                       Written to allow smaller outer radius.
!                       Code now assumes R=Ro+V(r).t
!
	INTEGER ND
	INTEGER LU
!
! The following vectors and variables are returned.
!
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
	REAL*8 RMAX
	REAL*8 RCORE
!
! The following must be input
!
	REAL*8 RMAX_ON_RCORE
	REAL*8 SN_AGE_DAYS
	INTEGER N_IB_INS
	INTEGER N_OB_INS
	LOGICAL RDINR
	LOGICAL PURE_HUBBLE_FLOW
!
! Local arrays and variables.
!
	REAL*8, ALLOCATABLE :: R_HYDRO(:)
	REAL*8, ALLOCATABLE :: LOG_R_HYDRO(:)
	REAL*8, ALLOCATABLE :: V_HYDRO(:)
	REAL*8, ALLOCATABLE :: SIGMA_HYDRO(:)
	REAL*8, ALLOCATABLE :: DENSITY_HYDRO(:)
	REAL*8, ALLOCATABLE :: KAPPA_HYDRO(:)
	REAL*8, ALLOCATABLE :: TAU_HYDRO(:)
!
	REAL*8, ALLOCATABLE :: OLD_R(:)
	REAL*8, ALLOCATABLE :: LOG_OLD_R(:)
	REAL*8, ALLOCATABLE :: OLD_TAU(:)
	REAL*8, ALLOCATABLE :: COEF(:,:)
!
	REAL*8 TA(ND)
	REAL*8 XN(ND)
	REAL*8 KAPPA(ND)
	REAL*8 TAU(2*ND)
	REAL*8 LOG_R(2*ND)
	REAL*8 OLD_SN_AGE_DAYS
	REAL*8 SN_EXP_FACTOR
	REAL*8 dTAU
	REAL*8 dLOGR
	REAL*8 T1,T2
	REAL*8 NEXT_R
	REAL*8 DEN_SCL_FAC
	REAL*8 TAU_BEG,TAU_END
!
! These indicate the ratio of DTAU in the fine grid near the inner and outer boudary.
! They are currently hardwird, but could be passed as parameters.
!
	REAL*8, PARAMETER :: IB_FAC=2.5D0
	REAL*8, PARAMETER :: OB_FAC=2.0D0
!
	INTEGER ND_TMP
	INTEGER NX
	INTEGER NS
	INTEGER NSP
	INTEGER IOS
	INTEGER NOLD,NDOLD
	INTEGER I,J,J_SAV,L
	INTEGER N_INS_TOT
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER, PARAMETER :: IONE=1
	CHARACTER(LEN=200) STRING
!
	LUER=ERROR_LU()
	IF(.NOT. RDINR)THEN
	  IF(N_IB_INS .LT. 0 .OR. N_IB_INS .GT. ND/4)THEN
	    WRITE(LUER,'(A)')' '
	    WRITE(LUER,*)'Error is SET_RV_HYDRO_V3 - invalid N_IB_INS: N_IB_INS=',N_IB_INS
	    WRITE(LUER,*)'Valid range is 1 to ',ND/4
	    STOP
	   END IF
	  IF(N_OB_INS .LT. 0 .OR. N_OB_INS .GT. ND/4)THEN
	    WRITE(LUER,'(A)')' '
	    WRITE(LUER,*)'Error is SET_RV_HYDRO_V3 - invalid N_OB_INS: N_OB_INS=',N_OB_INS
	    WRITE(LUER,*)'Valid range is 1 to ',ND/4
	    STOP
	   END IF
	END IF
	 
	OPEN(UNIT=LU,FILE='SN_HYDRO_DATA',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL_V3 --- SN_HYDRO_DATA file not found'
	    WRITE(LUER,*)'Create file or EDIT option in VADAT'
	    STOP
	  END IF
!
! Get the number of data points in the HYDRO model, and the number 
! of species.
!
	NX=0; NSP=0; OLD_SN_AGE_DAYS=0.0D0
	DO WHILE (1 .EQ. 1)
	  STRING=' '
	  DO WHILE (STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING 
	  END DO
	  IF(INDEX(STRING,'Number of data points:') .NE. 0)THEN
	    I=INDEX(STRING,'points:')+7
	    READ(STRING(I:),*)NX
	  ELSE IF(INDEX(STRING,'Number of mass fractions:') .NE. 0)THEN
	    I=INDEX(STRING,'fractions:')+10
	    READ(STRING(I:),*)NSP
	  ELSE IF(INDEX(STRING,'Time(days) since explosion:') .NE. 0)THEN
	    I=INDEX(STRING,'explosion:')+10
	    READ(STRING(I:),*)OLD_SN_AGE_DAYS
	  ELSE IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	     IF(NSP .EQ. 0 .OR. NX .EQ. 0 .OR. OLD_SN_AGE_DAYS .EQ. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in RD_SN_DATA'
	      WRITE(LUER,*)'NSP or NX undefined'
	      STOP
	    ELSE
	      EXIT
	    END IF
	  END IF
	END DO
!
	SN_EXP_FACTOR=SN_AGE_DAYS/OLD_SN_AGE_DAYS
!
	ALLOCATE (V_HYDRO(NX));          V_HYDRO=0.0D0
	ALLOCATE (SIGMA_HYDRO(NX));      SIGMA_HYDRO=0.0D0
	ALLOCATE (DENSITY_HYDRO(NX));    DENSITY_HYDRO=0.0D0
	ALLOCATE (KAPPA_HYDRO(NX));      KAPPA_HYDRO=0.0D0
	ALLOCATE (R_HYDRO(NX));          R_HYDRO=0.0D0
	ALLOCATE (LOG_R_HYDRO(NX));      LOG_R_HYDRO=0.0D0
	ALLOCATE (TAU_HYDRO(NX));        TAU_HYDRO=0.0D0
!
	ALLOCATE (OLD_R(NX));		 OLD_R=0.0D0
	ALLOCATE (LOG_OLD_R(NX));        LOG_OLD_R=0.0D0
	ALLOCATE (OLD_TAU(NX));		 OLD_TAU=0.0D0
!
! Get basic HYDRO grid vectors.
!
	DO WHILE (1 .EQ. 1)
	  DO WHILE (STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING 
	  END DO
	  IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	    READ(LU,*)R_HYDRO
	  ELSE IF(INDEX(STRING,'Velocity') .NE. 0)THEN
	    READ(LU,*)V_HYDRO
	  ELSE IF(INDEX(STRING,'Sigma') .NE. 0)THEN
	    READ(LU,*)SIGMA_HYDRO
	  ELSE IF(INDEX(STRING,'Density') .NE. 0)THEN
	    READ(LU,*)DENSITY_HYDRO
	  ELSE IF(INDEX(STRING,'Kappa') .NE. 0)THEN
	    READ(LU,*)KAPPA_HYDRO
	  ELSE IF(INDEX(STRING,'ass fraction') .NE. 0)THEN
	    IF(R_HYDRO(1) .EQ. 0.0D0 .OR. 
	1             V_HYDRO(1) .EQ. 0.0D0 .OR.
	1             DENSITY_HYDRO(1) .EQ. 0.0D0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error reading SN data in SET_RV_HYDRO_MODEL'
	      WRITE(LUER,*)'R, V, the DENSITY is zero'
	      STOP
	    ELSE IF(KAPPA_HYDRO(1) .EQ. 0.0D0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Warning: reading SN data in SET_RV_HYDRO_MODEL'
	      WRITE(LUER,*)'KAPPA is zero'
	      EXIT
	    ELSE
	      EXIT
	    END IF
	  END IF
	  STRING=' '
	END DO
	WRITE(LUER,*)'Successfuly read SN data in SET_RV_HYDRO_V3'
	CLOSE(LU)
!
        IF(PURE_HUBBLE_FLOW)THEN
	  WRITE(6,*)'Setting Velocity so pure Hubble law in SET_RV_HYDRO_V3'
          T1=24.0D0*3600.0D0*1.0D+05*OLD_SN_AGE_DAYS/1.0D+10
          DO I=1,NX
            V_HYDRO(I)=R_HYDRO(I)/T1
            SIGMA_HYDRO(I)=0.0D0
          END DO
        END IF
!
! Adjust R grid of SN for expansion.
!	R_HYDRO=R_HYDRO*SN_EXP_FACTOR (Valid only if Hubble flow).
!	DENSITY_HYDRO=DENSITY_HYDRO/(SN_EXP_FACTOR**3)
!
! Valid for any flow where V, at a given mass coordinate, is constnat.
! Thus valid for a velocity law of the form R = Ro + V.t
!
! In the constant T1 we convert from days to seconds, and allow for the units of
! V (km/s) and R (10^10 cm).
!
	T1=24.0D0*3600.0D0*1.0D+05*(SN_AGE_DAYS-OLD_SN_AGE_DAYS)/1.0D+10
	DO I=1,NX
	  DEN_SCL_FAC=1.0D0/(1.0D0+T1*V_HYDRO(I)/R_HYDRO(I)*(SIGMA_HYDRO(I)+1.0D0))/
	1                   (1.0D0+T1*V_HYDRO(I)/R_HYDRO(I))**2
	  DENSITY_HYDRO(I)=DEN_SCL_FAC*DENSITY_HYDRO(I)
	  R_HYDRO(I)=R_HYDRO(I)+T1*V_HYDRO(I)
	END DO
!
	OPEN(UNIT=LU,FILE='OLD_SN_R_GRID',STATUS='UNKNOWN')
	  DO I=1,NX-1
	    WRITE(LU,'(I5,4ES16.8)')I,R_HYDRO(I),V_HYDRO(I),SIGMA_HYDRO(I),R_HYDRO(I)/R_HYDRO(I+1)
	  END DO
	  WRITE(LU,'(I5,3ES16.8)')I,R_HYDRO(I),V_HYDRO(I),SIGMA_HYDRO(I)
	CLOSE(LU)
!
	IF(RDINR)THEN
	  OPEN(UNIT=LU,STATUS='OLD',FILE='RDINR',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL_V3 --- File with R grid not found'
	    WRITE(LUER,*)'Create file or EDIT option in VADAT'
	    STOP
	  END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file.
!
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(LU,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU)
!
	  READ(LU,*,IOSTAT=IOS)TA(1),TA(1),NOLD,NDOLD
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL_V3 --- unable to read header in file with R grid'
	    STOP
	  END IF
!
! Check relative values.
!
	  IF(ND .NE. NDOLD)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error-NDOLD and ND are not equal in RDINR'
	    WRITE(LUER,*)'NDOLD=',NDOLD,' ND=',ND
	    STOP
	  END IF
!
! TA is used for everything but R which is all we want.
!
	  DO I=1,ND
	    READ(LU,*,IOSTAT=IOS)R(I),TA(I),TA(I),TA(I),TA(I),V(I)
	    IF(IOS .EQ. 0)READ(LU,*,IOSTAT=IOS)(TA(J),J=1,NOLD)
	    IF(IOS .NE. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL_V3 --- unable to R grid in file with R grid'
	      STOP
	    END IF
          END DO
!
! Make sure the R grid that is read in is compatible with the R grid found in SN_HYDRO_DATA.
! If R at the inner boundary differs signifcantly (currently 1 part in 10^6) from that in the
! SN_HYDRO_DATA file, we assume that the supplied R grid was for the previous time step.
!
	  T2=1.0D0-R(ND)/R_HYDRO(NX)
	  IF(T2 .GT. 1.0D-06)THEN
	    WRITE(LUER,*)'In SET_RV_HYDRO_MODEL_V3 I am assuming that the supplied R grid is from'
	    WRITE(LUER,*)'the previous time step.'
!
! In the constant T1 we convert from days to seconds, and allow for the units of
! V (km/s) and R (10^10 cm).
!
	    T1=24.0D0*3600.0D0*1.0D+05*(SN_AGE_DAYS-OLD_SN_AGE_DAYS)/1.0D+10
	    DO I=1,ND
	      R(I)=R(I)+T1*V(I)
	    END DO
	    T2=ABS(R(ND)/R_HYDRO(NX)-1.0D0)
	    IF(T2 .GT. 1.0D-06)THEN
	      WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL_V3'
	      WRITE(LUER,*)'Supplied R grid does not match current or previous time step'
	      WRITE(LUER,*)'      R(ND)=',R(ND)
	      WRITE(LUER,*)'R_HYDRO(NX)=',R_HYDRO(NX)
	      STOP
	    END IF
	  END IF
!
! Ensure we have an exact match between R_HYDRO and new R grid.
! Currently we assume that the inner boundary of model is not changing.
! The outer boundary can be moved inwards with time.
!
!	  WRITE(6,*)'Scale factor in SET_RV_HYDRO_V3 is',R_HYDRO(NX)/R(ND)
	  DO I=1,ND               	!need to issue a warning here
	    R(I)=R_HYDRO(NX)*(R(I)/R(ND))
	  END DO
!
! Check that RMAX is compatible.
!
	  IF(R(1) .GT. R_HYDRO(1))THEN
	    R(1)=R_HYDRO(1)
	    IF(R(1) .LT. R(2))THEN
	      WRITE(LUER,*)'Error setting RMAX at outer boundary in SET_RV_HYDRO_MODEL_V3'
	      WRITE(LUER,*)'We now have non-monotonic R grid'
	      WRITE(LUER,*)'R(1:4)=',R(1:4)
	      STOP
	    END IF
	  END IF
	  RMAX=R(1); RCORE=R(ND)
	  dLOGR=0.0D0; dTAU=0.0D0		!Set for diagnostic output file.
	  WRITE(LUER,*)'Read in R grid from RDINR in SET_RV_HYDRO_MODEL'
!
	ELSE
!
! Set RMAX and RCORE. These are returned and are not take from VADAT.
!
	  RCORE=R_HYDRO(NX)
	  IF(RMAX_ON_RCORE .GT. 1.0D0)THEN
	    RMAX=RMAX_ON_RCORE*RCORE
	    IF(RMAX .GT. (1.0D0+1.0D-07)*R_HYDRO(1))THEN
	      WRITE(LUER,'(A)')' Error is SET_RV_HYDRO_MODEL'
	      WRITE(LUER,'(A)')' RMAX_ON_RCORE is too large'
	      WRITE(LUER,'(A,F20.10)')' RMAX_ON_RCORE must be =< ',R_HYDRO(1)/R_HYDRO(NX)
	      STOP
	    ELSE IF(RMAX .GT. (1.0D0-1.0D-07)*R_HYDRO(1))THEN
	      RMAX=R_HYDRO(1)
	    END IF
	  ELSE
	    RMAX=R_HYDRO(1)
	  END IF
!
! Compute TAU. Since we are only setting the grid, a simple trapzoidal
! rule integration will suffice. We assume KAPPA does not change during
! the expansion.
!
	  KAPPA_HYDRO=1.0D+10*KAPPA_HYDRO*DENSITY_HYDRO
	  TAU_HYDRO(1)=KAPPA_HYDRO(1)*R_HYDRO(1)
	  T1=LOG(KAPPA_HYDRO(2)/KAPPA_HYDRO(1))/LOG(R_HYDRO(1)/R_HYDRO(2))
	  IF(T1 .GT. 1.5D0)TAU_HYDRO(1)=TAU_HYDRO(1)/(T1-1)
	  DO I=2,NX
	    TAU_HYDRO(I)=TAU_HYDRO(I-1)+0.5D0*(KAPPA_HYDRO(I-1)+KAPPA_HYDRO(I))*
	1                  (R_HYDRO(I-1)-R_HYDRO(I))
	  END DO
!
! Reduce RMAX if necessary. We store the reduced HYDRO grid in OLD_...
!
	  IF(RMAX .LT. R_HYDRO(1))THEN
	    J=1
	    DO WHILE (RMAX .LT. R_HYDRO(J))
	      J=J+1
	    END DO
	    NS=NX-J+2
	    OLD_R(2:NS)=R_HYDRO(J:NX)
	    OLD_R(1)=RMAX
	    OLD_TAU(2:NS)=TAU_HYDRO(J:NX)
	    T1=(RMAX-R_HYDRO(J))/(R_HYDRO(J-1)-R_HYDRO(J))
	    OLD_TAU(1)=T1*TAU_HYDRO(J-1)+(1.0D0-T1)*TAU_HYDRO(J)
	  ELSE
	    NS=NX
	    OLD_R(1:NS)=R_HYDRO(1:NX)
	    OLD_TAU(1:NS)=TAU_HYDRO(1:NX)
	  END IF
!
! Estimate spacing to get required grid spacing.
!
	  OLD_TAU=LOG(OLD_TAU)
	  J=ND-1-N_OB_INS-N_IB_INS
	  dTAU=(OLD_TAU(NS)-OLD_TAU(1))/J
	  LOG_OLD_R=LOG(OLD_R)
	  dLOGR=(LOG_OLD_R(1)-LOG_OLD_R(NS))/J  		!(ND-7)
!
! Define the new radius grid. The step size in R corresponds to the smaller of
! dLOGR and dLOG_TAU. We create a "uniform" grid. The finer grid at the inner and
! uter boudaries is now only set when we set the FINAL grid.
!
! N_INS_TOTS is a little larger than necessary.
!
	  J=1; I=1
	  LOG_R(1)=LOG_OLD_R(1)
	  DO WHILE(1 .EQ. 1)
	    I=I+1
	    IF(I .GT. 2*ND)THEN
	      WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL_V3 --- LOG_R and TAU vectors too small'
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
	    TAU_BEG=T1*OLD_TAU(J)+(1.0D0-T1)*OLD_TAU(J+1)
!
! Compute dTAU step, and check to see if we need to use that.
!
	    J_SAV=J
	    NEXT_R=LOG_R(I-1)-dLOGR
	    DO WHILE(LOG_OLD_R(J+1) .GT. NEXT_R)
	      J=J+1
	    END DO
	    T1=(NEXT_R-LOG_OLD_R(J+1))/(LOG_OLD_R(J)-LOG_OLD_R(J+1))
	    TAU_END=T1*OLD_TAU(J)+(1.0D0-T1)*OLD_TAU(J+1)
!
	    IF(TAU_END-TAU_BEG .GT. dTAU)THEN
	      NEXT_R=LOG_R(I-1)-dLOGR*dTAU/(TAU_END-TAU_BEG)
	    END IF
	    J=J_SAV
	    LOG_R(I)=NEXT_R
!
! Check whether close enough to Outer Bondary.
!
	    IF(LOG_R(I)-1.5*(LOG_R(I-1)-LOG_R(I)) .LT. LOG_OLD_R(NS))EXIT
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
	  DO I=1,ND_TMP; TAU(I)=I; END DO
	  DO I=1,J
	    XN(I)=1.0D0+(I-1.0D0)*(ND_TMP-1.0D0)/(J-1.0D0)
	  END DO
	  CALL MON_INTERP(R,J,IONE,XN,ND,LOG_R,ND_TMP,TAU,ND_TMP)
	  LOG_R=R
	  ND_TMP=J
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
	      T2=IB_FAC*T2+1.0D0
	    END DO
	    T1=T1/T2
	    DO J=N_IB_INS,1,-1
	      LOG_R(I+J-1)=LOG_OLD_R(NS)+T1
	      T1=T1*IB_FAC
	    END DO
	    ND_TMP=I+N_IB_INS
	  END IF
	  LOG_R(ND_TMP)=LOG_OLD_R(NS)
!
! Add finer grid at outer boundary.
!
! Shift grid to allow for insertion of extra ponts
!
	  DO I=ND_TMP,2,-1
	     LOG_R(I+N_OB_INS)=LOG_R(I)
	  END DO
!
	  T1=LOG_R(1)-LOG_R(2)
	  IF(N_OB_INS .GE. 5)THEN
!
! The ration of successive step size is assumed to be FAC. The last step size is
! taken to be ~1/20 of the second last step size.
!
	    T2=1.0D0
	    DO I=1,N_OB_INS-1
	      T2=T2*OB_FAC+1.0D0
	    END DO
	    T1=T1/T2
	    LOG_R(2)=LOG_R(1)-0.05D0*T1
	    LOG_R(3)=LOG_R(1)-T1
	    DO I=3,N_OB_INS
	      T1=T1*OB_FAC
	      LOG_R(1+I)=LOG_R(I)-T1
	    END DO
	    ND_TMP=ND_TMP+N_OB_INS
!
	  ELSE IF(N_OB_INS .EQ. 4)THEN
	    LOG_R(2)=LOG_OLD_R(1)-0.003D0*T1
	    LOG_R(3)=LOG_OLD_R(1)-0.03D0*T1
	    LOG_R(4)=LOG_OLD_R(1)-0.1D0*T1
	    LOG_R(5)=LOG_OLD_R(1)-0.4D0*T1
	    ND_TMP=ND_TMP+4
!
	  ELSE IF(N_OB_INS .EQ. 3)THEN
	    LOG_R(2)=LOG_OLD_R(1)-0.01D0*T1
	    LOG_R(3)=LOG_OLD_R(1)-0.1D0*T1
	    LOG_R(4)=LOG_OLD_R(1)-0.3D0*T1
	    ND_TMP=ND_TMP+3
!
	  ELSE IF(N_OB_INS .EQ. 2)THEN
	    LOG_R(2)=LOG_OLD_R(1)-0.01D0*T1
	    LOG_R(3)=LOG_OLD_R(1)-0.3D0*T1
	    ND_TMP=ND_TMP+2
!
	  ELSE IF(N_OB_INS .EQ. 1)THEN
	    LOG_R(2)=LOG_OLD_R(1)-0.03D0*T1
	    ND_TMP=ND_TMP+1
	  END IF
!
	  R=EXP(LOG_R)
	  R(1)=OLD_R(1)
	  R(ND)=R_HYDRO(NX)
	  WRITE(LUER,*)R(1:ND)
	  WRITE(LUER,*)'Computed R grid in SET_RV_HYDRO_MODEL_V3'
!
	END IF
!
! Compute V and SIGMA.
!
	ALLOCATE (COEF(NX,4))
	CALL MON_INT_FUNS_V2(COEF,V_HYDRO,R_HYDRO,NX)
	J=1
	DO I=1,ND
500	  IF(R(I) .GE. R_HYDRO(J+1))THEN
	    T1=R(I)-R_HYDRO(J)
	    V(I)=((COEF(J,1)*T1+COEF(J,2))*T1+COEF(J,3))*T1+COEF(J,4)
	    SIGMA(I)=((3.0D0*COEF(J,1)*T1+2.0D0*COEF(J,2))*T1+COEF(J,3))*R(I)/V(I)-1.0D0
	  ELSE
	    J=J+1
	    GOTO 500
	  END IF
	END DO
!
! If very close to Hubble law, we set SIGMA to zero.
!
	IF(MAXVAL(ABS(SIGMA)) .LT. 1.0D-04)SIGMA(1:ND)=0.0D0
!
	LOG_R(1:ND)=LOG(R(1:ND))
	LOG_R_HYDRO(1:NX)=LOG(R_HYDRO(1:NX))
	CALL MON_INTERP(KAPPA,ND,IONE,LOG_R,ND,KAPPA_HYDRO,NX,LOG_R_HYDRO,NX)
	TAU(1)=KAPPA(1)*R(1)
	T1=LOG(KAPPA(2)/KAPPA(1))/LOG(R(1)/R(2))
	IF(T1 .GT. 1.5D0)TAU(1)=TAU(1)/(T1-1)
	DO I=2,ND
	  TAU(I)=TAU(I-1)+0.5D0*(KAPPA(I-1)+KAPPA(I))*(R(I-1)-R(I))
	END DO
!
	OPEN(UNIT=LU,FILE='NEW_SN_R_GRID',STATUS='UNKNOWN')
	  WRITE(LU,'(/,3X,A,F6.3,A,F6.3,A)')'  dLOGR was ',0.4343*dLOGR,'(',EXP(dLOGR),')'
	  WRITE(LU,'(3X,A,F6.3,A,F6.3,A)')'dLOGTAU was ',0.4343*dTAU,'(',EXP(dTAU),')'
	  WRITE(LU,'(3X,A,F8.4)')'Rmax=',R(1)/R(ND)
	  WRITE(LU,'(/,A,15X,A,3(5X,A))')'    I','R','V(km/s)','  SIGMA','    TAU'
	  DO I=1,ND-1
	    WRITE(LU,'(I5,ES16.8,7ES12.4)')I,R(I),V(I),SIGMA(I),TAU(I),R(I)/R(I+1),LOG10(TAU(I+1)/TAU(I))
	  END DO
	  WRITE(LU,'(I5,ES16.8,3ES12.4)')ND,R(ND),V(ND),SIGMA(ND),TAU(ND)
	CLOSE(LU)
!
	DEALLOCATE (R_HYDRO, LOG_R_HYDRO, V_HYDRO, SIGMA_HYDRO)
	DEALLOCATE (DENSITY_HYDRO, KAPPA_HYDRO, TAU_HYDRO)
	DEALLOCATE (OLD_R, LOG_OLD_R, OLD_TAU, COEF)
!
	RETURN
	END
