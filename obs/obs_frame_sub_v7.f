!
! Subroutine to compute line profile in the observers frame.
!
	SUBROUTINE OBS_FRAME_SUB_V7(ETA_CMF,CHI_CMF,FREQ_CMF,
	1            R,V,T,ED,ND,NCF,
	1            P,MU_AT_RMAX,HQW_AT_RMAX,NC,NP,
	1            OBS_FREQ,OBS_FLUX,NOS,
	1            MAX_DEL_V_RES_ZONE,
	1            TAU_MAX,ES_DTAU,N_INS_OBS,INT_METHOD,
	1            WRITE_IP,DO_REL_CORRECTIONS,PLANE_PARALLEL)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 04-May-2009 : Changed to handle a hollow core. This is important
!                         for SN.
! Altered 26-Aug-2007 : Added additional deallocation of ETA_RAY_CMF etc.
! Altered 11-Feb-2006 : Minor bug fix which caused some frequency acces to go
!                         outside arrays. Error detecting inserted.
! Altered 29-Aug-2003 : Relativistic corrections included. Corrections were
!                         inserted in a manner so as not to effect the operation
!                         of the code when relativistic corrections were
!                         switched off.
! Altered 06-Jun-2003 : Evaluation of Z(I) changed to avoid bad FP numbers
!                         with the INTEL compiler.
! Altered 08-Jul-2002 : N_INS_OBS inserted in call
!                         Changed to V4
! Altered 14-Aug-2000 : WRITE_IP option installed.
! Altered 26-Jan-2000 : NUM_BANDS made variable, with a maximum of 100.
!                         DONE to save allocation of space.
! Altered 03-Jan-2000 : Extra points inserted along ray when DTAU_ES >
!                         0.25D0 when TAU_ES < 20. In V83 O star models
!                         there were large increments (DTAU > 1) in
!                         the continuum optical depth near TAU=1.
!                         We now check that Z_RAY is big enough for the
!                         the point insertions.
! Altered 22-Jan-1998 : Treatment of optically thick lines at the outer
!                         boundary improved.
!
	INTEGER ND
	INTEGER NCF
!
! Contain J, ETA and CHI from the Comoving-frame computations. Must be
! supplied by calling routine. At present it is assumed that these are
! on the same depth grid, and on the same frequency grid.
!
	REAL(KIND=LDP) ETA_CMF(ND,NCF)
	REAL(KIND=LDP) CHI_CMF(ND,NCF)
	REAL(KIND=LDP) FREQ_CMF(NCF)
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) ED(ND)
!
! Impact parameters: For rays not striking the core, P must be defined by
! the R grid. NP should be ND+NC
!
	INTEGER NC
	INTEGER NP
	REAL(KIND=LDP) P(NP)
	REAL(KIND=LDP) MU_AT_RMAX(NP)
	REAL(KIND=LDP) HQW_AT_RMAX(NP)
!
! Observer's frame frequencies in units of 10^15 Hz. Should be monotonically
! decreasing.
!
	INTEGER NOS
	REAL(KIND=LDP) OBS_FREQ(NOS)
	REAL(KIND=LDP) OBS_FLUX(NOS)
!
! Maximum velocity spacing across consequitive grid points. Ideally this
! should be least than 0.25 of a Doppler width.
!
	REAL(KIND=LDP) MAX_DEL_V_RES_ZONE(ND)
!
	REAL(KIND=LDP) TAU_MAX
	REAL(KIND=LDP) ES_DTAU
	INTEGER N_INS_OBS
	CHARACTER*(*) INT_METHOD
!
! Local vectors and arrays.
!
	REAL(KIND=LDP) VMU(ND)		!Velocity*(Z/R) on R grid.
	REAL(KIND=LDP) Z(ND)		!Distance along ray on R grid.
	REAL(KIND=LDP) ESEC(ND)		!Electron scattering opacity
	REAL(KIND=LDP) TAU_ES(ND)	!Rough TAU_ES along ray
	REAL(KIND=LDP) DTAU_ES(ND)	!Rough DTAU_ES along ray
!
! Emissivities and Opacities in the CMF frame. These are created in this
! routing from ETA_CMF and CHI_CMF, and are on the fine grid used to solve
! the transfer equation.
!
	REAL(KIND=LDP), ALLOCATABLE ::  ETA_CMF_RAY(:,:)
	REAL(KIND=LDP), ALLOCATABLE ::  CHI_CMF_RAY(:,:)
	REAL(KIND=LDP), ALLOCATABLE ::  SM_FREQ_CMF(:)
!
! Intensity as a function of impact parameter and frequency in the observers
! frame. Need for interpreting observations of Eta Car, and for understanding
! interferometric observations.
!
	LOGICAL WRITE_IP
	REAL(KIND=LDP), ALLOCATABLE ::  IP_OBS(:,:)
!
	LOGICAL DO_REL_CORRECTIONS
	LOGICAL PLANE_PARALLEL
!
! Vectors defined along a ray. HALF_DZ, DZSQ_ON_12, and RECIP_DEL_Z are used
! to minimize computations.
!
	REAL(KIND=LDP), ALLOCATABLE ::  ETA_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE ::  CHI_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE ::  R_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE ::  V_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE ::  VMU_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE ::  GAM_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE ::  NU_ON_NU0_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE ::  Z_RAY(:)
	REAL(KIND=LDP), ALLOCATABLE ::  TAU(:)
	REAL(KIND=LDP), ALLOCATABLE ::  dZ(:)
	REAL(KIND=LDP), ALLOCATABLE ::  HALF_DZ(:)
	REAL(KIND=LDP), ALLOCATABLE ::  DZSQ_ON_12(:)
	REAL(KIND=LDP), ALLOCATABLE ::  RECIP_DEL_Z(:)
!
	REAL(KIND=LDP), ALLOCATABLE ::  A0(:)
	REAL(KIND=LDP), ALLOCATABLE ::  A1(:)
	REAL(KIND=LDP), ALLOCATABLE ::  A2(:)
	REAL(KIND=LDP), ALLOCATABLE ::  A3(:)
	REAL(KIND=LDP), ALLOCATABLE ::  A4(:)
!
	REAL(KIND=LDP), ALLOCATABLE ::  EE(:)
	REAL(KIND=LDP), ALLOCATABLE ::  E0(:)
	REAL(KIND=LDP), ALLOCATABLE ::  E1(:)
	REAL(KIND=LDP), ALLOCATABLE ::  E2(:)
	REAL(KIND=LDP), ALLOCATABLE ::  E3(:)
!
	REAL(KIND=LDP), ALLOCATABLE ::  DTAU(:)
	REAL(KIND=LDP), ALLOCATABLE ::  S(:)
	REAL(KIND=LDP), ALLOCATABLE ::  dS(:)
!
	INTEGER NOS_INC
	INTEGER OUT_ML
	INTEGER SM_NCF
	INTEGER NUM_BANDS
	INTEGER, ALLOCATABLE :: STRT_INDX_CMF(:)       !NUM_BANDS
	INTEGER, ALLOCATABLE :: END_INDX_CMF(:)        !NUM_BANDS
	INTEGER, ALLOCATABLE :: INDX(:)                !NRAY
!
! NR is the number of points along a ray for Z .GE. 0. It is a function
! of LS.
!
	INTEGER NR
!
! NRAY total number of points along ray. For rays stringing the core,
! NRAY=NR, otherwise NRAY=2*NR-1
!
	INTEGER NRAY
	INTEGER SM_NRAY
!
	INTEGER I,J,K,ML,LS
	INTEGER NP_TMP
	INTEGER IOS
	INTEGER NINS
	INTEGER NINS_CONT
!
	INTEGER LUER
	INTEGER NI
	INTEGER NI_MAX
	INTEGER ERR_COUNT
	INTEGER REC_SIZE,UNIT_SIZE
	INTEGER WORD_SIZE,N_PER_REC,ACCESS_F
	INTEGER IST,IEND
!
	REAL(KIND=LDP) MAX_VMU
	REAL(KIND=LDP) T1,T2,T3,T4
	REAL(KIND=LDP) PSQ
	REAL(KIND=LDP) C_KMS
	REAL(KIND=LDP) PAR_FLUX
!
	LOGICAL BACK_SIDE
	LOGICAL HOLLOW_CORE
	INTEGER ERROR_LU
	REAL(KIND=LDP) SPEED_OF_LIGHT,FUN_PI,PI
	EXTERNAL SPEED_OF_LIGHT,FUN_PI,ERROR_LU
!
	REAL(KIND=LDP), PARAMETER :: ONE=1
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
	CALL TUNE(1,'OBS_FRAME_SUB')
!
	C_KMS=SPEED_OF_LIGHT()*1.0E-05_LDP
	PI=FUN_PI()
	LUER=ERROR_LU()
	OBS_FLUX(1:NOS)=0.0_LDP
	ERR_COUNT=0
!
	IF(INT_METHOD .NE. 'STAU' .AND.
	1         INT_METHOD .NE. 'ETAZ')THEN
	  WRITE(LUER,*)'Error in OBS_FRAME_SUB_V4'
	  WRITE(LUER,*)'Invalid integration method'
	  WRITE(LUER,*)'Integration Method =',INT_METHOD
	  STOP
	END IF
	IF(TAU_MAX .LT. 10 .OR. TAU_MAX .GT. 30.0_LDP)THEN
	  WRITE(LUER,*)'Warning in OBS_FRAME_SUB_V4'
	  WRITE(LUER,*)'TAU_MAX is outside the recommended range of (10,20)'
	  WRITE(LUER,*)'TAU_MAX=',TAU_MAX
	END IF
	IF(TAU_MAX .LE. 0)THEN
	  WRITE(LUER,*)'Error in OBS_FRAME_SUB_V4'
	  WRITE(LUER,*)'TAU_MAX cannot be zero or -ve'
	  STOP
	END IF
!
! Ideally ES_DTAU should be 0.1 to 0.25
!
	IF(ES_DTAU .LT. 0.001_LDP .OR. ES_DTAU .GT. 1.0_LDP)THEN
	  WRITE(LUER,*)'Warning in OBS_FRAME_SUB_V4'
	  WRITE(LUER,*)'ES_DTAU is outside the recommended range'
	  WRITE(LUER,*)'ES_DTAU=',ES_DTAU
	END IF
	IF(ES_DTAU .LE. 0)THEN
	  WRITE(LUER,*)'Error in OBS_FRAME_SUB_V4.'
	  WRITE(LUER,*)'ES_DTAU cannot be zero or -ve.'
	  STOP
	END IF
!
	IF(N_INS_OBS .LT. 0)THEN
	  WRITE(LUER,*)'Warning in OBS_FRAME_SUB_V4'
	  WRITE(LUER,*)'N_INS_OBS is outside the recommended range'
	  WRITE(LUER,*)'N_INS_OBS=',N_INS_OBS
	END IF
!
	IF(WRITE_IP)THEN
	  ALLOCATE (IP_OBS(NP,NOS),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	    WRITE(LUER,*)'Unable to allocate memory for IP_OBS'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
	END IF
!
! Check that frequencies are monotonically decreasing.
!
	DO ML=2,NCF
	  IF(FREQ_CMF(ML) .GE. FREQ_CMF(ML-1))THEN
	    WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	    WRITE(LUER,*)'CMF frequencies not monotonically decreasing'
	    WRITE(LUER,*)'ML=',ML
	    WRITE(LUER,*)'NU(ML-1)=',FREQ_CMF(ML-1)
	    WRITE(LUER,*)'NU(ML)=',FREQ_CMF(ML)
	    STOP
	  END IF
	END DO
!
	DO ML=2,NOS
	  IF(OBS_FREQ(ML) .GE. OBS_FREQ(ML-1))THEN
	    WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	    WRITE(LUER,*)'OBS frequencies not monotonically decreasing'
	    WRITE(LUER,*)'ML=',ML
	    WRITE(LUER,*)'NU(ML-1)=',OBS_FREQ(ML-1)
	    WRITE(LUER,*)'NU(ML)=',OBS_FREQ(ML)
	    STOP
	  END IF
	END DO
!
! We split the computation into NUM_BANDS loops. In each inner loop, we compute
! the fluxes for NOS_INC observer's frame frequencies. Here we determine
! NUM_BANDS, and allocate the vectors which specify the frequency range
! need in the CMF for each set of observer's frequencies.
!
	NUM_BANDS=100
	IF(NOS/NUM_BANDS .LT. 1000)NUM_BANDS=NOS/1000
	IF(NUM_BANDS .LT. 1)NUM_BANDS=1
	IF(ALLOCATED(STRT_INDX_CMF))DEALLOCATE(STRT_INDX_CMF)
	IF(ALLOCATED(END_INDX_CMF))DEALLOCATE(END_INDX_CMF)
	ALLOCATE (STRT_INDX_CMF(NUM_BANDS),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (END_INDX_CMF(NUM_BANDS),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for STRT_INDX_CMF'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
!
! Compute electron scattering opacity.
!
	ESEC(1:ND)=6.65E-15_LDP*ED(1:ND)
!
! If we have a SN model, we determine if we have a hollow core.
! We use TAU_MAX.
!
	TAU_ES(1)=0.0_LDP
	DO I=2,ND
	   TAU_ES(I)=TAU_ES(I-1)+0.5_LDP*(ESEC(I-1)+ESEC(I))*(R(I-1)-R(I))
	END DO
	HOLLOW_CORE=.FALSE.
	IF(V(ND) .GT. 100.0_LDP .AND. TAU_ES(ND) .LT. TAU_MAX)HOLLOW_CORE=.TRUE.
!
! 2* since ray can extend across whole atmosphere. The 2TAU_MAX/ES_DTAU is
! for the insertion of extra continuum points if DTAU_ES > ES_DTAU (for
! TAU_ES < TAU_MAX).
!
	T1=MINVAL(MAX_DEL_V_RES_ZONE)
	NI_MAX=2*( (V(1)/T1)+ (N_INS_OBS+1)*ND + 2.0_LDP*TAU_MAX/ES_DTAU )
	ALLOCATE (Z_RAY(NI_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for Z_RAY'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
!
! Only do for NP-1 since last ray of a single point contains no flux.
!
	NP_TMP=NP-1
	IF(PLANE_PARALLEL)NP_TMP=NC
	DO 50000 LS=1,NP_TMP
!
	  IF(LS .GT. NC)THEN
	    NI=ND-(LS-NC-1)
	  ELSE
	    NI=ND
	  END IF
!
	  IF(PLANE_PARALLEL)THEN
	    DO I=1,NI
	      Z(I)=R(I)/MU_AT_RMAX(LS)
	      VMU(I)=V(I)*R(I)/Z(I)
	    END DO
	  ELSE
	    DO I=1,NI
	      Z(I)=SQRT( (R(I)+P(LS))*(R(I)-P(LS)) )
	      VMU(I)=V(I)*Z(I)/R(I)
	    END DO
	    IF(LS .GT. NC)THEN
	      Z(NI)=0.0_LDP
	      VMU(NI)=0.0_LDP
	    END IF
	  END IF
!
! Compute optical depth increments in Electron scattering optical depth.
!
	  TAU_ES(1)=0.0_LDP
	  DO I=2,NI
	   DTAU_ES(I-1)=0.5_LDP*(ESEC(I-1)+ESEC(I))*(Z(I-1)-Z(I))
	   TAU_ES(I)=TAU_ES(I-1)+DTAU_ES(I-1)
	  END DO
	  DTAU_ES(NI)=0.0_LDP
!
	  Z_RAY(1)=Z(1)
	  J=1
	  DO I=2,NI
	    NINS=ABS(VMU(I-1)-VMU(I))/MAX_DEL_V_RES_ZONE(I-1)+1
            IF(TAU_ES(I) .LE. TAU_MAX)THEN
	      NINS_CONT=DTAU_ES(I-1)/ES_DTAU+1
	    ELSE
	      NINS_CONT=1
	    END IF
	    NINS=MAX(NINS,NINS_CONT)
	    NINS=MAX(NINS,N_INS_OBS+1)	
	    IF(J+NINS .GT. NI_MAX)THEN
	      WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	      WRITE(LUER,*)'NI_MAX is too small: NI_MAX=',NI_MAX
	      WRITE(LUER,*)'I=',I
	      WRITE(LUER,*)'VMU(I-1)=',VMU(I-1),'VMU(I)=',VMU(I)
	      WRITE(LUER,*)'MAX_DEL_V_RES_ZONE(I-1)=',MAX_DEL_V_RES_ZONE(I-1)
	      STOP
	    END IF
	    IF(NINS .GT. 1)THEN
	      T1=(Z(I)-Z(I-1))/NINS
	      DO K=2,NINS
	        J=J+1
	        Z_RAY(J)=Z(I-1)+(K-1)*T1
	      END DO
	    END IF
	    J=J+1
	    Z_RAY(J)=Z(I)
	  END DO
	  NR=J
!
	  IF(LS .NE. 1)THEN
	    DEALLOCATE (R_RAY)
	    DEALLOCATE (V_RAY)
	    DEALLOCATE (VMU_RAY)
	    DEALLOCATE (GAM_RAY)
	    DEALLOCATE (NU_ON_NU0_RAY)
	    DEALLOCATE (dZ)
	    DEALLOCATE (HALF_DZ)
	    DEALLOCATE (DZSQ_ON_12)
	    DEALLOCATE (RECIP_DEL_Z)
	    DEALLOCATE (TAU)
	    DEALLOCATE (ETA_VEC)
	    DEALLOCATE (CHI_VEC)
	    DEALLOCATE (INDX)
	    DEALLOCATE (A0)
	    DEALLOCATE (A1)
	    DEALLOCATE (A2)
	    DEALLOCATE (A3)
	    DEALLOCATE (A4)
	    DEALLOCATE (EE)
	    DEALLOCATE (E0)
	    DEALLOCATE (E1)
	    DEALLOCATE (E2)
	    DEALLOCATE (E3)
	    DEALLOCATE (DTAU)
	    DEALLOCATE (S)
	    DEALLOCATE (dS)
!
! We deallocate ETA_CMF_RAY etc here to keep allocated data contiguous.
!
	    IF(ALLOCATED(ETA_CMF_RAY))THEN
	      DEALLOCATE (ETA_CMF_RAY)
	      DEALLOCATE (CHI_CMF_RAY)
	      DEALLOCATE (SM_FREQ_CMF)
	    END IF
	  END IF
!
	  NRAY=2*NR-1
	  IF(LS .LE. NC)NRAY=NR
	  IF(LS .LE. NC .AND. HOLLOW_CORE)NRAY=2*NR
!
	  ALLOCATE (R_RAY(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (V_RAY(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VMU_RAY(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_RAY(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NU_ON_NU0_RAY(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dZ(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HALF_DZ(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DZSQ_ON_12(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RECIP_DEL_Z(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TAU(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ETA_VEC(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_VEC(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (INDX(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (A0(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (A1(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (A2(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (A3(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (A4(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EE(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (E0(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (E1(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (E2(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (E3(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAU(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (S(NRAY),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dS(NRAY),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	    WRITE(LUER,*)'Unable to allocate R_RAY etc'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
!
	  IF(PLANE_PARALLEL)THEN
	    DO I=1,NR
	      R_RAY(I)=Z_RAY(I)*MU_AT_RMAX(LS)
	    END DO
	  ELSE
	    PSQ=P(LS)*P(LS)
	    DO I=2,NR-1
	      R_RAY(I)=SQRT(PSQ+Z_RAY(I)*Z_RAY(I))
	    END DO
	  END IF
	  R_RAY(1)=R(1)
	  R_RAY(NR)=R(NI)
!
	  CALL MON_INTERP_FAST(V_RAY,NRAY,IONE,R_RAY,NR,V,ND,R,ND)
	  IF(PLANE_PARALLEL)THEN
	    VMU_RAY(1:NR)=V_RAY(1:NR)*MU_AT_RMAX(LS)/C_KMS
	  ELSE
	    DO I=1,NR
	      VMU_RAY(I)=V_RAY(I)*Z_RAY(I)/R_RAY(I)/C_KMS
	    END DO
	  END IF
!
	  IF(LS .GT. NC)THEN
!
! Since LS>NC extend trajectory to extend across the whole of the atmosphere.
!
	    DO I=1,NR-1
	      J=2*NR-I
	      R_RAY(J)=R_RAY(I)
	      V_RAY(J)=V_RAY(I)
	      Z_RAY(J)=-Z_RAY(I)
	      VMU_RAY(J)=-VMU_RAY(I)
	    END DO
	  ELSE IF(LS .LE. NC .AND. HOLLOW_CORE)THEN
	    DO I=1,NR
	      J=2*NR+1-I
	      R_RAY(J)=R_RAY(I)
	      V_RAY(J)=V_RAY(I)
	      Z_RAY(J)=-Z_RAY(I)
	      VMU_RAY(J)=-VMU_RAY(I)
	    END DO
	  END IF
!
! NB: VMU_RAY has already been normalized by C_KMS. V_RAY has not.
!
	  IF(DO_REL_CORRECTIONS)THEN
	    DO I=1,NRAY
	      GAM_RAY(I)=1.0_LDP/SQRT(1.0_LDP-V_RAY(I)*V_RAY(I)/C_KMS/C_KMS)
	      NU_ON_NU0_RAY(I)=1.0_LDP/GAM_RAY(I)/(1.0_LDP-VMU_RAY(I))
	    END DO
	  ELSE
	     GAM_RAY(:)=1.0_LDP
	     NU_ON_NU0_RAY(:)=1.0_LDP
	  END IF
!
	  DO J=1,NRAY-1
	    dZ(J)=Z_RAY(J)-Z_RAY(J+1)
	    HALF_DZ(J)=0.5_LDP*DZ(J)
	    DZSQ_ON_12(J)=HALF_DZ(J)*HALF_DZ(J)/3.0_LDP		!ie. DZ/12.0D0
	  END DO
	  DO J=2,NRAY-1
	    RECIP_DEL_Z(J)=1.0_LDP/(Z_RAY(J-1)-Z_RAY(J+1))
	  END DO
!
! 
!
! We split the loop over observer's frame frequencies into two loops.
! While this makes the code logic more complicated, it can be used to minimize
! the amount of storage use by ETA_CMF_RAY and CHI_CMF_RAY. Without this
! loop, ETA_CMF_RY would have to be dimensioned ETA_CMF_RAY(NR,CMF) which
! is a very big array for NCF large, since NR generally will be > 100.
!
	  NOS_INC=(NOS+NUM_BANDS-1)/NUM_BANDS
!
! Determine the sections of ETA_CMF and CHI_CMF need for each band of
! observer's frame frequencies.
!
	  SM_NCF=0
	  DO OUT_ML=1,NUM_BANDS
!
! Determine range of frequencies in CMF required for determination of the
! observer's frame fluxes. We do this for all bands. The maximum length
! is then used to declare ETA_CMF_RAY and CHI_CMF_RAY so that we only
! allocate them once for each band. Factor of 1.05 is used to expand the
! window slightly.
!
	    MAX_VMU=MAXVAL(VMU_RAY)
	    T2=1.0_LDP+1.2_LDP*MAX_VMU
	    IF(DO_REL_CORRECTIONS)THEN
!	       T2=(MAXVAL((1.0D0-VMU_RAY)*GAM_RAY) -1.0D0)*1.05D0+1.0D0
	       T2=MAXVAL((1.0_LDP-VMU_RAY)*GAM_RAY)
	    END IF
	    I=MIN(NOS,1+(OUT_ML-1)*NOS_INC)
	    T1=OBS_FREQ(I)*T2
	    J=2
	    DO WHILE (T1 .LT. FREQ_CMF(J))
	      J=J+1
	    END DO
	    STRT_INDX_CMF(OUT_ML)=MAX(1,J-10)
!
	    T2=1.0_LDP-1.2_LDP*MAX_VMU
	    IF(DO_REL_CORRECTIONS)THEN
!	      T2=(MINVAL((1.0D0-VMU_RAY)*GAM_RAY) -1.0D0)*1.05D0+1.0D0
	      T2=MINVAL((1.0_LDP-VMU_RAY)*GAM_RAY)
	    END IF
	    I=MIN(NOS,OUT_ML*NOS_INC)
	    T1=OBS_FREQ(I)*T2
	    J=NCF
	    DO WHILE (T1 .GT. FREQ_CMF(J))
	      J=J-1
	    END DO
	    END_INDX_CMF(OUT_ML)=MIN(NCF,J+10)
	    WRITE(301,'(5(A,I7))')'LS=',LS,'OUT_ML=',OUT_ML,
	1           'IST=',STRT_INDX_CMF(OUT_ML),'IND=',END_INDX_CMF(OUT_ML),
	1           'SM_NCF=',END_INDX_CMF(OUT_ML)-STRT_INDX_CMF(OUT_ML)
	    SM_NCF=MAX(SM_NCF,END_INDX_CMF(OUT_ML)-STRT_INDX_CMF(OUT_ML)+1)
	  END DO
!
! We only declare ETA_CMF_RAY and CHI_CMF_RAY using NR (and not NRAY) since
! in the CMF these quanties are symmetric on a ray. We utilize this symmetry
! when accessing ETA_CMF_RAY and CHI_CMF_RAY.
!
	  IF(ALLOCATED(ETA_CMF_RAY))THEN
	    DEALLOCATE (ETA_CMF_RAY)
	    DEALLOCATE (CHI_CMF_RAY)
	    DEALLOCATE (SM_FREQ_CMF)
	  END IF
	  ALLOCATE (SM_FREQ_CMF(SM_NCF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ETA_CMF_RAY(NR,SM_NCF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_CMF_RAY(NR,SM_NCF),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in OBS_FRAME_SUB'
	    WRITE(LUER,*)'Unable to allocate ETA_CMF_RAY etc'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
!
	  DO 45000 OUT_ML=1,NUM_BANDS
	    SM_NCF=END_INDX_CMF(OUT_ML)-STRT_INDX_CMF(OUT_ML)+1
!
	    CALL TUNE(1,'MON_INTERP')
	      ETA_CMF_RAY(:,:)=0.0_LDP
	      CHI_CMF_RAY(:,:)=0.0_LDP
	      CALL MON_INTERP_FAST(ETA_CMF_RAY,NR,SM_NCF,
	1            R_RAY,NR,ETA_CMF(1,STRT_INDX_CMF(OUT_ML)),ND,R,ND)
	      CALL MON_INTERP_FAST(CHI_CMF_RAY,NR,SM_NCF,
	1            R_RAY,NR,CHI_CMF(1,STRT_INDX_CMF(OUT_ML)),ND,R,ND)
	    CALL TUNE(2,'MON_INTERP')
	    SM_FREQ_CMF(1:SM_NCF)=
	1       FREQ_CMF(STRT_INDX_CMF(OUT_ML):END_INDX_CMF(OUT_ML))
!
	    INDX(1:NRAY)=1
!
! NB: In this loop we corrupt NRAY because we cease our integration at TAU_MAX.
!     Thus at the beginning of each loop, we restore NRAY.
!
	  CALL TUNE(1,'ML')
!
!$OMP PARALLEL PRIVATE(ETA_VEC, CHI_VEC, TAU, DTAU, S, dS, EE, E0, E1, E2, E3,
!$OMP1   A0,A1,A2,A3,A4, PAR_FLUX, SM_NRAY, INDX, I, J, K, T1, T2, T3, T4,
!$OMP1   IST, IEND, BACK_SIDE )
!$OMP DO
	  DO 40000 ML=1+(OUT_ML-1)*NOS_INC,MIN(NOS,OUT_ML*NOS_INC)
!
	    ETA_VEC(1:NRAY)=0.0_LDP
	    CHI_VEC(1:NRAY)=0.0_LDP
!
!	    K=INDX(1)
	    K=1   !INDX(1)
	    INDX(1)=1
	    T1=OBS_FREQ(ML)*(1.0_LDP-VMU_RAY(1))*GAM_RAY(1)
	    DO WHILE(T1 .LT. SM_FREQ_CMF(K))
	      K=K+1
	      IF(K .GT. SM_NCF)THEN
	        WRITE(LUER,*)'Error in OBS_FRAME_SUB_V5'
	        WRITE(LUER,*)'Bad K indx: K=',K
	        WRITE(LUER,*)T1,SM_FREQ_CMF(1),SM_FREQ_CMF(SM_NCF)
	        WRITE(LUER,*)'LS, OUT_ML, ML=',LS,OUT_ML,ML
	      END IF
	    END DO
	    INDX(1)=K
	    T2=(T1-SM_FREQ_CMF(K))/(SM_FREQ_CMF(K-1)-SM_FREQ_CMF(K))
	    ETA_VEC(1)=
	1      T2*(ETA_CMF_RAY(1,K-1)-ETA_CMF_RAY(1,K))+ETA_CMF_RAY(1,K)
	    CHI_VEC(1)=
	1      T2*(CHI_CMF_RAY(1,K-1)-CHI_CMF_RAY(1,K))+CHI_CMF_RAY(1,K)
!
	    IF(DO_REL_CORRECTIONS)THEN
	      ETA_VEC(1)=ETA_VEC(1)*NU_ON_NU0_RAY(1)**2
	      CHI_VEC(1)=CHI_VEC(1)/NU_ON_NU0_RAY(1)
	    END IF
!
! Set TAU(1) to a reasonable value at the outer boundary. This avoids
! having TAU(1)=0 and TAU2(300) which can cause wild oscillations in the
! Flux for lines which are optically thick at the outer boundary
! (eg CIII 976 in WCE stars, MGII 2780 in LBVs).
!
	    TAU(1)=CHI_VEC(1)*HALF_DZ(2)
	    SM_NRAY=NRAY
!	    CALL TUNE(1,'FREQ INTERP')
	    DO I=2,NRAY
	      J=I
	      IF(I .GT. NR)J=2*NR-I
!
! Interpolate in frequency to transform ETA and CHI from the comoving frame to
! the observers frame for this ray and frequency.
!
	      K=INDX(I-1)      !INDX(I)
!	      K=INDX(I)
	      T1=OBS_FREQ(ML)*(1.0_LDP-VMU_RAY(I))*GAM_RAY(I)
!	      DO WHILE(T1 .LT. SM_FREQ_CMF(K))
	      DO WHILE(T1 .GT. SM_FREQ_CMF(K-1))
	        K=K-1
!	        IF(K .GT. SM_NCF)THEN
	        IF(K .LE. 1)THEN
	          WRITE(LUER,*)'Error in OBS_FRAME_SUB_V5'
	          WRITE(LUER,*)'Bad K indx: K=',K
	          WRITE(LUER,*)T1,SM_FREQ_CMF(1),SM_FREQ_CMF(SM_NCF)
	          WRITE(LUER,*)'LS, OUT_ML, ML, I=',LS,OUT_ML,ML,I
	        END IF
	      END DO
	      INDX(I)=K
	      T2=(T1-SM_FREQ_CMF(K))/(SM_FREQ_CMF(K-1)-SM_FREQ_CMF(K))
	      ETA_VEC(I)=
	1       T2*(ETA_CMF_RAY(J,K-1)-ETA_CMF_RAY(J,K))+ETA_CMF_RAY(J,K)
	      CHI_VEC(I)=
	1       T2*(CHI_CMF_RAY(J,K-1)-CHI_CMF_RAY(J,K))+CHI_CMF_RAY(J,K)
!
	      IF(DO_REL_CORRECTIONS)THEN
	        ETA_VEC(I)=ETA_VEC(I)*NU_ON_NU0_RAY(I)**2
	        CHI_VEC(I)=CHI_VEC(I)/NU_ON_NU0_RAY(I)
	      END IF
!
! With a hollow core, we need to make sure we don't include
! opacity due to the hollow core (located between NR and NR+1).
!
	      DTAU(I-1)=HALF_DZ(I-1)*(CHI_VEC(I-1)+CHI_VEC(I))
	      IF(I .EQ. NR+1 .AND. LS .LE. NC)THEN
	        TAU(I)=TAU(I-1)
	      ELSE
	        TAU(I)=TAU(I-1)+DTAU(I-1)
	      END IF
	      IF(LS .LE. NC .AND. HOLLOW_CORE)THEN
	        IF(TAU(I) .GT. TAU_MAX .AND. (I .LE. NR .OR. I .GT. NR+5))THEN
	          SM_NRAY=I
	          GOTO 2500
	        END IF
	      ELSE IF(TAU(I) .GT. TAU_MAX .AND. I .GT. 3)THEN
	        SM_NRAY=I
	        GOTO 2500
	      END IF
	    END DO
!
2500	    CONTINUE
!	    CALL TUNE(2,'FREQ INTERP')
!
! Use the first derivatives in conjunction with the Euler-MaucLarin summation
! formula to increase the accuracy of the optical depth integration. No
! correction is made for the integration over the end points. Note that TAU(I)
! already contains the trapezoidal estimate for the optical depth scale.
! The derivative at each node is estimated using the function values at the
! two adjacent data points.
!
!	    CALL TUNE(1,'TAU')
!	    T3=0.0D0
!	    T2=(CHI_VEC(1)-CHI_VEC(3))*RECIP_DEL_Z(2)
!	    DO I=3,SM_NRAY-1
!	      T1=T2
!	      T2=(CHI_VEC(I-1)-CHI_VEC(I+1))*RECIP_DEL_Z(I)
!	      T4=DZSQ_ON_12(I)*(T2-T1)
!	      IF(ABS(T4) .LT. 0.5*(TAU(I)-TAU(I-1)))T3=T3+T4
!             TAU(I)=TAU(I)+T3
!	    END DO
!	    TAU(SM_NRAY)=TAU(SM_NRAY)+T3
!	    CALL TUNE(2,'TAU')
!
! If we have a hollow core, and ray passes through core, we need to
! do the ray integration in two parts. This avoids issues with the
! "empty" interval inside the hollow core.
!
! To proceed, we first compute the flux coming from the backside
! of the hollow core. The flux from this integration then provides
! the incident flux for the NEAR-SIDE integration. To save uneccesary
! programming we loop over the integration section (35000 loop) twice
! when passing through the hollow core. Instead of 1 and NRAY_SM, we
! now use IST and IEND to define the integration limits.
!
	    IST=1; IEND=SM_NRAY
	    BACK_SIDE=.FALSE.
	    IF(LS .LE. NC .AND. HOLLOW_CORE)BACK_SIDE=.TRUE.
35000	    CONTINUE
!
! Check whether we need to do the BACK_SIDE loop.
!
	      IF(HOLLOW_CORE .AND. LS .LE. NC)THEN
	        IST=1; IEND=MIN(NR,SM_NRAY)
	        IF(BACK_SIDE .AND. SM_NRAY .LE. NR)THEN
	          BACK_SIDE=.FALSE.
	          PAR_FLUX=0.0_LDP
	        ELSE IF(BACK_SIDE)THEN
	          IST=NR+1; IEND=SM_NRAY
	        END IF
	      END IF
!
! Compute the slopes in each interval.
!
	      DO I=IST,IEND-1
                S(I)=(CHI_VEC(I)-CHI_VEC(I+1))/dZ(I)
	      END DO
!
! Now compute the derivatives of CHI at node I.
!
	      dS(IST)=S(IST) +(S(IST)-S(IST+1))*DZ(IST)/(DZ(IST)+DZ(IST+1))
	      DO I=IST+1,IEND-1
	        dS(I)=(S(I-1)*DZ(I)+S(I)*DZ(I-1))/
	1                  (DZ(I-1)+DZ(I))
	      END DO
	      dS(IEND)=S(IEND-1)+
	1                  (S(IEND-1)-S(IEND-2))*DZ(IEND-1)/
	1                  (DZ(IEND-2)+DZ(IEND-1))
!
! Adjust the first derivatives so that function is monotonic in each interval.
!
	      dS(IST)=( SIGN(ONE,S(IST))+SIGN(ONE,dS(IST)) )*
	1                      MIN(ABS(S(IST)),0.5_LDP*ABS(dS(IST)))
	      DO I=IST+1,IEND-1
	        dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	      END DO
	      dS(IEND)=( SIGN(ONE,S(IEND-1))+SIGN(ONE,dS(IEND)) )*
	1            MIN(ABS(S(IEND-1)),0.5_LDP*ABS(dS(IEND)))
!
! Compute the revised optical depth scale.
!
	    DO I=IST,IEND-1
	      DTAU(I)=DTAU(I)+DZSQ_ON_12(I)*(dS(I+1)-dS(I))
	    END DO
!
!
!
! We now use ETA_VEC for the source function.
!
	    T1=MINVAL(CHI_VEC(IST:IEND))
	    IF(T1 .GT. 0 .AND. INT_METHOD .EQ. 'STAU')THEN
!	      CALL TUNE(1,'ST FLUX')
	      DO I=IST,IEND
                ETA_VEC(I)=ETA_VEC(I)/CHI_VEC(I)
	      END DO
	      DO I=IST,IEND-1
                S(I)=(ETA_VEC(I+1)-ETA_VEC(I))/DTAU(I)
	      END DO
!
! Compute the functions used to evaluate the weighted integral over the
! polynomial fit to the source function.
!
! NB:  En(I)= EXP(-DTAU) {Integral[0 to DTAU] t^n EXP(t) dt }/ DTAU^n
!
	      DO I=IST,IEND-1
	        T1=DTAU(I)
	        IF(T1 .GT. 0.5_LDP)THEN
	          EE(I)=EXP(-T1)
	          E0(I)=1.0_LDP-EE(I)
	          E1(I)=1.0_LDP-E0(I)/T1
	          E2(I)=1.0_LDP-2.0_LDP*E1(I)/T1
	          E3(I)=1.0_LDP-3.0_LDP*E2(I)/T1
	        ELSE IF(T1 .GT. 0.1_LDP)THEN
	          E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1               (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1               (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP*
	1               (1.0_LDP-T1/10.0_LDP*(1.0_LDP-T1/11.0_LDP*
	1               (1.0_LDP-T1/12.0_LDP*(1.0_LDP-T1/13.0_LDP)))))))) )
	          E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	          E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	          E0(I)=T1*( 1.0_LDP-E1(I) )
	          EE(I)=1.0_LDP-E0(I)
	        ELSE
	          E3(I)=0.25_LDP*T1*( 1.0_LDP-0.20_LDP*T1*
	1               (1.0_LDP-T1/6.0_LDP*(1.0_LDP-T1/7.0_LDP*
	1               (1.0_LDP-T1/8.0_LDP*(1.0_LDP-T1/9.0_LDP) ))))
	          E2(I)=T1*( 1.0_LDP-E3(I) )/3.0_LDP
	          E1(I)=T1*( 1.0_LDP-E2(I) )/2.0_LDP
	          E0(I)=T1*( 1.0_LDP-E1(I) )
	          EE(I)=1.0_LDP-E0(I)
	        END IF
	      END DO
!
              DO I=IST,IEND-1
	        A0(I)=EE(I)
	        A1(I)=E0(I)-3.0_LDP*E2(I)+2.0_LDP*E3(I)
	        A2(I)=3.0_LDP*E2(I)-2.0_LDP*E3(I)
	        A3(I)=DTAU(I)*(E1(I)-2.0_LDP*E2(I)+E3(I))
	        A4(I)=DTAU(I)*(E3(I)-E2(I))
	      END DO
!
! Now compute the derivatives at node I.
!
	      dS(IST)=S(IST) +(S(IST)-S(IST+1))*DTAU(IST)/(DTAU(IST)+DTAU(IST+1))
	      DO I=IST+1,IEND-1
	        dS(I)=(S(I-1)*DTAU(I)+S(I)*DTAU(I-1))/
	1                  (DTAU(I-1)+DTAU(I))
	      END DO
	      dS(IEND)=S(IEND-1)+
	1                  (S(IEND-1)-S(IEND-2))*DTAU(IEND-1)/
	1                  (DTAU(IEND-2)+DTAU(IEND-1))
!
! Adjust the first derivatives so that function is monotonic in each interval.
!
	      dS(IST)=( SIGN(ONE,S(IST))+SIGN(ONE,dS(IST)) )*
	1                      MIN(ABS(S(IST)),0.5_LDP*ABS(dS(IST)))
	      DO I=IST+1,IEND-1
	        dS(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5_LDP*ABS(dS(I)))
	      END DO
	      dS(IEND)=( SIGN(ONE,S(IEND-1))+SIGN(ONE,dS(IEND)) )*
	1            MIN(ABS(S(IEND-1)),0.5_LDP*ABS(dS(IEND)))
!
	      IF(LS .GT. NC .OR. BACK_SIDE)THEN
	        PAR_FLUX=0.0_LDP			!No incident radiation
	      ELSE IF(LS .LE. NC .AND. HOLLOW_CORE)THEN
!
! PAR_FLUX was computed for this option of the previous loop.
! NB: ETA_VEC is actually the source function (done to save space).
!
	      ELSE
	        PAR_FLUX=ETA_VEC(SM_NRAY)	!Ray is thick!
	      END IF
	      DO I=IEND-1,IST,-1
	        PAR_FLUX=PAR_FLUX*A0(I)+ (
	1              ETA_VEC(I+1)*A1(I)
	1        +     ETA_VEC(I)*A2(I)
	1        -     dS(I+1)*A3(I)
	1        -     dS(I)*A4(I) )
	      END DO
!	      CALL TUNE(2,'ST FLUX')
	    ELSE
!
	      TAU(IST)=CHI_VEC(IST)*HALF_DZ(IST+1)
	      DO I=IST,IEND-1
	        TAU(I+1)=TAU(I)+DTAU(I)
	      END DO
!
! Compute the integrand = ETA*EXP(-TAU)
!
!	      CALL TUNE(1,'FLUX INTEG')
	      DO I=IST,IEND
	        ETA_VEC(I)=ETA_VEC(I)*EXP(-TAU(I))
	      END DO
!
! Use the Euler-Mclaurin Summation formula to compute the integrand.
! NB: We integrate with Z as the independent variable, since TAU may not
! be monotonic.
!
	      IF(HOLLOW_CORE .AND. LS .LE. NC .AND. .NOT. BACK_SIDE)THEN
	        PAR_FLUX=PAR_FLUX*EXP(-TAU(IEND))
	      ELSE
	        PAR_FLUX=0.0_LDP
	      END IF
	      T2=(ETA_VEC(IST)-ETA_VEC(IST+1))/dZ(IST)
	      DO I=IST+1,IEND-1
	        T1=T2
	        T2=(ETA_VEC(I-1)-ETA_VEC(I+1))*RECIP_DEL_Z(I)
	        T3=HALF_DZ(I-1)*(ETA_VEC(I-1)+ETA_VEC(I))
	        T4=DZSQ_ON_12(I-1)*(T2-T1)
	        IF(ABS(T4) .GT. 0.9_LDP*T3)T4=SIGN(0.9_LDP*T3,T4)
	        PAR_FLUX=PAR_FLUX + T3 + T4
	      END DO
	      I=IEND
	      T1=T2
	      T2=(ETA_VEC(I-1)-ETA_VEC(I))/dZ(I-1)
	      T3=HALF_DZ(I-1)*(ETA_VEC(I-1)+ETA_VEC(I))
	      T4=DZSQ_ON_12(I-1)*(T2-T1)
	      IF(ABS(T4) .GT. 0.9_LDP*T3)T4=SIGN(0.9_LDP*T3,T4)
	      PAR_FLUX=PAR_FLUX + (T3+T4)
!
!	      CALL TUNE(2,'FLUX INTEG')
	    END IF
!	
! Integrate this solution on to line profile. If we have a hollow core,
! we update the observed profile for the NEAR-SIDE integration --- this
! will also contain the back-ide contribution corrected for optical depth
! effects in the near side envelope.
!
	    IF(BACK_SIDE)THEN
	      BACK_SIDE=.FALSE.
	      GOTO 35000
	    ELSE
	      OBS_FLUX(ML)=OBS_FLUX(ML)+HQW_AT_RMAX(LS)*PAR_FLUX
	      IF(WRITE_IP)IP_OBS(LS,ML)=PAR_FLUX
	    END IF
!
40000	  CONTINUE
!$OMP END DO
!$OMP END PARALLEL
!
	  CALL TUNE(2,'ML')
45000	  CONTINUE
!
	  WRITE(LUER,'(A,I3,A,I5)')
	1            ' LS loop',LS,' is finished. Number of points along ray=',NRAY
50000	CONTINUE
!
! Convert the fluxes to Janskies for an object at 1kpc.
! The constant is 1E23*2*pi*1E20/(3.0856E21**2).
!     :the factor 10^23 is the conversion to Janskies.
!     :the factor 2pi is from the integration over a circular annulus.
!     :the factor 10^20 arises because R and P are in units of 10^10 cm.
!     :The factor 3.0856E21 is 1 kpc.
!
! NB: HQW_AT_RMAX is MU dMU and that PdP is R^2 MU dMU.
!
	OBS_FLUX=OBS_FLUX*6.59934_LDP*R(1)*R(1)          !Jansky's (1kpc)
!
! Output IP_OBS, if saved.
!
	IF(WRITE_IP)THEN
	  CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	  ACCESS_F=5
	  I=WORD_SIZE*(NP+1)/UNIT_SIZE; J=82
	  CALL OPEN_DIR_ACC_V1(NP,I,'10-Apr-2007','IP_DATA',J)
	  WRITE(82,REC=3)ACCESS_F,NOS,NP
	  WRITE(82,REC=ACCESS_F)(P(I),I=1,NP)
	  WRITE(82,REC=ACCESS_F+1)(MU_AT_RMAX(I),I=1,NP)
	  DO ML=1,NOS
	    WRITE(82,REC=ACCESS_F+1+ML)(IP_OBS(I,ML),I=1,NP),OBS_FREQ(ML)
	  END DO
	  CLOSE(UNIT=82)
	END IF
!
	CALL TUNE(2,'OBS_FRAME_SUB')
!
	RETURN
	END
