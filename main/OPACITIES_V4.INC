C
C Compute opacity and emissivity. This is a general include file
C provided program uses exactly the same variables. Can be achieved
C by copying declaration statements from CMFGEN. Always advisable
C to use ``IMPLICIT NONE''.
C
	IF(COMPUTE_NEW_CROSS)THEN
C
C Compute EXP(-hv/kT) and zero CHI, and ETA.
C
	  T1=-HDKT*CONT_FREQ
	  DO I=1,ND
	    EMHNUKT(I)=EXP(T1/T(I))
	    CHI(I)=0.0D0
	    ETA(I)=0.0D0
	  END DO
C
C Compute continuum intensity incident from the core assuming a TSTAR
C blackbody.
C
	  IC=TWOHCSQ*(CONT_FREQ**3)/(EXP(HDKT*CONT_FREQ/TSTAR)-1.0D0)
C
C Compute opacity and emissivity. ESOPAC must be call first since
C CHI is effectively zeroed in that routine.
C
	  CALL ESOPAC(ESEC,ED,ND)		!Electron scattering emission factor.
	  CALL ESOPAC(CHI,ED,ND)		!Electron scattering opacity.
C
C Free-free and bound-free opacities.
C
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        PHOT_ID=J
	        CALL GENOPAETA_V7(ID,CHI,ETA,CONT_FREQ,
	1           ATM(ID)%XzV_F,  ATM(ID)%XzVLTE_F, ATM(ID)%EDGEXzV_F,
	1           ATM(ID)%ZXzV,   ATM(ID)%NXzV_F,
	1           ATM(ID+1)%XzV,  ATM(ID+1)%XzVLTE, ATM(ID+1)%NXzV, 
	1           PHOT_ID,        ATM(ID)%XzV_ION_LEV_ID(J),
	1           ED,T,EMHNUKT,L_TRUE,ND,LST_DEPTH_ONLY)
	      END DO
	    END IF
	  END DO
C
C 
C
C Add in 2-photon emissivity and opacity.
C
	  CALL TWO_PHOT_OPAC(ETA,CHI,POPS,T,CONT_FREQ,ND,NT)
C
C Compute X-ray opacities and emissivities (at local electron temp).
C
	  IF(XRAYS)THEN
	    INCLUDE 'XOPAC_V4.INC'
	  END IF
C
	  CHI_C_EVAL(:)=CHI(:)
	  ETA_C_EVAL(:)=ETA(:)
C
	END IF
C
C 
C
C Section to revise continuum opacities etc so that they are computed at
C the correct frequency. We have stored the orginal continuum opacity and
C emissivity in CHI_C_EVAL and ETA_C_EVAL, which were computed at CONT_FREQ.
C
	IF(FL .NE. CONT_FREQ)THEN
C
	  DO J=1,ND
	    EMHNUKT_CONT(J)=EXP(-HDKT*CONT_FREQ/T(J))
	    EMHNUKT(J)=EXP(-HDKT*FL/T(J))
	  END DO
C
C Compute continuum intensity incident from the core assuming a TSTAR
C blackbody.
C
	  IC=TWOHCSQ*(FL**3)/(EXP(HDKT*FL/TSTAR)-1.0D0)
C
C We assume that the photoionization cross-section has not changed since the
C last iteration. Using the result that the stimulated emission occurs in
C LTE and is given by
C                     ETA/(2hv^3/c^2)
C we can adjust CHI and ETA so that the condition of constant photoionization
C cross-section is met. This adjustment automatically ensures that ETA/CHI 
C gives the Planck function in LTE. 
C
	  T1=(FL/CONT_FREQ)**3
	  T2=TWOHCSQ*(CONT_FREQ**3)
	  T3=TWOHCSQ*(FL**3)
	  DO J=1,ND
	    T4=ETA_C_EVAL(J)*T1*EXP(-HDKT*(FL-CONT_FREQ)/T(J))
	    CHI(J)=CHI_C_EVAL(J)+(ETA_C_EVAL(J)/T2-T4/T3)
	    ETA(J)=T4
	  END DO
	ELSE
C
C We reset CHI and ETA in case shock X-ray emission has been added to ETA.
C
	  CHI(1:ND)=CHI_C_EVAL(1:ND)
	  ETA(1:ND)=ETA_C_EVAL(1:ND)
	END IF
C
C 
C
C The shock emission is added separately since it does not occur at the
C local electron temperature.
C
	IF(XRAYS)THEN
!i
	  IF(FF_XRAYS)THEN
!
! Since T_SHOCK is depth indpendent, Z^2 * (the free-free Gaunt factors) 
! are depth independent.
!
! We use T3 for the Electron density. We asume H, He, and C are fully ionized
! in the X-ray emitting plasma. All other species are assumed to have Z=6.0
!
	    IF(T_SHOCK_1 .NE. 0.0D0)THEN
	      T1=CHIFF*TWOHCSQ*(FILL_FAC_XRAYS_1**2)*
	1                    EXP(-HDKT*CONT_FREQ/T_SHOCK_1)/SQRT(T_SHOCK_1)
	      T2=1.0D0 ; TA(1)=GFF(CONT_FREQ,T_SHOCK_1,T2)
	      T2=2.0D0 ; TA(2)=4.0D0*GFF(CONT_FREQ,T_SHOCK_1,T2)
	      T2=6.0D0 ; TA(3)=36.0D0*GFF(CONT_FREQ,T_SHOCK_1,T2)
	      DO I=1,ND
	        T2=TA(1)*POP_SPECIES(I,1)+TA(2)*POP_SPECIES(I,2) +
	1         TA(3)*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        T3=POP_SPECIES(I,1)+POP_SPECIES(I,2) +
	1          6.0D0*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        ZETA(I)=T1*T2*T3*EXP(-V_SHOCK_1/V(I))	 !Zeta is temporary
              END DO
	    END IF
	    IF(T_SHOCK_2 .NE. 0.0D0)THEN
	      T1=CHIFF*TWOHCSQ*(FILL_FAC_XRAYS_2**2)*
	1                    EXP(-HDKT*CONT_FREQ/T_SHOCK_2)/SQRT(T_SHOCK_2)
	      T2=1.0D0 ; TA(1)=GFF(CONT_FREQ,T_SHOCK_2,T2)
	      T2=2.0D0 ; TA(2)=4.0D0*GFF(CONT_FREQ,T_SHOCK_2,T2)
	      T2=6.0D0 ; TA(3)=36.0D0*GFF(CONT_FREQ,T_SHOCK_2,T2)
	      DO I=1,ND
	        T2=TA(1)*POP_SPECIES(I,1)+TA(2)*POP_SPECIES(I,2) +
	1       TA(3)*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        T3=POP_SPECIES(I,1)+POP_SPECIES(I,2) +
	1          6.0D0*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        ZETA(I)=T1*T2*T3*EXP(-V_SHOCK_2/V(I))	 !Zeta is temporary
              END DO
	    END IF
	  ELSE
!
! Use X-ray emission as tubulated by a PLASMA code. Emission should  be
! tabulated per electron and per ion.
!
	    CALL GET_SCL_XRAY_FLUXES_V1(CONT_FREQ,
	1              XRAY_EMISS_1,XRAY_EMISS_2,
	1              NU_EVAL_CONT,NCF,FREQ_INDX,
	1              VSMOOTH_XRAYS,SECTION)
!
! We use T3 for the Electron density. We asume H, He, and C are fully ionized
! in the X-ray emitting plasma. All other species are assumed have Z=6.0
!
	    DO I=1,ND
	      T1=EXP(-V_SHOCK_1/V(I))*(FILL_FAC_XRAYS_1)**2
	      T2=EXP(-V_SHOCK_2/V(I))*(FILL_FAC_XRAYS_2)**2
	      T3=POP_SPECIES(I,1)+2.0D0*POP_SPECIES(I,2)+
	1              6.0D0*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	      ZETA(I)=(T1*XRAY_EMISS_1+T2*XRAY_EMISS_2)*T3*POP_ATOM(I)
	    END DO
	  END IF
!
	  IF(XRAY_SMOOTH_WIND)ZETA(1:ND)=ZETA(1:ND)*CLUMP_FAC(1:ND)
          ETA(1:ND)=ETA(1:ND)+ZETA(1:ND)
!
	  T1=0.241838D0		!eV to 10^15Hz
	  IF(SECTION .EQ. 'CONTINUUM')THEN
	    IF(ML .EQ. 1)THEN
	       XRAY_LUM_TOT(1:ND)=0.0D0
	       XRAY_LUM_0P1(1:ND)=0.0D0
	       XRAY_LUM_1KEV(1:ND)=0.0D0
	    END IF
	    XRAY_LUM_TOT(1:ND)=XRAY_LUM_TOT(1:ND)+ZETA(1:ND)*FQW(ML)
	    IF(NU(ML) .GE. 100.0D0*T1)XRAY_LUM_0P1(1:ND)=
	1                XRAY_LUM_0P1(1:ND)+ZETA(1:ND)*FQW(ML)
	    IF(NU(ML) .GE. 1000.0D0*T1)XRAY_LUM_1KEV(1:ND)=
	1                XRAY_LUM_1KEV(1:ND)+ZETA(1:ND)*FQW(ML)
	  END IF
	END IF
C
C The continuum source function is defined by:
C                                              S= ZETA + THETA.J
	DO I=1,ND
	  ZETA(I)=ETA(I)/CHI(I)
	  THETA(I)=ESEC(I)/CHI(I)
	END DO
C
C Store TOTAL continuum line emissivity and opacity.
C
	ETA_CONT(:)=ETA(:)
	CHI_CONT(:)=CHI(:)
