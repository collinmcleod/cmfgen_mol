C
C Compute the opacity AND emissivity variation as a function of the changes
C in population levels. As the variation is linear, they can be added
C independently. Note that the VCHI and VETA are not zeroed when the
C variation routines are called.
C
	IF(COMPUTE_NEW_CROSS)THEN
C
C Zero VCHI and VETA, then allow for the variation of the electron
C scattering opacity and emissivity with ED.
C
C TA is used as a work vector.
C
	  CALL DP_ZERO(VCHI,NT*ND)
	  CALL DP_ZERO(VETA,NT*ND)
	  CALL DP_ZERO(VCHI_ALL,NT*ND)
	  CALL DP_ZERO(VETA_ALL,NT*ND)
C
C We use EMHNUKT_CONT since we are evaluating the variation in the opacities
C and emissivities at CONT_FREQ, not FL
C
	  DO L=1,ND
	    EMHNUKT_CONT(L)=EXP(-HDKT*CONT_FREQ/T(L))
	  END DO
C
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        PHOT_ID=J
	        CALL VAR_OP_V8(VCHI,VETA,VCHI_ALL,VETA_ALL,
	1         ATM(ID)%XzV, ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT, ATM(ID)%NXzV,
	1         ATM(ID)%XzVLTE_F, ATM(ID)%EDGEXzV_F, ATM(ID)%NXzV_F, ATM(ID)%F_TO_S_XzV,
	1         ATM(ID+1)%XzV,    ATM(ID+1)%XzVLTE,  ATM(ID+1)%dlnXzVLTE_dlnT,
	1         ATM(ID+1)%NXzV,   PHOT_ID,           ATM(ID)%XzV_ION_LEV_ID(J),
	1         ED,T,EMHNUKT_CONT,IMP_VAR,
	1         CONT_FREQ,ATM(ID)%ZXzV,ID,L_TRUE,
	1         ATM(ID)%EQXzV,ATM(ID+1)%EQXzV,NT,ND,LST_DEPTH_ONLY)
	      END DO
	    END IF
	 END DO
C $END_ION$
C $END_EXPAND$
C
C 
C
C Add in 2-photon opacities and emissivities.
C
	  CALL TWO_PHOT_VAR_OPAC(VETA,VCHI,POPS,T,CONT_FREQ,ND,NT)
C
	  IF(XRAYS)THEN
	    INCLUDE 'VAR_XOPACH_V4.INC'
	  END IF
C
	  VCHI_SAV(:,:)=VCHI(:,:)
	  VETA_SAV(:,:)=VETA(:,:)
!
! _ALL will contain the variation of ETA/CHI for ALL species.
!
	  VCHI_ALL(:,:)=VCHI_ALL(:,:)+VCHI(:,:)
	  VETA_ALL(:,:)=VETA_ALL(:,:)+VETA(:,:)
C
	  VCHI_ALL_SAV(:,:)=VCHI_ALL(:,:)
	  VETA_ALL_SAV(:,:)=VETA_ALL(:,:)
C
	END IF
C
C 
C
C This section of code can be utilized by both the DTDR and CONTINUUM code
C sections.
C
	IF(CONT_FREQ .NE. FL)THEN
C
C Scale dCHI and dETA to allow for the slight variation in frequency.
C Correction should generally be small.
C
! NB: Previously T4 was defined by T4= EMHNUKT(L)/EMHNUKT_CONT(L). In the
! presence of X-rays, and with NU=1000 or larger, EMHNUKT_CONT could be
! zero, thus causing a divide by zero.
!
	  T1=(FL/CONT_FREQ)**3
	  T2=TWOHCSQ*(CONT_FREQ**3)
	  T3=TWOHCSQ*(FL**3)
	  J=1
	  IF(LST_DEPTH_ONLY)J=ND
!
	  DO L=J,ND
	    T4=T1*EXP(-HDKT*(FL-CONT_FREQ)/T(L))
            DO K=1,NT
              VETA(K,L)=VETA_SAV(K,L)*T4
              VETA_ALL(K,L)=VETA_ALL_SAV(K,L)*T4
            END DO
C
C Need to correct VETA for the variation in T in the factor T4.
C Over correction at present because impurity species included.
C
	    VETA(NT,L)=VETA(NT,L) + 
	1                ETA_C_EVAL(L)*T4*HDKT*(FL-CONT_FREQ)/T(L)/T(L)
	    VETA_ALL(NT,L)=VETA_ALL(NT,L) + 
	1                ETA_C_EVAL(L)*T4*HDKT*(FL-CONT_FREQ)/T(L)/T(L)
C
C Can now correct the opacity.
C
	   DO K=1,NT
              VCHI(K,L)=VCHI_SAV(K,L)+(VETA_SAV(K,L)/T2-VETA(K,L)/T3)
              VCHI_ALL(K,L)=VCHI_ALL_SAV(K,L)+(VETA_ALL_SAV(K,L)/T2-VETA_ALL(K,L)/T3)
	   END DO
	    
          END DO
C
	END IF			!Correct cross-section?
C
C 
C
C
C In the CONTINUUM code section (with blanketing) the variation of the opacity 
C and emissivity due to the electron scattering term is handled separately.
C It is required in the DTDR section, DIELECTRONIC and LINE sections.
C In DTDR the cross-sections can also be held fixed. In order that our
C correction process implemented above works, we must keep the electron
C scattering cross-section separate (primarily in the emissivity).
C
	IF(SECTION .NE. 'CONTINUUM')THEN
	  DO K=1,ND
	    VCHI(NT-1,K)=VCHI(NT-1,K)+ESEC(K)/ED(K)
	    VETA(NT-1,K)=VETA(NT-1,K)+ESEC(K)*RJ(K)/ED(K)
	  END DO
	END IF
 
