C**************************************************************************
C
C Include file to calculate the variation of J as a function of
C the population values. For use with the Sobolev approximation,
C and the Dielectronic recombination section.
C
C Altered 11-Dec-1997 - Section for computing J in CONT_VEL section  with
C                         additional data points inserted. Designed to give
C                         higher accuracy, especially near ionization fronts.
C                         Additional points are inserted for ALL frequencies.
C                         For this option, ACCURATE must be set to TRUE.
C
C Altered 12-Jun-91 - Call to PERTJFEAUNEW replaced by call to
C                     PERTJFEAU_IBC. This has improved handling of
C                     the outer boundary condition.
C Altered 13-Feb-88 - Changed to allow Banded-Newton Raphson Operators.
C                     The number of bands is given by NUM_BNDS and must
C                     be odd. If NUM_BNDS=ND, the matrix is assumed
C                     to represent the full Newton-Raphson Operator.
C                     The index, DIAG_INDX must also be specified, but is
C                     only used for the BANDED case. Note 
C                     DIAG_INDX=( NUM_BNDS/2+1 )
C Altered 6-May-86 -  PERTJFEAU installed. (now varcontfeau)
C Created 26-Jan-88 - (Final version).
C
C BA and VJ have the following dimensions.
C
C                     BA(NT,NT,NUM_BNDS,ND)
C                     VJ(NT,NUM_BNDS,ND)
C
C                     In addition, the variables HYD_PRES, HEI_PRES
C                     and HE2_PRES have been installed. These have
C                     been included so the routine is completely general.
C                     These LOGICAL variables should be specified by
C                     parameter statements - they can then be optimized
C                     out of the code at compile time.
C
C Altered 21-Mar-1989 - VETA routine installed. VETA and VCHI are now
C                       computed simultaneously.
C
C Altered 04-Apr-89 - EDDINGTON logical variable installed.
C                     THIS_FREQ_EXT restored (was THIS_FREQ).
C                     ACCURATE variable removed.
C
C Altered 31-May-1989 - THICK variable replaced by THK_CONT.
C
C 
C**************************************************************************
C
C Solve for the perturbations to J in terms of the perturbations
C to CHI and ETA. F2DA is the dCHI matrix ; FC the dETA matrix
C and FA is the d(diffusion approx) vector for the boundary
C condition at the core.
C
	  IF(THIS_FREQ_EXT .AND. .NOT. CONT_VEL)THEN
	    CALL TUNE(1,'DJFEAUEXT')
	      CALL PERTJFEAU_IBC(F2DAEXT,FCEXT,FAEXT,
	1            DTAU,CHIEXT,REXT,ZETAEXT,
	1            THETAEXT,RJEXT,QEXT,FEXT,dCHIdR,
	1            TA,TB,TC,HBC_J,HBC_S,INBC,DBB,DIF,
	1            THK_CONT,NDEXT,METHOD)
C
C Put variation matrices on old grid. Note that FA is used for the
C diffusion approximation.
C
	      CALL REGRID_dCHI(F2DA,CHI,ND,POS_IN_NEW_GRID,
	1                          F2DAEXT,CHIEXT,NDEXT,COEF,INDX)
	      CALL REGRID_dCHI(FC,ETA,ND,POS_IN_NEW_GRID,
	1                          FCEXT,ETAEXT,NDEXT,COEF,INDX)
	      DO I=1,ND
	        FA(I)=FAEXT(POS_IN_NEW_GRID(I))
	      END DO
	    CALL TUNE(2,'DJFEAUEXT')
C
C 
C
	  ELSE IF(CONT_VEL .AND. .NOT. ACCURATE)THEN
	    IF(FIRST_FREQ)THEN
	      TX(:,:,:)=0.0D0
	      TVX(:,:,:)=0.0D0
	    ELSE
	      RAT_TOO_BIG=.FALSE.
	      DO L=1,ND
	        TA(L)=(CHI_PREV(L)-ESEC(L))/(CHI_CONT(L)-ESEC(L))
	        IF(ETA_CONT(L) .EQ. 0)THEN
                  TB(L)=1.0D0
	        ELSE
	          TB(L)=ETA_PREV(L)/ETA_CONT(L)
	        END IF
	        IF(TA(L) .GT. 1.5)RAT_TOO_BIG=.TRUE.
	      END DO
	      IF(RAT_TOO_BIG)THEN
	        DO L=1,ND
	          TA(L)=0.0D0
	          TB(L)=0.0D0
	        END DO
	      END IF
	      DO J=1,ND
	        DO K=1,ND
	          TX(K,J,3)=TX(K,J,3) * TA(J)
	          TX(K,J,4)=TX(K,J,4) * TB(J)
	        END DO
	      END DO
	      DO J=1,ND
	        DO K=1,ND-1
	          TVX(K,J,3)=TVX(K,J,3) * TA(J)
	          TVX(K,J,4)=TVX(K,J,4) * TB(J)
	        END DO
	      END DO
	    END IF
	    DO I=1,NM
	      DO_THIS_TX_MATRIX(I)=.TRUE.
	    END DO
	    DO I=TX_OFFSET+1,NM
	      IF(VAR_IN_USE_CNT(I) .EQ. 0)THEN
	        DO_THIS_TX_MATRIX(I)=.FALSE.
	      END IF
	    END DO
C
C Use TA as temporary storage for the emissivity.
C
	    IF(COHERENT_ES)THEN
	      TA(1:ND)=ETA_CLUMP(1:ND)
	      ES_COH_VEC(1:ND)=ESEC_CLUMP(1:ND)/CHI_CLUMP(1:ND)
	    ELSE
C
C Two scenarios: 
C    (i) We use a lambda iteration to allow for the variation of J in
C           the electron scattering term. ES_COH_VEC (==THETA) is zero,
C           and the electron scattering emissivity (related to RJ_ES)
C           is included with the emissivity directly.
C   (ii) We assume that the variation of J in the electron scattering
C           term can be treated coherently. ES_COH_VEC is non-zero.
C           A correction to the emissivity is made to allow for the
C           difference between RJ and RJ_ES.
C
C	      ES_COH_VEC(1:ND)=0.0D0
C	      IF(MIXED_ES_VAR)THEN
C	        DO I=ND,1,-1
C	          IF( ABS(RJ(I)-RJ_ES(I))/(RJ(I)+RJ_ES(I)) .GT.
C	1                          0.5*ES_VAR_FAC)EXIT
C	          ES_COH_VEC(I)=ESEC_CLUMP(I)/CHI_CLUMP(I)
C	        END DO
C	      END IF
C	      DO I=1,ND
C	        IF(ES_COH_VEC(I) .EQ. 0)THEN
C	          TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*RJ_ES(I)
C	        ELSE
C	          TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*(RJ_ES(I)-RJ(I))
C	        END IF
C	      END DO
C	    END IF
C
	      IF(MIXED_ES_VAR)THEN
	        T1=2.0D0
	        DO I=1,ND
	          IF(RJ_ES(I) .GT. RJ(I))THEN
	            TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*
	1                   (RJ_ES(I)-RJ(I)/T1)
	            ES_COH_VEC(I)=ESEC_CLUMP(I)/CHI_CLUMP(I)/T1
	          ELSE
	            TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*RJ_ES(I)/T1
	            ES_COH_VEC(I)=ESEC_CLUMP(I)/CHI_CLUMP(I)*
	1                    (RJ_ES(I)/RJ(I))/T1
	          END IF
	        END DO
	      ELSE
	        ES_COH_VEC(1:ND)=0.0D0
	        TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND)
	      END IF
	   END IF
C
	   CALL VAR_MOM_J_CMF_V7(TA,CHI_CLUMP,ESEC_CLUMP,
	1                  ES_COH_VEC,V,SIGMA,R,
	1                  TX,TVX,dJ_DIF_d_T,dJ_DIF_d_dTdR, 
	1                  dRSQH_DIF_d_T,dRSQH_DIF_d_dTdR, 
	1                  KI,WM,RHS_dHdCHI,
	1                  FEDD,GEDD,N_ON_J,HBC_CMF,INBC,NBC_CMF,
	1                  FEDD_PREV,GEDD_PREV,N_ON_J_PREV,
	1                  HBC_PREV,INBC_PREV,NBC_PREV,
	1                  FIRST_FREQ,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1                  DO_THIS_TX_MATRIX,METHOD,ND,NM,NM_KI)
C
C Correcting for clumping this way does it for both the continuum and lines.
C
	    IF(DO_CLUMP_MODEL)THEN
	      DO J=1,ND
	        DO K=1,ND
	          TX(K,J,1)=TX(K,J,1)*CLUMP_FAC(J)
	          TX(K,J,2)=TX(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
	      DO J=1,ND
	        DO K=1,ND-1
	          TVX(K,J,1)=TVX(K,J,1)*CLUMP_FAC(J)
	          TVX(K,J,2)=TVX(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
	    END IF
C
C We use TB as a temporary vector for J in the electron scattering emissivity.
C Its value depends on whether we have coherent or incoherent e.s.
C
	    IF(COHERENT_ES)THEN
	      TB(1:ND)=RJ(1:ND)
	    ELSE
	      TB(1:ND)=RJ_ES(1:ND)
	    END IF
	    DO J=1,ND
	      DO K=1,ND
	        TX(K,J,3)=TX(K,J,3) + TX(K,J,1)
	        TX(K,J,4)=TX(K,J,4) + TX(K,J,2)
	        TX(K,J,5)=TX(K,J,5) + TX(K,J,1) + TX(K,J,2)*TB(J)
	      END DO
	    END DO
	    DO J=1,ND
	      DO K=1,ND-1
	        TVX(K,J,3)=TVX(K,J,3) + TVX(K,J,1)
	        TVX(K,J,4)=TVX(K,J,4) + TVX(K,J,2)
	        TVX(K,J,5)=TVX(K,J,5) + TVX(K,J,1) + TVX(K,J,2)*TB(J)
	      END DO
	    END DO
C
C Update line variation matrices. Note that the matrices now refer to the
C variation with respect to levels (e.g. the lower and upper level) and
C not CHIL and ETAL. This should save arrays when we are dealing with
C overlapping transitions between the same level(s).
C
C For simplicity we have ignored the T dependance of L_STAR_RATIO and
C U_STAR_RATIO.
C
	    DO SIM_INDX=1,MAX_SIM
	      NL=LOW_POINTER(SIM_INDX)
	      NUP=UP_POINTER(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
	        DO J=1,ND
	          OPAC_FAC=lINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      U_STAR_RATIO(J,SIM_INDX)
	          STIM_FAC=LINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          DO K=1,ND
	            TX(K,J,NL)=TX(K,J,NL) + 
	1             OPAC_FAC*TX(K,J,1)
	            TX(K,J,NUP)=TX(K,J,NUP) +
	1             ( EMIS_FAC*TX(K,J,2) -
	1               STIM_FAC*TX(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
C
C We now do the update for dH (i.e. TVX)
C
	    DO SIM_INDX=1,MAX_SIM
	      NL=LOW_POINTER(SIM_INDX)
	      NUP=UP_POINTER(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
	        DO J=1,ND
	          OPAC_FAC=lINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      U_STAR_RATIO(J,SIM_INDX)
	          STIM_FAC=LINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          DO K=1,ND-1
	            TVX(K,J,NL)=TVX(K,J,NL) + 
	1             OPAC_FAC*TVX(K,J,1)
	            TVX(K,J,NUP)=TVX(K,J,NUP) +
	1             ( EMIS_FAC*TVX(K,J,2) -
	1               STIM_FAC*TVX(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
C
C Now zero dCHI and dETA storage locations. These refer to the TOTAL 
C opacity and emissivity.
C
	    TX(:,:,1:2)=0.0D0
	    TVX(:,:,1:2)=0.0D0
C
C 
C
	  ELSE IF(CONT_VEL .AND. ACCURATE)THEN
	    IF(FIRST_FREQ)THEN
	      TX_EXT(:,:,:)=0.0D0
	      TVX_EXT(:,:,:)=0.0D0
	    ELSE
	      RAT_TOO_BIG=.FALSE.
	      DO L=1,ND
	        TA(L)=(CHI_PREV(L)-ESEC(L))/(CHI_CONT(L)-ESEC(L))
	        TB(L)=ETA_PREV(L)/ETA_CONT(L)
	        IF(TA(L) .GT. 1.5)RAT_TOO_BIG=.TRUE.
	      END DO
	      IF(RAT_TOO_BIG)THEN
	        DO L=1,ND
	          TA(L)=0.0D0
	          TB(L)=0.0D0
	        END DO
	      END IF
	      DO J=1,ND
	        DO K=1,NDEXT
	          TX_EXT(K,J,3)=TX_EXT(K,J,3) * TA(J)
	          TX_EXT(K,J,4)=TX_EXT(K,J,4) * TB(J)
	        END DO
	      END DO
	      DO J=1,ND
	        DO K=1,NDEXT-1
	          TVX_EXT(K,J,3)=TVX_EXT(K,J,3) * TA(J)
	          TVX_EXT(K,J,4)=TVX_EXT(K,J,4) * TB(J)
	        END DO
	      END DO
	    END IF
	    DO I=1,NM
	      DO_THIS_TX_MATRIX(I)=.TRUE.
	    END DO
	    DO I=TX_OFFSET+1,NM
	      IF(VAR_IN_USE_CNT(I) .EQ. 0)THEN
	        DO_THIS_TX_MATRIX(I)=.FALSE.
	      END IF
	    END DO
C
C Use TA as temporary storage for the emissivity.
C
	    IF(COHERENT_ES)THEN
	      TA(1:NDEXT)=ETAEXT(1:NDEXT)
	    ELSE
	      TA(1:NDEXT)=ETAEXT(1:NDEXT)+
	1                 ESECEXT(1:NDEXT)*RJ_ES(1:NDEXT)
	    END IF
	    CALL VAR_MOM_JEXT_CMF_V2(TA,CHIEXT,ESECEXT,
	1                  VEXT,SIGMAEXT,REXT,
	1                  ETA_CLUMP,CHI_CLUMP,ESEC_CLUMP,
	1                  INDX,COEF,INTERP_TYPE,ND,
	1                  TX_EXT,TVX_EXT,
	1                  dJ_DIF_d_T_EXT,dJ_DIF_d_dTdR_EXT,
	1                  dRSQH_DIF_d_T,dRSQH_DIF_d_dTdR, 
	1                  KI,RHS_dHdCHI,
	1                  FEDD,GEDD,N_ON_J,HBC_CMF,INBC,NBC_CMF,
	1                  FEDD_PREV,GEDD_PREV,N_ON_J_PREV,
	1                  HBC_PREV,INBC_PREV,NBC_PREV,
	1                  FIRST_FREQ,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1                  DO_THIS_TX_MATRIX,METHOD,COHERENT_ES,NDEXT,NM,NM_KI)
C
C Correcting for clumping this way does it for both the continuum and lines.
C
	    IF(DO_CLUMP_MODEL)THEN
	      DO J=1,ND
	        DO K=1,NDEXT
	          TX_EXT(K,J,1)=TX_EXT(K,J,1)*CLUMP_FAC(J)
	          TX_EXT(K,J,2)=TX_EXT(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
	      DO J=1,ND
	        DO K=1,NDEXT-1
	          TVX_EXT(K,J,1)=TVX_EXT(K,J,1)*CLUMP_FAC(J)
	          TVX_EXT(K,J,2)=TVX_EXT(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
	    END IF
C
C We use TB as a temporary vector for J in the electron scattering emissivity.
C Its value depends on whether we have coherent or incoherent e.s.
C
	    IF(COHERENT_ES)THEN
	      TB(1:ND)=RJ(1:ND)
	    ELSE
	      TB(1:ND)=RJ_ES(1:ND)
	    END IF
	    DO J=1,ND
	      DO K=1,NDEXT
	        TX_EXT(K,J,3)=TX_EXT(K,J,3) + TX_EXT(K,J,1)
	        TX_EXT(K,J,4)=TX_EXT(K,J,4) + TX_EXT(K,J,2)
	        TX_EXT(K,J,5)=TX_EXT(K,J,5) + TX_EXT(K,J,1) + 
	1                                         TX_EXT(K,J,2)*TB(J)
	      END DO
	    END DO
	    DO J=1,ND
	      DO K=1,NDEXT-1
	        TVX_EXT(K,J,3)=TVX_EXT(K,J,3) + TVX_EXT(K,J,1)
	        TVX_EXT(K,J,4)=TVX_EXT(K,J,4) + TVX_EXT(K,J,2)
	        TVX_EXT(K,J,5)=TVX_EXT(K,J,5) + TVX_EXT(K,J,1) + 
	1                                        TVX_EXT(K,J,2)*TB(J)
	      END DO
	    END DO
C
C Update line variation matrices. Note that the matrices now refer to the
C variation with respect to levels (e.g. the lower and upper level) and
C not CHIL and ETAL. This should save arrays when we are dealing with
C overlapping transitions between the same level(s).
C
C For simplicity we have ignored the T dependance of L_STAR_RATIO and
C U_STAR_RATIO.
C
	    DO SIM_INDX=1,MAX_SIM
	      NL=LOW_POINTER(SIM_INDX)
	      NUP=UP_POINTER(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
	        DO J=1,ND
	          OPAC_FAC=lINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      U_STAR_RATIO(J,SIM_INDX)
	          STIM_FAC=LINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          DO K=1,NDEXT
	            TX_EXT(K,J,NL)=TX_EXT(K,J,NL) + 
	1             OPAC_FAC*TX_EXT(K,J,1)
	            TX_EXT(K,J,NUP)=TX_EXT(K,J,NUP) +
	1             ( EMIS_FAC*TX_EXT(K,J,2) -
	1               STIM_FAC*TX_EXT(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
C
C We now do the update for dH (i.e. TVX_EXT)
C
	    DO SIM_INDX=1,MAX_SIM
	      NL=LOW_POINTER(SIM_INDX)
	      NUP=UP_POINTER(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
	        DO J=1,ND
	          OPAC_FAC=lINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      U_STAR_RATIO(J,SIM_INDX)
	          STIM_FAC=LINE_OPAC_CON(SIM_INDX)*
	1                      LINE_PROF_SIM(SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*
	1                      U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          DO K=1,NDEXT-1
	            TVX_EXT(K,J,NL)=TVX_EXT(K,J,NL) + 
	1             OPAC_FAC*TVX_EXT(K,J,1)
	            TVX_EXT(K,J,NUP)=TVX_EXT(K,J,NUP) +
	1             ( EMIS_FAC*TVX_EXT(K,J,2) -
	1               STIM_FAC*TVX_EXT(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
C
C Now zero dCHI and dETA storage locations. These refer to the TOTAL
C opacity and emissivity.
C
	    TX_EXT(:,:,1:2)=0.0D0		!dCHI,dETA
	    TVX_EXT(:,:,1:2)=0.0D0		!dCHI,dETA
C
C Compute dJ/d? just on the small grid. We don't need dH/d? since the
C statistical and radiative equilibrium equations depend only on J.
C
	    DO I=3,NM
	      DO J=1,ND
	        DO K=1,ND
	          TX(K,J,I)=TX_EXT(POS_IN_NEW_GRID(K),J,I)
	        END DO
	      END DO
	    END DO
C
	    DO K=1,ND
	      dJ_DIF_d_T(K)=dJ_DIF_d_T_EXT(POS_IN_NEW_GRID(K))
	      dJ_DIF_d_dTdR(K)=dJ_DIF_d_dTdR_EXT(POS_IN_NEW_GRID(K))
	    END DO
C
C 
C
	  ELSE IF(EDDINGTON)THEN
	    CALL TUNE(1,'PERTJFEAU')
	      CALL PERTJFEAU_IBC(F2DA,FC,FA,
	1            DTAU,CHI,R,ZETA,
	1            THETA,RJ,QEDD,FEDD,dCHIdR,
	1            TA,TB,TC,HBC_J,HBC_S,INBC,DBB,DIF,
	1            THK_CONT,ND,METHOD)
	    CALL TUNE(2,'PERTJFEAU')
	  ELSE
	    CALL TUNE(1,'PERTJD')
	      CALL MULTVEC(SOURCE,ZETA,THETA,RJ,ND)
	      CALL NEWPERTJD(F2DA,FC,FA,FB,VK,WM,AQW
	1       ,DTAU,CHI,dCHIdR,R,Z,P,THETA,SOURCE,TA,TB,TC,XM
	1       ,DIF,DBB,IC,ESEC,THK_CONT,NC,ND,NP,METHOD)
	    CALL TUNE(2,'PERTJD')
	  END IF
C 
C
C Compute the opacity AND emissivity variation as a function of the changes
C in population levels. 
C
C
	INCLUDE 'VAROPAC_V4.INC'
C 
C
C Zero VJ array.
C
	CALL DP_ZERO(VJ,NT*NUM_BNDS*ND) 
C
	IF(CONT_VEL)THEN
C
C NB: We no longer include the variation ESEC variation with CHI, but treat it
C separately. this correction has now been specifically included in VAROPAC.
C
C
C The matrix TX gives dJ ( J depth, X depth, X) where X represent some 
C fundamental parameter such as CHI_C, ETA_C, ESEC, ETAL_ etc
C
C We now convert TX  to a smaller matrix, taking into account only the 
C variation of X at J's depth, and neighbouring depths.
C
C dJ_LOC( X , X depth [1=3:NUM_BNDS], J depth )
C
C We only treat Z_INDX from 3 onwards as the first refer to CHI and ETA which
C have been zeroed.
C
	  DO X_INDX=3,NM
	    IF(DO_THIS_TX_MATRIX(X_INDX))THEN
	      DO K=1,ND
	        DO J=BNDST(K),BNDEND(K)
	          L=BND_TO_FULL(J,K)
	          dJ_LOC(X_INDX,J,K)=TX(K,L,X_INDX)
	        END DO
	      END DO
	    END IF
	  END DO
C
C Compute VJ which gives the variation of J with respect to the atomic
C populations. NB: Electron scattering cross-section (6.65D-15) was replaced 
C by ESEC(L)/ED(L) 24-Sep-1997.
C
C NOPS= ND*NUM_BNDS*( 6NT + 2 + 7NUM_SIM )
C
	  CALL TUNE(1,'VC_VCHI')
	  DO K=1,ND
	    DO J=BNDST(K),BNDEND(K)
	      L=BND_TO_FULL(J,K)
	      DO I=1,NT
	        VJ(I,J,K)=VJ(I,J,K) +
	1              ( VCHI(I,L)*dJ_LOC(3,J,K) + VETA(I,L)*dJ_LOC(4,J,K) )
	      END DO
	      VJ(NT-1,J,K)=VJ(NT-1,J,K) + ESEC(L)*dj_LOC(5,J,K)/ED(L)
C                    
C Now must do line terms.
C
	      DO I=TX_OFFSET+1,NM
	        IF(VAR_IN_USE_CNT(I) .GT. 0)THEN
	          NL=VAR_LEV_ID(I)
	          VJ(NL,J,K)=VJ(NL,J,K) + dj_LOC(I,J,K)
	        END IF
	      END DO
	    END DO	!Over Variable depth (1:NUM_BNDS)
	  END DO		!Over J depth.
	  CALL TUNE(2,'VC_VCHI')
C
C	  DO L=1,ND
C	    VCHI(EQNE,L)=VCHI(EQNE,L)+ESEC(L)/ED(L)
C	    VETA(EQNE,L)=VETA(EQNE,L)+ESEC(L)/ED(L)*RJ(L)
C	  END DO
C
C Update VJ for perturbations in diffusion approximation. The case ND=NUM_BNDS
C is no longer treated.
C
	  IF(DIF)THEN
	    T1=DBB/DTDR
	    DO J=DIAG_INDX,NUM_BNDS
	      K=ND+DIAG_INDX-J
	      DO I=1,NT-1
	        VJ(I,J,K)=VJ(I,J,K) + dJ_DIF_d_dTdR(K)*DIFFW(I)
	      END DO  
	      VJ(NT,J,K)=VJ(NT,J,K) +
	1                       dJ_DIF_d_T(K)+dJ_DIF_d_dTdR(K)*DIFFW(NT)
	    END DO
	  END IF
	ELSE
C
C Section to increment the variation matrices.
C
	  DO K=1,ND			!Depth of intensity
	    DO J=BNDST(K),BNDEND(K)		!Depth of variable index
	      L=BND_TO_FULL(J,K)
	      DO I=1,NT
	        VJ(I,J,K)=VJ(I,J,K)+
	1          (F2DA(K,L)*VCHI(I,L)+FC(K,L)*VETA(I,L))
	      END DO
	    END DO
	  END DO
C
C Update VJ for perturbations in diffusion approximation.
C 18-Dec-1991 replaced ND in VJ( ,ND,K) by VJ( ,NUM_BNDS,K) to avoid
C             compilations errors when NUM_BNDS .NE. ND
C             (only in first clause)
C
	  IF(DIF .AND. ND .EQ. NUM_BNDS)THEN
	    T1=DBB/DTDR
	    DO K=1,ND
	      DO I=1,NT-1
	        VJ(I,NUM_BNDS,K)=VJ(I,NUM_BNDS,K)+FA(K)*T1*DIFFW(I)
	      END DO
	      VJ(NT,NUM_BNDS,K)=VJ(NT,NUM_BNDS,K)+ 
	1                       FA(K)*(DDBBDT+T1*DIFFW(NT))
	    END DO
	  ELSE IF(DIF)THEN
	    T1=DBB/DTDR
	    DO J=DIAG_INDX,NUM_BNDS
	      K=ND+DIAG_INDX-J
	      DO I=1,NT-1
	        VJ(I,J,K)=VJ(I,J,K)+FA(K)*T1*DIFFW(I)
	      END DO
	      VJ(NT,J,K)=VJ(NT,J,K)+FA(K)*(DDBBDT+T1*DIFFW(NT))
	    END DO
	  END IF
	END IF
C
