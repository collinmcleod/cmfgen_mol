	SUBROUTINE COMP_JBAR(ETA_CONT,CHI_CONT,ESEC,ETAL,CHIL,
	1                T,V,SIGMA,R,P,AQW,HMIDQW,KQW,NMIDQW,
	1                JBAR,DIF,DTDR,IC,THK_CONT,
	1                VTURB,VDOP_FG_FRAC,VDOP_MOM_FRAC,
	1                RED_EXT,AMASS,LINE_FREQ,
	1                METHOD,FG_SOL_OPTIONS,N_TYPE,
	1                LUM,PLT_J,PLT_H,PLT_LF,PLT_FM,
	1                LEV,NLEV,NLF,NC,NP,ND)
	USE SET_KIND_MODULE
!
	USE MOD_WR_STRING
	USE MOD_USR_OPTION
!
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NC
	INTEGER NP
	INTEGER NLF
	INTEGER NLEV
!	
	REAL(KIND=LDP) ETA_CONT(ND)
	REAL(KIND=LDP) CHI_CONT(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) ETAL(ND)
	REAL(KIND=LDP) CHIL(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) P(NP)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) JBAR(ND)
!
	REAL(KIND=LDP) AQW(ND,NP)
	REAL(KIND=LDP) HMIDQW(ND,NP)
	REAL(KIND=LDP) KQW(ND,NP)
	REAL(KIND=LDP) NMIDQW(ND,NP)
!
	REAL(KIND=LDP) DTDR
	REAL(KIND=LDP) IC
	REAL(KIND=LDP) LUM
	REAL(KIND=LDP) RED_EXT
	REAL(KIND=LDP) AMASS
	REAL(KIND=LDP) VTURB
	REAL(KIND=LDP) VDOP_FG_FRAC
	REAL(KIND=LDP) VDOP_MOM_FRAC
	REAL(KIND=LDP) LINE_FREQ
!
	LOGICAL DIF
	LOGICAL THK_CONT
	LOGICAL PLT_J
	LOGICAL PLT_H
	LOGICAL PLT_LF
	LOGICAL PLT_FM
!
	CHARACTER*(*) METHOD
	CHARACTER*(*) FG_SOL_OPTIONS
	CHARACTER*(*) N_TYPE
!
	INTEGER LEV(NLEV)
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	EXTERNAL SPEED_OF_LIGHT
	REAL(KIND=LDP) SPEED_OF_LIGHT
!
! Local arrays and variables
!
	REAL(KIND=LDP) RJ(ND)
	REAL(KIND=LDP) RJ_ES(ND)
	REAL(KIND=LDP) RJ_STORE(ND,NLF)
	REAL(KIND=LDP) RJ_ES_STORE(ND,NLF)
!
	REAL(KIND=LDP) LAM(NLF)
	REAL(KIND=LDP) RSQ_LF_STORE(NLF,ND)
	REAL(KIND=LDP) RSQH_STORE(NLF,ND)
	REAL(KIND=LDP) RSQJ_STORE(NLF,ND)
!
	REAL(KIND=LDP) ETA(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) RH(ND)
	REAL(KIND=LDP) RSQHNU(ND)
	REAL(KIND=LDP) FEDD(ND)
	REAL(KIND=LDP) GEDD(ND)
	REAL(KIND=LDP) N_ON_J(ND)
	REAL(KIND=LDP) SOB(ND)
	REAL(KIND=LDP) FORCE_MULT(ND)
	REAL(KIND=LDP) FORCE_MULT_FG(ND)
!
	REAL(KIND=LDP) CHI_PREV(ND)
	REAL(KIND=LDP) ETA_PREV(ND)
	REAL(KIND=LDP) FEDD_PREV(ND)
	REAL(KIND=LDP) GEDD_PREV(ND)
	REAL(KIND=LDP) JNU_PREV(ND)
	REAL(KIND=LDP) RSQHNU_PREV(ND)
	REAL(KIND=LDP) N_ON_J_PREV(ND)
!
	REAL(KIND=LDP) TA(ND)
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) TC(ND)
	REAL(KIND=LDP) FQW(ND)
	REAL(KIND=LDP) LINE_QW(ND)
	REAL(KIND=LDP) NORM_CHK(ND)
	REAL(KIND=LDP) FOLD(ND)
	REAL(KIND=LDP) DOP_NU_VEC(ND)
	REAL(KIND=LDP) VDOP_VEC(ND)
!
	REAL(KIND=LDP) IPLUS(NP)
!
	REAL(KIND=LDP) HBC_CMF(3)
	REAL(KIND=LDP) INBC(3)
	REAL(KIND=LDP) NBC_CMF(3)
!
	REAL(KIND=LDP) HBC_PREV(3)
	REAL(KIND=LDP) INBC_PREV(3)
	REAL(KIND=LDP) NBC_PREV(3)
!
	REAL(KIND=LDP) FL
	REAL(KIND=LDP) FL_OLD
	REAL(KIND=LDP) VDOP_MAX
	REAL(KIND=LDP) VDOP_MIN
	REAL(KIND=LDP) VDOP
	REAL(KIND=LDP) DOP_NU
	REAL(KIND=LDP) DEL_NU
	REAL(KIND=LDP) NU_ST
	REAL(KIND=LDP) EMHNUKT
	REAL(KIND=LDP) DBB
	REAL(KIND=LDP) dLOG_NU
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) T2
	REAL(KIND=LDP) H_IN
	REAL(KIND=LDP) H_OUT
!
	REAL(KIND=LDP) AMASS_LOC
	REAL(KIND=LDP) VTURB_MIN
	REAL(KIND=LDP) VTURB_IM
	REAL(KIND=LDP) VTURB_MAX
!
	INTEGER I
	INTEGER L
	INTEGER ML
	INTEGER FG_COUNT
!
	LOGICAL CONT_VEL
	LOGICAL FIRST_FREQ
	LOGICAL NEW_FREQ
	LOGICAL INACCURATE
	LOGICAL COHERENT_ES
!
	CHARACTER*10 TEMP_STR
	CHARACTER*132 DEFAULT
!
	CONT_VEL=.TRUE.
	COHERENT_ES=.TRUE.
!
	JBAR(:)=0.0_LDP
	RJ_STORE(:,:)=0.0_LDP
	RJ_ES_STORE(:,:)=0.0_LDP
	FOLD(:)=0.0_LDP
	NORM_CHK(:)=0.0_LDP
	FORCE_MULT(:)=0.0_LDP
	FORCE_MULT_FG(:)=0.0_LDP
!
!	DO I=1,20
!	  ESEC(I)=ESEC(11)
!	  CHI_CONT(I)=CHI_CONT(11)
!	END DO
!	CHI_CONT(1)=CHI_CONT(20)
!	ESEC(1)=ESEC(20)
!
!	CHI_CONT(1)=CHI_CONT(1)+999.0D0*ESEC(1)
!	ESEC(1)=1000.0D0*ESEC(1)
!
!	DO I=1,ND
!	  IF(ESEC(I)*(R(I)-R(ND)) .LT. 0.001D0)THEN
!	     CHI_CONT(I)=CHI_CONT(I)-ESEC(I)
!	     T1=ESEC(I)
!	     ESEC(I)=0.001D0/(R(I)-R(ND))
!	     CHI_CONT(I)=CHI_CONT(I)+ESEC(I)
!	     WRITE(6,*)I,T1,ESEC(I)
!	   ELSE
!	     EXIT
!	   END IF
!	END DO
!
	DEFAULT=WR_STRING(AMASS)
	CALL USR_OPTION(AMASS_LOC,'AMASS',DEFAULT,'Atomic mass in AMU')
	DEFAULT=WR_STRING(VTURB)
	CALL USR_OPTION(VTURB_MIN,'VT_MIN',DEFAULT,'Minimum turbulent velocity in km/s')
	DEFAULT=WR_STRING(VTURB_MIN)
	CALL USR_OPTION(VTURB_IM,'VT_IM',DEFAULT,'Intermediate turbulent velocity in km/s')
	DEFAULT=WR_STRING(VTURB_IM)
	CALL USR_OPTION(VTURB_MAX,'VT_MAX',DEFAULT,'Maximum turbulent velocity in km/s')

	VDOP_MAX=0.0_LDP
	VDOP_MIN=1000.0_LDP
	DO I=1,ND
	  T1=VTURB_MIN+MIN((VTURB_IM-VTURB_MIN),V(I))+ (VTURB_MAX-VTURB_IM)*V(I)/V(1)
	  T2=12.85_LDP*SQRT( T(I)/AMASS_LOC + (T1/12.85_LDP)**2 )
	  VDOP_MAX=MAX(VDOP_MAX,T2)
	  VDOP_MIN=MAX(VDOP_MIN,T2)
	  VDOP_VEC(I)=T2
	  DOP_NU_VEC(I)=LINE_FREQ*T2/2.998E+05_LDP
	END DO
	VDOP=VDOP_MAX
!
	DOP_NU=LINE_FREQ*VDOP/2.998E+05_LDP
	DEL_NU=(6.0_LDP+RED_EXT)*DOP_NU/(NLF-1)
	FIRST_FREQ=.TRUE.
	NU_ST=LINE_FREQ+6.0_LDP*DOP_NU
!
	WRITE(6,*)'ND=',ND
	WRITE(6,*)'NLF=',NLF
	WRITE(6,*)'VTURB=',VTURB
	WRITE(6,*)'AMASS=',AMASS_LOC
	WRITE(6,*)'LINE_FREQ=',LINE_FREQ
	WRITE(6,*)'MEANT=',T1
	WRITE(6,*)'DOP_NU_MIN=',MINVAL(DOP_NU_VEC)
	WRITE(6,*)'DOP_NU_MAX=',MAXVAL(DOP_NU_VEC)
	WRITE(6,*)'DEL_NU=',DEL_NU
	WRITE(6,*)'DEL_NU/DOP_NU=',DEL_NU/DOP_NU
!
	DO ML=1,NLF
!
	  FL= NU_ST-(ML-1)*DEL_NU
!
! Compute DBB for diffusion approximation. DBB=dB/dR
! and DDBBDT= dB/dTR .
!
	  T1=HDKT*FL/T(ND)
	  EMHNUKT=EXP(-T1)
	  T2=1.0_LDP-EMHNUKT
	  DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT/(T2**2)
!
	  DO I=1,ND
	    T1=1.0E-15_LDP/1.77245385095516_LDP/DOP_NU_VEC(I)
            T2=T1*EXP(- ((FL-LINE_FREQ)/DOP_NU_VEC(I))**2 )
	    CHI(I)=CHI_CONT(I)+T2*CHIL(I)
	    ETA(I)=ETA_CONT(I)+T2*ETAL(I)
	    FQW(I)=1.0E+15_LDP*DEL_NU
	    LINE_QW(I)=1.0E+15_LDP*DEL_NU*T2
	    NORM_CHK(I)=NORM_CHK(I)+LINE_QW(I)
	  END DO
C
C NB: CHI_PREV is used to refer to the continuum opacity at the previous
C frequency. Is does not need to be multiplied by CLUMP_FAC, as it is compared
C directly to CHI_CONT.
C
C For HBC and NBC only the first vector element is used.
C
	  NEW_FREQ=.TRUE.
	  IF(FIRST_FREQ)THEN
	    dLOG_NU=0.0_LDP
	    DO I=1,ND
	      RJ(I)=0.0_LDP
	      CHI_PREV(I)=CHI(I)
	      ETA_PREV(I)=ETA(I)
	      FEDD_PREV(I)=0.0_LDP		!Not required.
	      GEDD_PREV(I)=0.0_LDP
	      JNU_PREV(I)=0.0_LDP
	      N_ON_J_PREV(I)=0.0_LDP
	      RSQHNU_PREV(I)=0.0_LDP
	    END DO
	    HBC_PREV(:)=0.0_LDP		!1:3
	    NBC_PREV(:)=0.0_LDP		!1:3
	    HBC_CMF(:)=0.0_LDP		!1:3
	    NBC_CMF(:)=0.0_LDP		!1:3
	    FG_COUNT=0
	  ELSE
	    dLOG_NU=LOG(FL_OLD/FL)
	    DO I=1,ND
	      FEDD_PREV(I)=FEDD(I)
	      GEDD_PREV(I)=GEDD(I)
	      N_ON_J_PREV(I)=N_ON_J(I)
	      JNU_PREV(I)=RJ(I)
	      RSQHNU_PREV(I)=RSQHNU(I)
	    END DO
	    HBC_PREV(:)=HBC_CMF(:)
	    NBC_PREV(:)=NBC_CMF(:)
	  END IF
!
!	 RJ(:,ML)=RJ_STORE(:,ML)
!	 IF(.NOT. COHERENT_ES)THEN
!	    RJ_ES(:,ML)=RJ_ES_STORE(:,ML)
!	 END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  L=0
	  DO WHILE(INACCURATE)
C
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA(1:ND)+ESEC(1:ND)*RJ(1:ND)
	     ELSE
	       TA(1:ND)=ETA(1:ND)+ESEC(1:ND)*RJ_ES(1:ND)
	     END IF
C
C NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
C
	     IF(VDOP_FG_FRAC .NE. 0)THEN
	       CALL FG_J_CMF_V10(TA,CHI,ESEC,V,SIGMA,R,P,
	1                  TC,TB,FEDD,GEDD,N_ON_J,
	1                  AQW,HMIDQW,KQW,NMIDQW,
	1                  INBC,HBC_CMF,NBC_CMF,
	1                  IPLUS,FL,dLOG_NU,DIF,DBB,IC,
	1                  VDOP_VEC,VDOP_FG_FRAC,
	1                  METHOD,FG_SOL_OPTIONS,THK_CONT,
	1                  FIRST_FREQ,NEW_FREQ,N_TYPE,NC,NP,ND)
	     ELSE
	       CALL FG_J_CMF_V9(TA,CHI,ESEC,V,SIGMA,R,P,
	1                  TC,TB,FEDD,GEDD,N_ON_J,
	1                  AQW,HMIDQW,KQW,NMIDQW,
	1                  INBC,HBC_CMF,NBC_CMF,
	1                  IPLUS,FL,dLOG_NU,DIF,DBB,IC,
	1                  METHOD,FG_SOL_OPTIONS,THK_CONT,
	1                  FIRST_FREQ,NEW_FREQ,N_TYPE,NC,NP,ND)
	     END IF
	     FG_COUNT=FG_COUNT+1
C
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA(1:ND)
	     ELSE
	       TA(1:ND)=ETA(1:ND)+ESEC(1:ND)*RJ_ES(1:ND)
	     END IF
	     IF(VDOP_MOM_FRAC .NE. 0)THEN
	       CALL MOM_J_CMF_V6(TA,CHI,ESEC,V,SIGMA,R,
	1                  FEDD,GEDD,N_ON_J,RJ,RSQHNU,
	1                  VDOP_VEC,VDOP_MOM_FRAC,
	1                  HBC_CMF,INBC,NBC_CMF,
	1                  FL,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,COHERENT_ES,
	1                  FIRST_FREQ,NEW_FREQ,NC,NP,ND)
	     ELSE
	       CALL MOM_J_CMF_V5(TA,CHI,ESEC,V,SIGMA,R,
	1     	           FEDD,GEDD,N_ON_J,
	1                  FEDD_PREV,GEDD_PREV,N_ON_J_PREV,
	1                  RJ,RSQHNU,JNU_PREV,RSQHNU_PREV,
	1                  HBC_CMF,INBC,NBC_CMF,
	1                  HBC_PREV,INBC_PREV,NBC_PREV,
	1                  FL,dLOG_NU,DIF,DBB,IC,METHOD,COHERENT_ES,
	1                  FIRST_FREQ,NEW_FREQ,NC,NP,ND)
	     END IF
	     WRITE(163,'(7ES14.6)')FL,TC(1)*R(1)*R(1),RJ(1)*R(1)*R(1),TB(1)*0.25*(R(1)+R(2))**2,
	1                           RSQHNU(1),HBC_CMF(1),NBC_CMF(1)
C
C We set NEW_FREQ to false so that FG_J_CMF continues to use the same
C AV_PREV and CV_PREV. NEW_FREQ must be set to true again outside the
C F iteration loop.
C
	     NEW_FREQ=.FALSE.
C
C Update "inaccurate" iteration counter
C
	      L=L+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(L .LT. 20)THEN
	        T1=0.0_LDP
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. 1.0E-07_LDP)INACCURATE=.TRUE.
	      END IF
C
	      IF(L .GT. 10)THEN
	         WRITE(6,*)'Possible error converging f - T1 is',T1
	      	 INACCURATE=.FALSE.
	      END IF
	  END DO
!
	  H_IN=DBB/CHI(ND)/3.0_LDP
	  H_OUT=HBC_CMF(1)*RJ(1)
	  CALL REGRID_H(SOB,R,RSQHNU,H_OUT,H_IN,ND,TA)
	  DO I=1,ND
	    FORCE_MULT(I)=FORCE_MULT(I)+SOB(I)*FQW(I)*(CHI(I)-CHI_CONT(I))
	  END DO
	  WRITE(72,'(I3,3X,6ES12.4)')ML,SOB(11),SOB(15),(CHI(10)-CHI_CONT(10)),(CHI(15)-CHI_CONT(15)),FORCE_MULT(10),FORCE_MULT(1)
	  DO I=1,ND-1
	    TB(I)=TB(I)*0.25_LDP*(R(I)+R(I+1))**2
	  END DO
	  CALL REGRID_H(SOB,R,TB,H_OUT,H_IN,ND,TA)
	  DO I=1,ND
	    FORCE_MULT_FG(I)=FORCE_MULT_FG(I)+SOB(I)*FQW(I)*(CHI(I)-CHI_CONT(I))
	  END DO
C
	  LAM(ML)=1.0E-07_LDP*SPEED_OF_LIGHT()/FL
	  RSQ_LF_STORE(ML,1:ND)=RSQHNU(1:ND)*(CHI(1:ND)-CHI_CONT(1:ND))
	  RSQH_STORE(ML,1:ND)=RSQHNU(1:ND)
	  RSQJ_STORE(ML,1:ND)=RJ(1:ND)*R(1:ND)*R(1:ND)
!
	  FIRST_FREQ=.FALSE.
	  FL_OLD=FL
	  JBAR(:)=JBAR(:)+LINE_QW(:)*RJ(:)
	  TEMP_STR='RSQH'
	  WRITE(TEMP_STR(5:7),'(I3.3)')ML
	  CALL WRITV(RSQHNU,ND-1,TEMP_STR,27)
	  TEMP_STR='RJ'
	  WRITE(TEMP_STR(3:5),'(I3.3)')ML
	  CALL WRITV(RJ,ND,TEMP_STR,28)
	END DO
	TEMP_STR='JBAR'
	CALL WRITV(JBAR,ND,TEMP_STR,28)
	TEMP_STR='NORM_CHK'
	CALL WRITV(NORM_CHK,ND,TEMP_STR,28)
!
	IF(PLT_J)THEN
	  I=1
	  DO WHILE(LEV(I) .GE. 1 .AND. LEV(I) .LE. ND)
	    CALL DP_CURVE(NLF,LAM,RSQJ_STORE(:,LEV(I)))
	    I=I+1
	    IF(I .GT. NLEV)EXIT
	  END DO
	END IF
!
	IF(PLT_H)THEN
	  I=1
	  DO WHILE(LEV(I) .GE. 1 .AND. LEV(I) .LE. ND)
	    CALL DP_CURVE(NLF,LAM,RSQH_STORE(:,LEV(I)))
	    I=I+1
	    IF(I .GT. NLEV)EXIT
	  END DO
	END IF
!
	IF(PLT_LF)THEN
	  I=1
	  DO WHILE(LEV(I) .GE. 1 .AND. LEV(I) .LE. ND)
	    CALL DP_CURVE(NLF,LAM,RSQ_LF_STORE(:,LEV(I)))
	    I=I+1
	    IF(I .GT. NLEV)EXIT
	  END DO
	END IF
!
	IF(PLT_FM)THEN
	  DO I=1,ND
	    TA(I)=LOG10(R(I)/R(ND))
	    FORCE_MULT(I)=FORCE_MULT(I)*4.1274E-12_LDP/LUM/ESEC(I)
	    FORCE_MULT_FG(I)=FORCE_MULT_FG(I)*4.1274E-12_LDP/LUM/ESEC(I)
	  END DO
	  CALL DP_CURVE(ND,V,FORCE_MULT)
	  CALL DP_CURVE(ND,V,FORCE_MULT_FG)
!	  CALL DP_CURVE(ND,TA,FORCE_MULT)
	END IF
!
	RETURN
	END
