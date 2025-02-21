C
C Routine to compute approximate STARK profiles for any Hydrogenic species
C of charge Z.
C
C Uses the GRIEM theory as modified by AUER and MIHALAS AP J S 24 1972.
C
C Output:
C        PROF( ,I) - Log of STARK profile convolved with a Doppler profile.
C                    The Doppler profile can have a turbulent contribution
C                    given by VDOP. Originally output as Log10, now
C                    output as Loge (8-Aug-91).
C Input:
C        NU        - Velocity in units of C from line center [= (V-Vo)/Vo ].
C        ED_IN     - Electron density (/cm^3) (Vector, length ND)
C        TEMP_IN   - Temperature (10^4 K)     (Vector, length ND)
C        AMASS_IN  - Atomic mass for Doppler profile in amu's.
C        Z_IN      - Ion charge.
C        NL        - Lower level of transition.
C        NUP       - Upper level of transition.
C        ND        - Nuber of (Ne,T) values profile is to be computed for.
C        VDOP_IN   - Turbulent velocity in km/s (function of depth).
C
C Altered  9-Aug-1991 : Call order changed. VDOP_IN is now an array.
C
C Finalized 8-Aug-1991 : Based on original Auer and Mihalas routine
C                        as supplied by Keith Butler (GRIEm.ORIG).
C                        original routine allowed a maximum of only 80
C                        frequencies, and a maximum of 32 levels.
C
	SUBROUTINE SET_GRIEM(PROF,NU,NFREQ,
	1                    ED_IN,TEMP_IN,VDOP_IN,ND,
	1                    NL,NUP,AMASS_IN,Z_IN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER NFREQ,ND,NL,NUP
	REAL(KIND=LDP) PROF(NFREQ,ND),NU(NFREQ)
	REAL(KIND=LDP) ED_IN(ND),TEMP_IN(ND),VDOP_IN(ND)
	REAL(KIND=LDP) AMASS_IN,Z_IN
C
	INTEGER I,J,ML
	REAL(KIND=LDP) LAM_ZERO,NU_ZERO,RHYD
	REAL(KIND=LDP) RYDBERG_HYDROGEN
C
	INTEGER NWS,NF_MAX
	PARAMETER (NF_MAX=100)
	REAL(KIND=LDP) DWS(NF_MAX),PR(NF_MAX)
C
C Check validity of passed parameters.
C
	IF(NUP .LE. NL)THEN
	  WRITE(2,*)'Error : NUP < NL'
	  STOP
	END IF
	IF(NFREQ .GT. NF_MAX)THEN
	  WRITE(2,'(A,I3,A)')' Error in SET_GRIEM - only allowed a '//
	1                    'maximum of ',NF_MAX,' frequencies'
	  STOP
	END IF
C
C Convert from FREQ space to Angstroms. DWS is the distance from
C line center in Angstroms.
C
	RHYD=RYDBERG_HYDROGEN()		!cm^{-1}
	NU_ZERO=RHYD*Z_IN*Z_IN*(NUP-NL)*(NUP+NL)/NUP/NUP/NL/NL	!cm^-1
	LAM_ZERO=1.0E+08_LDP/NU_ZERO				!Angstroms.
	DO I=1,NFREQ
	  DWS(I)=ABS( NU(I)/(1.0_LDP+NU(I)) )*LAM_ZERO
	END DO
C
C Now compute, and store the Doppler broadened Stark profile for each
C depth.
C
	DO I=1,ND
          CALL GRIEM(PR,DWS,NFREQ,ED_IN(I),TEMP_IN(I),
	1             NL,NUP,Z_IN,AMASS_IN,VDOP_IN(I))
	  DO ML=1,NFREQ
	    PROF(ML,I)=PR(ML)
	  END DO
	END DO
C
	RETURN
	END
C
C 
C
      SUBROUTINE GRIEM(PR,DWS,NWS,ED_IN,TEMP_IN,IL,IU,ZZ,AMASS,VDOP)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
C
C Create a profile table using the griem theory as used by A and M
C in their AP J S24 paper.
C
C NWS is the number of wvaelengths profile is to be computed at.
C
      DIMENSION DWS(NWS),PR(NWS)
      COMMON /GRIEM_COM/ SS(109),SX(109),AS,PS,ODOP,NS
C
      DIMENSION BET(55)
      DATA BET/ 0._LDP  ,0.1_LDP  ,0.2_LDP  ,0.3_LDP  ,0.4_LDP  ,0.5_LDP  ,0.6_LDP  ,0.7_LDP  ,0.8_LDP  ,
     ._LDP 0.9_LDP  ,1._LDP  ,1.2_LDP  ,1.4_LDP  ,1.6_LDP  ,1.8_LDP  ,2.0_LDP  ,2.2_LDP  ,2.4_LDP  ,2.6_LDP  ,2.8_LDP  ,
     ._LDP 3.0_LDP  ,3.25_LDP  ,3.5_LDP  ,3.75_LDP  ,4.0_LDP  ,4.5_LDP  ,5.0_LDP  ,6.0_LDP  ,8.0_LDP  ,10._LDP  ,
     ._LDP 12.5_LDP  ,15._LDP  ,17.5_LDP  ,20._LDP  ,25._LDP  ,30._LDP  ,40._LDP  ,50._LDP  ,60._LDP  ,80._LDP  ,
     ._LDP 1.E2_LDP,1.25E2_LDP,1.5E2_LDP,1.75E2_LDP,2.E2_LDP,2.5E2_LDP,3.E2_LDP,4.E2_LDP,5.E2_LDP,6.E2_LDP,8.E2_LDP,
     ._LDP 1.E3_LDP,1.5E3_LDP,2.E3_LDP,2.5E3_LDP /
C
      DATA CLIGHT/2.997925E18_LDP/,  VLIGHT/2.997925E10_LDP/
      DATA CBOLTZ/1.3805E-16_LDP/, AMU/1.67333E-24_LDP/
      DATA THIRD/0.333333333_LDP/
      DATA NBET/55/
      DATA GAMIN,GAMAX/0.01_LDP,25.0_LDP/
      DATA PI/3.1415926536_LDP/
C
C RYD need not be this accurate.
C
      RYD=1.09737312E+05_LDP*VLIGHT*ZZ*ZZ*AMASS/(AMASS+5.48593E-04_LDP)
C
C Set up the constant terms for evaluation of S(BET,GAM) and DONV
      TCON=WAVE*SQRT(2.0_LDP*CBOLTZ/(AMASS*AMU))/VLIGHT
C
      AA=IL
      ASQ=IL*IL
      BSQ=IU*IU
      FRQ=RYD/ASQ-RYD/BSQ               !Line frequency in Hz.
      WAVE=CLIGHT/FRQ			!Wavelength in Angstroms.
      CON=WAVE/FRQ
      OBA=1.0_LDP/(BSQ-ASQ)
      WSQ=WAVE*WAVE
      CCKAB=5.5E-5_LDP*(ASQ*BSQ)**2*OBA/ZZ**5
      BETP1=2.995E-15_LDP*WSQ
      BETW1=6.951E-9_LDP*WSQ*ZZ/BSQ
      GAM01=3.78E-5_LDP*(BSQ-3.0_LDP*AA)*BSQ*OBA/ZZ
      GAM02=4.0E6_LDP*ZZ/BSQ
C
C      TCON=WAVE*SQRT(2.0E0*CBOLTZ/(AMASS*AMU))/VLIGHT
C
C The loops over electron density, and temperature have been removed.
C
C Set Electron density parameters. The electron density is no longer
C logarithmic.
C
      EE=ED_IN
      SRE=SQRT(EE)
      CUE=EE**THIRD
C
C Normal field strength
C
      F0=1.25E-9_LDP*CUE*CUE
C
C Set temperature parameters.
C
      TT=1.0E+04_LDP*TEMP_IN  		!As TEMP_IN in units of 10^4 K.
      SRT=SQRT(TT)
C
C Doppler width (A) (old expression was TCON*SRT). 10^10 arrises as we need
C to convert VDOP from km/s to cm/s.
C
      DOP=WAVE*SQRT( 2.0_LDP*CBOLTZ*TT/(AMASS*AMU) +
	1               1.0E+10_LDP*VDOP*VDOP )/VLIGHT
      ODOP=1.0_LDP/DOP
C
C Ratio STARK/DOPPLER width.
C
      CKABD=CCKAB*F0*ODOP
C
C 1/STARK WIDTH     (A)**-U
C
      OCKAB=ODOP/CKABD
      BETAP=BETP1*SRE*OCKAB
      BETAW=BETW1*TT*OCKAB
      BETAP=MIN(0.99_LDP*BETAW,BETAP)
      GAM0=GAM01*CUE*LOG10(GAM02*TT/SRE)/SRT
      GAML=LOG(BETAW/BETAP)
      GAMW=4.712389_LDP/SQRT(BETAW)
      FPW=GAM0/GAML
C
C Initialize GAMA for ZERO,BETAP
C
      GAMA=GAMW+GAM0
C
C Generate S(BETA) on a fixed basis.
C
      JM=NBET
      JP=JM
      DO 6 I=1,NBET
        BETA=BET(I)
        IF(BETA.LT.BETAP) GOTO 11
        IF(BETA.GT.BETAW) GOTO 12
C
C ASSIGN GAMA FOR BETAP,BETAW
C
        GAMA=GAMW+FPW*LOG(BETAW/BETA)
C
C Look for regions of the GAMA,BETA plane where fast methods are ok.
C This is not just to save time as the T(BETA,GAMA) analysis
C Is numerically unstable in these same regions.
C
   11   IF(GAMA.GT.GAMAX) GOTO 13
        IF(GAMA.LT.GAMIN*BETA) GOTO 14
        IF(GAMA.LT.GAMIN) GOTO 14
C FULL CASE
        SSS=TBG(BETA,GAMA)
        GOTO 7
C   BETA GT 100*GAMA FOR BETA LT 1
   14   SSS=TH(BETA)
        IF(BETA.LT.1.0_LDP) GOTO 7
C GAMA LT GAMIN OR BETA GT 100*GAMA FOR BETA GT 1
        SSS=SSS+GAMA/(PI*BETA*BETA)
        GOTO 7
C   BETA GT BETAW
   12   SSS=2.0_LDP*TH(BETA)
        GOTO 7
C   GAMA GT GAMAX
   13   SSS=GAMA/(PI*(GAMA*GAMA+BETA*BETA))
C FILL THE SYMETRIC SS,SX VECTORS
    7   SS(JM)=SSS
        SX(JM)=-BETA*CKABD
        SS(JP)=SSS
        SX(JP)= BETA*CKABD
        JM=JM-1
        JP=JP+1
    6   CONTINUE
      NS=2*NBET-1
C Normalize to unit wavelength integral.
      REN=0.0_LDP
      DO 8 I=2,NS
        REN=REN+(SS(I)+SS(I-1))*(SX(I)-SX(I-1))
    8 CONTINUE
      FAC=CON*ODOP*2.0_LDP/REN
C Make the asymtotic power law constants.
      PS=LOG(SS(NS)/SS(NS-1)) / LOG(SX(NS)/SX(NS-1))
      PS=MIN(PS,-2.0_LDP)
      AS=SS(NS)/SX(NS)**PS
C
C Map the profile onto the DLAM set by convolution with the doppler
C profile.
C
      DO 10 I=1,NWS
        PR(I)=LOG(FAC*DCONV(DWS(I)))
   10 CONTINUE
      RETURN
      END
C
C 
C
      FUNCTION DCONV(DLAM)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
C
C Convolution of gaussian profile with S function  6 AUG 77
C
      COMMON /GRIEM_COM/ SS(109),SX(109),AS,PS,ODOP,NS
      DIMENSION X(109),ERX(109),EX(109)
      DATA HALF/0.5_LDP/, ZERO/0.0_LDP/, SRTPI/5.6418958E-1_LDP/
      DATA RANGE/6.0_LDP/
      DB=ABS(DLAM)*ODOP
      X1=SX(1)+DB
      IF(X1.LT.RANGE) GOTO 10
C ASINT -RANGE,RANGE
      DCONV=ASINT(-RANGE,RANGE,DB)
      RETURN
   10 DCONV=ZERO
C ASINT -RANGE,X1
      IF(X1.GT.-RANGE) DCONV=DCONV+ASINT(-RANGE,X1,DB)
C ASINT XN,RANGE
      XN=SX(NS)+DB
      IF(XN.LT.RANGE) DCONV=DCONV+ASINT(XN,RANGE,DB)
C SET UP X IN THE GAUSSIAN FRAME
C STORE ALL EXPONENTIAL AND ERROR FUNCTIONS ONCE
      N0=2
      CON=HALF*SRTPI
      DO 20 I=1,NS
         X(I)=SX(I)+DB
         IF(X(I).LT.-RANGE) N0=I+1
         EX(I)=DOPLER(X(I))
         ERX(I)=ERR(X(I))
   20 CONTINUE
C BY S SEGMENT
      TERX=ZERO
      TEX=ZERO
      DO 50 J=N0,NS
         JM=J-1
         DX=X(J)-X(JM)
         IF(DX.LT.0.1_LDP) GOTO 30
         GR=(SS(J)-SS(JM))/DX
         CR=SS(J)-GR*X(J)
         TERX=TERX+CR*(ERX(J)-ERX(JM))
         TEX=TEX+GR*(EX(J)-EX(JM))
         GOTO 40
   30    DCONV=DCONV +  CON*(SS(J)*EX(J)+SS(JM)*EX(JM))*DX
   40    CONTINUE
         IF(X(J).GT.RANGE) GOTO 60
   50 CONTINUE
   60 DCONV=DCONV+HALF*(TERX-SRTPI*TEX)
      RETURN
      END
      FUNCTION TBG(BET,GAM)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
      DATA PI/3.14159265358979_LDP  /
      IF(GAM.GT.1.E-2_LDP) GOTO 10
C     SMALL GAMMA
      TBG=TH(BET)
      IF(BET.GT.1._LDP) TBG=TBG+GAM/PI/BET**2
      RETURN
C     NORMAL CASE
   10 IF(BET/GAM.GT.1.E2_LDP) GOTO 20
      T=GAM*(TF(BET,GAM)+TF(-BET,GAM)+TG(BET,GAM)+TG(-BET,GAM))/PI
      TBG=T
      RETURN
   20 TBG=TH(BET)+GAM/PI/BET**2
      RETURN
      END
C
C 
C
      FUNCTION TG(BET,GAM)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
      DATA C4,C5,C6/1.5_LDP  ,3.4636008E1_LDP,-1.3253986E2_LDP/
      DATA PI2/1.57079632679489_LDP  /
      G2=GAM*GAM
      B=2.0_LDP  *BET
      C=BET*BET+G2
      BET2=BET*BET
      CC=C*C
      R= SQRT(C)
      Q= SQRT(R)
      Q2=Q*Q
      P=1._LDP  /Q
      P2=P*P
      SINA=GAM/R
      COSA2= SQRT(0.5_LDP  *(1._LDP  -BET/R))
      SINA2= SQRT(0.5_LDP  *(1._LDP  +BET/R))
      X5=BET+4.
      Y2=LOG((X5*X5+G2)/16._LDP  )
      IF(X5/GAM.GT.1.0_LDP  ) GOTO 10
      Y1=PI2- ATAN(X5/GAM)
      GOTO 20
   10 Y1= ATAN(GAM/X5)
   20 D2=C/192.
      D3=-B/32.
      D4=0.25_LDP  *(3._LDP  *BET2-G2)/C
      D5=B*(G2-BET2)/CC
      D6=(BET2*(BET2-6.0_LDP  *G2)+G2*G2)/CC
      SUM1=((D6*Y1)/GAM+(D4+D3))+D2
      SUM1=(SUM1+D5*Y2)*C5/CC
      D1=C/1024.
      D2=-B/192.
      D3=(3.0_LDP  *BET2-G2)/(32._LDP  *C)
      D4=BET*(G2-BET2)/CC
      D5=0.5_LDP  *(G2*G2+5._LDP  *BET2*(BET2-2.0_LDP  *G2))/(C*CC)
      D6=-BET*(BET2*BET2+5.0_LDP  *G2*(G2-2.0_LDP  *BET2))/(CC*C)
      SUM2=(D5*Y2+D3)+D1
      SUM2=(SUM2+(D2+D4+(D6*Y1)/GAM))*C6/CC
      D7=C4/C
      D8=D7*(B*B/C-1._LDP  )/(2._LDP  *Q*Q2*SINA)
      D9=D7*(B/(C*C))/(2._LDP  *P*P2*SINA)
      X1=(4.0_LDP  -Q2)/(4.0_LDP  *Q*SINA2)
      X2=(Q*(Q+4._LDP  *COSA2)+4._LDP  )/(Q*(Q-4._LDP  *COSA2)+4._LDP  )
      X3=(0.25_LDP  -P2)/(P*SINA2)
      X4=(P*(P+COSA2)+0.25_LDP  )/(P*(P-COSA2)+0.25_LDP  )
      IF(X1.GT.1._LDP  ) GOTO 30
      Y1=PI2- ATAN(X1)
      GOTO 40
   30 Y1= ATAN(1._LDP  /X1)
   40 IF(X3.GT.-1._LDP  ) GOTO 50
      Y2=- ATAN(1._LDP  /X3)
      GOTO 60
   50 Y2=PI2+ ATAN(X3)
   60 SUM3=D8*(2._LDP  *COSA2*Y1-SINA2*LOG(X2))
      SUM4=D9*(2._LDP  *COSA2*Y2+SINA2*LOG(X4))
      TG=(SUM4+D7*(1._LDP  /12._LDP  -B/C))+SUM3
      TG=TG+SUM1+SUM2
      RETURN
      END
C
C 
C
      FUNCTION TF(B,G)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
      DATA C0,C1,C2,C3/1.0007744E-1_LDP,4.93208719E-3_LDP,-7.09873526E-3_LDP,
     ._LDP 7.11559325E-4_LDP/
      D3=C3
      D2=C2-3.0_LDP  *B*C3
      D1=(3.0_LDP  *C3*B-2.0_LDP  *C2)*B+C1
      D =((-C3*B+C2)*B-C1)*B+C0
      G2=G*G
      X1=B+4.
      TF=4.0_LDP  *D2+4.0_LDP  *D3*(B+2.0_LDP  )+0.5_LDP*(D1-G2*D3)*LOG((X1*X1+G2)/(B*B
     ._LDP +G2)) +(D -G2*D2)*( ATAN(X1/G)- ATAN(B/G))/G
      RETURN
      END
C
C 
C
      FUNCTION TH(X)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
C GRIEM MICROFIELD FUDGE - NORMALIZED TO UNIT B INTEGRAL
      DATA C0,C1,C2,C3/1.0007744E-1_LDP,4.93208719E-3_LDP,-7.09873526E-3_LDP,
     ._LDP 7.11559325E-4_LDP/
      DATA C4,C5,C6/1.5_LDP  ,3.4636008E1_LDP,-1.3253986E2_LDP/
      B=ABS(X)
      IF(B.GT.4._LDP  ) GOTO 10
         TH=((C3*B+C2)*B+C1)*B+C0
      RETURN
   10    SRT= SQRT(B)
         B2=B*B
         TH=((C6/B+C5)/B2+C4/SRT)/B2
      RETURN
      END
C
C 
C
      FUNCTION ASINT(X0,X1,DB)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
C
C Integral over asymtotic region by SIMPSONS rule.
C
      COMMON /GRIEM_COM/ SS(109),SX(109),AS,PS,ODOP,NS
      DATA SRTPI/5.6418958E-1_LDP/,  C23/6.666666667E-1_LDP/
      DATA HALF,ONE,TWO/0.5_LDP,1.0_LDP,2.0_LDP/
      DATA STEP/0.2_LDP/
      F(X)=DOPLER(X)*ABS(X-DB)**PS
C
C Chose H to be sept
C
      N=(X1-X0)/STEP+ONE
      DX=(X1-X0)/FLOAT(N)
      H=HALF*DX
      X=X0-H
      ASINT=F(X0)
      DO 10 I=1,N
        X=X+DX
        ASINT=ASINT+  TWO*F(X)+F(X+H)
   10 CONTINUE
      ASINT=ASINT*AS*SRTPI*H*C23
      RETURN
      END
C
C 
C
      FUNCTION DOPLER(X)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
      SQ=X*X
      IF(SQ.GT.64._LDP) THEN
        DOPLER=0.
      ELSE
        DOPLER=EXP(-SQ)
      END IF
      RETURN
      END
C
C 
C
C Error function.
C
      FUNCTION ERR(X)
	USE SET_KIND_MODULE
      IMPLICIT REAL(KIND=LDP) (A-H,O-Z)
      ONE=1.0_LDP
      T=ABS(X)
      IF(T.GT.5._LDP)THEN
        ERR=SIGN(ONE,X)
      ELSE
        T=1._LDP/(1._LDP+0.3275911_LDP*T)
        ERR=((((1.061405429_LDP  *T-1.453152027_LDP  )*T +1.421413741_LDP  )*T-
     ._LDP 0.284496736_LDP  )*T +0.254829592_LDP  )*T
        EX=DOPLER(X)
        ERR=SIGN(ONE-EX*ERR,X)
      END IF
      RETURN
      END
