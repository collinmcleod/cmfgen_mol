C
C This routine solves for the mean intensity as a function of depth using the
C Feautrier Technique. A Schuster or diffusion approaximation is used for the
C lower boundary condition. This routine may need to be in a loop so that
C the f values are iterated to convergence.
C
C Created 24-FEB-1987
C Altered 21-Sep-1990 --- Bug fix: IBOUND was not being zeroed when THK
C                         option was off. Only effected direct computation of
C                         HBCNEW.
C Altered 04-Oct-1990 --- INBCNEW was being evaluated incorrectly.
C                         Not used if DIFF is true.
C Altered 15-Oct-2023 - Fixed bug -- DACOS was relaced by COS instead of ACOS.
C
	SUBROUTINE FQCOMP(TA,TB,TC,XM,DTAU,R,Z,P,NEWRJ,NEWRK
	1  ,SOURCE,CHI,dCHIdr,AQW,AQW3,DBB,HBCNEW
	1  ,INBCNEW,IC,S1,THK,DIFF,NC,ND,NP,METHOD)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
	INTEGER NC,ND,NP,I,NI,LS
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND),XM(ND),R(ND),Z(ND)
	REAL(KIND=LDP) NEWRJ(ND),NEWRK(ND),DTAU(ND),dCHIdR(ND)
	REAL(KIND=LDP) SOURCE(ND),CHI(ND),AQW(ND,NP),AQW3(ND,NP),P(NP)
C
	REAL(KIND=LDP) S1,DBB,DBC,IBOUND,TOR,HBCNEW,INBCNEW,IC,E1,E2,E3
	LOGICAL DIFF,THK
	CHARACTER*6 METHOD
C
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
C
C Zero all parameters.
C
	CALL DP_ZERO(NEWRJ,ND)
	CALL DP_ZERO(NEWRK,ND)
	HBCNEW=0.0_LDP
	INBCNEW=0.0_LDP
C
C ENTER LOOP FOR EACH IMPACT PARAMETER P
C
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*SQRT( (R(ND)-P(LS))*(R(ND)+P(LS)) )
	1   /R(ND)/CHI(ND)
	  END IF
C
	  IF(THK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796_LDP-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=S1*(1.0_LDP-EXP(-TOR))
	  ELSE
	    IBOUND=0.0_LDP
	  END IF
C
C Compute Z and the optical depth scale DTAU for this imapct parameter.
C
	  IF(NI .GT. 1)THEN
C	    CALL ZALONGP(R,Z,P(LS),NI)
C	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdr,NI)
	    DO I=1,NI
	      Z(I)=SQRT( (R(I)-P(LS))*(R(I)+P(LS)) )
	    END DO
	    DO I=1,NI-1
	      DTAU(I)=0.5_LDP*(Z(I)-Z(I+1))*(CHI(I)+CHI(I+1)+(Z(I)-Z(I+1))
	1     *(dCHIdR(I+1)*Z(I+1)/R(I+1)-dCHIdR(I)*Z(I)/R(I))/6.0_LDP)
	    END DO
	  END IF
C
C Compute XM. Compute T ( a tridiagonal matrix) and store it as three vectors
C TA,TB and TC . This code is a combined version of XVECD and TCOMPD.
C
C	    CALL XVECD(DTAU,SOURCE,XM,DIFF,DBC,IC,LS,NC,ND,NI)
C	    IF(THK)XM(1)=-IBOUND
C	    CALL TCOMPD(TA,TB,TC,DTAU,DIFF,LS,NC,ND,NI)
C
	  IF(NI .GT. 2)THEN
	    XM(1)=-IBOUND
	    TA(1)=0.0
	    TC(1)=1._LDP/DTAU(1)
	    TB(1)=-1.0_LDP-TC(1)
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TC(I)=1.0_LDP/DTAU(I)
	      TB(I)=-0.5_LDP*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5
	    END DO
C
	    IF(LS .LE. NC .AND. DIFF)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=TC(NI-1)
	      XM(NI)=DBC
	    ELSE IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-TA(NI)+DTAU(NI-1)/2.
	      XM(NI)=0.5_LDP*DTAU(NI-1)*SOURCE(NI)
	    ELSE
	      TA(NI)=-TC(NI-1)
	      TB(NI)=1-TA(NI)
	      XM(NI)=IC
	    END IF
	    TC(NI)=0.0
C
C Solve the tridiagonal system of equations.
C
	    CALL THOMAS(TA,TB,TC,XM,NI,1)
C
	  ELSE IF(NI .EQ. 1)THEN
	    XM(1)=IBOUND
	  ELSE IF(NI .EQ. 2)THEN
	    E1=EXP(-DTAU(1))
	    E2=1.0_LDP-(1.0_LDP-E1)/DTAU(1)
	    E3=(1.0_LDP-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0E-03_LDP)THEN
	      E2=DTAU(1)*0.5_LDP+DTAU(1)*DTAU(1)/6.0_LDP
	      E3=DTAU(1)*0.5_LDP-DTAU(1)*DTAU(1)/3.0_LDP
	    END IF
C
	    XM(2)=IBOUND*E1+SOURCE(2)*E2+SOURCE(1)*E3
            XM(1)=0.5_LDP*(IBOUND+XM(2)*E1+SOURCE(1)*E2+SOURCE(2)*E3)
	  END IF
C
C UPDATE THE FA AND FB MATRICES . (SEE NOTES)
C
	  DO I=1,NI
	    NEWRJ(I)=NEWRJ(I)+AQW(I,LS)*XM(I)
	    NEWRK(I)=NEWRK(I)+AQW3(I,LS)*XM(I)
	  END DO
C
	  HBCNEW=HBCNEW+AQW(1,LS)*(XM(1)-IBOUND)*Z(1)/R(1)
	  IF(NI .EQ. ND)THEN
	    INBCNEW=INBCNEW + AQW(ND,LS)*
	1    (XM(ND)-(XM(ND)-XM(ND-1))/DTAU(ND-1))*Z(ND)/R(ND)
	  END IF
C
2000	CONTINUE
C
C Compute the factor for the outer boundary condition.
C
	HBCNEW=HBCNEW/NEWRJ(1)
	INBCNEW=INBCNEW/(2.0_LDP*NEWRJ(ND)-IC)
C
C Compute the new Feautrier factors.
C
	DO I=1,ND
	  NEWRK(I)=NEWRK(I)/NEWRJ(I)
	END DO
C
C Compute the Q factors from F. The Q values are stored in NEWRJ.
C We now compute Q in JFEAUNEW so that we nolonger need to save Q (only f).
C We still compute Q here because
C    (i) Q is passed (as NEWRJ vector) and hence is corrupted by this
C        routine. We still need Q, however, for the varaition
C        routines.
C   (ii) Q is now completely consistent with the updated f value.
C
	CALL QFROMF(NEWRK,NEWRJ,R,TA,TB,ND)	!TA work vector
C
	RETURN
	END
