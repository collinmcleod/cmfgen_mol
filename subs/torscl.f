C
C Routine to compute the radial optical depth scale. DTAU and dCHI_dR
C are work vectors.
C
	SUBROUTINE TORSCL(TOR,CHI,R,DTAU,dCHI_dR,ND,METHOD,TYPE_ATM)
	IMPLICIT NONE
C
C Altered 02-Jul-1998 : Length of TYPE_ATM checked to avoid bounds problem
C                        when TYPE_ATM is passed as a single (blank) character.
C Altered 28-May-1996 : Removed for [jdh.disp]SETVEC routine.
C                       DOUBLE PRECISION declarations removed.
C                       ERROR_LU installed.
C
C Altered 11-Nov-88 - TYPE_ATM store as option. Alows better
C                     estimate of optical depth to infinity.
C
	INTEGER*4 ND,I
	REAL*8 TOR(ND),CHI(ND),R(ND),DTAU(ND),dCHI_dR(ND),T1
	CHARACTER*(*) METHOD,TYPE_ATM
C
	INTEGER*4 ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	CALL DERIVCHI(dCHI_dR,CHI,R,ND,METHOD)
	CALL NORDTAU(DTAU,CHI,R,R,dCHI_dR,ND)
C
C Determine optical depth from outer boundary to infinity. 
C
	LUER=ERROR_LU()
	IF(TYPE_ATM(1:MIN(3,LEN(TYPE_ATM))) .EQ. 'EXP')THEN
	  IF(CHI(1) .GT. 0 .AND. CHI(3) .GT. CHI(1))THEN
	   TOR(1)=CHI(1)*(R(1)-R(3))/LOG(CHI(3)/CHI(1))
	  ELSE
	    TOR(1)=0.00001
	    WRITE(LUER,*)'Warning - optical depth at boundary set to 10^{-5}'
	  END IF
	ELSE IF (CHI(1) .GT. 0 .AND. CHI(3) .GT. CHI(1))THEN
	  T1=DLOG(CHI(1)/CHI(3))/DLOG(R(1)/R(3))
	  IF(T1 .GT. 1)THEN
	   TOR(1)=CHI(1)*R(1)/(T1-1.0D0)
	  ELSE
	    TOR(1)=CHI(1)*R(1)
	    WRITE(LUER,*)'Warning - opacity assumed to be r**(-2) '
	  END IF
	ELSE
	  TOR(1)=CHI(1)*R(1)
	  WRITE(LUER,*)'Warning - opacity assumed to be r**(-2) '
	END IF
C
	DO I=2,ND
	  TOR(I)=TOR(I-1)+DTAU(I-1)
	END DO
C
	RETURN
	END
