!
! Routine to compute the radial optical depth scale. DTAU and dCHI_dR
! are work vectors.
!
	SUBROUTINE TORSCL_V3(TOR,CHI,R,DTAU,dCHI_dR,ND,METHOD,
	1             TYPE_ATM,INDX,WR_ASSUMED_ATM)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 27-Jan-2014 : Changed to V3, INDEX inserted in call. Modified
!                         power law options. PCOMP added.
! Altered 14-May-2010 : WR_ASSUMED_ATM inserted into call.
!                         CHanged to V2.
! Altered 21-Dec-2004 : Bug fix. For TYPE_ATM .NE. 'EXP', the routine was always
!                         returning TAU(1) = CHI(1)*R(1). Routine now computes
!                         optical depth over depth indices 1 to 5, and limits the
!                         power law variation (CHI propto r^{-n}) to n > 1.5.
! Altered 02-Jul-1998 : Length of TYPE_ATM checked to avoid bounds problem
!                        when TYPE_ATM is passed as a single (blank) character.
! Altered 28-May-1996 : Removed for [jdh.disp]SETVEC routine.
!                       DOUBLE PRECISION declarations removed.
!                       ERROR_LU installed.
!
! Altered 11-Nov-88 - TYPE_ATM store as option. Alows better
!                     estimate of optical depth to infinity.
!
	INTEGER ND
	INTEGER INDX
	REAL(KIND=LDP) TOR(ND),CHI(ND),R(ND),DTAU(ND),dCHI_dR(ND)
	LOGICAL WR_ASSUMED_ATM
	CHARACTER*(*) METHOD,TYPE_ATM
!
	INTEGER I
	REAL(KIND=LDP) T1
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	CALL DERIVCHI(dCHI_dR,CHI,R,ND,METHOD)
	CALL NORDTAU(DTAU,CHI,R,R,dCHI_dR,ND)
!
! Determine optical depth from outer boundary to infinity.
! For the extrapolations, we use depths 1 to INDX.
!
	LUER=ERROR_LU()
	IF(TYPE_ATM(1:MIN(3,LEN(TYPE_ATM))) .EQ. 'EXP')THEN
	  IF(CHI(1) .GT. 0 .AND. CHI(INDX) .GT. CHI(1))THEN
	   TOR(1)=CHI(1)*(R(1)-R(INDX))/LOG(CHI(INDX)/CHI(1))
	  ELSE
	    TOR(1)=CHI(1)/10.0D0
	    WRITE(LUER,*)'Warning - optical depth at boundary set to',TOR(1),' in TORSCL'
	  END IF
!
	ELSE IF (TYPE_ATM(1:5) .EQ. 'PCOMP')THEN
	  IF(CHI(1) .GT. 0 .AND. CHI(INDX) .GT. CHI(1))THEN
	    T1=LOG(CHI(1)/CHI(INDX))/LOG(R(INDX)/R(1))
	    IF(T1 .GT. 1.5D0)THEN
	      TOR(1)=CHI(1)*R(1)/(T1-1.0D0)
	    ELSE
	      TOR(1)=CHI(1)*R(1)
	      IF(WR_ASSUMED_ATM)WRITE(LUER,*)'Warning - opacity assumed to be r**(-2) in TORSCL'
	    END IF
	  ELSE
	    TOR(1)=ABS(CHI(1))*R(1)
	    WRITE(LUER,*)'Negative opacity encountered in TORSCL_V3'
	    WRITE(LUER,*)'Warning - opacity assumed to be r**(-2) in TORSCL'
	    WRITE(LUER,*)'CHI(1)=',CHI(1)
	    WRITE(LUER,*)'CHI(INDX)=',CHI(INDX),INDX
	  END IF
!
! For this option, we assume the pwoer-law exponent is passed.
!
        ELSE IF(TYPE_ATM(1:1) .EQ. 'P')THEN
          READ(TYPE_ATM(2:),*)T1
          TOR(1)=CHI(1)*R(1)/(T1-1.0D0)
!
	ELSE IF(TYPE_ATM .EQ. ' ')THEN
	  TOR(1)=CHI(1)*R(1)
	  IF(WR_ASSUMED_ATM)WRITE(LUER,*)'Warning - opacity assumed to be r**(-2) in TORSCL'
	ELSE
	  WRITE(LUER,*)'Atmosphere type:',TRIM(TYPE_ATM),' not recognized in TOR_SCL_V3'
	  STOP
	END IF
!
	DO I=2,ND
	  TOR(I)=TOR(I-1)+DTAU(I-1)
	END DO
!
	RETURN
	END
