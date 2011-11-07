	SUBROUTINE ARNAUD_CROSS_V3()
	USE MOD_NON_THERM
	USE MOD_CMFGEN
	IMPLICIT NONE
!
!	INTEGER NKT
!	REAL*8 XKT(NKT)
!
	INTEGER IT
	INTEGER IKT
	REAL*8 U1
	REAL*8 T1,T2,T3
	INTEGER ID
!	INTEGER MAX_NUM_IONS
	LOGICAL NEGATIVE_CROSS_SEC
!	CHARACTER*12 ION_ID(MAX_NUM_IONS)
!
! Fitting formula from Arnaud & Rothenflug 1985:
! http://adsabs.harvard.edu/abs/1985A%26AS...60..425A
! Gives the cross section for direct ionisation.
!
	DO IT=1,NUM_THD
	  IF(THD(IT)%PRES)THEN
	    NEGATIVE_CROSS_SEC=.FALSE.
	    THD(IT)%CROSS_SEC=0.0D0
	    ID=THD(IT)%LNK_TO_ION
	    DO IKT=1,NKT
	      U1 = XKT(IKT) / THD(IT)%ION_POT
	      IF (U1 .GT. 1.0D0) THEN
	        T1=(1.0D0-1.0D0/U1)
	        T2=DLOG(U1)
	        T3 = 1.0D-14 * ( THD(IT)%A_COL*T1 + THD(IT)%B_COL*T1*T1 + &
	              THD(IT)%C_COL*T2 + THD(IT)%D_COL*T2/U1  ) / U1 / THD(IT)%ION_POT**2
	        T3 = T3*ATM(ID)%ION_CROSEC_NTFAC
	        IF(T3 .GE. 0.0D0)THEN
	          THD(IT)%CROSS_SEC(IKT)=T3
	        ELSE
	          THD(IT)%CROSS_SEC(IKT)=0.0D0
	          NEGATIVE_CROSS_SEC=.TRUE.
	        END IF
	      END IF
	    END DO
	    IF(NEGATIVE_CROSS_SEC)THEN
	      ID=THD(IT)%LNK_TO_ION
	      WRITE(6,*)'Negative cross section(s) for ',ION_ID(ID),' in ARNAUD_CROSS_V2!'
	      WRITE(6,*)'Forcing them to be zeros.'
	    END IF
	  END IF
	END DO
!	
	RETURN
	END
