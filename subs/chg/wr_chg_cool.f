!
! Routine to write out cooling/heating for charge exchange reactions.
!
	SUBROUTINE WR_CHG_COOL(NETCR,TOTCR,COUNTER,ND,LU)
	USE SET_KIND_MODULE
	USE CHG_EXCH_MOD
	IMPLICIT NONE
!
! Created 26-Aug-1998 : based on WR_AD_COOL
!
	INTEGER*4 COUNTER
	INTEGER*4 ND,LU
	REAL(KIND=LDP) NETCR(ND)	!Accumlated net cooling rate
	REAL(KIND=LDP) TOTCR(ND)	!Accumulated sum of absoulte cooling rates.
!
	INTEGER*4 I,J
	INTEGER*4 MS,MF
!
	MS=(COUNTER-1)*10+1	
	MF=COUNTER*10
	IF(MF .GT. ND)MF=ND
!
	IF(.NOT. DO_CHG_EXCH)RETURN
!
	WRITE(LU,'(/3X,A)')'Charge exchange cooling rate [ergs/cm**3/sec]'
	DO J=1,N_CHG
	  IF(CHG_REACTION_AVAILABLE(J))THEN
	    DO I=MS,MF
	      NETCR(I)=NETCR(I)+COOL_CHG(J,I)
	      TOTCR(I)=TOTCR(I)+ABS(COOL_CHG(J,I))
	    END DO
	    WRITE(LU,999)(COOL_CHG(J,I),I=MS,MF)
	  END IF
	END DO
!
999	FORMAT(X,1P,10E12.4)
	RETURN
	END
