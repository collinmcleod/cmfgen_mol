	SUBROUTINE SET_DC_OR_POP_OR_TX_V2(YV,LEV,HE2,LOG_HE2LTE,EDGEHE2,NHE2,ND,T,X,DESC,FLAG)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 11-Apr-2023 : Added DCRGS option.
! Altered 10-Dec-2010 : Changed to V2
!                       Pass LOG_He2LTE rather than He2LTE.
! Altered 20-Mar-1997 : YV is now REAL(KIND=LDP)
!
	INTEGER NHE2,ND,LEV,I,J
!
	REAL(KIND=LDP) YV(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) EDGEHE2(NHE2)
	REAL(KIND=LDP) HE2(NHE2,ND)
	REAL(KIND=LDP) LOG_HE2LTE(NHE2,ND)
!
	CHARACTER*(*) X,DESC
	CHARACTER*20 XSPEC
	LOGICAL FLAG,LOCAL
	INTEGER, PARAMETER :: T_OUT=6
!
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) DELTA_T
	REAL(KIND=LDP) T_EXCITE
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
	LOCAL=.FALSE.
	I=LEN_TRIM(X)
	J=INDEX(X,'_')
	IF(J .EQ. 0)THEN
	    WRITE(T_OUT,*)'Error in SET_DC_OR_POP - missing _ from option'
	    WRITE(T_OUT,*)'NSPEC=',X
	    RETURN
	ELSE
	  XSPEC=X(J+1:)
	END IF
!	
	IF(X(1:5) .EQ. 'DCRGS' .AND. XSPEC .EQ. DESC)THEN
	  IF(LEV .LE. 0 .OR. LEV .GT. NHE2)THEN
	    WRITE(T_OUT,*)'Error in SET_DC_OR_POP_OR_TX - invalid level'
	    WRITE(T_OUT,*)'Lev=',LEV
	    WRITE(T_OUT,*)'NSPEC=',NHE2
	    RETURN
	  END IF
	  DO J=1,ND
	    YV(J)=LOG10(HE2(LEV,J)/HE2(1,J)) -
	1             (LOG_HE2LTE(LEV,J)-LOG_HE2LTE(1,J))/LOG(10.0_LDP)
	  END DO
	  LOCAL=.TRUE.
	ELSE IF(X(1:2) .EQ. 'DC' .AND. XSPEC .EQ. DESC)THEN
	  IF(LEV .LE. 0 .OR. LEV .GT. NHE2)THEN
	    WRITE(T_OUT,*)'Error in SET_DC_OR_POP_OR_TX - invalid level'
	    WRITE(T_OUT,*)'Lev=',LEV
	    WRITE(T_OUT,*)'NSPEC=',NHE2
	    RETURN
	  END IF
	  DO J=1,ND
	    YV(J)=LOG10(HE2(LEV,J))-LOG_HE2LTE(LEV,J)/LOG(10.0_LDP)
	  END DO
	  LOCAL=.TRUE.
	ELSE IF(X(1:2) .EQ. 'TX' .AND. XSPEC .EQ. DESC)THEN
	  IF(LEV .LE. 0 .OR. LEV .GT. NHE2)THEN
	    WRITE(T_OUT,*)'Error in SET_DC_OR_POP_OR_TX - invalid level'
	    WRITE(T_OUT,*)'Lev=',LEV
	    WRITE(T_OUT,*)'NSPEC=',NHE2
	    RETURN
	  END IF
	  T_EXCITE=T(ND)
	  I=LEV
	  DO J=ND,1,-1
	    DELTA_T=100
	    DO WHILE(DELTA_T .GT. 1.0E-08_LDP)
	      T1=EXP( LOG(He2(I,J))-LOG_He2LTE(I,J)+HDKT*EDGEHE2(I)*(1.0_LDP/T(J)-1.0_LDP/T_EXCITE) )*
	1           (T_EXCITE/T(J))**1.5_LDP
	      DELTA_T=(T1-1.0_LDP)*T_EXCITE/T1/(1.5_LDP+HDKT*EDGEHE2(I)/T_EXCITE)
	      T_EXCITE=T_EXCITE-DELTA_T
	    END DO
	    YV(J)=T_EXCITE
	  END DO
	  LOCAL=.TRUE.
	ELSE IF(X(1:3) .EQ. 'RAT' .AND. LEV .GT. NHE2
	1      .AND. XSPEC .EQ. DESC)THEN
	  DO J=1,ND
	    T1=0.0_LDP		!By using T1 advoids log(ZERO).
	    DO I=1,NHE2
	      T1=T1+HE2(I,J)
	    END DO
	    YV(J)=LOG10(T1)
	  END DO
	  LOCAL=.TRUE.
	ELSE IF(XSPEC .EQ. DESC)THEN  		!RAT or POP
	  IF(LEV .LE. 0 .OR. LEV .GT. NHE2)THEN
	    WRITE(T_OUT,*)'Error in SET_DC_OR_POP_OR_TX - invalid level'
	    WRITE(T_OUT,*)'Lev=',LEV
	    WRITE(T_OUT,*)'NSPEC=',NHE2
	    RETURN
	  END IF
	  DO J=1,ND
	    YV(J)=LOG10(HE2(LEV,J))
	  END DO
	  LOCAL=.TRUE.
	END IF
!
	IF(LOCAL .AND. FLAG)THEN
	  WRITE(T_OUT,*)'Error - POP, DC or RAT has been called twice'
	ELSE IF(LOCAL)THEN
	  FLAG=LOCAL
	END IF
!
	RETURN
	END
