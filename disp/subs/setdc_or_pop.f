	SUBROUTINE SETDC_OR_POP(YV,LEV,HE2,HE2LTE,NHE2,ND,X,DESC,FLAG)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Aletered 20-Mar-1997 : YV is now REAL(KIND=LDP)
C
	INTEGER NHE2,ND,LEV,I,J
	REAL(KIND=LDP) YV(ND)
	REAL(KIND=LDP) HE2(NHE2,ND),HE2LTE(NHE2,ND),T1
	CHARACTER*(*) X,DESC
	CHARACTER*20 XSPEC
	LOGICAL FLAG,LOCAL
	INTEGER, PARAMETER :: T_OUT=6
C
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
C	
	IF(X(1:2) .EQ. 'DC' .AND. XSPEC .EQ. DESC)THEN
	  IF(LEV .LE. 0 .OR. LEV .GT. NHE2)THEN
	    WRITE(T_OUT,*)'Error in SET_DC_OR_POP - invalid level'
	    WRITE(T_OUT,*)'Lev=',LEV
	    WRITE(T_OUT,*)'NSPEC=',NHE2
	    RETURN
	  END IF
	  DO J=1,ND
	    YV(J)=LOG10(HE2(LEV,J)/HE2LTE(LEV,J))
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
	ELSE IF(XSPEC .EQ. DESC )THEN  		!RAT or POP
	  IF(LEV .LE. 0 .OR. LEV .GT. NHE2)THEN
	    WRITE(T_OUT,*)'Error in SET_DC_OR_POP - invalid level'
	    WRITE(T_OUT,*)'Lev=',LEV
	    WRITE(T_OUT,*)'NSPEC=',NHE2
	    RETURN
	  END IF
	  DO J=1,ND
	    YV(J)=LOG10(HE2(LEV,J))
	  END DO
	  LOCAL=.TRUE.
	END IF
C
	IF(LOCAL .AND. FLAG)THEN
	  WRITE(T_OUT,*)'Error - POP, DC or RAT has been called twice'
	ELSE IF(LOCAL)THEN
	  FLAG=LOCAL
	END IF
C
	RETURN
	END
