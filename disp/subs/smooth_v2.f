	SUBROUTINE SMOOTH_V2(SPEC,LTE,VAL,N,TYPE,DEPTH_LIM,
	1                     ND,OPT,IDENT,FLAG)
	IMPLICIT NONE
C
	INTEGER*4 N,ND,DEPTH_LIM
	REAL*8 SPEC(N,ND),LTE(N,ND),VAL
	CHARACTER*(*) IDENT,OPT,TYPE
	LOGICAL FLAG
C
C Local variables
C
	INTEGER*4 I,J
	REAL*8 T1,T2
C
	IF( (OPT .NE. IDENT) .AND. (OPT .NE. 'ALL') )RETURN
C
	IF(TYPE .EQ. 'POP')THEN
	  FLAG=.TRUE.
	  DO J=DEPTH_LIM-1,1,-1
	    DO I=1,N
	      IF(SPEC(I,J) .GT. SPEC(I,J+1)*VAL)THEN
                SPEC(I,J)=MAX(SPEC(I,J)*0.01D0,SPEC(I,J+1)*VAL)
	      ELSE IF(SPEC(I,J) .LT. SPEC(I,J+1)/VAL)THEN
                SPEC(I,J)=MIN(SPEC(I,J)*100.0D0,SPEC(I,J+1)/VAL)
	      END IF
	    END DO
	  END DO
	ELSE IF(TYPE .EQ. 'DC')THEN
	  IF(LTE(1,ND) .EQ. 0)RETURN
	  FLAG=.TRUE.
	  DO J=DEPTH_LIM-1,1,-1
	    DO I=1,N
	      T1=SPEC(I,J)/LTE(I,J)
	      T2=SPEC(I,J+1)/LTE(I,J+1)
	      IF(T1 .GT. T2*VAL)THEN
                  SPEC(I,J)=LTE(I,J)*VAL
	      ELSE IF(T1 .GT. T2/VAL)THEN
                  SPEC(I,J)=LTE(I,J)*VAL
	      END IF
	    END DO
	  END DO
	ELSE
	  WRITE(6,*)'Invalid type in SMOOTH_V2'
	END IF
C
	RETURN
	END
