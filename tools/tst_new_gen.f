	PROGRAM TST_NEW_GEN_IN
	USE SET_KIND_MODULE
	USE NEW_GEN_IN_INTERFACE
	IMPLICIT NONE
C
	INTEGER, PARAMETER :: NMAX=10
	REAL(KIND=LDP) VEC(NMAX)
	REAL*4 SP(NMAX)
	INTEGER INTV(NMAX)
	INTEGER NSP,NDP,NI,I
	CHARACTER*20 HOPE
	LOGICAL LOG_VAL
C
	DO I=1,NMAX
	  SP(I)=I*I+0.1
	  VEC(I)=I*I+0.3
	  INTV(I)=I*I
	END DO
	HOPE='  TEST'
C
	CALL NEW_GEN_IN_OPTS('OPEN_LOG_FILE','LG_FILE',11)
	LOG_VAL=.TRUE.
	CALL NEW_GEN_IN(LOG_VAL,'Logical value')
	CALL NEW_GEN_IN(SP,NSP,NMAX,'Single precision Values')
	CALL NEW_GEN_IN(VEC,NDP,NMAX,'Double precision Values')
	CALL NEW_GEN_IN(INTV,NI,NMAX,'Integer Values')
	CALL NEW_GEN_IN(HOPE,'String value')
	CALL NEW_GEN_IN_OPTS('CLOSE_LOG_FILE','LG_FILE',11)
C
	WRITE(6,*)'Logical value=',LOG_VAL
	WRITE(6,*)'String value=',TRIM(HOPE)
	DO I=1,NSP
	  WRITE(6,*)I,SP(I)
	END DO
	DO I=1,NDP
	  WRITE(6,*)I,VEC(I)
	END DO
	DO I=1,NI
	  WRITE(6,*)I,INTV(I)
	END DO
!
	CALL NEW_GEN_IN_OPTS('OPEN_FILE_INPUT','LG_FILE',12)
	LOG_VAL=.TRUE.
	CALL NEW_GEN_IN(LOG_VAL,'Logical value')
	CALL NEW_GEN_IN(SP,NSP,NMAX,'Single precision Values')
	CALL NEW_GEN_IN(VEC,NDP,NMAX,'Double precision Values')
	CALL NEW_GEN_IN(INTV,NI,NMAX,'Integer Values')
!
	WRITE(6,*)'Logical value=',LOG_VAL
	WRITE(6,*)'String value=',TRIM(HOPE)
	DO I=1,NSP
	  WRITE(6,*)I,SP(I)
	END DO
	DO I=1,NDP
	  WRITE(6,*)I,VEC(I)
	END DO
	DO I=1,NI
	  WRITE(6,*)I,INTV(I)
	END DO
!
	STOP
	END
