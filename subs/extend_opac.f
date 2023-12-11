C
C
	SUBROUTINE EXTEND_OPAC(CHIEXT,ETAEXT,ESECEXT,RJEXT,COEF,INDX,
	1                     NX,CHI,ETA,ESEC,RJ,ND)
	USE SET_KIND_MODULE
C
	IMPLICIT NONE
	INTEGER NX,ND,I,J,INDX(NX)
	REAL(KIND=LDP) CHIEXT(NX),ETAEXT(NX),ESECEXT(NX),RJEXT(NX),COEF(0:3,NX)
	REAL(KIND=LDP) CHI(ND),ETA(ND),ESEC(ND),RJ(ND)
C
	DO I=1,NX
	  CHIEXT(I)=0.0_LDP
	  ETAEXT(I)=0.0_LDP
	  ESECEXT(I)=0.0_LDP
	  RJEXT(I)=0.0_LDP
	  DO J=0,3
	    CHIEXT(I)=CHIEXT(I)+COEF(J,I)*LOG( CHI(J+INDX(I)) )
	    ETAEXT(I)=ETAEXT(I)+COEF(J,I)*LOG( ETA(J+INDX(I)) )
	    ESECEXT(I)=ESECEXT(I)+COEF(J,I)*LOG( ESEC(J+INDX(I)) )
	    RJEXT(I)=RJEXT(I)+COEF(J,I)*LOG( RJ(J+INDX(I)) )
	  END DO
	  CHIEXT(I)=EXP(CHIEXT(I))
	  ETAEXT(I)=EXP(ETAEXT(I))
	  ESECEXT(I)=EXP(ESECEXT(I))
	  RJEXT(I)=EXP(RJEXT(I))
	END DO
C
	RETURN	
	END
