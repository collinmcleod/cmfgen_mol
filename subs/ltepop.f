C
C Subroutine to compute the LTE populations (at NR depth points)
C given ED (electron density) and DI,GU (density and statistical
C weight of the ground state of the next ionization state 
C respectively).
C
	SUBROUTINE LTEPOP(HNST,ED,DI,G,NUION,T,GU,N,NR)
	IMPLICIT NONE
C
C Altered 05-Dec-1996 : END DO used to terminate DO LOOPS.
C Altered 28-May-1996 : Generica calls for EXP and LOG
C Altered 10-Apr-1989 - Implicit none installed.
C Altered 14-AUG-1984
C
	INTEGER*4 N,NR
	REAL*8 HNST(N,NR),ED(NR),DI(NR),T(NR),G(N),NUION(N),GU
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local varaiables.
C
	INTEGER*4 I,J
	REAL*8 X,Y,RGU
C
	RGU=2.07078D-22/GU
	RGU=LOG(RGU)
	DO I=1,NR
	  X=HDKT/T(I)
	  Y=ED(I)*DI(I)*( T(I)**(-1.5) )
	  DO J=1,N
	    HNST(J,I)=G(J)*Y*EXP(NUION(J)*X+RGU)
	  END DO
	END DO
C
	RETURN
	END
