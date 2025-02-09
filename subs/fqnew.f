C
C Subroutine to update dCHI and dETA matrices to variation of
C mean continuum intensity in the electron scattering emmission
C coefficient.
C
C
	SUBROUTINE FQNEW(FQAF,F2DA,FC,W,ESEC,ND,NM)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 22-Jun-2022 Removed shared DO end.
C Altered 05-Dec-1996 DO 100 I(J) replaced by DO, END DO
C Altered 24-May-1996  IMPLICIT NONE installed.
C Altered  7-SEP-1982 (%CHI changed to 3 and %ETA changed to 4)
C Altered 31-AUG-1982 (Major bug fixed)
C Altered 23-AUG-1982 (FC changed to pure dETA matrix)
C Created 27-JUL-1982
C
	INTEGER NM,ND
	REAL(KIND=LDP) FQAF(ND,ND,NM),F2DA(ND,ND),FC(ND,ND)
	REAL(KIND=LDP) W(ND,ND),ESEC(ND)
C
	INTEGER I,J,K
C
	DO J=1,ND
	  DO I=1,ND
	    W(I,J)=FQAF(I,J,4)*ESEC(J)
	  END DO
	END DO
C
	DO K=1,ND
	  DO J=1,ND
	    DO I=1,ND
	      FQAF(J,K,3)=FQAF(J,K,3)+W(J,I)*F2DA(I,K)
	      FQAF(J,K,4)=FQAF(J,K,4)+W(J,I)*FC(I,K)
	    END DO
	  END DO
	END DO
C
	RETURN
	END
