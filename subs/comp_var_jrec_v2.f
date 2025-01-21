	SUBROUTINE COMP_VAR_JREC_V2(JREC,dJRECdT,JPHOT,
	1                  JREC_CR,dJREC_CRdT,JPHOT_CR,BPHOT_CR,
	1                  RJ,EMHNUKT,T,NU,FQW,TWOHCSQ,HDKT,ND,COMP_NEW_CONT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 14-Jul-2019 : Added dJREC_CRdT to call and hanged to V2
!
	INTEGER ND
	REAL(KIND=LDP) JREC(ND)
	REAL(KIND=LDP) dJRECdT(ND)
	REAL(KIND=LDP) JPHOT(ND)
	REAL(KIND=LDP) JREC_CR(ND)
	REAL(KIND=LDP) dJREC_CRdT(ND)
	REAL(KIND=LDP) JPHOT_CR(ND)
	REAL(KIND=LDP) BPHOT_CR(ND)
	REAL(KIND=LDP) RJ(ND)
	REAL(KIND=LDP) EMHNUKT(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) NU
	REAL(KIND=LDP) FQW
	REAL(KIND=LDP) HDKT,TWOHCSQ
	LOGICAL COMP_NEW_CONT
C
	REAL(KIND=LDP) T1,T2
	INTEGER I
C
	IF(COMP_NEW_CONT)THEN
	  JREC(:)=0.0_LDP
	  dJRECdT(:)=0.0_LDP
	  JPHOT(:)=0.0_LDP
	  JREC_CR(:)=0.0_LDP
	  dJREC_CRdT(:)=0.0_LDP
	  JPHOT_CR(:)=0.0_LDP
	  BPHOT_CR(:)=0.0_LDP
	END IF
C
	T1=TWOHCSQ*(NU**3)
	DO I=1,ND
	  T2=(T1+RJ(I))*EMHNUKT(I)*FQW/NU
	  JREC(I)=JREC(I)+T2
	  dJRECdT(I)=dJRECdT(I)+T2*HDKT*NU/T(I)/T(I)
	  JPHOT(I)=JPHOT(I)+RJ(I)*FQW/NU
	  JREC_CR(I)=JREC_CR(I)+T2*NU
	  dJREC_CRdT(I)=dJREC_CRdT(I)+T2*NU*HDKT*NU/T(I)/T(I)
	  JPHOT_CR(I)=JPHOT_CR(I)+RJ(I)*FQW
	  BPHOT_CR(I)=BPHOT_CR(I)+T1*FQW*EMHNUKT(I)/(1.0_LDP-EMHNUKT(I))
	END DO
C
	RETURN
	END
