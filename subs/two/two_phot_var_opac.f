!
! Routine to compute the variation of CHI and ETA due to 2-photon
! processes. We ignore any temperature dependence (which may arise
! from the use of super levels).
!
	SUBROUTINE TWO_PHOT_VAR_OPAC(VETA,VCHI,POPS,T,FREQ,ND,NT)
	USE SET_KIND_MODULE
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Created 26-Jun-1998
!
	INTEGER*4 NT,ND
	REAL(KIND=LDP) VETA(NT,ND)
	REAL(KIND=LDP) VCHI(NT,ND)
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) T(ND)
!
	REAL(KIND=LDP) FREQ
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local constants.
!
	REAL(KIND=LDP) h
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) CONST
	REAL(KIND=LDP) ETA_CONST
	REAL(KIND=LDP) CHI_CONST
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) FREQ_B
!
! See TWO_PHOT_OPAC for definitions
!
	REAL(KIND=LDP) AY,Y,U,FU
!
	INTEGER*4 LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER*4 J,L
	INTEGER*4 NL,NUP
	
	h=6.626D-27			!cgs units
	PI=4.0D0*ATAN(1.0D0)
	CONST=1.0D+10*H/4.0D0/PI
!
	DO J=1,N_TWO
	  IF(TWO_PHOT_AVAILABLE(J) .AND. FREQ .LT. FREQ_TWO(J))THEN
!
	    FREQ_B=FREQ_TWO(J)-FREQ
	    NL=LOW_LEV_TWO(J)
	    NUP=UP_LEV_TWO(J)
	    IF(TYPE_TWO(J) .EQ. 1)THEN
	      Y=FREQ/FREQ_TWO(J)
	      U=Y*(1.0D0-Y)
	      FU=4.0D0*U
	      AY=24.56D0*COEF_TWO(J,1)*( U*(1-FU**0.8) +
	1                 0.88D0*(U**1.53)*(FU**0.8) )
	    END IF
!
	    ETA_CONST=CONST*FREQ/FREQ_TWO(J)
	    CHI_CONST=G_UP_TWO(J)*ETA_CONST/TWOHCSQ/FREQ**3
	    DO L=1,ND
	      VETA(NUP,L)=VETA(NUP,L) + ETA_CONST*AY*FS_RAT_UP(L,J)
	      T1=EXP(-HDKT*FREQ_B/T(L))
	      VCHI(NL,L)=VCHI(NL,L) + CHI_CONST*AY*
	1                   FS_RAT_LOW(L,J)*T1/G_LOW_TWO(J)
	      VCHI(NUP,L)=VCHI(NUP,L) - CHI_CONST*AY*
	1                   *FS_RAT_UP(L,J)/G_UP_TWO(J)
	    END DO
	  END IF
	END DO
!
	RETURN
	END
