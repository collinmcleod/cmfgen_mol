!
! Routine to increment the Emissivity and Opacity for
! 2-photon processes.
!
	SUBROUTINE TWO_PHOT_OPAC(ETA,CHI,POPS,T,FREQ,ND,NT)
	USE SET_KIND_MODULE
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Created 26-Jun-1998
!
	INTEGER NT,ND
	REAL(KIND=LDP) ETA(ND)		!Emisivity
	REAL(KIND=LDP) CHI(ND)		!Opacity [in (10^10 cm)^-1]
	REAL(KIND=LDP) POPS(NT,ND)	!Atomic populations
	REAL(KIND=LDP) T(ND)		!in 10^4 K
!
	REAL(KIND=LDP) FREQ		!Current frequency (10^15 Hz)
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local constants.
!
	REAL(KIND=LDP) PLANKS_CONSTANT	!cgs units
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) CONST
	REAL(KIND=LDP) ETA_CONST	!Used to evaluate ETA
	REAL(KIND=LDP) CHI_CONST	!Used to evaluate CHI
	REAL(KIND=LDP) FREQ_B		!Frequency of other photon
	REAL(KIND=LDP) T1
!
! The 2-photon distribution functions, AY, are usually in wrtten in terms
! of the variable y=FREQ/MAX_FREQ, and which extends from 0 to 1.
! As the distribution, AY, is symmetric about y/2 we have defined the
! variables  U=Y(1-Y) and FU=4*U
!
	REAL(KIND=LDP) AY,Y,U,FU
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER J,L
	INTEGER NL,NUP
	
	PLANKS_CONSTANT=6.626E-27_LDP			!cgs units
	PI=4.0_LDP*ATAN(1.0_LDP)
!
! The factor 10^10 arises since R is units of 10^10cm, and we
! scale CHI (and ETA) so that R.CHI is dimensionless, and
! ETA/CHI has the dimensions of the Planck Function.
! We don't have to worry about the frequency units, as it multilplied
! by a ratio of 2 frequencies.
!
	CONST=1.0E+10_LDP*PLANKS_CONSTANT/4.0_LDP/PI
!
	DO J=1,N_TWO
	  IF(TWO_PHOT_AVAILABLE(J) .AND. FREQ .LT. FREQ_TWO(J))THEN
!
	    FREQ_B=FREQ_TWO(J)-FREQ
	    NL=LOW_LEV_TWO(J)
	    NUP=UP_LEV_TWO(J)
!
	    IF(TYPE_TWO(J) .EQ. 1)THEN
	      Y=FREQ/FREQ_TWO(J)
	      U=Y*(1.0_LDP-Y)
	      FU=4.0_LDP*U
	      AY=24.56_LDP*COEF_TWO(J,1)*( U*(1.0_LDP-FU**0.8_LDP) +
	1                 0.88_LDP*(U**1.53_LDP)*(FU**0.8_LDP) )
	    ELSE
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in TWO_PHOT_OPAC'
	      WRITE(LUER,*)'Unrecognized two photon transition'
	      STOP
	    END IF
!
! The expressions for CHI and ETA are approximate. Their ratio gives the
! Planck function at depth, and the expression for ETA is exact at low
! densities and for small dilution factors.
!
	    ETA_CONST=CONST*FREQ/FREQ_TWO(J)
	    CHI_CONST=G_UP_TWO(J)*ETA_CONST/TWOHCSQ/FREQ**3
	    DO L=1,ND
	      ETA(L)=ETA(L) + ETA_CONST*AY*POPS(NUP,L)*FS_RAT_UP(L,J)
	      T1=EXP(-HDKT*FREQ_B/T(L))
	      CHI(L)=CHI(L) + CHI_CONST*AY*(
	1               POPS(NL,L)*FS_RAT_LOW(L,J)*T1/G_LOW_TWO(J)-
	1               POPS(NUP,L)*FS_RAT_UP(L,J)/G_UP_TWO(J) )
	    END DO
	  END IF
	END DO
!
	RETURN
	END
