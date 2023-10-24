!
! Routine to compute the variation of CHI and ETA due to 2-photon
! processes. We ignore any temperature dependence (which may arise
! from the use of super levels).
!
	SUBROUTINE TWO_PHOT_VAR_OPAC(VETA,VCHI,POPS,T,FREQ,ND,NT)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Altered 24-Sep-2023: Updated to use consistent physical constants (24-Sep-223).
! Altered 01-Oct-2015 : Added TWO_METHOD option (passed by module). i
! Created 26-Jun-1998
!
	INTEGER NT,ND
	REAL(10) VETA(NT,ND)
	REAL(10) VCHI(NT,ND)
	REAL(10) POPS(NT,ND)
	REAL(10) T(ND)
!
	REAL(10) FREQ
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL(10) CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local constants.
!
	REAL(10) h
	REAL(10) PI
	REAL(10) CONST
	REAL(10) ETA_CONST
	REAL(10) CHI_CONST
	REAL(10) T1
	REAL(10) FREQ_B
	REAL(10) PLANCKS_CONSTANT
	EXTERNAL PLANCKS_CONSTANT
!
! See TWO_PHOT_OPAC for definitions
!
	REAL(10) AY,Y,U,FU
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER J,L
	INTEGER NL,NUP
	
	h=PLANCKS_CONSTANT()			!cgs units
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
	      AY=24.56D0*COEF_TWO(J,1)*( U*(1.0D0-FU**0.8D0) +
	1                 0.88D0*(U**1.53D0)*(FU**0.8D0) )
	    ELSE
	      WRITE(6,'(/,1X,A)')'Error in TWO_PHOT_VAR_OPAC -- unrecognized type for TWO_PHOTON transition'
	      STOP
	    END IF
!
	    ETA_CONST=CONST*FREQ/FREQ_TWO(J)
	    CHI_CONST=G_UP_TWO(J)*ETA_CONST/TWOHCSQ/FREQ**3
	    IF(TWO_METHOD .EQ. 'OLD_DEFAULT')THEN
	      DO L=1,ND
	        VETA(NUP,L)=VETA(NUP,L) + ETA_CONST*AY*FS_RAT_UP(L,J)
	        T1=EXP(-HDKT*FREQ_B/T(L))      
	        VCHI(NL,L)=VCHI(NL,L)   + CHI_CONST*AY*FS_RAT_LOW(L,J)*T1/G_LOW_TWO(J)
	        VCHI(NUP,L)=VCHI(NUP,L) - CHI_CONST*AY*FS_RAT_UP(L,J)/G_UP_TWO(J)
	      END DO
	    ELSE
	      DO L=1,ND
	        VETA(NUP,L)=VETA(NUP,L) + ETA_CONST*AY*FS_RAT_UP(L,J)*(1.0D0+PHOT_OC_TWO(L,J))
	        VCHI(NL,L)=VCHI(NL,L)   + CHI_CONST*AY*FS_RAT_LOW(L,J)*PHOT_OC_TWO(L,J)/G_LOW_TWO(J)
	        VCHI(NUP,L)=VCHI(NUP,L) - CHI_CONST*AY*FS_RAT_UP(L,J)/G_UP_TWO(J)*(1.0D0+PHOT_OC_TWO(L,J))
	      END DO
	    END IF
	  END IF
	END DO
!
	RETURN
	END
