!
! Subroutine to fits a straigt line, and a set of Gaussians,
! to a set of data. Routine is designed to be called in GRAMON_PGPLOT.
!
	SUBROUTINE DO_GAUS_FIT(XST_PASSED,XEND_PASSED)
	USE MOD_CURVE_DATA
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	REAL*4 XST_PASSED,XEND_PASSED
	REAL*8 GAUS_FIT_FUNC
	EXTERNAL GAUS_FIT_FUNC
!
	REAL*8 XST,XEND
	REAL*8 YST,YEND
!
	REAL*8 TOL
	REAL*8 T1
	INTEGER IP
	INTEGER NVALS
	INTEGER ITER
	INTEGER, SAVE :: N_PAR
	INTEGER, SAVE :: NUM_GAUS
	INTEGER, SAVE :: N_PAR_OLD
!
	INTEGER I,J,K
!
	XST=XST_PASSED
	XEND=XEND_PASSED
	CALL GEN_IN(XST,'Start wavlength for fitting')
	CALL GEN_IN(XEND,'End wavelength for fitting')
	IP=1
	CALL GEN_IN(IP,'Plot for fitting')
!
	CALL GEN_IN(NUM_GAUS,'Number of gaussians to fit:')
!
	IF(ALLOCATED(SIM))THEN
	  DEALLOCATE(SIM,SUM_SQ)
	  N_PAR_OLD=N_PAR
	END IF
	N_PAR=2+3*NUM_GAUS
	ALLOCATE (SIM(N_PAR+1,N_PAR))
	ALLOCATE (SUM_SQ(N_PAR+1))
!
! If possible, we use previous fit parameters as default.
!
	SIM=0.0D0
	IF(ALLOCATED(PAR))THEN
	  SIM(1,1:MIN(N_PAR,N_PAR_OLD))=PAR(1:MIN(N_PAR,N_PAR_OLD))
	  DEALLOCATE(PAR)
	ELSE
	  SIM(1,1)=1.0D0
	  SIM(1,2)=0.0D0	
	ENDIF
	ALLOCATE (PAR(N_PAR))
	CALL GEN_IN(SIM(1,1),'Mean value')
	CALL GEN_IN(SIM(1,2),'Continuum slope')
	DO J=1,NUM_GAUS
	  K=2+(J-1)*3+1
	  CALL GEN_IN(SIM(1,K),'Central wavlength of Gaussian')
          IF(SIM(1,K+1) .EQ. 0.0D0 .AND. J .GT. 1)SIM(1,K+1)=SIM(1,K-2)
	  CALL GEN_IN(SIM(1,K+1),'Sigma of Gassian')
	  IF(SIM(1,K+2) .EQ. 0.0D0)SIM(1,K+2)=-0.2D0
	  CALL GEN_IN(SIM(1,K+2),'Offset from continuum (-ve for absorption)')
	END DO
!
	WRITE(6,*)CD(IP)%XVEC(1),CD(IP)%XVEC(NPTS(IP))
	WRITE(6,*)CD(IP)%DATA(1),CD(IP)%DATA(NPTS(IP))
	CALL SET_GAUS_DATA(CD(IP)%XVEC,CD(IP)%DATA,XST,XEND,
	1                    NPTS(IP),N_PAR,YST,YEND,NVALS)
!
! Set other vertices.
!
	DO I=2,N_PAR+1
	  SIM(I,1:N_PAR)=SIM(1,1:N_PAR)
	END DO
	SIM(2,1)=0.9D0*SIM(1,1)
	SIM(3,2)=(YEND-SIM(I,1))/(XEND-XST)
	DO J=1,NUM_GAUS
	  K=3+(J-1)*3
	  SIM(K+1,K)=SIM(1,K)+SIM(1,K+1)
	  SIM(K+2,K+1)=SIM(1,K+1)*0.9D0
	  SIM(K+3,K+2)=SIM(1,K+2)*0.9D0
	END DO
	WRITE(6,*)'Simples vertices have been set'
	DO J=1,N_PAR
	  WRITE(6,'(20ES12.4)')(SIM(I,J),I=1,N_PAR+1)
	END DO
! 
        SUM_SQ(:)=0.0D0
        DO I=1,N_PAR+1
          PAR(1:N_PAR)=SIM(I,1:N_PAR)
          SUM_SQ(I)=GAUS_FIT_FUNC(PAR)
        END DO
	WRITE(6,*)'Evaluated SUM_SQ'
!
        TOL=1.0D-08
	I=N_PAR+1
        CALL AMOEBA(SIM,SUM_SQ,I,N_PAR,N_PAR,TOL,GAUS_FIT_FUNC,ITER)
	WRITE(6,*)'Called AMOEBA'
!
	WRITE(6,*)'Fit parameters are:'
	WRITE(6,*)SIM(1,1),SIM(1,2)
	WRITE(6,'(7X,A,3X,5X,A,6X,6X,A,7X,3X,A,4X,A)')
	1          'Height','Lam','a','Sigma(km/s)','FWHM(km/s)'
	DO I=1,NUM_GAUS
	  K=2+(I-1)*3+1
	  T1=2.99794D+05*SIM(1,K+1)/SIM(1,K)
	  WRITE(6,'(2X,5ES14.4)')SIM(I,K+2),SIM(1,K),SIM(1,K+1),
	1               T1/SQRT(2.0D0),1.6651*T1
	END DO
!
        SUM_SQ(:)=0.0D0
        PAR(1:N_PAR)=SIM(1,1:N_PAR)
        SUM_SQ(1)=GAUS_FIT_FUNC(PAR)
	WRITE(6,*)'Final SUM_SQ=',SUM_SQ(1)
	I=1
	CALL PGSCI(I)
	CALL PGLINE(NG_DATA,XFIT,YFIT)
!	
	RETURN
	END
!
! Routine to draw the fit on top of an existing plot. At present
! plot is drawn in black only.
!
	SUBROUTINE DRAW_GAUS
	USE MOD_CURVE_DATA
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	REAL*8 GAUS_FIT_FUNC
	EXTERNAL GAUS_FIT_FUNC
	INTEGER I
!
        SUM_SQ(:)=0.0D0
        PAR(1:NG_PAR)=SIM(1,1:NG_PAR)
        SUM_SQ(1)=GAUS_FIT_FUNC(PAR)
	I=1
	CALL PGSCI(I)
	CALL PGLINE(NG_DATA,XFIT,YFIT)
!
	RETURN
	END
