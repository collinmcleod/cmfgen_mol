!
! Subroutine to compute the temperature structure for a grey atmosphere
! with a Rosseland mean opacity of ROSSMEAN. When passed, ROSSMEAN must include
! the effects (if any) of clumping.
!
! At present, two different options are available.
!
! (a) Spherical atmosphere with out velocity terms.
! (b) If JGREY_WITH_V_TERMS =.TRUE. we use a spherical atmosphere,
!       and include the effect of the first order Doppler shift.
!
	SUBROUTINE COMP_GREY_V4(POPS,TGREY,TAU_ROSS,ROSSMEAN,PLANCKMEAN,COMPUTED,LUER,NC,ND,NP,NT)
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered 27-Jan-2014 : Now use TORSCL_V3.
! Altered 03-Apr-2009 : Inserted OUT_JH call in REL section.
! Altered 20-Jan-2008 : POPS & NT inserted into call. Changed to V2.
!                       Now call JGREY_HUB_DDT_V2
! 
! Altered 30_Aug-2007 : Bug fix. Invalid lower boundary condition flux for
!                         MOM_JREL and for FCOMP_PP_V2i. Forgot factor of CHI. 
!                         
! Created 21-Dec-2004
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
	INTEGER NT
	INTEGER LUER
	LOGICAL COMPUTED
!
	REAL*8 POPS(NT,ND)
	REAL*8 TGREY(ND)
	REAL*8 TAU_ROSS(ND)
	REAL*8 ROSSMEAN(ND)
	REAL*8 PLANCKMEAN(ND)
!
	REAL*8 RJ(ND)
	REAL*8 DTAU(ND)
	REAL*8 Z(ND)
	REAL*8 CHI(ND)
	REAL*8 dCHIdR(ND)
!
	REAL*8 TA(ND)
	REAL*8 TB(ND)
	REAL*8 TC(ND)
	REAL*8 Q(ND)
	REAL*8 GAM(ND)
	REAL*8 GAMH(ND)
	REAL*8 SOB(ND)
	REAL*8 XM(ND)
	REAL*8 SOURCE(ND)
	REAL*8 FEDD(ND)
	REAL*8 H_ON_J(ND)
	REAL*8 RSQHNU(ND)
	REAL*8 dlnJdlnR(ND)
!
	REAL*8 T1,T2
	REAL*8 HBC_J
	REAL*8 HBC_CMF
	REAL*8 NBC_CMF
	REAL*8 INBC
	REAL*8 DBB
	REAL*8 FL
	REAL*8 PI
	REAL*8 HFLUX
!
	INTEGER I,J,K,L
	INTEGER LUWARN,WARNING_LU
	EXTERNAL WARNING_LU
!
	LOGICAL LST_DEPTH_ONLY
	LOGICAL FIRST_FREQ
	LOGICAL NEW_FREQ
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
	COMPUTED=.TRUE.
	LUWARN=WARNING_LU()
	WRITE(LUWARN,'(/,A,3I6)')' Begin COMP_GREY_V4',ND,NC,NP
	WRITE(LUER,'(/,A,3I6)')' Begin COMP_GREY_V4',ND,NC,NP
	PI=4.0D0*ATAN(1.0D0)
	CHI(1:ND)=ROSSMEAN(1:ND)
!
! Compute the Rosseland optical depth scale. J is used to indicate the second
! depth used to compute exponential or power law slope.
!
	J=5
	DO K=5,10
	  IF(CHI(J)/CHI(1) .GT. 5.0D0)EXIT
	  J=K
	END DO
	IF(PLANE_PARALLEL_NO_V)THEN
	  CALL TORSCL_V3(TAU_ROSS,CHI,R,TB,TC,ND,METHOD,'EXP',J,L_FALSE)
	ELSE
	  CALL TORSCL_V3(TAU_ROSS,CHI,R,TB,TC,ND,METHOD,'PCOMP',J,L_FALSE)
	END IF
!
	IF(PLANE_PARALLEL_NO_V)THEN
           FEDD=0.3333D0
	   HBC_CMF=0.7D0; NBC_CMF=0.0D0; FL=1.0D0; INBC=0.1D0
	   NEW_FREQ=.TRUE.; FIRST_FREQ=.TRUE.
!
! Note HFLUX=LUM*Lsun/16/(PI*PI)/10**2 (10**2 for 1/R**2).
! DBB - dBdR = 3. Chi. L/16(piR)**2 and is used for the lower boundary
! diffusion approximation. Since we are dealing with a plane-parallel
! atmopshere, we divide HFLUX by R*^2.
!
	   HFLUX=3.826D+13*LUM/16.0D0/PI**2/R(ND)/R(ND)
           DBB=3.0D0*CHI(ND)*HFLUX
	   T1=1000.0D0
!
! Compute radial (vertical) optical depth increments.
!
	   CALL DERIVCHI(TB,CHI,R,ND,METHOD)
           CALL NORDTAU(DTAU,CHI,R,R,TB,ND)
!
!
! Compute the solution vector. Note that the units need to be
! eventually included. The following follows direcly from d2K/d2Tau=0.
!
	   DO WHILE(T1 .GT. 1.0D-08)
	     T2=FEDD(1)/HBC_CMF
	     RJ(1)=HFLUX/HBC_CMF
	     DO I=2,ND
	       T2=T2+DTAU(I-1)
	       RJ(I)=T2*HFLUX/FEDD(I)
	     END DO
!
	     SOURCE(1:ND)=RJ
	     CALL FCOMP_PP_V2(R,TC,GAMH,SOURCE,CHI,IPLUS,HBC_CMF,
	1               NBC_CMF,INBC,DBB,IC,THK_CONT,DIF,ND,NC,METHOD)
	     T1=0.0D0
	     DO I=1,ND
	       T1=MAX(ABS(FEDD(I)-GAMH(I)),T1)
	       FEDD(I)=GAMH(I)
	     END DO
	     NEW_FREQ=.FALSE.
	   END DO
!
	   OPEN(UNIT=7,FILE='GREY_CHK',STATUS='UNKNOWN')
	     WRITE(7,'(A)')
	     WRITE(7,'(A)')'Check for plane-parallel grey atmosphere calculation'
	     WRITE(7,'(A)')
	     WRITE(7,'(4X,A,5(11X,A))')'I','    R','    J','Ray J','    f','  Tau'
	     DO I=1,ND
	       WRITE(7,'(1X,I4,5(ES16.6))')I,R(I),RJ(I),TC(I),GAMH(I),TAU_ROSS(I)
	     END DO
	   CLOSE(UNIT=7)
!
	ELSE IF(JGREY_WITH_V_TERMS)THEN
!
! This routine will supersede the one above, and included to zero
! order the effect of the velocity field.
!
	   T2=1.0D-05		!Accuracy to converge f
	   CALL JGREY_WITH_FVT(RJ,SOB,CHI,R,V,SIGMA,
	1           P,AQW,HMIDQW,KQW,NMIDQW,
	1           LUM,METHOD,DIF,IC,T2,ND,NC,NP)
!
	ELSE IF(USE_J_REL)THEN
!
	  FEDD(1:ND)=0.3D0		!Initial guess
	  H_ON_J(1:ND)=0.0D0		!Initial guess
	  GAMH(1:ND)=0.0D0		!Old FEDD
	  XM(1:ND)=0.0D0		!As grey solution, not needed (ETA)
	  dlnJdlnR=0.0D0
	  NEW_FREQ=.TRUE.
	  WRITE(LUER,*)'Using MOM_JREL_GREY_V2 for grey solution'
!
! Note
!   HFLUX=LUM*Lsun/16/(PI*PI)/10**2/R**2 (10**2 for 1/R**2).
!   DBB = dBdR = 3.Chi.L/16(piR)**2 
!   DBB is used for the lower boundary diffusion approximation. 
!
	  HFLUX=3.826D+13*LUM/(4.0D0*PI*R(ND))**2
          DBB=3.0D0*HFLUX*CHI(ND)
	  T1=1.0D0
	  DO WHILE(T1 .GT. 1.0D-05)
            CALL MOM_JREL_GREY_V2(XM,CHI,CHI,V,SIGMA,R,
	1              H_ON_J,FEDD,dlnJdlnR,
	1              RJ,RSQHNU,HBC_CMF,INBC,
	1              DIF,DBB,IC,METHOD,
	1              L_TRUE,L_TRUE,NEW_FREQ,COMPUTED,ND)
	    IF(.NOT. COMPUTED)RETURN
!
	    CALL FGREY_NOREL_V1(FEDD,H_ON_J,RJ,CHI,R,V,SIGMA,
	1              P,AQW,HMIDQW,KQW,LUM,IC,METHOD,
	1              HBC_CMF,INBC,DIF,ND,NC,NP)
	    T1=0.0D0
	    DO I=1,ND
	      T1=MAX(ABS(FEDD(I)-GAMH(I)),T1)
	      GAMH(I)=FEDD(I)
	    END DO
	    NEW_FREQ=.FALSE.
	    WRITE(LUER,*)'Current grey iteration accuracy is',T1
	  END DO
!
! Output gray fluxes and boundary Eddington factors for next model in time sequence.
!     T1  is used as a dummy argument for NU.
!      I  is used as a dummy argument for NCF. Neither value is used when the
!  'GREY' option is set. The logical variable L_TRUE is also not accessed.
!
	DO I=1,ND
	  TA(I)=RJ(I)*R(I)*R(I)
	END DO
	T2=HFLUX*R(ND)*R(ND)      !R^2.H at inner boundary.
	CALL OUT_JH(TA,RSQHNU,T2,HBC_CMF,T1,I,R,V,ND,L_TRUE,'GREY')
!
	ELSE IF(USE_DJDT_RTE)THEN
!
! TA is a temporary vector with the change in enthalpy.
!
	   WRITE(LUER,*)'Calling JGREY_HUB_DDT_V4'
	   WRITE(LUWARN,*)'Calling JGREY_HUB_DDT_V4'
	   T2=1.0D-06		!Accuracy to converge f
	   CALL JGREY_HUB_DDT_V4(RJ,SOB,CHI,PLANCKMEAN,R,V,SIGMA,POPS,
	1              P,AQW,HMIDQW,KQW,LUM,METHOD,DIF,IC,
	1              T2,INCL_DJDT_TERMS,TIME_SEQ_NO,COMPUTED,ND,NC,NP,NT)
	   IF(.NOT. COMPUTED)RETURN
!
	ELSE
!
! Will use FEDD for F, GAM for NEWRJ, GAMH for NEWRK, and T2 for NEWHBC.
! Will use HBC_J for HBC. No need to modify JGREY, as outer boundary
! will always be optically thin.
!
	    DO I=1,ND
	      FEDD(I)=1.0D0/3.0D0
	    END DO
	    HBC_J=1.0D0
	    T1=1000.0D0
	    DO WHILE(T1 .GT. 1.0D-06)
	      CALL JGREY(TA,TB,TC,XM,DTAU,R,Z,P,RJ,
	1        GAM,GAMH,Q,FEDD,CHI,dCHIdR,
	1        AQW,KQW,LUM,HBC_J,T2,NC,ND,NP,METHOD)
	      T1=0.0D0
	      DO I=1,ND
	        T1=MAX(ABS(FEDD(I)-GAMH(I)),T1)
	        FEDD(I)=GAMH(I)
	      END DO
	      T1=MAX(ABS(HBC_J-T2),T1)
	      HBC_J=T2
!	      WRITE(LUER,'('' Maximum change in Grey F is '',1P,E11.4)')T1
	    END DO
	END IF
	WRITE(LUWARN,'(A)')' '
!
! Compute the temperature distribution, and the Rosseland optical depth scale.
! NB sigma=5.67D-05 and the factor of 1.0D-04 is to convert T from units of 
! K to units of 10^4 K. The ' ' in TORSCL indicates TYPE of atmosphere,
! and here is set to ' ' so as TORSCL assumes a 1/r^2 density dependence
! at boundary.
! 
	DO I=1,ND
	  TGREY(I)=((PI/5.67D-05*RJ(I))**0.25D0)*1.0D-04
	END DO
!
	RETURN
	END
