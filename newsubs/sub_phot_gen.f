!
! Subroutine to compute the photoionization cross-section for XzV. The
! cross-sections have been scaled by dex(10) to keep chi*z constant.
!
! The cross-sections are stored in the data module PHOT_DATA_MOD. They can have
! a functional form, or be a tabulation of cross-section versus frequency
! (normalized by threshold value).
!
!***************************
!***************************
! Currently implemented fits.
!
!     1 - Seaton formula fit [sigma_o,alpha,beta]
!     2 - Hydrogenic: Split l (z states, n > 11)
!     3 - Hydrogenic: Pure n level (all l, n >= 13)
!     4 - Used for CIV rates from Leobowitz (JQSRT 1972,12,299)
!     5 - Opacity project fits (from Peach, Sraph, and Seaton (1988)
!     6 - Hummer fits to the opacity cross-sections for HeI
!     7 - Modifed Seaton fit --- cross-section zero until offset edge.
!     8 - Modifed Hydrogenic split l: cross-section zero until offset edge.
!     9 - Ground-state photoioinization cross-sections from Werner et al.
!    20 - Opacity Project: smoothed [number of data pairs]
!    21 - Opacity Project: scaled, smoothed [number of data pairs]
!    30 - Call XCROSS_V2 (For Lithium-like ions).
!    31 - Call XCROSS_V2 (For inner shell cross-sections: n,l specified)
!
!******************************************************************************
!******************************************************************************
!
	SUBROUTINE SUB_PHOT_GEN(ID,PHOT,FREQ,GS_EDGE,NLEVS,PHOT_ID,SET_TO_EDGE)
	USE SET_KIND_MODULE
!
! IF SET_TO_EDGE is TRUE, the routine will return the
! bound-free cross-section at threshold for those levels which have
! NU < EDGE. If SET_TO_EDGE if false, ZERO will be returned for those
! levels with NU < EDGE. To return all of the THRESHOLD cross-sections,
! set NU=0, SET_TO_EDGE=.TRUE.
!
! If FREQ=0, and SET_EDGE_FREQ is true the EDGE frequecnies can be obtained
!   by passing -PHOT_ID as PHOT_ID.
! If FREQ=0, and SET_EDGE_FREQ is true the EDGE purge edge cross-sections
!   can be obtained by passing PHOT_ID+100 as PHOT_ID.
!
! Contains cross-section data. Also include in RD_PHOT_XzV and
!
	USE PHOT_DATA_MOD
	USE HYD_BF_PHOT_DATA
	USE XRAY_DATA_MOD
!
	IMPLICIT NONE
!
! Altered 25-Oct-2023 : Added error messages concerning MAX_N_PQN and MAX_L_PQN (20-Aug-2023)
! Altered 15-Nov-2021 : Fixed call to XCROSS_V2 for type 31. Added STOP type 30.
! Altered 09-Oct-2018 : Now set the b-f gaunt factor to unity for n>30  (previously crashed).
! Altered 07-Oct-2015 : Bug fix for Type 7 (modified Seaton formula).
!                         Offset was beeing added to the current frequency instead
!                            of the ionization edge.
! Altered 21-Apr-2011 : Bug fix for Types 30, 31 (inadvertantly applying to all levels)
! Altered 11-Mar-2011 : 2s,2p opacity added for SiIV, PV
!                       3s,3p opacity added for states with one 3d electron.
! Altered 09-Nov-2010 : Modified fit type 9.
! Altered 17-Apr-2008 : Installed fit Type 8:
! Altered 30-Apr-2004 : Modified Seaton fit (Type 7) installed.
! Altered 07-Dec-2001 : EQUAL_CORFAC installed for AMD processors and PGI compiler.
!                       Bug fix. Correct edge cross-section was not being returned for
!                           when PHOT>100.
! Altered 08-Dec-2000 : XCROSS_V2 installed.
! Altered 19-DEC-1997 : Altered to allow pure edge photoionization cross-
!                         sections B, C etc. to be computed. The actual
!                         B and C edge frequencies can also be returned.
!                         A trick is used which avoids altering routines that
!                         call SUB_PHOT_XzV.
!                       Bug fix for option '5' for NU outside fit range.
! Altered 16-Dec-1996 : Several bug fixes for handling new hydrogenic-
!                         cross section computation.
!                       Edge cross-section returned if nu < EDGE and
!                         SET_TO_EDGE for all PHOT_ID values.
!                         (needed by SET_EDGE_FREQ_V2)
! Altered 14-Nov-1996 : Bug fixed for Cross-section TYPE 20,21 when NU
!                         is greater than tabulated value.
! Created 18-Sep-1996 : Based on PHOT_XzV
!
	INTEGER ID
	INTEGER NLEVS			!Number of levels
	INTEGER PHOT_ID		!Which photoionization route.
	REAL(KIND=LDP) PHOT(NLEVS)		!Cross-section
	REAL(KIND=LDP) GS_EDGE(NLEVS)		!Energy for ionization to Ground State!
	REAL(KIND=LDP) FREQ
	LOGICAL SET_TO_EDGE
!
! External functions.
!
	INTEGER ERROR_LU
	REAL(KIND=LDP) VOIGT,HYDCROSSL,XCROSS_V2,GBF
	EXTERNAL VOIGT,HYDCROSSL,XCROSS_V2,GBF,ERROR_LU
!
! Common block with opacity/emissivity constants.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
! Local vectors
!
	REAL(KIND=LDP) FREQ_VEC(NLEVS)
	REAL(KIND=LDP) INDX(NLEVS)
!
! Local variables.
!
	INTEGER I,J,K,l,N,TERM
	INTEGER IZ,INE
	INTEGER LMIN,LMAX,LMID
	INTEGER INC,INC_SAV
	INTEGER GRID_ID
	INTEGER LST,LEND
!
	REAL(KIND=LDP) U			!Defined as FREQ/EDGE
	REAL(KIND=LDP) RU			!Defined as EDGE/FREQ
!
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) DOP_NU
	REAL(KIND=LDP) EDGE			!Ionization energy
	REAL(KIND=LDP) A_VOIGT
	REAL(KIND=LDP) V_VOIGT
	REAL(KIND=LDP) X
	REAL(KIND=LDP) RJ
	REAL(KIND=LDP) SUM
!
	REAL(KIND=LDP), PARAMETER :: EV_TO_HZ=0.241798840766_LDP
	REAL(KIND=LDP), PARAMETER :: EQUAL_COR_FAC=1.0_LDP+1.0E-14_LDP
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ITHREE=3
	LOGICAL ALL_DONE
	LOGICAL DO_PURE_EDGE
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! If PHOT_ID is < 0, FREQ=0, and SET_TO_EDGE is true, we return the
! ionization edges in PHOT. Works for all photoionization routes. NO
! photoionization cross-sections are computed. For use with SET_TO_EDGE.
!
	IF(PHOT_ID .LT. 0)THEN
	  IF(SET_TO_EDGE .AND. FREQ .EQ. 0.0_LDP)THEN
	    PHOT(1:NLEVS)=GS_EDGE(1:NLEVS)+PD(ID)%EXC_FREQ(ABS(PHOT_ID))
	  ELSE
	    WRITE(ERROR_LU(),*)'Invalid option list in SUB_PHOT_XzV[1]'
	    STOP
	  END IF
	  RETURN
	END IF
!
! IF PHOT_ID > 100, we adjust DO_PURE_EDGE so that a pure B level-cross
! section is returned, independent of whether the B levels is present.
! For use with SET_TO_EDGE.
!
	IF(PHOT_ID .GT. 100)THEN
	  IF(SET_TO_EDGE .AND. FREQ .EQ. 0)THEN
	    DO_PURE_EDGE=.TRUE.
	    PHOT_ID=PHOT_ID-100
	  ELSE
	    WRITE(ERROR_LU(),*)'Invalid option list in SUB_PHOT_XzV[2]'
	    STOP
          END IF
	ELSE
	  DO_PURE_EDGE=.FALSE.
	END IF	
!
	PHOT(1:NLEVS)=0.0_LDP
	FREQ_VEC(1:NLEVS)=FREQ
	DO I=1,NLEVS
	  EDGE=GS_EDGE(I)+PD(ID)%EXC_FREQ(PHOT_ID)
	  IF(FREQ .LT. EDGE .AND. SET_TO_EDGE)FREQ_VEC(I)=EDGE
	END DO
!
! Check if cross-sections were just computed. If so, we use those value.
!
	ALL_DONE=.TRUE.
	DO I=1,NLEVS
	  IF(FREQ_VEC(I) .NE. PD(ID)%LST_FREQ(I,PHOT_ID))ALL_DONE=.FALSE.
	END DO
	IF(ALL_DONE)THEN
	  PHOT(1:NLEVS)=PD(ID)%LST_CROSS(1:NLEVS,PHOT_ID)
	  RETURN
	END IF
!
! 
!
! Dam --- we have to do some work!
!
!
! NB: If one of DO_PHOT(PHOT_ID,..) is true, we return a cross-section.
! Usually if PHOT_ID=1 we always return a photoionization cross section.
! This rate is modified to include the other photoionization cross-sections
! when their destination level in the ion is unavailable.
!
	DO K=1,PD(ID)%NUM_PHOT_ROUTES
	  GRID_ID=PD(ID)%PHOT_GRID(K)
	  IF( PD(ID)%DO_PHOT(PHOT_ID,K) .OR.
	1            (DO_PURE_EDGE .AND. K .EQ. PHOT_ID) )THEN
C
C First do those parts not readily vectorized.
C
	    DO I=1,NLEVS
	      EDGE = (GS_EDGE(I) + PD(ID)%EXC_FREQ(K))/EQUAL_COR_FAC
	      IF(FREQ_VEC(I) .GE. EDGE)THEN
	        TERM=PD(ID)%A_ID(I,GRID_ID)
!
! Determine cross-section formulation, and do appropriate fitting procedure.
! We check whether tabulated format first, as this will become the usual
! fitting procedure.
!
! NB: ST_LOC and END_LOC refer to the storage, and hence are referenced
! according to the term. LST_LOC refers to some location between
! ST_LOC and END_LOC and are referenced according to the particular
! level under consideration.
!
	        U=FREQ_VEC(I)/EDGE
	        LMIN=PD(ID)%ST_LOC(TERM,GRID_ID)
	        IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 20 .OR.
	1              PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 21)THEN
!
! Check if in range of tabulation. If not, we extrapolate according to
! 1/nu^3.
!
	          N=PD(ID)%END_LOC(TERM,GRID_ID)
	          IF(U .LT. PD(ID)%NU_NORM(N))THEN
!
! Most cross-sections are computed sequentially. Therefore we start from
! that point, and search from there.
!
! The first step is to bracket the location of U in the table.
!
	            INC=1
	            INC_SAV=0	  		!Some value other than inc
	            LMIN=PD(ID)%LST_LOC(I,GRID_ID)
	            LMAX=LMIN+1
	            DO WHILE(INC .NE. INC_SAV)
	              INC_SAV=INC
	              IF(U .LT. PD(ID)%NU_NORM(LMIN))THEN
	                LMIN=LMIN-INC
	                LMIN=MAX(LMIN,PD(ID)%ST_LOC(TERM,GRID_ID))
	                INC=INC+INC
	              ELSE IF(U .GT. PD(ID)%NU_NORM(LMAX))THEN
	                LMAX=LMAX+INC
	                LMAX=MIN(LMAX,PD(ID)%END_LOC(TERM,GRID_ID))
	                INC=INC+INC
	              END IF
	            END DO
!
! We have located the region. We now narrow in on the required region.
!
	            DO WHILE(LMAX-LMIN .GT. 1)
	              LMID=(LMAX+LMIN)/2
	              IF(U .GT. PD(ID)%NU_NORM(LMID))THEN
	                LMIN=LMID
	              ELSE
	                LMAX=LMID
	              END IF
	            END DO
	            INDX(I)=LMIN
	          END IF
!
! Hydrogenic transitions.
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 2)THEN
!
! N, LST and LEND must be integer - all other values are double precision.
! If it wasn't for the loop over l, this could be done in the main loop.
!
! BF_CROSS contains the LOG10 hydrogenic cross-section.
!
	          N=NINT( PD(ID)%CROSS_A(LMIN) )
	          LST=NINT( PD(ID)%CROSS_A(LMIN+1) )
	          LEND=NINT( PD(ID)%CROSS_A(LMIN+2) )
	          IF(LEND+1 .GT. N)THEN
	            WRITE(6,*)'Error in SUB_PHOT_GEN -- invalid PHOT_PARAMS for TYPE 2'
	            WRITE(6,*)'ID=',ID
	            WRITE(6,*)'N,LST,LEND=',N,LST,LEND
	            STOP
	          END IF
!
	          X=LOG10(U)
	          RJ=X/L_DEL_U
	          J=RJ
	          T1=RJ-J
	          J=J+1
	          SUM=0.0_LDP
	          DO L=LST,LEND
	            IF(J .LT. N_PER_L)THEN
	              J=J+BF_L_INDX(N,L)-1
	              T2=T1*BF_L_CROSS(J+1)+(1.0_LDP-T1)*BF_L_CROSS(J)
	            ELSE
	              J=BF_L_INDX(N,L)+N_PER_L-1
	              T2=(BF_L_CROSS(J)-BF_L_CROSS(J-1))*
	1                   (RJ-N_PER_L)+BF_L_CROSS(J)
	            END IF
	            SUM=SUM+(2*L+1)*(10.0_LDP**T2)
		    J=RJ+1 		!Restore as corrupted	
	          END DO
	          SUM=SUM*( PD(ID)%NEF(I,GRID_ID)/(N*PD(ID)%ZION) )**2
	          PHOT(I)=PHOT(I) + SUM/( (LEND-LST+1)*(LEND+LST+1) )
!
! Hydrogenic transitions.
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 8)THEN
!
! N, LST and LEND must be integer - all other values are double precision.
! If it wasn't for the loop over l, this could be done in the main loop.
!
! BF_CROSS contains the LOG10 hydrogenic cross-section.
!
	          IF(FREQ_VEC(I) .GE. EDGE+PD(ID)%CROSS_A(LMIN+3))THEN
	            U=FREQ_VEC(I)/(EDGE+PD(ID)%CROSS_A(LMIN+3))
	            N=NINT( PD(ID)%CROSS_A(LMIN) )
	            LST=NINT( PD(ID)%CROSS_A(LMIN+1) )
	            LEND=NINT( PD(ID)%CROSS_A(LMIN+2) )
!
	            IF(N .GT. MAX_N_PQN)THEN
	              WRITE(6,*)'Error in SUB_PHOT_GEN -- invalid N value for cross-section type 8'
	              WRITE(6,*)'Maximum N value is ',MAX_N_PQN
	              WRITE(6,*)'N,LST,LEND=',N,LST,LEND
	              WRITE(6,*)'ID=',ID,'NLEVS=',NLEVS
	              WRITE(6,*)'PHOT_ID=',PHOT_ID,'SET_TO_EDGE=',SET_TO_EDGE
	              WRITE(6,*)'TERM=',TERM,K,PD(ID)%CROSS_TYPE(TERM,GRID_ID)
	              FLUSH(UNIT=6)
	              STOP
	            END IF
!
	            IF(LST .GT. MAX_L_PQN .OR. LEND .GT. MAX_L_PQN)THEN
	              WRITE(6,*)'Error in SUB_PHOT_GEN -- invalid L value for cross-section type 8'
	              WRITE(6,*)'Maximum L value is ',MAX_L_PQN
	              WRITE(6,*)'N,LST,LEND=',N,LST,LEND
	              WRITE(6,*)'ID=',ID,'NLEVS=',NLEVS
	              WRITE(6,*)'PHOT_ID=',PHOT_ID,'SET_TO_EDGE=',SET_TO_EDGE
	              WRITE(6,*)'TERM=',TERM,K,PD(ID)%CROSS_TYPE(TERM,GRID_ID)
	              FLUSH(UNIT=6)
	              STOP
	            END IF
!
	            X=LOG10(U)
	            RJ=X/L_DEL_U
	            J=RJ
	            T1=RJ-J
	            J=J+1
	            SUM=0.0_LDP
	            DO L=LST,LEND
	              IF(J .LT. N_PER_L)THEN
	                J=J+BF_L_INDX(N,L)-1
	                T2=T1*BF_L_CROSS(J+1)+(1.0_LDP-T1)*BF_L_CROSS(J)
	              ELSE
	                J=BF_L_INDX(N,L)+N_PER_L-1
	                T2=(BF_L_CROSS(J)-BF_L_CROSS(J-1))*
	1                   (RJ-N_PER_L)+BF_L_CROSS(J)
	              END IF
	              SUM=SUM+(2*L+1)*(10.0_LDP**T2)
		      J=RJ+1 		!Restore as corrupted	
	            END DO
	            SUM=SUM/PD(ID)%ZION/PD(ID)%ZION	!Ignore neff correction :( PD(ID)%NEF(I,GRID_ID)/(N*PD(ID)%ZION) )**2
	            PHOT(I)=PHOT(I) + SUM/( (LEND-LST+1)*(LEND+LST+1) )
	          END IF	!Above threshold
!
	        END IF		!Type of cross-section
	      END IF		!Non zero
	    END DO		!Loop over level.
!
!
! Do the loop over N for the code that should be vectorizable.
! The Function GBF should be inlined.
!
	    DO I=1,NLEVS
	      EDGE = (GS_EDGE(I) + PD(ID)%EXC_FREQ(K))/EQUAL_COR_FAC
	      TERM=PD(ID)%A_ID(I,GRID_ID)
	      LMIN=PD(ID)%ST_LOC(TERM,GRID_ID)
	      IF(FREQ_VEC(I) .GE. EDGE)THEN
!
! Cross-sections are tabulated in terms of v/v_o
!
	        U=FREQ_VEC(I)/EDGE
!
! Determine cross-section formulation, and do appropriate fitting procedure.
! We check whether tabulated format first, as this will become the usual
! fitting procedure.
!
	        IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 20 .OR.
	1                   PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 21)THEN
!
! Check if in range of tabulation. If not, we extrapolate according to
! 1/nu^3.
!
	          N=PD(ID)%END_LOC(TERM,GRID_ID)
	          IF(U .GE. PD(ID)%NU_NORM(N))THEN
	            PHOT(I)=PHOT(I)+CONV_FAC*PD(ID)%CROSS_A(N)*(PD(ID)%NU_NORM(N)/U)**3
	            PD(ID)%LST_LOC(I,GRID_ID)=N-1
	          ELSE
!
! Do linear interpolation.
!
	            LMIN=INDX(I)
	            LMAX=LMIN+1
	            T1=(PD(ID)%NU_NORM(LMAX)-U)/(PD(ID)%NU_NORM(LMAX)-PD(ID)%NU_NORM(LMIN))
	            PHOT(I)=PHOT(I)+ CONV_FAC*(
	1              (1.0_LDP-T1)*PD(ID)%CROSS_A(LMAX)+T1*PD(ID)%CROSS_A(LMIN)  )
	            PD(ID)%LST_LOC(I,GRID_ID)=LMIN
	          END IF		!Outside table range?
!
!
!
! Seaton fit.
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 1 .AND.
	1                             PD(ID)%CROSS_A(LMIN) .NE. 0)THEN
	          RU=EDGE/FREQ_VEC(I)
	          PHOT(I)=PHOT(I) + CONV_FAC*
	1            PD(ID)%CROSS_A(LMIN)*( PD(ID)%CROSS_A(LMIN+1) +
	1             (1.0_LDP-PD(ID)%CROSS_A(LMIN+1))*RU )*( RU**PD(ID)%CROSS_A(LMIN+2) )
!
!
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 3)THEN
!
! HYD_N_DATA contains the Bound-free gaunt factor.
!
	          N=PD(ID)%CROSS_A(LMIN+1)
	          IF(N .GT. 30)THEN
	            T1=1.0_LDP
	          ELSE
	            X=LOG10(U)
	            RJ=X/N_DEL_U
	            J=RJ
	            T1=RJ-J
	            J=J+1
	            IF(J .LT. N_PER_N)THEN
	              J=J+BF_N_INDX(N)-1
	              T1=T1*BF_N_GAUNT(J+1)+(1.0_LDP-T1)*BF_N_GAUNT(J)
	            ELSE
!
! Power law extrapolation.
!
	              J=BF_N_INDX(N)+N_PER_N-1
	              T1=LOG10(BF_N_GAUNT(J-1)/BF_N_GAUNT(J))
1                     T1=BF_N_GAUNT(J)*( 10.0_LDP**(T1*(N_PER_N-RJ)) )
	            END IF
	          END IF
!
! NB: ZION is already include in ALPHA_BF
!
	          PHOT(I)=PHOT(I) + PD(ID)%ALPHA_BF*T1*PD(ID)%CROSS_A(LMIN)/
	1                PD(ID)%NEF(I,GRID_ID)/N/
	1                ( (FREQ_VEC(I)*PD(ID)%NEF(I,GRID_ID))**3 )
!
!
! Used for CIV recombination rates for s and p states.
!(ref -Leibowitz J.Q.S.R.T 1972,12,299)
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 4)THEN
	          RU=EDGE/FREQ_VEC(I)
	          T1=CONV_FAC*(  PD(ID)%CROSS_A(LMIN)+RU*( PD(ID)%CROSS_A(LMIN+1) +
	1                   RU*(PD(ID)%CROSS_A(LMIN+2) + RU*(PD(ID)%CROSS_A(LMIN+3)+
	1                   RU*(PD(ID)%CROSS_A(LMIN+4)+RU*PD(ID)%CROSS_A(LMIN+5)))) )  )
	          IF(T1 .GT. 0.0_LDP)PHOT(I)=PHOT(I)+T1
!
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 5)THEN
!
! These fits are fits to the opacity cross section. Data taken from opacity
! project - Peach, Saraph, and Seaton (1988, C J. Phys. B: 21, 3669-3683)
!
! If the cross-section lies outside the range given by the fits, we assume
! that it scales as nu^{-2}. The values in this region should be unimportant.
!
	          X=MIN( U,PD(ID)%CROSS_A(LMIN+4) )
	          X=LOG10(X)
	          T1=10**(  PD(ID)%CROSS_A(LMIN)+X*( PD(ID)%CROSS_A(LMIN+1) +
	1                 X*(PD(ID)%CROSS_A(LMIN+2) +
	1                 X*PD(ID)%CROSS_A(LMIN+3)) ) + LG10_CONV_FAC  )
	          IF(U .GT. PD(ID)%CROSS_A(LMIN+4))T1=T1*(PD(ID)%CROSS_A(LMIN+4)/U)**2
	          PHOT(I)=PHOT(I)+T1
!
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 6)THEN
!
! This type is for the Hummer fits to the Opacity cross sections of HeI.
! See HEI_PHOT_OPAC.
!
! We issue the SPACING command to ensure that rounding error does not
! cause X to be < 0 at the bound-free edge.
!
! U+SPACING(U) is the smallest number different from U (and lager).
!
	          X=LOG10(U+3.0_LDP*SPACING(U))
	          IF(X .GE. 0.0_LDP)THEN
                    IF(X .LT. PD(ID)%CROSS_A(LMIN+4))THEN
	              T1=((PD(ID)%CROSS_A(LMIN+3)*X+PD(ID)%CROSS_A(LMIN+2))*X +
	1                  PD(ID)%CROSS_A(LMIN+1))*X+PD(ID)%CROSS_A(LMIN)
	            ELSE
                      T1=PD(ID)%CROSS_A(LMIN+5)+PD(ID)%CROSS_A(LMIN+6)*X
	            END IF
	            PHOT(I)=PHOT(I)+10.0_LDP**(T1+LG10_CONV_FAC)
	          END IF
!
! Modified seaton fit.
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 7 .AND.
	1                             PD(ID)%CROSS_A(LMIN) .NE. 0)THEN
!	          RU=EDGE/(FREQ_VEC(I)+PD(ID)%CROSS_A(LMIN+3))
	          RU=(EDGE+PD(ID)%CROSS_A(LMIN+3))/FREQ_VEC(I)
	          IF(RU .LE. 1.0_LDP)THEN
	            PHOT(I)=PHOT(I) + CONV_FAC*
	1              PD(ID)%CROSS_A(LMIN)*( PD(ID)%CROSS_A(LMIN+1) +
	1               (1.0_LDP-PD(ID)%CROSS_A(LMIN+1))*RU )*( RU**PD(ID)%CROSS_A(LMIN+2) )
	          END IF
!
! We assume all ionizations are to the ground state. Previously we read in
! an additional offset, but this has now been superceded.
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 9)THEN
	          LMIN=LMIN-8
	          DO J=1,(PD(ID)%END_LOC(TERM,GRID_ID)-PD(ID)%ST_LOC(TERM,GRID_ID)+1)/8
	            LMIN=LMIN+8
	            IF(J .NE. 1)EDGE=0.241798840766_LDP*PD(ID)%CROSS_A(LMIN+2)
	            IF(FREQ_VEC(I) .GE. EDGE)THEN
                      U=FREQ_VEC(I)/PD(ID)%CROSS_A(LMIN+3)/EV_TO_HZ
                      T1=(U-1.0_LDP)**2 + PD(ID)%CROSS_A(LMIN+7)**2
                      T2=U**( 5.5_LDP+PD(ID)%CROSS_A(LMIN+1)-0.5_LDP*PD(ID)%CROSS_A(LMIN+6) )
                      T3=( 1.0_LDP+SQRT(U/PD(ID)%CROSS_A(LMIN+5)) )**PD(ID)%CROSS_A(LMIN+6)
                      PHOT(I)=PHOT(I)+1.0E-08_LDP*T1*PD(ID)%CROSS_A(LMIN+4)/T2/T3
                    END IF
	          END DO
!
!
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 30)THEN
	          T1=PD(ID)%AT_NO+1.0_LDP-PD(ID)%ZION	!# of elec. in species.
	          T1=XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,IZERO,IZERO,L_FALSE,L_TRUE)
	          IF(T1 .GT. 0)PHOT(I)=PHOT(I)+T1
	          WRITE(6,*)'Error in SUB_PHOT_GEN -- calling CROSS TYPE 30 - 31'
	          WRITE(6,*)'Loop   needs fixing'
	          WRITE(6,*)GS_EDGE,NLEVS,PHOT_ID
	          STOP
!
	        ELSE IF(PD(ID)%CROSS_TYPE(TERM,GRID_ID) .EQ. 31)THEN
	          T1=PD(ID)%AT_NO+1.0_LDP-PD(ID)%ZION	!# of elec. in species.
	          N=NINT(PD(ID)%CROSS_A(LMIN))
	          L=NINT(PD(ID)%CROSS_A(LMIN))
	          T1=XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,N,L,L_FALSE,L_FALSE)
	          IF(T1 .GT. 0)PHOT(I)=PHOT(I)+T1
	          WRITE(6,*)'Error in SUB_PHOT_GEN -- calling CROSS TYPE 30 - 31'
	          WRITE(6,*)'Loop   needs fixing'
	          WRITE(6,*)GS_EDGE,NLEVS,PHOT_ID
	          STOP
!
! More cross-section types can be added in here.
!
	      	END IF		!Type of cross-section
	      END IF		!FREQ .GE. EDGE
	    END DO		!Loop over level
	  END IF		!Do this final state for this PHOT_ID
	END DO			!Loop over PHOT_ID
!
! Include dielectronic! lines if present. The factor of 10^{-15} in the
! expression for the opacity arises in the conversion of DOP_NU to Hz.
!
	IF(PHOT_ID .EQ. 1 .AND. PD(ID)%NUM_DIE .NE. 0)THEN
	  DO I=1,NLEVS
	    DO J=PD(ID)%ST_INDEX(I),PD(ID)%END_INDEX(I)
	      IF( FREQ_VEC(I) .GE. PD(ID)%NU_MIN(J) .AND.
	1                           FREQ_VEC(I) .LE. PD(ID)%NU_MAX(J) )THEN
	        DOP_NU=PD(ID)%NU_ZERO(J)*PD(ID)%VSM_KMS/2.998E+05_LDP
	        A_VOIGT=PD(ID)%GAMMA(J)/DOP_NU
	        V_VOIGT=(FREQ-PD(ID)%NU_ZERO(J))/DOP_NU
	        PHOT(I)=PHOT(I)+
	1             1.0E-15_LDP*OPLIN*PD(ID)%OSC(J)*VOIGT(A_VOIGT,V_VOIGT)/DOP_NU
	      END IF
	    END DO
	  END DO
	END IF
!
! Include X-rays opacity with regular cross-sections if requested. This is
! usually done when only the ground state of the next ioization stage is
! included.
!
	IF(PD(ID)%DO_KSHELL_W_GS .AND. PHOT_ID .EQ. 1 .AND. FREQ .GT. 0)THEN
	  T1=PD(ID)%AT_NO+1.0_LDP-PD(ID)%ZION			!Number of electrons in species.
	  T2=XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,IZERO,IZERO,L_FALSE,L_TRUE)
	  IF(NINT(T1) .EQ. 3)T2=T2+XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,IONE,IZERO,L_FALSE,L_TRUE)
	  IF(T2 .GT. 0)THEN
	    PHOT(1:NLEVS)=PHOT(1:NLEVS)+T2
	  END IF
	END IF
!
! Add in L-shell X-ray absorption for species such as SiIV, PV. We assume that only
! one electron is ejected.
!
	T1=PD(ID)%AT_NO+1.0_LDP-PD(ID)%ZION			!Number of electrons in species.
	IZ=NINT(PD(ID)%AT_NO)
        INE=NINT(T1)
	IF(INE .EQ. 11 .AND. PHOT_ID .EQ. 1 .AND. FREQ .GT. 0)THEN
	  T2=XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,ITWO,IZERO,L_FALSE,L_TRUE)
	  T2=T2+XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,ITWO,IONE,L_FALSE,L_TRUE)
	  IF(T2 .GT. 0.0_LDP)THEN
	    PHOT(1:NLEVS)=PHOT(1:NLEVS)+T2
	  END IF
	END IF
!
! Add in 3s, 3p X-ray absorption for species in which only 1 electron is ejected.
!
	T1=PD(ID)%AT_NO+1.0_LDP-PD(ID)%ZION			!Number of electrons in species.
	IZ=NINT(PD(ID)%AT_NO)
        INE=NINT(T1)
	IF(INE .GE. 19 .AND. PHOT_ID .EQ. 1 .AND. FREQ .GT. 0)THEN
	  T2=0.0_LDP
	  IF(N_ED_EJ(IZ,INE,ITHREE,IZERO) .EQ. 1)THEN
	    T2=XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,ITHREE,IZERO,L_FALSE,L_TRUE)
	  END IF
	  IF(N_ED_EJ(IZ,INE,ITHREE,IONE) .EQ. 1)THEN
	    T2=T2+XCROSS_V2(FREQ,PD(ID)%AT_NO,T1,ITHREE,IONE,L_FALSE,L_TRUE)
	  END IF
	  IF(T2 .GT. 0.0_LDP)THEN
	    PHOT(1:NLEVS)=PHOT(1:NLEVS)+T2
	  END IF
	END IF
!
	PHOT(1:NLEVS)=PHOT(1:NLEVS)*PD(ID)%SCALE_FAC(PHOT_ID)
!
! Store values in case they required again. This will speed up the computation.
!
	PD(ID)%LST_FREQ(1:NLEVS,PHOT_ID)=FREQ_VEC(1:NLEVS)
	PD(ID)%LST_CROSS(1:NLEVS,PHOT_ID)=PHOT(1:NLEVS)
!
	RETURN
	END
