!
! Subroutine to compute the collisional excitation and ionization cross
! sections for an arbitrary species taking into account super levels.
!
! The collison rates among the levels in the full atom are first computed.
! These are then used to compute the rates amongst the super levels.
!
! The collison strengths (OMEGA) must be supplied by OMEGA_COL ---
! this subroutine name is passed in the call. OMEGA_COL uses data in
! the file COL_FILE. SUB_PHOT is also passed to OMEGA_COL.
!
! Single routine as technique is the same for ALL species.
!
	SUBROUTINE SUBCOL_MULTI_V6(
	1                 OMEGA_F,dln_OMEGA_dlnT,COL_S,dCOL_S,
	1                 HN_S,HNST_S,dlnHNST_S_dlnT,N_S,
	1                 HN_F,HNST_F_ON_S,W_F,EDGE_F,AHYD_F,GHYD_F,
	1                 LEVNAME_F,N_F,
	1                 ZION,ID,COL_FILE,OMEGA_COL,
	1                 F_TO_S_MAPPING,COOL,T,ED,ND,
	1                 COMPUTE_BA,FIXED_T,LAST_ITERATION)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 23-Oct-2016 - Bug fix with level dissolution.
!                          Fixed bug introduced in creating V6.
! Altered 17-Oct-2016 - Changed to V6. Updated to improve speed - based on "fast"
!                          version created earlier on OSPREY.
! Altered 05-Apr-2011 - Changed to V5.
!                       HNST_F_ON_S (rather than HNST_F) is passed in call.
!                       HNST_F/HNST_S replaced by HNST_F_ON_S - done to faciliate
!                         modifications allowing lower temperaturs.
!                       Most of editing done early 2011
!
! Altered 16-Jun-1996 : Call to ZERO removed. COL and dCOL now initializd
!                         outside depth loop.
!                       Bug FIX: COL_FILE was declared REAL(KIND=LDP), NOW declared
!                           as CHARACTER. No efffect on VAX, important on CRAY.
! Altered 27-May-1996 : Generic calls used for EXP, SQRT.
! Altered 03-Jan-1995 - HN_F  inserted in call (_V2 changed to _V3)
!                       Adjustment made to allow for variable departure
!                          coeficients in a given super level (collisional
!                          excitation nad dexciation only).
!                       Collisional ionization rates still assume constant
!                          departure coeficients among levels in a given super
!                          levl.
!
! Altered 24-Nov-1995 - Bug fixed: Incorrect dimension for HN_S
! Altered 10-Nov-1995 - Bug fixed in DCOL.
!                       HN_S inserted in call, and HN_F deleted.
! Created 06-Sep-1995 - Based on SUBCOL_HYD version.
!                       Created so that all species with super-levels use
!                         the same collisional routine. Only the routine to
!                         compute OMEGA in the FULL atom is distinct.
!                       Vector COOL was installed to enable checking of
!                         collisonal cooling rates.
!
	INTEGER N_S,N_F,ND
	REAL(KIND=LDP) OMEGA_F(N_F,N_F),dln_OMEGA_dlnT(N_F,N_F)
	REAL(KIND=LDP) COL_F(N_F,N_F),DCOL_F(N_F,N_F)
	REAL(KIND=LDP) COL_S(N_S,N_S,ND),DCOL_S(N_S,N_S,ND)
!
	REAL(KIND=LDP) HN_S(N_S,ND)		!Population of atom with super levels.
	REAL(KIND=LDP) HNST_S(N_S,ND)		!LTE pop. of atom with super levels.
	REAL(KIND=LDP) dlnHNST_S_dlnT(N_S,ND)
!
	REAL(KIND=LDP) HN_F(N_F,ND)		!Population of FULL atom
	REAL(KIND=LDP) HNST_F(N_F,ND)		!LTE population of FULL atom
	REAL(KIND=LDP) HNST_F_ON_S(N_F,ND)	!LTE population of FULL atom
	REAL(KIND=LDP) W_F(N_F,ND)		!Occupation probability
	REAL(KIND=LDP) AHYD_F(N_F,N_F)		!Einstein A coefficient
	REAL(KIND=LDP) EDGE_F(N_F)		!Ionization frequency (10^15 Hz)
	REAL(KIND=LDP) GHYD_F(N_F)		!Statistical weight.
	REAL(KIND=LDP) ZION			!Charge on ion (i.e. 1 for H)
	CHARACTER*(*) COL_FILE		!Name of file with collisonal data.
	CHARACTER*(*) LEVNAME_F(N_F)	!Level names in FULL ATOM.
	INTEGER ID			!Specifies ident. for photiozation data
!
	INTEGER F_TO_S_MAPPING(N_F)
!
	REAL(KIND=LDP) T(ND)			!Temperature (10^4 K)
	REAL(KIND=LDP) ED(ND)			!Electron density
!
! NB: On exit COOL needs to be multiplied by COOL.
!
	REAL(KIND=LDP) I_COOL,K_COOL
	REAL(KIND=LDP) COOL(ND)			!Net collisional cooling
	LOGICAL COMPUTE_BA,FIXED_T,LAST_ITERATION
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	EXTERNAL OMEGA_COL
!
	INTEGER I,J,K
	INTEGER L,U
	REAL(KIND=LDP) X
	REAL(KIND=LDP) BRAT
	REAL(KIND=LDP) CIJ,CJI,CII
	LOGICAL DO_dCOLdT
!
!$OMP PARALLEL WORKSHARE
	COL_S(:,:,:)=0.0_LDP  		!N_S, N_S, ND
	dCOL_S(:,:,:)=0.0_LDP  		!N_S, N_S, ND
!$OMP END PARALLEL WORKSHARE
!
! Loop over all depths. This is provided for consistency with other collisonal
! routines. We operate on a single depth at a time to minimize the size of
! OMEGA and dln_OMEGA_dlnT.
!
	DO K=1,ND
!
! The following correspond to the collision rates between levels I and J in
! the FULL atom.
!
!$OMP PARALLEL WORKSHARE
	  OMEGA_F(:,:)=0.0_LDP
	  dln_OMEGA_dlnT(:,:)=0.0_LDP
	  COL_F(:,:)=0.0_LDP
	  DCOL_F(:,:)=0.0_LDP
!$OMP END PARALLEL WORKSHARE
!
	  DO_dCOLdT=.TRUE.
	  IF(.NOT. COMPUTE_BA .OR. FIXED_T)DO_dCOLdT=.FALSE.
!
	  CALL TUNE(1,'OMEGA')
	  CALL OMEGA_COL(OMEGA_F,dln_OMEGA_dlnT,EDGE_F,AHYD_F,GHYD_F,
	1                  LEVNAME_F,ZION,N_F,T(K),
	1                  ID,COL_FILE)
	  CALL TUNE(2,'OMEGA')
!
! Allow for collisional ionization. OMEGA for the collisional ionization
! is assumed to be define in exactly the same way as for collisional
! excitation. NB: The variation due to a change in X cancels out with the
! change in HNST_F.
!
	  DO I=1,N_F
	    L=F_TO_S_MAPPING(I)
	    CII=8.63E-08_LDP*ED(K)*OMEGA_F(I,I)/GHYD_F(I)/SQRT(T(K))
	    X=HDKT*EDGE_F(I)/T(K)
	    CII=CII*EXP(-X)*HNST_F_ON_S(I,K)
	    COL_S(L,L,K)=COL_S(L,L,K)+CII
	    DCOL_S(L,L,K)=DCOL_S(L,L,K)+
	1       CII*( dln_OMEGA_dlnT(I,I) - 2.0_LDP - dlnHNST_S_dlnT(L,K) )/T(K)
	    COOL(K)=COOL(K)+EDGE_F(I)*CII*(HN_S(L,K)-HNST_S(L,K))
	  END DO
	  K_COOL=0.0_LDP
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) REDUCTION(+:K_COOL)
!$OMP1 PRIVATE(I,J,I_COOL,CII,CIJ,CJI,X,L,U,BRAT)
	  DO I=1,N_F-1
	    I_COOL=0.0_LDP
	    DO J=I+1,N_F
!
	      CIJ=8.63E-08_LDP*ED(K)*OMEGA_F(I,J)/SQRT(T(K))
	      CJI=CIJ/GHYD_F(J)
	      X=HDKT*(EDGE_F(I)-EDGE_F(J))/T(K)
	      CIJ=CIJ*EXP(-X)/GHYD_F(I)
!
! Allow for collisional ionization through level dissolution.
!
	      L=F_TO_S_MAPPING(I)
	      CII=CIJ*HNST_F_ON_S(I,K)*(1.0_LDP-W_F(J,K)/W_F(I,K))
	      COL_F(I,I)=COL_F(I,I)+CII
	      IF(DO_dCOLdT)THEN
	        DCOL_F(I,I)=CII*( dln_OMEGA_dlnT(I,J) +
	1        X -2.0_LDP - HDKT*EDGE_F(I)/T(K)-dlnHNST_S_dlnT(L,K) )/T(K)
	      END IF
!
! We use the full ionization energy because of the following argument.
! We think of the ionization as a 3 body process. The extra ionization
! energy comes from the third electron and hence is lost from the electron
! thermal pool.
!
	      I_COOL=I_COOL+EDGE_F(I)*CII*(HN_S(L,K)-HNST_S(L,K))
!
! The following section is independent of the atomic structure. We use
! F_TO_S_MAPPING to describe how the FULL atom is mapped on to the smaller
! atom. Note that the collsion rate is only important when L is not equal to
! U.
!
	      L=F_TO_S_MAPPING(I)
	      U=F_TO_S_MAPPING(J)
!
! BRAT is the ration of b (level in full atom) to b (in the super level).
! Rates are written using BRAT since BRAT is assumed to remain constant
! during the linearization.
!
	      IF(L .NE. U)THEN
	        BRAT=(HN_F(I,K)/HN_S(L,K))/HNST_F_ON_S(I,K)
	        CIJ=CIJ*BRAT*HNST_F_ON_S(I,K)*(W_F(J,K)/W_F(I,K))
!
	        BRAT=(HN_F(J,K)/HN_S(U,K))/HNST_F_ON_S(J,K)
	        CJI=CJI*BRAT*HNST_F_ON_S(J,K)
!
	        COL_F(I,J)=CIJ
	        COL_F(J,I)=CJI
!
! We now compute the derivatives of the collision rates with respect to T.
!
	        IF(DO_dCOLdT)THEN
	          DCOL_F(I,J)=CIJ*( dln_OMEGA_dlnT(I,J) + X -
	1             2.0_LDP - HDKT*EDGE_F(I)/T(K)-dlnHNST_S_dlnT(L,K) )/T(K)
	          DCOL_F(J,I)=CJI*( dln_OMEGA_dlnT(I,J) -
	1           2.0_LDP - HDKT*EDGE_F(J)/T(K)-dlnHNST_S_dlnT(U,K) )/T(K)
	        END IF
!
	        I_COOL=I_COOL+(HN_S(L,K)*CIJ-HN_S(U,K)*CJI)*(EDGE_F(I)-EDGE_F(J))
	      END IF
	    END DO		!J
	    K_COOL=K_COOL+I_COOL
	  END DO		!I
	  COOL(K)=COOL(K)+K_COOL
!
	  DO J=1,N_F
	    U=F_TO_S_MAPPING(J)
	    DO I=1,N_F
	      L=F_TO_S_MAPPING(I)
	      COL_S(L,U,K)=COL_S(L,U,K)+COL_F(I,J)
	    END DO
	  END DO
!
	  IF(DO_dCOLdT)THEN
	    DO J=1,N_F
	      U=F_TO_S_MAPPING(J)
	      DO I=1,N_F
	        L=F_TO_S_MAPPING(I)
	        DCOL_S(L,U,K)=DCOL_S(L,U,K)+DCOL_F(I,J)
	      END DO
	    END DO
	  END IF
!
	END DO
!
	RETURN
	END
