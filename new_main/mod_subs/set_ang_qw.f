!
! Declares arrays required for the angular quadrature weights, and
! computes their vaule. Routine can be called more than once, although
! the size of thei declared arrays must not change.
!
	SUBROUTINE SET_ANG_QW(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,
	1                 TRAPFORJ,ACCURATE)
	USE SET_KIND_MODULE
	USE ANG_QW_MOD
	IMPLICIT NONE
!
! Created 02-May-2004
!
	INTEGER NC,ND,NP
	INTEGER NCEXT,NDEXT,NPEXT
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) REXT(NDEXT)
	LOGICAL TRAPFORJ
	LOGICAL ACCURATE
!
! External function calls.
!
	EXTERNAL JWEIGHT,HWEIGHT,KWEIGHT,NWEIGHT
	EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
	EXTERNAL ERROR_LU
	INTEGER ERROR_LU
!
	REAL(KIND=LDP) W1(NPEXT),W2(NPEXT),W3(NPEXT)		!Work vectors
!
	INTEGER LU_ER
	INTEGER LS
	INTEGER IOS
	LOGICAL MID
!
! Angular quadrature weights.
!
	IF(.NOT. ALLOCATED(P))THEN
!	
	  ALLOCATE(P(NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(AQW(ND,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(HQW(ND,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(KQW(ND,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(NQW(ND,NP),STAT=IOS)
!
! Defined at the mid points of the mesh.
!
	  IF(IOS .EQ. 0)ALLOCATE(HMIDQW(ND,NP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(NMIDQW(ND,NP),STAT=IOS)
!
! If required, these arrays shoukd have size NDEXT*NPEXT
!
	  IF(ACCURATE)THEN
	    IF(IOS .EQ. 0)ALLOCATE(PEXT(NPEXT),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE(AQWEXT(NDEXT,NPEXT),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE(HQWEXT(NDEXT,NPEXT),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE(KQWEXT(NDEXT,NPEXT),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE(NQWEXT(NDEXT,NPEXT),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE(HMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE(NMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	  END IF
!
	  IF(IOS .NE. 0)THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Unable to allocate memory for JQW in SET_ANG_QW'
	    STOP
	  END IF
	END IF
!
! Compute impact parameter values P
!
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
!
! Compute the angular quadrature weights
!
	IF(TRAPFORJ)THEN
	    CALL NORDANGQW(AQW,R,P,W1,W2,W3,NC,ND,NP,JTRPWGT)
	    CALL NORDANGQW(HQW,R,P,W1,W2,W3,NC,ND,NP,HTRPWGT)
	    CALL NORDANGQW(KQW,R,P,W1,W2,W3,NC,ND,NP,KTRPWGT)
	    CALL NORDANGQW(NQW,R,P,W1,W2,W3,NC,ND,NP,NTRPWGT)
	    MID=.TRUE.
	    CALL GENANGQW(HMIDQW,R,P,W1,W2,W3,NC,ND,NP,HTRPWGT,MID)
	    CALL GENANGQW(NMIDQW,R,P,W1,W2,W3,NC,ND,NP,NTRPWGT,MID)
	  ELSE
	    CALL NORDANGQW(AQW,R,P,W1,W2,W3,NC,ND,NP,JWEIGHT)
	    CALL NORDANGQW(HQW,R,P,W1,W2,W3,NC,ND,NP,HWEIGHT)
	    CALL NORDANGQW(KQW,R,P,W1,W2,W3,NC,ND,NP,KWEIGHT)
	    CALL NORDANGQW(NQW,R,P,W1,W2,W3,NC,ND,NP,NWEIGHT)
	    MID=.TRUE.
	    CALL GENANGQW(HMIDQW,R,P,W1,W2,W3,NC,ND,NP,HWEIGHT,MID)
	    CALL GENANGQW(NMIDQW,R,P,W1,W2,W3,NC,ND,NP,NWEIGHT,MID)
	  END IF
!
! Compute impact parameter values P
!
	IF(ACCURATE)THEN
	  CALL IMPAR(PEXT,REXT,REXT(NDEXT),NCEXT,NDEXT,NPEXT)
	  IF(TRAPFORJ)THEN
	    CALL NORDANGQW(AQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,JTRPWGT)
	    CALL NORDANGQW(HQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,HTRPWGT)
	    CALL NORDANGQW(KQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,KTRPWGT)
	    CALL NORDANGQW(NQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,NTRPWGT)
	    MID=.TRUE.
	    CALL GENANGQW(HMIDQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,HTRPWGT,MID)
	    CALL GENANGQW(NMIDQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,NTRPWGT,MID)
	  ELSE IF(ACCURATE)THEN
	    CALL NORDANGQW(AQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,JWEIGHT)
	    CALL NORDANGQW(HQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,HWEIGHT)
	    CALL NORDANGQW(KQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,KWEIGHT)
	    CALL NORDANGQW(NQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,NWEIGHT)
	    MID=.TRUE.
	    CALL GENANGQW(HMIDQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,HWEIGHT,MID)
	    CALL GENANGQW(NMIDQWEXT,REXT,PEXT,W1,W2,W3,NCEXT,NDEXT,NPEXT,NWEIGHT,MID)
	  END IF
	END IF
!
! Allocate arrays and vectors for computing observed fluxes.
!
	IF(ACCURATE)THEN
	  NP_OBS_MAX=NPEXT+12
	ELSE
	  NP_OBS_MAX=NP+12
	END IF
!
	IF(.NOT. ALLOCATED(P_OBS))THEN
	  ALLOCATE (P_OBS(NP_OBS_MAX),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (IPLUS_STORE(NST_CMF,NP_OBS_MAX),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (IPLUS(NP_OBS_MAX),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (MU_AT_RMAX(NP_OBS_MAX),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HQW_AT_RMAX(NP_OBS_MAX),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Unable to allocate memory for P_OBS etc in SET_ANG_QW'
	    STOP
	  END IF
	END IF
!
! Used when computing the observed fluxes. These will get overwritten
! if we do an accurate comoving frame soluton using CMF_FORM_SOL.
!
	IF(ACCURATE)THEN
	  DO LS=1,NPEXT
	    MU_AT_RMAX(LS)=SQRT( 1.0_LDP -(PEXT(LS)/REXT(1))**2 )
	    HQW_AT_RMAX(LS)=HQWEXT(1,LS)
	  END DO
	ELSE
	  DO LS=1,NP
	    MU_AT_RMAX(LS)=SQRT( 1.0_LDP -(P(LS)/R(1))**2 )
	    HQW_AT_RMAX(LS)=HQW(1,LS)
	  END DO
	END IF
!
	RETURN
	END
