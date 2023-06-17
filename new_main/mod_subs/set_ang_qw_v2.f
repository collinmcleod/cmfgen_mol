!
! Declares arrays required for the angular quadrature weights, and
! computes their vaule. Routine can be called more than once, although
! the size of thei declared arrays must not change.
!
	SUBROUTINE SET_ANG_QW_V2(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,
	1              R_PNT_SRCE,NC_PNT_SRCE,TRAPFORJ,ACCURATE)
	USE ANG_QW_MOD
	IMPLICIT NONE
!
! Altered 16-Jun-2023 : Changed angular quadrature weight routines to V2.
! Alteerd 16-Jun-2023 - Changed quadrature weight routines to V2.
! Altered 17-AUg-2109 - Placed in IBIS
! Created 06-Jun-2019 (on OSIRIS) - added PNT SRCE options to call.
! Created 02-May-2004
!
	INTEGER NC,ND,NP
	INTEGER NCEXT,NDEXT,NPEXT
	INTEGER NC_PNT_SRCE
	REAL*8 R(ND)
	REAL*8 REXT(NDEXT)
	REAL*8 R_PNT_SRCE
	LOGICAL TRAPFORJ
	LOGICAL ACCURATE 
!
! External function calls.
!
	EXTERNAL JWEIGHT_V2,HWEIGHT_V2,KWEIGHT_V2,NWEIGHT_V2
	EXTERNAL JTRPWGT_V2,HTRPWGT_V2,KTRPWGT_V2,NTRPWGT_V2
	EXTERNAL ERROR_LU
	INTEGER ERROR_LU
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
	CALL IMPAR_WITH_PNT_SRCE(P,R,R(ND),NC,ND,NP,R_PNT_SRCE,NC_PNT_SRCE)
!
! Compute the angular quadrature weights
!
	IF(TRAPFORJ)THEN
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,JTRPWGT_V2)
	    CALL NORDANGQW_V2(HQW,R,P,NC,ND,NP,HTRPWGT_V2)
	    CALL NORDANGQW_V2(KQW,R,P,NC,ND,NP,KTRPWGT_V2)
	    CALL NORDANGQW_V2(NQW,R,P,NC,ND,NP,NTRPWGT_V2)
	    MID=.TRUE.
	    CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HTRPWGT_V2,MID)
	    CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NTRPWGT_V2,MID)
	  ELSE
	    CALL NORDANGQW_V2(AQW,R,P,NC,ND,NP,JWEIGHT_V2)
	    CALL NORDANGQW_V2(HQW,R,P,NC,ND,NP,HWEIGHT_V2)
	    CALL NORDANGQW_V2(KQW,R,P,NC,ND,NP,KWEIGHT_V2)
	    CALL NORDANGQW_V2(NQW,R,P,NC,ND,NP,NWEIGHT_V2)
	    MID=.TRUE.
	    CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HWEIGHT_V2,MID)
	    CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NWEIGHT_V2,MID)
	  END IF
!
! Compute impact parameter values P
!
	IF(ACCURATE)THEN
	  CALL IMPAR_WITH_PNT_SRCE(PEXT,REXT,REXT(NDEXT),NCEXT,NDEXT,NPEXT,R_PNT_SRCE,NC_PNT_SRCE)
	  IF(TRAPFORJ)THEN
	    CALL NORDANGQW_V2(AQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,JTRPWGT_V2)
	    CALL NORDANGQW_V2(HQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HTRPWGT_V2)
	    CALL NORDANGQW_V2(KQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,KTRPWGT_V2)
	    CALL NORDANGQW_V2(NQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NTRPWGT_V2)
	    MID=.TRUE.
	    CALL GENANGQW_V2(HMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HTRPWGT_V2,MID)
	    CALL GENANGQW_V2(NMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NTRPWGT_V2,MID)
	  ELSE IF(ACCURATE)THEN
	    CALL NORDANGQW_V2(AQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,JWEIGHT_V2)
	    CALL NORDANGQW_V2(HQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HWEIGHT_V2)
	    CALL NORDANGQW_V2(KQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,KWEIGHT_V2)
	    CALL NORDANGQW_V2(NQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NWEIGHT_V2)
	    MID=.TRUE.
	    CALL GENANGQW_V2(HMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,HWEIGHT_V2,MID)
	    CALL GENANGQW_V2(NMIDQWEXT,REXT,PEXT,NCEXT,NDEXT,NPEXT,NWEIGHT_V2,MID)
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
	    MU_AT_RMAX(LS)=SQRT( 1.0D0 -(PEXT(LS)/REXT(1))**2 )
	    HQW_AT_RMAX(LS)=HQWEXT(1,LS)
	  END DO
	ELSE
	  DO LS=1,NP
	    MU_AT_RMAX(LS)=SQRT( 1.0D0 -(P(LS)/R(1))**2 )
	    HQW_AT_RMAX(LS)=HQW(1,LS)
	  END DO
	END IF
!
	RETURN
	END
