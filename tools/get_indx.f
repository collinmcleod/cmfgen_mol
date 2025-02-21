C
C Returns index value J such that RVAL lies between R(J) and
C R(J+1).
C
	INTEGER FUNCTION GET_INDX_DP(RVAL,R,NW)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 17-Jan-2018 - Error reporting improved (added to IBIS 17-Mar-2018)
C Altered 19-May-2014 - Changed error unit form 2 to 6 (25-May-2014).
C Altered 02-May-1990 - Bug fix: If on two consecutive calls, RVAL was the
C                       same but the R array different, the returned index
C                       was incorrect (referring to the first array). No
C                       longer use RVAL_SAV, but directly check R array.
C                       Routine now correctly returns 1 or NW-1 if RVAL is
C                       outside R range. Previously only J and not GET_INDX_DP
C                       was being set. (Minor correction on 8/May/91).
C
C Altered 20-Nov-1990 - Call to SEED_RITE changed.
C Altered 12-Nov-1990 = RVAL check altered. Routine nolonger STOPS program.
C Altered 23-Apr-1990 - Check that RVAL is inside range given by R included.
C Created 11-Apr-1990
C
	INTEGER NW
	REAL(KIND=LDP) R(NW),RVAL
C
	INTEGER J,ILOW,IHIGH,ILOW_SAV
	SAVE ILOW_SAV
	DATA ILOW_SAV/1/
C
C
C NW may have changed bewteen this and previous call, hence ILOW_SAV may
C lie outside the valid range (1:NW).
C
	IF(ILOW_SAV .GE. NW)ILOW_SAV=NW-1	
C
C Have two cases to find RVAL for, depending on whether the array R
C increases or decreases with its index.
C
	IF(R(1) .GT. R(NW))THEN
C
C Need t compare RVAl with R array (and not RSAV) in case value is
C same, but array is different.
C
	  IF( (RVAL .GE.  R(ILOW_SAV+1)) .AND.
	1       (RVAL .LE. R(ILOW_SAV)) )THEN
            GET_INDX_DP=ILOW_SAV
	    RETURN
	  END IF
	  IF(RVAL .GT. R(1))THEN
	    WRITE(6,*)'Error in GET_INDX_DP - RVAL too large.'
	    WRITE(6,70)'R(1)=',R(1),'RVAL=',RVAL
	    WRITE(6,70)'RVAL=',RVAL
	    WRITE(6,*)'       NW=',NW
70	    FORMAT(1X,1P,(3X,A,E12.4))
	    RVAL=R(1)
	    GET_INDX_DP=1
	    RETURN
	  END IF
	  IF(RVAL .LT. R(NW))THEN
	    WRITE(6,*)'Error in GET_INDX_DP - RVAL too small.'
	    WRITE(6,70)'R(NW)=',R(NW),'RVAL=',RVAL
	    WRITE(6,70)'RVAL=',RVAL
	    WRITE(6,*)'       NW=',NW
	    RVAL=R(NW)
	    GET_INDX_DP=NW-1
	    RETURN
	  END IF
C
C Here LOW and HIGH refer to the size of the INDEX - not R which
C gets smaller as the index decreases.
C
	  ILOW=1
	  IHIGH=NW
	  DO WHILE( (IHIGH-ILOW) .GT. 1)
	    J=(ILOW+IHIGH)/2
	    IF(RVAL .LT. R(J))THEN
	      ILOW=J
	    ELSE
	      IHIGH=J
	    END IF
	  END DO
	  ILOW_SAV=ILOW
	  GET_INDX_DP=ILOW
	  RETURN
	ELSE
C
C Need t compare RVAl with R array (and not RSAV) in case value is
C same, but array is different.
C
	  IF( (RVAL .LE.  R(ILOW_SAV+1)) .AND.
	1       (RVAL .GE. R(ILOW_SAV)) )THEN
            GET_INDX_DP=ILOW_SAV
	    RETURN
	  END IF
	  IF(RVAL .LT. R(1))THEN
	    WRITE(6,*)'Error in GET_INDX_DP - RVAL too small.'
	    WRITE(6,70)'R(1)=',R(1),'RVAL=',RVAL
	    RVAL=R(1)
	    GET_INDX_DP=1
	    RETURN
	  END IF
	  IF(RVAL .GT. R(NW))THEN
	    WRITE(6,*)'Warning in GET_INDX_DP - RVAL too large.'
	    WRITE(6,70)'R(NW)=',R(NW),'RVAL=',RVAL
	    RVAL=R(NW)
	    GET_INDX_DP=NW-1
	    RETURN
	  END IF
C
	  ILOW=1
	  IHIGH=NW
	  DO WHILE( (IHIGH-ILOW) .GT. 1)
	    J=(ILOW+IHIGH)/2
	    IF(RVAL .LT. R(J))THEN
	      IHIGH=J
	    ELSE
	      ILOW=J
	    END IF
	  END DO
	  ILOW_SAV=ILOW
	  GET_INDX_DP=ILOW
	  RETURN
	END IF
C
	END FUNCTION GET_INDX_DP
C
C
C Returns index value J such that RVAL lies between R(J) and
C R(J+1).
C
	INTEGER FUNCTION GET_INDX_SP(RVAL,R,NW)
	USE SET_KIND_MODULE
	IMPLICIT NONE
C
C Altered 02-May-1990 - Bug fix: If on two consecutive calls, RVAL was the
C                       same but the R array different, the returned index
C                       was incorrect (referring to the first array). No
C                       longer use RVAL_SAV, but directly check R array.
C                       Routine now correctly returns 1 or NW-1 if RVAL is
C                       outside R range. Previously only J and not GET_INDX
C                       was being set. (Minor correction on 8/May/91).
C
C Altered 20-Nov-1990 - Call to SEED_RITE changed.
C Altered 12-Nov-1990 = RVAL check altered. Routine nolonger STOPS program.
C Altered 23-Apr-1990 - Check that RVAL is inside range given by R included.
C Created 11-Apr-1990
C
	INTEGER NW
	REAL*4 R(NW),RVAL
C
	INTEGER J,ILOW,IHIGH,ILOW_SAV
	SAVE ILOW_SAV
	DATA ILOW_SAV/1/
C
C
C NW may have changed bewteen this and previous call, hence ILOW_SAV may
C lie outside the valid range (1:NW).
C
	IF(ILOW_SAV .GE. NW)ILOW_SAV=NW-1	
C
C Have two cases to find RVAL for, depending on whether the array R
C increases or decreases with its index.
C
	IF(R(1) .GT. R(NW))THEN
C
C Need t compare RVAl with R array (and not RSAV) in case value is
C same, but array is different.
C
	  IF( (RVAL .GE.  R(ILOW_SAV+1)) .AND.
	1       (RVAL .LE. R(ILOW_SAV)) )THEN
            GET_INDX_SP=ILOW_SAV
	    RETURN
	  END IF
	  IF(RVAL .GT. R(1))THEN
	    WRITE(6,*)'Error in GET_INDX_SP - RVAL too large.'
	    WRITE(6,70)'R(1)=',R(1),'RVAL=',RVAL
70	    FORMAT(1X,1P,(3X,A,E12.4))
	    RVAL=R(1)
	    GET_INDX_SP=1
	    RETURN
	  END IF
	  IF(RVAL .LT. R(NW))THEN
	    WRITE(6,*)'Error in GET_INDX_SP - RVAL too small.'
	    WRITE(6,70)'R(NW)=',R(NW),'RVAL=',RVAL
	    RVAL=R(NW)
	    GET_INDX_SP=NW-1
	    RETURN
	  END IF
C
C Here LOW and HIGH refer to the size of the INDEX - not R which
C gets smaller as the index decreases.
C
	  ILOW=1
	  IHIGH=NW
	  DO WHILE( (IHIGH-ILOW) .GT. 1)
	    J=(ILOW+IHIGH)/2
	    IF(RVAL .LT. R(J))THEN
	      ILOW=J
	    ELSE
	      IHIGH=J
	    END IF
	  END DO
	  ILOW_SAV=ILOW
	  GET_INDX_SP=ILOW
	  RETURN
	ELSE
C
C Need t compare RVAl with R array (and not RSAV) in case value is
C same, but array is different.
C
	  IF( (RVAL .LE.  R(ILOW_SAV+1)) .AND.
	1       (RVAL .GE. R(ILOW_SAV)) )THEN
            GET_INDX_SP=ILOW_SAV
	    RETURN
	  END IF
	  IF(RVAL .LT. R(1))THEN
	    WRITE(6,*)'Error in GET_INDX_SP - RVAL too small.'
	    WRITE(6,70)'R(1)=',R(1),'RVAL=',RVAL
	    RVAL=R(1)
	    GET_INDX_SP=1
	    RETURN
	  END IF
	  IF(RVAL .GT. R(NW))THEN
	    WRITE(6,*)'Warning in GET_INDX_SP - RVAL too large.'
	    WRITE(6,70)'R(NW)=',R(NW),'RVAL=',RVAL
	    RVAL=R(NW)
	    GET_INDX_SP=NW-1
	    RETURN
	  END IF
C
	  ILOW=1
	  IHIGH=NW
	  DO WHILE( (IHIGH-ILOW) .GT. 1)
	    J=(ILOW+IHIGH)/2
	    IF(RVAL .LT. R(J))THEN
	      IHIGH=J
	    ELSE
	      ILOW=J
	    END IF
	  END DO
	  ILOW_SAV=ILOW
	  GET_INDX_SP=ILOW
	  RETURN
	END IF
C
	END FUNCTION GET_INDX_SP
