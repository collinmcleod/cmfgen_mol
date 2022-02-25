	SUBROUTINE REP_DATA(IP_IN,IP_OUT)
	USE MOD_CURVE_DATA
        USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	INTEGER IP_IN
	INTEGER IP_OUT
!
	REAL*4 XV(2*NPTS(IP_IN))
	REAL*4 YV(2*NPTS(IP_IN))
	REAL*4 TMP_XV(2*NPTS(IP_IN))
	REAL*4 WRK_YV(2*NPTS(IP_IN))
!
	INTEGER, PARAMETER :: IONE=1
!
	REAL*4 XVAL,YVAL,SYMB_EXP_FAC
	REAL*4 T1
	INTEGER ND
	INTEGER I,J
	INTEGER PGCURS
	INTEGER CURSERR
	CHARACTER(LEN=1) CURSVAL
!
	ND=NPTS(IP_IN)
	XV(1:ND)=CD(IP_IN)%XVEC
	YV(1:ND)=CD(IP_IN)%DATA
!
	WRITE(6,'(A)')' '
	WRITE(6,*)'Use ''r'' to replace a data point'
	WRITE(6,*)'Use ''a'' to add a data point'
	WRITE(6,*)'Use ''d'' to delete a data point'
	WRITE(6,*)'Use ''m'' move next, and points further out/inward by a fixed amount'
	WRITE(6,*)'Use ''v'' move next, and points further  up/down by a fixed amount in V'
	WRITE(6,*)'Use ''e'' to exit'
	WRITE(6,'(A)')DEF_PEN
!
	XVAL=5.0; CALL PGSCH(XVAL)
	DO WHILE(1 .EQ. 1)				!Multiple cursor entries
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  WRITE(6,*)'Cursor values are:',XVAL,YVAL
	  IF(CURSVAL .EQ. 'e' .OR. CURSVAL .EQ. 'E')EXIT
!
	  IF(CURSVAL .EQ. 'r' .OR. CURSVAL .EQ. 'R')THEN
	    T1=XVAL
	    TMP_XV(1:ND)=ABS(XV(1:ND)-XVAL)
	    I=MINLOC(TMP_XV(1:ND),IONE)
	    YV(I)=YVAL
	    J=1; CALL PGSCI(J)
	    CALL PGPT(IONE,XVAL,YVAL,IONE)
	    WRITE(6,*)'Replaced YV at  XV=',XV(I)
!
	  ELSE IF(CURSVAL .EQ. 'd' .OR. CURSVAL .EQ. 'D')THEN
	    T1=XVAL
	    TMP_XV(1:ND)=ABS(XV(1:ND)-T1)
	    I=MINLOC(TMP_XV(1:ND),IONE)
	    XV(I:ND-1)=XV(I+1:ND)
	    YV(I:ND-1)=YV(I+1:ND)
	    ND=ND-1
	    WRITE(6,*)'Deleted XV=',XV(I),' from grid'
!
	  ELSE IF(CURSVAL .EQ. 'a' .OR. CURSVAL .EQ. 'A')THEN
	    DO I=1,ND-1
	      IF( (XV(I)-XVAL)*(XV(I+1)-XVAL) .LT. 0 )THEN
	        DO J=ND,I+1,-1
                  XV(J+1)=XV(J)
                  YV(J+1)=YV(J)
	        END DO
	        XV(I+1)=XVAL; XV(I+1)=YVAL
	        ND=ND+1
	        J=1; CALL PGSCI(J)
	        CALL PGPT(IONE,XVAL,YVAL,IONE)
	        WRITE(6,*)'Add XV=',XV(I),' to grid'
	        EXIT
	      END IF
	    END DO
!
	  ELSE IF(CURSVAL .EQ. 'm' .OR. CURSVAL .EQ. 'M')THEN
	     WRK_YV(1:ND)=ABS(YV(1:ND)-YVAL)
	     I=MINLOC(WRK_YV(1:ND),IONE)
	     J=1; CALL PGSCI(J)
	     T1=XVAL-XV(I)
	     DO J=1,I
	       XV(J)=XV(J)+T1
	       XVAL=XV(J); YVAL=YV(J)
	       CALL PGPT(IONE,XVAL,YVAL,IONE)
	     END DO	
	     WRITE(6,*)'Shited XV grid by',T1
!
	  ELSE IF(CURSVAL .EQ. 'v' .OR. CURSVAL .EQ. 'V')THEN
	     TMP_XV(1:ND)=ABS(XV(1:ND)-XVAL)
	     I=MINLOC(TMP_XV(1:ND),IONE)
	     J=1; CALL PGSCI(J)
	     T1=YVAL-YV(I)
	     DO J=1,I
	       YV(J)=YV(J)+T1
	       XVAL=XV(J); YVAL=YV(J)
	       CALL PGPT(IONE,XVAL,YVAL,IONE)
	     END DO	
	     WRITE(6,*)'Shifted Y grid by',T1
!
	   ELSE
	      WRITE(6,*)RED_PEN
	      WRITE(6,*)'Error - use r(eplace), a(dd), d(elete), e'
	      WRITE(6,*)DEF_PEN
	   END IF
	 END DO
!
	 IF(IP_OUT .GT. NPLTS+1)IP_OUT=NPLTS+1
	 IF(IP_OUT .GT. NPLTS)NPLTS=NPLTS+1
	 IF(ALLOCATED(CD(IP_OUT)%XVEC))THEN
	   DEALLOCATE(CD(IP_OUT)%XVEC,CD(IP_OUT)%DATA)
	 END IF
	 IF(ALLOCATED(CD(IP_OUT)%EMIN))THEN
	   DEALLOCATE(CD(IP_OUT)%EMIN,CD(IP_OUT)%EMIN)
	 END IF
	 NPTS(IP_OUT)=ND
	 ALLOCATE(CD(IP_OUT)%XVEC(ND))
	 ALLOCATE(CD(IP_OUT)%DATA(ND))
	 CD(IP_OUT)%XVEC(1:ND)=XV(1:ND)
	 CD(IP_OUT)%DATA(1:ND)=YV(1:ND)
	 ERR(IP_OUT)=.FALSE.
	 CD(IP_OUT)%CURVE_ID=' '
!
	RETURN
	END
