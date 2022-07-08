!
! Subroutine to replace a section of data by a straight line. Average of aound current
! location, oy y cursor vlaue can be used. Controls allow plot to moved
! and shifted (in X).
!
! Options:
!      (1) Cursor control
!            Using cursors the user specifies the range limits (X coordinates).
!            Two chocies for Y.
!                 Y location indicates the continuum
!             or  Continuum defined by average (over NPIX) of the data at the X coordinates.
!
	SUBROUTINE PG_REPLACE_DATA_INT(
	1             XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE,
	1             PEN_COL,PEN_OFFSET,REVERSE_PLOT_ORDER)
	USE MOD_CURVE_DATA
	USE MOD_EW_VARIABLES
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Created: 30-Jun-2022 - Can no append transition name to EW file.
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ISIX=6
!
! Most of the passed variables are simply required for redrawing the
! plot.
!
	REAL*4 XPAR(2),XT(2)
	REAL*4 YPAR(2),YT(2)
        REAL*4 XINC,XNUMST
        REAL*4 YINC,YNUMST
	REAL*4 dY
        REAL*4 EXPCHAR,TICK_FAC
        INTEGER IDX,IXTICK
        INTEGER IDY,IYTICK
!
	LOGICAL TITONRHS
	LOGICAL NORMAL_R_Y_AXIS
!
	CHARACTER(LEN=*) XLABEL,YLABEL
        CHARACTER(LEN=*) LOG_AXIS
	CHARACTER(LEN=*) OPTION
	CHARACTER(LEN=*) XLAB_FILE,YLAB_FILE
!
! Local vairables.
!
	INTEGER, SAVE :: NPIX=3  	!Integraton band pass around line limits
	INTEGER, SAVE :: IP=1
	REAL*4,  SAVE :: LOW_CONT=0.998
	REAL*4,  SAVE :: HIGH_CONT=1.002
!
	INTEGER PEN_COL(0:MAX_PLTS)
	INTEGER PEN_OFFSET
	LOGICAL REVERSE_PLOT_ORDER
!
	INTEGER IST,IEND		!Line limits in pixel space
	REAL*4 XLOC,YLOC		!Used to read cursor location
	REAL*4 YST,YEND  		!Continuum flux at line limits
!
! In the following I use Y and F interchangably.
! Also X will normally be Lambda.
!
	REAL*4 dX               !X spacing (pixel centered)
	REAL*4 YMEAN  		!Average value of (F-Fc)
	REAL*4 XVAL  		!Current X value
	REAL*4 YVAL  		!Current Y value
	REAL*4 YINT
	REAL*4 SLOPE  		!Continuum slope
	REAL*4 T1		!Work variable
!
	REAL*4, SAVE :: CONT_ACC=0.2		!percentage error in continuum
!
        INTEGER PGCURS
        INTEGER CURSERR
	INTEGER CUR_LW
	INTEGER, SAVE :: IPEN=1
	CHARACTER(LEN=1) CURSVAL
!
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, SAVE :: FIRST_WRITE=.TRUE.
	INTEGER I,J,K
	INTEGER, SAVE :: LUIN=0
	INTEGER, SAVE :: LUOUT=0
!
	INTEGER GET_INDX_SP
	EXTERNAL GET_INDX_SP
!
	LOGICAL RESET_DEFAULTS
	LOGICAL END_FILE
	LOGICAL FILE_PRES
	CHARACTER(LEN=80) LOC_FILE_WITH_LINE_LIMS
	CHARACTER(LEN=80) OUT_FILE 
	CHARACTER(LEN=80) STRING
!
! XLOC,YLOC will set inital cursor location.
!
	XLOC=0.5D0*(XPAR(1)+XPAR(2))
	YLOC=0.5D0*(YPAR(1)+YPAR(2))
	dY=0.05*(YPAR(2)-YPAR(1))
!
	RESET_DEFAULTS=.FALSE.
	CALL GEN_IN(RESET_DEFAULTS,'Reset default parameters/settings?')
	IF(RESET_DEFAULTS)THEN
	  IF(NPLTS .EQ. 1)THEN
	    IP=1
	  ELSE
	    CALL GEN_IN(IP,'Plot for fitting')
	  END IF
!
! Used if average data within X-limits is used to define the level.
!
	  CALL GEN_IN(NPIX,'Number of pixels at X location to average data (must be odd)')
	  IF(MOD(NPIX,2) .EQ. 0)THEN
	    NPIX=NPIX+1
	    WRITE(6,*)'NPIX increase by 1 to make odd; NPIX=',NPIX
	  END IF
	END IF
!
! Open output file. Data is appended if it alread exists.
!
	IF(LUOUT .EQ. 0)THEN
	  OUT_FILE='INT_REG'
	  CALL GET_LU(LUOUT,'LUOUT in DO_CURSOR_EW_V2')
	  CALL GEN_IN(OUT_FILE,'File to OUTPUT EWs etc')
	  INQUIRE(FILE=OUT_FILE,EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)'File already exists -- appending new data'
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	  ELSE
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='NEW',ACTION='WRITE')
	    CALL WRITE_EW_HEADER(LUOUT)
	  END IF
	END IF
!
	WRITE(6,'(A)')' '
	I=6; CALL WRITE_EW_HEADER(I)
	WRITE(6,'(A)')' '
!
	DO WHILE(1 .EQ. 1)
!
! Change pen colors for defining limits of the next line.
!
	  IF(IPEN .EQ. 6)THEN
	    IPEN=1
	  ELSE
	    IPEN=6
	  END IF
	  CALL PGSCI(IPEN)  
!
! Get start of line.
!
	  CURSERR = PGCURS(XLOC,YLOC,CURSVAL)
	  IST=GET_INDX_SP(XLOC,CD(IP)%XVEC,NPTS(IP))
	  XST=XLOC; YST=YLOC
	  IF(XST .LT. XPAR(1))THEN
	    WRITE(6,*)'Invalid X start location -- start again'
	    GOTO 1000
	  END IF
	  IF(CURSVAL .EQ. 'e')THEN
	    EXIT
	  ELSE IF(CURSVAL .EQ. 's')THEN
	    GOTO 1000
	  ELSE IF(CURSVAL .EQ. 'y')THEN
	    CALL PGPT(IONE,XLOC,YLOC,ITWO)
!
! Y value set by average at cursor location
!
	  ELSE IF(CURSVAL .EQ. 'a')THEN
	    YST=0.0
	    DO J=IST-NPIX/2,IST+NPIX/2
	      YST=YST+CD(IP)%DATA(J)
	    END DO
	    YST=YST/NPIX
	    YT(1)=YST-dY; YT(2)=YST+dY
	    XT(1)=XLOC; XT(2)=XLOC
	    CALL PGLINE(2,XT,YT)
!
! Options to move/expand/contract plot in X direction.
!
	  ELSE IF(CURSVAL .EQ. 'n' .OR.
	1           CURSVAL .EQ. 'p' .OR.
	1           CURSVAL .EQ. 'x' .OR.
	1           CURSVAL .EQ. 'c' .OR.
	1           CURSVAL .EQ. 'r')THEN
	    T1=XPAR(2)-XPAR(1)
	    IF(CURSVAL .EQ. 'n')THEN
	      XPAR(1)=XPAR(2)-XINC
	      XPAR(2)=XPAR(1)+T1
	    ELSE IF(CURSVAL .EQ. 'p')THEN
	      XPAR(1)=XPAR(1)+XINC-(XPAR(2)-XPAR(1))
	      XPAR(2)=XPAR(1)+T1
	    ELSE IF(CURSVAL .EQ. 'x')THEN
	      XPAR(2)=XPAR(2)+XINC
	    ELSE IF(CURSVAL .EQ. 'c')THEN
	      XPAR(2)=XPAR(2)-XINC
	    END IF
	    XLOC=XPAR(1)+XINC
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
            CALL PGERAS
	    IF(REVERSE_PLOT_ORDER)THEN
	      DO J=NPLTS,1,-1
	        K=PEN_COL(J)+PEN_OFFSET
	        CALL PGSCI(K)
                CALL PGLINE(NPTS(J),CD(J)%XVEC,CD(J)%DATA)
	      END DO
	    ELSE
	      DO J=1,NPLTS
	        K=PEN_COL(J)+PEN_OFFSET
	        CALL PGSCI(K)
                CALL PGLINE(NPTS(J),CD(J)%XVEC,CD(J)%DATA)
	      END DO
	    END IF
!
! Add line IDs
!
	    CALL DRAW_LINE_IDS_V2(XPAR,YPAR,EXPCHAR,IONE,ISIX)
            CALL PGSCH(EXPCHAR)		!Reset character size
!
	    CALL PGSCI(IONE)
            CALL MONBORD_V4(XPAR,XINC,XNUMST,IXTICK,IDX,
	1           YPAR,YINC,YNUMST,IYTICK,IDY,
	1           TICK_FAC,EXPCHAR,
	1           XLABEL,YLABEL,TITLE,N_TITLE,TITONRHS,
	1           LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1           XLAB_FILE,YLAB_FILE)
	    GOTO 1000
	  ELSE
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)' Cursor value not recognized'	
	    CALL PRINT_CURSOR_DESC
	    GOTO 1000
	  END IF
!
! Get end of line.
!
	  CURSERR = PGCURS(XLOC,YLOC,CURSVAL)
	  IF(XLOC .GT. XPAR(2))THEN
	    WRITE(6,*)'Invalid X end location -- start again'
	    GOTO 1000
	  END IF
	  IEND=GET_INDX_SP(XLOC,CD(IP)%XVEC,NPTS(IP))
	  XEND=XLOC; YEND=YLOC
	  IF(CURSVAL .EQ. 'e')THEN
	    EXIT
	  ELSE IF(CURSVAL .EQ. 's')THEN
	    GOTO 1000
	  ELSE IF(CURSVAL .EQ. 'y')THEN
	    CALL PGPT(IONE,XLOC,YLOC,ITWO)
	  ELSE IF(CURSVAL .EQ. 'a')THEN
	    YEND=0.0
	    DO J=IEND-NPIX/2,IEND+NPIX/2
	      YEND=YEND+CD(IP)%DATA(J)
	    END DO
	    YEND=YEND/NPIX
	    YT(1)=YEND-dY; YT(2)=YEND+dY
	    XT(1)=XLOC; XT(2)=XLOC
	    CALL PGLINE(2,XT,YT)
	  ELSE
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)' Cursor value not recognized'	
	    CALL PRINT_CURSOR_DESC
	    GOTO 1000
	  END IF
!
	  IF(IST .EQ. IEND)THEN
	    WRITE(6,*)'Error -- limits on line are identical'
	    WRITE(6,*)'Need to redefine the line'
	    GOTO 1000
	  END IF
!
! We now determine the line parameters
!
	  EW=0.0;  XMEAN=0.0
	  SLOPE=(YEND-YST)/(CD(IP)%XVEC(IEND)-CD(IP)%XVEC(IST))
	  DO I=IST,IEND
	    YVAL=YST+(CD(IP)%XVEC(I)-CD(IP)%XVEC(IST))*SLOPE
	    CD(IP)%DATA(I)=YVAL
	  END DO
	  WRITE(6,*)XST,XEND
	  WRITE(LUOUT,*)XST,XEND
!
1000	  CONTINUE
	END DO
!
2000	CONTINUE
	RETURN
	CONTAINS
!
	SUBROUTINE PRINT_CURSOR_DESC
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')' Program assumes X axis increases with X index '
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')' Use cursor to define left and right side of line:'
	  WRITE(6,'(A)')'   y - uses cursor Y value for location'
	  WRITE(6,'(A)')'   s - reset and start selection again'
	  WRITE(6,'(A)')'   a - Y value is average over data at X location'
	  WRITE(6,'(A)')'   e - exit line selection'
	  WRITE(6,'(A)')' '
	  WRITE(6,'(A)')'Plot movement'
	  WRITE(6,'(A)')'   n(ext)      -- shift window to left by X increment'
	  WRITE(6,'(A)')'   p(revious)  -- shift window to right by X increment'
	  WRITE(6,'(A)')'   x(extend)   -- expand window by increment but no shift'
 	  WRITE(6,'(A)')'   c(ontract)  -- contract window by increment but no shift'
 	  WRITE(6,'(A)')'   r(edraw)    -- redraw display'
	  WRITE(6,'(A)')DEF_PEN
	  RETURN
	END SUBROUTINE PRINT_CURSOR_DESC
	END SUBROUTINE PG_REPLACE_DATA_INT
