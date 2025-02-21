!
! General purpose routine to generate the folowing 2D maps:
!
!         (a) CONTOUR plots
!         (b) P,THETA vector plots
!         (c) Gray-scale plots
!         (d) Color images
!
! Based on GRAMON_PGPLOT: Some GRAMON options are obsolete for this routine.
!
	SUBROUTINE MAP_PLOT(DP_MAP,DP_THETA,DP_XVEC,DP_YVEC,XDIM,YDIM,
	1                   XLAB,YLAB,TITL,PASSED_OPT)
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created 24-Mar-1998
!
        INTEGER, PARAMETER :: MAXSTR=20
	INTEGER, PARAMETER :: MAXVEC=20
	INTEGER, PARAMETER :: MAXPEN=30
	INTEGER, PARAMETER :: MAX_PLTS=1
	INTEGER, PARAMETER :: NPLTS=1
	INTEGER, PARAMETER :: NSCL=100
!
	INTEGER NDEC
	EXTERNAL SPACING
	REAL*4 SPACING
!
	INTEGER XDIM,YDIM
	REAL(KIND=LDP) DP_MAP(XDIM,YDIM)
	REAL(KIND=LDP) DP_THETA(XDIM,YDIM)
	REAL(KIND=LDP) DP_XVEC(XDIM)
	REAL(KIND=LDP) DP_YVEC(YDIM)
!
	REAL*4 THETA(XDIM,YDIM)
	REAL*4, ALLOCATABLE :: XVEC(:)
	REAL*4, ALLOCATABLE :: YVEC(:)
	REAL*4, ALLOCATABLE :: TMP_VEC(:)
	REAL*4, ALLOCATABLE :: MAP(:,:)
	REAL*4, ALLOCATABLE :: TMP_MAP(:,:)
!
	REAL*4 X1,X2,Y1,Y2

	CHARACTER*1 TYPE_CURVE(MAX_PLTS)
!
	LOGICAL DO_ERROR
	CHARACTER*5 LOG_AXIS
!
	REAL*4 XINC,XNUMST,YNUMST
	REAL*4 YINC
	INTEGER IDX,IXTICK
	INTEGER IDY,IYTICK
	REAL*4 XPAR(2),YPAR(2),ZPAR(2)
	REAL*4 XT(2),YT(2)
	REAL*4 XMIN,XMAX,YMIN,YMAX
	REAL*4 XPAR_SAV(2)
	REAL*4 MAPMAX,MAPMIN
	REAL*4 SCALEPOL
!
	REAL*4 YPAR_R_AX(2)
	REAL*4 YNUMST_R_AX
	REAL*4 YINC_R_AX
	LOGICAL NORMAL_R_Y_AXIS
	INTEGER IDY_R_AX,IYTICK_R_AX
	CHARACTER*1 WHICH_Y_AX(MAX_PLTS)
	CHARACTER*80 YLABEL_R_AX
!
	REAL*4 SCALE_VEC(4,NSCL)
	LOGICAL, SAVE :: DO_SCALE
!
	INTEGER, PARAMETER :: N_TITLE=10
	CHARACTER(LEN=50), SAVE :: XLABEL,YLABEL,TITLE(N_TITLE)
	CHARACTER(LEN=*) XLAB,YLAB,TITL,PASSED_OPT
        CHARACTER(LEN=35) FILNAME
	CHARACTER(LEN=80) WK_STR,OPTION
	CHARACTER(LEN=6) TO_TEK
	CHARACTER(LEN=2) TO_VT
	CHARACTER(LEN=50) PLT_ID,RD_PLT_ID
	CHARACTER(LEN=80) XLAB_FILE,YLAB_FILE
!
	LOGICAL HAN
	INTEGER XHAN,YHAN,KST,KEND,LST,LEND
	REAL(KIND=LDP) SM_ARRAY(-10:10,-10:10)			!Array for smoothing
	REAL(KIND=LDP) FAC
	REAL(KIND=LDP) MAXSIZEY
	REAL(KIND=LDP) MAXSIZEX
	EXTERNAL FAC
!
	INTEGER, PARAMETER :: MAX_N_CONT=20
	INTEGER NUM_CONTOUR_LEVELS
	REAL*4 CONTOUR_LEVELS(MAX_N_CONT)
	REAL*4 MAX_CONT_LEV
	REAL*4 MIN_CONT_LEV
	LOGICAL LOG_CONT_LEVS
!
	CHARACTER*20 TYPE_OF_IMAGE
	REAL*4 TR_VEC(6)
	REAL*4 TR_VEC_SCL(6)
	REAL*4 XPIX_SIZE
	REAL*4 YPIX_SIZE
!
! Vector arrays
!
	REAL*4 LINEXST(MAXVEC),LINEYST(MAXVEC)
	REAL*4 LINEXEND(MAXVEC),LINEYEND(MAXVEC)
	INTEGER NVEC
	LOGICAL VEC,FLAGLINE(MAXVEC),INIT
	LOGICAL QUERYFLAG
	INTEGER VECPEN(MAXVEC)
!
	REAL(KIND=LDP) LCI(256)
	REAL(KIND=LDP) RCI(256)
	REAL(KIND=LDP) BCI(256)
	REAL(KIND=LDP) GCI(256)
	REAL*4 CONTRAST
	REAL*4 BRIGHTNESS
!
! String arrays (not labels or titles)
!
	REAL*4 ORIENTATION(MAXSTR),XSTR(MAXSTR),YSTR(MAXSTR)
	REAL*4 STR_EXP(MAXSTR),LOC_PG(MAXSTR)
	REAL*4 XSTRPOS(MAXSTR),YSTRPOS(MAXSTR)
	CHARACTER*50 STRING(MAXSTR)
	INTEGER NSTR,ISTR
	INTEGER LOC(MAXSTR)
	INTEGER STR_COL(MAXVEC)
	LOGICAL FLAGSTR(MAXSTR),STR
!
	REAL*4 DXST,DXEND,DYST,DYEND
!
! Pen type, Device ID's
!
! LINE_STYLE   : Type of curve (solid, dsahed etc drawn through data)
! MARKER_STYLE : Type of marker drawn at each datum. If negative, no
!                   curve is drawn between the markers.
!
        INTEGER LINE_STYLE(MAX_PLTS)
        INTEGER LINE_WGT(MAX_PLTS)
	INTEGER MARKER_STYLE(MAX_PLTS)
	LOGICAL HARD,TITONRHS,FIRST,FSTOPEN,DASH,MARK,CON,HIST
	INTEGER, SAVE :: IPALLET
!
! E, cursor, and continuum parameters.
!
	REAL*4 X,Y
	REAL*4 EW,CENTROID,CONTIN,DWL,XCUR(2),YCUR(2),SLOPE
	INTEGER CURSERR
	CHARACTER*1 CURSVAL
!
! Functions
!
	LOGICAL END_CURS
	INTEGER PGBEG, PGCURS
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
! Parameters for indicating the size of the plot.
!
	REAL*4 XCM,DXLENGTH,DYLENGTH,ASR,TEMPASR,DASR
!
! Miscellaneous
	CHARACTER*3 ANS
	INTEGER LXL,LYL,LT1,LT2,ISTAT,LENGTH
	REAL*4 V1,IT,RL1,HT1,EXPCHAR,TICK_FAC,EXPMARK
	REAL*4 MARGINX(2),MARGINY(2),MARGINX_SCL(2)
	REAL*4 T1,T2,T3
	REAL*4 XVAL,YVAL
	REAL*4 XVAL_SAV,YVAL_SAV,ZVAL_SAV
!
	INTEGER BEG,Q,NPLTS2
	CHARACTER*8 ITEM
	CHARACTER*50 VALUE,FILE,FILE_CUT
	CHARACTER*10 VALUE_CUT
!
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
!
! Loop variables
	INTEGER I,J,K,L,IOS,CNT
!
! Color variables
	REAL*4 RED(0:15),BLUE(0:15),GREEN(0:15)
	INTEGER PEN_COL(0:MAXPEN)
	REAL*4 RTEMP,BTEMP,GTEMP
!
! Printer variables.
!
	CHARACTER(LEN=80) PRINTER
	CHARACTER(LEN=80), SAVE :: HARD_FILE
	CHARACTER(LEN=20), SAVE :: HARD_TYPE
	CHARACTER(LEN=20), SAVE :: HARD_FMT_STR
	INTEGER, SAVE :: HARD_CNT
!
	INTEGER PGOPEN,PLT_LINE_WGT
	INTEGER ID
	REAL*4 TOTXMM,TOTYMM,SCALEFACY,SCALEFAC
	REAL*4 PRINTX1,PRINTX2,PRINTY1,PRINTY2
	REAL*4 PRINTX1_SCL,PRINTX2_SCL
!
! Variables for options 'WP' (i.e. write plot) and 'RP' (i.e. read plot).
!
	INTEGER IST,IEND,N_REC_SIZE
!
	SAVE RED,BLUE,GREEN
	SAVE FSTOPEN,PEN_COL,DASH,PRINTER
	SAVE MARGINX,MARGINY,MARGINX_SCL
	SAVE PLT_LINE_WGT
	DATA FSTOPEN,DASH,PRINTER/.TRUE.,.FALSE.,'FIRST'/
	DATA PLT_LINE_WGT/1/
	DATA HARD_FMT_STR/'_1'/
	DATA HARD_CNT/1/
	DATA XCM/20/
	DATA FIRST/.TRUE./
	N_REC_SIZE=1000
!
	IF(NPLTS .EQ. 0)THEN
	  WRITE(T_OUT,*)'Error - No calls made to curve'
	  RETURN
	END IF
!
! Define character strings for switching between VT and TEK modes.
!
	NPLTS2=NPLTS
	TO_TEK(1:1)=CHAR(27)
	TO_TEK(2:6)='[?38h'
	TO_VT(1:1)=CHAR(27)
	TO_VT(2:2)=CHAR(3)
!
! Set the Aspect ratio of the plot to the default value. No longer ask
! for a new value, since can be set with N option.
!
	ASR=0.0			!Leave as is
!
	LOG_AXIS=' '
	XLABEL=XLAB
	YLABEL=YLAB
	XLAB_FILE=' '
	YLAB_FILE=' '
	YLABEL_R_AX=' '
	LXL=LEN(XLABEL)
	LYL=LEN(YLABEL)
	STR=.FALSE.
	VEC=.FALSE.
	NORMAL_R_Y_AXIS=.TRUE.
	DO I=1,MAXSTR
	  FLAGSTR(I)=.FALSE.
	  FLAGLINE(I)=.FALSE.
	  STR_EXP(I)=1.0	   !Default (changed using SE option only)
	  STR_COL(I)=1		   !Default (changed using SE option only)
	END DO
	DO_ERROR=.TRUE.
	OPTION=PASSED_OPT
	CALL SET_CASE_UP(OPTION,1,0)
	SCALEPOL=1.0D0
!
	IF(PASSED_OPT .EQ. 'COL')THEN
	  TYPE_OF_IMAGE='COL'
	ELSE IF(PASSED_OPT .EQ. 'VEC')THEN
	  TYPE_OF_IMAGE='VEC'
	ELSE IF(PASSED_OPT .EQ. 'GREY')THEN
	  TYPE_OF_IMAGE='GREY'
	ELSE IF(PASSED_OPT .EQ. 'GRAY')THEN
	  TYPE_OF_IMAGE='GREY'
	ELSE IF(PASSED_OPT .EQ. 'CONT')THEN
 	  TYPE_OF_IMAGE='CONT'
	ELSE
	  TYPE_OF_IMAGE='COL'
	END IF
!
	YPAR_R_AX(1)=0.0D0 ; YPAR_R_AX(2)=0.0D0
!
! Plts connected
!
	DO I=1,NPLTS
	  TYPE_CURVE(I)='L'
	END DO
	MARK=.FALSE.
	WHICH_Y_AX(1:NPLTS)='L'
!
! Define line representations (initially not dashed)
!
	DO I=1,NPLTS
	  LINE_STYLE(I)=1
	  LINE_WGT(I)=1
	END DO
!
! Assign a color index to each pen. Keep previus assignments if they have been
! made.
!
	IF (FSTOPEN) THEN
	  DO I=0,MAXPEN
	    PEN_COL(I)=MIN(I,15)
	  END DO
	  TITLE(1:N_TITLE)=' '
	END IF
	TITLE(1)=TITL
!
! Assign a color index to each vecpen
!
	IF (FSTOPEN) THEN
	  DO I=1,10
	    VECPEN(I)=I+5
	  END DO
	  DO I=11,MAXVEC
	    VECPEN(I)=14
	  END DO
	  DO_SCALE=.TRUE.
	END IF
!
! Define default type for "Data markers"
!
	DO I=1,MIN(12,MAX_PLTS)
	  MARKER_STYLE(I)=I+1           !1 is not appropriate
	END DO
!
	EXPCHAR=1.0   			!Was 1.01 to prevent crashes.
	EXPMARK=1.0   			!Symbol size on plots.
	TICK_FAC=1.0D0
	PLT_LINE_WGT=1
!
	ANS='P'
	IF(ALLOCATED(MAP))THEN
	  DEALLOCATE(MAP)
	  DEALLOCATE(XVEC)
	  DEALLOCATE(YVEC)
	END IF
	ALLOCATE(MAP(XDIM,YDIM))
	ALLOCATE(XVEC(XDIM))
	ALLOCATE(YVEC(YDIM))
!
	MAP=DP_MAP
	THETA=DP_THETA
	XVEC=DP_XVEC
	YVEC=DP_YVEC
!
! Set paparmeter for CONTOUR and GREY scale plots.
!
	XPIX_SIZE=(XVEC(XDIM)-XVEC(1))/(XDIM-1)
	YPIX_SIZE=(YVEC(YDIM)-YVEC(1))/(YDIM-1)
	TR_VEC(1)=XVEC(1)-XPIX_SIZE
	TR_VEC(2)=XPIX_SIZE
	TR_VEC(3)=0.0D0
	TR_VEC(4)=YVEC(1)-YPIX_SIZE
	TR_VEC(5)=0.0D0
	TR_VEC(6)=YPIX_SIZE
!
	T1=1.0/3.0D0
	T2=1.0/(NSCL-1)
	TR_VEC_SCL(1)=0.0-T1
	TR_VEC_SCL(2)=T1
	TR_VEC_SCL(3)=0.0D0
	TR_VEC_SCL(4)=0.0-T2
	TR_VEC_SCL(5)=0.0D0
	TR_VEC_SCL(6)=T2
!
! Look for absica limits
!
	XMIN=XVEC(1)
	XMAX=XVEC(1)
	DO I=1,XDIM
	  IF(XVEC(I) .LT. XMIN)XMIN=XVEC(I)
	  IF(XVEC(I) .GT. XMAX)XMAX=XVEC(I)
	END DO
!
! Look for ordinate limits
	YMIN=YVEC(1)
	YMAX=YVEC(1)
	DO I=1,YDIM
	  IF(YVEC(I) .GT. YMAX) YMAX=YVEC(I)
	  IF(YVEC(I) .LT. YMIN) YMIN=YVEC(I)
	END DO
	IF(ABS(YMAX-YMIN) .LT. 1.0D-08)THEN
	  YMIN=0.0
	  YMAX=1.0
	  WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	END IF
	IF(ABS(XMAX-XMIN) .LT. 1.0D-08)THEN
	  XMIN=0.0
	  XMAX=1.0
	  WRITE(T_OUT,*)'Invalid X limits - setting default values'
	END IF
!
	MAPMAX=MAXVAL(MAP)
	MAPMIN=MINVAL(MAP)
	MIN_CONT_LEV=MAPMIN
	MAX_CONT_LEV=MAPMAX
	LOG_CONT_LEVS=.FALSE.
!
	XPAR(1)=XMIN
	XPAR(2)=XMAX
	YPAR(1)=YMIN
	YPAR(2)=YMAX
	ZPAR(1)=MAPMIN
	ZPAR(2)=MAPMAX
	XPAR_SAV(1)=XPAR(1)		!Indicate value limits evaluated for.
	XPAR_SAV(2)=XPAR(2)
!
! Open user set workstation (default is set into the system).
! We only ask
! for work-station name if first call to GRAMON.
!
	WRITE(T_OUT,4) XMIN,XMAX
4	FORMAT(' Abisca limits  :',1P2E14.4)
	WRITE(T_OUT,11) YMIN,YMAX
11	FORMAT(' Ordinate limits:',1P2E14.4)
	WRITE(T_OUT,'(A)')' '
!
	IF (FSTOPEN) THEN
	  WRITE(T_OUT,*)'NOTE: /XWINDOW supports color'
	  WRITE(T_OUT,*)'    while /XTERM allows auto-switching of the cursor'
!	  BEG = PGBEG(IZERO,'?',IONE,IONE)
	  BEG=PGOPEN('?')
	  CALL PGQID(ID)
!
! Initialize viewport.
!
	  CALL PGASK(.FALSE.)                      !TURNS OFF PAGE PROMPT
	  CALL PGVSTD
	  CALL PGENV(XPAR(1),XPAR(2),YPAR(1),YPAR(2),IZERO,IZERO)
	  CALL PGQVP(IZERO,MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	  IF (MARGINX(1) .LT. (.9)) MARGINX(1)=MARGINX(1)+.05
!
!
! Increase default border size.
!
          MARGINX(1)=0.15; MARGINX(2)=0.9
          MARGINX_SCL(1)=0.96; MARGINX_SCL(2)=0.98
          MARGINY(1)=0.15; MARGINY(2)=0.9
!
! Set preferred defaults for pen colors.
!
	  CALL PGSCR(0,.7,.7,.7)     !set color representations
	  CALL PGSCR(2,1.0,0.0,0.0)
	  CALL PGSCR(3,0.0,0.0,1.0)
	  CALL PGSCR(4,0.0,0.6,0.3)
	  CALL PGSCR(5,.5,0.0,.7)
	  CALL PGSCR(13,0.0,1.0,1.0)
	  CALL PGSCR(15,.95,.95,.95)
          DO I=0,15                  !Get these + def color representations.
            CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
	END IF
	FSTOPEN=.FALSE.
	HARD=.FALSE.
1000	ANS='P'

!
!	WRITE(T_OUT,*)TO_VT
!
	WRITE(T_OUT,*)' H,P,Z(ZN),A(F),L,CC,CP,N, M,W,D,C,B,'//
	1              ' VC,VF,VE, SC,SF,SE, E'
	CALL GEN_IN(ANS,'ANS')
	L=LEN_TRIM(ANS)
	CALL SET_CASE_UP(ANS,IONE,IZERO)
        IF(ANS .EQ. 'H')THEN
	  WRITE(T_OUT,*)'P=Plot - default'
          WRITE(T_OUT,*)'E=EXIT from PLOT package'
	  WRITE(T_OUT,*)'Z=Hardcopy (ZN=Asks for new hard device)'
          WRITE(T_OUT,*)'A=Define Axis Parameters'
          WRITE(T_OUT,*)'F=Change default axis parameters'
	  WRITE(T_OUT,*)'L=Modify Axis Labels and Titles'
 	  WRITE(T_OUT,*)'D=Switch dashed lines on/off'
	  WRITE(T_OUT,*)'M=Switch marking data points on/off'
	  WRITE(T_OUT,*)'W=Change thickness of curves'
	  WRITE(T_OUT,*)'C=Indicate how curves are to be connected (L,H,A)'
	  WRITE(T_OUT,*)'B=Switch error bars on/off'
	  WRITE(T_OUT,*)'NC=Set number and placement of contour levels'
	  WRITE(T_OUT,*)'IM=Set image type'
	  WRITE(T_OUT,*)'SM=Smooth image'
	  WRITE(T_OUT,*)'CC=Change Color setting'
	  WRITE(T_OUT,*)'CP=Change Pen (Color Index)'
          WRITE(T_OUT,*)'N=Define IDX,HT etc'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)'IM=Set image display mode'
	  WRITE(T_OUT,*)'FL=Flip map in X'
	  WRITE(T_OUT,*)'RO=Rotate map by -90"'
	  WRITE(T_OUT,*)'LOG=Take log of image'
c
	  WRITE(T_OUT,*)'LX=Switch between LINEAR/LOG X axis'
	  WRITE(T_OUT,*)'LY=Switch between LINEAR/LOG Y axis'
	  WRITE(T_OUT,*)'LXY=Switch between LINEAR/LOG for X and Y axes'
	  WRITE(T_OUT,*)'VC=Define line vectors using cursor'
	  WRITE(T_OUT,*)'VF=Define line vectors using file input'
	  WRITE(T_OUT,*)'VF=Online edit of vectors'
          WRITE(T_OUT,*)'SC=Define strings using cursor'
          WRITE(T_OUT,*)'SF=Define strings using file input'
          WRITE(T_OUT,*)'SE=Online edit of strings'
	  WRITE(T_OUT,*)'CL=Clear Graphics Screen'
	  WRITE(T_OUT,*)'H=Help'
          GOTO 1000
	ELSE IF(ANS .EQ. 'E')THEN
          IF(STR .OR. VEC)THEN
	    FILNAME='PGOUT'
	    CALL GEN_IN(FILNAME,'Filname for STRING/VEC storage (no ext)')
	    IF(FILNAME .EQ. ' ')FILNAME='JNK'	!Force output to jnk file as
	    L=LEN_TRIM(FILNAME)			!precaution.
	    IF(FILNAME .NE. ' ')THEN
	      IF(STR)THEN
	        FILNAME=FILNAME(:L)//'.STR'
                L=L+4
                OPEN(UNIT=33,FILE=FILNAME(:L),STATUS='NEW')
	          DO I=1,MAXSTR
                    IF(FLAGSTR(I))THEN
                      WRITE(33,17)LOC(I),XSTR(I),YSTR(I),ORIENTATION(I),
	1                   STR_EXP(I),TRIM(STRING(I))
17	              FORMAT(1X,I1,',',4(F9.4,','),1H',A,1H')
	           END IF
	         END DO
                CLOSE(UNIT=33)
	      END IF
	      IF(VEC)THEN
	        IF(STR)THEN
	          FILNAME(L-3:L)='.VEC'
	        ELSE
	          FILNAME=FILNAME(:L)//'.VEC'
                  L=L+4
	        END IF
                OPEN(UNIT=33,FILE=FILNAME(:L),STATUS='NEW')
	          DO I=1,MAXVEC
                    IF(FLAGLINE(I))THEN
	              WRITE(33,18)LINEXST(I),LINEYST(I),
	1                         LINEXEND(I),LINEYEND(I)
18	              FORMAT(1X,1P,4E18.8)
	            END IF
	          END DO
                CLOSE(UNIT=33)
	      END IF		!VEC end
            END IF		!Filename check
          END IF		!STR or VEC present
	  RETURN
	END IF			!ANS = 'E'
!
	IF(FIRST)THEN
! Define xinc
	  XINC=SPACING(XPAR(1),XPAR(2))
	  YINC=SPACING(YPAR(1),YPAR(2))
	  IXTICK=2
	  IYTICK=2
! Compute Xmin,Xmax,Ymin,Ymax if not input
	  V1=1000
	  XPAR(1)=ABS(XINC)*(AINT(XMIN/ABS(XINC)+V1)-V1)
	  XPAR(2)=ABS(XINC)*(AINT(XMAX/ABS(XINC)-V1)+V1)
	  YPAR(1)=ABS(YINC)*(AINT(YMIN/ABS(YINC)+V1)-V1)
	  YPAR(2)=ABS(YINC)*(AINT(YMAX/ABS(YINC)-V1)+V1)
!
	  XNUMST=XPAR(1)
	  YNUMST=YPAR(1)
!
! Determine number of decimals
	  IDX=NDEC(XINC)
	  IT=NDEC(XPAR(1))
	  IF(IT .GT. IDX)IDX=IT
	  IDY=NDEC(YINC)
	  IT=NDEC(YPAR(1))
	  IF(IT .GT. IDY)IDY=IT
!
	  IPALLET=1
	  FIRST=.FALSE.
	END IF
!
	IF(ANS .EQ. 'A' .OR. ANS .EQ. 'F')THEN
	  WRITE(T_OUT,4) XMIN,XMAX
	  CALL GEN_IN(XPAR,I,ITWO,'XST,XEND')
!
! Look for ordinate limits in X range. This will now generate a plot with
! roughly the correct scaling even though the oridinates values may be vastly
! different outside the plot window.
!
	  IF(XPAR(2) .NE. XPAR_SAV(2) .OR. XPAR(1) .NE. XPAR_SAV(1))THEN
	    YMIN=1.0E+32
	    YMAX=-1.0E+32
	    DO I=1,XDIM
	      IF( (XVEC(I) .GE. XPAR(1) .AND. XVEC(I) .LE. XPAR(2)) .OR.
	1         (XVEC(I) .GE. XPAR(2) .AND. XVEC(I) .LE. XPAR(1)) )THEN
	        DO J=1,YDIM
	          YMAX=MAX(YMAX,YVEC(J))
	          YMIN=MIN(YMIN,YVEC(J))
	        END DO
	      END IF
	    END DO
	    IF(ABS(YMAX-YMIN) .LT. 1.0D-08 .OR.
	1        YMIN .EQ. 1.0E+32 .OR. YMAX .EQ. -1.0E+32)THEN
	      YMIN=0.0
	      YMAX=1.0
	    END IF
	    YINC=SPACING(YMIN,YMAX)
	    IYTICK=2
	    V1=1000
	    YPAR(1)=ABS(YINC)*(AINT(YMIN/ABS(YINC)+V1)-V1)
	    YPAR(2)=ABS(YINC)*(AINT(YMAX/ABS(YINC)-V1)+V1)
	    XPAR_SAV(1)=XPAR(1)		!Indicate value limits evaluated for.
	    XPAR_SAV(2)=XPAR(2)
	  END IF
!
	  WRITE(T_OUT,11)YMIN,YMAX
	  CALL GEN_IN(YPAR,I,ITWO,'YST,YEND')
!
! Check that limits inserted are not absurd.
!
	  IF(ABS(YPAR(2)-YPAR(1)) .LT. 1.0E-08)THEN
	    YPAR(1)=0
	    YPAR(2)=1.00
	    WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	  END IF
	  IF(ABS(XPAR(2)-XPAR(1)) .LT. 1.0E-08)THEN
	    XPAR(1)=0
	    XPAR(2)=1.00
	    WRITE(T_OUT,*)'Invalid X limits - setting default values'
	  END IF
!
	  WRITE(T_OUT,11)MAPMIN,MAPMAX
	  CALL GEN_IN(ZPAR,I,ITWO,'MAP_ST,MAP_END')
!
	  XNUMST=XPAR(1)
	  YNUMST=YPAR(1)
!
! Set defaults for XINC and YINC. IXTICK-1 is the number of tick
! marks between numbered ticks.
!
	  XINC=SPACING(XPAR(1),XPAR(2))
	  IXTICK=2
	  YINC=SPACING(YPAR(1),YPAR(2))
	  IYTICK=2
	  WRITE(T_OUT,*)'Indicate spacing between numberd tickmarks'
	  CALL GEN_IN(XINC,'For X tickmarks')
	  CALL GEN_IN(YINC,'For Y tickmarks')
!
! Determine number of decimals
!
	  IDX=NDEC(XINC)
	  IT=NDEC(XPAR(1))
	  IF(IT .GT. IDX)IDX=IT
	  IDY=NDEC(YINC)
	  IT=NDEC(YPAR(1))
	  IF(IT .GT. IDY)IDY=IT
!
	  IF(ANS .EQ. 'A')GOTO 1000
!
! This section is for additional axis fiddling. Can change number
! of digits after decimal point, and numbers at which axis labeling
! begins. These parameters are restored to their original values
! everytime 'A' option is called.
!
	  CALL GEN_IN(IXTICK,'X minor divisions')
	  IXTICK=MAX(IXTICK,1)
	  CALL GEN_IN(IYTICK,'Y minor divisions')
	  IYTICK=MAX(IYTICK,1)
	  CALL GEN_IN(IDX,'# of X Dec digits')
	  CALL GEN_IN(IDY,'# of Y Dec digits')
	  CALL GEN_IN(XNUMST,'Number beginning X Axis')
	  CALL GEN_IN(YNUMST,'Number beginning Y Axis')
!
	  FIRST=.FALSE.
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'L')THEN
	  CALL GEN_IN(XLABEL,'XLAB')
	  CALL GEN_IN(YLABEL,'YLAB')
	  DO J=1,N_TITLE
	    CALL GEN_IN(TITLE(J),'TITLE')
	    IF(TITLE(J) .EQ. ' ')THEN
	      TITLE(J+1:N_TITLE)=' '
	      EXIT
	    END IF
	  END DO
	  CALL GEN_IN(TITONRHS,'TITONRHS')
	  GOTO 1000
!
! 
	ELSE IF(ANS .EQ. 'RO')THEN
	  ALLOCATE (TMP_MAP(YDIM,XDIM))
	  DO I=1,XDIM
	    DO J=1,YDIM
	       TMP_MAP(YDIM+1-J,I)=MAP(I,J)
	    END DO
	  END DO
!
	  ALLOCATE(TMP_VEC(XDIM))
	  TMP_VEC(1:XDIM)=XVEC(1:XDIM)
	  DEALLOCATE(XVEC) ; ALLOCATE(XVEC(YDIM))
	  XVEC(1:YDIM)=YVEC(1:YDIM)
	  DEALLOCATE(YVEC) ;ALLOCATE(YVEC(XDIM))
	  YVEC(1:XDIM)=TMP_VEC(1:XDIM)
!
	  I=XDIM ; XDIM=YDIM ;	  YDIM=I
	  DEALLOCATE(MAP)
	  ALLOCATE(MAP(XDIM,YDIM))
	  MAP=TMP_MAP
	  DEALLOCATE(TMP_MAP)
	  DEALLOCATE(TMP_VEC)
	ELSE IF(ANS .EQ. 'FL')THEN
	  DO J=1,YDIM
	    DO I=1,XDIM/2
	      T1=MAP(I,J)
	      MAP(I,J)=MAP(XDIM+1-I,J)
	      MAP(XDIM+1-I,J)=T1
	    END DO
	  END DO
!
	  DO I=1,XDIM/2
	    T1=XVEC(I)
	    XVEC(I)=-XVEC(XDIM+1-I)
	    XVEC(XDIM+1-I)=-T1
	  END DO
	  T1=XPAR(1)
	  XPAR(1)=-XPAR(2)
	  XPAR(2)=-T1
!
	ELSE IF(ANS .EQ. 'LOG')THEN
	  T1=MINVAL(MAP, MASK=MAP .GT. 0.0)
	  WHERE (MAP .GE. T1)
	    MAP=LOG10(MAP)
	  ELSEWHERE
	    MAP=LOG10(T1)-1.0
	  END WHERE
	  MAPMAX=MAXVAL(MAP)
	  MAPMIN=MINVAL(MAP)
	  WRITE(6,*)'Maximum value in map is',MAPMAX
	  WRITE(6,*)'Useful ninimum value in map is',MAPMIN+1.0D0
!
	ELSE IF(ANS .EQ. 'NC')THEN
	  WRITE(T_OUT,'(A,1P2E14.4)')' Map limits: ',MAPMIN,MAPMAX
	  CALL GEN_IN(NUM_CONTOUR_LEVELS,'Number of contour levels')
	  CALL GEN_IN(LOG_CONT_LEVS,'Logarithmically spaced contours?')
	  CALL GEN_IN(MIN_CONT_LEV,'Minimum contour level')
	  CALL GEN_IN(MAX_CONT_LEV,'Maximum contour level')
	  IF(LOG_CONT_LEVS)THEN
	    IF(MIN_CONT_LEV .LE. 0)MIN_CONT_LEV=MAX_CONT_LEV/1.0E+05
	    T1=ALOG10(MAX_CONT_LEV/MIN_CONT_LEV)/(NUM_CONTOUR_LEVELS-1)
	    CONTOUR_LEVELS(1)=MIN_CONT_LEV
	    DO I=2,NUM_CONTOUR_LEVELS
	      CONTOUR_LEVELS(I)=CONTOUR_LEVELS(I-1)*T1
	    END DO
	  ELSE
	    T1=(MAX_CONT_LEV-MIN_CONT_LEV)/(NUM_CONTOUR_LEVELS-1)
	    DO I=1,NUM_CONTOUR_LEVELS
	      CONTOUR_LEVELS(I)=T1*(I-1)
	    END DO
	  END IF
!
	ELSE IF(ANS .EQ. 'IM')THEN
	  CALL GEN_IN(TYPE_OF_IMAGE,'Image type: VEC, GREY, COL, CONT')
	  CALL SET_CASE_UP(TYPE_OF_IMAGE,IZERO,IZERO)
	  IF(TYPE_OF_IMAGE .EQ. 'GRAY')TYPE_OF_IMAGE='GREY'
	  IF(TYPE_OF_IMAGE .NE. 'VEC' .AND.
	1           TYPE_OF_IMAGE .NE. 'GREY' .AND.
	1           TYPE_OF_IMAGE .NE. 'COL' .AND.
	1           TYPE_OF_IMAGE .NE. 'CONT')THEN
	     WRITE(T_OUT,*)'Invalid image type'
	     GOTO 1000
	  END IF
	ELSE IF(ANS .EQ. 'SM')THEN
	  HAN=.TRUE.
	  XHAN=3
	  CALL GEN_IN(XHAN,'# of points to HAN over in X direction')
	  YHAN=XHAN
	  CALL GEN_IN(YHAN,'# of points to HAN over in Y direction')
	  IF(XHAN .EQ. 1 .AND. YHAN .EQ. 1)THEN
	    HAN=.FALSE.
	    GOTO 1000
	  END IF
!
! Do smoothing in Y direction. It is assumed that data outside
! array is zero (reasonable approximation if spatial array
! is bigger than image).
!
	  T1=FAC(XHAN-1)/(2.0**(XHAN-1))
	  DO K=0,XHAN/2
	   SM_ARRAY(K,0)=T1/FAC(XHAN-1-K)/FAC(K)
	  END DO
	  DO J=1,YHAN/2
	   SM_ARRAY(K,J)=SM_ARRAY(K,0)
	  END DO
!
	  T1=FAC(YHAN-1)/(2.0**(YHAN-1))
	  DO K=0,XHAN/2
	    DO J=0,YHAN/2
	      SM_ARRAY(K,J)=SM_ARRAY(K,J)/T1/FAC(YHAN-1-J)/FAC(J)
	      SM_ARRAY(-K,-J)=SM_ARRAY(K,J)
	    END DO
	  END DO
	ELSE IF(ANS .EQ. 'CC')THEN
	  CALL CHANGE_COLOR(RED,BLUE,GREEN)
	  GOTO 1000
	ELSE IF(ANS .EQ. 'CP')THEN
	  CALL CHANGE_PEN(PEN_COL,MAXPEN,NPLTS2)
	  GOTO 1000
	ELSE IF(ANS .EQ. 'N')THEN
!
	  CALL GEN_IN(EXPCHAR,'Expand Char.')
	  TICK_FAC=EXPCHAR
	  CALL GEN_IN(TICK_FAC,'Expand TICK')
!
! Define the Aspect ratio of the plot.
!
	  WRITE(T_OUT,*)'Y/X > 0, 0(Device default) , X/Y < 0'
	  CALL GEN_IN(ASR,'Aspect Ratio')
!
700	  CALL GEN_IN(MARGINX,I,ITWO,'Left and Right Margins (0:1)')
	  IF((MARGINX(1) .GT. 1) .OR. (MARGINX(1) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Left-Margin is incorrect, try again.'
	    GOTO 700
	  END IF
	  IF((MARGINX(2) .GT. 1) .OR. (MARGINX(2).LT. 0)) THEN
	    WRITE(T_OUT,*)'Right-Margin is incorrect, try again.'
	    GOTO 700
	  END IF
!
710	  CALL GEN_IN(MARGINY,I,ITWO,'Bottom and Top Margins (0:1)')
	  IF((MARGINY(1) .GT. 1.0) .OR. (MARGINY(1) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Bottom-Margin is incorrect, try again.'
	    GOTO 710
	  END IF
	  IF((MARGINY(2) .GT. 1.0) .OR. (MARGINY(2) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Top-Margin is incorrect, try again.'
	    GOTO 710
	  END IF
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'W')THEN
	    CALL GEN_IN(LINE_WGT,I,NPLTS,
	1                'Line Weight:: 1,2 etc')
	ELSE IF( ANS .EQ. 'M')THEN
	  IF(MARK)THEN
            WRITE(T_OUT,*)'Will not mark data points'
	    MARK=.FALSE.
	  ELSE
	    MARK=.TRUE.
            WRITE(T_OUT,*)'Wil now mark data points'
	    CALL GEN_IN(EXPMARK,'Expand Symbol Size')
	    CALL GEN_IN(MARKER_STYLE,I,NPLTS,
	1                'Marker style {+/-}1,...,31)')
	  END IF
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'C')THEN
	  DO I=1,NPLTS
850	    WRITE(T_OUT,'(I3,'' : '',$)')I
	    CALL GEN_IN(TYPE_CURVE(I),'TC=')
	    CALL SET_CASE_UP(TYPE_CURVE(I),1,1)
	    IF( TYPE_CURVE(I) .NE. ' ' .AND.
	1       TYPE_CURVE(I) .NE. 'L' .AND.              !Normal line
	1       TYPE_CURVE(I) .NE. 'B' .AND.              !Broken
	1       TYPE_CURVE(I) .NE. 'I' .AND.              !Invisible
	1       TYPE_CURVE(I) .NE. 'V' .AND.              !Verticle lines
	1       TYPE_CURVE(I) .NE. 'A' .AND.              !Hist - X vert
	1       TYPE_CURVE(I) .NE. 'H' )THEN              !Histogram
	      WRITE(T_OUT,*)'Invalid connection specifier (L,H,A,I OR B)'
	      GOTO 850
	    END IF
	  END DO
	  GOTO 1000
!	
 	ELSE IF( ANS .EQ. 'D')THEN
	  IF(DASH)THEN
            WRITE(T_OUT,*)'Now all solid line plots '
	    DO I=1,NPLTS
	      LINE_STYLE(I)=1
	    END DO
	    DASH=.FALSE.
	  ELSE
	    WRITE(T_OUT,*)'Input dashed pens for ',NPLTS
	    DO I=1,NPLTS,5
	      DO J=1,5
	        IF(I+J-1 .LE. NPLTS)LINE_STYLE(I+J-1)=J
	      END DO
	    END DO
	    CALL GEN_IN(LINE_STYLE,I,NPLTS,
	1      'Dashed pen types (1,...,5)')
	    DASH=.TRUE.
	  END IF
	  GOTO 1000
!
! Edit each dashed pen separately.
!
 	ELSE IF( ANS .EQ. 'DE')THEN
	  DO J=1,NPLTS
	    CALL GEN_IN(LINE_STYLE(J),'Dashed pen types (1,...,5)')
	  END DO
	  DASH=.TRUE.
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'B')THEN
	  IF(DO_ERROR)THEN
            WRITE(T_OUT,*)'Swithcing off error bar plotting'
	    DO_ERROR=.FALSE.
	  ELSE
	    DO_ERROR=.TRUE.
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'LXY')THEN
	  IF(LOG_AXIS(1:3) .EQ. 'LOG')THEN
            WRITE(T_OUT,*)'Swithcing off ALL LOG axes'
	    LOG_AXIS=' '
	  ELSE
	    LOG_AXIS='LOGXY'
	  END IF
	  GOTO 1000
	ELSE IF( ANS .EQ. 'LX')THEN
	  IF(LOG_AXIS(1:5) .EQ. 'LOGX')THEN
            WRITE(T_OUT,*)'Swithcing off X LOG axes'
	    LOG_AXIS=' '
	  ELSE IF(LOG_AXIS(1:5) .EQ. 'LOGXY')THEN
            WRITE(T_OUT,*)'Swithcing off X LOG axes'
	    LOG_AXIS='LOGY'
	  ELSE IF(LOG_AXIS(1:4) .EQ. 'LOGY')THEN
	    LOG_AXIS='LOGXY'
	  ELSE
	    LOG_AXIS='LOGX'
	  END IF
	  GOTO 1000
	ELSE IF( ANS .EQ. 'LY')THEN
	  IF(LOG_AXIS(1:4) .EQ. 'LOGY')THEN
            WRITE(T_OUT,*)'Swithcing off Y LOG axes'
	    LOG_AXIS=' '
	  ELSE IF(LOG_AXIS(1:5) .EQ. 'LOGXY')THEN
            WRITE(T_OUT,*)'Swithcing off Y LOG axes'
	    LOG_AXIS='LOGX'
	  ELSE IF(LOG_AXIS(1:4) .EQ. 'LOGX')THEN
	    LOG_AXIS='LOGXY'
	  ELSE
	    LOG_AXIS='LOGY'
	  END IF
	  GOTO 1000
!
!
! String handling section.
!
	ELSE IF (ANS .EQ. 'SF')THEN
	  INIT=.FALSE.          !Strings automatically initialized first time.
	  IF(STR)CALL GEN_IN(INIT,'Initialize string list ?')
          STR=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXSTR
	      FLAGSTR(I)=.FALSE.
	      STR_EXP(I)=1.0	  !Default (changed at editing stage only)
	    END DO
	  END IF
	  ISTR=0
923	  CALL GEN_IN(FILNAME,'Input string file to be read')
	  CALL SET_CASE_UP(FILNAME,1,0)
	  L=LEN_TRIM(FILNAME)
	  IF(L .NE. 0)THEN
	    IF(INDEX(FILNAME,'.') .EQ. 0)FILNAME=FILNAME(:L)//'.STR'
            OPEN(UNIT=33,FILE=TRIM(FILNAME),STATUS='OLD',ERR=923)
	    L=0
            DO ISTR=1,MAXSTR
	      IF( .NOT. FLAGSTR(ISTR))THEN
	        L=L+1
	        READ(33,*,END=930,ERR=925)LOC(ISTR),XSTR(ISTR),YSTR(ISTR),
	1                   ORIENTATION(ISTR),STR_EXP(ISTR),STRING(ISTR)
	        FLAGSTR(ISTR)=.TRUE.
                GOTO 928
925             WRITE(T_OUT,*)' Error reading string',ISTR
928	        CONTINUE
	      END IF
	      IF(L .NE. 0)
	1             WRITE(T_OUT,932)L,'STRINGS read in from file'
            END DO
930	    CLOSE(UNIT=33)
          END IF
932	  FORMAT(1X,I2,2X,(A))
!
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SE')THEN
!
! First output strings to user, before doing online editing.
!
	  DO I=1,MAXSTR
	    IF(FLAGSTR(I))THEN
	      WRITE(T_OUT,955)I,TRIM(STRING(I))
955	      FORMAT(1X,I2,3X,A)
	    END IF
	  END DO
!
! Allow online editing of strings. GEN_IN is used so user can
! see current value.
!
	  DO WHILE (0 .EQ. 0)
	    ISTR=0			!Default terminates input
	    CALL GEN_IN(ISTR,'String number - 0 exits')
	    IF(ISTR .LE. 0 .OR. ISTR .GT. MAXSTR)GOTO 1000
	    CALL GEN_IN(FLAGSTR(ISTR),'FLAG')
	    IF(FLAGSTR(ISTR))THEN
	      CALL GEN_IN(LOC(ISTR),'LOC')
	      CALL GEN_IN(XSTR(ISTR),'XPOS')
	      CALL GEN_IN(YSTR(ISTR),'YPOS')
	      CALL GEN_IN(ORIENTATION(ISTR),'ORI')
	      CALL GEN_IN(STR_EXP(ISTR),'EXP')
	      CALL GEN_IN(STR_COL(ISTR),'COL')
	      CALL GEN_IN(STRING(ISTR),'STRING')
	    END IF
	  END DO
	  STR=.TRUE.
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'SC')THEN
	  INIT=.FALSE.          !Strings automatically initialized first time.
	  IF(STR)CALL GEN_IN(INIT,'Initialize string list ?')
	  STR=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXSTR
	      FLAGSTR(I)=.FALSE.
	    END DO
	  END IF
	      WRITE(T_OUT,820)
 820	      FORMAT(1X,'LOC: 1-9, 1=lef bot, 2=cent bot, 3=rit bot')
	      WRITE(T_OUT,821)
 821	      FORMAT(1X,'          4=lef mid, 5=cent mid, 6=rit mid')
	      WRITE(T_OUT,822)
 822	      FORMAT(1X,'          7=lef top, 8=cent top, 9=rit top')
	      WRITE(T_OUT,840)
 840	      FORMAT(1X,'ORI= AN ANGLE FROM 0.0 TO 360')
	  ISTR=1
	  DO WHILE(ISTR .LE. MAXSTR)
	    IF(.NOT. FLAGSTR(ISTR))THEN
	      CURSERR = PGCURS(XSTR(ISTR),YSTR(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
830	      WRITE(T_OUT,810,ADVANCE='NO')
810	      FORMAT(1X,'Desc(LOC,ORI,''STRING''): ')
	      READ(T_IN,'(A)',ERR=830)WK_STR
	      READ(WK_STR,*,ERR=830)LOC(ISTR)
	      IF(LOC(ISTR) .EQ. 0)GOTO 1000
	      READ(WK_STR,*,ERR=830)LOC(ISTR),ORIENTATION(ISTR),STRING(ISTR)
	      FLAGSTR(ISTR)=.TRUE.
	    END IF
	    ISTR=ISTR+1
	  END DO
	  WRITE(T_OUT,*)'Maximum number of strings is',MAXSTR
	  GOTO 1000
!
!
	ELSE IF(ANS .EQ. 'LOC')THEN
	  CALL PGQCS(4,T1,T3)			!T3=Character Height
	  T3=T3*1.5
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  IF(END_CURS(CURSVAL))GOTO 70
!
	  I=NINT( (XVAL-XVEC(1))/XPIX_SIZE )+1
	  XVAL=XVEC(I)
	  J=NINT( (YVAL-YVEC(1))/YPIX_SIZE )+1
	  YVAL=YVEC(J)
!
	  T1=XPAR(1)-(XPAR(2)-XPAR(1))/15.0
	  T2=YPAR(2)+(YPAR(2)-YPAR(1))/15.0
	  CALL MON_NUM(T1,T2,XVAL,1,IDX+2)
	  CALL MON_NUM(T1,T2-T3,YVAL,1,IDY+2)
	  CALL MON_NUM(T1,T2-2*T3,MAP(I,J),1,IDY+4)
!
	  DO WHILE(1 .EQ. 1)
	    XVAL_SAV=XVAL
	    YVAL_SAV=YVAL
	    ZVAL_SAV=MAP(I,J)
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    CALL PGSCI(0)     ! Background
	    CALL MON_NUM(T1,T2,XVAL_SAV,1,IDX+2)
	    CALL MON_NUM(T1,T2-T3,YVAL_SAV,1,IDY+2)
	    CALL MON_NUM(T1,T2-2*T3,ZVAL_SAV,1,IDY+4)
	    CALL PGSCI(1)     ! Axes
	    IF(END_CURS(CURSVAL))GOTO 70
!
	    I=NINT( (XVAL-XVEC(1))/XPIX_SIZE )+1
	    XVAL=XVEC(I)
	    J=NINT( (YVAL-YVEC(1))/YPIX_SIZE )+1
	    YVAL=YVEC(J)
	    CALL MON_NUM(T1,T2,XVAL,1,IDX+2)
	    CALL MON_NUM(T1,T2-T3,YVAL,1,IDY+2)
	    CALL MON_NUM(T1,T2-2*T3,MAP(I,J),1,IDY+4)
	  END DO
70	  CONTINUE
!
!
! Set Line vectors
!
	ELSE IF(ANS .EQ. 'VF')THEN
	  INIT=.FALSE.
	  IF(VEC)CALL GEN_IN(INIT,'Initialize line drawing list ?')
	  VEC=.TRUE.
	  IF(INIT)THEN
	    NVEC=0
	    DO I=1,MAXVEC
	      FLAGLINE(I)=.FALSE.
	    END DO
	  END IF
!
1100	  CALL GEN_IN(FILNAME,'Line vector file forinput')
	  CALL SET_CASE_UP(FILNAME,1,0)
	  L=LEN_TRIM(FILNAME)
	  ISTR=1
	  IF(L .NE. 0)THEN
	    IF(INDEX(FILNAME,'.') .EQ. 0)FILNAME=FILNAME(:L)//'.VEC'
            OPEN(UNIT=33,FILE=TRIM(FILNAME),STATUS='OLD',ERR=1100)
	    DO ISTR=1,MAXVEC
	      READ(33,*,ERR=1100)LINEXST(ISTR),LINEYST(ISTR),
	1              LINEXEND(ISTR),LINEYEND(ISTR)
	      FLAGLINE(ISTR)=.TRUE.
	    END DO
	    CLOSE(UNIT=33)
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'VE')THEN
!
! First output vectors to user, before doing online editing.
! GEN_IN is used so user can see current value. This section
! also used to input vectors by hand.
!
	  DO I=1,MAXVEC
	    IF(FLAGLINE(I))THEN
	      WRITE(T_OUT,'(1X,I2,4(3X,E14.6))')I,
	1                        LINEXST(I),LINEXEND(I),
	1                        LINEYST(I),LINEYEND(I)
	    END IF
	  END DO
!
	  DO WHILE (0 .EQ. 0)
	    ISTR=0			!Default terminates input
	    CALL GEN_IN(ISTR,'Vector number - 0 exits')
	    IF(ISTR .LE. 0 .OR. ISTR .GT. MAXVEC)GOTO 1150
	    CALL GEN_IN(FLAGLINE(ISTR),'FLAG')
	    IF(FLAGLINE(ISTR))THEN
	      CALL GEN_IN(LINEXST(ISTR),'XST')
	      CALL GEN_IN(LINEYST(ISTR),'YST')
	      CALL GEN_IN(LINEXEND(ISTR),'XEND')
	      CALL GEN_IN(LINEYEND(ISTR),'YEND')
	    END IF
	  END DO
 1150	  VEC=.TRUE.
	  QUERYFLAG=.FALSE.
	  CALL GEN_IN(QUERYFLAG,'Change vector pens?')
	  CALL VECTORPEN(VECPEN,MAXVEC,FLAGLINE)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'VC')THEN
	  INIT=.FALSE.
	  IF(VEC)CALL GEN_IN(INIT,'Initialize line drawing list?')
	  VEC=.TRUE.
	  IF(INIT)THEN
	    NVEC=0
	    DO I=1,MAXVEC
	      FLAGLINE(I)=.FALSE.
	    END DO
	  END IF
	  ISTR=1
	  DO WHILE (ISTR .LT. MAXVEC)
	    IF( .NOT. FLAGLINE(ISTR) )THEN
	      CURSERR = PGCURS(LINEXST(ISTR),LINEYST(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
	      CURSERR = PGCURS(LINEXEND(ISTR),LINEYEND(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
	      IF(LINEXST(ISTR) .EQ. LINEXEND(ISTR) .AND.
	1        LINEYST(ISTR) .EQ. LINEYEND(ISTR) )THEN
	         ISTR=1000
	      ELSE
	        FLAGLINE(ISTR)=.TRUE.
	      END IF
	    END IF
	    ISTR=ISTR+1
	  END DO
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'Z' .OR. ANS .EQ. 'ZN')THEN
	  IF (PRINTER .EQ. 'FIRST' .OR. ANS .EQ. 'ZN') THEN
	    WRITE(T_OUT,*)'Choose a post-script device and file',
	1              ' for printing [file/dev] '
	    BEG=PGOPEN('?')
	  ELSE
	    CALL GEN_IN(HARD_FILE,'Plot file')
	    PRINTER=TRIM(HARD_FILE)//'/'//HARD_TYPE
	    BEG=PGOPEN(PRINTER)
	  END IF
!
! Save hard device and set default file name for net plot.
!
	  CALL PGQINF('TYPE',HARD_TYPE,LENGTH)
	  CALL PGQINF('FILE',HARD_FILE,LENGTH)
!
! Want last occurrence of HARD_FMT_STR in filename to avoid possible
! occurences of HARD_FMT_STR in directory names.
!
	  K=-10
	  J=0
	  DO WHILE(J .NE. K)
	   K=J
	   J=J+INDEX(HARD_FILE(J+1:),TRIM(HARD_FMT_STR))
	  END DO
	  IF(J .NE. 0)HARD_FILE=HARD_FILE(1:J-1)
	  J=INDEX(HARD_FILE,'.ps')
	  IF(J .NE. 0)HARD_FILE=HARD_FILE(1:J-1)
!
	  HARD_CNT=HARD_CNT+1
	  WRITE(HARD_FMT_STR,'(I10)')HARD_CNT
	  DO WHILE(HARD_FMT_STR(1:1) .EQ. ' ')
	    HARD_FMT_STR(1:)=HARD_FMT_STR(2:)
	  END DO
	  HARD_FMT_STR='_'//TRIM(HARD_FMT_STR)
	  HARD_FILE=TRIM(HARD_FILE)//TRIM(HARD_FMT_STR)//'.ps'
	  PRINTER=TRIM(HARD_FILE)//'/'//HARD_TYPE
!
	  HARD=.TRUE.
!
1200	  CALL GEN_IN(PLT_LINE_WGT,'line weight')
	  IF (PLT_LINE_WGT .GT. 201) THEN
	    WRITE(T_OUT,*)'value too large'
	    GOTO 1200
	  END IF
	  GOTO 5000
!
!
	ELSE IF(ANS .EQ. 'CL')THEN
	  IF(.NOT. HARD)THEN		 	!Terminal
	    CALL PGERAS
	    CALL PGETXT
	    CALL PGPAGE
	    READ(T_IN,'(A)')ANS
	  END IF
	  GOTO 1000
!
	ELSE IF (ANS .NE. 'P' .AND. ANS .NE. 'R')THEN
	  WRITE(T_OUT,*)' Invalid option - try again'
	  GOTO 1000
	END IF
!
! Actual plotting section. All parameters except the plot size have been set.
!
5000	CONTINUE
!
! Open the appropriate plotting device
!
 6000	IF (HARD) CALL PGBBUF
!
! Initialize plotting parameters to the correct values
!
	CALL PGENV(XPAR(1),XPAR(2),YPAR(1),YPAR(2),0,-1)
	CALL PGSCH(EXPCHAR)
	IF(HARD)THEN
	  CALL PGSCR(0,1.0,1.0,1.0)  !change background to white
          DO I=1,15                  !Make sure color are set for hard device.
            CALL PGSCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
	  CALL PGSLW(PLT_LINE_WGT)
	ELSE
	  CALL PGQCR(0,RTEMP,GTEMP,BTEMP)
	  CALL PGSLW(IONE)
	END IF
	CALL PGSVP(MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
	CALL PGERAS
!
! Set the desired aspect ratio of the plot (useful so that
! aspect ratio is preserved on the screen and hardcopy device)
! 0 implies leave aspect ratio at default value
!
	CALL PGQVP(ITWO,DXST,DXEND,DYST,DYEND)
	DASR=(DYEND-DYST)/(DXEND-DXST)
	IF(ASR .LT. 0) TEMPASR=-1.0/ASR
	IF(ASR .GT. 0) TEMPASR=ASR
	IF(ASR .EQ. 0) TEMPASR=DASR
!
	TOTXMM=(DXEND-DXST)
	TOTYMM=(DYEND-DYST)
 350	IF(HARD) CALL GEN_IN(XCM,'Plot size')
!
! Commented out to test the box squareing routine
! Gregson Vaux December 1996
!
	IF(HARD .AND. (XCM .NE. 0)) THEN
	  SCALEFAC=XCM*10/TOTXMM
	  SCALEFACY=XCM*10*TEMPASR/TOTYMM
	ELSE IF (DASR .LE. 1.0) THEN
	  SCALEFAC=1
	  SCALEFACY=TEMPASR/DASR
          IF(DASR .LT. TEMPASR)THEN
	    T1=DASR/TEMPASR
	    SCALEFAC=SCALEFAC*T1
	    SCALEFACY=SCALEFACY*T1
	  END IF
	ELSE
	  SCALEFAC=TEMPASR/DASR
	  SCALEFACY=1
          IF(DASR .GT. TEMPASR)THEN
	    T1=TEMPASR/DASR
	    SCALEFAC=SCALEFAC*T1
	    SCALEFACY=SCALEFACY*T1
	  END IF
	ENDIF
!
	PRINTX1=MARGINX(1)
	PRINTX2=MARGINX(1)+(MARGINX(2)-MARGINX(1))*SCALEFAC
!	PRINTX1_SCL=MARGINX_SCL(1)
!	PRINTX2_SCL=MARGINX_SCL(2)
	PRINTX1_SCL=MARGINX(1)+(MARGINX_SCL(1)-MARGINX(1))*SCALEFAC
	PRINTX2_SCL=MARGINX(1)+(MARGINX_SCL(2)-MARGINX(1))*SCALEFAC
	PRINTY1=MARGINY(1)
	PRINTY2=MARGINY(1)+(MARGINY(2)-MARGINY(1))*SCALEFACY
        WRITE(6,*)PRINTX1,PRINTX2,PRINTX1_SCL,PRINTX2_SCL,PRINTY1,PRINTY2
	IF(HARD) THEN
	  IF(PRINTX2 .GT. 1.0 .OR. PRINTY2 .GT. 1.0)THEN
	    WRITE(T_OUT,*)' Error - plot too large - try again'
	    GOTO 350
	  ENDIF
	ENDIF
	CALL PGSVP(PRINTX1,PRINTX2,PRINTY1,PRINTY2)
!
! 
	IF(TYPE_OF_IMAGE .EQ. 'VEC')THEN
!
! Draw vectorplot.
!
! Draw Polarization vectors. Factor of 0.5 included since
! we draw the vector either side of the location.
!
	  CALL PGSLW(LINE_WGT(1))
	  CALL PGSLS(LINE_STYLE(1))
	  Q=PEN_COL(2)
	  CALL PGSCI(Q)     ! START WITH COLOR INDEX 2
!
          MAXSIZEX=0.5*SCALEPOL*(XPAR(2)-XPAR(1))/( MAPMAX*XDIM )
	  MAXSIZEY=MAXSIZEX
!
! NB - Theta is defined relative to the Y axis in a clockwise
! sense (i.e. towards +X axis).
!
	DO J=1,YDIM
	  DO I=1,XDIM
	    IF(HAN)THEN
	      KST=-XHAN/2
	      KEND=XHAN/2
	      IF(I+KST .LT. 1)KST=1-I
	      IF(I+KEND .GT. XDIM)KEND=XDIM-I
	      LST=-YHAN/2
	      LEND=YHAN/2
	      IF(I+LST .LT. 1)LST=1-I
	      IF(I+LEND .GT. YDIM)LEND=YDIM-I
	      T1=0
	      T2=0
	      DO L=LST,LEND
	        DO K=KST,KEND
	          T1=T1+MAP(I+K,J+L)*SIN(THETA(I+K,J+L))*SM_ARRAY(K,J)
	          T2=T2+MAP(I+K,J+L)*COS(THETA(I+K,J+L))*SM_ARRAY(K,J)
	        END DO
	      END DO
	      X1=X1-T1*0.5*MAXSIZEX
	      Y1=Y1-T2*0.5*MAXSIZEX
	      X2=X1+T1*0.5*MAXSIZEY
	      Y2=Y1+T2*0.5*MAXSIZEY
	    ELSE
	      X1=XVEC(I)+MAXSIZEX*MAP(I,J)*SIN(THETA(I,J))
	      Y1=YVEC(J)+MAXSIZEY*MAP(I,J)*COS(THETA(I,J))
	      X2=XVEC(I)-MAXSIZEX*MAP(I,J)*SIN(THETA(I,J))
	      Y2=YVEC(J)-MAXSIZEY*MAP(I,J)*COS(THETA(I,J))
	    END IF
	    CALL PGMOVE(X1,Y1)
	    CALL PGDRAW(X2,Y2)
	  END DO
	END DO
!
! 
!
	ELSE IF (TYPE_OF_IMAGE .EQ. 'COL')THEN
	  CALL PGPAGE
	  CALL PGQCIR(I,J)
	  CALL PGSCIR(I,J)
	  CONTRAST=1.0
	  BRIGHTNESS=0.5
	  CALL GEN_IN(IPALLET,'Pallet?')
	  CALL PALETT(IPALLET,CONTRAST,BRIGHTNESS)
	  CALL PGIMAG(MAP,XDIM,YDIM,IONE,XDIM,IONE,YDIM,ZPAR(1),ZPAR(2),TR_VEC)
	  IF(DO_SCALE)THEN
	    CALL PGSVP(MARGINX_SCL(1),MARGINX_SCL(2),MARGINY(1),MARGINY(2))
	    T1=0.0; T2=1.0
	    CALL PGSWIN(T1,T2,T1,T2)
	    DO I=1,NSCL
	      SCALE_VEC(1:4,I)=ZPAR(1)+(I-1)*(ZPAR(2)-ZPAR(1))/(NSCL-1)
	    END DO
	    I=4
	    WRITE(6,*)'Doing scale'
	    CALL PGIMAG(SCALE_VEC,I,NSCL,IONE,I,IONE,NSCL,ZPAR(1),ZPAR(2),TR_VEC_SCL)
!
! Reset to old values.
!
!	    CALL PGSVP(MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	    CALL PGSVP(PRINTX1,PRINTX2,PRINTY1,PRINTY2)
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
	  END IF
	ELSE IF (TYPE_OF_IMAGE .EQ. 'GREY')THEN
	  WRITE(6,*)'Ploting grey image'
	  T1=0.0
	  CALL PGGRAY(MAP,XDIM,YDIM,IONE,XDIM,IONE,YDIM,
	1                    ZPAR(1),ZPAR(2),TR_VEC)
	  IF(DO_SCALE)THEN
!	    CALL PGSVP(MARGINX_SCL(1),MARGINX_SCL(2),MARGINY(1),MARGINY(2))
	    CALL PGSVP(PRINTX1_SCL,PRINTX2_SCL,PRINTY1,PRINTY2)
	    T1=0.0; T2=1.0
	    CALL PGSWIN(T1,T2,T1,T2)
	    DO I=1,NSCL
	      SCALE_VEC(1:4,I)=ZPAR(1)+(I-1)*(ZPAR(2)-ZPAR(1))/(NSCL-1)
	    END DO
	    I=4
	    CALL PGGRAY(SCALE_VEC,I,NSCL,IONE,I,IONE,NSCL,ZPAR(1),ZPAR(2),TR_VEC_SCL)
!
! Reset to old values.
!
!	    CALL PGSVP(MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	    CALL PGSVP(PRINTX1,PRINTX2,PRINTY1,PRINTY2)
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
	  END IF
	ELSE IF (TYPE_OF_IMAGE .EQ. 'CONT')THEN
	  CALL PGCONT(MAP,XDIM,YDIM,IONE,XDIM,IONE,YDIM,CONTOUR_LEVELS,
	1                      NUM_CONTOUR_LEVELS,TR_VEC)
	END IF
!
	CALL PGSLW(PLT_LINE_WGT)
!
! Draw borders.
!
	I=1
	CALL PGSLS(I)
	Q=PEN_COL(1)
	CALL PGSCI(Q)   ! Draw box in color index 1
!
! *** PGBOX is commented out while MON_BORD is being tested
! September 1996 Gregson Vaux
!
!	CALL PGBOX('ABCNT',0.0,0,'ABCNT',0.0,0)
!
	Q=PEN_COL(2)
	CALL PGSCI(Q)     ! START WITH COLOR INDEX 2
	CALL MONBORD_V4(XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITLE,N_TITLE,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS,
	1             XLAB_FILE,YLAB_FILE)
	IF(.NOT. NORMAL_R_Y_AXIS)THEN
	  CALL DRAW_RIGHT_Y_AXIS(YPAR_R_AX,YINC_R_AX,YNUMST_R_AX,
	1        IYTICK_R_AX,IDY_R_AX,TICK_FAC,
	1        EXPCHAR,YLABEL_R_AX,LOG_AXIS)
	END IF
!
! Draw strings on graph
!
	IF(STR)THEN
	  CALL JUSTIFY_CONVERT_V2(XSTR,YSTR,LOC,LOC_PG,ORIENTATION,FLAGSTR,
     *    XSTRPOS,YSTRPOS,STRING,MAXSTR)
	  DO I=1,MAXSTR
	    IF(FLAGSTR(I))THEN
	      CALL PGSCI(STR_COL(I))
	      CALL PGSCH(EXPCHAR*STR_EXP(I))
	      CALL PGPTXT(XSTRPOS(I),YSTRPOS(I),ORIENTATION(I),LOC_PG(I),
	1                    STRING(I))
	    END IF
	  END DO
	END IF
!
! Draw vectors on graphs.
!
	IF(VEC)THEN
	  DO I=1,MAXVEC
	    Q=VECPEN(I)
	    CALL PGSCI(Q)
	    IF(FLAGLINE(I))THEN
	      CALL PGMOVE(LINEXST(I),LINEYST(I))
	      CALL PGDRAW(LINEXEND(I),LINEYEND(I))
	    END IF
	  END DO
	END IF
!
! Centre text output on terminal
!
	CALL PGMOVE((XPAR(2)-XPAR(1))/2.0,(YPAR(2)-YPAR(1))/2.0)
	IF(HARD) THEN
	  CALL PGEBUF		!send plot to printer file
	  CALL PGCLOS
	  CALL PGSLCT(ID)
	  HARD=.FALSE.
	END IF
	CALL PGSCR(0,RTEMP,GTEMP,BTEMP)
!
	GO TO 1000
	END
