!
! Subroutine to measure line EWs using direct numerical integration of the data.
! Routine is designed to be called in GRAMON_PGPLOT.
!
	SUBROUTINE DO_FILE_EW_V1(FILE_WITH_LINE_LIMS)
	USE MOD_CURVE_DATA
	USE GEN_IN_INTERFACE
	USE MOD_EW_VARIABLES
	IMPLICIT NONE
!
! Altered: 30-Jun-2022 - Can no append transition name to EW file.
! Created: 27-FEb-2022
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
	CHARACTER(LEN=80) FILE_WITH_LINE_LIMS
!
!
	INTEGER, SAVE :: NPIX=3		!Integraton band pass around line limits
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
	REAL*4, SAVE :: CONT_ACC=0.2		!percentage
	REAL*4, SAVE :: LOW_CONT=0.998
	REAL*4, SAVE :: HIGH_CONT=1.002
!
	INTEGER IPEN
	CHARACTER(LEN=1) CURSVAL
!
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, SAVE :: FIRST_WRITE=.TRUE.
	INTEGER I,J,K
	INTEGER IP
	INTEGER, SAVE :: LUIN=0
	INTEGER, SAVE :: LUOUT=0
!
	INTEGER GET_INDX_SP
	EXTERNAL GET_INDX_SP
!
	LOGICAL END_FILE
	LOGICAL FILE_PRES
	CHARACTER(LEN=80) LOC_FILE_WITH_LINE_LIMS
	CHARACTER(LEN=80) OUT_FILE 
	CHARACTER(LEN=80) STRING
!
! Used if average data on X-limits to defined the continuum level.
!
	NPIX=1
	CALL GEN_IN(NPIX,'Number of pixels at X location to average continuum (must be odd)')
	IF(MOD(NPIX,2) .EQ. 0)THEN
	  NPIX=NPIX+1
	  WRITE(6,*)'NPIX increase by 1 to make odd; NPIX=',NPIX
	END IF
	CALL GEN_IN(USE_MILLI_ANG,'Output EWs in milli-Angstroms?')
	CALL GEN_IN(CONT_ACC,'Measure accuracy for continuum -- percentage?')
	LOW_CONT=1.0-CONT_ACC/100.0
	HIGH_CONT=1.0+CONT_ACC/100.0
!
	IF(LUIN .EQ. 0)CALL GET_LU(LUIN,'LUIN in DO_CURSOR_EW_V2')
	LOC_FILE_WITH_LINE_LIMS='LINE_LIMS'
	CALL GEN_IN(LOC_FILE_WITH_LINE_LIMS,'File with XSTART, XEND defining lines')
	OPEN(UNIT=LUIN,FILE=LOC_FILE_WITH_LINE_LIMS,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error the follwing file cannot be opened',TRIM(LOC_FILE_WITH_LINE_LIMS)
	  RETURN
	END IF
!
! Open output file. Data is appended if it alread exists.
!
	IF(LUOUT .EQ. 0)THEN
	  OUT_FILE='EW_DATA'
	  CALL GET_LU(LUOUT,'LUOUT in DO_CURSOR_EW_V2')
	  CALL GEN_IN(OUT_FILE,'File to OUTPUT EWs etc')
	  INQUIRE(FILE=OUT_FILE,EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)'File already exists -- appending new data'
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	  ELSE
	    OPEN(UNIT=LUOUT,FILE=OUT_FILE,STATUS='NEW',ACTION='WRITE')
	    CALL WRITE_EW_HEADER(LUOUT)
	    CALL WRITE_EW_HEADER(6)
	  END IF
	END IF
!
! File input
!
	DO WHILE(1 .EQ. 1)
	  STRING(1:1)='!'
	  DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    END_FILE=.TRUE.
	    READ(LUIN,'(A)',END=2000)STRING
	    END_FILE=.FALSE.
	  END DO
	  READ(STRING,*)XST,XEND
!
	  DO IP=1,NPLTS
!
! Get continuum fluxes at line limits.
!
	    IST=GET_INDX_SP(XST,CD(IP)%XVEC,NPTS(IP))
	    YST=0.0
	    DO J=IST-NPIX/2,IST+NPIX/2
	      YST=YST+CD(IP)%DATA(J)
	    END DO
	    YST=YST/NPIX
!
	    IEND=GET_INDX_SP(XEND,CD(IP)%XVEC,NPTS(IP))
	    YEND=0.0
	    DO J=IEND-NPIX/2,IEND+NPIX/2
	      YEND=YEND+CD(IP)%DATA(J)
	    END DO
	    YEND=YEND/NPIX
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
	      dX=(CD(IP)%XVEC(MIN(I+1,IEND))-CD(IP)%XVEC(MAX(IST,I-1)))/2
	      EW=EW+dX*(YVAL-CD(IP)%DATA(I))
	      EWL=EW+dX*(LOW_CONT*YVAL-CD(IP)%DATA(I))
	      EWH=EW+dX*(HIGH_CONT*YVAL-CD(IP)%DATA(I))
	      XMEAN=XMEAN+CD(IP)%XVEC(I)*(YVAL-CD(IP)%DATA(I))*dX
	    END DO
!
	    IF(ABS(EW)/YST .LT. 1.0D-10)THEN
	      WRITE(6,*)'Possibe error EW is close to zero', EW/YST
	      WRITE(6,*)'Skipping this line. XST,END=',XST,XEND
	      GOTO 1000
	    END IF
	    XMEAN=XMEAN/EW
	    YCONT=YST+(XMEAN-CD(IP)%XVEC(IST))*SLOPE
	    YINT=EW
	    YMEAN=EW/(CD(IP)%XVEC(IEND)-CD(IP)%XVEC(IST))
	    EW=EW/YCONT
	    EWL=EWL/(YCONT*LOW_CONT)
	    EWH=EWH/(YCONT*HIGH_CONT)
!
! These parameters are used to provide information on whether a line is blended.
!
	    SIGMA=0.0D0; SKEWNESS=0.0D0; KURTOSIS=0.0D0
	    DO I=IST,IEND
	      XVAL=CD(IP)%XVEC(I)
	      YVAL=YST+(CD(IP)%XVEC(I)-CD(IP)%XVEC(IST))*SLOPE - CD(IP)%DATA(I)
	      dX=(CD(IP)%XVEC(MIN(I+1,IEND))-CD(IP)%XVEC(MAX(IST,I-1)))/2
	      SIGMA=SIGMA+dX*YVAL*(XVAL-XMEAN)**2
	      SKEWNESS=SKEWNESS+dX*YVAL*(XVAL-XMEAN)**3
	      KURTOSIS=KURTOSIS+dX*YVAL*(XVAL-XMEAN)**4
	    END DO
	    IF(SIGMA*YINT .LE. 0)THEN
	      SIGMA=-1.0; SKEWNESS=-1.0; KURTOSIS=-1.0
	    ELSE
	      SIGMA=SQRT(SIGMA/YINT)
	      SKEWNESS=SKEWNESS/YINT/SIGMA**3
	      KURTOSIS=KURTOSIS/YINT/SIGMA**4
	    END IF
	    CALL GET_LINE_ID_PG(TRANS_NAME,EW,XMEAN,T1)
!
	    CALL WR_EW_VARIABLES(IP,LUOUT)
	    I=6; CALL WR_EW_VARIABLES(IP,I)
	  END DO
!
1000	  CONTINUE
	END DO
2000	CONTINUE
!
	RETURN
	END
