!
! Routine to perform an NG accleration for a Comoving-Frame Model. Progam
! uses the last 4 iterations which are stored in the last "4 records" 
! (effectively) of SCRTEMP.
!
! The NG acceleration is perfomed separately on each depth, or over a band
! of depths.
!
! Output is to the last record of SCRTEMP.
!
! Input files required:
!		     MODEL(.DAT)
!                    POINT1(.DAT)
!                    SCRTEMP(.DAT)
! Output files:
!                    POINT1(.DAT)
!                    POINT2(.DAT)
!                    SCRTEMP(.DAT)
!
	PROGRAM DO_NG_V2
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 19-May-2004 : Bug fix: NG_DONE set to true, even when no NG done.
! Altered 28-Mar-2004 : Changed to handle new format SCRTEMP files.
!                       Now choice of 3 options: NG, AV and SOR.
! Altered 01-Jul-2003 : IT_STEP inserted as option.  
! Altered 15-Apr-2003 : IFLAG in RD_4_ITS initialized.
!                       NUM_BAD_NG in NG_MIT_OPST initialized.
! Altered 18-Feb-2001 : Extensive rewrite.
!                       NBAND option installed.
!
! Altered 04-Jan-1998 : Cleaned. ND, NT now read from MODEL file.
!                       Based on a very old version in [JDH.BASOL] and 
!                        REWRITE_SCR.
!
	REAL*8, ALLOCATABLE :: RDPOPS(:,:,:)		!NT+3,ND
	REAL*8, ALLOCATABLE :: BIG_POPS(:,:)		!NT+3,ND
	REAL*8, ALLOCATABLE :: POPS(:,:)		!NT,ND
	REAL*8, ALLOCATABLE :: R(:)			!ND
	REAL*8, ALLOCATABLE :: V(:)			!ND
	REAL*8, ALLOCATABLE :: SIGMA(:)			!ND
!
! Local variables which are adjusted to match the particular model under
! consideration.
!
	REAL*8 T1,T2,T3
	REAL*8 SCALE_FAC
	REAL*8 BIG_FAC
!
	INTEGER*4 ND,NT
	INTEGER*4 NBAND
!
	INTEGER*4 IOS
	INTEGER*4 IREC
	INTEGER*4 NITSF
	INTEGER*4 LST_NG
	INTEGER*4 IFLAG
	INTEGER*4 N_ITS_TO_RD
	INTEGER*4 IT_STEP
	INTEGER*4 I,J
!
	INTEGER*4, PARAMETER :: RITE_N_TIMES=1
	INTEGER*4, PARAMETER :: T_OUT=6
	INTEGER*4, PARAMETER :: LUSCR=10
!
	LOGICAL NEWMOD
	LOGICAL NG_DONE
	LOGICAL DO_REGARDLESS
	LOGICAL WRITE_RVSIG
	LOGICAL SCALE_INDIVIDUALLY
	CHARACTER*3 OPTION
	CHARACTER*132 STRING
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       MODEL(.DAT)'
	WRITE(T_OUT,*)'                                       POINT1(.DAT)'
	WRITE(T_OUT,*)'                                       SCRTEMP(.DAT)'
	WRITE(T_OUT,*)' '
!
! Get basic Model data (i.e., NT and ND).
!
	OPEN(UNIT=12,FILE='MODEL',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of depth') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)ND
	  DO WHILE(INDEX(STRING,'!Total number of variables') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)NT
!
! Iw we couldn't successfully red the MODEL file, get NT and ND from terminal.
!
100	CONTINUE
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read MODEL file'
	  CALL GEN_IN(NT,'Total number of levels')
	  CALL GEN_IN(ND,'Number of depth points')
	END IF 
	CLOSE(UNIT=12)
!
	ALLOCATE (RDPOPS(NT+3,ND,4))
	ALLOCATE (BIG_POPS(NT+3,ND))
	ALLOCATE (POPS(NT,ND))
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
!
	OPTION='NG'
	CALL GEN_IN(OPTION,'Acceleration method: NG, AV(ERAGE), SOR')
	CALL SET_CASE_UP(OPTION,0,0)
	IF(OPTION(1:2) .EQ.'NG')THEN
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')' This option was inserted mainly for testing purposes.'
	  WRITE(T_OUT,'(A)')' The default value of 1 should generally be used.'
	  IT_STEP=1
	  CALL GEN_IN(IT_STEP,'Iteration step size for NG acceleration')
	  NBAND=1
	  DO_REGARDLESS=.FALSE.
	  CALL GEN_IN(NBAND,'Band width for NG acceleration')
	  CALL GEN_IN(DO_REGARDLESS,'Do acceleration independent of corection size?')
	  SCALE_INDIVIDUALLY=.FALSE.
	  IF(DO_REGARDLESS)THEN
	    CALL GEN_IN(SCALE_INDIVIDUALLY,'Scale each population individually at each depth to limit change')
	  END IF
	  N_ITS_TO_RD=4
	ELSE IF(OPTION(1:2) .EQ.'AV')THEN
	  IT_STEP=1
	  N_ITS_TO_RD=2
	ELSE IF(OPTION(1:3) .EQ. 'SOR')THEN
	  IT_STEP=1
	  N_ITS_TO_RD=2
	  CALL GEN_IN(SCALE_FAC,'Factor to scale last correction by')
          CALL GEN_IN(BIG_FAC,'Maximum correction to any depth')
	ELSE
	   WRITE(6,*)'Invalid acceleration option'
	   STOP
	END IF
!
! Read POPULATIONS that were output on last iteration. This is primarily
! done to get NITSF etc.
! 
	NEWMOD=.FALSE.
	IREC=0			!Get last iteration
	CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	IF(NEWMOD)THEN
	  WRITE(T_OUT,*)'Unable to read last iteration in scratch file'
	  STOP
	END IF
	BIG_POPS(1:NT,:)=POPS(1:NT,:)
	BIG_POPS(NT+1,:)=R(:)
	BIG_POPS(NT+2,:)=V(:)
	BIG_POPS(NT+3,:)=SIGMA(:)+1.0D0
	WRITE(6,*)'Read in last iteration'
!
! Read in the last N_ITS_TO_RD estimates of the populations, as output to SCRTEMP.
!
	NEWMOD=.FALSE.
	DO I=1,N_ITS_TO_RD
	  J=IREC-(I-1)*IT_STEP
	  WRITE(6,*)'Reading record', J
	  CALL SCR_READ_V2(R,V,SIGMA,POPS,J,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  RDPOPS(1:NT,:,I)=POPS(1:NT,1:ND)
	  RDPOPS(NT+1,:,I)=R(:)
	  RDPOPS(NT+2,:,I)=V(:)
	  RDPOPS(NT+3,:,I)=SIGMA(:)+1.0D0
	  IF(NEWMOD)THEN
	    WRITE(T_OUT,*)'Unable to read scratch file'
	    STOP
	  END IF
	END DO

	WRITE(6,*)'Starting the accleration'
	IF(OPTION(1:2) .EQ. 'NG')THEN
!
! Perform the NG acceleration.
!
	  J=NT+3
	  CALL NG_MIT_OPTS(BIG_POPS,RDPOPS,ND,J,NBAND,DO_REGARDLESS,
	1                     SCALE_INDIVIDUALLY,NG_DONE,T_OUT)
!
	  WRITE(6,*)'Finished NG accleration'
	ELSE IF(OPTION(1:2) .EQ. 'AV')THEN
	  DO J=1,ND
	    DO I=1,NT+3
	      BIG_POPS(I,J)=0.5D0*(RDPOPS(I,J,1)+RDPOPS(I,J,2))
	    END DO
	  END DO
	  NG_DONE=.TRUE.
	ELSE IF(OPTION(1:3) .EQ. 'SOR')THEN
	  DO J=1,ND
	    T1=-1000.0D0
	    T2=1000.0D0
	    DO I=1,NT+3
	      T1=MAX(T1,(RDPOPS(I,J,1)-RDPOPS(I,J,2))/RDPOPS(I,J,1))
	      T2=MIN(T2,(RDPOPS(I,J,1)-RDPOPS(I,J,2))/RDPOPS(I,J,1))
	    END DO
	    T2=ABS(T2)+1.0D0
	    T3=SCALE_FAC
	    IF(T3*T1 .GT. BIG_FAC)T3=BIG_FAC/T1
	    IF(T3*T2 .GT. (1.0D0-1.0D0/BIG_FAC))T3=BIG_FAC/T2
	    DO I=1,NT+3
	      BIG_POPS(I,J)=RDPOPS(I,J,1)+T3*(RDPOPS(I,J,1)-RDPOPS(I,J,2))
	    END DO
	  END DO
	  NG_DONE=.TRUE.
	END IF
!
	POPS(1:NT,:)=BIG_POPS(1:NT,:)
	R(2:ND-1)=BIG_POPS(NT+1,2:ND-1)
	V(2:ND-1)=BIG_POPS(NT+2,2:ND-1)
	SIGMA(2:ND-1)=BIG_POPS(NT+3,2:ND-1)-1.0D0
!
	IF(NG_DONE)THEN
	  NITSF=NITSF+1
	  LST_NG=NITSF
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,
	1                 LST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  IF(OPTION(1:2) .EQ.'NG')THEN
	    WRITE(T_OUT,*)'Results of successful NG acceleration ',
	1                 'output to SCRTEMP.'
	  ELSE IF(OPTION(1:3) .EQ.'SOR')THEN
	    WRITE(T_OUT,*)'Results of successful SOR output to SCRTEMP.'
	  ELSE
	    WRITE(T_OUT,*)'Last 2 iterations averaged and output to SCRTEMP'
	  END IF
	END IF
!
	STOP
	END
!
!
!
! Altered 3-April-1989. -Call changed. MAXDEC, MAXINC and IFLAG installed,
!                        ACCELERATE flag removed.
!
! If IFLAG is returned with zero, the NG acceleration was successful.
! Otherwise an error condition occurred.
!
	SUBROUTINE GENACCEL_V2(NEWPOP,RDPOPS,LST,LEND,NT,ND,NS)
	IMPLICIT NONE
!
	INTEGER*4 LST,LEND
	INTEGER*4 NT
	INTEGER*4 ND
	INTEGER*4 NS
	REAL*8 NEWPOP(NS)
	REAL*8 RDPOPS(NT,ND,4)
!
	REAL*8 TEMP(4,NS)
	INTEGER*4 I,J,K,L
	LOGICAL WEIGHT
!
! Use the NG acceleration method to improve the population estimates. 
! We can perform the NG acceleration at each depth individually, or
! over a range of depths. Weight is used to indicate that we
! are to minimize the percentage errors - not the absolute magnitudec
! of the errors. The absolute maximum percentage change is also
! determined.
!
! Rewrite the relevant poulations in a form suitable for NGACCEL.
!
	DO J=1,4
	  DO L=LST,LEND
	    DO I=1,NT
	      K=I+NT*(L-LST)
	      TEMP(J,K)=RDPOPS(I,L,J)
	    END DO
	  END DO
	END DO
	WEIGHT=.TRUE.
!
	WRITE(6,*)'Calling NGACCEL'
	CALL NGACCEL(NEWPOP,TEMP,NS,WEIGHT)
!
	RETURN
	END
!
!
	SUBROUTINE NG_MIT_OPTS(POPS,RDPOPS,ND,NT,NBAND,DO_REGARDLESS,
	1                           SCALE_INDIVIDUALLY,NG_DONE,LUER)
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 30-May-2003: Bug fix: When SCALE_INDIVIDUALLY was true, not all
!                         populations were being accelerated.
!
	INTEGER*4 ND
	INTEGER*4 NT
	INTEGER*4 NBAND
	INTEGER*4 LUER
	REAL*8 POPS(NT,ND)
	REAL*8 RDPOPS(NT,ND,4)
	LOGICAL DO_REGARDLESS
	LOGICAL SCALE_INDIVIDUALLY
	LOGICAL NG_DONE
!
	REAL*8 NEWPOP(NT,ND)
	REAL*8 VEC_INC(ND)
	REAL*8 VEC_DEC(ND)
	REAL*8 INT_ARRAY(ND)
!
	REAL*8 RINDX(NT)
	REAL*8 RAT(NT)
!
	REAL*8 T1
	REAL*8 LOCINC,LOCDEC
	REAL*8 MAXINC,MAXDEC
!
	INTEGER*4 NS
	INTEGER*4 NUM_BAD_NG
	INTEGER*4 DEC_LOC,INC_LOC
	INTEGER*4 I,K,L
	INTEGER*4 LST,LEND
	INTEGER*4, PARAMETER :: IONE=1
	LOGICAL DO_PLTS
!
	NUM_BAD_NG=0
	VEC_INC(1:ND)=0.0D0
	VEC_DEC(1:ND)=0.0D0
	IF(NBAND .GE. ND)THEN
	  NS=NT*ND
	  CALL GENACCEL_V2(NEWPOP,RDPOPS,IONE,ND,NT,ND,NS)
	ELSE
	  DO K=ND,1,-NBAND
	    LST=MAX(K-NBAND+1,1)
	    LEND=LST+NBAND-1
	    NS=(LEND-LST+1)*NT
	    CALL GENACCEL_V2(NEWPOP(1,LST),RDPOPS,LST,LEND,NT,ND,NS)
	  END DO
	END IF
!
	WRITE(6,*)' '
	WRITE(6,'(X,70A)')('*',I=1,70)
	WRITE(6,*)' '
        WRITE(6,*)'Plotting Y(NG accel)/Y(old) for each depth'
        WRITE(6,*)'Ten depths plotted each time.'
        WRITE(6,*)'These plots are for diagnostic purposes only'
        WRITE(6,*)' '
	WRITE(6,'(X,70A)')('*',I=1,70)
        WRITE(6,*)' '
	DO_PLTS=.FALSE.
        CALL GEN_IN(DO_PLTS,'Plot the above diagnostic information?')
	IF(DO_PLTS)THEN
	  DO L=1,ND
	    DO K=1,NT
	      RINDX(K)=K
	      RAT(K)=NEWPOP(K,L)/POPS(K,L)
	    END DO
	    CALL DP_CURVE(NT,RINDX,RAT)
	    IF(MOD(L,10) .EQ. 0 .OR. L .EQ. ND)THEN
	      WRITE(6,*)' ** Plotting depths',L-9+MOD(L,10),' to ',L
	      CALL GRAMON_PGPLOT('Index','Ratio',' ',' ')
	    END IF
	  END DO
	END IF
	WRITE(6,*)' '
!
! Now check whether the NG acceleation has been reasonable.
!
	MAXINC=-1000.0
	MAXDEC=1000.0
	DO L=1,ND 
	  LOCINC=-1000.0
	  LOCDEC=1000.0
	  DO K=1,NT
	    T1=NEWPOP(K,L)/POPS(K,L)
	    LOCINC=MAX(LOCINC,T1)
	    LOCDEC=MIN(LOCDEC,T1)
	  END DO
	  VEC_INC(L)=100.0D0*(LOCINC-1.0D0)
	  VEC_DEC(L)=100.0D0*(LOCDEC-1.0D0)
!
	  IF(DO_REGARDLESS)THEN
	    IF(SCALE_INDIVIDUALLY)THEN
	      DO K=1,NT
	        IF(NEWPOP(K,L) .LT. 0.1D0*POPS(K,L))THEN
	          POPS(K,L)=0.1D0*POPS(K,L)
	        ELSE IF(NEWPOP(K,L) .GT. 10.0D0*POPS(K,L))THEN
	          POPS(K,L)=10.0D0*POPS(K,L)
	        ELSE
	          POPS(K,L)=NEWPOP(K,L)
	        END IF
	      END DO
	    ELSE
	      T1=1.0D0
	      DO K=1,NT
                IF(NEWPOP(K,L) .LT. 0.1D0*POPS(K,L))THEN
                  T1=MIN( T1, 0.9D0*POPS(K,L)/(POPS(K,L)-NEWPOP(K,L)) )
                ELSE IF(NEWPOP(K,L) .GT. 10.0D0*POPS(K,L))THEN
                  T1=MIN( T1, 9.0D0*POPS(K,L)/(NEWPOP(K,L)-POPS(K,L)) )
                END IF
	      END DO
	      DO K=1,NT
                POPS(K,L)=POPS(K,L)+T1*(NEWPOP(K,L)-POPS(K,L))
              END DO
	      IF(T1 .NE. 1.0D0)WRITE(LUER,8000)L,LOCINC,LOCDEC
8000	      FORMAT(X,'Scaling performed at depth ',I3,':',
	1          X,'Biggest increase/decrease was ',2ES12.2)
	    END IF
!
	    IF(LOCINC .GT. MAXINC)THEN
	      MAXINC=LOCINC
	      INC_LOC=L
	    END IF
	    IF(LOCDEC .LT. MAXDEC)THEN
	      MAXDEC=LOCDEC
	      DEC_LOC=L
	    END IF
!
! Before storing the NG acceleration at this depth, we check to
! see whether the predicted corrections are "reasonable".
!
	  ELSE IF(LOCINC .GT. 10.1D0 .OR. LOCDEC .LT. 0.09D0)THEN
	    NUM_BAD_NG=NUM_BAD_NG+1
	    WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	    WRITE(LUER,9000)L,LOCINC,LOCDEC
9000	    FORMAT(X,'No NG acceleration performed at depth ',I3,/,
	1          X,'Biggest increase was ',1PE10.2,/,
	1          X,'Biggest decrease was ',E10.2)
	
	  ELSE
	    DO K=1,NT
	      POPS(K,L)=NEWPOP(K,L)
	    END DO
	    IF(LOCINC .GT. MAXINC)THEN
	      MAXINC=LOCINC
	      INC_LOC=L
	    END IF
	    IF(LOCDEC .LT. MAXDEC)THEN
	      MAXDEC=LOCDEC
	      DEC_LOC=L
	    END IF
	  END IF
	END DO
!
	WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	IF(NUM_BAD_NG .GT. 3)THEN
	  WRITE(LUER,*)'Too many bad NG accelerations - '//
	1              'NG acceleration cancelled'
	  NG_DONE=.FALSE. 
!
! Restore old population estimates for all depths.
!
	  DO L=1,ND
	    DO K=1,NT
	      POPS(K,L)=RDPOPS(K,L,1)
	    END DO
	  END DO
	  RETURN
	END IF
!
	WRITE(6,*)' '
	WRITE(6,'(X,70A)')('*',I=1,70)
	WRITE(6,*)' '
	WRITE(6,*)'Plotting largest % decrease and increase at each depth.'
	WRITE(6,*)' '
	WRITE(6,'(X,70A)')('*',I=1,70)
	WRITE(6,*)' '
	DO I=1,ND
	  INT_ARRAY(I)=I
	END DO
	CALL DP_CURVE(ND,INT_ARRAY,VEC_DEC)
	CALL DP_CURVE(ND,INT_ARRAY,VEC_INC)
	CALL GRAMON_PGPLOT('Depth Index','100(N/O-1)',' ',' ')
!
! By definition MAXINC and MININC are both positive. These are only
! defined for successful NG accelerations (does not include any depths
! at which NG acceleration did not work).
!
	MAXINC=100.0D0*(MAXINC-1.0D0)
	MAXDEC=100.0D0*(1.0D0/MAXDEC-1.0D0)
	WRITE(LUER,9800)INC_LOC,MAXINC
9800	FORMAT(X,'Max NG % increase at depth ',I3,' is',1PE10.2)
	WRITE(LUER,9900)DEC_LOC,MAXDEC
9900	FORMAT(X,'Max NG % decrease at depth ',I3,' is',1PE10.2)
C
	NG_DONE=.TRUE.  		!Successful acceleration
C
	RETURN
	END
