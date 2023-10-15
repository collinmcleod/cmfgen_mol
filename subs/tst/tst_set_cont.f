	PROGRAM TST_SET_CONT
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
C
	INTEGER*4 NCF_CONT
	INTEGER*4 NFREQ
	INTEGER*4 N_EDGE
	INTEGER*4 N_LINES
	INTEGER*4, PARAMETER :: NFREQ_MAX=200000
	INTEGER*4, PARAMETER :: NLINE_MAX=200000
	INTEGER*4, PARAMETER :: LUOUT=40
!
	REAL(10) NU_EVAL(NFREQ_MAX)		!Frequencies at which continuum is evaluated.
	REAL(10) NU_CONT(NFREQ_MAX)		!Pure continuum frequencies
	REAL(10) FREQ(NFREQ_MAX)			!Continuum frequencies
	REAL(10) EDGE_FREQ(NFREQ_MAX) 		
	CHARACTER*6 EDGE_TYPE(NFREQ_MAX)	!End index for the line 
	INTEGER*4 INDX(NFREQ_MAX)
!
	REAL(10) NU_LINE(NLINE_MAX)
	REAL(10) NU_STRT_LINE(NLINE_MAX)
	REAL(10) VEC_MIN_VDOP(NLINE_MAX)
	INTEGER*4 LINE_ST_INDX(NLINE_MAX)
	INTEGER*4 LINE_END_INDX(NLINE_MAX)
	CHARACTER*6 TRANS_TYPE(NLINE_MAX)
	CHARACTER*1 CONT_TYPE(NFREQ_MAX)
	CHARACTER*80 FILENAME
	CHARACTER*80 STRING
!
	INTEGER*4 LINES_THIS_FREQ(NFREQ_MAX)
!
! Passed constants:
!
        REAL(10) VINF             !Terminal velocity of wind.
        REAL(10) FRAC_DOP         !Indicates dNU across line in Doppler widths.
C
        REAL(10) MAX_DOP
        REAL(10) dV_CMF_PROF
        REAL(10) dV_CMF_WING
        REAL(10) ES_WING_ExT
        REAL(10) R_CMF_WING_EXT
!
	REAL(10) MAX_FREQ
	REAL(10) MIN_FREQ
	REAL(10) SMALL_RAT
	REAL(10) BIG_AMP
	REAL(10) DNU_MAX
	REAL(10) dV_LEV
	REAL(10) AMP_DIS
	REAL(10) MIN_FREQ_LEV_DIS
!
	REAL(10) dV_DOP
	REAL(10) dV_CONT
	REAL(10) MIN_dV_CONT
!
	REAL(10) C_KMS,T1
	REAL(10) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER *4 I,ML,K,J,ISEED
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
!
	FILENAME='EDGE_FREQ'
	CALL GEN_IN(FILENAME,'Continuum edge file')
	OPEN(UNIT=10,FILE=FILENAME,STATUS='OLD',READONLY)
	ML=0
	DO WHILE(1 .EQ. 1)
	 ML=ML+1
	  READ(10,*,END=100)EDGE_FREQ(ML),EDGE_TYPE(ML)
	  N_EDGE=ML
	END DO
100	CONTINUE
	CLOSE(UNIT=10)
	EDGE_TYPE(1:N_EDGE)='D'
	ISEED=-145623
	DO I=1,N_EDGE
	  J=1+3.0*RAN(ISEED)
	  IF(J .EQ. 1 .OR. J .EQ. 2)EDGE_TYPE(I)='S'
	END DO
	DO I=1,N_EDGE
	  WRITE(51,'(X,ES16.6,2X,A)')EDGE_FREQ(I),EDGE_TYPE(I)
	END DO
!
	MAX_FREQ=80.0D0
	CALL GEN_IN(MAX_FREQ,'Maximum frequency')
	MIN_FREQ=0.0005
	CALL GEN_IN(MIN_FREQ,'Minimum frequency')
!
	SMALL_RAT=1.10
	BIG_AMP=1.05
	DNU_MAX=0.1
!
	dV_LEV=200.0D0
	AMP_DIS=1.4D0
	MIN_FREQ_LEV_DIS=0.1D0
	dV_DOP=10.0D0
	dV_CONT=1000.0D0
!
	CALL GEN_IN(dV_CONT,'Maximum continuum spacing in km/s')
	CALL GEN_IN(dV_LEV,'Spacing on red side of edge for level dissolution')
	CALL GEN_IN(AMP_DIS,'Amplification factor on red side of edge for level dissolution')
!
	CALL SET_CONT_FREQ_V2(NU_CONT,CONT_TYPE, EDGE_FREQ,EDGE_TYPE,INDX,
	1                        SMALL_RAT,BIG_AMP,DNU_MAX,MAX_FREQ,MIN_FREQ,
	1                        dV_LEV,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        dV_CONT,dV_DOP,
	1                        N_EDGE,NCF_CONT,NFREQ_MAX,LUOUT)
C
	WRITE(6,*)'Returned: NCF_CONT=',NCF_CONT
	DO ML=1,NCF_CONT-1
	  WRITE(17,'(I6,2X,3F15.6)')
	1         ML,NU_CONT(ML), NU_CONT(ML)-NU_CONT(ML+1), 
	1         C_KMS*(NU_CONT(ML)-NU_CONT(ML+1))/NU_CONT(ML)
	END DO
!
! We will now test the line insertion routine.
!
	FILENAME(1:4)='LINE'
	CALL GEN_IN(FILENAME,'Line frequency file')
	OPEN(UNIT=10,FILE=FILENAME,STATUS='OLD',READONLY)
	  READ(10,'(A)')STRING
	  IF(INDEX(STRING,'Lam') .NE. 0)THEN
            ML=1
	    DO WHILE(ML .LE. NLINE_MAX)
	      READ(10,*,END=200)I,I,I,NU_LINE(ML)
              ML=ML+1
	    END DO
	  ELSE
	    REWIND(10)
            ML=1
	    DO WHILE(ML .LE. NLINE_MAX)
	      READ(10,*,END=200)NU_LINE(ML)
              ML=ML+1
	    END DO
	  END IF
200     N_LINES=ML-1
	WRITE(6,*)'Nubmer of lines read in is ',N_LINES
    	VEC_MIN_VDOP(1:N_LINES)=20.0D0            !km/s
        TRANS_TYPE(1:N_LINES)='BLANK'
!
        MIN_dV_CONT=10.0              !km/s
        FRAC_DOP=1.0D0
        MAX_DOP=6.0D0
        VINF=840.0D0            !km/s
!
	DO ML=1,N_LINES
	  NU_STRT_LINE(ML)=NU_LINE(ML)*(1.0D0+6.0D0*20.0D0/C_KMS)
	END DO
!	WRITE(6,*)NU_LINE(1:N_LINES),NU_STRT_LINE(1:N_LINES)
!
        dV_CMF_PROF=100.0D0               !km/s
	CALL GEN_IN(dV_CMF_PROF,'Spacing in km/s across CMF profile')
        dV_CMF_WING=300.0D0
	CALL GEN_IN(dV_CMF_WING,'Spacing in km/s across e.s. wings')
        ES_WING_ExT=2500.0              !km/s
        R_CMF_WING_EXT=3.0D0
!
        CALL INS_LINE_V6(FREQ,LINES_THIS_FREQ,NFREQ,NFREQ_MAX,
	1              NU_LINE,NU_STRT_LINE,VEC_MIN_VDOP,TRANS_TYPE,
	1              LINE_ST_INDX,LINE_END_INDX,N_LINES,
	1              NU_CONT,CONT_TYPE,NCF_CONT,MIN_dV_CONT,
	1              FRAC_DOP,VINF,dV_CMF_PROF,
	1              dV_CMF_WING,ES_WING_EXT,R_CMF_WING_EXT )
!
	WRITE(6,*)'Exited INS_LINE'
!
!	K=1
!	FREQ(K)=NU_CONT(1)
!	DO ML=2,NCF_CONT
!	  T1=0.25D0*(NU_CONT(ML)-NU_CONT(ML-1))
!          DO J=1,4
!            K=K+1
!	    FREQ(K)=FREQ(K-1)+T1
!	  END DO
!	END DO
!	NFREQ=K
!
	CALL DET_MAIN_CONT_FREQ(FREQ,NFREQ,NU_CONT,NCF_CONT,
	1            NU_EVAL,.TRUE.,.FALSE.)
C
	STOP
	END
