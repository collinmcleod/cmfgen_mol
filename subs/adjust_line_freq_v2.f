!
! Subroutine to replace program wavelength with an accurate wavelength.
! Matching is done by level name.
!
	SUBROUTINE ADJUST_LINE_FREQ_V2(LINE_FREQ,LINE_STRT_FREQ,VEC_TRANS_NAME,N_LINE_FREQ,LUIN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Atered  29-June-2011  : Changed to V2. LINE_STRT_FREQ inserted in file.
! Created 30-March-2008
!
	INTEGER N_LINE_FREQ
	INTEGER LUIN
	REAL(KIND=LDP) LINE_FREQ(N_LINE_FREQ)
	REAL(KIND=LDP) LINE_STRT_FREQ(N_LINE_FREQ)
	CHARACTER(LEN=*) VEC_TRANS_NAME(N_LINE_FREQ)
!
	REAL(KIND=LDP) WAVE,FREQ,ANG_TO_HZ,T1
	REAL(KIND=LDP) LAM_VAC,SPEED_OF_LIGHT
	CHARACTER(LEN=10) SPECIES,OLD_SPECIES
	CHARACTER(LEN=40) NAME
	CHARACTER(LEN=80) STRING
	INTEGER LUER,ERROR_LU,GET_INDX_DP
	EXTERNAL ERROR_LU,GET_INDX_DP,LAM_VAC,SPEED_OF_LIGHT
!
	INTEGER IOS
	INTEGER J,JST,K,L
!
	OPEN(UNIT=LUIN,FILE='REVISED_LAMBDAS',STATUS='OLD',IOSTAT=IOS,ACTION='READ')
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,'(/,1X,A)')'WARNING: file with accurate wavelengths not found'
	  RETURN
	END IF
!
	OLD_SPECIES='JJJ'
	ANG_TO_HZ=1.0D-07*SPEED_OF_LIGHT()
	DO WHILE(1 .EQ. 1)
	  STRING=' '
	  DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LUIN,'(A)',END=100)STRING
	  END DO
	  STRING=ADJUSTL(STRING)
	  K=INDEX(STRING,'  ')-1
	  NAME=STRING(1:K)
	  READ(STRING(K+1:),*)WAVE
	  IF(INDEX(STRING,'AIR') .NE. 0)WAVE=LAM_VAC(WAVE)
	  FREQ=ANG_TO_HZ/WAVE
!
	  L=INDEX(STRING,'(')-1
	  SPECIES=STRING(1:L)
!
	  IF(SPECIES .EQ. OLD_SPECIES)THEN
	  ELSE
	    DO J=1,N_LINE_FREQ
	      IF(SPECIES .EQ. VEC_TRANS_NAME(J)(1:L))THEN
	        JST=J
	        OLD_SPECIES=SPECIES
	        EXIT
	      END IF
	    END DO
	  END IF
!
! Now find the line. These are not ordered.
!
	  IF(SPECIES .EQ. OLD_SPECIES)THEN
	    DO J=JST,N_LINE_FREQ
	      IF(NAME(1:K) .EQ. VEC_TRANS_NAME(J))THEN
	        T1=LINE_STRT_FREQ(J)/LINE_FREQ(J)
	        LINE_FREQ(J)=FREQ
	        LINE_STRT_FREQ(J)=T1*FREQ
	        EXIT
	      END IF
	    END DO
	  ELSE
	    WRITE(6,*)'Warning -- unrecognized species in REVISED_LAMBDAS: SPECIES=',SPECIES
	  END IF
!
	END DO
!
100	CONTINUE
	CLOSE(LUIN)
!
	RETURN
	END
