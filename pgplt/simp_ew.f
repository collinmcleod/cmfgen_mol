	SUBROUTINE SIMP_EW(NORM,WAVE,N,LAM_ST,LAM_END,TOLERANCE)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered: 15-Jun-2014: DOING_LINE now initialized on entry into subroutine.
!
	INTEGER N
	REAL*4 NORM(N)
	REAL*4 WAVE(N)
!
	REAL*4 LAM_ST
	REAL*4 LAM_END
	REAL*4 TOLERANCE
!
	REAL*4 LINE_ST,LINE_END
	REAL(KIND=LDP) SUM,LAM
	INTEGER I,IST,IEND
	INTEGER GET_INDX_SP
	EXTERNAL GET_INDX_SP
	LOGICAL DOING_LINE
!
	WRITE(6,*)'N=',N
	WRITE(6,*)'WAVE(1)=',WAVE(1)
	WRITE(6,*)'WAVE(N)=',WAVE(N)
!
	IST=GET_INDX_SP(LAM_ST,WAVE,N)
	IEND=GET_INDX_SP(LAM_END,WAVE,N)
	WRITE(6,*)'IST=',IST,WAVE(IST)
	WRITE(6,*)'IEND=',IEND,WAVE(IEND)
!
	DOING_LINE=.FALSE.
	DO I=IST,IEND-1
	  IF( ABS(NORM(I)-1.0D0) .GT. TOLERANCE)THEN
	     IF(.NOT. DOING_LINE)THEN
	       SUM=0.0D0
	       LAM=0.0D0
	       DOING_LINE=.TRUE.
	       LINE_ST=WAVE(I)
	     END IF
	     SUM=SUM+0.5D0*(WAVE(I+1)-WAVE(I))*(NORM(I)+NORM(I+1)-2.0D0)
	     LAM=LAM+0.5D0*(WAVE(I+1)-WAVE(I))*(WAVE(I)*(NORM(I)-1.0D0) +
	1                                     WAVE(I+1)*(NORM(I+1)-1.0D0) )
	  ELSE IF(DOING_LINE)THEN
	     LAM=LAM/SUM
	     SUM=SUM*1000.0D0			!mill-Ang
!	     WRITE(6,*)'Line EQ is',SUM
!	     WRITE(6,*)'Line centroid is',LAM
!	     WRITE(6,*)'Line start lambda',LINE_ST
!	     WRITE(6,*)'Line end lambda',WAVE(I)
	     IF( ABS(SUM) .GT. 10.0D0)THEN
	       WRITE(16,'(F12.2,2X,3F12.4)')SUM,LAM,LINE_ST,WAVE(I)
	     END IF
	     DOING_LINE=.FALSE.
	  END IF
	END DO
!
	RETURN
	END
