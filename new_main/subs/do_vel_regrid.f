	SUBROUTINE DO_VEL_REGRID(POPS,R,V,POP_ATOM,DONE_R_REV,ND,NT)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NM=5
!
	INTEGER ND
	INTEGER NT
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) POP_ATOM(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
!
	LOGICAL DONE_R_REV
!
	REAL(KIND=LDP) CUR_LOG_POP(NT)
	REAL(KIND=LDP) OLD_LOG_POP(NT)
	REAL(KIND=LDP) LOG_ATOM_DEN
!
	REAL(KIND=LDP) XNEW(NM*ND)
	REAL(KIND=LDP) VNEW(NM*ND)
	REAL(KIND=LDP) RNEW(NM*ND)
	REAL(KIND=LDP) REVR(ND)
	REAL(KIND=LDP) REVV(ND)
	REAL(KIND=LDP) REVX(ND)
	REAL(KIND=LDP) ROLD(ND)
!
	REAL(KIND=LDP) TA(ND),RTA(NM*ND)
	REAL(KIND=LDP) TB(ND),RTB(NM*ND)
!
	REAL(KIND=LDP) VLOW
	REAL(KIND=LDP) VHIGH
	REAL(KIND=LDP) dDEN
!
	REAL(KIND=LDP) dV
	REAL(KIND=LDP) dV_MAX
	REAL(KIND=LDP) T1,T2
!
	INTEGER, PARAMETER :: NVAR_MAX=20
	INTEGER IVAR(NVAR_MAX)
	INTEGER NVAR
!
	INTEGER NIB
	INTEGER NOB
	REAL(KIND=LDP) OBND_PARAMS(5)           !Parameters specifying grid placement at outer boundary.
	REAL(KIND=LDP) IBND_PARAMS(5)           !Parameters specifying grid placement at nner boundary.
	CHARACTER(LEN=10) OUT_BND_OPT   !Outer boundary option
	CHARACTER(LEN=10) IN_BND_OPT
	CHARACTER(LEN=6)  KEY
!
	INTEGER NI
	INTEGER IV1,IV2
	INTEGER JST,JEND
	INTEGER I,J,K,L
	INTEGER NMAX
	INTEGER LU
	INTEGER, PARAMETER :: LUER=6
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	WRITE(LUER,*)'Entered DO_VEL_REGRID'
	DONE_R_REV=.FALSE.
	OUT_BND_OPT='DEFAULT'
	IN_BND_OPT='DEFAULT'
	IV1=1; IV2=21
	NMAX=NM*ND
	dV_MAX=500.0_LDP
!
	CALL RD_STORE_DBLE(VLOW,'VLOW',L_TRUE,'Minimum velocity (km/s) of zone to be regridded')
	CALL RD_STORE_DBLE(VHIGH,'VHIGH',L_TRUE,'Maximum velocity (km/s) of zone to be regridded')
	CALL RD_STORE_DBLE(dDEN,'dDEN',L_TRUE,'Maximum fractional change in density')
!
! IVAR determins which variables are used to define the R grid. IVAR is defined using the linearization
! variables, and should be 1 to NT-1. Use MODEL to check linkage between species and linearization
! variable. The first 2 variables are used for diagnostic output. A good choice is IV1=1, and IV2=NHYD+1
!
	NVAR=0
	CALL RD_STORE_INT(NVAR,'NVAR',L_FALSE,'Number of variables used to help constrain new R grid')
	IF(NVAR .EQ. 0)THEN
	  NVAR=2
	  IVAR(1)=IV1; IVAR(2)=IV2
	ELSE
	  IF(NVAR .GT. NVAR_MAX)THEN
	    NVAR=NVAR_MAX
	    WRITE(LUER,*)'Warning: NVAR is limited to ',NVAR_MAX,'variables'
	  END IF
	  DO I=1,MIN(NVAR,NVAR_MAX)
	    WRITE(KEY,'(I3)')I
	    KEY='IV'//ADJUSTL(KEY)
	    CALL RD_STORE_INT(IVAR(I),TRIM(KEY),L_TRUE,'Level variable for constraining R')
	    IF(IVAR(I) .LT. 1 .OR. IVAR(I) .GT. NT-1)THEN
	      WRITE(LUER,*)'Error in DO_VEL_REGRID'
	      WRITE(LUER,*)'Invalid IVAR value: IVAR(I)=',IVAR(I)
	      RETURN
	    END IF
	    IF(I .EQ. 1)IV1=IVAR(1)
	    IF(I .EQ. 2)IV2=IVAR(2)
	  END DO
	END IF
!
	I=10; CALL RD_STORE_NCHAR(OUT_BND_OPT,'OB_OPT',I,L_FALSE,'Outer boundary option: SPECIFY or DEFAULT')
	IF(OUT_BND_OPT .EQ. 'DEFAULT')THEN
	  NOB=2
	ELSE IF(OUT_BND_OPT .EQ. 'SPECIFY')THEN
	  CALL RD_STORE_INT(NOB,'NOB_PARS',L_TRUE,'Number of outer boundary parameters')
	  DO I=1,NOB
	    WRITE(KEY,'(I3)')I
	    KEY='OB_P'//ADJUSTL(KEY)
	    CALL RD_STORE_DBLE(OBND_PARAMS(I),TRIM(KEY),L_TRUE,'Parameters for outer boundary condition')
	  END DO
	ELSE
	  WRITE(LUER,*)'Invalid outer boundary option in DO_VEL_REGRID'
	  WRITE(LUER,*)'OUT_BND_OPT = ',TRIM(OUT_BND_OPT)
	  RETURN
	END IF
!
	I=10; CALL RD_STORE_NCHAR(IN_BND_OPT,'IB_OPT',I,L_FALSE,'Inner boundary option: SPECIFY or DEFAULT')
	IF(IN_BND_OPT .EQ. 'DEFAULT')THEN
	  NIB=2
	ELSE IF(IN_BND_OPT .EQ. 'SPECIFY')THEN
	  CALL RD_STORE_INT(NIB,'NIB_PARS',L_TRUE,'Number of inner boundary parameters')
	  DO I=1,NIB
	    WRITE(KEY,'(I3)')I
	    KEY='IB_P'//ADJUSTL(KEY)
	    CALL RD_STORE_DBLE(IBND_PARAMS(I),TRIM(KEY),L_TRUE,'Paremeters for inner boundary condition')
	  END DO
	ELSE
	  WRITE(LUER,*)'Invalid inner boundary option in DO_VEL_REGRID'
	  WRITE(LUER,*)'IN_BND_OPT = ',TRIM(IN_BND_OPT)
	  RETURN
	END IF
!
	CALL GET_LU(LU)
	OPEN(UNIT=LU,FILE='R_REGRIDDING_LOG',STATUS='UNKNOWN',ACTION='WRITE')
!
	TA(1:ND)=POPS(IV1,1:ND)
	TB(1:ND)=POPS(IV2,1:ND)
!
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A)')'! Current grid'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A,7X,A,17X,A,16X,A,4X,2(8X,A,7X,A,X))')'! Index','R','V','dV','N(V1)','dN(V1)/N(V1)','N(V2)','dN(V2)/N(V2)'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(I6,3X,ES12.6,6X,ES12.6)')1,R(1),V(1)
	DO I=2,ND
	  WRITE(LU,'(I6,3X,ES12.6,3(6X,ES12.6,5X,F10.3))')
	1             I,R(I),V(I),V(I-1)-V(I),
	1             TA(I),TA(I)/TA(I-1)-1.0_LDP,
	1             TB(I),TB(I)/TB(I-1)-1.0_LDP
	END DO
!
! Find the interval over which we will adjust the grid.
! Remember that V(1) > V(ND), and R(1) > R(ND).
!
	JST=1
	DO I=2,ND
	  IF(V(I) .LT. VHIGH)THEN
	    JST=I-1
	    EXIT
	  END IF
	END DO
!
	JEND=ND
	DO I=ND-1,1,-1
	  IF(V(I) .GT. VLOW)THEN
	    JEND=I+1
	    EXIT
	  END IF
	END DO
!
! Adjust NI to allow for finer grid at boundaries. We initially make  no special adjustments
! for the boundaries when defining the revised grid. After the grid is defined, we add
! NIB points at the inner boundary (if necessary) and NOB points at the outer boundary.
! The code may not define a good boudary grid if JST=2 and/or JEND=ND-1. In practice, this
! should not occur.
!
! We use a factor of 2 to allow for extra point insertion due to the density criterion.
!
	NI=JEND-JST+1
	IF(JST .EQ. 1)NI=NI-NOB
	IF(JEND .EQ. ND)NI=NI-NIB
	dV=EXP(2.0_LDP*LOG(V(JST)/V(JEND))/(NI-1))
	dV=1.04_LDP
!
	WRITE(LUER,'(A)')' '
	WRITE(LUER,'(2(A,I4,14X))')    '    JST=',   JST,'   JEND=',JEND
	WRITE(LUER,'(2(A,F10.4,8X))') ' V(JST)=',V(JST),'V(JEND)=',V(JEND)
	WRITE(LUER,'(2(A,ES14.8,4X))')' R(JST)=',R(JST),'R(JEND)=',R(JEND)
	WRITE(LUER,'(A)')' '
!
100	CONTINUE
	WRITE(LUER,*)' Fractional change in Velocity (=dV) is',dV
!
! Define the new velocity grid. We limit the change in V to dV, or to a
! size such that the % change an important population is smaller than dDEN.
!
! Set V(L) using only dV. We stop when close enough to the end of the
! revised grid.
!
	L=1
	VNEW(1)=V(JST)
	CUR_LOG_POP=LOG(POPS(:,JST))
	DO WHILE (VNEW(L) .GT. V(JEND))
 	  L=L+1
	  IF(L .GT. NMAX)THEN
	    WRITE(LUER,*)'Error in DO_VEL_REGRID -- NMAX not sufficently large'
	    WRITE(LUER,*)'Adjusting dV and starting again'
	    WRITE(LUER,*)'L=',L
	    dV=dV*1.02_LDP
	    GOTO 100
	  END IF
!
	  VNEW(L)=VNEW(L-1)/dV
	  IF(VNEW(L)-0.2_LDP*(VNEW(L-1)-VNEW(L)) .LE. V(JEND))THEN
	    VNEW(L)=V(JEND)
	    EXIT
	  END IF
	  IF(VNEW(L-1)-VNEW(L) .GE. dV_MAX)THEN
	    VNEW(L)=VNEW(L-1)-dV_MAX
	  END IF
!
! This section checks whether the density criteria is met.
!
	  OLD_LOG_POP=CUR_LOG_POP
1000	  CONTINUE
	  K=JST
	  DO WHILE(VNEW(L) .LT. V(K+1))
	    K=K+1
	  END DO
	  T1=(VNEW(L)-V(K))/(V(K+1)-V(K))
	  LOG_ATOM_DEN=T1*LOG(POP_ATOM(K+1))+(1.0_LDP-T1)*LOG(POP_ATOM(K))
	  T2=-1000.0_LDP
	  DO J=1,NVAR
	    I=IVAR(J)
	    CUR_LOG_POP(I)=T1*LOG(POPS(I,K+1))+(1.0_LDP-T1)*LOG(POPS(I,K))
	    IF(CUR_LOG_POP(I) .GT. LOG_ATOM_DEN-12.0_LDP)THEN
	      T2=MAX(T2,ABS(CUR_LOG_POP(I)-OLD_LOG_POP(I)))
	      WRITE(70,*)I,L,LOG_ATOM_DEN,CUR_LOG_POP(I)
	    END IF
	  END DO
!
! We don't use simple ineterpolation to estimate a revised V, since the
! increment in velocity may take us from an interval where the density
! varies slowly to an interval where it varies rapidly.
!
	  IF(T2 .GT. dDEN)THEN
	    VNEW(L)=VNEW(L-1)+0.95_LDP*(VNEW(L)-VNEW(L-1))
	    GOTO 1000
	  END IF
	END DO
!
! For best regidding results, it would be best if these two numbers
! were similar.
!
! This section is only done for testing and verification purposes.
!
	WRITE(LUER,*)'Number of points in old interval inclusive was',JEND-JST+1
	WRITE(LUER,*)' Number of points in new interval inclusive is ',L
!
	TA=0.0_LDP; TB=0.0_LDP
	TA(1:ND)=LOG(POPS(IV1,1:ND))
	TB(1:ND)=LOG(POPS(IV2,1:ND))
	CALL LIN_INTERP(VNEW,RTA,L,V,TA,ND)
	CALL LIN_INTERP(VNEW,RTB,L,V,TB,ND)
	CALL LIN_INTERP(VNEW,RNEW,L,V,R,ND)
	RTA=EXP(RTA); RTB=EXP(RTB)
!
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A)')'! Revised grid --- regridded interval only.'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A,7X,A,17X,A,16X,A,4X,2(8X,A,7X,A,X))')'! Index','R','V','dV','N(V1)','dN(V1)/N(V1)','N(V2)','dN(V2)/N(V2)'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(I6,3X,ES12.6,6X,ES12.6)')1,RNEW(1),VNEW(1)
	DO I=2,L
	  WRITE(LU,'(I6,3X,ES12.6,3(6X,ES12.6,5X,F10.3))')
	1            I,RNEW(I),VNEW(I),VNEW(I-1)-VNEW(I),
	1             RTA(I),RTA(I)/RTA(I-1)-1.0_LDP,
	1             RTB(I),RTB(I)/RTB(I-1)-1.0_LDP
	END DO
!
! Adjust the number of points so that it matches the number of points in the original
! interval.
!
	DO I=1,NI
	  REVX(I)=1.0_LDP+FLOAT((I-1)*(L-1))/(NI-1)
	END DO
	DO I=1,L
	  XNEW(I)=I
	END DO
	CALL LIN_INTERP(REVX,REVV,NI,XNEW,VNEW,L)
	WRITE(LUER,*)'Computed revised R grid'
	CALL LIN_INTERP(REVV,REVR,NI,V,R,ND)
!
	ROLD=R
	IF(JST .EQ. 1 .AND. JEND .EQ. ND)THEN
	  R(2+NOB:ND-NIB-1)=REVR(2:NI-1)
	ELSE IF(JST .EQ. 1)THEN
	  R(2+NOB:JEND)=REVR(2:NI)
	ELSE IF(JEND .EQ. ND)THEN
	  R(JST:ND-NIB-1)=REVR(1:NI-1)
	ELSE
	  R(JST:JEND)=REVR(1:NI)
	END IF
	R(1)=ROLD(1); R(ND)=ROLD(ND)
!
! If necessary, refine the grid at the boundaries.
!
	IF(JST .EQ. 1)THEN
	  IF(OUT_BND_OPT .EQ. 'DEFAULT')THEN
	    R(2)=R(1)-0.02_LDP*(REVR(1)-REVR(2))
	    R(3)=R(1)-0.30_LDP*(REVR(1)-REVR(2))
	  ELSE IF(OUT_BND_OPT .EQ. 'SPECIFY')THEN
	    DO J=1,NOB
	      R(J+1)=R(1)-(REVR(1)-REVR(2))/OBND_PARAMS(J)
	    END DO
	  END IF
	END IF
!
	IF(JEND  .EQ. ND)THEN
	  IF(IN_BND_OPT .EQ. 'DEFAULT')THEN
	    R(ND-2)=R(ND)+0.4_LDP*(REVR(NI-1)-REVR(NI))
	    R(ND-1)=R(ND)+0.1_LDP*(REVR(NI-1)-REVR(NI))
	  ELSE IF(IN_BND_OPT .EQ. 'SPECIFY')THEN
	    DO J=1,NOB
	      R(ND-J)=R(ND)+(REVR(NI-1)-REVR(NI))/IBND_PARAMS(J)
	    END DO
	  ELSE
	    WRITE(LUER,*)'Invalid outer boundary option in DO_VEL_REGRID'
	    WRITE(LUER,*)'IN_BND_OPT = ',TRIM(IN_BND_OPT)
	    RETURN
	  END IF
	END IF
!
! For testing purposes.
!	
	TA(1:ND)=LOG(POPS(IV1,1:ND))
	TB(1:ND)=LOG(POPS(IV2,1:ND))
	CALL LIN_INTERP(R,RTA,ND,ROLD,TA,ND)
	CALL LIN_INTERP(R,RTB,ND,ROLD,TB,ND)
	CALL LIN_INTERP(R,VNEW,ND,ROLD,V,ND)
	RTA=EXP(RTA); RTB=EXP(RTB)
!
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A)')'! Revised R grid'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(A,7X,A,17X,A,16X,A,4X,2(8X,A,7X,A,X))')'! Index','R','V','dV','N(V1)','dN(V1)/N(V1)','N(V2)','dN(V2)/N(V2)'
	WRITE(LU,'(A)')'!'
	WRITE(LU,'(I6,3X,ES12.6,6X,ES12.6)')1,R(1),VNEW(1)
	DO I=2,ND
	  WRITE(LU,'(I6,3X,ES12.6,3(6X,ES12.6,5X,F10.3))')
	1            I,R(I),VNEW(I),VNEW(I-1)-VNEW(I),
	1             RTA(I),RTA(I)/RTA(I-1)-1.0_LDP,
	1             RTB(I),RTB(I)/RTB(I-1)-1.0_LDP
	END DO
	CLOSE(UNIT=LU)
!
	OPEN(UNIT=LU,FILE='NEW_R_GRID',STATUS='UNKNOWN')
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A,T40,A)')' 24-Feb-2004','!Format date'
	WRITE(LU,'(A)')' '
	WRITE(LU,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4,5X,I4)')R(ND),0.0D0,1,ND
	DO I=1,ND
          WRITE(LU,'(A)')' '
          WRITE(LU,'(1X,1P,E15.7,6E15.5,2X,I4,A1)')R(I),1.0_LDP,1.0_LDP,1.0_LDP,1.0_LDP,VNEW(I),1.0_LDP,I
          WRITE(LU,'(F7.1)')1.0_LDP
	END DO
	CLOSE(UNIT=LU)
	DONE_R_REV=.TRUE.
!
	RETURN
	END
