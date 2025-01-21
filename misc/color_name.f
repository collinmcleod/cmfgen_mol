!
! Returned size is LEN of function in calling program.
! BLUE_PEN etc are 5 characters long
!
	function color_lev_name(trans_name,level) result(res)
	USE MOD_COLOR_PEN_DEF
	character(len=*) res
	character(len=*), intent(in) :: trans_name
	character(len=*), intent(in) :: level 
	character(len=100) string
	integer ist,iend
	integer lev_len
!
	string=trans_name
	lev_len=len_trim(level)
	ist=index(string,trim(level))
	iend=index(string,'-')
	if(ist .lt.iend)then
	  iend=index(string,'-')
	  string=BLUE_PEN//string(1:iend-1)//DEF_PEN//string(iend:)
	else if(ist .ne. 0)then
	  ist=iend
	  string=string(1:ist)//BLUE_PEN//string(ist+1:)
	  iend=len_trim(string)
	  string=string(1:iend)//DEF_PEN
	end if
	res = trim(string)
	end function
