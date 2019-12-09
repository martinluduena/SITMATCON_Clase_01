MODULE mod_drive_files
contains
! SUBROUTINE open_file(jj,files)                    
! IMPLICIT NONE                                     
! INTEGER,intent(in) ::jj                           
! INTEGER :: ierror 
! INTEGER :: uini
! CHARACTER(len=60),dimension(5),intent(in):: files
 !OPEN(UNIT=123,file=trim(files(jj)),IOSTAT=ierror)
! uini=find_io(20)

! print*, "uini",uini
! OPEN(uini,file=trim(files(jj)),IOSTAT=ierror) 
! END SUBROUTINE open_file                          
                                                   
integer function find_io(start)

!  find an unused unit number for input or output. unit n=start is used
!  if available; otherwise n is incremented until an unused unit is found.
!  unit numbers are limited to the range 1-100; if n reaches 100 the
!  search starts again at 1.

implicit none
integer, intent(in) :: start
logical :: in_use, exists
integer :: n, n0
integer, parameter :: max_unit=99

n0=start
if (n0 <= 1 .or. n0 > max_unit) n0=1
n=n0
in_use=.true.
do while (in_use)
  inquire(n,opened=in_use,exist=exists)
  if (exists) then
    if (.not. in_use) exit
  else
    !FIXME write (unit=string,fmt="(a,i3,a)") "unit number", n, " out of range"
    !call report (string)
  endif
  n=n+1
  if (n > max_unit) n=1
  if (n == n0) then
    !FIXME call report ("no i/o unit available")
  end if
end do
find_io=n

end function find_io

END MODULE mod_drive_files 
                      
                      
