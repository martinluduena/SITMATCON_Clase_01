
! This is a src file for GeMS program and is under subversion control.
! Here there are info about the last (relevant) commit of this file:
! $Rev: 43 $ $Author: alexis $ 
! $Date: 2013-03-07 11:50:21 -0300 (jue 07 de mar de 2013) $
 
module strings
  use constants, only:zmax,z_sym,dp 

  implicit none

  private

  interface operator(.ich.)
    module procedure int2char,float2char
  end interface

  public    :: upcase,locase,operator(.ich.),int2char0
  public    :: inq_z,contar_cifras

  character(len=26), parameter ::                                        &
    upper_case="ABCDEFGHIJKLMNOPQRSTUVWXYZ",                           &
    lower_case="abcdefghijklmnopqrstuvwxyz"

 contains
 
  function contar_cifras(num) result(c)
   integer,intent(in)          :: num
   integer                     :: a,c
   a=num
   c=0
   if(a==0) then
     c=1
     return
   endif
   do while (a>0)
     a=int(a/10)
     c=c+1
   enddo
  end function  
 
  function int2char(num)
   ! Convert integer to string
   integer,intent(in)           :: num
   character(20)                :: int2char
   write(int2char,'(I0)') num 
  end function  

  function float2char(num)
   ! Convert float to string
   real(dp),intent(in)          :: num
   character(20)                :: float2char
   write(float2char,'(f10.6)') num 
  end function   

  subroutine upcase(word)
  ! Change the word to upercase
  character(len=*), intent(inout) :: word
  integer :: i,k

  do i=1,len(word)
    k=index(lower_case,word(i:i))
    if (k .ne. 0) word(i:i)=upper_case(k:k)
  end do

  end subroutine upcase

  subroutine locase(word)
  ! Change the word to lowercase
  character(len=*), intent(inout) :: word
  integer :: i,k

  do i=1,len(word)
    k=index(upper_case,word(i:i))
    if (k .ne. 0) word(i:i)=lower_case(k:k)
  end do
  end subroutine locase
 

  function inq_z(sym)
    ! Devuelve el id del elemento correspondiente al simbolo sym
    integer                   :: inq_z
    character(2)              :: sym,l1,l2

    l1 = sym
    call locase(l1)
    do inq_z = 0, zmax
      l2 = z_sym(inq_z)
      call locase(l2)
      if (l1==l2)  return
    enddo  
    inq_z = -1
  end function inq_z


  function int2char0(num1,num2) result(c)
   ! Like integer_char but add ceros on the left until 
   ! fill the same digit number of num2.
   integer,intent(in)          :: num1,num2
   integer                     :: a,i
   character(20)               :: c
   a=contar_cifras(num2)-contar_cifras(num1)
   c = ''
   do i=1,a
     c(i:i) = '0' 
   enddo
   write(c,'(a)') trim(c) // trim(.ich.num1)
  end function  
 


end module strings 


