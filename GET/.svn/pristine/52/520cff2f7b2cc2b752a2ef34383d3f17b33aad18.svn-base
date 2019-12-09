 MODULE biseccion
 contains
 
!_---------------------------------------------------
 subroutine fraccion_gas(CDN,f)
 
implicit none
!---------------------------------------------------
integer,parameter :: dp=kind(0.d0)
real(dp) :: CDN
real(dp) :: f,d_f
!---------------------------------------------------
integer :: iter
real(dp),parameter :: tol=1.e-8_dp
real(dp) :: fold
!---------------------------------------------------   
f=0.00001_dp      
fold=fun(CDN,f)
d_f=0.0001_dp

iter=0
do
iter=iter+1
f=f+d_f

if(fold*fun(CDN,f)<0.)then
f=f-d_f
d_f=0.5_dp*d_f
endif

if(abs(d_f)<tol)EXIT
enddo

end subroutine

!----------------------------------------------------------------
function fun(CDN,f)     
implicit none
integer,parameter :: dp=kind(0.d0)

real(dp) :: CDN,f
real(dp) :: fun
real(dp) :: c1,c2,c3,c4

 c1= 2._dp*CDN**(-4.5)
 c2=-6._dp*CDN**(-3)
 c3=-1._dp*CDN**(-1.5)
 c4=-6._dp*c3

fun=c1*(f**7.5)+c2*(f**5)+c3*(f**3.5)+c4*(f**2.5)+2._dp*f-2._dp
end function 

!---------------------------------------------------------------
end MODULE biseccion
