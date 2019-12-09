MODULE weighting_functions 
 contains 
!_-------------------------------------- Weighting_gas_Entropia (Wg_S)
function Wg_S(m,k,h,T,f,p,Z,Y)
implicit none
integer,parameter :: dp=kind(0.d0)
real(dp),parameter :: pi=3.141592654_dp
real(dp) :: m,k,h,T,f,p,Z,Y
real(dp) :: Wg_S
real(dp) :: bT,bT3
real(dp) :: c1,c2
real(dp) :: SHs
 bT=h/dsqrt(2._dp*pi*m*k*T)   !Thermal de Broglie wavelength


 write(*,*)
 write(*,*)'Thermal de Broglie wavelength',bT
 bT3=bT*bT*bT 
 c1=Z/(bT3*f*p)
 c2=(Y*(3._dp*Y-4._dp))/((1._dp-Y)*(1._dp-Y))
 SHs=2.5_dp + log(c1) + c2    !Entropia de Gas de Esferas Duras
 Wg_S=SHs/(3._dp)           !Weighting_gas_Entropia (Wg_S)


!  write(*,*)'ln',log(c1)
!  write(*,*)'3er miembro',c2
  
!  write(*,*)'kb',k,T
end function
!_---------------------------------------  Weighting_solido_Entropia (Wg_S)
function Ws_S(bhc,vnu)
implicit none
integer,parameter :: dp=kind(0.d0)
real(dp) :: bhc,vnu,Ws_S
real(dp) :: arg
 arg=bhc*vnu
 Ws_S= (arg/(dexp(arg)-1._dp))-dlog(1._dp-dexp(-arg))
end function
 
!_---------------------------------------  Weighting_solido_Helmholtz (Ws_A)

function Ws_A(bhc,vnu) !Chequear esto, porque creo que le falta la temperatura...
implicit none
integer,parameter :: dp=kind(0.d0)
real(dp) :: bhc,vnu,Ws_A
real(dp) :: arg
 arg=bhc*vnu
 Ws_A= 0.5_dp*arg + dlog(1._dp-dexp(-arg))
end function

!_---------------------------------------
!_---------------------------------------- Funcion Coeficiente de Difusion Normalizado
function fun_CDN(T,p,N,m,So,kb,v_luz)
integer,parameter :: dp=kind(0.d0) 
integer :: natom
real(dp),parameter :: pi=3.141592654_dp
real(dp) :: fun_CDN,T,p,m,So,v_luz,kb
real(dp) ::c1,c2,c3,c4
 c1 = 2._dp*So/(9._dp*dble(N))
 c2 = dsqrt(pi*kb*T/m)
 c3 = p**(1._dp/3._dp)
 c4 = (6._dp/pi)**(2._dp/3._dp)
 fun_CDN = c1*c2*c3*c4/v_luz                ![Adimensonal]
end function
!_--------------------------------------------
!_---------------------------------------------------------Ec. Estado Carnahan-Starliing
 function fun_ZCS(Y)
 implicit none
 integer,parameter :: dp=kind(0.d0) 
 real(dp) ::fun_ZCS,Y
 fun_ZCS=(1._dp+ Y +Y*Y - Y*Y*Y)/((1._dp-Y)**3) 
 end function
!_-------------------------------------------

 
end MODULE weighting_functions
