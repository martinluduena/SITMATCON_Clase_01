MODULE Simpson
contains
!_------------------------------------ Cuadratura Simpson
subroutine Cuad_Simps(Ndat,h,g,SUM)
implicit none
integer,parameter :: dp=kind(0.d0)
integer  :: j,Ndat
real(dp),dimension(1:Ndat) :: g !funcion a integrar
real(dp) ::h                    !paso de cuadratura
real(dp) :: SUM,FAC
if(mod(real(Ndat),2.)/=0) Ndat=Ndat-1
SUM=g(2)
FAC=2._dp
do j=3,Ndat-1
IF(FAC==2._dp)THEN
FAC=4._dp
ELSE
FAC=2._dp
ENDIF
SUM=SUM+FAC*g(j)
enddo
SUM=SUM+g(Ndat)
SUM=SUM*h/3._dp
end subroutine
!_-----------------------------------------
end MODULE Simpson
