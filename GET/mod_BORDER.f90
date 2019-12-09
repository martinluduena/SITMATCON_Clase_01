!Este Módulo calcula el orden de enlace de cada átomo
MODULE mod_BORDER
contains
SUBROUTINE doBORDER()
implicit none
integer,parameter :: dp=kind(0.0d0)
real(dp),dimension(:),allocatable,save :: OE_Au,OE_S,OE_N,OE_C

real(dp) :: rx,ry,rz !Posiciones de los atomos (auxiliar)
real(dp),dimension(:),allocatable :: rxAu,ryAu,rzAu
real(dp),dimension(:),allocatable :: rxS,ryS,rzS
real(dp),dimension(:),allocatable :: rxN,ryN,rzN
real(dp),dimension(:),allocatable :: rxC,ryC,rzC

integer :: ntot=0             !numero de atomos totales
integer :: natom_Au=0         !numero de atomos de Au
integer :: natom_S=0          !numero de atomos de S
integer :: natom_N=0          !numero de atomos de N
integer :: natom_C=0          !numero de atomos de C
integer :: nerror=0

character(len=3) :: Mm
character(len=4) :: lline
integer :: nline=0
integer :: nframe=0
integer :: nblanco=0
integer :: nconfig=0

integer :: ierror
integer :: i,j,k
!-----------------------------------------------------------------------
!---- Cuento el Numero total de lineas ---------------------------------
DO
READ(123, *, IOSTAT = ierror) 
IF(ierror  > 0) nline=nline + 1   !Error en la lectura
IF(ierror == 0) nline=nline + 1   !No hay error en la lectura
IF(ierror  < 0) EXIT              !Encontro el final del archivo
ENDDO
REWIND(123)
!.......................................................................
READ(123,*) ntot !Numero de atomos totales por config
REWIND(123) 
nframe=nline/(ntot+2)
WRITE(*,*) "Transitorio Nconfig="
READ(*,*) nblanco ; nconfig=nframe - nblanco

ALLOCATE(rxAu(nconfig*ntot),ryAu(nconfig*ntot),rzAu(nconfig*ntot))
ALLOCATE(rxS(nconfig*ntot),ryS(nconfig*ntot),rzS(nconfig*ntot))
ALLOCATE(rxN(nconfig*ntot),ryN(nconfig*ntot),rzN(nconfig*ntot))
ALLOCATE(rxC(nconfig*ntot),ryC(nconfig*ntot),rzC(nconfig*ntot))

DO j=1,nblanco  !Adelanto las Configuraciones que no quiero que entren en el promedio
 READ(123,*, IOSTAT = ierror)                                                        
 READ(123,*, IOSTAT = ierror)                                                        
  DO i=1,ntot                                                                        
  READ(123,*, IOSTAT = ierror)                                                       
  ENDDO                                                                              
ENDDO                                                                                
WRITE(*,*)"Paso el Blanco",nblanco                                                   

!...... Cuenta el Numero de atomos de cada especie .............   
DO 
 READ(123,'(a2,x,3(f10.6,x))',IOSTAT = ierror) Mm ,rx,ry,rz
 IF(ierror == 0)then
  IF(trim(adjustl(Mm))=="Au")then ; natom_Au = natom_Au  + 1 ; rxAu(natom_Au)=rx ;ryAu(natom_Au)=ry ;rzAu(natom_Au)=rz ;endif
  IF(trim(adjustl(Mm))=="S")then  ; natom_S  = natom_S   + 1 ; rxS(natom_S)  =rx ;ryS(natom_S)=ry   ;rzS(natom_S)=rz   ;endif
  IF(trim(adjustl(Mm))=="N")then  ; natom_N  = natom_N   + 1 ; rxN(natom_N)  =rx ;ryN(natom_N)=ry   ;rzN(natom_N)=rz   ;endif
  IF(trim(adjustl(Mm))=="C")then  ; natom_C  = natom_C   + 1 ; rxC(natom_C)  =rx ;ryC(natom_C)=ry   ;rzC(natom_C)=rz   ;endif
 ENDIF
 IF(ierror > 0)then; nerror=nerror + 1 ; ENDIF      !Cuenta posibles errores de lectura.
 IF(ierror < 0)EXIT                                 !Encuentra el final del archivo.
ENDDO 

write(*,*)"nframes",nframe
write(*,*)"nconfig",nconfig
write(*,*)"nerror",nerror
write(*,*)"ntot",ntot        !Numero total de atomos por config
natom_Au=natom_Au/nconfig
natom_S=natom_S/nconfig
natom_N=natom_N/nconfig
natom_C=natom_C/nconfig

write(*,*)"natom_Au",natom_Au
write(*,*)"natom_S",natom_S
write(*,*)"natom_N",natom_N
write(*,*)"natom_C",natom_C

if(nframe<nblanco)then;print*,"nframe<nblanco";stop;endif


!..................... Archivo de SALIDA ..........
open(unit=400,file="BOND_ORDER.dat")
 
ALLOCATE(OE_Au(natom_Au),OE_S(natom_S),OE_N(natom_N),OE_C(natom_C))
REWIND(123)
!..........................................................
DO j=1,nconfig
 call dist_OE(j,nconfig,natom_Au,ntot,rxAu,ryAu,rzAu,OE_Au)
! call dist_r(j,nconfig,natom_S,ntot , rxS, ryS, rzS,OE_S)
! call dist_r(j,nconfig,natom_N,ntot , rxN, ryN, rzN,OE_N)
! call dist_r(j,nconfig,natom_C,ntot , rxC, ryC, rzC,OE_C)  
ENDDO


DO i=1,natom_Au
write(400,*)i,OE_Au(i)/dfloat(nconfig)
ENDDO

close(400)
END SUBROUTINE doBORDER

SUBROUTINE dist_OE(k,nconfig,natom,ntot,rx,ry,rz,OE)
implicit none
integer,parameter  :: dp=kind(0.0d0)
integer,intent(in) :: k,nconfig,natom,ntot
real(dp),dimension(nconfig*ntot),intent(in) :: rx,ry,rz
real(dp),dimension(natom),intent(inout) :: OE                     
real(dp) :: OEaux,rr
real(dp) :: rrx,rry,rrz
integer :: i,j,ii
ii=0                                   
do i=(k-1)*natom+1,k*natom 
 ii=ii+1 
 do j=(k-1)*natom+1,k*natom 
   if(i/=j)then     
   rrx=rx(j)-rx(i)
   rry=ry(j)-ry(i)                     
   rrz=rz(j)-rz(i)                     
   rr=dsqrt(rrx*rrx+rry*rry+rrz*rrz)
   call fun_OE(rr,OEaux)
   OE(ii)=OE(ii)+OEaux
   endif
  enddo
enddo                                
END SUBROUTINE

SUBROUTINE fun_OE(rr,OE)
IMPLICIT NONE
integer,parameter  :: dp=kind(0.0d0)
real(dp) :: rr,rc1,rc2
real(dp) ::OE
real(dp) ::PI,RMID2

rc1=2.90d0
rc2=4.06d0
RMID2=(rc1+rc2)/2.d0
PI=acos(-1.d0)

if(rr.le.rc1)then                        
OE=1.0d0                                      
else
 if(rr.gt.rc1.and.rr.lt.rc2)then                     
 OE=0.5D0-0.5D0*DSIN(PI*(rr-RMID2)/(rc2-rc1))  
 else                                    
 OE=0.0d0
 endif
end if 
END SUBROUTINE fun_OE 

 END MODULE mod_BORDER























