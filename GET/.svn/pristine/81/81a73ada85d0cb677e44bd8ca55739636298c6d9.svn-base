!Este Módulo genera los "feff.inp"
!que son las tarjetas de entrada para 
!el codigo FEFF6 que calcula los expectros 
!EXAFS. 
MODULE mod_FEFFinp
USE mod_ran2, only : ran2
contains
SUBROUTINE doFEFFinp()
implicit none
integer,parameter :: dp=kind(0.0d0)
real,dimension(:),allocatable :: rxyz,rrxyz,rr
real,dimension(3) :: rcm      !posicion centro de masa
integer :: ntot=0             !numero de atomos totales
integer :: natom=0            !numero de atomos metálicos
character(len=2) :: Mm
character(len=4) :: lline
integer :: nline=0
integer :: nframe
integer :: nblanco=0
!.........................
integer  :: nshell=1
real,dimension(:),allocatable :: rshell
integer,dimension(:),allocatable ::abyshell
integer,dimension(:,:),allocatable :: indexshell
integer,dimension(:,:),allocatable:: phabs
!.........................
character(len=5) :: chi,chj
!.........................
integer :: ierror
!-------------------------
integer :: i,ii,j,jj,k  !Contadores
!...................................
integer :: idum
!...................................
logical :: flag=.true.
!-----------------------------------------------------------------------
!---- Cuento el Numero total de lineas ---------------------------------
DO
READ(123,*, IOSTAT = ierror) 
IF(ierror  > 0) nline=nline + 1   !Error en la lectura
IF(ierror == 0) nline=nline+1     !No hay error en la lectura
IF(ierror  < 0) EXIT              !Encontro el final del archivo
ENDDO
REWIND(123)
!.........................................................
READ(123,*) ntot !Numero de atomos totales
READ(123,*)
nframe=nline/(ntot+2)
!...... Cuenta el Numero de atomos Metalicos .............
DO 
READ(123,*) Mm;IF(Mm/="Au") EXIT;IF(Mm=="Au") natom=natom +1
ENDDO 
write(*,*)"nlineas",nline
write(*,*)"nframes",nframe
write(*,*)"natom",natom
if(nframe<nblanco)then;print*,"nframe<nblanco";stop;endif
!........................................................
!********************************************************
write(*,*)" ............................. " 
write(*,*) "nshell ? " ; read(*,*) nshell
!....... Allocateo posiciones,distancias,shells......... 
ALLOCATE(rxyz(3*natom));ALLOCATE(rrxyz(3*natom))
ALLOCATE(rr(natom))  
ALLOCATE(rshell(0:nshell))
ALLOCATE(abyshell(0:nshell))
ALLOCATE(indexshell(0:nshell,natom))
ALLOCATE(phabs(0:nshell,natom))
!.......................................................
OPEN(301,file="abyshell.dat")
!.......................................................
write(*,*) "Ditancias de las Shells ? "
rshell(0)= 0.0
DO i=1,nshell
read(*,*) rshell(i)
if(i==nshell) write(*,*) "end charging shells"
ENDDO
REWIND(123)
!..........................................................
DO j=1,nblanco !Adelanto los Configuraciones que no quiero que entren en el promedio
 READ(123,*, IOSTAT = ierror) 
 READ(123,*, IOSTAT = ierror) 
  DO i=1,ntot 
  READ(123,*, IOSTAT = ierror) 
  ENDDO
ENDDO  
DO j=nblanco+1,nframe
 READ(123,*) ntot
 READ(123,*)
  DO i=1,3*natom,3
  READ(123,*) Mm,rxyz(i),rxyz(i+1),rxyz(i+2)   
  ENDDO
  call cm(natom,rxyz,rr,rcm,rrxyz)
  call atom_by_shell(rcm,natom,rr,nshell,rshell,abyshell,indexshell)      
  write(chj,"(i4)")j-nblanco   ! Configuración
  do jj=0,nshell
   write(chi,"(i4)")jj !Shell
   write(301,*)j-nblanco,jj,abyshell(jj) 
    OPEN(300,file="feff_inputs/feff.inp_"//trim(adjustl(chj))//"_"//trim(adjustl(chi)))         
    call choose_atom_phabs(jj,idum,nshell,natom,abyshell,indexshell,phabs)
    call write_feffinp()
     do ii=0,nshell
      do i=1,abyshell(ii)

       if(flag)then
       flag=.false. 
         write(300,'(3(x,f9.4),i3,x,a10,5x,i2,2x,f9.4)') &
             rrxyz(3*indexshell(ii,i)-2), &
             rrxyz(3*indexshell(ii,i)-1), &
             rrxyz(3*indexshell(ii,i)), &
             phabs(ii,i),"shell",ii,rshell(ii) 
       else
        write(300,'(3(x,f9.4),i3)') &
             rrxyz(3*indexshell(ii,i)-2), &                  
             rrxyz(3*indexshell(ii,i)-1), &                
             rrxyz(3*indexshell(ii,i)), &                
             phabs(ii,i)          
       endif

      enddo
      flag=.true.
     enddo
    close(300)
  enddo
  DO i=natom+1,ntot !Adelanto lo que no me interesa
  READ(123,*, IOSTAT = ierror) 
  ENDDO
ENDDO 

1050 FORMAT(a2,x,3(F10.6,x))
DEALLOCATE(rxyz)
DEALLOCATE(rr)
DEALLOCATE(rshell)          
DEALLOCATE(abyshell)
DEALLOCATE(indexshell)
ENDSUBROUTINE doFEFFinp
SUBROUTINE choose_atom_phabs(ii,idum,nshell,natom,abyshell,indexshell,phabs)
implicit none                                                                             
integer,parameter  :: dp=kind(0.0d0)
integer :: ii
integer :: natom,nshell
integer,dimension(0:nshell,natom),intent(inout) :: indexshell
integer,dimension(0:nshell,natom),intent(inout) :: phabs
integer,dimension(0:nshell),intent(in) :: abyshell
!---------------
integer  ::idum
integer  :: i,x
!---------------
phabs=1
x=int(ran2(idum)*abyshell(ii))+1
phabs(ii,x)=0
ENDSUBROUTINE choose_atom_phabs
SUBROUTINE atom_by_shell(rcm,natom,rr,nshell,rshell,abyshell,indexshell)
implicit none                                                                             
integer,parameter  :: dp=kind(0.0d0)
integer :: natom,nshell
real,dimension(3),intent(in) :: rcm 
real,dimension(natom),intent(in)  :: rr
real,dimension(0:nshell),intent(in)  :: rshell
integer,dimension(0:nshell),intent(inout) :: abyshell
integer,dimension(0:nshell,natom),intent(out) :: indexshell
integer :: i,ii
abyshell=0
indexshell=0
do i=1,natom
 do ii=0,nshell-1

 if(rr(i)==rshell(ii))then
  abyshell(ii)=abyshell(ii)+1
  indexshell(ii,abyshell(ii)) = i
 else
  if(rr(i)>rshell(ii).and.rr(i)<=rshell(ii+1))then
  abyshell(ii+1) = abyshell(ii+1) + 1
  indexshell(ii+1,abyshell(ii+1)) = i
  endif 
endif
 enddo
enddo
ENDSUBROUTINE atom_by_shell
SUBROUTINE cm(natom,rxyz,rr,rcm,rrxyz)
 implicit none 
 integer,parameter  :: dp=kind(0.0d0)
 integer,intent(in) :: natom
 real,dimension(3*natom),intent(in) :: rxyz 
 real,dimension(3*natom),intent(out) :: rrxyz 
 real,dimension(natom),intent(out)  :: rr
 real,dimension(3),intent(inout) :: rcm 
 real  :: xcm,ycm,zcm
 real :: rrx,rry,rrz
 integer :: i,j

 real :: rclose
 real,dimension(3) :: rxyzclose

 rclose=100.0
  xcm=0.0;ycm=0.0;zcm=0.0   
  do i=1,3*natom,3
   xcm=xcm+rxyz(i)
   ycm=ycm+rxyz(i+1)
   zcm=zcm+rxyz(i+2)
  enddo
  rcm(1)=xcm/real(natom)
  rcm(2)=ycm/real(natom)
  rcm(3)=zcm/real(natom)    
  j=0
  do i=1,3*natom,3
  j=j+1
  rrx=rxyz(i)-rcm(1)
  rry=rxyz(i+1)-rcm(2)
  rrz=rxyz(i+2)-rcm(3)
  rrxyz(i)=rrx
  rrxyz(i+1)=rry
  rrxyz(i+2)=rrz

  rr(j)=sqrt(rrx*rrx+rry*rry+rrz*rrz)
  if(rr(j)<rclose)then
    rclose=rr(j)
    rxyzclose(1)=rrxyz(i)
    rxyzclose(2)=rrxyz(i+1)
    rxyzclose(3)=rrxyz(i+2)
  endif
  enddo

  j=0
  do i=1,3*natom,3
  j=j+1
  rrxyz(i)=rrxyz(i)-rxyzclose(1)
  rrxyz(i+1)=rrxyz(i+1)-rxyzclose(2)
  rrxyz(i+2)=rrxyz(i+2)-rxyzclose(3)

  rrx=rrxyz(i)
  rry=rrxyz(i+1)
  rrz=rrxyz(i+2)

  rr(j)=sqrt(rrx*rrx+rry*rry+rrz*rrz)

  !write(147,*)"rr",j,rr(j)
  enddo
  rclose=100.0
 END SUBROUTINE cm
 SUBROUTINE write_feffinp()
 write(300,*) " * Sample input file "                                          
 write(300,*) "TITLE Au,Amina"
 write(300,*) " "
 write(300,*) "DEBYE  20  160  Au at 90K, Debye temp 16K (Ashcroft & Mermin)"
 write(300,*) " "                                                              
 write(300,*) "POTENTIALS"                                                     
 write(300,*) " 0 79 "                                                         
 write(300,*) " 1 79 "                                                         
 write(300,*) " "                                                              
 write(300,*) "ATOMS" 
 ENDSUBROUTINE write_feffinp   
 END MODULE mod_FEFFinp
