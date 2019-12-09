!Este Módulo calcula el orden de enlace de cada átomo
MODULE mod_PDF
contains
SUBROUTINE doPDF(unit_file,name_files)
implicit none
integer,parameter :: dp=kind(0.0d0)

CHARACTER(len=*),dimension(5),intent(in) :: name_files 
integer,intent(in) :: unit_file

integer,dimension(:,:),allocatable,save :: Hist_Au,Hist_S,Hist_N,Hist_C
real(dp),dimension(:),allocatable,save :: PDF_Au,PDF_S,PDF_N,PDF_C
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
integer :: ierror
integer :: i,j,k   !..........Contadores

character(len=3) :: Mm
character(len=4) :: lline
integer :: nline=0
integer :: nframe=0
integer :: nblanco=0
integer :: nconfig=0
!.....................................................
integer :: nbin=1000

real(dp) ::rhist_Au=150.d0         !radio máximo para Au
real(dp) ::rhist_N =150.d0         !radio máximo para N
real(dp) ::rhist_S =150.d0         !radio máximo para S
real(dp) ::rhist_C =150.d0         !radio máximo para C
real(dp) :: delr_Au
real(dp) :: delr_N
real(dp) :: delr_S
real(dp) :: delr_C

delr_Au=rhist_Au/real(nbin)
delr_N=rhist_N/real(nbin)
delr_S=rhist_S/real(nbin)
delr_C=rhist_C/real(nbin)  

 CALL count_lines(nline)   
!.......................................................................
READ(123,*) ntot !Numero de atomos totales por config
REWIND(123) 
nframe=nline/(ntot+2)
WRITE(*,*)"Transitorio Nconfig="
READ(*,*) nblanco ; nconfig=nframe - nblanco

ALLOCATE(rxAu(nconfig*ntot),ryAu(nconfig*ntot),rzAu(nconfig*ntot))
ALLOCATE(rxS(nconfig*ntot),ryS(nconfig*ntot),rzS(nconfig*ntot))
ALLOCATE(rxN(nconfig*ntot),ryN(nconfig*ntot),rzN(nconfig*ntot))
ALLOCATE(rxC(nconfig*ntot),ryC(nconfig*ntot),rzC(nconfig*ntot))

CALL jump_lines(nblanco,ntot)

WRITE(*,*)"Paso el Blanco",nblanco                                                   
!...... Cuenta el Numero de atomos y leo sus posiciones........   
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
natom_Au=natom_Au/nconfig  ; write(*,*)"natom_Au", natom_Au
natom_S=natom_S/nconfig    ; write(*,*)"natom_S" , natom_S
natom_N=natom_N/nconfig    ; write(*,*)"natom_N" , natom_N
natom_C=natom_C/nconfig    ; write(*,*)"natom_C" , natom_C

if(nframe<nblanco)then;print*,"nframe<nblanco";stop;endif

!........ Allocateo distancias e histogramas
ALLOCATE(Hist_Au(natom_Au,nbin),Hist_S(natom_S,nbin),Hist_N(natom_N,nbin),Hist_C(natom_C,nbin))
ALLOCATE(PDF_Au(nbin),PDF_S(nbin),PDF_N(nbin),PDF_C(nbin))
REWIND(123)
!..........................................................
DO j=1,nconfig
call dist_p(j,nconfig,natom_Au,ntot,nbin,rxAu,ryAu,rzAu,delr_Au,rhist_Au,Hist_Au)
call dist_p(j,nconfig,natom_S,ntot,nbin , rxS,ryS , rzS,delr_S,rhist_S, Hist_S) 
call dist_p(j,nconfig,natom_N,ntot,nbin , rxN,ryN , rzN,delr_N,rhist_N, Hist_N) 
call dist_p(j,nconfig,natom_C,ntot,nbin , rxC,ryC , rzC,delr_C,rhist_C, Hist_C) 
ENDDO
!............ Escritura del Histograma .... PDF
OPEN(301,file=name_files(unit_file)//'/PDF_Au.'//name_files(unit_file)//'.dat')
OPEN(302,file=name_files(unit_file)//'/PDF_S.'//name_files(unit_file)//'.dat')
OPEN(303,file=name_files(unit_file)//'/PDF_N.'//name_files(unit_file)//'.dat')
OPEN(304,file=name_files(unit_file)//'/PDF_C.'//name_files(unit_file)//'.dat')


 DO i=1,natom_Au
 CALL prom_over_i(natom_Au,nbin,Hist_Au,PDF_Au)
 ENDDO 
 
 DO i=1,natom_S
 CALL prom_over_i(natom_S,nbin,Hist_S,PDF_S)
 ENDDO 
 
 DO i=1,natom_N
 CALL prom_over_i(natom_N,nbin,Hist_N,PDF_N)
 ENDDO 
 
 DO i=1,natom_C
 CALL prom_over_i(natom_C,nbin,Hist_C,PDF_C)
 ENDDO

 DO j=1,nbin
  write(301,*) dfloat(j)*delr_Au,PDF_Au(j)/nconfig
  write(302,*) dfloat(j)*delr_S,PDF_S(j)/nconfig
  write(303,*) dfloat(j)*delr_N,PDF_N(j)/nconfig
  write(304,*) dfloat(j)*delr_C,PDF_C(j)/nconfig
 ENDDO
!.............................................
DEALLOCATE(Hist_Au,Hist_S,Hist_N,Hist_C)
CLOSE(301);CLOSE(302);CLOSE(303);CLOSE(304)  
END SUBROUTINE doPDF

SUBROUTINE dist_p(k,nconfig,natom,ntot,nbin,rx,ry,rz,delr,rhist,Hist)
implicit none
integer,parameter     :: dp=kind(0.0d0)
integer,intent(in)    :: k,nconfig,natom,ntot,nbin
real(dp),dimension(nconfig*ntot),intent(in) :: rx,ry,rz
real(dp) :: rr                  
real(dp),intent(in)   :: rhist
real(dp) :: rrx,rry,rrz
integer,intent(inout) :: Hist(natom,nbin)
real(dp),intent(in)  :: delr
integer :: inb
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
   inb=int(rr/delr)+1
   Hist(ii,inb)=Hist(ii,inb)+1
   endif
  enddo
enddo   
END SUBROUTINE

SUBROUTINE prom_over_i(max_i,max_j,vector,vec_sum) 
integer,parameter     :: dp=kind(0.0d0)
integer,intent(in) :: max_i,max_j
integer :: i,j
integer,dimension(max_i,max_j) :: vector
real(dp),dimension(max_j) :: vec_sum
real(dp) :: suma

do j=1,max_j
suma=0.d0
 do i=1,max_i
 suma=suma + vector(i,j)
 enddo
vec_sum(j)= suma/dfloat(max_i)
enddo
 ENDSUBROUTINE

SUBROUTINE count_lines(nline)
IMPLICIT NONE
INTEGER :: ierror,nline
DO
READ(123, *, IOSTAT = ierror) 
IF(ierror  > 0) nline=nline + 1   !Error en la lectura
IF(ierror == 0) nline=nline+1     !No hay error en la lectura
IF(ierror  < 0) EXIT              !Encontro el final del archivo
ENDDO
REWIND(123)
END SUBROUTINE count_lines

SUBROUTINE jump_lines(nblanco,ntot)
IMPLICIT NONE
INTEGER :: i,j,ierror
INTEGER :: nblanco,ntot
DO j=1,nblanco  !Adelanto las Configuraciones que no quiero que entren en el promedio
 READ(123,*, IOSTAT = ierror)                                                        
 READ(123,*, IOSTAT = ierror)                                                        
  DO i=1,ntot                                                                        
  READ(123,*, IOSTAT = ierror)                                                       
  ENDDO                                                                              
ENDDO                                                                                
WRITE(*,*)"Paso el Blanco",nblanco 
END SUBROUTINE jump_lines

END MODULE mod_PDF

