!Este Módulo calcula la función de distribución radial
!respecto del centro de masas.
!Esto es lo que uno debe esperar como resultado.

!   6 ++-------+--------+--------+--------+-------+--------+--------+-------++
!    +        +        +        +       './RDF_AuCoor-S147-300K.dat' ****** +
!    |                                                                      |
!  5 ++                               *                                    ++
!    |                                *                                     |
!    |                                *                                     |
!  4 ++                               *                                    ++
!    |                                *                                     |
!    |                                *                                     |
!  3 ++                               *                                    ++
!    |                             ****                                     |
!    |                    *        ****                                     |
!    |                    *        ****  *                                  |
!  2 ++                   * **     ***** *                                 ++
!    |             *      * **    ****** *                                  |
!    |             *     ******   **  *****                                 |
!  1 ++         * ***    *******  *    ****                                ++
!    |  *      ******    ******** *    *****                                |
!    +  *     +******* +*   *  ****    * ***     +        +        +        +
!  0 *************-******------***-------+-******************************--++
!    0        2        4        6        8       10       12       14       16

MODULE mod_RDF
USE mod_count_line
integer,parameter,private :: dp=kind(0.0d0)

contains
SUBROUTINE doRDF(unit_file,name_files)
implicit none
integer,dimension(:),allocatable,save :: Hist_Au,Hist_S,Hist_N,Hist_C,Hist_Cu
integer,dimension(:),allocatable,save :: Hist_Ag,Hist_Ni,Hist_Co
integer,intent(in) :: unit_file
CHARACTER(len=*),dimension(5),intent(in) :: name_files ! poniendo len=* entra
                                                       !el character entra con la longitud 
                                                       !trimiada que le entre
                                                       !como variable
                                                       !dummy(entra por la
                                                       !subroutina)
                                                    


real(dp),dimension(:),allocatable :: rrAu,rrS,rrN,rrC,rrCu
real(dp),dimension(:),allocatable :: rrAg,rrNi,rrCo
real(dp) :: rx,ry,rz !Posiciones de los atomos (auxiliar)
real(dp),dimension(:),allocatable :: rxAu,ryAu,rzAu
real(dp),dimension(:),allocatable :: rxS,ryS,rzS
real(dp),dimension(:),allocatable :: rxN,ryN,rzN
real(dp),dimension(:),allocatable :: rxC,ryC,rzC
real(dp),dimension(:),allocatable :: rxCu,ryCu,rzCu
real(dp),dimension(:),allocatable :: rxAg,ryAg,rzAg
real(dp),dimension(:),allocatable :: rxCo,ryCo,rzCo
real(dp),dimension(:),allocatable :: rxNi,ryNi,rzNi
real(dp),dimension(:),allocatable :: rxtot,rytot,rztot

real(dp),dimension(3) :: rcm  !posicion centro de masa

integer :: ntot=0             !numero de atomos totales
integer :: natom=0
integer :: natom_Au=0         !numero de atomos de Au
integer :: natom_S=0          !numero de atomos de S
integer :: natom_N=0          !numero de atomos de N
integer :: natom_C=0          !numero de atomos de C
integer :: natom_Cu=0         !numero de atomos de Cu
integer :: natom_Ag=0         !numero de atomos de Ag
integer :: natom_Co=0         !numero de atomos de Co
integer :: natom_Ni=0         !numero de atomos de Ni
integer :: nerror=0
integer :: ierror
integer :: i,j,k !............Contadores  

character(len=3) :: Mm  !Identidad de los elementos
character(len=4) :: lline
!integer :: nline=0,nbas=0
integer :: nframe=0
integer :: nblanco=0
integer :: nconfig=0
!.................. RDF .......................................
integer :: nbin=200               !nbin (barras del histograma)

real(dp) ::rhist_Au=8.d0         !radio máximo para Au
real(dp) ::rhist_N =10.d0         !radio máximo para N
real(dp) ::rhist_S =10.d0         !radio máximo para S
real(dp) ::rhist_C =12.d0         !radio máximo para C
real(dp) ::rhist_Cu=8.d0         !radio máximo para Cu
real(dp) ::rhist_Ag=8.d0         !radio máximo para Ag
real(dp) ::rhist_Ni=8.d0         !radio máximo para Ni
real(dp) ::rhist_Co=8.d0         !radio máximo para Co
real(dp) :: delr_Au
real(dp) :: delr_N
real(dp) :: delr_S
real(dp) :: delr_C
real(dp) :: delr_Cu
real(dp) :: delr_Ag
real(dp) :: delr_Ni
real(dp) :: delr_Co

delr_Au=rhist_Au/real(nbin)
delr_N=rhist_N/real(nbin)
delr_S=rhist_S/real(nbin)
delr_C=rhist_C/real(nbin)
delr_Cu=rhist_Cu/real(nbin)
delr_Ag=rhist_Ag/real(nbin)
delr_Ni=rhist_Ni/real(nbin)
delr_Co=rhist_Co/real(nbin)
!................................................................
 CALL count_lines()
!.......................
READ(123,*) ntot !Numero de atomos totales por config
REWIND(123) 
nframe=nline/(ntot+2)
WRITE(*,*)"Transitorio Nconfig= "
READ(*,*) nblanco ; nconfig=nframe - nblanco



ALLOCATE(rxtot(nconfig*ntot),rytot(nconfig*ntot),rztot(nconfig*ntot))

ALLOCATE(rxAu(nconfig*ntot),ryAu(nconfig*ntot),rzAu(nconfig*ntot))
ALLOCATE( rxS(nconfig*ntot), ryS(nconfig*ntot), rzS(nconfig*ntot))
ALLOCATE( rxN(nconfig*ntot), ryN(nconfig*ntot), rzN(nconfig*ntot))
ALLOCATE( rxC(nconfig*ntot), ryC(nconfig*ntot), rzC(nconfig*ntot))
ALLOCATE( rxCu(nconfig*ntot), ryCu(nconfig*ntot), rzCu(nconfig*ntot))
ALLOCATE(rxAg(nconfig*ntot),ryAg(nconfig*ntot),rzAg(nconfig*ntot))
ALLOCATE(rxNi(nconfig*ntot),ryNi(nconfig*ntot),rzNi(nconfig*ntot))
ALLOCATE(rxCo(nconfig*ntot),ryCo(nconfig*ntot),rzCo(nconfig*ntot))

!............................
 CALL jump_lines(nblanco,ntot)
!...... Cuento Numero de atomos y leo sus posiciones .............    

! Error frecuente :
! Subscript #1 of the array HIST has value 9962 which is greater than the upper bound of 300
! Tener MUCHO cuidado con el formato con el que está escrito el archivo.

DO 
! READ(123,'(a2,x,3(e16.9,x))',IOSTAT = ierror) Mm , rx , ry , rz
 READ(123,*,IOSTAT = ierror) Mm , rx , ry , rz

 IF(ierror == 0)then

   IF( (trim(adjustl(Mm))=="Au").or.&
       (trim(adjustl(Mm))=="Cu").or.&
       (trim(adjustl(Mm))=="S").or.&
       (trim(adjustl(Mm))=="N").or.&
       (trim(adjustl(Mm))=="Ag").or.&
       (trim(adjustl(Mm))=="Ni").or.&
       (trim(adjustl(Mm))=="Co").or.&
       (trim(adjustl(Mm))=="C")     )THEN

        WRITE(550,'(a2,x,3(f16.9,x))',IOSTAT = ierror) Mm , rx , ry , rz
        natom = natom + 1 
        rxtot(natom)=rx ; rytot(natom)=ry ; rztot(natom)=rz
   ENDIF  
                                    
                                                                   
  IF(trim(adjustl(Mm))=="Au")then  ; natom_Au = natom_Au  + 1 ; rxAu(natom_Au)=rx ;ryAu(natom_Au)=ry ;rzAu(natom_Au)=rz ;endif
  IF(trim(adjustl(Mm))=="S")then   ; natom_S  = natom_S   + 1 ; rxS(natom_S)  =rx ;ryS(natom_S)=ry   ;rzS(natom_S)=rz   ;endif
  IF(trim(adjustl(Mm))=="N")then   ; natom_N  = natom_N   + 1 ; rxN(natom_N)  =rx ;ryN(natom_N)=ry   ;rzN(natom_N)=rz   ;endif
  IF(trim(adjustl(Mm))=="C")then   ; natom_C  = natom_C   + 1 ; rxC(natom_C)  =rx ;ryC(natom_C)=ry   ;rzC(natom_C)=rz   ;endif
  IF(trim(adjustl(Mm))=="Cu")then  ; natom_Cu = natom_Cu  + 1 ; rxCu(natom_Cu)=rx ;ryCu(natom_Cu)=ry ;rzCu(natom_Cu)=rz ;endif
  IF(trim(adjustl(Mm))=="Ag")then  ; natom_Ag = natom_Ag  + 1 ; rxAg(natom_Ag)=rx ;ryAg(natom_Ag)=ry ;rzAg(natom_Ag)=rz ;endif
  IF(trim(adjustl(Mm))=="Ni")then  ; natom_Ni = natom_Ni  + 1 ; rxNi(natom_Ni)=rx ;ryNi(natom_Ni)=ry ;rzNi(natom_Ni)=rz ;endif
  IF(trim(adjustl(Mm))=="Co")then  ; natom_Co = natom_Co  + 1 ; rxCo(natom_Co)=rx ;ryCo(natom_Co)=ry ;rzCo(natom_Co)=rz ;endif

 ENDIF
 IF(ierror > 0)then; nerror=nerror + 1 ; ENDIF      !Cuenta posibles errores de lectura.
 IF(ierror < 0)EXIT                                 !Encuentra el final del archivo.

ENDDO 

write(*,*)"nframes",nframe
write(*,*)"nconfig",nconfig
write(*,*)"nerror",nerror
write(*,*)"ntot",ntot        !Numero total de atomos por config


natom=natom/nconfig ; write(*,*)"natom",natom,"ntot",ntot

natom_Au=natom_Au/nconfig  ; write(*,*)"natom_Au",natom_Au
natom_S=natom_S/nconfig    ; write(*,*)"natom_S",natom_S
natom_N=natom_N/nconfig    ; write(*,*)"natom_N",natom_N
natom_C=natom_C/nconfig    ; write(*,*)"natom_C",natom_C
natom_Cu=natom_Cu/nconfig  ; write(*,*)"natom_Cu",natom_Cu
natom_Ag=natom_Ag/nconfig  ; write(*,*)"natom_Ag",natom_Ag
natom_Ni=natom_Ni/nconfig  ; write(*,*)"natom_Ni",natom_Ni
natom_Co=natom_Co/nconfig  ; write(*,*)"natom_Co",natom_Co

if(nframe<nblanco)then;print*,"nframe<nblanco";stop;endif
!....... Allocateo distancias e histograma  
ALLOCATE(rrAu(natom_Au),rrS(natom_S),rrN(natom_N),rrC(natom_C),rrCu(natom_Cu))
ALLOCATE(rrAg(natom_Ag),rrNi(natom_Ni),rrCo(natom_Co))
ALLOCATE(Hist_Au(nbin),Hist_S(nbin),Hist_N(nbin),Hist_C(nbin),Hist_Cu(nbin))
ALLOCATE(Hist_Ag(nbin),Hist_Ni(nbin),Hist_Co(nbin))
REWIND(123)
!..........................................................
DO j=1,nconfig
 !call cm(j,natom_Cu,nconfig,ntot,rxCu,ryCu,rzCu,rcm)
 
 call cm(j,natom,nconfig,ntot,rxtot,rytot,rztot,rcm)

 call dist_r(j,nconfig,natom_Au,ntot,nbin,rxAu,ryAu,rzAu,rrAu,rcm,delr_Au,rhist_Au,Hist_Au)
 call dist_r(j,nconfig,natom_Cu,ntot,nbin, rxCu, ryCu, rzCu, rrCu,rcm,delr_Cu,rhist_Cu, Hist_Cu)  

 call dist_r(j,nconfig,natom_Ag,ntot,nbin,rxAg,ryAg,rzAg,rrAg,rcm,delr_Ag,rhist_Ag,Hist_Ag)
 call dist_r(j,nconfig,natom_Ni,ntot,nbin,rxNi,ryNi,rzNi,rrNi,rcm,delr_Ni,rhist_Ni,Hist_Ni)
 call dist_r(j,nconfig,natom_Co,ntot,nbin,rxCo,ryCo,rzCo,rrCo,rcm,delr_Co,rhist_Co,Hist_Co)

 
! call dist_r(j,nconfig,natom_S,ntot,nbin, rxS, ryS, rzS, rrS,rcm,delr_S,rhist_S, Hist_S)
! call dist_r(j,nconfig,natom_N,ntot,nbin, rxN, ryN, rzN, rrN,rcm,delr_N,rhist_N, Hist_N)
! call dist_r(j,nconfig,natom_C,ntot,nbin, rxC, ryC, rzC, rrC,rcm,delr_C,rhist_C, Hist_C)  
ENDDO 
!............ Escritura del Histograma .... RDF
 
OPEN(201,file=name_files(unit_file)//'/RDF_Au.'//name_files(unit_file)//'.dat')
OPEN(202,file=name_files(unit_file)//'/RDF_S.'//name_files(unit_file)//'.dat')
OPEN(203,file=name_files(unit_file)//'/RDF_N.'//name_files(unit_file)//'.dat')
OPEN(204,file=name_files(unit_file)//'/RDF_C.'//name_files(unit_file)//'.dat')
OPEN(205,file=name_files(unit_file)//'/RDF_Cu.'//name_files(unit_file)//'.dat')
OPEN(206,file=name_files(unit_file)//'/RDF_Ag.'//name_files(unit_file)//'.dat')
OPEN(207,file=name_files(unit_file)//'/RDF_Ni.'//name_files(unit_file)//'.dat')
OPEN(208,file=name_files(unit_file)//'/RDF_Co.'//name_files(unit_file)//'.dat')

 DO i=1,nbin
 WRITE(201,*)real(i)*delr_Au,Hist_Au(i)/real(nconfig)
 WRITE(202,*)real(i)*delr_S,Hist_S(i)/real(nconfig)
 WRITE(203,*)real(i)*delr_N,Hist_N(i)/real(nconfig)          
 WRITE(204,*)real(i)*delr_C,Hist_C(i)/real(nconfig)
 WRITE(205,*)real(i)*delr_Cu,Hist_Cu(i)/real(nconfig)
 WRITE(206,*)real(i)*delr_Ag,Hist_Ag(i)/real(nconfig)
 WRITE(207,*)real(i)*delr_Ni,Hist_Ni(i)/real(nconfig)
 WRITE(208,*)real(i)*delr_Co,Hist_Co(i)/real(nconfig)
ENDDO
!............................................................
DEALLOCATE(rrAu,rrS,rrN,rrC,rrCu)
DEALLOCATE(rrAg,rrNi,rrCo)
DEALLOCATE(Hist_Au,Hist_S,Hist_N,Hist_C,Hist_Cu)
DEALLOCATE(Hist_Ag,Hist_Ni,Hist_Co)
CLOSE(201);CLOSE(202);CLOSE(203);CLOSE(204);CLOSE(205);CLOSE(206);CLOSE(207);CLOSE(208) 
END SUBROUTINE doRDF
!............................................................................................

SUBROUTINE cm(k,natom,nconfig,ntot,rx,ry,rz,rcm) 
 implicit none 
 integer,intent(in) :: k,natom,nconfig,ntot
 real(dp),dimension(nconfig*ntot),intent(in) :: rx,ry,rz    
 real(dp),dimension(3),intent(inout) :: rcm 
 real(dp) :: xcm,ycm,zcm
 integer :: i,j


 xcm=0.0_dp ; ycm=0.d0 ; zcm=0.d0   
 do i=(k-1)*natom+1,k*natom
  xcm=xcm+rx(i)
  ycm=ycm+ry(i)
  zcm=zcm+rz(i)
 enddo
 rcm(1)=xcm/real(natom);rcm(2)=ycm/real(natom);rcm(3)=zcm/real(natom)   
 END SUBROUTINE
!............................................................................................                 
 
SUBROUTINE dist_r(k,nconfig,natom,ntot,nbin,rx,ry,rz,rr,rcm,delr,rhist,Hist)
 implicit none
 integer,intent(in) :: k,nconfig,natom,ntot,nbin
 real(dp),dimension(nconfig*ntot),intent(in) :: rx,ry,rz
 real(dp),dimension(natom),intent(out)  :: rr          
 real(dp),dimension(3),intent(inout) :: rcm            
 real(dp),intent(in)   :: rhist
 real(dp) :: rrx,rry,rrz
 integer,intent(inout) :: Hist(nbin)
 real(dp),intent(in)  :: delr
 integer :: inb
 integer :: i,j
 
 if(natom==0) return


 write(*,*) rcm


 j=0                                  
 do i=(k-1)*natom+1,k*natom           
 j=j+1                                
 rrx=rx(i)-rcm(1)
 rry=ry(i)-rcm(2)                     
 rrz=rz(i)-rcm(3)                     
 rr(j)=dsqrt(rrx*rrx+rry*rry+rrz*rrz)
 enddo                                
  
 do i=1,natom
 inb=int(rr(i)/delr)+1
 Hist(inb)=Hist(inb)+1 
 enddo
 END SUBROUTINE
!............................................................................................ 

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

 END MODULE mod_RDF
