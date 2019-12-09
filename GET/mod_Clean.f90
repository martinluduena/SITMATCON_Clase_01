!Este Módulo 
MODULE mod_Clean
USE mod_count_line
integer,parameter,private :: dp=kind(0.0d0)
contains
SUBROUTINE doClean(unit_file,name_files)
implicit none
integer,intent(in) :: unit_file
CHARACTER(len=*),dimension(5),intent(in) :: name_files 

 
real(dp),dimension(:),allocatable :: rrAu,rrS,rrN,rrC
real(dp) :: rx,ry,rz !Posiciones de los atomos (auxiliar)
real(dp),dimension(:),allocatable :: rxAu,ryAu,rzAu
real(dp),dimension(:),allocatable :: rxS,ryS,rzS
real(dp),dimension(:),allocatable :: rxN,ryN,rzN
real(dp),dimension(:),allocatable :: rxC,ryC,rzC

integer :: ntot=0             !numero de atomos totales
integer :: natom=0
integer :: natom_Au=0         !numero de atomos de Au
integer :: natom_S=0          !numero de atomos de S
integer :: natom_N=0          !numero de atomos de N
integer :: natom_C=0          !numero de atomos de C
integer :: largo_cadena=0       !Cantidad de C en la Cadena
integer :: mol_ads_N=0,mol_ads_S=0
integer :: nerror=0
integer :: ierror
integer :: i,j,k,ii,jj !............Contadores 
integer :: iiAu
integer :: iiN
integer :: iiS
integer :: iiC


character(len=3) :: Mm  !Identidad de los elementos
character(len=4) :: lline
integer :: nframe=0
!.................. RDF .......................................

real(dp) :: dist_N_Au=3.5d0         !Distancia de Maxima Adsorcion N-Au
real(dp) :: dist_S_Au=3.5d0         !Distancia de Maxima Adsorcion S-Au

real(dp) :: dist2_N_Au
real(dp) :: dist2_S_Au


real(dp) :: drx,dry,drz,dr2

dist2_N_Au=dist_N_Au*dist_N_Au
dist2_S_Au=dist_S_Au*dist_S_Au
!................................................................
 CALL count_lines()
!.......................
READ(123,*) ntot !Numero de atomos totales por config
REWIND(123) 
nframe=nline/(ntot+2)
!...... Cuento Numero de atomos .............    
 DO 
 READ(123,'(a2,x,3(f10.6,x))',IOSTAT = ierror) Mm ,rx,ry,rz

 IF(ierror == 0)then !No hay error en la lectura
 
  IF(trim(adjustl(Mm))=="Au")then ; natom_Au = natom_Au  + 1 ;endif
  IF(trim(adjustl(Mm))=="S")then  ; natom_S  = natom_S   + 1 ;endif
  IF(trim(adjustl(Mm))=="N")then  ; natom_N  = natom_N   + 1 ;endif
  IF(trim(adjustl(Mm))=="C")then  ; natom_C  = natom_C   + 1 ;endif
 
 ENDIF
 
 IF(ierror > 0)then; nerror=nerror + 1 ; ENDIF      !Cuenta posibles errores de lectura.
 IF(ierror < 0)EXIT                                 !Encuentra el final del archivo.

 ENDDO 

write(*,*)"nframes",nframe
write(*,*)"nerror",nerror
write(*,*)"ntot",ntot        !Numero total de atomos por config

natom_Au=natom_Au/nframe  ; write(*,*)"natom_Au",natom_Au
natom_S=natom_S/nframe    ; write(*,*)"natom_S",natom_S
natom_N=natom_N/nframe    ; write(*,*)"natom_N",natom_N
natom_C=natom_C/nframe    ; write(*,*)"natom_C",natom_C

 if(natom_S==0)then 
 largo_cadena= natom_C/natom_N ; write(*,*)"largo_cadena",largo_cadena
  elseif(natom_N==0)then
 largo_cadena= natom_C/natom_S ; write(*,*)"largo_cadena",largo_cadena
  elseif(natom_S/=0.and.natom_N/=0)then
 write(*,*)"cuantos atomos tiene la cadena de carbonos?"
 read(*,*) largo_cadena
 endif

 
ALLOCATE(rxAu(nframe*natom_Au),ryAu(nframe*natom_Au),rzAu(nframe*natom_Au)) 
ALLOCATE( rxS(nframe*natom_S), ryS(nframe*natom_S), rzS(nframe*natom_S))    
ALLOCATE( rxN(nframe*natom_N), ryN(nframe*natom_N), rzN(nframe*natom_N))    
ALLOCATE( rxC(nframe*natom_C), ryC(nframe*natom_C), rzC(nframe*natom_C))    

natom_Au=0 ; natom_S=0 ; natom_N=0 ; natom_C=0 

REWIND(123)

DO 
 READ(123,'(a2,x,3(f10.6,x))',IOSTAT = ierror) Mm ,rx,ry,rz

 IF(ierror == 0)then !No hay error en la lectura
 
  IF(trim(adjustl(Mm))=="Au")then ; natom_Au = natom_Au  + 1 ; rxAu(natom_Au)=rx ;ryAu(natom_Au)=ry ;rzAu(natom_Au)=rz ;endif
  IF(trim(adjustl(Mm))=="S")then  ; natom_S  = natom_S   + 1 ; rxS(natom_S)  =rx ;ryS(natom_S)=ry   ;rzS(natom_S)=rz   ;endif
  IF(trim(adjustl(Mm))=="N")then  ; natom_N  = natom_N   + 1 ; rxN(natom_N)  =rx ;ryN(natom_N)=ry   ;rzN(natom_N)=rz   ;endif
  IF(trim(adjustl(Mm))=="C")then  ; natom_C  = natom_C   + 1 ; rxC(natom_C)  =rx ;ryC(natom_C)=ry   ;rzC(natom_C)=rz   ;endif
 
 ENDIF
 
 IF(ierror > 0)then; nerror=nerror + 1 ; ENDIF      !Cuenta posibles errores de lectura.
 IF(ierror < 0)EXIT                                 !Encuentra el final del archivo.

ENDDO 

REWIND(123)

OPEN(201,file=name_files(unit_file)//'/Clean.'//name_files(unit_file)//'.xyz')
OPEN(202,file=name_files(unit_file)//'/Adorvidos_NH2.dat')
OPEN(204,file=name_files(unit_file)//'/Adorvidos_S.dat')
OPEN(203,file=name_files(unit_file)//'/configuraciones.'//name_files(unit_file))


 natom_Au=natom_Au/nframe
 natom_S=natom_S/nframe  
 natom_N=natom_N/nframe  
 natom_C=natom_C/nframe  



DO i=1,nframe

 write(201,*)ntot ; write(203,*)
 write(201,*)     ; write(203,*)
 
  DO iiAu=1,natom_Au !Escrivo la configuracion de la NP
  j=(i-1)*natom_Au+iiAu
  write(201,'(a2,x,3(f10.6,x))') "Au",rxAu(j),ryAu(j),rzAu(j)
  write(203,'(a2,x,3(f10.6,x))') "Au",rxAu(j),ryAu(j),rzAu(j)
  ENDDO

               mol_ads_N=0
   DO iiN=1,natom_N ;if(size(rxN)==0)EXIT
  
   ii=(i-1)*natom_N+iiN         !Index de N
     
       DO iiAu=1,natom_Au
       jj=(i-1)*natom_Au+iiAu       !INDEX DE Au
        drx=rxN(ii)-rxAu(jj) ; dry=ryN(ii)-ryAu(jj) ; drz=rzN(ii)-rzAu(jj)

        dr2= drx*drx + dry*dry + drz*drz

        if(dr2<dist2_N_Au)then !Si encuentro 1, salgo...! y escribo la
       !coordenada de la molecula.

        mol_ads_N=mol_ads_N+1  
       
        write(201,'(a2,x,3(f10.6,x))') "N",rxN(ii),ryN(ii),rzN(ii)
        write(203,'(a2,x,3(f10.6,x))') "N",rxN(ii),ryN(ii),rzN(ii)
       
        DO iiC=1,largo_cadena !Largo de Cadena ...!
        k=(ii-1)*largo_cadena+iiC
         write(201,'(a2,x,3(f10.6,x))')"C",rxC(k),ryC(k),rzC(k)
         write(203,'(a2,x,3(f10.6,x))')"C",rxC(k),ryC(k),rzC(k) !en 203 escribo
                                                               !solo lo que de adsrobio 
        ENDDO 
       
        EXIT 

        elseif(iiAu==natom_Au)then

       !En VMD, no se puede modificar el numero de atomos en la corrida,..! Hay
       !que hacer un truco. poner el maximo numero de atomos que se va a tener,
       !pero si hay moleculas, tambien hay que tener en cuenta la conectividad
       !de los atomos desde el comienzo. por lo tanto, los atomos de moleculas,
       !que aun no deben aparacer, deben ser atomos Dummy, pero respetando la
       !conectividad que luego van a tener cuando aparescan en escena. Por esto
       !es importante escribirlos respetando no solo la naturaleza de los atomos
       !sino tambien la conectividad que en el futuro vban a tener ...! 


        write(201,'(a2,x,3(f10.6,x))') "N",100.d0,100.d0,100.d0 ! Deben ir en un
       !lugar donde no sean visto...! puede ser el centro de masa de la NP. 
       
        DO iiC=1,largo_cadena !Largo de Cadena ...!
         write(201,'(a2,x,3(f10.6,x))')"C",100.d0,100.d0,100.d0
        ENDDO 

      endif

      ENDDO   
   ENDDO
   write(202,*)i, mol_ads_N

               mol_ads_S=0
   DO iiS=1,natom_S ; if(size(rxS)==0)EXIT
    ii=(i-1)*natom_S+iiS         !Index de S
      DO iiAu=1,natom_Au
       jj=(i-1)*natom_Au+iiAu       !INDEX DE Au
        drx=rxS(ii)-rxAu(jj) ; dry=ryS(ii)-ryAu(jj) ; drz=rzS(ii)-rzAu(jj)

        dr2= drx*drx + dry*dry + drz*drz

        if(dr2<dist2_S_Au)then !Si encuentro 1, salgo...! y escribo la
        !coordenada de la molecula.

        mol_ads_S=mol_ads_S+1  
        
        write(201,'(a2,x,3(f10.6,x))') "S",rxS(ii),ryS(ii),rzS(ii)
        write(203,'(a2,x,3(f10.6,x))') "S",rxS(ii),ryS(ii),rzS(ii)
        
        DO iiC=1,largo_cadena !Largo de Cadena ...!
        k=(ii-1)*largo_cadena+iiC
         write(201,'(a2,x,3(f10.6,x))')"C",rxC(k),ryC(k),rzC(k)
         write(203,'(a2,x,3(f10.6,x))')"C",rxC(k),ryC(k),rzC(k) !en 203 escribo
                                                                !solo lo que de adsrobio 
        ENDDO 
        
        EXIT 

        elseif(iiAu==natom_Au)then

       !En VMD, no se puede modificar el numero de atomos en la corrida,..! Hay
       !que hacer un truco. poner el maximo numero de atomos que se va a tener,
       !pero si hay moleculas, tambien hay que tener en cuenta la conectividad
       !de los atomos desde el comienzo. por lo tanto, los atomos de moleculas,
       !que aun no deben aparacer, deben ser atomos Dummy, pero respetando la
       !conectividad que luego van a tener cuando aparescan en escena. Por esto
       !es importante escribirlos respetando no solo la naturaleza de los atomos
       !sino tambien la conectividad que en el futuro vban a tener ...! 


        write(201,'(a2,x,3(f10.6,x))') "S",100.d0,100.d0,100.d0 ! Deben ir en un
        !lugar donde no sean visto...! puede ser el centro de masa de la NP. 
        
        DO iiC=1,largo_cadena !Largo de Cadena ...!
         write(201,'(a2,x,3(f10.6,x))')"C",100.d0,100.d0,100.d0
        ENDDO 

      endif

      ENDDO   
   ENDDO
write(204,*)i, mol_ads_S


ENDDO




CLOSE(201);CLOSE(202);CLOSE(203)
 END SUBROUTINE doClean
END MODULE mod_Clean
