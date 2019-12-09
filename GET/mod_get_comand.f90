MODULE mod_get_comand
USE mod_RDF
USE mod_PDF
USE mod_Clean
USE mod_SMTB
USE mod_BORDER
USE mod_FFT
USE mod_DoS_2PT
USE mod_DoS_low
USE mod_FEFFinp

CHARACTER(len=60),dimension(5)    :: files,name_files
CHARACTER(len=25),dimension(5,10) :: makes
CONTAINS

SUBROUTINE read_comand()                
IMPLICIT NONE                                                   
!INTEGER :: i
INTEGER :: jj=0 !Unit_file
INTEGER :: ii=0 !numero de comandos

INTEGER :: arglen
INTEGER :: inunit,logunit
CHARACTER(len=32) :: arg                                        
CHARACTER(len=25) :: begin_arg,make_arg                 
CHARACTER(len= 1) :: point_arg
CHARACTER(len=4)  :: ext_get
CHARACTER(len=32) :: get_file


LOGICAL :: eof
CHARACTER(LEN=16) :: w
DOUBLE PRECISION :: T1, T2, DT, R1, R2, DR, lfact=1d0, Tzero=0d0, &
    value(10)
INTEGER :: i=1, int(10)


i= 1 ;DO                                                       
 CALL get_command_argument(i,arg,arglen) ;IF(arglen==0)EXIT 
! ..................... doCOMANDs por Linea de Comandos ......... ! 

 !Averiguo si es un archivo                                                              
  point_arg=trim(arg(arglen-3:arglen-3))                        
   if(point_arg==".")then  !Es un Archivo                                      
    jj=jj+1;ii=1                                                  
    write(files(jj),*) trim(arg)  !name.ext                                 
    write(name_files(jj),*) trim(arg(1:arglen-4))  !name                                                   
   
    write(*,*) "name=" , name_files(jj)
  endif                                                         
  !----------------------------------                            
 !Averiguo comienzo de arg.                                     
  begin_arg=trim(arg(1:2))                                      
  if(begin_arg=="do")then                                       
  make_arg=trim(arg(3:arglen))                                  
  write(makes(jj,ii),*)trim(make_arg) !Escribe adentro de makes, el nombre del comando                           
  
  write(*,*) "comand=",makes(jj,ii)
  
  ii=ii+1                                                       
  endif


 !----------------------------------                            
 i=i+1;ENDDO                                                 
END SUBROUTINE read_comand                                  
!............................................................
SUBROUTINE open_file(jj)                    
 IMPLICIT NONE                                                                
 INTEGER :: ierror ,jj
 OPEN(UNIT=123,file=trim(adjustl(files(jj))),IOSTAT=ierror) 
 END SUBROUTINE open_file                          
!............................................................
SUBROUTINE open_directory(jj)
!ESTA SUBROUTINE ABRE UNA CARPETA CON EL NOMBRE DEL ARCHIVO CON EL QUE ESTAS
!TRABAJANDO.
!VER SI SE PUEDE HACER QUE HABRA UNA SUB-CARPETA ADEMAS, CON EL NOMBRE DEL COMANDO
!QUE VOY A EJECUTAR. TENDRIA QUE PASARLE EL MAKES, Y QUE HABRA ADENTRO OTRA
!CARPETA CON ESE NOMBRE.

!Con el 2>/dev/null , (2=error estandar), mando el error al dev/null , 
!asi evito el mensaje de alerta de error, en caso de querer crear 
!una carpeta que ya existia..! Hay que tener ojo de todos modos, 
!porque se sobreescriben los archivos dentro de la carpeta...!
IMPLICIT NONE                                                                
INTEGER :: jj
character(66) :: mkdir_directory
mkdir_directory='mkdir '//trim(adjustl(name_files(jj)))//' 2>/dev/null'
CALL SYSTEM (mkdir_directory)
END SUBROUTINE open_directory 
!.............................................................
SUBROUTINE make_comand(jj,ii)                        
IMPLICIT NONE                                              
INTEGER,intent(in) :: jj,ii
CHARACTER(len=25) :: aux                                   
aux=makes(jj,ii) 

SELECT CASE (TRIM(ADJUSTL(aux)))
 
 CASE("RDF")        ; write(*,*) "doing",aux ;call doRDF(jj,trim(adjustl(name_files(jj))))
 CASE("PDF")        ; write(*,*) "doing",aux ;call doPDF(jj,trim(adjustl(name_files(jj))))
 CASE("Clean")      ; write(*,*) "doing",aux ;call doClean(jj,trim(adjustl(name_files(jj))))
 CASE("tb")         ; write(*,*) "doing",aux ;call doSMTB(jj,trim(adjustl(name_files(jj))))
 CASE("BORDER")     ; write(*,*) "doing",aux ;call doBORDER()  
 CASE("FFT")        ; write(*,*) "doing",aux ;call doFFT(jj,trim(adjustl(name_files(jj))))
 CASE("DoS")        ; write(*,*) "doing",aux ;call doDoS(jj,trim(adjustl(name_files(jj))))
 CASE("DoSlow")     ; write(*,*) "doing",aux ;call doDoSlow(jj,trim(adjustl(name_files(jj))))
 CASE("FEFFinp")    ; write(*,*) "doing",aux ;call doFEFFinp()    
 CASE("EXAFS")      ; write(*,*) "doing",aux
 
 CASE DEFAULT        ;call write_lista()
END SELECT
RETURN
ENDSUBROUTINE make_comand
!............................................................                                                    
 SUBROUTINE rewind_file()
 rewind(123)
 END SUBROUTINE rewind_file                                    
!............................................................                                                   
 SUBROUTINE close_file()                           
 close(123)                                        
 END SUBROUTINE close_file  
!............................................................
SUBROUTINE write_coments()                                  

 IMPLICIT NONE                                                   
 write(*,*)"**************************************************"  
 write(*,*)"Se ingresan comandos por linea de comandos"          
 write(*,*)"Ej: ./get velocidades.xyz doDoSlow"               
 write(*,*)"Los archivos con datos deben llamarse : nombre.xxx"  
 write(*,*)"Ej: Coor.xyz ; Velocidades.xys; energias_pot.txt"    
 write(*,*)"Los comandos de acción empiezan con do---"           
 write(*,*)"Ej: doDoSlow; realiza el calculos de Densidad de Estados Vibracioanles (vDoS)"                
 write(*,*)"**************************************************"  
 END SUBROUTINE
!..............................................................- 
 SUBROUTINE write_lista()
  IMPLICIT NONE
  write(*,*) "NO elegiste un make_comand correcto" 
  write(*,*) "Lista de make comands es:"           
  write(*,*)  " doRDF"
  write(*,*)  " doPDF"
  write(*,*)  " doClean"
  write(*,*)  " doBORDER"                          
  write(*,*)  " doFFT"                             
  write(*,*)  " doDoS"                             
  write(*,*)  " doDoSlow"                             
  write(*,*)  " doFEFFinp"                         
  write(*,*)  " doEXAFS"     
  call write_coments()
 END SUBROUTINE


END MODULE mod_get_comand 
