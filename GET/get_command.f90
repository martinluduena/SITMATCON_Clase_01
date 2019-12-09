PROGRAM GET
USE mod_get_comand    
!...................................................
IMPLICIT NONE
INTEGER :: jj,ii

!...................................................
CALL read_comand()
!...................................................
DO jj=1,5  !Se leen hasta 5 archivos.
 IF(LEN_TRIM(files(jj))==60)EXIT
  CALL open_file(jj) ; CALL open_directory(jj)
  DO ii=1,10
   IF(LEN_TRIM(makes(jj,ii))==25)EXIT
   CALL make_comand(jj,ii)
   CALL rewind_file()
  ENDDO
  CALL close_file()
ENDDO
END PROGRAM
!************       Fin del PROGRAMA        ****************** 
