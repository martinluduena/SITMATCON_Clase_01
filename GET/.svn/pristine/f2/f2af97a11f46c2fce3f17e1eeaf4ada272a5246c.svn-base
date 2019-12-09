MODULE mod_count_line
public
integer :: nline,nbas
CONTAINS
SUBROUTINE count_lines()
IMPLICIT NONE
integer :: ierror
nline=0;nbas=0
DO
READ(123, *, IOSTAT = ierror) 
IF(ierror  > 0) nbas=nbas + 1     !Error en la lectura
IF(ierror == 0) nline=nline+1     !No hay error en la lectura
IF(ierror  < 0) EXIT              !Encontro el final del archivo
ENDDO;REWIND(123)

write(*,*)"nline = " ,nline
write(*,*)"nerrors = ",nbas
RETURN
END SUBROUTINE count_lines
END MODULE mod_count_line
