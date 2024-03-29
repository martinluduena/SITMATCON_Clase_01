!OJO QUE HAY UN "ERROR", NO SE CUAL ES ; EL TEMA ESTA CUANDO
!SALGO DE LA SUBROUTINE ,... HAY ALGUN PROBLEMA DE MEMORIA
!O ALGO POR EL ESTILO. DE TODOS MODOS ASI COMO ESTA, ESTA FUNCIONANDO, 
!PARA ESTO PUSE UN PARCHE AL FINAL DE LA SUBROUTINE, "STOP", PARA QUE CORTE 
!AL FINAL Y NO SALGA DE SUBROUTINA, LOS RESULTADOS ESTAN TESTEADOS CON OTRO 
!CODIGO QUE NO DABA ESTE TIPO DE ERROR , ASI QUE NO HAY PROBLEMA EN 
!CUANTO A LA TRANSFORMADA DE FOURIER SE REFIERE.....!!

!Este MODULO realiza una transformada de fourier (FFTW3) de una dimension
!de una funcion REAL(fx) a una COMPLEJA(Fx_t).
!-----------------------------------------------------------------------------------
!Este codigo Fortran 90 usa la libreria FFTW-3.2.2 previamente instaladas en usr/local/lib:
!esta libreria fue bajada de la siguiente direccion :
!PLAN_R2C.0.2      www.fftw.org/download.html
!e instalada ejecutando los sig. comandos en el directorio donde fue guardado el paquete.
! ./configure
! make
! make install
!------------------------------------------------------------------------------- 
!--------------------------------------------------------------------------------
!Esta parte del programa realiza una transformada de fourier (FFTW3) de una dimension
!de una funcion REAL(fx) a una COMPLEJA(Fx_t).
!
!To transform a one-dimensional real array in Fortran, you might do:

             !double precision in
             !dimension in(N)
             !double complex out
             !dimension out(N/2 + 1)
             !integer*8 plan
     
             !call dfftw_plan_dft_r2c_1d(plan,N,in,out,FFTW_ESTIMATE)
             !call dfftw_execute_dft_r2c(plan, in, out)
             !call dfftw_destroy_plan(plan)
!---------------------------------------------------------------------------------	     
                                       !DATOS

!     1 **------------+------***----+------------+***----------+--------***-++
!      +**           +     ** **   +            ** **        cos(4t) ****** +
!  0.8 ++**               **   *                *   **               **   **+
!      |  *               *    **              **    *               *     *|
!  0.6 ++ **             **     *              *     **              *     **
!      |   *             *      **            **      *             **      *
!  0.4 ++  *            **       *            *       **            *      +*
!      |   **           *        *            *        *           **       |
!  0.2 ++   *           *        **          **        *           *       ++
!    0 ++   *          **         *          *         **          *       ++
!      |    **         *          *          *          *         **        |
! -0.2 ++    *         *          **        **          *         *        ++
!      |     **       **           *        *            *        *         |
! -0.4 ++     *       *            **      **            *       *         ++
!      |      *      **             *      *             **      *          |
! -0.6 ++     **     *              *      *              *     **         ++
!      |       *    **              **    **              **    *           |
! -0.8 ++      **   *                **  **                *   **          ++
!      +        ** **+             +  *  *      +          ** **            +
!   -1 ++--------***-+-------------+--****------+-----------***------------++
!      0             1             2            3             4             5

                                !Ejemplo
!Si,..      fx=a*sen(2pi*f1*t) + b*cos(2pi*f2*t)
!Entonces.. Fx_t(f,y)= (f1,a) ; (f2,b) 

!Puede ser que la altura no de justo la amplitud , puesto que puede
!ser que no caigas justo en la frecuencia , y eso hace que caiga un poquito por
!debajo de lo esperado para esa frecuancia.


                                    !RESULTADO
!  .7 ++---------------+-----------------+----------------+---------------++
!      +      *         +                 +       'FT[cos(4t)]' ******      +
!      |      *                                                             |
!  0.6 ++     *                                                            ++
!      |      *                                                             |
!  0.5 ++     *                                                            ++
!      |      *                                                             |
!      |      *                                                             |
!  0.4 ++     *                                                            ++
!      |      *                                                             |
!      |      *                                                             |
!  0.3 ++     *                                                            ++
!      |      *                                                             |
!      |      *                                                             |
!  0.2 ++     **                                                           ++
!      |     * *                                                            |
!  0.1 ++    * *                                                           ++
!      |     * *                                                            |
!      +   **   *****   +                 +                +                +
!    0 *******|***************************************************************
!      0      2          5                 10               15               20
!                                  

!__________________________________________________________________ Martin Ludue�a 

MODULE mod_FFT
USE mod_count_line
integer,parameter :: dp=kind(0.d0)
CONTAINS

SUBROUTINE doFFT(unit_file,name_files)
IMPLICIT NONE 
include "fftw3.f"
 integer,intent(in) :: unit_file
 CHARACTER(len=*),dimension(5),intent(in) :: name_files 
!------------------------------------------------- 
real(dp),dimension(:),allocatable :: fx,fy
complex(dp),dimension(:),allocatable :: Fx_t,fxx
!--------------------------------------------------
integer :: N_fftw   !grillado de la f a transformar
integer :: N_fftwi 
!................... Planes ....................
integer :: plan_r2c
integer :: plan_c2c
!-----------------------------------------------
real(dp) :: pi  ! N�mero Pi 
integer  :: i , j , n
real(dp) :: ft, dt , t1 , t2   
real(dp) :: t,spot
real(dp) :: v_fftw,w_fftw,fact_fftw
!-----------------------------------------------
real(dp) :: xreal,yreal
!-----------------------------------------------
CALL count_lines()    !Cuenta la cantidad de lineas 

! write(*,*)"n_lineas................",nline 
 !write(*,*)"WARNING : errores de lectura = ",nbas

 N_fftw=nline      
 N_fftwi=N_fftw/2 

 ALLOCATE(fx(1:N_fftw),fy(1:N_fftw))
 ALLOCATE(Fx_t(1:N_fftwi),fxx(1:N_fftwi))
 
 !Lectura del Archivo de Datos
  read(123,*)t1,fx(1)
   do i=2,N_fftw-1
    read(123,*)t,fx(i)
   enddo          
  read(123,*)t2,fx(N_fftw)  
  dt=(t2-t1)/dble(N_fftw-1)  !Paso de tiempo

 pi=dacos(-1._dp)
 v_fftw=1._dp/(dt*dble(N_fftw))    !frecuencia  (v)
 w_fftw=pi/(dt*dble(N_fftw))       !frecuencia angular (w)
 fact_fftw= 2._dp/dble(N_fftw)     !Con este valor normalizo para que me de los coeficientes de la fx

call dfftw_plan_dft_r2c_1d(plan_r2c,N_fftw,fx,Fx_t,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan_r2c,fx,Fx_t)                                     

!------------------------------------------------------------------------------------------------------  
!OUTPUT : Transformada de Fourier
open(200,file=name_files(unit_file)//'/FFT.'//name_files(unit_file)//'.dat')  !'ABS(Trasformada de la funcion)')

 do i=1,N_fftwi
  spot=abs(Fx_t(i))*fact_fftw 
  write(200,*) v_fftw*dble(i-1),spot !real(Fx_t(i))*fact_fftw,aimag(Fx_t(i))*fact_fftw    
 enddo

 call dfftw_destroy_plan(plan_r2c)
!---------------------------------------------------------------------------------------------------
!TRANSFORMADA INVERSA
!Esta parte del programa realiza una transformada de fourier inversa (FFTW3) de una dimension
!de una funcion COMPLEJA(Fx_t) a una COMPLEJA(fxx)

!OUTPUT : Transformada Inversa
 open(300,file=name_files(unit_file)//'/anti_FFT.'//name_files(unit_file)//'.dat') !inversa de la transformada
 
 call dfftw_plan_dft_1d(plan_c2c,N_fftwi,Fx_t,fxx,FFTW_BACKWARD,FFTW_ESTIMATE)
 call dfftw_execute_dft(plan_c2c,Fx_t,fxx) 
  t=0._dp
  do i=1,N_fftwi
   write(300,*)t,real(fxx(i))*fact_fftw
   t=t + 2._dp*dt
  enddo
   
 call dfftw_destroy_plan(plan_c2c)
 
 
 close(200)
 close(300)
stop 
END SUBROUTINE doFFT

END MODULE mod_FFT
