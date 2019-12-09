
!--------------------------------------------------------------------------------------------------
!Este codigo Fortran 90 usa la libreria FFTW-3.2.2 previamente instaladas en usr/local/lib:
!esta libreria fue bajada de la siguiente direccion :
!      www.fftw.org/download.html    
!e instalada ejecutando los sig. comandos en el directorio donde fue guardado el paquete.
! ./configure
! make
! make install 
!---------------------------------------------------------------------------------------------------	     
!____________________________________________________________________________________ Martin Ludueña 

!     4 ++--------+---------+---------+--------+---------+---------+--------++
!      +         *         +         +        +  './DoS_Co_cm-1.dat' ****** +
!  3.5 ++        *                                                         ++
!      |         * *                                                        |
!      |         ***                                                        |
!    3 ++        ***                                                       ++
!      |         ***                                                        |
!  2.5 ++        ****                                                      ++
!      |       * *****                                                      |
!    2 ++      ******* *                                                   ++
!      |      ******** *   *                                                |
!      |      ******** *  **                                                |
!  1.5 ++    **** ****** ***                                               ++
!      |     **** ****** ****                                               |
!    1 ++    ****   ******  *                                              ++
!      |     ****   ******  ****** * *   *                                  |
!      |     ****     * **  *********** *** **       *                      |
!  0.5 ++    ***             **************************                    ++
!      +    ***  +         +        *+******************+        *+ ***     +
!    0 ******----+---------+---------+--------+--------**********************
!      0         50       100       150      200       250       300       350
MODULE mod_DoS_2PT
USE mod_count_line
USE biseccion
USE Simpson
USE weighting_functions

integer,parameter :: dp=kind(0.d0)

!--------------------------------------------------
integer :: N_fftw   !grillado de la f a transformar
integer :: N_fftwi
integer :: Natom
!--------------------------------------------------
 real(dp),dimension(:),allocatable :: fx,fy,fz       !Funcion a transformar
 complex(dp),dimension(:),allocatable :: Fx_t        !transformada de la fx
 complex(dp),dimension(:),allocatable :: Fy_t        !transformada de la fy
 complex(dp),dimension(:),allocatable :: Fz_t        !transformada de la fz
!----------------------------------------------------------------------------------------------------
 complex(dp),dimension(:),allocatable :: Rxyz       !Funcion de Autocorrelacion Total
 complex(dp),dimension(:),allocatable :: Rxx        !Funcion de Autocorrelacion de Vx 
 complex(dp),dimension(:),allocatable :: Ryy        !Funcion de Autocorrelacion de Vy
 complex(dp),dimension(:),allocatable :: Rzz        !Funcion de Autocorrelacion de Vz
!----------------------------------------------------------------------------------------------------
 complex(dp),dimension(:,:),allocatable :: spot_x,spot_y,spot_z  !Espectro de Energia (/Fx_t/**2)
 complex(dp),dimension(:),allocatable :: DoS,dos_x,dos_y,dos_z   !Densidad de Estados
 complex(dp),dimension(:),allocatable :: dos_core,dos_shell   
 !---------------------------------------------------------------------------------------------------
 real(dp),dimension(:,:),allocatable :: vx,vy,vz     !Velocidades en x,y,z
 !----------------------------------------------------------------------------------------------------
 real(dp),dimension(:),allocatable :: mas_v          !Vector masa de atomos
 real(dp),dimension(:),allocatable :: OE             !Orden de Enlace de los Atomos
 !----------------------------------------------------------------------------------------------------
 !----------------  MODELO 2PT ------------------------------------
 real(dp),dimension(:),allocatable :: DoS_gas,Dos_S_Debye
 real(dp),dimension(:),allocatable :: g_SWGS      !DoS(Tipo Gas) * WGS
 real(dp),dimension(:),allocatable:: g_SWSS       !DoS(Tipo Solido) * Ws_S (Ws_S = funcion de peso Entropia solido)
 real(dp),dimension(:),allocatable :: g_SDoS      !DoS(Total)*Ws_S 
  
 real(dp),dimension(:),allocatable :: g_SWGA       !DoS(Tipo Gas) * WGA
 real(dp),dimension(:),allocatable :: g_SWSA       !DoS(Tipo Solido) * Ws_A (Ws_A = funcion de peso E. Helmholtz solido)
 real(dp),dimension(:),allocatable :: g_ADoS       !DoS(Total)*Ws_A 
!----------------------------------------------------------
 real(dp)  :: SWGS             !Entropia Fase tipo Gas
 real(dp)  :: SWSS             !Entropia Fase tipo Solido
 real(dp)  :: S2PT             !Entropia del Sistema (con 2PT model) 
 real(dp)  :: SDoS             !Entropia del Sistema (sin 2PT model)    
  
 real(dp)  :: SWGA             !E. Helmholtz Fase tipo Gas
 real(dp)  :: SWSA             !E. Helmholtz Fase tipo Solido
 real(dp)  :: A2PT             !E. Helmholtz del Sistema (con 2PT model) 
 real(dp)  :: ADoS             !E. Helmholtz del Sistema (sin 2PT model)    
!-----------------------------------------------------------------------------------------------------
CONTAINS

SUBROUTINE ALLOCATER()

ALLOCATE(fx(1:N_fftw),fy(1:N_fftw),fz(1:N_fftw))

ALLOCATE(Fx_t(1:N_fftwi),Fy_t(1:N_fftwi),Fz_t(1:N_fftwi))

ALLOCATE(Rxyz(1:N_fftwi),Rxx(1:N_fftwi),Ryy(1:N_fftwi),Rzz(1:N_fftwi))

ALLOCATE(spot_x(1:Natom,1:N_fftwi),spot_y(1:Natom,1:N_fftwi),spot_z(1:Natom,1:N_fftwi))

ALLOCATE(DoS(1:N_fftwi),dos_x(1:N_fftwi),dos_y(1:N_fftwi),dos_z(1:N_fftwi))

ALLOCATE(vx(1:Natom,1:N_fftw),vy(1:Natom,1:N_fftw),vz(1:Natom,1:N_fftw))

ALLOCATE(OE(1:Natom),mas_v(1:Natom))

ALLOCATE(DoS_gas(1:N_fftwi),DoS_S_Debye(1:N_fftwi))


ALLOCATE(g_SWGS(1:N_fftwi))
ALLOCATE(g_SWSS(1:N_fftwi))
ALLOCATE(g_SDoS(1:N_fftwi))     
ALLOCATE(g_SWGA(1:N_fftwi))          
ALLOCATE(g_SWSA(1:N_fftwi))
ALLOCATE(g_ADoS(1:N_fftwi))
               
RETURN
END SUBROUTINE

SUBROUTINE doDoS(unit_file,name_files)
IMPLICIT NONE 
include "fftw3.f"
 integer,intent(in) :: unit_file
 CHARACTER(len=*),dimension(5),intent(in) :: name_files  !asi entra con la dimension que venga de la variable Dummy
 integer :: plan_r2c
 integer :: plan_c2c
!-----------------------------------------------------------------------------------------------------
! MASAS ATOMICAS
!----------------- Metales de Transicion---------------------------------------------
real(dp),parameter :: masa_Pt= 3.239422119e-25_dp         ! (input) Masa del Platino [kg]
real(dp),parameter :: masa_Au= 3.271338426e-25_dp         ! (input) Masa del Oro [kg]                                                                                     
real(dp),parameter :: masa_Co= 9.786086422e-26_dp         ! (input) Masa del Oro [kg]
!..................GASES NOBLES........................................................
real(dp),parameter :: masa_Ar=6.63e-26_dp                 ! (input) Masa del Argon [kg]
!-----------------------------------------------------------------------------------------------------
real(dp) :: masa                                       !real auxiliar de mas_v
!-----------------------------------------------------------------------------------------------------
                                !Constantes Universales
      
real(dp),parameter :: kb= 1.38e-23_dp           !Constante de Boltzman  [J/k]   		       
real(dp),parameter :: v_luz=29979245800._dp     !velocidad de la luz [cm/s]
real(dp),parameter :: hp=6.62606896e-34_dp      !Constante de planck (cuanto elemental de acción) [Js]
real(dp),parameter :: pi=3.141592654_dp         !Numero pi
real(dp),parameter :: N_av=6.02214179e23_dp     !Numero de Avogadro
real(dp) :: bhc                                 !hp*v_luz/KT [cm]
real(dp) :: So
real(dp) :: CDN
real(dp) :: Ycs
real(dp) :: ZCS
real(dp) :: f
real(dp) :: WGS
real(dp) :: WGE
real(dp) :: WGA

!-----------------------------------------------------------------------------------------------------
!Propiedades Mecanicas
real(dp) :: TT                                           !Temperatura (input) 
real(dp) :: Temp=1.1_dp
real(dp) :: Etotal=-2295.998_dp                          !Energia total reducida   
real(dp) :: p=7.0125e28_dp                               !Densidad    (1/m3)
!----------------------------------------------------------------------------------------------------- 
 integer  :: i,j                               !Contadores i=>(Natom) ,j=>(N_fftw)
 real(dp) :: dt                                !paso de tiempo (tiempo_total /Nª_datos)  (input)                                               
 real(dp) :: vnu                               !numero de onda [1/cm]
 real(dp) :: t                                 !tiempo
 integer  :: np,num
 !....................................................................................................
 character(2) :: Id                             !identidad de los atomos
 character(2) :: Au,Co,Pt                       !Metales de Transicion
 character(2) :: Ar                             !Gases Nobles
 !....................................................................................................
 real(dp) :: vel2                              !acumulador de velocidad (Ecinetica) 
!------------------------------------------------------------------------------
                                !Factores para convertir unidades
  real(dp),parameter :: Afs_ms=1.0e5_dp 

  real(dp),parameter :: Aps_ms=1.0e2_dp !Unidades q usa GeMS
!------------------------------------------------------------------------------------------------------
                                 !Factores de Escaleo 
 
  real(dp) :: w_fftw             ! (escalea a frecuncia (Hz))
  real(dp) :: nu_fftw            ! (escalea a numero de onda (cm-1))
  real(dp) :: fact_fftw          ! Con este valor normalizo para que me de los coeficientes de la fx
  real(dp) :: f_fftw_cm          ! (Dimensionaliza DoS a cm)
!_________________________--------OUTPUT-------_______________________________
open (70,file=name_files(unit_file)//'/VAC.'//name_files(unit_file)//'.dat')             !Output; escribe :: t,Rxyz,Rxx,Ryy,Rfzz
open (100,file=name_files(unit_file)//'/DoS_'//name_files(unit_file)//'_Hz.dat')         !Output; escribe :: vnu,DoS,dos_x,dos_y,dos_z (Hz) 
open (101,file=name_files(unit_file)//'/DoS_'//name_files(unit_file)//'_cm-1.dat')                       !Output; escribe :: 
!--------OUTPUTS (Densidad de Estados proyectada sobre CORE/SHELL)---------
open (102,file=name_files(unit_file)//'/core_dos_'//name_files(unit_file)//'_Hz.dat')
open (103,file=name_files(unit_file)//'/core_dos_'//name_files(unit_file)//'_cm-1.dat')
open (104,file=name_files(unit_file)//'/shell_dos_'//name_files(unit_file)//'_Hz.dat')
open (105,file=name_files(unit_file)//'/shell_dos_'//name_files(unit_file)//'_cm-1.dat')
!________________________________________________________________________________________
!----- 2PT ----
open (106,file=name_files(unit_file)//'/DoS_2PT'//name_files(unit_file)//'_cm-1.dat')
open (107,file=name_files(unit_file)//'/Prop_Termo'//name_files(unit_file)//'.dat')
!________________________________________________________________________________________

! Lectura de las velocidades de la Dinamica (unit = 123)


call count_lines() ; read(123,*) Natom ; rewind(123)
N_fftw=nline/(Natom+2) ; N_fftwi=N_fftw/2 


!...............
CALL ALLOCATER()   !Allocateo todos lo arrays  
!...............


!dt=2.0e-14_dp
!dt=30.0e-15_dp
write(*,*) "dt/fs"
read(*,*) dt
dt=dt*1.0e-15_dp !paso a segundos

  w_fftw    = 1._dp/(dt*dble(N_fftw))
  
  nu_fftw   = w_fftw/v_luz
      
  fact_fftw = 2._dp/dble(N_fftw) 
   
  f_fftw_cm = 1._dp/(2._dp*nu_fftw)  

TT=0._dp
 do j=1,N_fftw
 vel2=0._dp 
  read(123,*)         !linea correspondiente al numero de atomos
  read(123,*)         !linea correspondiente a la leyenda
  do i=1,Natom
  !read(123,*)Id,OE(i),vx(i,j),vy(i,j),vz(i,j)  !Velocidades 
  read(123,*)Id,vx(i,j),vy(i,j),vz(i,j)  !Velocidades 
  
 vx(i,j)=Aps_ms*vx(i,j)     !Conversion de velocidades a m/s
 vy(i,j)=Aps_ms*vy(i,j)     !Conversion de velocidades a m/s
 vz(i,j)=Aps_ms*vz(i,j)     !Conversion de velocidades a m/s
 
 selectcase(Id)
 case('Au')
 mas_v(i)=masa_Au
 vel2=vel2 + masa_Au*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)+vz(i,j)*vz(i,j))
 case('Co')
 mas_v(i)=masa_Co             
 vel2=vel2 + masa_Co*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)+vz(i,j)*vz(i,j))
 case('Pt')
 mas_v(i)=masa_Pt
 vel2=vel2 + masa_Pt*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)+vz(i,j)*vz(i,j))
 case('Ar')
 mas_v(i)=masa_Ar
 vel2=vel2 + masa_Ar*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)+vz(i,j)*vz(i,j))
 case default
   write(*,*)'1)_ WARNING: Hay elementos que no conosco' 
 end select

 enddo
 vel2=vel2/real(Natom)
 
 TT=TT + vel2/(3._dp*kb)
 enddo 
 TT=TT/real(N_fftw)
 bhc=hp*v_luz/(kb*TT)               !Factor Q particion de solido de Debye [cm] 
 
 write(*,*)'Temperatura',TT  !Temperatura Promedio del Sistema
!______________________________________________________________________________________________________
                             
                     !Transformada Directa   (Calculo de DoS)
     
!Esta parte del programa realiza una transformada de fourier directa (FFTW3) de una dimension
!de una funcion REAL(fx) a una COMPLEJA(Fx_t)			     
     

 call dfftw_plan_dft_r2c_1d(plan_r2c,N_fftw,fx,Fx_t,FFTW_ESTIMATE)  !Genero el plan para la T directa
  
do i=1,Natom     
  
  do j=1,N_fftw
  fx(j)=vx(i,j)     !Asigno las velocidades a fx ,que va a ser la que se transforma a Fx_t
  fy(j)=vy(i,j)
  fz(j)=vz(i,j)
  enddo 
  
  call dfftw_execute_dft_r2c(plan_r2c,fx,Fx_t)  !transformacion directa
  call dfftw_execute_dft_r2c(plan_r2c,fy,Fy_t)
  call dfftw_execute_dft_r2c(plan_r2c,fz,Fz_t)
 
  do j=1,N_fftwi
  spot_x(i,j)=Fx_t(j)*conjg(Fx_t(j))*fact_fftw*fact_fftw     !Calculo del espectro de Energia 
  spot_y(i,j)=Fy_t(j)*conjg(Fy_t(j))*fact_fftw*fact_fftw
  spot_z(i,j)=Fz_t(j)*conjg(Fz_t(j))*fact_fftw*fact_fftw    
  end do
 
 
enddo
 
 DoS=0._dp;dos_x=0._dp;dos_y=0._dp;dos_z=0._dp               !Calculo de la Densidad de Estados 
 dos_core=0._dp;dos_shell=0._dp                              !DoS PROYECTADA sobre Core y Shell  

OE=12.0_dp   !Esto es un parche porque no leo OREDENES DE ENLACE, asi no se el programa cuando pide OE

do j=1,N_fftwi

  do i=1,Natom
   masa=mas_v(i)
   DoS(j)=DoS(j) + masa*( spot_x(i,j) + spot_y(i,j) + spot_z(i,j) )        !Densidad de Estados totales
   dos_x(j)=dos_x(j) + masa*spot_x(i,j)                                    !Densidad de Estados en x
   dos_y(j)=dos_y(j) + masa*spot_y(i,j)                                    !Densidad de Estados en y
   dos_z(j)=dos_z(j) + masa*spot_z(i,j)                                    !Densidad de Estados en z
 !..... Densidad de estado proyectada Sobre Core y Shell 
   if(OE(i)<10.5)then     
   dos_shell(j)=dos_shell(j)+masa*( spot_x(i,j) + spot_y(i,j) + spot_z(i,j) )
   elseif(10.5<=OE(i))then
   dos_core(j)=dos_core(j)+masa*( spot_x(i,j) + spot_y(i,j) + spot_z(i,j) )
   endif   
 !.....................................................   
  enddo  
  
  dos_x(j)=( dos_x(j)/(kb*TT) )*f_fftw_cm
  dos_y(j)=( dos_y(j)/(kb*TT) )*f_fftw_cm
  dos_z(j)=( dos_z(j)/(kb*TT) )*f_fftw_cm 
  DoS(j)=( DoS(j)/(kb*TT) )*f_fftw_cm          ! [cm] (.....  " DENSIDAD DE ESTADOS " .......... )
 !------------Densidad de Estado proyectada sobre Core y Shell---------------------
  dos_core(j) =( dos_core(j)/(kb*TT) )*f_fftw_cm
  dos_shell(j)=( dos_shell(j)/(kb*TT) )*f_fftw_cm
 !---------------------------------------------------------------------------------- 
  write(100,13) w_fftw*dble(j-1),real(DoS(j)),real(dos_x(j)),real(dos_y(j)),real(dos_z(j))
  write(101,13) nu_fftw*dble(j-1),real(DoS(j)),real(dos_x(j)),real(dos_y(j)),real(dos_z(j))
 !----------Densidad de Estados proyectada sobre Core----------- 
  write(102,*) w_fftw*dble(j-1) ,real(dos_core(j))
  write(103,*) nu_fftw*dble(j-1),real(dos_core(j))
 !---------Densidad de Estados proyectada sobre Shell-----------
  write(104,*) w_fftw*dble(j-1) ,real(dos_shell(j))
  write(105,*) nu_fftw*dble(j-1),real(dos_shell(j))
 !-------------------------------------------------------------

 enddo 

 call dfftw_destroy_plan(plan_r2c)

 !______________________________________________________________________________________________

                          !Transformada Inversa  (Calculo de Autocorrelacion)

!Esta parte del programa realiza una transformada de fourier inversa (FFTW3) de una dimension
!de una funcion COMPLEJA(DoS) a una COMPLEJA(Rxyz)


 call dfftw_plan_dft_1d(plan_c2c,N_fftwi,DoS,Rxyz,FFTW_BACKWARD,FFTW_ESTIMATE)!Genero el plan para la T inversa
 
 call dfftw_execute_dft(plan_c2c,DoS,Rxyz)
 call dfftw_execute_dft(plan_c2c,dos_x,Rxx)
 call dfftw_execute_dft(plan_c2c,dos_y,Ryy)
 call dfftw_execute_dft(plan_c2c,dos_z,Rzz)
 
 t=0._dp
 do j=1,N_fftwi
  write(70,*)t,0.5_dp*real(Rxyz(j)) !,0.5_dp*real(Rxx(j)),0.5_dp*real(Ryy(j)),0.5_dp*real(Rzz(j))
  t=t + 2._dp*dt
 enddo
 
 call dfftw_destroy_plan(plan_c2c) 

!---------------------------------------------------------
!___________________   2PT Model    ______________________

So=DoS(1) !Obtencion de So=S(0)

 CDN=fun_CDN(TT,p,Natom,masa,So,kb,v_luz) !Obtencion de Constante de Difusividad Normalizada (CDN)
 
 call fraccion_gas(CDN,f)  !Calculo fraccion de Fase tipo Gas (f)

 write(*,*)'Coef. Dif. Normalizada',CDN ; write(107,*)'Coef. Dif. Normalizada',CDN
 write(*,*)'fraccion de fase tipo gas',f ; write(107,*)'fraccion de fase tipo gas',f
 
 Ycs=(f**2.5_dp)/(CDN**1.5_dp)       !Factor de Compresibilidad de Carnahan-Starling								   								   
 ZCS=fun_ZCS(Ycs)                    !Ec. Estado Carnahan-Starliing

!______________________________________ Cálculo de Factores de Peso (Fase tipo GAS) (NO dependen de la frecuencia)

 WGS=Wg_S(masa,kb,hp,TT,f,p,ZCS,Ycs) !Factor de peso de Entropia de Gas (Weighting_gas_Entropia)(Wg_S)
 WGE=0.5_dp                          !Factor de peso de Energia de  Gas
 WGA=WGE-WGS                         !Factor de peso de Energia libre de Gas (Helmholtz)

!______________________________________ DoS (GAS) ____________________________

do j=1,N_fftwi
 
 vnu=nu_fftw*dble(j-1)   !Frecuenicas (numero de onda (cm-1)) 

 DoS_gas(j)=So/(1._dp + (pi*So*vnu/(6._dp*f*dble(Natom)))**2)  !DoS de la fase Gaseosa 
 
 DoS_S_Debye(j)=dabs( real(DoS(j))-DoS_gas(j) )                !DoS de la fase Solida

 write(106,57)vnu,DoS_gas(j),DoS_S_Debye(j),real(DoS(j)) 

  !---- DoS pesadas por factor de Entropia 
 g_SWSS(j)=DoS_S_Debye(j)*Ws_S(bhc,vnu)    !DoS(fase solida)*factor de peso de Entropia del solido
 g_SWGS(j)= DoS_gas(j)*WGS                 !DoS(fase gas)*factor de peso de Entropia del gas
 g_SDoS(j)=real(DoS(j))*Ws_S(bhc,vnu)      !DoS(TOTAL)*factor de peso de Entropia del solido
 
 !---- DoS pesadas por factor de Helmholtz
 g_SWSA(j)=DoS_S_Debye(j)*Ws_A(bhc,vnu)    !DoS(fase solida)*factor de peso de Energia de Helmholtz
 g_SWGA(j)=DoS_gas(j)*WGA                  !DoS(fase gas)*factor de peso de Energia de Helmholtz 
 g_ADoS(j)=real(DoS(j))*Ws_A(bhc,vnu)      !DoS(TOTAL)*factor de peso de Energia de Helmholtz
 
enddo

 call Cuad_Simps(N_fftwi,nu_fftw,g_SWGS,SWGS)   !calculo de cudratura para el calculo de la entropia fase gas 
 call Cuad_Simps(N_fftwi,nu_fftw,g_SWSS,SWSS)   !calculo de cudratura para el calculo de la entropia de fase solido
 call Cuad_Simps(N_fftwi,nu_fftw,g_SDoS,SDoS)   !calculo de cuadratura para el calculo de entropia  (sin 2PT)

 call Cuad_Simps(N_fftwi,nu_fftw,g_SWGA,SWGA)   !calculo de cudratura para el calculo de Energia de Helmholtz fase gas
 call Cuad_Simps(N_fftwi,nu_fftw,g_SWSA,SWSA)   !calculo de cudratura para el calculo de Energia de Helmholtz fase solido
 call Cuad_Simps(N_fftwi,nu_fftw,g_ADoS,ADoS)   !calculo de cuadratura para el calculo de Energia de Helmholtz (sin 2PT)

S2PT=(SWSS+SWGS)/real(Natom)                                        !SWSS(entropia fase solida) ; SWGS(entropia fase gaseosa)
SDoS=SDoS/real(Natom)                                               !SDoS(entropia desde DoS Total, sin 2PT)

A2PT= (Etotal-3._dp*real(Natom)*(1._dp-0.5_dp*f)*Temp+Temp*(SWSA+SWGA))/real(Natom)  !SWSA(E.Helmholtz fase solida) ; SWGA(E. Helmholtz fase gaseosa)
ADoS= (Etotal-3._dp*real(Natom)*Temp+Temp*ADoS)/real(Natom)                          !ADoS(Helmholtz desde DoS Total, sin 2PT)


 !! Ver las unidades de esto de A2PT  y ADoS, y como usar la temperatura con sus
 !unidades


write(*,*)'Entropia',S2PT        !Entropia total del sistema (normalizada por el numero de atomos)
write(107,*)'Entropia',S2PT      !Entropia total del sistema (normalizada por el numero de atomos)
write(*,*)'Entropia_DOS',SDoS
write(107,*)'Entropia_DOS',SDoS  !Entropia sin usar el modelo 2PT


write(*,*)'E. Helmholtz',A2PT        !E. Helmholtz del sistema (normalizada por el numero de atomos)
write(107,*)'E. Helmholtz',A2PT      !E. Helmholtz del sistema (normalizada por el numero de atomos)
write(*,*)'E. Helmholtz',ADoS
write(107,*)'E. Helmholtz',ADoS      !E. Helmholtz sin usa el modelo 2PT

!-------------------------------------------
13  format(f19.4,x,4(f16.10,x))
57  format(4(f14.8,x))
!-------------------------------------------

stop
END SUBROUTINE doDoS

subroutine bond_order(border,n_oe)
implicit none
integer,parameter :: sp=kind(0.)
real(sp),intent(in) :: border
integer,intent(out) :: n_oe

!if(0.5 <=border.and.border< 1.5) n_oe=1
!if(1.5 <=border.and.border< 2.5) n_oe=2
!if(2.5 <=border.and.border< 3.5) n_oe=3     
!if(3.5 <=border.and.border< 4.5) n_oe=4        
!if(4.5 <=border.and.border< 5.2) n_oe=5        
if(4.0 <=border.and.border< 7.0)  n_oe=6       
!if(6.8 <=border.and.border< 7.2) n_oe=7        
if(7.0 <=border.and.border< 10.0) n_oe=8        
!if(8.9 <=border.and.border< 9.5) n_oe=9       
!if(9.5 <=border.and.border<10.5) n_oe=10        
!if(10.5<=border.and.border<11.2) n_oe=11      
if(10.0<=border.and.border<14.0) n_oe=12       
end subroutine  

END MODULE mod_DoS_2PT
