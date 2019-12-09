!----------------------------------------------------------------------------------------------------------
!Este codigo Fortran 90 usa la libreria FFTW-3.2.2 previamente instaladas en usr/local/lib:
!esta libreria fue bajada de la siguiente direccion :
!      www.fftw.org/download.html
!e instalada ejecutando los sig. comandos en el directorio donde fue guardado el paquete.
! ./configure
! make
! make install 
!---------------------------------------------------------------------------------------------------------- 
!Se compila de la siguiente forma :
!ifort -O0 -traceback biseccion.f90 weighting_functions.f90 prueba_fftw_1.f90 
!-L/usr/local/lib -lfftw3 -lm -o ex.autocorr

! ./compiler   !Compilacion (compiler es un script de compilacion)
!./ex.autocorr !Ejecucion

!----------------------------------------------------------------------------------------------------------	     
!____________________________________________________________________________________        Martin Ludueña 

                                program Autocorrelacion
!__________________________________________________________________________________________________________                                
!----------------------------------------------------------------------------------------------------------
!                               Modulos que Utiliza el Programa
 USE biseccion
 USE Simpson
 USE weighting_functions
!----------------------------------------------------------------------------------------------------------
implicit none
!----------------------------------------------------------------------------------------------------------
 include "fftw3.f" 
!---------------------------------------------------------------------------------------------------------- 
 integer,parameter :: dp=kind(0.d0)
!---------------------------------------------------------------------------------------------------------- 
 integer,parameter :: N_fftw=734                                  !Grillado de la f a transformar
 integer,parameter :: N_fftwi=N_fftw/2 +1 
 integer,parameter :: Natom=500                                    !Numero de particulas 
!----------------------------------------------------------------------------------------------------------
 integer :: plan_r2c                                               !Plan Real a Complejo
 integer :: plan_c2c                                               !Plan Complejo a Complejo
!----------------------------------------------------------------------------------------------------------       
 real(dp),dimension(1:N_fftw) :: fx ,fy,fz                         !Funcion a transformar
 complex(dp),dimension (1:N_fftwi) :: Fx_t                         !transformada de la fx
 complex(dp),dimension (1:N_fftwi) :: Fy_t                         !transformada de la fy
 complex(dp),dimension (1:N_fftwi) :: Fz_t                         !transformada de la fz
!----------------------------------------------------------------------------------------------------------
 complex(dp),dimension (1:N_fftwi) :: Rxyz                         !Funcion de Autocorrelacion Total
 complex(dp),dimension (1:N_fftwi) :: Rxx                          !Funcion de Autocorrelacion de Vx 
 complex(dp),dimension (1:N_fftwi) :: Ryy                          !Funcion de Autocorrelacion de Vy
 complex(dp),dimension (1:N_fftwi) :: Rzz                          !Funcion de Autocorrelacion de Vz
!----------------------------------------------------------------------------------------------------------
 complex(dp),dimension(1:Natom,1:N_fftwi) :: spot_x,spot_y,spot_z  !Espectro de Energia (/Fx_t/**2)
 complex(dp),dimension(1:N_fftwi) :: DoS,dos_x,dos_y,dos_z         !Densidad de Estados
 real(dp),dimension(1:N_fftwi)    :: DoS_gas,DoS_S_Debye           !Densidad de Estados 2PT model [cm]
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------

                       !Two Phase Model (Cálculo de propiedades termodinamicas) 
             
 real(dp) :: CDN       !Coeficiente de difucion normalizado [Adimensional]
 real(dp) :: So        !DoS(0)  [cm] 
 real(dp) :: f         !fraccion de gas  (factor de fluidez)
 real(dp) :: Ycs       !Factor de Compresibilidad de Esferas Duras (Carnahan - Starling)
 real(dp) :: ZCS       !Funcion de estado de Esferas Duras (Carnahan - Starling)
!----------------------------------------------------------------------------------------------------------
 real(dp)  ::  WGS,WGE,WGA                           !Factor de Peso Fase GAS 
!----------------------------------------------------------------------------------------------------------      
 real(dp),dimension(1:N_fftwi) :: g_SWGS             !DoS(Tipo Gas) * WGS
 real(dp),dimension(1:N_fftwi) :: g_SWSS             !DoS(Tipo Solido) * Ws_S (Ws_S = funcion de peso Entropia solido)
 real(dp),dimension(1:N_fftwi) :: g_SDoS             !DoS(Total)*Ws_S 
  
 real(dp),dimension(1:N_fftwi) :: g_SWGA             !DoS(Tipo Gas) * WGA
 real(dp),dimension(1:N_fftwi) :: g_SWSA             !DoS(Tipo Solido) * Ws_A (Ws_A = funcion de peso E. Helmholtz solido)
 real(dp),dimension(1:N_fftwi) :: g_ADoS             !DoS(Total)*Ws_A 
  
 !----------------------------------------------------------------------------------------------------------
 real(dp)  :: SWGS                                   !Entropia Fase tipo Gas
 real(dp)  :: SWSS                                   !Entropia Fase tipo Solido
 real(dp)  :: S2PT                                   !Entropia del Sistema (con 2PT model) 
 real(dp)  :: SDoS                                   !Entropia del Sistema (sin 2PT model)    
  
 real(dp)  :: SWGA                                   !E. Helmholtz Fase tipo Gas
 real(dp)  :: SWSA                                   !E. Helmholtz Fase tipo Solido
 real(dp)  :: A2PT                                   !E. Helmholtz del Sistema (con 2PT model) 
 real(dp)  :: ADoS                                   !E. Helmholtz del Sistema (sin 2PT model)    
 !---------------------------------------------------------------------------------------------------------- 
 real(dp),dimension(1:Natom,1:N_fftw) :: vx,vy,vz                   !Velocidades Reducidas en x,y,z
!----------------------------------------------------------------------------------------------------------
     
                               !Parametros de L-J  (Argon) 

real(dp),parameter :: eps=1.65e-21_dp                             !Pozo de potencial  [J]
real(dp),parameter :: sig=3.405e-10_dp                            !Distancia de referencia del potencial [m]
real(dp),parameter :: masa=6.63e-26_dp                            !Masa del Argon [kg]
!-----------------------------------------------------------------------------------------------------------

                               !Constantes Universales
       
real(dp),parameter :: kb= 1.38e-23_dp           !Constante de Boltzman  [J/k]   		       
real(dp),parameter :: v_luz=29979245800._dp     !velocidad de la luz [cm/s]
real(dp),parameter :: hp=6.62606896e-34_dp      !Constante de planck (cuanto elemental de acción) [Js]
real(dp),parameter :: pi=3.141592654_dp         !Numero pi
real(dp) :: bhc                                 !hp*v_luz/KT [cm]
!-----------------------------------------------------------------------------------------------------

                              !Propiedades Mecanicas/Termodinamicas
 
real(dp),parameter :: densidad=1.1_dp       !Densidad Reducida
real(dp),parameter :: Temp=1.1_dp           !Temperatura Reducida 
real(dp) :: Etotal=-2295.998_dp             !Energia total reducida   
real(dp) :: p                               !Densidad    (1/m3)
real(dp) :: TT                              !Temperatura (K)

!----------------------------------------------------------------------------------------------------- 
 integer  :: i,j                                    !Contadores i=>(Natom) ,j=>(N_fftw)
 real(dp) :: dt=0.0175_dp                           !paso de tiempo reducido(tiempo_total/Nª_datos)           
 real(dp) :: t                                      !tiempo reducido
 real(dp) :: vnu                                    !numero de onda [1/cm]
!------------------------------------------------------------------------------------------------------

                                !Factores de Escaleo 
 
 real(dp) :: w_fftw             ! (escalea a frecuncia (Hz))
 real(dp) :: nu_fftw            ! (escalea a numero de onda (cm-1))
 real(dp) :: fact_fftw          ! Con este valor normalizo para que me de los coeficientes de la fx
 real(dp) :: f_fftw_cm          ! (Dimensionaliza DoS a cm)
  
  w_fftw    = 1._dp/(dt*dble(N_fftw))
      
  nu_fftw   = w_fftw*dsqrt(eps/(masa*sig*sig))/v_luz  !Al trabajar con unidades reducidas
    
  fact_fftw = 2._dp/dble(N_fftw) 
    
  f_fftw_cm = 1._dp/(2._dp*nu_fftw)  

!--------------------------------------------------------------------------------------------------------  
  
  p = densidad/(sig*sig*sig)         !Densidad del sistema particular [1/m3]
  TT = (eps/kb)*Temp                 !Temperatura del sistema particular [K]
    
  bhc=hp*v_luz/(kb*TT)               !Factor Q particion de solido de Debye [cm] 
  
  
  write(*,*)'Densidad', densidad   ; write(*,*)'T',Temp 
!______________________________________________________________________________________________________

open(40,file='Ac.xyz')          !Input (Velocidades de las particulas)
open(70,file='OUTPUT/VAC/Ac.dat')       !Output; escribe :: t,Rxyz,Rxx,Ryy,Rzz
open(100,file='OUTPUT/DoS/Ac.dat')      !Output; escribe :: vnu,DoS,dos_x,dos_y,dos_z  
open(250,file='OUTPUT/DoS_2PT/Ac.dat')  !Output; escribe :: vnu,DoS_gas
open(150,file='OUTPUT/Info/Ac.dat')          
!______________________________________________________________________________________________________
    
          ! Lectura de las velocidades de la Dinamica (unit = 40)
 
 do j=1,N_fftw
 do i=1,Natom
  read(40,*)t,vx(i,j),vy(i,j),vz(i,j)  
 enddo
 enddo 
 
!______________________________________________________________________________________________________
                             
                            !Transformada Directa   (Calculo de DoS)
                             
!Esta parte del programa realiza una transformada de fourier directa (FFTW3) de una dimension
!de una funcion REAL(fx) a una COMPLEJA(Fx_t)			     
     

 call dfftw_plan_dft_r2c_1d(plan_r2c,N_fftw,fx,Fx_t,FFTW_ESTIMATE)  !Genero el plan para la Trans directa
  
do i=1,Natom     
  
  do j=1,N_fftw
  fx(j)=vx(i,j)     !Asigno las velocidades a fx ,que va a ser la que se transforma a Fx_t
  fy(j)=vy(i,j)
  fz(j)=vz(i,j)
  enddo 
  
  call dfftw_execute_dft_r2c(plan_r2c,fx,Fx_t)  !transformacion directa
  call dfftw_execute_dft_r2c(plan_r2c,fy,Fy_t)
  call dfftw_execute_dft_r2c(plan_r2c,fz,Fz_t)
 
  do j=1,N_fftw/2 +1
  spot_x(i,j)=Fx_t(j)*conjg(Fx_t(j))*fact_fftw*fact_fftw     !Calculo del espectro de Energia 
  spot_y(i,j)=Fy_t(j)*conjg(Fy_t(j))*fact_fftw*fact_fftw
  spot_z(i,j)=Fz_t(j)*conjg(Fz_t(j))*fact_fftw*fact_fftw    
  end do
 
 
enddo
 
 DoS=0._dp;dos_x=0._dp;dos_y=0._dp;dos_z=0._dp                  !Calculo de la Densidad de Estados 
 do j=1,N_fftwi
 do i=1,Natom
  DoS(j)=DoS(j) + spot_x(i,j) + spot_y(i,j) + spot_z(i,j)       !Densidad de Estados totales
  dos_x(j)=dos_x(j) + spot_x(i,j)                               !Densidad de Estados en x
  dos_y(j)=dos_y(j) + spot_y(i,j)                               !Densidad de Estados en y
  dos_z(j)=dos_z(j) + spot_z(i,j)                               !Densidad de Estados en z
 enddo  
 
 DoS(j)=(DoS(j)/Temp)*f_fftw_cm      ! [cm] (.....  " DENSIDAD DE ESTADOS " .......... )
 
 write(100,13)nu_fftw*dble(j-1),real(DoS(j)),real(dos_x(j)),real(dos_y(j)),real(dos_z(j))
 enddo 

 
 call dfftw_destroy_plan(plan_r2c)
!______________________________________________________________________________________________

                  ! Transformada Inversa  (Calculo de Autocorrelacion)

!Esta parte del programa realiza una transformada de fourier inversa (FFTW3) de una dimension
!de una funcion COMPLEJA(DoS) a una COMPLEJA(Rxyz)


 call dfftw_plan_dft_1d(plan_c2c,N_fftwi,DoS,Rxyz,FFTW_BACKWARD,FFTW_ESTIMATE)!Genero el plan para la T inversa
 
 call dfftw_execute_dft(plan_c2c,DoS,Rxyz)
 call dfftw_execute_dft(plan_c2c,dos_x,Rxx)
 call dfftw_execute_dft(plan_c2c,dos_y,Ryy)
 call dfftw_execute_dft(plan_c2c,dos_z,Rzz)
 
 t=0._dp
 do j=1,N_fftwi
 write(70,13)t,0.5_dp*real(Rxyz(j)),0.5_dp*real(Rxx(j)),0.5_dp*real(Ryy(j)),0.5_dp*real(Rzz(j))
 t=t + 2._dp*dt
 enddo
 
 call dfftw_destroy_plan(plan_c2c) 

!_______________________________________________________________ 2PT Model ______________________

 So=DoS(1) !Obtencion de So=S(0)

 CDN=fun_CDN(TT,p,Natom,masa,So,kb,v_luz) !Obtencion de Constante de Difusividad Normalizada (CDN)
 
 call fraccion_gas(CDN,f)  !Calculo fraccion de Fase tipo Gas (f)

 write(*,*)'Coef. Dif. Normalizada',CDN ; write(150,*)'Coef. Dif. Normalizada',CDN
 write(*,*)'fraccion de fase tipo gas',f ; write(150,*)'fraccion de fase tipo gas',f
 
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

 write(250,57)vnu,DoS_gas(j),DoS_S_Debye(j),real(DoS(j)) 

 

   
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


write(*,*)'Entropia',S2PT        !Entropia total del sistema (normalizada por el numero de atomos)
write(150,*)'Entropia',S2PT      !Entropia total del sistema (normalizada por el numero de atomos)
write(*,*)'Entropia_DOS',SDoS
write(150,*)'Entropia_DOS',SDoS  !Entropia sin usar el modelo 2PT


write(*,*)'E. Helmholtz',A2PT        !E. Helmholtz del sistema (normalizada por el numero de atomos)
write(150,*)'E. Helmholtz',A2PT      !E. Helmholtz del sistema (normalizada por el numero de atomos)
write(*,*)'E. Helmholtz',ADoS
write(150,*)'E. Helmholtz',ADoS      !E. Helmholtz sin usa el modelo 2PT


13  format(5(f14.8,x))
57  format(4(f14.8,x))
end program 









