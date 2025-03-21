program main
use Domain
use particles
use ini
use properties
use kernal
use ptcl_matrix
use density
use pressure
use Compute_press_grad
use CompStress
use CompStrainRate
use CompGravity
use kernal_gradient_correction
use wavelength
use Damping_Zone
use step
implicit none

double precision::x!length of domain
double precision::z 
double precision::dp!particle spacing
double precision::dp_fly
double precision::wmaker_posx!x-position of 2d-wavemaker
double precision::fl_length !!length of fluid excluding damping zone
double precision::fl_height !fluid depth
double precision::wavemaker_height !wavemaker height
double precision::ini_posx,ini_posy,ini_posx_piston,ini_posy_piston
double precision::Xmin,Ymin
double precision::dt!time step
double precision::coeffs!coeff of sound
double precision,parameter::vis=1d-6
double precision::nu!kinematic viscosity of water
double precision::a!artificial viscosity of fluid
double precision::rho0 !density of fluid particles
double precision::rhob0!density of boundary particles
double precision::h,time,time_half,diff,diffx,diffy
integer::nx,nz,ntotal,ntotal_f,nx_fl,nz_fl,i,j
integer::nz_wavemaker
integer::n
type(particle),allocatable::p(:)
double precision,allocatable::drodt(:),drodt_old(:),dens_old(:)
integer::rank,ierr,status,ierror
integer::nt,t
double precision,allocatable::vel_halfx(:),vel_halfy(:),pos_half_x(:),pos_half_y(:),acc_oldx(:),acc_oldy(:)
double precision,allocatable::vel_oldx(:),vel_oldy(:),pos_oldx(:),pos_oldy(:),corr_vx(:),corr_vy(:)
double precision,allocatable::ini_coord_x(:),ini_coord_y(:)
type(vector_2d),allocatable::pres_vel(:)!prescribed wall velocity
type(vector_2d),allocatable::extra_vel(:)!extrapolated wall velocity
type(vector_2d),allocatable::fun(:)!reduction function
type(vector_2d)::damp_i,damp_f
double precision::simtime!simulation  time
double precision::var1,dp2,dp_res,s1,s2,s3
double precision::ramp_period 
double precision::ramp_func
character(13)::Wtype
double precision::cff
integer::dt_option

integer::type                !type 1=laminar +sps turbulnce viscosity                                      
                             !type 2=artificial viscosity
                             !type 3=orginal viscossity

 integer::k_type          !type 1=quintic spline kernal
                          !type 2=wendland kernal

integer::time_integrator=1  !1=predictor-corrector time integration
                            !2=velocity-verlet time integration

integer::wavemaker_theory       !1=1st order wave maker theory
                                  !2=2nd order wavemaker theory
                                  
integer::bound_type           !type 1=slip-boundary condition
                              !type 2=no-slip boundary condition                                            
character::err_msg
double precision::out_time
integer::t1,t2,clock_rate,clock_max
real::a1,b1,c1
integer::day,hr,mint,out_nt
real::sec







!properties of incident wave that is to be generated
double precision::wave_period 
double precision::wave_length 
double precision::wave_height 
double precision::dmp_size
double precision::dmp_ini 
double precision::dmp_final  !minimum 0.5 to 1.0 wavelngth size of damping zone is required to damp out incoming wave
double precision::beta !reduction parameter
double precision::f !angular frequency
double precision::k !wavenumber
double precision,parameter::pi=3.1416
double precision::S0,m1
double precision::fl_total_length



!minimum y position which is 0.0 for both flume bed and piston. Changing this parameter other than 0.0 might cause program to crash.
Ymin=0.0d0


write(*,20)"2d regular wave generation by piston type wavemaker using SPH numerical method copyright(c) 2025 by Md.Jahid Hasan,M.S.C at Water resources engineering,Bangladesh University of Engineering and Technology,student id-0419162024"
20 format(/,2x,a)
write(*,21)"#################################################################################################################################################################################################################################"
21 format(2x,a)
write(*,22) "This 2d fortran program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License"
write(*,22) "as published by the Free Software Foundation either version 3 of the License, or (at your option) any later version"
  
write(*,22) "This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE."
write(*,22) "See the General Public License for more details.You should have received a copy of the GNU General Public Licensealong with this program.If not, see <https://www.gnu.org/licenses/>."

   
22 format(20x,a)
write(*,23)"########################################################################################################################################################################################"
23 format(2x,a)
write(*,24)"All values must be given in SI units"
24 format(2x,a80)
write(*,25)"*************************************************************"
25 format(15x,a80)

write(*,43)"Incident Wave properties which is to be generated by 2d piston type wavemaker"
write(*,44)"**********************************************************************************"
43 format(/,1x,a)
44 format(1x,a)
write(*,*)"Order of wave generation=?(1=1st order,2=2nd order)"
read(*,*) wavemaker_theory
write(*,*)"Target wave period T=?"
read(*,*) wave_period
write(*,*)"Incident wave height=?"
read(*,*) wave_height
write(*,32)"Initial water depth=?"
32 format(/,1x,a)
read(*,*) fl_height
write(*,*)"Computing wavelength of the incident wave........."
call ComputeWavelength(wave_period,fl_height,wave_length,Wtype)
write(*,47)"Done"
47 format(1x,a30)
write(*,48)"********************"
48 format(1x,a38)
write(*,50)"This is a",Wtype,"water depth condition."
50 format(1x,a,1x,a,a)
write(*,*)"Length of main fluid domain along wave flume excluding damping zone=?"
read(*,*) fl_length
write(*,45)"Note:its assuming that size of damping zone is an integer multiple of wave length(i.e 1L,2L,3L......nL where L is wavelength)"
45 format(/,2x,a)
write(*,*)"Size of daming zone(only integer multiplier i.e 1,2,3 etc)"
read(*,*) n

write(*,54)"Note:Damping zone is assumed to be a rectangular shape"
54 format(/,1x,a)

dmp_ini=wmaker_posx+fl_length
dmp_size=n*wave_length
dmp_final=dmp_ini+dmp_size

damp_i%x(1)=dmp_ini
damp_f%x(1)=dmp_final


write(*,52)"Initial or minimum x-position of damping zone=",dmp_ini
52 format(/,1x,a,1x,f8.5)
write(*,53)"Final or maximum x-position of damping zone=",dmp_final
53 format(1x,a,1x,f8.5)
write(*,*)"Minimum y-position of damping zone=?"
read(*,*) damp_i%x(2)
write(*,*)"Maximum y position of damping zone=?"
read(*,*) damp_f%x(2)

write(*,49)"Damping reduction coefficient=?(coefficient value of 10 is suggested higher value can also be used)"
49 format(1x,a)
read(*,*) beta
write(*,46)"value of Ramp period=?(1.0T,2.0T,3.0T.....nT where T is incident wave period and n is real or integer positive number)"
46 format(1x,a)
read(*,*) ramp_period

!computing parameters related to wave generation
f=(2.0*pi)/wave_period !angular frequency
k=(2.0d0*pi)/wave_length !wave number
m1=(2.0d0*(sinh(k*fl_height))**2)/((sinh(k*fl_height)*cosh(k*fl_height))+(k*fl_height))
S0=wave_height/m1

   if(wavemaker_theory.gt.1)then !second order wavemaker theory(madsen 1971)
     s1=3.0*cosh(k*fl_height)
     s2=(sinh(k*fl_height))**3
     s3=wave_height**2/(32.0*fl_height)
   var1=s3*((s1/s2) -(2.0/m1))
   else
  var1=0.0d0
  end if


write(*,26)"spacing between particles(dp)=?"
26 format(/,1x,a32)
write(*,27)"Note:dp should be choosen in such a way that ratio of length of waveflume in x direction and dp is an integer"
27 format(/,1x,a110)
read(*,*) dp

write(*,500)"coefficient to calculate Smoothing length(h)=?(generally h=k*dp where k is real numberand value of k=1.0-1.3 is recomended)"
500 format(1x,a)
read(*,*) cff
  h=cff*dp


write(*,*)"Minimum x-position of domain,Xmin=?"
read(*,*) Xmin


write(*,28)"Note:Its assuming that bottom of flume bed is horizontal and minimum y-position of bottom bed is zero same as minimum y-position of wavemaker"
28 format(/,1x,a)

write(*,*)"length of wave flume in y-direction,y=?"
read(*,*) z
write(*,*)"Initial x-position of pistion=?"
read(*,*) wmaker_posx
if(wmaker_posX.le.Xmin)then
  write(*,*)"Error:X-position of piston must be greater than minumum X-position of domain"
  write(*,*)"Press any key to exit.."
  read(*,*)
  stop
  endif
write(*,*)"Height of wavemaker or piston=?"
read(*,*) wavemaker_height
 




write(*,544)"density of fluid particles=?"
544 format(/,1x,a)
read(*,*) rho0
write(*,29)"density of wall/boundary particles=?"
29 format(/,1x,a)
write(*,30)"Note:its suggessed that density of wall particles should be equal to that of fluid particles."
30 format(1x,a93)
write(*,31)"higher density cause higher repulsive force excerted by wall to fluid partlces which results in gap betn wall and fluid."
31 format(1x,a120)
read(*,*) rhob0
write(*,33)"type of kernal function that will be used =?(1=quintic spline kernal,2=wendland kernal)"
33 format(/,1x,a)
read(*,*) k_type
write(*,34)"coefficient to calculate artificial sound speed in fluid=?(10 to 20 gives good approximation of sound)"
34 format(/,1x,a)
write(*,158)"Note:Higher value of coefficent results in small time-step size which further increases no.of time-step."
158 format(1x,a)
read(*,*) coeffs
write(*,36)"which type of viscosity will be used=?(1=laminar+sps turbulence,2=artificial viscosity,3=orginal viscosity)"
36 format(/,1x,a)
read(*,*) type
if(type.eq.1.or.type.eq.3)then
  write(*,*)"kinematic viscosity of fluid=?"
  read(*,*) nu
elseif(type.eq.2)then
 write(*,37)"artificial viscosity coefficient(alpha)=?(0.01 to 0.04 give good approximation of fluid viscosity)"
 37 format(1x,a98)
 read(*,*) a
 end if
 
 write(*,38)"Note:wall boundary is modelled by using dummy particles(paper:A generalized wall boundary condition for smoothed particle hydrodynamics by Adami 2012)"
 38 format(/,1x,a)
write(*,39)"type of boundary condition=?(1-slip,2-no-slip)"
39 format(/,1x,a)
write(*,40)"Note:Its better to use slip bc when using artificial viscosity,otherwise particles sticking to wall/boundary can occur."
40 format(1x,a)
write(*,41)"In case of orginal viscosity as kgf-sph formulation is used,no-slip boundary condition causes program to crash after some timesteps(may be its related to  matrix formulation for boundary particles,and this problm is still under investigation)"
write(*,42)"so its advisable to use slip boundary condition when using orginal viscosity formulation"
 41 format(1x,a)
 42 format(1x,a)
read(*,*) bound_type





 fl_total_length=fl_length+dmp_size !total fluid length including damping zone
 x=(wmaker_posx-Xmin)+fl_total_length!total length of domain in x direction




 
dp2=dp !cutoff radius for boundary particle creation
  
nx=floor(x/dp)
nz=floor(z/dp2)
nz_wavemaker=floor(wavemaker_height/dp2)
nx_fl=nint(fl_total_length/dp)
nz_fl=ceiling(fl_height/dp)
ntotal_f=(nx_fl-5)*nz_fl
ntotal=(3*nx+3*nz+4*nz_wavemaker-9)+ntotal_f
allocate(p(ntotal))
allocate(vel_halfx(ntotal),vel_halfy(ntotal),pos_half_x(ntotal),pos_half_y(ntotal),stat=status,errmsg=err_msg)
allocate(vel_oldx(ntotal),vel_oldy(ntotal),pos_oldx(ntotal),pos_oldy(ntotal),corr_vx(ntotal),corr_vy(ntotal),stat=status,errmsg=err_msg)
allocate(ini_coord_x(3*nz_wavemaker),ini_coord_y(3*nz_wavemaker))
allocate(drodt(ntotal),drodt_old(ntotal),dens_old(ntotal),acc_oldx(ntotal),acc_oldy(ntotal),stat=status,errmsg=err_msg)
allocate(pres_vel(ntotal),extra_vel(ntotal))
allocate(fun(ntotal))


!measuring time
call system_clock(t1,clock_rate,clock_max)

write(*,147) "Generating particles coordinate(fluid,wall boundary.etc)...................................."
147 format(1x,a100)

!creation of problem domain
 call create_domain(p,ntotal,Xmin,x,fl_length,z,wmaker_posx,wavemaker_height,fl_height,dmp_size,dp,dp2)
 
write(*,148)"done"
148 format(50x,a5)
write(*,180)'*****************************************************'
180 format(1x,a80)




write(*,187)"total no. of fluid particles=",ntotal_f
187 format(/,1x,a,i7)
write(*,188)"total no. of boundary particles=",(ntotal-ntotal_f)
188 format(1x,a,i7)
write(*,191)"total no. of particles=",ntotal
191 format(1x,a,i7)
write(*,182)"writing initial particles coordinate to a text file..............................."
182 format(/,1x,a80)


!save geometry to txt file
open(unit=10,file='coordinate_ini.txt',action='write'&
   ,status='replace',iostat=ierror)
     
    write(10,102) 'coordx','coordy','coordz'
    102 format(1x,a12,2x,a12,2x,a12)
   do i=1,ntotal
    write(10,105) p(i)%coord%x(1),p(i)%coord%x(2),p(i)%coord%x(3)
      
    105 format(1x,f12.5,2x,f12.5,2x,f12.5)
  end do

write(*,149)"done"
149 format(50x,a5)
write(*,181)'*****************************************************'
181 format(1x,a80)




!checking if fluid particles are too close to boundary particles 
 do i=1,ntotal
   if(p(i)%id.eq.0)then!fluid particles
   do j=1,ntotal
   if(p(i)%id*p(j)%id.lt.0)then!boundary particles
    diffx=p(i)%coord%x(1)-p(j)%coord%x(1)
    diffy=p(i)%coord%x(2)-p(j)%coord%x(2)
    diff=diffx*diffx+diffy*diffy
    diff=sqrt(diff)
    
    if(diff.lt.0.25*dp)then
      write(*,*)"Error:fluid particles are too close to boundary particles"
      write(*,*)"Press any key to exit..."
      read(*,*)
      stop
     end if
      end if
      end do
       end if
         end do



write(*,150)"Assigning properties(mass,density) and initializing variables(velocity,pressure,etc)to all particles in the domain......................"
150 format(/,2x,a120)

!assign properties(mass,density) to each of particles in domain
 call PROP(p,ntotal,rho0,rhob0,dp,dp2)
 

 !initialize velocity,acceleration for all particles in the domain
 call initialize(p,ntotal)

write(*,151)"done"
151 format(50x,a5)
write(*,183)"**********************************************************"
183 format(1x,a80)

write(*,600)"choose the option to determine the time step size(dt)=?(1=value is given by the user,2=dt is automatically calculated by the program)"
600 format(1x,a)
read(*,*) dt_option
 if(dt_option .eq.1)then
  write(*,*)"Time-step size(dt)=?"
  read(*,*) dt
  end if 
                                                                  
write(*,*)"Simulation time(s)=?"
read(*,*) simtime
write(*,*)"At which time(s)output results(coordinates only) will be saved=?"
read(*,*) out_time
 
  if(dt_option.gt.1)then
  call timestep(1,p,ntotal,h,fl_height,coeffs,vis,dt) 
  end if
  write(*,130)"Initial timestep size=",dt,"sec"
  130 format(1x,a28,ES12.4,a5)
  
  write(*,*)"*****************************************************"
  nt=int(simtime/dt)
  write(*,152)"Total no. of timestep=",nt  
  152 format(1x,a30,1x,i8)
  write(*,*)"****************************************************"
  out_nt=int(out_time/dt)

 
 

  !initialize variables for all particles at time=0
  do i=1,ntotal
   if(p(i)%id.eq.0)then !only fluid particles
   vel_oldx(i)=0.0d0
   vel_oldy(i)=0.0d0
   pos_oldx(i)=p(i)%coord%x(1)
   pos_oldy(i)=p(i)%coord%x(2)
   drodt(i)=0.0d0
   drodt_old(i)=0.0
   dens_old(i)=p(i)%dens!initialize old_density of fluid particles
   corr_vx(i)=0.0
   corr_vy(i)=0.0 
   acc_oldx(i)=0.0d0
   acc_oldy(i)=0.0d0
   vel_halfx(i)=0.0d0
   vel_halfy(i)=0.0d0
    end if
   end do

  !initialize the coordinates of 2d wavemaker(moving boundary particles) at time=0.0s
   do i=1,ntotal
    if(p(i)%id.eq.-1)then!wavemaker      
      pos_oldx(i)=p(i)%coord%x(1)
     end if
       end do

  !initialize prescribed velocity of moving boundary particles
  do i=1,ntotal
   pres_vel(i)%x(1)=0.0d0
   pres_vel(i)%x(2)=0.0d0
   end do     
 
   
   
!start time integration
write(*,*)"Starting time integration............................"
write(*,153)"predictor-corrector time integration is used."
153 format(3x,a50)
write(*,*)"****************************************************************"

  write(*,110)"Timestep" ,"Time(sec)"             
  110 format(4x,a12,2x,a12)
  write(*,111)"========","=========="
  111 format(4x,a12,2x,a12)
  time=0.0
  do t=1,nt
   !time=t*dt
     
     !write(*,120) t,time
     !120 format(1x,i12,3x,f12.6)
 

 if(time_integrator.eq.1)then                           
  !predictor-corrector time integration
 
!predictor step

!compute acceleration of fluid particles at n timestep
   if(type.eq.3)then
  call matrix(type,k_type,Xmin,Ymin,x,z,p,ntotal,h)
  end if
  call compute_press_fluid(p,ntotal,rho0,coeffs,fl_height) !compute pressure between fluid  particles
  call Compute_press_and_vel_bound(p,ntotal,h,dp,pres_vel,extra_vel,Xmin,Ymin,x,z,k_type) !compute pressure and velocity of boundary particles by extrapolation
  call Compute_density_bound(p,ntotal,rho0,coeffs,fl_height)!compute density of boundary particles
  call KGC(p,ntotal,h,dp,Xmin,Ymin,x,z,k_type) !kernal gradient correction
  call gravity(time,p,ntotal)
  call Press_Gradient(p,ntotal,h,dp,Xmin,Ymin,x,z,k_type)
  
    if(bound_type.eq.1)then
       if(type.eq.1.or.type.eq.3)then
         call  Strain_Rate_Slip(k_type,type,p,ntotal,h,dp,Xmin,Ymin,x,z)
      end if
       end if
    
       if(bound_type.gt.1)then
         if(type.eq.1.or.type.eq.3)then          
       call Strain_Rate_NoSlip(k_type,type,p,ntotal,h,dp,extra_vel,Xmin,Ymin,x,z)
         end if
            end if
     
   if(bound_type.eq.1)then
     if(type.eq.1.or.type.eq.3)then
   call viscosity_Slip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,rho0,Xmin,Ymin,x,z,fl_height)
   end if
    end if
    
   if(bound_type.gt.1)then
     if(type.eq.1.or.type.eq.3)then
     call viscosity_NoSlip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,extra_vel,rho0,Xmin,Ymin,x,z,fl_height)
     end if
      end if


   
  if(dt_option.gt.1)then 
  call  timestep(t,p,ntotal,h,fl_height,coeffs,vis,dt)
  end if
     time=(t)*dt
     
     write(*,120) t,time
     120 format(1x,i12,3x,f12.6)
  
 !compute drodt at n time level    
  call Compute_density_fluid(p,ntotal,h,dp,pres_vel,Xmin,Ymin,x,z,k_type,coeffs,fl_height,drodt)

!update density at half time step(n+1/2)
   do i=1,ntotal
     if(p(i)%id.eq.0)then ! density of fluid only
        p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
     end if
   end do


!calculate velocity at half timestep(n+1/2)
  do i=1,ntotal
       if(p(i)%id.eq.0)then !onlyfluid particles
   
    p(i)%vel%x(1)=vel_oldx(i)+(dt/2.0)*p(i)%acc%x(1)
    p(i)%vel%x(2)=vel_oldy(i)+(dt/2.0)*p(i)%acc%x(2)                      
       end if
    end do

 
!calling damping function here
 call Damp(p,ntotal,damp_i,damp_f,beta,dt,fun)
  do i=1,ntotal
    if(p(i)%id.eq.0)then
      !if(p(i)%coord%x(1).gt.dmp_ini)then
       if(p(i)%coord%x(1).gt.damp_i%x(1))then
      p(i)%vel%x(1)=p(i)%vel%x(1)*fun(i)%x(1)
      p(i)%vel%x(2)=p(i)%vel%x(2)*fun(i)%x(2)
      end if
       end if
      end do

   

!update posotion at half timestep(n+1/2)
 do i=1,ntotal
   if(p(i)%id.eq.0)then !onlyfluid particles
     p(i)%coord%x(1)=pos_oldx(i)+(dt/2.0)*vel_oldx(i)
     p(i)%coord%x(2)=pos_oldy(i)+(dt/2.0)*vel_oldy(i)
     !p(i)%coord%x(1)=pos_oldx(i)+(dt/2.0)*corr_vx(i)
     !p(i)%coord%x(2)=pos_oldy(i)+(dt/2.0)*corr_vy(i)
     
  end if
    end do





!resetting particles acc to zero
   do i=1,ntotal
     p(i)%acc%x(1)=0.0
     p(i)%acc%x(2)=0.0
    end do


!corrector step
!compute acceleration of fluid particles at half timestep(n+1/2)

  if(type.eq.3)then
  call matrix(type,k_type,Xmin,Ymin,x,z,p,ntotal,h)
  end if
  call compute_press_fluid(p,ntotal,rho0,coeffs,fl_height) !compute pressure between fluid  particles
  call Compute_press_and_vel_bound(p,ntotal,h,dp,pres_vel,extra_vel,Xmin,Ymin,x,z,k_type) !compute pressure and velocity of boundary particles by extrapolation
  call Compute_density_bound(p,ntotal,rho0,coeffs,fl_height)!compute density of boundary particles
  call KGC(p,ntotal,h,dp,Xmin,Ymin,x,z,k_type) !kernal gradient correction
  call gravity(time,p,ntotal)
  call Press_Gradient(p,ntotal,h,dp,Xmin,Ymin,x,z,k_type)
  
    if(bound_type.eq.1)then
       if(type.eq.1.or.type.eq.3)then
         call  Strain_Rate_Slip(k_type,type,p,ntotal,h,dp,Xmin,Ymin,x,z)
      end if
       end if
    
       if(bound_type.gt.1)then
         if(type.eq.1.or.type.eq.3)then          
       call Strain_Rate_NoSlip(k_type,type,p,ntotal,h,dp,extra_vel,Xmin,Ymin,x,z)
         end if
            end if
     
   if(bound_type.eq.1)then
     if(type.eq.1.or.type.eq.3)then
   call viscosity_Slip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,rho0,Xmin,Ymin,x,z,fl_height)
   end if
    end if
    
   if(bound_type.gt.1)then
     if(type.eq.1.or.type.eq.3)then
     call viscosity_NoSlip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,extra_vel,rho0,Xmin,Ymin,x,z,fl_height)
     end if
      end if



  if(dt_option.gt.1)then
  call timestep(t,p,ntotal,h,fl_height,coeffs,vis,dt) 
  end if 
  time=t*dt



!compute corrected  drodt at half timestep(n+1/2)
call Compute_density_fluid(p,ntotal,h,dp,pres_vel,Xmin,Ymin,x,z,k_type,coeffs,fl_height,drodt)

!corrrect density of fluid particles at n+1/2 time level
 do i=1,ntotal
 if(p(i)%id.eq.0)then 
   p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
 end if
  end do
 
!correct velocity at n+1/2 timestep
  do i=1,ntotal
   if(p(i)%id.eq.0)then !onlyfluid particles
   
    p(i)%vel%x(1)=vel_oldx(i)+(dt/2.0)*p(i)%acc%x(1)
    p(i)%vel%x(2)=vel_oldy(i)+(dt/2.0)*p(i)%acc%x(2)                      
       end if
    end do



  !XSph correction
 !call xsph_corr(p,ntotal,h,corr_vx,corr_vy)
 
!correct position at n+1/2 timestep
do i=1,ntotal
   if(p(i)%id.eq.0)then !onlyfluid particles
    p(i)%coord%x(1)= pos_oldx(i)+(dt/2.0)*p(i)%vel%x(1)
    p(i)%coord%x(2)= pos_oldy(i)+(dt/2.0)*p(i)%vel%x(2)

     !p(i)%coord%x(1)= pos_oldx(i)+(dt/2.0)*corr_vx(i)
     !p(i)%coord%x(2)= pos_oldy(i)+(dt/2.0)*corr_vy(i)

  end if
    end do


!advance density,position and velocity for next time step(n+1)

do i=1,ntotal
 if(p(i)%id.eq.0)then !onlyfluid particles
p(i)%dens=2.0*p(i)%dens-dens_old(i)
p(i)%vel%x(1)=2.0*p(i)%vel%x(1)-vel_oldx(i)
p(i)%vel%x(2)=2.0*p(i)%vel%x(2)-vel_oldy(i)
end if
end do

!calling damping function here
 call Damp(p,ntotal,damp_i,damp_f,beta,dt,fun)
  do i=1,ntotal
    if(p(i)%id.eq.0)then
      !if(p(i)%coord%x(1).gt.dmp_ini)then
       if(p(i)%coord%x(1).gt.damp_i%x(1))then
      p(i)%vel%x(1)=p(i)%vel%x(1)*fun(i)%x(1)
      p(i)%vel%x(2)=p(i)%vel%x(2)*fun(i)%x(2)
      end if
       end if
      end do


!update position of fluid particles
do i=1,ntotal
 if(p(i)%id.eq.0)then !onlyfluid particles
p(i)%coord%x(1)=2.0*p(i)%coord%x(1)-pos_oldx(i)
p(i)%coord%x(2)=2.0*p(i)%coord%x(2)-pos_oldy(i)
end if
end do



!compute velocity of wavemaker/piston for next timestep(n+1)
  
   do i=1,ntotal
     
     if(p(i)%id.eq.-1)then
       if(time.le.ramp_period)then
        ramp_func=0.5d0*(1.0-cos(pi*(time/ramp_period))) !       
        pres_vel(i)%x(1)=((f*S0/2.0d0)*cos(f*time) +var1*2.0*f*cos(2.0*f*time)) 
        pres_vel(i)%x(1)= pres_vel(i)%x(1)*ramp_func
        pres_vel(i)%x(2)=0.0d0
        elseif(time.gt.ramp_period)then
       pres_vel(i)%x(1)=((f*S0/2.0d0)*cos(f*time) +var1*2.0*f*cos(2.0*f*time))
       pres_vel(i)%x(2)=0.0d0
       end if
       end if
       end do  
       
    !compute position of wavemaker/piston for next timestep(n+1)
   do i=1,ntotal
     if(p(i)%id.eq.-1)then         
       p(i)%coord%x(1)=pos_oldx(i)+pres_vel(i)%x(1)*dt               
       end if
      end do

!checking if particles go out of domain if then resetting velocities of these particles to zero
       do i=1,ntotal
      if(p(i)%id.eq.0)then
        if(p(i)%coord%x(2).lt.2.0*dp)then
           p(i)%coord%x(1)=Xmin
           p(i)%coord%x(2)=Ymin          
           p(i)%acc%x(1)=0.0
           p(i)%acc%x(2)=0.0
           p(i)%id=1
           end if
        if(p(i)%coord%x(1).gt.(x-2.0*dp))then
           p(i)%coord%x(1)=Xmin
           p(i)%coord%x(2)=Ymin
           p(i)%acc%x(1)=0.0
           p(i)%acc%x(2)=0.0
           p(i)%id=1
       end if
         end if
            end do


!XSph correction
   !call xsph_corr(p,ntotal,h,corr_vx,corr_vy)


 do i=1,ntotal
  if(p(i)%id.eq.0)then
   vel_oldx(i)=p(i)%vel%x(1)
   vel_oldy(i)=p(i)%vel%x(2)
   pos_oldx(i)=p(i)%coord%x(1)
   pos_oldy(i)=p(i)%coord%x(2)
   dens_old(i)=p(i)%dens
  end if
   if(p(i)%id.eq.-1)then
    pos_oldx(i)=p(i)%coord%x(1)
    end if

 end do

!resetting particles acc to zero
   do i=1,ntotal
     p(i)%acc%x(1)=0.0
     p(i)%acc%x(2)=0.0
    end do

      end if





  !velocity verlet algorithm(still buggy:())
  
  if(time_integrator.gt.1)then !velocity verlet algorithm
    !compute velocity and position at half timestep
    do i=1,ntotal
    if(p(i)%id.eq.0)then
      vel_halfx(i)=p(i)%vel%x(1)+(dt/2.0)*acc_oldx(i)
      vel_halfy(i)=p(i)%vel%x(2)+(dt/2.0)*acc_oldy(i)
      
     
      p(i)%coord%x(1)=p(i)%coord%x(1)+(dt/2.0)*vel_halfx(i)
      p(i)%coord%x(2)=p(i)%coord%x(2)+(dt/2.0)*vel_halfy(i)
     end if
       end do

   !compute dro/dt at n timestep
    call Compute_density_fluid(p,ntotal,h,dp,pres_vel,Xmin,Ymin,x,z,k_type,coeffs,fl_height,drodt)

     !compute density of fluid  at half timestep
     do i=1,ntotal
       if(p(i)%id.eq.0)then
       p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
       end if
       end do
 
   !compute dro/dt at halftime step(n+1/2)
       call Compute_density_fluid(p,ntotal,h,dp,pres_vel,Xmin,Ymin,x,z,k_type,coeffs,fl_height,drodt)
      
       !compute density of fluid  at next timestep(n+1)
         do i=1,ntotal
       if(p(i)%id.eq.0)then
       p(i)%dens=dens_old(i)+(dt/2.0)*drodt(i)
       end if
       end do

    !compute position at next timestep(n+1)
     do i=1,ntotal
    if(p(i)%id.eq.0)then           
      p(i)%coord%x(1)=p(i)%coord%x(1)+(dt/2.0)*p(i)%vel%x(1)
      p(i)%coord%x(2)=p(i)%coord%x(2)+(dt/2.0)*p(i)%vel%x(2)
     end if
       end do

  !compute velocity at next timestep(n+1)
      
        do i=1,ntotal
         if(p(i)%id.eq.0)then
      p(i)%vel%x(1)=vel_oldx(i) +dt*acc_oldx(i)
      p(i)%vel%x(2)=vel_oldy(i) +dt*acc_oldy(i)
      
          end if
            end do




   !compute acceleration of fluid particles for next time step(n+1)
  if(type.eq.3)then
  call matrix(type,k_type,Xmin,Ymin,x,z,p,ntotal,h)
  end if
  call compute_press_fluid(p,ntotal,rho0,coeffs,fl_height) !compute pressure between fluid  particles
  call Compute_press_and_vel_bound(p,ntotal,h,dp,pres_vel,extra_vel,Xmin,Ymin,x,z,k_type) !compute pressure and velocity of boundary particles by extrapolation
  call Compute_density_bound(p,ntotal,rho0,coeffs,fl_height)!compute density of boundary particles
  call KGC(p,ntotal,h,dp,Xmin,Ymin,x,z,k_type) !kernal gradient correction
  call gravity(time,p,ntotal)
  call Press_Gradient(p,ntotal,h,dp,Xmin,Ymin,x,z,k_type)
  
    if(bound_type.eq.1)then
       if(type.eq.1.or.type.eq.3)then
         call  Strain_Rate_Slip(k_type,type,p,ntotal,h,dp,Xmin,Ymin,x,z)
      end if
       end if
    
       if(bound_type.gt.1)then
         if(type.eq.1.or.type.eq.3)then          
       call Strain_Rate_NoSlip(k_type,type,p,ntotal,h,dp,extra_vel,Xmin,Ymin,x,z)
         end if
            end if
     
   if(bound_type.eq.1)then
     if(type.eq.1.or.type.eq.3)then
   call viscosity_Slip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,rho0,Xmin,Ymin,x,z,fl_height)
   end if
    end if
    
   if(bound_type.gt.1)then
     if(type.eq.1.or.type.eq.3)then
     call viscosity_NoSlip(type,k_type,p,ntotal,nu,a,coeffs,h,dp,extra_vel,rho0,Xmin,Ymin,x,z,fl_height)
     end if
      end if


     !correct velocity at next timestep(n+1)
      
         do i=1,ntotal
       if(p(i)%id.eq.0)then
        p(i)%vel%x(1)=vel_halfx(i)+(dt/2.0)*p(i)%acc%x(1)
        p(i)%vel%x(2)=vel_halfy(i)+(dt/2.0)*p(i)%acc%x(2)                 
       end if
         end do


 !checking if particles go out of domain if then resetting velocities of these particles to zero
    do i=1,ntotal
   if(p(i)%id.eq.0)then
  if(p(i)%coord%x(2).lt.2.0*dp)then
   p(i)%vel%x(1)=0.0d0
   p(i)%vel%x(2)=0.0
   p(i)%press=0.0
    end if
  if(p(i)%coord%x(1).gt.(x-2.0*dp))then
   p(i)%vel%x(1)=0.0d0
   p(i)%vel%x(2)=0.0
   p(i)%press=0.0
   end if
   end if
  end do




    do i=1,ntotal
    if(p(i)%id.eq.0)then
      vel_oldx(i)=p(i)%vel%x(1)
      vel_oldy(i)=p(i)%vel%x(2)
      acc_oldx(i)=p(i)%acc%x(1)
      acc_oldy(i)=p(i)%acc%x(2)   
      dens_old(i)=p(i)%dens
     end if
      end do

   !resetting particles acc to zero
    do i=1,ntotal
     p(i)%acc%x(1)=0.0
     p(i)%acc%x(2)=0.0
    end do

  
        end if
      


if(mod(t,out_nt).eq.0)then
  write(*,*)"writing out final position of particles in a text file......."
open(unit=10,file='coordinate_final.txt',action='write'&
   ,status='replace',iostat=ierror)
     
    write(10,101) 'coordx','coordy','coordz'
    101 format(1x,a12,2x,a12,2x,a12)
   do i=1,ntotal
    write(10,100) p(i)%coord%x(1),p(i)%coord%x(2),p(i)%coord%x(3)
      
    100 format(1x,f12.5,2x,f12.5,2x,f12.5)
  end do
  write(*,161)"done"
 161 format(15x,a)
 write(*,162)"**********************"
 162 format(6x,a)
 end if


   end do

   



   !write(*,*)"p%acc%x(3)=",p(:)%acc%x(3)
!$$$$$$    do i=1,ntotal
!$$$$$$      if(p(i)%id.lt.0)then
!$$$$$$    write(*,*) "p%mat=",p(i)%mat(1,1)
!$$$$$$    end if
!$$$$$$    end do

   

   
     !resetting the acceleration of particles to zeroes
  
     
!$$$$$$             do i=1,ntotal
!$$$$$$            if(p(i)%id.eq.1)then
!$$$$$$             write(*,*)"p(i)%dens=",p(i)%dens
!$$$$$$             end if
!$$$$$$             end do

           

            
            !do i=1,ntotal
            !if(p(i)%id.lt.2)then
            !write(*,*)"p(i)%drodt=",drodt(i)
            !end if
            !end do
                
               !do i=1,ntotal
                 !if(p(i)%id.gt.1)then
               !write(*,*)"p(i)%vel=",p(i)%vel%x(2)
              !end if
               !end do
 


   
     
     

          !do i=1,ntotal
            !if(p(i)%id.lt.2)then
            !write(*,*)"p(i)%press=",p(i)%press
            !end if
           !end do
                  
 !write(*,*)"ref_press=",ref_press


            !do i=1,ntotal
           !if(p(i)%id.gt.1)then
            !write(*,*)"p(i)%acc%x(2)=",p(i)%acc%x(2)
            !end if
            !end do



 !write(*,*)"p%coord%x(3)=",p(:)%coord%x(3)
 !write(*,*)"p%hxx=",p(:)%hxx
  !write(*,*)"p%acc=",p(:)%acc%x(1)
  
 




 call system_clock(t2,clock_rate,clock_max)
  var1=(real(t2-t1)/real(clock_rate))
  if(var1.le.60.0d0)then
  write(*,113)var1,"sec"
  113 format("time need to finish the simulation =",1x,f8.5,a4)
  end if
  if(var1.gt.60.0d0 .and. var1.lt.3600.0d0)then
    a1=dble(var1/60.0d0)
    mint=int(var1/60.0d0)
    sec=dble((a1-mint)*60.0)
    write(*,114)"Time need to finish the simulation =",mint,"min",sec,"sec"
    114 FORMAT(1x,a,i4,a5,f6.2,a4)
    end if
   if(var1.gt.3600.0d0 .and. var1.le.86400.0d0)then
     a1=dble(var1/3600.0d0)
     hr=int(var1/3600.0d0)
     b1=dble((a1-hr)*60.0)
     mint=int((a1-hr)*60.0)
     sec=dble((b1-mint)*60.0)
     write(*,115)"Time needs to finish the simulation =",hr,"hr",mint,"min",sec,"sec"
     115 format(1x,a,i3,a4,i3,a4,f6.2,a4)
    end if
   if(var1.gt.86400.0d0)then
     a1=dble(var1/86400.0d0)
     day=int(var1/86400.0d0)
     b1=dble((a1-day)*24)
     hr=int((a1-day)*24)
     mint=int((b1-hr)*60.0)
     c1=dble((b1-hr)*60.0)
     sec=dble((c1-mint)*60.0)     
     write(*,116)"Time need to finish the simulation =",day,"day",hr,"hr",mint,"min",sec,"sec"                 
     116 format(1x,a,i2,a5,i3,a5,i3,a5,f6.2,a4)
    end if   
    
 deallocate(p)
 write(*,*)"press enter key to exit..........."
 read(*,*)


end program
