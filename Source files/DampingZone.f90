module Damping_Zone
contains
subroutine Damp(p,ntotal,x0,xf,beta,dt,func)
use particles
implicit none
type(particle),intent(in)::p(:)!particles data
type(vector_2d),intent(in)::x0,xf
double precision,intent(in)::beta!reduction function
double precision,intent(in)::dt!time step size
type(vector_2d),intent(out)::func(:)
integer,intent(in)::ntotal
!local variable
integer::i
double precision::var1,var2
 do i=1,ntotal
 func(i)%x(1)=0.0d0
 func(i)%x(2)=0.0d0
 end do
 
  do i=1,ntotal
    if(p(i)%id.eq.0)then              
        var1=(p(i)%coord%x(1)-x0%x(1))/(xf%x(1)-x0%x(1))
        var2=(p(i)%coord%x(2)-x0%x(2))/(xf%x(2)-x0%x(2))        
        func(i)%x(1)=1.0-dt*beta*var1*var1
        func(i)%x(2)=1.0-dt*beta*var2*var2                         
         end if
          end do

   end subroutine
    end module
        
