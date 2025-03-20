module wavelength
contains
subroutine ComputeWavelength(T,d,L,DepthType)
implicit none
double precision,intent(in)::T
double precision,intent(in)::d
double precision,intent(out)::L
character(13),intent(out)::DepthType
double precision,parameter::pi=3.1416
double precision,parameter::g=9.81
double precision::xr,dL,length,Lhs,Rhs
double precision::l1=10.0d0! limit wave length for trials
integer::i
dL=1d-5 !step size for trials


 xr=(g*T*T)/(2.0*pi)!deep water condition

  if((d/xr).ge.0.50d0)then
   DepthType='Deep'
   elseif((d/xr).le.0.050d0)then
   DepthType='shallow'
   elseif((d/xr).gt.0.050 .and. (d/xr).lt.0.50d0)then
   DepthType='intermediate'
   end if
   
 Length=0.0d0
 Typ:select case(DepthType)
 case('Deep')
  L=xr
 case('shallow')
  L=sqrt(g*d*T)
 case('intermediate')
   length=length+dL
   Lhs=length/((g*T*T)/(2.0*pi))
   Rhs=tanh((2.0*pi*d)/length)
   do while((Lhs-Rhs).lt.1d-6)
    length=length+dL
    Lhs=length/((g*T*T)/(2.0*pi))
    Rhs=tanh((2.0*pi*d)/length)              
      end do
     L=length  
  end select Typ

  end subroutine
  end module
    