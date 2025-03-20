module Domain

contains

subroutine create_domain(p,np,Xmin,length,length_f,height,wavemaker_posX,wmaker_height,height_f,damp_size,dp,dp1)
use particles
implicit none
type(particle),intent(out)::p(:)
integer,intent(out)::np
double precision,intent(in)::length,height,dp,length_f,height_f,dp1,Xmin
double precision,intent(in)::wmaker_height,damp_size,wavemaker_posX
!local variable
integer::i,j,nn
integer::nx,ny,nz,nx_f,nz_f,nz_wavemaker
double precision::fl_total_length
double precision::trans

fl_total_length=length_f+damp_size !total length of fluid domain

nx=floor(length/dp)
nz=floor(height/dp1)
nz_wavemaker=floor(wmaker_height/dp1)
nx_f=nint(fl_total_length/dp)
nz_f=ceiling(height_f/dp)


p(:)%coord%x(1)=0.0!initialize all attributes
p(:)%coord%x(2)=0.0
nn=0
!3 layers of dummy fluid particles

!creating wavemaker particles

 !creating 1st layer
do i=1,nz_wavemaker
  nn=nn+1
   p(nn)%coord%x(1)=wavemaker_posX
   p(nn)%coord%x(2)=(i)*dp1+(2.0*dp1)
   p(nn)%id=-1
 end do

 !creating 2nd layer

 do i=1,nz_wavemaker
    nn=nn+1
     p(nn)%coord%x(1)=trans(p(nn-nz_wavemaker)%coord%x(1),dp)
     p(nn)%coord%x(2)=(i)*(dp1)+(2.0*dp1)
     p(nn)%id=-1
 end do
 
   ! creating 3rd layer
    do i=1,nz_wavemaker
    nn=nn+1
     p(nn)%coord%x(1)=trans(p(nn-nz_wavemaker)%coord%x(1),dp)
     p(nn)%coord%x(2)=i*dp1+(2.0*dp1)
     p(nn)%id=-1
 end do

  !creating 4th layer
  do i=1,nz_wavemaker
    nn=nn+1
     p(nn)%coord%x(1)=trans(p(nn-nz_wavemaker)%coord%x(1),dp)
     p(nn)%coord%x(2)=i*dp1+(2.0*dp1)
     p(nn)%id=-1
 end do



 !creation of right boundary

 
 !creating 1st layer
   do i=1,nz
    nn=nn+1
   p(nn)%coord%x(1)=wavemaker_posX+fl_total_length
   p(nn)%coord%x(2)=(i-1)*dp1
   p(nn)%id=-2
 end do
  
  

    
  !creation of second layer
    do i=1,nz-1
     nn=nn+1
   p(nn)%coord%x(1)=trans(p(nn-nz+1)%coord%x(1),-dp)
   p(nn)%coord%x(2)=p(nn-nz+1)%coord%x(2)
   p(nn)%id=-2
   end do
 
    !creation of 3rd layer
    do i=1,nz-2
     nn=nn+1
   p(nn)%coord%x(1)=trans(p(nn-nz+2)%coord%x(1),-dp)
   p(nn)%coord%x(2)=p(nn-nz+2)%coord%x(2)
   p(nn)%id=-2
   end do

  !create bottom boundary
  !creating 1st layer
    
     
      do j=1,nx-1
       nn=nn+1
      p(nn)%coord%x(1)=Xmin+(j)*dp
      p(nn)%coord%x(2)=0.0d0
      p(nn)%id=-2
    end do
     
!creating 2nd layer
   
    do j=1,nx-2
       nn=nn+1
    p(nn)%coord%x(1)=trans(p(nn-nx+1)%coord%x(1),0.0d0)
    p(nn)%coord%x(2)=trans(p(nn-nx+1)%coord%x(2),dp1)
    p(nn)%id=-2
      end do
     


 !creation of 3rd layer
      

    do j=1,nx-3
       nn=nn+1
    p(nn)%coord%x(1)=trans(p(nn-nx+2)%coord%x(1),0.0d0)
    p(nn)%coord%x(2)=trans(p(nn-nx+2)%coord%x(2),dp1)
    p(nn)%id=-2
      end do

!creation of fluid particles
 do i=1,nz_f
  do j=1,nx_f-5
  nn=nn+1
  p(nn)%coord%x(1)=wavemaker_posX+((j-1)*dp)+(3.0*dp)
  p(nn)%coord%x(2)=((i-1)*dp1)+(3.0*dp1)
  p(nn)%id=0
  end do
    end do
   
   np=nn

 end subroutine
 

  end module  
  

 function trans(x,dx)
 double precision,intent(in)::x
 double precision,intent(in)::dx
 double precision::trans
 trans=x+dx
 end function

     


