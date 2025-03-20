module Compute_press_grad

contains
subroutine Press_Gradient(p,ntotal,h,dp,Xmin,Ymin,length,height,kernal_type)
use particles
use kernal
!use ptcl_matrix
implicit none
type(particle),intent(inout)::p(:)
double precision,intent(in)::h
integer,intent(in)::ntotal
integer,intent(in)::kernal_type
double precision,intent(in)::Xmin,Ymin  
double precision,intent(in)::length,height!lenght and height of confined wall or simulation domain
double precision,intent(in)::dp
double precision::wij,dwdx,dwdy
double precision,allocatable::dwdx_corr(:),dwdy_corr(:)
integer::i,j,k,sc_index
double precision::rij_x,rij_y,rij
integer::cellx,celly,cell_total
double precision::xmax,ymax,ln_x,ln_y!max domain size
integer,allocatable::box(:),np(:),ll(:)
integer::ix,iy,neigh_cell,i1,i2
integer::ncellx(4),ncelly(4) !no.of neighbour cells
double precision,allocatable::sum1(:),sum2(:),sum3(:),sum4(:)
double precision::s1,s2,s3,s4
double precision::var1,var2,var3,var4,var5
double precision::sumv1,sumv2
integer::j1,j2,jj,nc
double precision::kf
double precision,allocatable::p_force_x(:),p_force_y(:)
double precision,allocatable::press_sum_x(:),press_sum_y(:)






if(kernal_type.eq.1)then
  kf=3.0d0
else
 kf=2.0d0
end if


ln_x=-(abs(Xmin)+kf*h) !minimum domain extension in both x and y direction
ln_y=-(abs(Ymin)+kf*h)

xmax=length+kf*h !assuming square cell (rectangular cell is also possible,in this case(xmax.ne.ymax))
ymax=height+kf*h


cellx=nint((xmax)/(kf*h)) !no. of cell in x direction
celly=nint((ymax)/(kf*h)) !no. of cell in y direction
cell_total=cellx*celly  !total no. of cell

  allocate(box(cell_total),np(cell_total),ll(ntotal),p_force_x(ntotal),p_force_y(ntotal))
  allocate(sum1(ntotal),sum2(ntotal),sum3(ntotal),sum4(ntotal))
  allocate(dwdx_corr(ntotal),dwdy_corr(ntotal))
  allocate(press_sum_x(ntotal),press_sum_y(ntotal))
  
    do i=1,ntotal
    dwdx_corr(i)=0.0d0
    dwdy_corr(i)=0.0d0
    end do
 


  box=0

  do i=1,cell_total
   np(i)=0
  end do
  
  do i=1,ntotal
    ll(i)=0
   end do  

 do i=1,ntotal
  p_force_x(i)=0.0
  p_force_y(i)=0.0
  sum1(i)=0.0
  sum2(i)=0.0
  sum3(i)=0.0
  sum4(i)=0.0
  press_sum_x(i)=0.0
  press_sum_y(i)=0.0
 end do

!call thee matrix here
  !call matrix(p,ntotal,h) 


  
 !sorting particles into neach cell
    
     !sorting particles into cells
     do k=1,ntotal      
      ix=int((p(k)%coord%x(1)-Xmin)/(kf*h))+1
      !ix=(int((p(k)%coord%x(1))/(kf*h)))+1
      iy=int(abs(p(k)%coord%x(2))/(kf*h))+1              
      sc_index=ix+(iy-1)*cellx
      ll(k)=box(sc_index)
      box(sc_index)=k
      np(sc_index)=np(sc_index)+1 !no.of particles in each cell      
      end do


!$$$$$$      do k=1,ntotal      
!$$$$$$       icell=int(p(k)%coord%x(1)/(3.0*h))+1
!$$$$$$       kcell=int(p(k)%coord%x(2)/(3.0*h))+1              
!$$$$$$       ii=icell+(kcell-1)*ncx            
!$$$$$$       np(ii)=np(ii)+1 !no.of particles in each cell 
!$$$$$$       ibox(ii,np(ii))=k   
!$$$$$$       end do



   
 
       do i2=1,celly
          do i1=1,cellx
         sc_index=i1+(i2-1)*cellx
       if(np(sc_index).gt.0)then!there are particles in the cell
         i=box(sc_index)
                                  
         do while(i.gt.0)
        
           j=ll(i)
           do while(j.gt.0) !loop in same cell
           if(p(i)%id*p(j)%id.le.0)then
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)
          
          !kernal gradient correction
!$$$$$$        dwdx_corr(i)=p(i)%mat(1,1)*dwdx+p(i)%mat(1,2)*dwdy
!$$$$$$        dwdy_corr(i)=p(i)%mat(2,1)*dwdx+p(i)%mat(2,2)*dwdy
!$$$$$$        dwdx_corr(j)=p(j)%mat(1,1)*dwdx+p(j)%mat(1,2)*dwdy
!$$$$$$        dwdy_corr(j)=p(j)%mat(2,1)*dwdx+p(j)%mat(2,2)*dwdy
!$$$$$$         
!$$$$$$    
!$$$$$$        !Computation of pressure gradients
!$$$$$$       sum1(i)=sum1(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx_corr(i))
!$$$$$$       sum2(i)=sum2(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy_corr(i))
!$$$$$$       sum1(j)=sum1(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx_corr(j))
!$$$$$$       sum2(j)=sum2(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy_corr(j))






      sum1(i)=sum1(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx)
      sum2(i)=sum2(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy)
      sum1(j)=sum1(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx)
      sum2(j)=sum2(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy)
      
 
!$$$$$$       sum3(i)=sum3(i)+(dwdx*(-rij_x)*(p(j)%mass/p(j)%dens))
!$$$$$$       sum4(i)=sum4(i)+(dwdy*(-rij_y)*(p(j)%mass/p(j)%dens))
!$$$$$$       sum3(j)=sum3(j)+(dwdx*(-rij_x)*(p(i)%mass/p(i)%dens))
!$$$$$$       sum4(j)=sum4(j)+(dwdy*(-rij_y)*(p(i)%mass/p(i)%dens))


    !computation of pressure gradients by KGF-SPH method
!$$$$$$       sum1(i)=sum1(i)+(p(j)%press*wij*p(j)%mass/p(j)%dens)
!$$$$$$       sum1(j)=sum1(j)+(p(i)%press*wij*p(i)%mass/p(i)%dens)
!$$$$$$       sum2(i)=sum2(i)+(p(j)%press*(-rij_x)*wij*p(j)%mass/p(j)%dens)
!$$$$$$       sum2(j)=sum2(j)+(p(i)%press*(rij_x)*wij*p(i)%mass/p(i)%dens)
!$$$$$$       sum3(i)=sum3(i)+(p(j)%press*(-rij_y)*wij*p(j)%mass/p(j)%dens)
!$$$$$$       sum3(j)=sum3(j)+(p(i)%press*(rij_y)*wij*p(i)%mass/p(i)%dens)

      

     
                                

             

               
            end if
             end if
           j=ll(j)
         end do

       !now for neighbour cells

          ncellx=(/1,0,1,-1/)
          ncelly=(/0,1,1, 1/)

        do k=1,4
        ix=i1+ncellx(k)
        iy=i2+ncelly(k)

        !PERIOD bc

       

        if(ix.lt.1) cycle
        
        if(ix.gt.cellx) cycle

        if(iy.lt.1) cycle
      
        if(iy.gt.celly) cycle

       neigh_cell=ix+(iy-1)*cellx
 
        if(np(neigh_cell).gt.0)then
           j=box(neigh_cell)
         do while(j.gt.0)
         if(p(i)%id*p(j)%id.le.0)then   
       rij_x=p(i)%coord%x(1)-p(j)%coord%x(1)
       rij_y=p(i)%coord%x(2)-p(j)%coord%x(2)
       rij=rij_x**2+rij_y**2

        if(rij.le.kf*kf*h*h)then
          call ComKer(kernal_type,rij,rij_x,rij_y,h,wij,dwdx,dwdy)

         !kernal gradient correction
!$$$$$$        dwdx_corr(i)=p(i)%mat(1,1)*dwdx+p(i)%mat(1,2)*dwdy
!$$$$$$        dwdy_corr(i)=p(i)%mat(2,1)*dwdx+p(i)%mat(2,2)*dwdy
!$$$$$$        dwdx_corr(j)=p(j)%mat(1,1)*dwdx+p(j)%mat(1,2)*dwdy
!$$$$$$        dwdy_corr(j)=p(j)%mat(2,1)*dwdx+p(j)%mat(2,2)*dwdy
!$$$$$$ 
!$$$$$$       !Computation of pressure gradients
!$$$$$$       sum1(i)=sum1(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx_corr(i))
!$$$$$$       sum2(i)=sum2(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy_corr(i))
!$$$$$$       sum1(j)=sum1(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx_corr(j))
!$$$$$$       sum2(j)=sum2(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy_corr(j))









      sum1(i)=sum1(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx)
      sum2(i)=sum2(i)+(p(j)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy)
      sum1(j)=sum1(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdx)
      sum2(j)=sum2(j)-(p(i)%mass*((p(i)%press/p(i)%dens**2)+(p(j)%press/p(j)%dens**2))*dwdy) 

!$$$$$$       sum3(i)=sum3(i)+(dwdx*(-rij_x)*(p(j)%mass/p(j)%dens))
!$$$$$$       sum4(i)=sum4(i)+(dwdy*(-rij_y)*(p(j)%mass/p(j)%dens))
!$$$$$$       sum3(j)=sum3(j)+(dwdx*(-rij_x)*(p(i)%mass/p(i)%dens))
!$$$$$$       sum4(j)=sum4(j)+(dwdy*(-rij_y)*(p(i)%mass/p(i)%dens))
     

     !computation of presure gradient by  KGF-SPH method

!$$$$$$       sum1(i)=sum1(i)+(p(j)%press*wij*p(j)%mass/p(j)%dens)
!$$$$$$       sum1(j)=sum1(j)+(p(i)%press*wij*p(i)%mass/p(i)%dens)
!$$$$$$       sum2(i)=sum2(i)+(p(j)%press*(-rij_x)*wij*p(j)%mass/p(j)%dens)
!$$$$$$       sum2(j)=sum2(j)+(p(i)%press*(rij_x)*wij*p(i)%mass/p(i)%dens)
!$$$$$$       sum3(i)=sum3(i)+(p(j)%press*(-rij_y)*wij*p(j)%mass/p(j)%dens)
!$$$$$$       sum3(j)=sum3(j)+(p(i)%press*(rij_y)*wij*p(i)%mass/p(i)%dens)



     
   

     
          end if
           end if

       j=ll(j)
   
     end do

       end if

         end do              

         i=ll(i)
          end do 
         
              end if
               end do
                     end do


     do i=1,ntotal
        if(p(i)%id.eq.0)then
      p(i)%acc%x(1)=p(i)%acc%x(1)-sum1(i)!/sum3(i)
      p(i)%acc%x(2)=p(i)%acc%x(2)-sum2(i)!/sum4(i)
       else
         p(i)%acc%x(1)=0.0
         p(i)%acc%x(2)=0.0
        end if
     end do

!$$$$$$    do i=1,ntotal
!$$$$$$      if(p(i)%id.eq.1)then
!$$$$$$      !acceleration due to pressure gradient(KGF-SPH method)
!$$$$$$        p(i)%acc%x(1)=p(i)%acc%x(1)-(1.0/p(i)%dens)*(p(i)%mat(2,1)*sum1(i)+p(i)%mat(2,2)*sum2(i)+p(i)%mat(2,3)*sum3(i))
!$$$$$$        p(i)%acc%x(2)=p(i)%acc%x(2)-(1.0/p(i)%dens)*(p(i)%mat(3,1)*sum1(i)+p(i)%mat(3,2)*sum2(i)+p(i)%mat(3,3)*sum3(i))
!$$$$$$     end if
!$$$$$$      end do


      

      end subroutine
        end module
