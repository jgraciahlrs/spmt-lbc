!!
!! Copyright (c) 2010        HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


module lb_bc
   use mpl_lib
   use lbmodel
   use function_set
   use timing
   use lb_init
   use lb_geo 
#include "include/replace.h"
contains





   subroutine set_periodic(lb_dom,fIn,s_par,prc,meas)
   !------------------------------------------------------------------------
   !
   ! Tested and ok
   !
      implicit none
      type(lb_block)       :: lb_dom
      type(sim_parameter)  :: s_par 
      type(mpl_var)        :: prc
      type(measure)        :: meas
      integer              :: ii,ll
      real(R8B),dimension(NDX(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
      real(R8B)    :: diff


      if(s_par%periodic(1).eqv..true.) then
      if(prc%size > 1) then; write(*,*) "Periodic boundaries not yet implemented for np>1"; stop; endif 
      !write(*,*) "x periodic"
      fIn(NDX(:,1,:,:)) = fIn(NDX(:,lb_dom%lx(1)-1,:,:)) 
      fIn(NDX(:,lb_dom%lx(1),:,:)) = fIn(NDX(:,2,:,:)) 
      endif      

      if(s_par%periodic(2).eqv..true.) then
      if(prc%size > 1) then; write(*,*) "Periodic boundaries not yet implemented for np>1"; stop; endif 
      !write(*,*) "y periodic"
      fIn(NDX(:,:,1,:)) = fIn(NDX(:,:,lb_dom%lx(2)-1,:)) 
      fIn(NDX(:,:,lb_dom%lx(2),:)) = fIn(NDX(:,:,2,:)) 
      endif     

      if(s_par%periodic(3).eqv..true.) then
      if(prc%size > 1) then; write(*,*) "Periodic boundaries not yet implemented for np>1"; stop; endif 
      !write(*,*) "z periodic"
      fIn(NDX(:,:,:,1)) = fIn(NDX(:,:,:,lb_dom%lx(3)-1)) 
      fIn(NDX(:,:,:,lb_dom%lx(3))) = fIn(NDX(:,:,:,2)) 
      endif     

   end subroutine set_periodic
   !------------------------------------------------------------------------






subroutine set_bnd(lb_dom,state,fIn,rho,u,s_par,prc,meas,tStep)

!------------------------------------------------------------------------
!
! set the boundary conditions in each timestep.
!


   implicit none
   type(lb_block)       :: lb_dom
   type(sim_parameter)  :: s_par 
   type(mpl_var)        :: prc
   type(measure)        :: meas
   integer              :: tStep
   integer              :: state(lb_dom%lx(1),lb_dom%lx(2),lb_dom%lx(3))
   integer              :: i,j,k,l,x_start,x_end,y_start,y_end,z_start,z_end
   integer              :: incr_x,incr_y,incr_z,lx,ly,lz
   real(R8B),dimension(NDX(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
   real(R8B),dimension(NDX(1:NDIM,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: u
   real(R8B),dimension(0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1)    :: rho

   intent(inout) :: fIn,rho,u


   call cpu_time_measure(meas%tSt_comm) 


   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)


   ! Set the different starting coordinates:
   ! If crd == 1     ->  first cell is boundary!   -> start bc on 2 
   ! If crd == np-1  ->  last  cell is boundary!   -> start bc on lx-1 

   if(prc%crd(1) == 0)          then; x_start = 2;    else; x_start=1;  end if
   if(prc%crd(1) == prc%np(1)-1) then; x_end   = lx-1; else; x_end  =lx; end if
   if(prc%crd(2) == 0)          then; y_start = 2;    else; y_start=1;  end if
   if(prc%crd(2) == prc%np(2)-1) then; y_end   = ly-1; else; y_end  =ly; end if
   if(prc%crd(3) == 0)          then; z_start = 2;    else; z_start=1;  end if
   if(prc%crd(3) == prc%np(3)-1) then; z_end   = lz-1; else; z_end  =lz; end if

   ! Set the delta increase of coordinates for process coordinates > 1 
   incr_x = prc%crd(1)*lx
   incr_y = prc%crd(2)*ly
   incr_z = prc%crd(3)*lz

!--------------------------------------------------------------------
! Treat non-reflecting and periodic boundaries 


 
! Treat inlet / outlet 
 

!--------------------------------------------------------------------
! Select defined problem  

      select case(s_par%problem)


!----------------------------------------
! Lid-driven cavity

     case(cavity)
        ! Reset boundary 
        ! macroscopic Variables
        j = ly
        if(prc%crd(2) == prc%np(2)-1) then !wenn Block in oberster Reihe
         do k=z_start,z_end
            do i=x_start,x_end
              ! top lid
               rho(i,j,k)   = 1._R8B
               u(NDX(1,i,j,k))   = s_par%umax/sqrt(2.d0)
               u(NDX(2,i,j,k))   = 0._R8B
               u(NDX(3,i,j,k))   = s_par%umax/sqrt(2.d0)
            end do
         end do
         end if

         if(prc%crd(2) == prc%np(2)-1) then
            call calc_fEq(lb_dom,rho,u,fIn,x_start,x_end,ly,ly,z_start,z_end,s_par) !top lid
         end if

   end select


   call cpu_time_measure(meas%tEnd_comm)

   meas%sbnd_duration = meas%tEnd_comm - meas%tSt_comm + meas%sbnd_duration



end  subroutine set_bnd
!------------------------------------------------------------------------








end module lb_bc        


