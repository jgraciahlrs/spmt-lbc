!!
!! Copyright (c) 2010        HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


module lbm_step
  use lbm_functions    
  
#include "include/replace.h"
contains
  
  subroutine lb_step(lb_dom, pnt_nxt, pnt_cur, s_par, prc, meas, tStep)	
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_nxt
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_cur
    type(sim_parameter),intent(in)  :: s_par
    type(mpl_var),intent(inout)     :: prc
    type(measure)                   :: meas  
    integer, intent(in)             :: tStep

    real(R8B),pointer,dimension(:,:,:,:)                :: swapper

#ifdef COMBINE_STREAM_COLLIDE
      !---------------------------
      ! Combined Stream-collide Routine
    call stream_collide_bgk( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#else
    !---------------------------
    ! Propagation fIn => fOut 
    call propagate(lb_dom,pnt_nxt,pnt_cur,lb_dom%state,meas)
    
    ! swap work array pointers; pnt_nxt: post-propagation, pnt_cur: post-collision
    swapper => pnt_cur; pnt_cur => pnt_nxt; pnt_nxt => swapper

    !---------------------------
    ! Collision routine fOut => fIn, excluding solid nodes 
    call collide_bgk( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#endif

    ! do boundary conditions after stream and collide steps
    call set_nrbc_antibb(lb_dom,pnt_nxt,pnt_cur,s_par,prc,meas)
    call set_bnd(lb_dom,lb_dom%state,pnt_nxt,lb_dom%rho,lb_dom%u,s_par,prc,meas,tStep)
    
  end subroutine lb_step

end module lbm_step
