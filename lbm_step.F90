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

#ifdef COMBINE_STREAM_COLLIDE
      !---------------------------
      ! Combined Stream-collide Routine
#ifdef TRT
    call stream_collide_trt(lb_dom,pnt_nxt,pnt_cur,s_par%omega)
#else
    call stream_collide_bgk( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#endif
#else
    !---------------------------
    ! Propagation fIn => fOut 
    call propagate(lb_dom,pnt_nxt,pnt_cur,lb_dom%state,meas)
    
    ! pnt_nxt : fOut: post-propagation, pnt_cur: fIn: post-collision
    
    pnt_nxt => lb_dom%fIn
    pnt_cur => lb_dom%fOut

    !---------------------------
    ! Collision routine fOut => fIn, excluding solid nodes 
#ifdef TRT
    call collide_trt(lb_dom,pnt_nxt,pnt_cur,s_par%omega)
#else
#ifdef MRT
    call collide_mrt( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#else
    call collide_bgk( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#endif
#endif
#endif

    ! do boundary conditions after stream and collide steps
    call set_nrbc_antibb(lb_dom,pnt_nxt,pnt_cur,s_par,prc,meas)
    call set_bnd(lb_dom,lb_dom%state,pnt_nxt,lb_dom%rho,lb_dom%u,s_par,prc,meas,tStep)
    
  end subroutine lb_step

end module lbm_step
