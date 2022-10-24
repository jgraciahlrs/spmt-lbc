!!
!! Copyright (c) 2010-2011   HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


! ------------------------------------------------------------------------
!
!         oooo   .o8
!          888   888
!          888   888oooo.   .ooooo.
!          888   d88   88b d88    Y8
!          888   888   888 888
!          888   888   888 888   .o8
!         o888o  `Y8bod8P  `Y8bod8P
!
! This is the Lattice Boltzmann Kernel LBC for testing purposes only.
! Loosely based on J. Bernsdorfs anb-Lattice Boltzmann-Code
! 2D Lattice Boltzmann calculations with boundary and obstacle implementation.
!
!------------------------------------------------------------------------
!
!
!========================================================================

PROGRAM lbmain
   use mpl_lib
   use lbmodel
   use lbm_functions
   use lb_bc
   use mpl_set
   use function_set
   use tools
   use timing
#ifdef USE_TILES
   use lbm_step_tiled
#else
   use lbm_step
#endif /* USE_TILES */

#include "include/replace.h"
   implicit none
   type(measure)        :: meas
   type(lb_block)       :: lb_dom
   type(sim_parameter)  :: s_par
   type(mpl_var)        :: prc
   real(R8B)            :: total_start,total_end
   integer              :: tStep,countOut,result_state
   integer              :: i,j,k
   integer,pointer      :: tstep_cnt
   logical              :: stop_request
   real(R8B),pointer,dimension(:,:,:,:) :: swapper

   countOut             = 0
   result_state         = 0
   tstep_cnt => gtstep_cur

   stop_request         = .false.

!------------------------------------------------------------------------
! Initialize environment

   call mpl_init(meas,prc)
   call cpu_time_measure(total_start)
   call cpu_time_measure(meas%tStart)
   call read_params(s_par,'lbc.params',prc)

   call mpl_init_domain(lb_dom,s_par,prc)

   call alloc_mem(lb_dom,s_par,prc)


!------------------------------------------------------------------------
! Initialize data and distribute to processes
!
! does this fix the MPI serialization?
if (mod(prc%rk,2)>0) then
   prc%nghb(1:2) = prc%nghb(2:1:-1)
endif

#ifdef USE_TILES
   call init_tile(lb_dom, s_par, prc)
!$OMP TASKWAIT ON(lb_dom)
#endif /* USE_TILES */

#ifdef USE_OMP_FINE
!$omp parallel default(none) shared(lb_dom)
!$omp workshare
#endif /* USE_OMP_FINE */
   lb_dom%fIn(NDX(1,:,:,:))=0.33333333333333d0
   lb_dom%fIn(NDX(2:7,:,:,:))=0.0555555555555d0
   lb_dom%fIn(NDX(8:19,:,:,:))=0.02777777777777d0
!   lb_dom%fIn  = -10.d0
   lb_dom%fOut = lb_dom%fIn
#ifdef USE_OMP_FINE
!$omp end workshare nowait
!$omp end parallel
#endif /* USE_OMP_FINE */


   !---------------------------
   ! Master thread initializes initial field and geometry

#ifdef INIT_WITH_ROOT
   if(prc%rk==0) then
      call get_geo(lb_dom,s_par,prc)
      call init_field(lb_dom%gfIn,lb_dom%gu,lb_dom%grho,lb_dom%gstate,s_par,lb_dom,prc)
      call calc_fEq_global(lb_dom%grho,lb_dom%gu,lb_dom%gfIn,s_par) ! reset complete domain to eq
   end if

   !---------------------------
   ! Decompose domain
   call mpl_decompose_domain(lb_dom,s_par,prc)

#else
   call get_geo_each(lb_dom,s_par,prc)
   call init_field_each(lb_dom,s_par,prc)
   call calc_fEq(lb_dom,lb_dom%rho,lb_dom%u,lb_dom%fIn,0,0,0,0,0,0,s_par)
#endif
   !---------------------------
   ! Initialize the obstacle and boundary vectors
   call init_vector(lb_dom,s_par,prc)

   call ramp_u(s_par,0)

   call set_bnd(lb_dom,lb_dom%state,lb_dom%fIn,lb_dom%rho,lb_dom%u,s_par,prc,meas,0)

! for timing
!$OMP TASKWAIT
   call cpu_time_measure(meas%tEnd)
   meas%init_duration = real(meas%tEnd,8) - real(meas%tStart,8)
!------------------------------------------------------------------------
! Main Loop
!
   if(prc%rk==0) then
     write(*,*); write(*,*) "  --------------------------------------"
     write(*,*) "  Starting main time loop"; write(*,*)
   endif

#ifdef USE_SMPSS
   ! time main loop
   call cpu_time_measure(meas%tStart)
#endif /* USE_SMPSS */

   tstep_cnt = 0
   lb_dom%pnt_nxt => lb_dom%fIn; lb_dom%pnt_cur => lb_dom%fOut    ! init work arrays
   do tStep = 1,s_par%tMax

#ifndef USE_SMPSS
      call cpu_time_measure(meas%tStart)
#endif /* Not USE_SMPSS */
      !---------------------------
      ! Set Boundaries on fIn
      call treat_initial(lb_dom,s_par,tstep,prc)

      tstep_cnt = tstep_cnt+1

      call ramp_u(s_par,tstep)
      ! swap work array pointers
#ifdef USE_SMPSS
      !$OMP TASK INOUT(lb_dom) label(swap_lb_dom_pnts_smpss) PRIORITY(50) no_copy_deps
      call swap_lb_dom_pnts_task(lb_dom)
      !$OMP END TASK
#else
      swapper => lb_dom%pnt_cur; lb_dom%pnt_cur => lb_dom%pnt_nxt; lb_dom%pnt_nxt => swapper
#endif /* USE_SMPSS */

      !---------------------------
      ! Exchange between domains
      if(prc%size > 1) then
#ifdef USE_SMPSS
        !$OMP TASK INOUT(lb_dom) label(mpl_exchange) PRIORITY(50)  no_copy_deps
        call mpl_exchange(lb_dom,lb_dom%pnt_cur,prc,meas)
        !$OMP END TASK
#else
        call mpl_exchange(lb_dom,lb_dom%pnt_cur,prc,meas)
#endif
       end if


      ! do a single Lattice Boltzmann step
      call lb_step(lb_dom, lb_dom%pnt_nxt, lb_dom%pnt_cur, s_par, prc, meas, tStep)

#ifndef USE_SMPSS
      call cpu_time_measure(meas%tEnd)
      meas%duration = meas%duration + real(meas%tEnd,8) - real(meas%tStart,8)
#endif /* Not USE_SMPSS */

      !---------------------------
      ! Output and error check
      if(MOD(tStep,s_par%tOut)==0 .and. tStep >= s_par%tOutBeg) then
#ifdef USE_SMPSS
         !$OMP TASKWAIT !!ON(lb_dom)
         ! time main loop at barrier, excluding output
         call cpu_time_measure(meas%tEnd)
         meas%duration = meas%duration + real(meas%tEnd,8) - real(meas%tStart,8)
#endif /* USE_SMPSS */
         ! calc macroscopic values (is also done in do_kinetic)
         call calc_macr_vals(lb_dom,lb_dom%u,lb_dom%rho,lb_dom%pnt_nxt,meas)
         countOut=countOut+1
         call calc_rho_ges(lb_dom,lb_dom%pnt_nxt,s_par%rhoges,tStep,prc)
         if(tStep > s_par%tMax-s_par%tOut  ) s_par%goend=.true.
         if(result_state == -1) s_par%goend=.true.
         if(stop_request .eqv. .true.) s_par%goend=.true.

         inquire(FILE='stop',exist=stop_request)
         if(s_par%save_output) then
            call gather_output(lb_dom,countOut,tStep,s_par,prc)
         endif

         ! Check for errors. If errors true, then exit after next iteration
         if(result_state == 0) then
            if(s_par%initial .eqv. .false.) &
                 call check_density(lb_dom,tStep,result_state,prc%rk)
            call comm_res(result_state,prc)
         elseif(result_state ==-1 .and. s_par%goend .eqv. .true.) then
            write(*,*) "Result is errorous. Exiting..."
            exit
         endif

         if(stop_request .and. s_par%goend ) then
            write(*,*) "Stop has been requested by stopfile 'Stop'"
            exit
         endif
#ifdef USE_SMPSS
         call cpu_time_measure(meas%tStart)
#endif /* USE_SMPSS */

      endif

   end do

#ifdef USE_SMPSS
!$OMP TASKWAIT
! finish timing main loop
call cpu_time_measure(meas%tEnd)
meas%duration = meas%duration + real(meas%tEnd,8) - real(meas%tStart,8)
#endif /* USE_SMPSS */

!------------------------------------------------------------------------
! Done with Main loop

   if(prc%rk==0) then
     write(*,*); write(*,*) "  --------------------------------------"
     write(*,*) "  "; write(*,*)
   endif
!------------------------------------------------------------------------
! Finish up
!
   call cpu_time_measure(total_end)
   meas%total_duration = total_end - total_start

   call write_status(meas,result_state,tstep_cnt,s_par,prc)

   call dealloc_mem(lb_dom,s_par)

   call mpl_finish(prc)





END PROGRAM lbmain
!========================================================================
