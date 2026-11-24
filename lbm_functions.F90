!!
!! Copyright (c) 2010-2012   HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


module lbm_functions
#ifdef USE_OMP_FINE
   use omp_lib
#endif /* USE_OMP_FINE */
#ifdef USE_OMP_TARGET_FINE
   use omp_lib
#endif
   use mpl_lib
   use lbmodel
   use function_set
   use timing
   use lb_init
   use lb_bc
   use lb_geo 
#include "include/replace.h"
contains



#ifndef COMBINE_STREAM_COLLIDE
   subroutine propagate(lb_dom,fOut,fIn,state,meas)

!------------------------------------------------------------------------
! propagation step: move the densities to the neighbor nodes
! performs PULL of densities


   implicit none
   type(lb_block),intent(in)    :: lb_dom
   type(measure),intent(inout)  :: meas  
   real(R8B),pointer,dimension(:,:,:,:),intent(in)        :: fIn
   real(R8B),pointer,dimension(:,:,:,:),intent(out)       :: fOut
   integer,pointer,dimension(:,:,:),intent(in)            :: state
   integer                      :: i,j,k,l,count_err
   real(R8B)                    :: ftmp(nnod)

      call cpu_time_measure(meas%tSt_comm) 
      count_err=0

#ifdef USE_OMP_TARGET_FINE
      !$omp  target teams distribute parallel do collapse(2) &
#endif
#ifdef USE_OMP_FINE
      !$omp  parallel do collapse(2) &
#endif
#if defined(USE_OMP_TARGET_FINE) || defined (USE_OMP_FINE)
      !$omp& default(none) private(ftmp,i,j,k) &
      !$omp& shared(fIn,fOut,opp,lb_dom,state,cx)
#endif /* USE_OMP_TARGET_FINE || USE_OMP_FINE */
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
                  fOut(NDX(I__0,i,j,k))=fIn(NDX(I__0,i-cx(I__0,1),j-cx(I__0,2),k-cx(I__0,3)))
                  fOut(NDX(I__E,i,j,k))=fIn(NDX(I__E,i-cx(I__E,1),j-cx(I__E,2),k-cx(I__E,3)))
                  fOut(NDX(I__W,i,j,k))=fIn(NDX(I__W,i-cx(I__W,1),j-cx(I__W,2),k-cx(I__W,3)))
                  fOut(NDX(I__N,i,j,k))=fIn(NDX(I__N,i-cx(I__N,1),j-cx(I__N,2),k-cx(I__N,3)))
                  fOut(NDX(I__S,i,j,k))=fIn(NDX(I__S,i-cx(I__S,1),j-cx(I__S,2),k-cx(I__S,3)))
                  fOut(NDX(I_NE,i,j,k))=fIn(NDX(I_NE,i-cx(I_NE,1),j-cx(I_NE,2),k-cx(I_NE,3)))
                  fOut(NDX(I_SW,i,j,k))=fIn(NDX(I_SW,i-cx(I_SW,1),j-cx(I_SW,2),k-cx(I_SW,3)))
                  fOut(NDX(I_SE,i,j,k))=fIn(NDX(I_SE,i-cx(I_SE,1),j-cx(I_SE,2),k-cx(I_SE,3)))
                  fOut(NDX(I_NW,i,j,k))=fIn(NDX(I_NW,i-cx(I_NW,1),j-cx(I_NW,2),k-cx(I_NW,3)))
                  fOut(NDX(I__T,i,j,k))=fIn(NDX(I__T,i-cx(I__T,1),j-cx(I__T,2),k-cx(I__T,3)))
                  fOut(NDX(I__B,i,j,k))=fIn(NDX(I__B,i-cx(I__B,1),j-cx(I__B,2),k-cx(I__B,3)))
                  fOut(NDX(I_TE,i,j,k))=fIn(NDX(I_TE,i-cx(I_TE,1),j-cx(I_TE,2),k-cx(I_TE,3)))
                  fOut(NDX(I_BW,i,j,k))=fIn(NDX(I_BW,i-cx(I_BW,1),j-cx(I_BW,2),k-cx(I_BW,3)))
                  fOut(NDX(I_BE,i,j,k))=fIn(NDX(I_BE,i-cx(I_BE,1),j-cx(I_BE,2),k-cx(I_BE,3)))
                  fOut(NDX(I_TW,i,j,k))=fIn(NDX(I_TW,i-cx(I_TW,1),j-cx(I_TW,2),k-cx(I_TW,3)))
                  fOut(NDX(I_TN,i,j,k))=fIn(NDX(I_TN,i-cx(I_TN,1),j-cx(I_TN,2),k-cx(I_TN,3)))
                  fOut(NDX(I_BS,i,j,k))=fIn(NDX(I_BS,i-cx(I_BS,1),j-cx(I_BS,2),k-cx(I_BS,3)))
                  fOut(NDX(I_BN,i,j,k))=fIn(NDX(I_BN,i-cx(I_BN,1),j-cx(I_BN,2),k-cx(I_BN,3)))
                  fOut(NDX(I_TS,i,j,k))=fIn(NDX(I_TS,i-cx(I_TS,1),j-cx(I_TS,2),k-cx(I_TS,3)))
   if(btest(state(i,j,k),wall) .eqv. .true.) then
                  ftmp = fOut(NDX(:,i,j,k)) 
                  fOut(NDX(opp(I__0),i,j,k))=ftmp(I__0)
                  fOut(NDX(opp(I__E),i,j,k))=ftmp(I__E)
                  fOut(NDX(opp(I__W),i,j,k))=ftmp(I__W)
                  fOut(NDX(opp(I__N),i,j,k))=ftmp(I__N)
                  fOut(NDX(opp(I__S),i,j,k))=ftmp(I__S)
                  fOut(NDX(opp(I_NE),i,j,k))=ftmp(I_NE)
                  fOut(NDX(opp(I_SW),i,j,k))=ftmp(I_SW)
                  fOut(NDX(opp(I_SE),i,j,k))=ftmp(I_SE)
                  fOut(NDX(opp(I_NW),i,j,k))=ftmp(I_NW)
                  fOut(NDX(opp(I__T),i,j,k))=ftmp(I__T)
                  fOut(NDX(opp(I__B),i,j,k))=ftmp(I__B)
                  fOut(NDX(opp(I_TE),i,j,k))=ftmp(I_TE)
                  fOut(NDX(opp(I_BW),i,j,k))=ftmp(I_BW)
                  fOut(NDX(opp(I_BE),i,j,k))=ftmp(I_BE)
                  fOut(NDX(opp(I_TW),i,j,k))=ftmp(I_TW)
                  fOut(NDX(opp(I_TN),i,j,k))=ftmp(I_TN)
                  fOut(NDX(opp(I_BS),i,j,k))=ftmp(I_BS)
                  fOut(NDX(opp(I_BN),i,j,k))=ftmp(I_BN)
                  fOut(NDX(opp(I_TS),i,j,k))=ftmp(I_TS)
endif
            enddo
         enddo
      enddo
#ifdef USE_OMP_FINE
      !$omp end parallel do
#endif /* USE_OMP_FINE */
#ifdef USE_OMP_TARGET_FINE
      !$omp end target teams distribute parallel do
#endif /* USE_OMP_TARGET_FINE */

      call cpu_time_measure(meas%tEnd_comm)
      meas%prop_duration = meas%tEnd_comm - meas%tSt_comm + meas%prop_duration

   end subroutine propagate
  !------------------------------------------------------------------------
#endif /*COMBINE_STREAM_COLLIDE */






#ifdef COMBINE_STREAM_COLLIDE
#endif /* COMBINE_STREAM_COLLIDE */
















#ifndef COMBINE_STREAM_COLLIDE
   subroutine collide_bgk(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the BGK collision step
! used by default



      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: u(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv
      real(R8B)                  :: rho0 
      real(R8B),pointer,dimension(:,:,:,:),intent(out)     :: fOut
      real(R8B),pointer,dimension(:,:,:,:),intent(inout)   :: fIn
      integer                    :: i,j,k,l


      call cpu_time_measure(meas%tSt_comm) 
      omega = s_par%omega
      rho0=s_par%rho0

      t2cs4inv = 1._R8B/(2._R8B*cs*cs*cs*cs)
      t2cs2inv = 1._R8B/(2._R8B*cs*cs)


    ! Single Relaxation time Collision routine
#ifdef USE_OMP_TARGET_FINE
      !$omp  target teams distribute parallel do collapse(2) &
#endif
#ifdef USE_OMP_FINE
      !$omp  parallel do collapse(2) &
#endif
#if defined(USE_OMP_TARGET_FINE) || defined (USE_OMP_FINE)
      !$omp& default(none) private(rholoc,fEq,u,usq,i,j,k) &
      !$omp& shared(fIn,fOut,t2cs2inv,t2cs4inv,s_par,lb_dom,cx,opp) &	
      !$omp& firstprivate( omega,rho0 )
#endif /* USE_OMP_TARGET_FINE || USE_OMP_FINE */
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
        !       if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then



! Calculate current Density 
       rholoc = fIn(NDX(I__0,i,j,k)) + &
                fIn(NDX(I__E,i,j,k)) + fIn(NDX(I__W,i,j,k)) + &
                fIn(NDX(I__N,i,j,k)) + fIn(NDX(I__S,i,j,k)) + &
                fIn(NDX(I__T,i,j,k)) + fIn(NDX(I__B,i,j,k)) + &
                fIn(NDX(I_NE,i,j,k)) + fIn(NDX(I_SW,i,j,k)) + &
                fIn(NDX(I_SE,i,j,k)) + fIn(NDX(I_NW,i,j,k)) + &
                fIn(NDX(I_TE,i,j,k)) + fIn(NDX(I_BW,i,j,k)) + &
                fIn(NDX(I_BE,i,j,k)) + fIn(NDX(I_TW,i,j,k)) + &
                fIn(NDX(I_TN,i,j,k)) + fIn(NDX(I_BS,i,j,k)) + &
                fIn(NDX(I_BN,i,j,k)) + fIn(NDX(I_TS,i,j,k))


! Choose which reference density is used. 
! Constant density rho0 for incompressible model

#ifdef INCOMPRESSIBLE
   rho0 = s_par%rho0 
#else
   rho0 = rholoc
#endif
  
! Calculate Momentum
    u(1) = fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k)) + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k)) + &
                      fIn(NDX(I_BE,i,j,k)) - fIn(NDX(I_TW,i,j,k))  
    u(2) = fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k)) + & 
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k)) + &
                      fIn(NDX(I_BN,i,j,k)) - fIn(NDX(I_TS,i,j,k))
    u(3) = fIn(NDX(I__T,i,j,k)) - fIn(NDX(I__B,i,j,k))  + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k))  + &
                      fIn(NDX(I_TW,i,j,k)) - fIn(NDX(I_BE,i,j,k))  + &
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k))  + &
                      fIn(NDX(I_TS,i,j,k)) - fIn(NDX(I_BN,i,j,k))
   u=u/rho0
 
               usq =  (u(1)**2 + u(2)**2 + u(3)**2)*t2cs2inv

! Calculate equilibrium distributions

               fEq(1) = t(1)*(rholoc -rho0*usq)
 
               fEq(2) = t(2)*(rholoc + rho0*((u(1))*cs2inv &
                      + (u(1))**2*t2cs4inv  &
                      - usq))
               fEq(3) = t(3)*(rholoc + rho0*((-u(1))*cs2inv &
                      + (-u(1))**2*t2cs4inv  &
                      - usq) )
               fEq(4) = t(4)*(rholoc + rho0*((u(2))*cs2inv &
                      + (u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(5) = t(5)*(rholoc + rho0*((-u(2))*cs2inv &
                      + (-u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(6) = t(6)*(rholoc + rho0*((u(3))*cs2inv &
                      + (u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(7) = t(7)*(rholoc + rho0*((-u(3))*cs2inv &
                      + (-u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(8) = t(8)*(rholoc + rho0*((u(1) + u(2))*cs2inv &
                      + (u(1) + u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(9) = t(9)*(rholoc + rho0*((-u(1) -u(2))*cs2inv &
                      + (-u(1) -u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(10) = t(10)*(rholoc + rho0*((u(1) -u(2))*cs2inv &
                      + (u(1) -u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(11) = t(11)*(rholoc + rho0*((-u(1) + u(2))*cs2inv &
                      + (-u(1) + u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(12) = t(12)*(rholoc + rho0*((u(1) + u(3))*cs2inv &
                      + (u(1) + u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(13) = t(13)*(rholoc + rho0*((-u(1) -u(3))*cs2inv &
                      + (-u(1) -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(14) = t(14)*(rholoc + rho0*((u(1) -u(3))*cs2inv &
                      + (u(1) -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(15) = t(15)*(rholoc + rho0*((-u(1) + u(3))*cs2inv &
                      + (-u(1) + u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(16) = t(16)*(rholoc + rho0*((u(2) + u(3))*cs2inv &
                      + (u(2) + u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(17) = t(17)*(rholoc + rho0*((-u(2) -u(3))*cs2inv &
                      + (-u(2) -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(18) = t(18)*(rholoc + rho0*((u(2) -u(3))*cs2inv &
                      + (u(2)  -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(19) = t(19)*(rholoc + rho0*((-u(2) + u(3))*cs2inv &
                      + (-u(2) + u(3))**2*t2cs4inv  &
                      - usq)) 

! Finally, relax to local equilibrium

               fOut(NDX(1,i,j,k)) = fIn(NDX(1,i,j,k)) - omega*(fIn(NDX(1,i,j,k)) - fEq(1))
               fOut(NDX(2,i,j,k)) = fIn(NDX(2,i,j,k)) - omega*(fIn(NDX(2,i,j,k)) - fEq(2))
               fOut(NDX(3,i,j,k)) = fIn(NDX(3,i,j,k)) - omega*(fIn(NDX(3,i,j,k)) - fEq(3))
               fOut(NDX(4,i,j,k)) = fIn(NDX(4,i,j,k)) - omega*(fIn(NDX(4,i,j,k)) - fEq(4))
               fOut(NDX(5,i,j,k)) = fIn(NDX(5,i,j,k)) - omega*(fIn(NDX(5,i,j,k)) - fEq(5))
               fOut(NDX(6,i,j,k)) = fIn(NDX(6,i,j,k)) - omega*(fIn(NDX(6,i,j,k)) - fEq(6))
               fOut(NDX(7,i,j,k)) = fIn(NDX(7,i,j,k)) - omega*(fIn(NDX(7,i,j,k)) - fEq(7))
               fOut(NDX(8,i,j,k)) = fIn(NDX(8,i,j,k)) - omega*(fIn(NDX(8,i,j,k)) - fEq(8))
               fOut(NDX(9,i,j,k)) = fIn(NDX(9,i,j,k)) - omega*(fIn(NDX(9,i,j,k)) - fEq(9))
               fOut(NDX(10,i,j,k)) = fIn(NDX(10,i,j,k)) - omega*(fIn(NDX(10,i,j,k)) - fEq(10))
               fOut(NDX(11,i,j,k)) = fIn(NDX(11,i,j,k)) - omega*(fIn(NDX(11,i,j,k)) - fEq(11))
               fOut(NDX(12,i,j,k)) = fIn(NDX(12,i,j,k)) - omega*(fIn(NDX(12,i,j,k)) - fEq(12))
               fOut(NDX(13,i,j,k)) = fIn(NDX(13,i,j,k)) - omega*(fIn(NDX(13,i,j,k)) - fEq(13))
               fOut(NDX(14,i,j,k)) = fIn(NDX(14,i,j,k)) - omega*(fIn(NDX(14,i,j,k)) - fEq(14))
               fOut(NDX(15,i,j,k)) = fIn(NDX(15,i,j,k)) - omega*(fIn(NDX(15,i,j,k)) - fEq(15))
               fOut(NDX(16,i,j,k)) = fIn(NDX(16,i,j,k)) - omega*(fIn(NDX(16,i,j,k)) - fEq(16))
               fOut(NDX(17,i,j,k)) = fIn(NDX(17,i,j,k)) - omega*(fIn(NDX(17,i,j,k)) - fEq(17))
               fOut(NDX(18,i,j,k)) = fIn(NDX(18,i,j,k)) - omega*(fIn(NDX(18,i,j,k)) - fEq(18))
               fOut(NDX(19,i,j,k)) = fIn(NDX(19,i,j,k)) - omega*(fIn(NDX(19,i,j,k)) - fEq(19))

#ifdef DEBUG
if(i==5 .and. j==8 .and. k==5) then
write(*,*) "u",u
write(*,*) "rholoc",rholoc
write(*,*) "fEq"
do l=1,nnod
write(*,*) fEq(l) 
enddo
write(*,*) "fIn"
do l=1,nnod
   write(*,*) fIn(NDX(l,5,8,5)) 
enddo
write(*,*) "fOut"
do l=1,nnod
   write(*,*) fOut(NDX(l,5,8,5)) 
enddo
endif
#endif /* DEBUG */

!else 
   ! simply copy the unrelaxed variables
!   fOut(NDX(:,i,j,k)) = fIn(NDX(:,i,j,k))
!endif
            end do
         end do
      end do
#ifdef USE_OMP_FINE
      !$omp end parallel do
#endif /* USE_OMP_FINE */
#ifdef USE_OMP_TARGET_FINE
      !$omp end target teams distribute parallel do
#endif /* USE_OMP_TARGET_FINE */
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_bgk 
  !------------------------------------------------------------------------
#endif /* COMBINE_STREAM_COLLIDE */







   subroutine stream_collide_bgk(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the BGK collision step
! with combined stream collide 


      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real (R8B)                 :: t2cs4inv,t2cs2inv
      real(R8B),pointer,dimension(:,:,:,:),intent(out)     :: fOut
      real(R8B),pointer,dimension(:,:,:,:),intent(inout)   ::  fIn
      integer                    :: i,j,k,l
#ifdef USE_OMP_FINE
      integer                    :: chunk
#endif /* USE_OMP_FINE */
#ifdef SPLITLOOPS
#define INDEXVAR i
#define SPLITDIMENSION lb_dom%lx(1)
#else
#define INDEXVAR 1
#define SPLITDIMENSION 1 
#endif
      real(R8B),dimension(SPLITDIMENSION) :: rholoc,u1,u2,u3,omega,usq
      real(R8B),dimension(SPLITDIMENSION) :: ftmp1,ftmp2,ftmp3,ftmp4,ftmp5,ftmp6,ftmp7,ftmp8,ftmp9
      real(R8B)                           :: fEq1,fEq2,fEq3,fEq4,fEq5,fEq6,fEq7,fEq8,fEq9
      real(R8B),dimension(SPLITDIMENSION) :: ftmp10,ftmp11,ftmp12,ftmp13,ftmp14,ftmp15,ftmp16,ftmp17,ftmp18,ftmp19
      real(R8B)                           :: fEq10,fEq11,fEq12,fEq13,fEq14,fEq15,fEq16,fEq17,fEq18,fEq19
      call cpu_time_measure(meas%tSt_comm) 

      omega = s_par%omega
      t2cs4inv = 1._R8B/(2._R8B*cs*cs*cs*cs)
      t2cs2inv = 1._R8B/(2._R8B*cs*cs)

      ! Single Relaxation time Collision routine

#ifdef USE_OMP_TARGET_FINE
      !$omp  target teams distribute parallel do collapse(2) &
#endif
#ifdef USE_OMP_FINE
      !$omp  parallel do collapse(2) &
#endif
#if defined(USE_OMP_TARGET_FINE) || defined (USE_OMP_FINE)
      !$omp& default(none) private(ftmp1,ftmp2,ftmp3,ftmp4,ftmp5,ftmp6,ftmp7,ftmp8,ftmp9,ftmp10,ftmp11, &
      !$omp& ftmp12,ftmp13,ftmp14,ftmp15,ftmp16,ftmp17,ftmp18,ftmp19, &
      !$omp& fEq1,fEq2,fEq3,fEq4,fEq5,fEq6,fEq7,fEq8,fEq9,fEq10,fEq11,fEq12,fEq13,fEq14,	&
      !$omp& fEq15,fEq16,fEq17,fEq18,fEq19,rholoc,u1,u2,u3,usq,i,j,k) &
      !$omp& shared(fIn,fOut,t2cs2inv,t2cs4inv,lb_dom,cx,opp) firstprivate( omega )
#endif /* USE_OMP_TARGET_FINE || USE_OMP_FINE */
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)

            !--------------------------------------------
            ! Do Streaming process ( Pull values) 

                  ftmp1(INDEXVAR)=fIn(NDX(1,i-cx(1,1),j-cx(1,2),k-cx(1,3)))
                  ftmp2(INDEXVAR)=fIn(NDX(2,i-cx(2,1),j-cx(2,2),k-cx(2,3)))
                  ftmp3(INDEXVAR)=fIn(NDX(3,i-cx(3,1),j-cx(3,2),k-cx(3,3)))
                  ftmp4(INDEXVAR)=fIn(NDX(4,i-cx(4,1),j-cx(4,2),k-cx(4,3)))
                  ftmp5(INDEXVAR)=fIn(NDX(5,i-cx(5,1),j-cx(5,2),k-cx(5,3)))
                  ftmp6(INDEXVAR)=fIn(NDX(6,i-cx(6,1),j-cx(6,2),k-cx(6,3)))
                  ftmp7(INDEXVAR)=fIn(NDX(7,i-cx(7,1),j-cx(7,2),k-cx(7,3)))
                  ftmp8(INDEXVAR)=fIn(NDX(8,i-cx(8,1),j-cx(8,2),k-cx(8,3)))
                  ftmp9(INDEXVAR)=fIn(NDX(9,i-cx(9,1),j-cx(9,2),k-cx(9,3)))
                  ftmp10(INDEXVAR)=fIn(NDX(10,i-cx(10,1),j-cx(10,2),k-cx(10,3)))
                  ftmp11(INDEXVAR)=fIn(NDX(11,i-cx(11,1),j-cx(11,2),k-cx(11,3)))
                  ftmp12(INDEXVAR)=fIn(NDX(12,i-cx(12,1),j-cx(12,2),k-cx(12,3)))
                  ftmp13(INDEXVAR)=fIn(NDX(13,i-cx(13,1),j-cx(13,2),k-cx(13,3)))
                  ftmp14(INDEXVAR)=fIn(NDX(14,i-cx(14,1),j-cx(14,2),k-cx(14,3)))
                  ftmp15(INDEXVAR)=fIn(NDX(15,i-cx(15,1),j-cx(15,2),k-cx(15,3)))
                  ftmp16(INDEXVAR)=fIn(NDX(16,i-cx(16,1),j-cx(16,2),k-cx(16,3)))
                  ftmp17(INDEXVAR)=fIn(NDX(17,i-cx(17,1),j-cx(17,2),k-cx(17,3)))
                  ftmp18(INDEXVAR)=fIn(NDX(18,i-cx(18,1),j-cx(18,2),k-cx(18,3)))
                  ftmp19(INDEXVAR)=fIn(NDX(19,i-cx(19,1),j-cx(19,2),k-cx(19,3)))

            ! Split loops. (Thomas: 3 und 5 splitloops)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then


         !--------------------------------------------
         ! Calculate macroscopic variables 

               rholoc(INDEXVAR)  = ftmp1(INDEXVAR)  + ftmp2(INDEXVAR)  + &
                         ftmp3(INDEXVAR)  + ftmp4(INDEXVAR)  + ftmp5(INDEXVAR) + &
                         ftmp6(INDEXVAR)  + ftmp7(INDEXVAR)  + &
                         ftmp8(INDEXVAR)  + ftmp9(INDEXVAR)  + ftmp10(INDEXVAR) + & 
                         ftmp11(INDEXVAR) + ftmp12(INDEXVAR) + &
                         ftmp13(INDEXVAR) + ftmp14(INDEXVAR) + ftmp15(INDEXVAR) + & 
                         ftmp16(INDEXVAR) + ftmp17(INDEXVAR) + &
                         ftmp18(INDEXVAR) + ftmp19(INDEXVAR)

               u1(INDEXVAR) = (  ftmp2(INDEXVAR)  - ftmp3(INDEXVAR)  + &
                         ftmp8(INDEXVAR)  - ftmp9(INDEXVAR)  + &
                         ftmp10(INDEXVAR) - ftmp11(INDEXVAR) + & 
                         ftmp12(INDEXVAR) - ftmp13(INDEXVAR) + &
                         ftmp14(INDEXVAR) - ftmp15(INDEXVAR)  & 
                             )/rholoc(INDEXVAR)    

               u2(INDEXVAR) = (  ftmp4(INDEXVAR)  - ftmp5(INDEXVAR) + &
                         ftmp8(INDEXVAR)  - ftmp9(INDEXVAR)  - ftmp10(INDEXVAR) + & 
                         ftmp11(INDEXVAR) + & 
                         ftmp16(INDEXVAR) - ftmp17(INDEXVAR) + &
                         ftmp18(INDEXVAR) - ftmp19(INDEXVAR))/rholoc(INDEXVAR)    

               u3(INDEXVAR) = (  ftmp6(INDEXVAR)  - ftmp7(INDEXVAR)  + &
                         ftmp12(INDEXVAR) - ftmp13(INDEXVAR)- ftmp14(INDEXVAR)  + & 
                         ftmp15(INDEXVAR) + ftmp16(INDEXVAR) - &
                         ftmp17(INDEXVAR) - ftmp18(INDEXVAR) + &
                         ftmp19(INDEXVAR))/rholoc(INDEXVAR)    
 
               usq(INDEXVAR)=  (u1(INDEXVAR)**2 + u2(INDEXVAR)**2 + u3(INDEXVAR)**2)*t2cs2inv


         !--------------------------------------------
         ! Calculate equilibrium distribution


               fEq1 = t(1)*rholoc(INDEXVAR)*(1._R8B  &
                      - usq(INDEXVAR)) 
               fEq2 = t(2)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq3 = t(3)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq4 = t(4)*rholoc(INDEXVAR)*(1._R8B + (u2(INDEXVAR))*cs2inv &
                      + (u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq5 = t(5)*rholoc(INDEXVAR)*(1._R8B + (-u2(INDEXVAR))*cs2inv &
                      + (-u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq6 = t(6)*rholoc(INDEXVAR)*(1._R8B + (u3(INDEXVAR))*cs2inv &
                      + (u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq7 = t(7)*rholoc(INDEXVAR)*(1._R8B + (-u3(INDEXVAR))*cs2inv &
                      + (-u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fOut(NDX(1,i,j,k)) = ftmp1(INDEXVAR) - omega(INDEXVAR)*(ftmp1(INDEXVAR) - fEq1)
               fOut(NDX(2,i,j,k)) = ftmp2(INDEXVAR) - omega(INDEXVAR)*(ftmp2(INDEXVAR) - fEq2)
               fOut(NDX(3,i,j,k)) = ftmp3(INDEXVAR) - omega(INDEXVAR)*(ftmp3(INDEXVAR) - fEq3)
               fOut(NDX(4,i,j,k)) = ftmp4(INDEXVAR) - omega(INDEXVAR)*(ftmp4(INDEXVAR) - fEq4)
               fOut(NDX(5,i,j,k)) = ftmp5(INDEXVAR) - omega(INDEXVAR)*(ftmp5(INDEXVAR) - fEq5)
               fOut(NDX(6,i,j,k)) = ftmp6(INDEXVAR) - omega(INDEXVAR)*(ftmp6(INDEXVAR) - fEq6)
               fOut(NDX(7,i,j,k)) = ftmp7(INDEXVAR) - omega(INDEXVAR)*(ftmp7(INDEXVAR) - fEq7)

#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#endif



               fEq8 = t(8)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) + u2(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) + u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq9 = t(9)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) -u2(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) -u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq10 = t(10)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) -u2(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) -u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq11 = t(11)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) + u2(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) + u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq12 = t(12)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq13 = t(13)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fOut(NDX(8,i,j,k))  = ftmp8(INDEXVAR)  - omega(INDEXVAR)*(ftmp8(INDEXVAR)  -  fEq8)
               fOut(NDX(9,i,j,k))  = ftmp9(INDEXVAR)  - omega(INDEXVAR)*(ftmp9(INDEXVAR)  -  fEq9)
               fOut(NDX(10,i,j,k)) = ftmp10(INDEXVAR) - omega(INDEXVAR)*(ftmp10(INDEXVAR) - fEq10)
               fOut(NDX(11,i,j,k)) = ftmp11(INDEXVAR) - omega(INDEXVAR)*(ftmp11(INDEXVAR) - fEq11)
               fOut(NDX(12,i,j,k)) = ftmp12(INDEXVAR) - omega(INDEXVAR)*(ftmp12(INDEXVAR) - fEq12)
               fOut(NDX(13,i,j,k)) = ftmp13(INDEXVAR) - omega(INDEXVAR)*(ftmp13(INDEXVAR) - fEq13)

#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#endif


               fEq14 = t(14)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq15 = t(15)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq16 = t(16)*rholoc(INDEXVAR)*(1._R8B + (u2(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (u2(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq17 = t(17)*rholoc(INDEXVAR)*(1._R8B + (-u2(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (-u2(INDEXVAR) -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq18 = t(18)*rholoc(INDEXVAR)*(1._R8B + (u2(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (u2(INDEXVAR)  -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq19 = t(19)*rholoc(INDEXVAR)*(1._R8B + (-u2(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (-u2(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fOut(NDX(14,i,j,k)) = ftmp14(INDEXVAR) - omega(INDEXVAR)*(ftmp14(INDEXVAR) - fEq14)
               fOut(NDX(15,i,j,k)) = ftmp15(INDEXVAR) - omega(INDEXVAR)*(ftmp15(INDEXVAR) - fEq15)
               fOut(NDX(16,i,j,k)) = ftmp16(INDEXVAR) - omega(INDEXVAR)*(ftmp16(INDEXVAR) - fEq16)
               fOut(NDX(17,i,j,k)) = ftmp17(INDEXVAR) - omega(INDEXVAR)*(ftmp17(INDEXVAR) - fEq17)
               fOut(NDX(18,i,j,k)) = ftmp18(INDEXVAR) - omega(INDEXVAR)*(ftmp18(INDEXVAR) - fEq18)
               fOut(NDX(19,i,j,k)) = ftmp19(INDEXVAR) - omega(INDEXVAR)*(ftmp19(INDEXVAR) - fEq19)



else

         !--------------------------------------------
         ! Do Bounce Back 

         fOut(NDX(opp(1),i,j,k)) = ftmp1(INDEXVAR)
         fOut(NDX(opp(2),i,j,k)) = ftmp2(INDEXVAR)
         fOut(NDX(opp(3),i,j,k)) = ftmp3(INDEXVAR)
         fOut(NDX(opp(4),i,j,k)) = ftmp4(INDEXVAR)
         fOut(NDX(opp(5),i,j,k)) = ftmp5(INDEXVAR)
         fOut(NDX(opp(6),i,j,k)) = ftmp6(INDEXVAR)
         fOut(NDX(opp(7),i,j,k)) = ftmp7(INDEXVAR)
#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif
         fOut(NDX(opp(8),i,j,k)) = ftmp8(INDEXVAR)
         fOut(NDX(opp(9),i,j,k)) = ftmp9(INDEXVAR)
         fOut(NDX(opp(10),i,j,k)) = ftmp10(INDEXVAR)
         fOut(NDX(opp(11),i,j,k)) = ftmp11(INDEXVAR)
         fOut(NDX(opp(12),i,j,k)) = ftmp12(INDEXVAR)
         fOut(NDX(opp(13),i,j,k)) = ftmp13(INDEXVAR)
#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif
         fOut(NDX(opp(14),i,j,k)) = ftmp14(INDEXVAR)
         fOut(NDX(opp(15),i,j,k)) = ftmp15(INDEXVAR)
         fOut(NDX(opp(16),i,j,k)) = ftmp16(INDEXVAR)
         fOut(NDX(opp(17),i,j,k)) = ftmp17(INDEXVAR)
         fOut(NDX(opp(18),i,j,k)) = ftmp18(INDEXVAR)
         fOut(NDX(opp(19),i,j,k)) = ftmp19(INDEXVAR)
      endif ! not wall
            end do
         end do
      end do
#ifdef USE_OMP_FINE
      !$omp end parallel do
#endif /* USE_OMP_FINE */
#ifdef USE_OMP_TARGET_FINE
      !$omp end target teams distribute parallel do
#endif /* USE_OMP_TARGET_FINE */

      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine stream_collide_bgk 
  !------------------------------------------------------------------------












end module lbm_functions


