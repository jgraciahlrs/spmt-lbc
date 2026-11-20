!!
!! Copyright (c) 2010        HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


module lb_init      
   use mpl_lib
   use lbmodel
   use function_set
   use timing
   use lb_geo
   use tools  
#include "include/replace.h"
contains
  !------------------------------------------------------------------------



   subroutine init_field_each(lb_dom,s_par,prc)
   !------------------------------------------------------------------------
   !
   ! initialization routine for initial conditions of domains
   ! this routine is intended for usage of more than one processes
   ! centralized initialization can be done by using init_field routine

      implicit none
      type(lb_block),intent(inout)  :: lb_dom
      type(sim_parameter)           :: s_par
      type(mpl_var)                 :: prc
      integer                 :: i,j,k,l
      integer                 :: lx(3),offset(3)

      lx(:) = lb_dom%lx(:)
      offset(:) = lb_dom%lx(:)*prc%crd(:) 


      ! set complete domain to density 1, u=0
       do k=1,lb_dom%lx(3)
          do j=1,lb_dom%lx(2)
             do i=1,lb_dom%lx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               lb_dom%rho(i,j,k)   = 1.0_R8B 

               lb_dom%u(NDX(1,i,j,k))   = 0._R8B
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
             end do
          end do
       end do


   select case(s_par%problem)


!----------------------------------------
! Default case 

   case default
      lb_dom%rho = 1.0d0
      lb_dom%u   = 0.0d0
      write(*,*) " Default case: Set complete domain to rho=1, u=0" 
   end select




  end subroutine init_field_each




#ifdef OLD_INIT_EACH


!--------------------------------------------------------------------
! now assign initial conditions depending on the problem 

      if(s_par%problem == cylinder) then
       ! Poisseuille Profile for Flow around Cylinder
      do k=1,lx(3)
        do j=1,lx(2)
          do i=1,lx(1)
               lb_dom%u(NDX(1,i,j,k)) = s_par%umax*get_parabolic3d( &
                  real(j+lb_dom%lx(2)*prc%crd(2),R8B),real(s_par%gx(2)+1,R8B), &
                  real(k+lb_dom%lx(3)*prc%crd(3),R8B),real(s_par%gx(3)+1,R8B))
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
               lb_dom%rho(i,j,k)   = 1._R8B
             end do
          end do
       end do
#ifdef MURKS
       do k=1,s_par%gx(3)
        do j=1,s_par%gx(2)
          do i=1,s_par%gx(1)
               lb_dom%u(NDX(1,i,j,k))   = 4._R8B*s_par%umax/real(s_par%gx(2)+1)**2  &
         & *(real(j*s_par%gx(2)+1)-real(j)*real(j))*(real(k*s_par%gx(3)+1)-real(k)*real(k))
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
               lb_dom%rho(i,j,k)   = 1._R8B

             end do
          end do
       end do
#endif
    else if(s_par%problem==cavity) then
       ! Driven Cavity
       do k=1,lb_dom%lx(3)
          do j=1,lb_dom%lx(2)
             do i=1,lb_dom%lx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               lb_dom%rho(i,j,k)   = 1.0_R8B 
               lb_dom%u(NDX(1,i,j,k))   = 0._R8B
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
             end do
          end do
       end do
      
   end if
!----------------------------------------
! Corotating Vortex 


   case(corotating_vortex)
      ! reset complete domain to eq
      rho      = 1.0d0
      u        = 0.0d0
      
      call init_vortex_macr(u,s_par,prc,lb_dom,(/0,0,0/))

  end subroutine init_field_each
  !------------------------------------------------------------------------
#endif /* OLD_INIT_EACH */


  subroutine nrbc_init(lb_dom,fIn,s_par,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
! this is non-diagonal for the corners

   implicit none
   type(lb_block     )    :: lb_dom
   type(sim_parameter)    :: s_par
   type(mpl_var      )    :: prc  
   real(R8B),intent(in)       :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   integer                :: i,ind
   integer                :: x,y,z ,x1,y1,z1,x2,y2,z2  

!   allocate(lb_dom%nrfplus(lb_dom%nnr_wall,nnod))
   do i=1,lb_dom%nnr_wall   
      x = lb_dom%obs_nrwall(i,1)
      y = lb_dom%obs_nrwall(i,2)
      z = lb_dom%obs_nrwall(i,3)
      x1= lb_dom%obs_nrwall(i,4)
      y1= lb_dom%obs_nrwall(i,5)
      z1= lb_dom%obs_nrwall(i,6)
      x2= lb_dom%obs_nrwall(i,7)
      y2= lb_dom%obs_nrwall(i,8)
      z2= lb_dom%obs_nrwall(i,9)
      
!  do ind=1,nnod
!     lb_dom%nrfplus(i,ind) = 0.5*(fIn(NDX(ind,x,y,z)) + fIn(NDX(opp(ind),x,y,z)))
!  enddo

   end do
   end subroutine nrbc_init
  !------------------------------------------------------------------------





  subroutine nrbc_find_neighbors(lb_dom,s_par,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
! this is non-diagonal for the corners

   implicit none
   type(lb_block     )    :: lb_dom
   type(sim_parameter)    :: s_par
   type(mpl_var      )    :: prc  
   integer                :: i
   integer                :: x,y,z ,x1,y1,z1,x2,y2,z2,x3,y3,z3  
   integer                :: nsp 
   real(R8B)              :: outer_omega,om_sponge(0:3)

   do i=1,lb_dom%nnr_wall   
      x = lb_dom%obs_nrwall(i,1)
      y = lb_dom%obs_nrwall(i,2)
      z = lb_dom%obs_nrwall(i,3)

      lb_dom%obs_nrwall(i,0) = 0  ! integer with directions (set by setting the bits)

      ! Check: Is Inlet Set?? If not, outlet. If yes, inlet
      if(btest(lb_dom%state(x,y,z),nr_wall_in) .eqv. .true.) then
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),0)
      endif

!#!define NRBC_DIAGONAL
#ifdef NRBC_DIAGONAL
      if(x == 1 .and. prc%crd(1)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,4) = x+1 
         lb_dom%obs_nrwall(i,7) = x+2
         lb_dom%obs_nrwall(i,10) = x+3
         lb_dom%obs_nrwall(i,13) = x+4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dWEST)
      elseif(x == lb_dom%lx(1) .and. prc%crd(1)==prc%np(1)-1) then ! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,4) = x-1
         lb_dom%obs_nrwall(i,7) = x-2
         lb_dom%obs_nrwall(i,10) = x-3
         lb_dom%obs_nrwall(i,13) = x-4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dEAST)
      else
         lb_dom%obs_nrwall(i,4) = x
         lb_dom%obs_nrwall(i,7) = x
         lb_dom%obs_nrwall(i,10) = x
         lb_dom%obs_nrwall(i,13) = x
      endif

      if(y == 1 .and. prc%crd(2)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,5) = y+1 
         lb_dom%obs_nrwall(i,8) = y+2 
         lb_dom%obs_nrwall(i,11) = y+3 
         lb_dom%obs_nrwall(i,14) = y+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dSOUTH)
      elseif(y == lb_dom%lx(2) .and. prc%crd(2)==prc%np(2)-1) then! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,5) = y-1
         lb_dom%obs_nrwall(i,8) = y-2
         lb_dom%obs_nrwall(i,11) = y-3 
         lb_dom%obs_nrwall(i,14) = y-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dNORTH)
      else
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
      endif 

      if(z == 1.and.  prc%crd(3)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,6) = z+1 
         lb_dom%obs_nrwall(i,9) = z+2 
         lb_dom%obs_nrwall(i,12) = z+3 
         lb_dom%obs_nrwall(i,15) = z+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dBOTTOM)
      elseif(z == lb_dom%lx(3).and. prc%crd(3)==prc%np(3)-1 ) then
         lb_dom%obs_nrwall(i,6) = z-1
         lb_dom%obs_nrwall(i,9) = z-2
         lb_dom%obs_nrwall(i,12) = z-3 
         lb_dom%obs_nrwall(i,15) = z-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dTOP)
      else ! if not, use same coordinate
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      endif  ! z


#else /* NRBC_DIAGONAL */
!----------------------------
! x: check for x border
      if(x == 1 .and. prc%crd(1)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,4) = x+1 
         lb_dom%obs_nrwall(i,7) = x+2
         lb_dom%obs_nrwall(i,10) = x+3
         lb_dom%obs_nrwall(i,13) = x+4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dWEST)
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      elseif(x == lb_dom%lx(1) .and. prc%crd(1)==prc%np(1)-1) then ! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,4) = x-1
         lb_dom%obs_nrwall(i,7) = x-2
         lb_dom%obs_nrwall(i,10) = x-3
         lb_dom%obs_nrwall(i,13) = x-4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dEAST)
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      else ! if not, use same coordinate

!----------------------------
! y if not x-border check y coordinate
         lb_dom%obs_nrwall(i,4) = x
         lb_dom%obs_nrwall(i,7) = x
         lb_dom%obs_nrwall(i,10) = x
         lb_dom%obs_nrwall(i,13) = x
      if(y == 1 .and. prc%crd(2)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,5) = y+1 
         lb_dom%obs_nrwall(i,8) = y+2 
         lb_dom%obs_nrwall(i,11) = y+3 
         lb_dom%obs_nrwall(i,14) = y+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dSOUTH)
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      elseif(y == lb_dom%lx(2) .and. prc%crd(2)==prc%np(2)-1) then! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,5) = y-1
         lb_dom%obs_nrwall(i,8) = y-2
         lb_dom%obs_nrwall(i,11) = y-3 
         lb_dom%obs_nrwall(i,14) = y-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dNORTH)
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      else

!----------------------------
! z:if neither on x- nor on y-border, check z coordinate
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
      if(z == 1.and.  prc%crd(3)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,6) = z+1 
         lb_dom%obs_nrwall(i,9) = z+2 
         lb_dom%obs_nrwall(i,12) = z+3 
         lb_dom%obs_nrwall(i,15) = z+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dBOTTOM)
      elseif(z == lb_dom%lx(3).and. prc%crd(3)==prc%np(3)-1 ) then
         lb_dom%obs_nrwall(i,6) = z-1
         lb_dom%obs_nrwall(i,9) = z-2
         lb_dom%obs_nrwall(i,12) = z-3 
         lb_dom%obs_nrwall(i,15) = z-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dTOP)
      else ! if not, use same coordinate
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      endif  ! z
      endif  ! y
      endif  ! x
#endif      


!      write(31,*) x,y,z ,"  >  ",lb_dom%obs_nrwall(i,4),lb_dom%obs_nrwall(i,5),lb_dom%obs_nrwall(i,6)
      x1= lb_dom%obs_nrwall(i,4)
      y1= lb_dom%obs_nrwall(i,5)
      z1= lb_dom%obs_nrwall(i,6)
      x2= lb_dom%obs_nrwall(i,7)
      y2= lb_dom%obs_nrwall(i,8)
      z2= lb_dom%obs_nrwall(i,9)
      x3= lb_dom%obs_nrwall(i,10)
      y3= lb_dom%obs_nrwall(i,11)
      z3= lb_dom%obs_nrwall(i,12)

      ! assign zero values
      lb_dom%nrwall_0val(i,0)          = lb_dom%rho(x,y,z)
      lb_dom%nrwall_0val(i,1:NDIM    ) = lb_dom%u(NDX(1:NDIM,x,y,z))
      lb_dom%nrwall_prev(i,0)          = lb_dom%rho(x,y,z)
      lb_dom%nrwall_prev(i,1:NDIM    ) = lb_dom%u(NDX(1:NDIM,x,y,z))
      lb_dom%nrwall_prev(i,4)          = lb_dom%rho(x1,y1,z1)
      lb_dom%nrwall_prev(i,5:4+NDIM  ) = lb_dom%u(NDX(1:NDIM,x1,y1,z1))
      lb_dom%nrwall_prev(i,8)          = lb_dom%rho(x2,y2,z2)
      lb_dom%nrwall_prev(i,9:8+NDIM  ) = lb_dom%u(NDX(1:NDIM,x2,y2,z2))
 
#ifdef SPONGE_NRBC
   !----------------------------------------
   ! Increase viscosity on border nodes to increase stability
   
   lb_dom%omega(x, y, z ) = om_sponge(0)
   lb_dom%omega(x1,y1,z1) = om_sponge(1)
   lb_dom%omega(x2,y2,z2) = om_sponge(2)
   lb_dom%omega(x3,y3,z3) = om_sponge(3)

#endif /* SPONGE_NRBC */
   end do
   end subroutine nrbc_find_neighbors
  !------------------------------------------------------------------------





  subroutine per_find_neighbors(lb_dom,s_par,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
!

   implicit none
   type(lb_block     )      :: lb_dom
   type(sim_parameter)      :: s_par
   type(mpl_var      )      :: prc  
   integer                  :: i,j,k,l
   integer                  :: x,y,z 

   do i=1,lb_dom%nper_wall   
      x = lb_dom%obs_per(i,1)
      y = lb_dom%obs_per(i,2)
      z = lb_dom%obs_per(i,3)
      lb_dom%obs_per(i,0) = 0  ! integer with directions (set by setting the bits)
      ! x coordinate
      if(x == 1) then ! if on left boarder, use right neighbors
         lb_dom%obs_per(i,4) = lb_dom%lx(1)-1 
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),3)
      elseif(x == lb_dom%lx(1)) then ! if on right border, use left neighbors
         lb_dom%obs_per(i,4) = 2
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),1)
      else ! if not, use same coordinate
         lb_dom%obs_per(i,4) = x
      endif
      ! y coordinate
      if(y == 1) then ! if on left boarder, use right neighbors
         lb_dom%obs_per(i,5) = lb_dom%lx(2)-1 
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),4)
      elseif(y == lb_dom%lx(2)) then! if on right border, use left neighbors
         lb_dom%obs_per(i,5) = 2  
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),2)
      else ! if not, use same coordinate
         lb_dom%obs_per(i,5) = y
      endif
      ! z coordinate
      if(z == 1) then ! if on left boarder, use right neighbors
         lb_dom%obs_per(i,6) = lb_dom%lx(3)-1 
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),6)
      elseif(z == lb_dom%lx(3)) then
         lb_dom%obs_per(i,6) = 2  
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),5)
      else ! if not, use same coordinate
         lb_dom%obs_per(i,6) = z
      endif
   end do
   end subroutine per_find_neighbors
  !------------------------------------------------------------------------







   subroutine init_vector(lb_dom,s_par,prc)

!------------------------------------------------------------------------
! create a vector for the boundary conditions containing the coordinates 
! of the boundary node (and additionally first and second neighbors for nrbc)
!


      implicit none
      type(lb_block)                :: lb_dom
      type(sim_parameter)           :: s_par
      type(mpl_var)                 :: prc  
      integer      :: i,j,k,l
      integer      :: sum_nnr,sum_nper
      integer      :: counter,counter_sponge,counter_reset,counter_nrwall,counter_perwall

      lb_dom%nobs = 0
!      lb_dom%ninlet= 0
!      lb_dom%noutlet= 0
!      lb_dom%nobs_sponge = 0
!      lb_dom%nobs_reset = 0
      lb_dom%nnr_wall = 0
      lb_dom%nper_wall = 0

      ! count how many of each boundaries are present
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
                 lb_dom%nobs = lb_dom%nobs + 1
               end if
!               if(btest(lb_dom%state(i,j,k),sponge_wall)) then
!                 lb_dom%nobs_sponge = lb_dom%nobs_sponge + 1
!               end if
!               if(btest(lb_dom%state(i,j,k),reset_wall)) then
!                 lb_dom%nobs_reset = lb_dom%nobs_reset + 1
!               end if
               if(btest(lb_dom%state(i,j,k),nr_wall) .or. &
                  btest(lb_dom%state(i,j,k),nr_wall_in) ) then
                  if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
                    lb_dom%nnr_wall = lb_dom%nnr_wall + 1
                  end if
               end if
               if(btest(lb_dom%state(i,j,k),per_wall)) then
                 lb_dom%nper_wall = lb_dom%nper_wall + 1
               end if
!               if(btest(lb_dom%state(i,j,k),inlet ) ) then
!                 lb_dom%ninlet  = lb_dom%ninlet   + 1
!               end if
!               if(btest(lb_dom%state(i,j,k),outlet) ) then
!                 lb_dom%noutlet  = lb_dom%noutlet  + 1
!               end if


            end do
         end do
      end do

      allocate(lb_dom%obs(lb_dom%nobs,3))
      allocate(lb_dom%obs_per(lb_dom%nper_wall,0:6))
      allocate(lb_dom%obs_nrwall(lb_dom%nnr_wall,0:15))
               ! 0     ... holds bitmap information, which border
               ! 1..3  ... coordinates of nrbc node
               ! 4..12 ... hold 1st, 2nd, 3rd neighbor node coordinates 
               ! 13.15 ... hold neighboring ghost node                  
      allocate(lb_dom%nrwall_0val(lb_dom%nnr_wall,0:3))  ! holds values of rho and u of bc
      allocate(lb_dom%nrwall_prev(lb_dom%nnr_wall,0:11))  ! holds lodi values of rho and u 
!      allocate(lb_dom%inlet(lb_dom%ninlet,3))
!      allocate(lb_dom%outlet(lb_dom%noutlet,3))

      counter        = 0
      counter_sponge = 0
      counter_reset  = 0
      counter_nrwall = 0
      counter_perwall= 0
!      counter_ninlet= 0
!      counter_noutlet=0
      
      ! then assign the boundaries with their coordinates
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
                  counter = counter+1
                  lb_dom%obs(counter,1:3) = (/i,j,k/)
               end if
!               if(btest(lb_dom%state(i,j,k),sponge_wall)) then
!                  counter_sponge = counter_sponge+1
!                  lb_dom%obs_sponge(counter_sponge,1:3) = (/i,j,k/)
!               end if
!               if(btest(lb_dom%state(i,j,k),reset_wall)) then
!                  counter_reset = counter_reset+1
!                  lb_dom%obs_reset(counter_reset,1:3) = (/i,j,k/)
!               end if
               if(btest(lb_dom%state(i,j,k),nr_wall) .or. &
                  btest(lb_dom%state(i,j,k),nr_wall_in)) then
                  if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
                     counter_nrwall = counter_nrwall +1
                     lb_dom%obs_nrwall(counter_nrwall,1:3) = (/i,j,k/)
                  end if
               end if
               if(btest(lb_dom%state(i,j,k),per_wall)) then
                  counter_perwall = counter_perwall +1
                  lb_dom%obs_per(counter_perwall,1:3) = (/i,j,k/)
               end if
!               if(btest(lb_dom%state(i,j,k),inlet ) ) then
!                  counter_ninlet= counter_ninlet+1
!                  lb_dom%inlet(counter_ninlet,1:3) = (/i,j,k/)
!               end if

            end do
         end do
      end do

#ifdef USE_MPI
       call MPI_REDUCE(lb_dom%nnr_wall, sum_nnr, 1,mpi_integer,MPI_SUM,prc%root_th,prc%cart_comm,prc%ierr)
       call MPI_REDUCE(lb_dom%nper_wall,sum_nper,1,mpi_integer,MPI_SUM,prc%root_th,prc%cart_comm,prc%ierr)
#else
       sum_nnr  = lb_dom%nnr_wall
       sum_nper = lb_dom%nper_wall
#endif



       if(prc%rk==0) then
          write(*,*) 
          write(*,*) "Setting Boundary conditions resulted in"
          write(*,'(a10,i6,a14,i6)') " # NRBC:  ",sum_nnr,"   # Periodic: ",sum_nper
          write(*,*) 
       endif

      if (lb_dom%nper_wall > 0) then
         if(prc%size > 1) then 
            write(*,*) "Periodic functions are only working with serial execution!"
            stop
         endif
         call per_find_neighbors(lb_dom,s_par,prc)
      end if




#ifdef SPARSE /* SPARSE */
! this is a preparation for the usage of indirect addressing. 
! not used until now
   ! allocate the number of (nodes + ghost nodes - solid nodes)*nnod
   allocate(lb_dom%fsIn(((lb_dom%lx(1)+2*prc%n_ghost_nodes)*&
      & (lb_dom%lx(1)+2*prc%n_ghost_nodes)*(lb_dom%lx(1)+2*prc%n_ghost_nodes)-lb_dom%nobs)*nnod)
   ! assign the densities one by one
      counter = 0
      do k=1-prc%n_ghost_nodes,lb_dom%lx(3)+prc%n_ghost_nodes
         do j=1-prc%n_ghost_nodes,lb_dom%lx(2)+prc%n_ghost_nodes
            do i=1-prc%n_ghost_nodes,lb_dom%lx(1)+prc%n_ghost_nodes
               if(btest(lb_dom%state(i,j,k),fluid)) then
                  counter = counter+1
                   
               end if
            end do
         end do
      end do

   ! assign adjacency list
!FIXME ghost nodes also have to be defined as fluid nodes
      ! create the target vector
      counter = 0
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               counter = counter+1
               do l=1,nnod
                  lb_dom%f_sparse((counter*nnod-1)+l) = lb_dom%fIn(NDX(l,i,j,k))
               end do
            end do
         end do
      end do
#endif

  end subroutine init_vector
  !------------------------------------------------------------------------
  





   subroutine calc_rho_ges(lb_dom,fIn,rho_ges,tStep,prc)
!------------------------------------------------------------------------
! calculate the total density of the domain
! This can be used for convergence checks



      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(mpl_var)              :: prc
      real(R8B),intent(in)       :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B),intent(out)      :: rho_ges
      integer,intent(in)         :: tStep  
      integer                    :: i,j,k,l
      real(R8B)                  :: sum_rho,uc(3),umax(3)
      sum_rho  = 0.0d0
      rho_ges  = 0.0d0
      umax     = 0.0d0


      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
               if(btest(lb_dom%state(i,j,k),lid) .eqv. .false.) then
               uc = 0.0d0
               do l=1,nnod
                  rho_ges=rho_ges+fIn(NDX(l,i,j,k))
               !   if(fIn(NDX(l,i,j,k)) .le. 0.) write(*,*) 'found l 0  at',i,j,k,l
               enddo
#ifdef CALC_MAX_VELOCITY
                  uc(1) = cx(l,1)*fIn(NDX(l,i,j,k)) 
                  uc(2) = cx(l,2)*fIn(NDX(l,i,j,k))
                  uc(3) = cx(l,3)*fIn(NDX(l,i,j,k)) 
               enddo
               if(abs(uc(1)) > umax(1)) then
                  umax(1)=abs(uc(1))
               endif
               if(abs(uc(2)) > umax(2)) then
                  umax(2)=abs(uc(2))
               endif
               if(abs(uc(3)) > umax(3)) then
                  umax(3)=abs(uc(3))
               endif
#endif /* CALC_MAX_VELOCITY */
               endif
               endif
            enddo
         enddo
      enddo

      if(prc%size > 0) then
#ifdef USE_MPI
         ! if parallel execution, sum up the calculated values and sen to process 0
         call MPI_REDUCE(rho_ges,sum_rho,1,mpi_double_precision,MPI_SUM,prc%root_th,prc%cart_comm,prc%ierr)
         rho_ges=sum_rho
#endif /* MPI */
      end if
      if(prc%rk == 0) then
#ifdef CALC_MAX_VELOCITY
         write(*,'(a16,f10.5,a4,f10.5,a4,f10.5)') "max velocity  x:",umax(1),"  y:",umax(2),"  z:",umax(3)
#endif /* CALC_MAX_VELOCITY */
         write(*,'(a16,f16.6,a12,i7)') "Total density:  ",rho_ges," at timestep ",tStep
         open(43,file='dens.inf',position='append')
         write(43,*) tStep,rho_ges
         close(43)
      end if


   end subroutine calc_rho_ges
  !------------------------------------------------------------------------
  



   subroutine calc_fEq(lb_dom,rho,u,fEq,xmin,xmax,ymin,ymax,zmin,zmax,s_par)
   !------------------------------------------------------------------------
   !
   ! calculate the Maxwellian / equlibrium distribution function
   !

      type(sim_parameter)           :: s_par
      type(lb_block) :: lb_dom
      real(R8B)      :: usq,cs2inv,t2cs4inv,t2cs2inv,rholoc,rho0
      real(R8B)      :: u(NDX(NDIM,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B)      :: rho(0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1))
      real(R8B)      :: fEq(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer        :: i,j,k,l,err
      integer        :: xmin,xmax,ymin,ymax,zmin,zmax
      integer        :: x_l,x_u,y_l,y_u,z_l,z_u

      intent(in)     :: u,rho
      intent(inout)  :: fEq,lb_dom

      !---------------------------------------------------
      ! define upper and lower limits
      ! if no limit given (0), use lower: 1 or upper: lb_dom%lx(i)
      if(xmin==0) then;  x_l=1;  else; x_l=xmin; end if
      if(xmax==0) then;  x_u=lb_dom%lx(1); else; x_u=xmax; end if
      if(ymin==0) then;  y_l=1;  else; y_l=ymin; end if
      if(ymax==0) then;  y_u=lb_dom%lx(2); else; y_u=ymax; end if
      if(zmin==0) then;  z_l=1;  else; z_l=zmin; end if
      if(zmax==0) then;  z_u=lb_dom%lx(3); else; z_u=zmax; end if
      err = 0

      !---------------------------------------------------
      ! define constants 
      cs2inv  = 1._R8B/cs**2
      t2cs4inv = 1._R8B/(2*cs**4)
      t2cs2inv = 1._R8B/(2*cs**2)

    do k=z_l,z_u
       do j=y_l,y_u
          do i=x_l,x_u

          usq =  (u(NDX(1,i,j,k))**2 + u(NDX(2,i,j,k))**2 + u(NDX(3,i,j,k))**2)*t2cs2inv
          rholoc   = rho(i,j,k) 
#ifdef INCOMPRESSIBLE
             rho0 = s_par%rho0
#else
             rho0 = rholoc
#endif 

           fEq(NDX(1,i,j,k)) = t(1)*rholoc*(1._R8B  &
                - usq) 
           fEq(NDX(2,i,j,k)) = t(2)*(rholoc + rho0*((u(NDX(1,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)))**2*t2cs4inv  &
                - usq))
           fEq(NDX(3,i,j,k)) = t(3)*(rholoc + rho0*((-u(NDX(1,i,j,k)))*cs2inv &
                + u(NDX(1,i,j,k))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(4,i,j,k)) = t(4)*(rholoc + rho0*((u(NDX(2,i,j,k)))*cs2inv &
                + (u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(5,i,j,k)) = t(5)*(rholoc + rho0*((-u(NDX(2,i,j,k)))*cs2inv &
                + (-u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(6,i,j,k)) = t(6)*(rholoc + rho0*((u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(7,i,j,k)) = t(7)*(rholoc + rho0*((-u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(8,i,j,k)) = t(8)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(9,i,j,k)) = t(9)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(10,i,j,k)) = t(10)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(11,i,j,k)) = t(11)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(12,i,j,k)) = t(12)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(13,i,j,k)) = t(13)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(14,i,j,k)) = t(14)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(15,i,j,k)) = t(15)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(16,i,j,k)) = t(16)*(rholoc  &
                + rho0*((u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(17,i,j,k)) = t(17)*(rholoc  &
                + rho0*((-u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(18,i,j,k)) = t(18)*(rholoc  &
                + rho0*((u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(19,i,j,k)) = t(19)*(rholoc  &
                + rho0*((-u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
               
          end do
       end do
    end do
  end subroutine calc_fEq
  !------------------------------------------------------------------------



  subroutine get_geo_each(lb_dom,s_par,prc)
    !------------------------------------------------------------------------
    !
    ! set up the geometry, each process by itself
    !
      implicit none
      type(lb_block),intent(inout)      :: lb_dom 
      type(sim_parameter),intent(inout) :: s_par
      type(mpl_var),intent(inout)       :: prc  
      integer        :: i,j,k,l
      integer        :: x_start,x_end,y_start,y_end,z_start,z_end
      integer        :: lx,ly,lz,gx,gy,gz
      integer        :: offset(3)
      ! Fill up the Coordinates for the Grid
      lx = lb_dom%lx(1)
      ly = lb_dom%lx(2)
      lz = lb_dom%lx(3)
      gx = s_par%gx(1)
      gy = s_par%gx(2)
      gz = s_par%gx(3)
      offset(:) = lb_dom%lx(:)*prc%crd(:) 

      if(prc%crd(1) == 0)          then; x_start = 2;    else; x_start=1;  end if
      if(prc%crd(1) == prc%np(1)-1) then; x_end   = lx-1; else; x_end  =lx; end if
      if(prc%crd(2) == 0)          then; y_start = 2;    else; y_start=1;  end if
      if(prc%crd(2) == prc%np(2)-1) then; y_end   = ly-1; else; y_end  =ly; end if
      if(prc%crd(3) == 0)          then; z_start = 2;    else; z_start=1;  end if
      if(prc%crd(3) == prc%np(3)-1) then; z_end   = lz-1; else; z_end  =lz; end if

      do i=1,prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,2)-prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,1)+1
         lb_dom%x(i)  = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,1)+real(i,8) - 2
      end do
      do j=1,prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,2)-prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,1)+1
         lb_dom%y(j)  = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,1)+real(j,8) - 2
      end do
      do k=1,prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,2)-prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,1)+1
         lb_dom%z(k)  = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,1)+real(k,8) - 2
      end do    

      ! Set all nodes to fluid
      lb_dom%state(:,:,:)   = 0 
      lb_dom%state = ibset(lb_dom%state,fluid) 




!--------------------------------------------------------------------
! Select defined problem  


      select case(s_par%problem)



!----------------------------------------
! Lid driven Cavity 

      case(cavity)
         ! where is the lid? y has to be ly, x and z have to be 2 or 1, depending on where the prc%crd is
         if(prc%crd(2) == prc%np(2) - 1) then 
             lb_dom%state(  x_start:x_end,  ly,  z_start:z_end  )  = &
         ibset(lb_dom%state(  x_start:x_end,  ly,  z_start:z_end  ),lid)
         end if
         ! left wall (at y- border)
!FIXME
         if(prc%crd(1) == 0)           lb_dom%state(1,:,:)     = ibset(lb_dom%state(1,:,:),wall)
         if(prc%crd(1) == prc%np(1)-1) lb_dom%state(lx,:,:)    = ibset(lb_dom%state(lx,:,:),wall)
         ! bottom wall
         if(prc%crd(2) == 0)           lb_dom%state(:,1,:)     = ibset(lb_dom%state(:,1,:) ,wall)
         if(NDIM==3) then
           ! set also front and back wall to wall condition
            if(prc%crd(3) == 0)        lb_dom%state(:,:,1)     = ibset(lb_dom%state(:,:,1),wall)
            if(prc%crd(3) == prc%np(3)-1) lb_dom%state(:,:,lz) = ibset(lb_dom%state(:,:,lz),wall )        
         end if


   end select

! get number of inlet nodes in directions
   end subroutine get_geo_each
   !------------------------------------------------------------------------

  

   subroutine calc_macr_vals(lb_dom,u,rho,fIn,meas)
   !------------------------------------------------------------------------
   !
   ! calculate the macroscopic values from the microscopic distributions 
   !
     implicit none
     type(lb_block),intent(inout)  :: lb_dom
     type(measure )                :: meas   
     real(R8B), intent(out)    :: u(NDX(NDIM,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
     real(R8B), intent(out)    :: rho(0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1))
     real(R8B), intent(in)     :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
     integer              :: i,j,k,l
     
         call cpu_time_measure(meas%tSt_comm) 

     do k=1,lb_dom%lx(3)
        do j=1,lb_dom%lx(2)
           do i=1,lb_dom%lx(1)
!               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
!if(lb_dom%state(i,j,k) /= wall) then
       rho(i,j,k)  =  fIn(NDX(I__0,i,j,k)) + &
                      fIn(NDX(I__E,i,j,k)) + fIn(NDX(I__W,i,j,k)) + &
                      fIn(NDX(I__N,i,j,k)) + fIn(NDX(I__S,i,j,k)) + &
                      fIn(NDX(I__T,i,j,k)) + fIn(NDX(I__B,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) + fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_SE,i,j,k)) + fIn(NDX(I_NW,i,j,k)) + &
                      fIn(NDX(I_TE,i,j,k)) + fIn(NDX(I_BW,i,j,k)) + &
                      fIn(NDX(I_BE,i,j,k)) + fIn(NDX(I_TW,i,j,k)) + &
                      fIn(NDX(I_TN,i,j,k)) + fIn(NDX(I_BS,i,j,k)) + &
                      fIn(NDX(I_BN,i,j,k)) + fIn(NDX(I_TS,i,j,k))
    u(NDX(1,i,j,k)) = fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k)) + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k)) + &
                      fIn(NDX(I_BE,i,j,k)) - fIn(NDX(I_TW,i,j,k))  
    u(NDX(2,i,j,k)) = fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k)) + & 
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k)) + &
                      fIn(NDX(I_BN,i,j,k)) - fIn(NDX(I_TS,i,j,k))
    u(NDX(3,i,j,k)) = fIn(NDX(I__T,i,j,k)) - fIn(NDX(I__B,i,j,k))  + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k))  + &
                      fIn(NDX(I_TW,i,j,k)) - fIn(NDX(I_BE,i,j,k))  + &
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k))  + &
                      fIn(NDX(I_TS,i,j,k)) - fIn(NDX(I_BN,i,j,k))


               u(NDX(:,i,j,k))=u(NDX(:,i,j,k))/rho(i,j,k)
!endif

            enddo
         enddo
      enddo
         call cpu_time_measure(meas%tEnd_comm)
         meas%ccmv_duration = meas%tEnd_comm - meas%tSt_comm + meas%ccmv_duration
   end subroutine calc_macr_vals
   !------------------------------------------------------------------------









  subroutine treat_initial(lb_dom,s_par,tstep,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
!

   implicit none
   type(lb_block)       :: lb_dom
   type(sim_parameter)  :: s_par
   integer, intent(in)  :: tstep 
   type(mpl_var)        :: prc  
   integer              :: i,j
   integer              :: x,y,z
   logical              :: current_status,status_changed


   current_status = s_par%initial 
   status_changed = .false.
   if(tstep <= s_par%tInitial ) then 
      s_par%initial=.true. 
   else 
      s_par%initial=.false.
   endif 

   
if(s_par%initial .eqv. .true.) then
if(tstep ==1) then
   write(*,*) "assigning zero field..."
   if( &
#ifdef F2003_ALLOC_DERIVED_TYPE
   allocated(  &
#else
   associated(  &
#endif
lb_dom%u0)) then
      lb_dom%u0   = lb_dom%u
   else
      write(*,*) "array u0 is not allocated! Stopping..."
      stop
   endif
   write(*,*) "Done."
endif
endif

   if(s_par%initial .eqv. .false.) then 
   if(current_status .eqv. .true.) then 
      status_changed = .true.
   endif; endif

if(status_changed) then
if(prc%rk == 0) then
   write(*,*) "------------------------ "
   write(*,*) "Initial Relaxation done. "
   write(*,*) "Here, a Norm calculation should be implemented to check, if pressure field has converged"
   write(*,*) "This routine must be called before the set_bnd of the current iteration."
   write(*,*) 
   write(*,*) "Setting Zero-values for nrbc to current values..."
endif
if(lb_dom%nnr_wall > 0) then
   do i=1,lb_dom%nnr_wall   
      x = lb_dom%obs_nrwall(i,1)
      y = lb_dom%obs_nrwall(i,2)
      z = lb_dom%obs_nrwall(i,3)

      ! assign zero values
      lb_dom%nrwall_0val(i,0)      = lb_dom%rho(x,y,z)
      lb_dom%nrwall_0val(i,1:NDIM) = lb_dom%u(NDX(1:NDIM,x,y,z))
   end do
   write(*,*) prc%rk," done."
   write(*,*) 
endif
endif
   end subroutine treat_initial 
  !------------------------------------------------------------------------





end module lb_init     


