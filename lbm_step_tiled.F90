!!
!! Copyright (c) 2010-2011   HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!

#ifdef USE_TILES    /* all the file */
#if defined(USE_OMP_COARSE) && defined(USE_OMP_TASK)
#error "USE_OMP_COARSE and USE_OMP_TASK are mutually exclusive - pick one tile strategy"
#endif
module lbm_step_tiled
  use lbm_functions

#ifdef USE_OMP_COARSE
  use omp_lib
#endif /* USE_OMP_COARSE */

#include "include/replace.h"
contains

  ! TODO:
  !
  ! lb_step() is called from the main time loop
  ! does one lbm step, ie. propagate & collide
  ! and possibly boundary conditions.
  !
  ! The lb_step() is not done on the block lb_dom owned by every MPI
  ! rank, but instead on smaller "tiles" lb_dom_tile as SMPSs tasks
  ! in the subroutine lb_step_tile().
  subroutine lb_step(lb_dom, pnt_nxt, pnt_cur, s_par, prc, meas, tStep)
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_nxt  ! intent(out) might be enough
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_cur  ! intent(out) might be enough
    type(sim_parameter),intent(in)  :: s_par
    type(mpl_var),intent(inout)     :: prc
    type(measure)                   :: meas
    integer, intent(in)             :: tStep

    integer                         :: nTiles, iTile      ! total # of tiles, tile index
    integer                         :: nOuterTiles        ! total # of outer tiles
    integer,dimension(:,:),pointer  :: bounds_tile        ! boundaries of given tile in parent block
    type(lb_block),dimension(:),pointer                 :: lb_dom_tiles    ! lb_blocks for tiles
    integer,dimension(:),pointer    :: sentinel_inner, sentinel_outer           ! sentinel to distinguish between, inner, outer, and all tiles

    !if (prc%rk==0) then
    !   write (*,*) ' .... Doing tiled lb_step!'
    !endif

    !---------------------------
    ! Work on SMPSs tiles rather than on the MPI lb_block

    !  number of tiles
    nTiles = prc%nTiles

    ! tile boundaries within parent block
    bounds_tile => prc%bounds_tile
    ! bounds_tile needs to be pointer, otherwise will be overwritten -> race-condition

    ! container of tile lb_blocks
    lb_dom_tiles => prc%lb_dom_tiles

    sentinel_outer => bounds_tile(:,prc%nOuterTiles)   ! set sentinel to "outer tiles"
    sentinel_inner => bounds_tile(:,0)   ! set sentinel to "inner tiles"

#ifdef USE_OMP_COARSE
!$omp parallel
!$omp do schedule(runtime)
#endif /* USE_OMP_COARSE */
    do iTile = 1, prc%nOuterTiles
       ! select a tile and
       ! do a single tiled Lattice Boltzmann step
#ifdef USE_OMP_TASK
       !$omp task &
       !$omp&     depend(inoutset: lb_dom) depend(inoutset: sentinel_outer)
#endif /* USE_OMP_TASK */
       call lb_step_tile_task(lb_dom_tiles(iTile), s_par, prc, meas, tStep)
#ifdef USE_OMP_TASK
       !$omp end task
#endif /* USE_OMP_TASK */
    end do
#ifdef USE_OMP_COARSE
!$omp end do
#endif /* USE_OMP_COARSE */

    ! the !$omp do above ends with an implicit barrier, so all outer tiles
    ! are complete before the inner tiles start; no extra single region needed

#ifdef USE_OMP_COARSE
!$omp do schedule(runtime)
#endif /* USE_OMP_COARSE */
    do iTile = prc%nOuterTiles+1, nTiles
       ! select a tile and
       ! do a single tiled Lattice Boltzmann step
#ifdef USE_OMP_TASK
       !$omp task &
       !$omp&     depend(inoutset: lb_dom) depend(inoutset: sentinel_inner)
#endif /* USE_OMP_TASK */
       call lb_step_tile_task(lb_dom_tiles(iTile), s_par, prc, meas, tStep)
#ifdef USE_OMP_TASK
       !$omp end task
#endif /* USE_OMP_TASK */
    end do
#ifdef USE_OMP_COARSE
!$omp end do
!$omp end parallel
#endif /* USE_OMP_COARSE */

! appply physical boundary conditions
#ifdef USE_OMPSS_FORK_JOIN
!$OMP TASK INOUT(lb_dom, sentinel_inner, sentinel_outer) PRIORITY(50)  label(set_bnds) no_copy_deps
#else
   ! pass
#endif /* USE_OMPSS_FORK_JOIN */
#ifdef USE_OMP_TASK
!$omp task &
!$omp&     depend(inoutset: lb_dom) depend(inout: sentinel_outer)
#endif /* USE_OMP_TASK */ 
    call set_bnds_task(lb_dom, s_par,prc,meas, tStep)
#ifdef USE_OMP_TASK
!$omp end task
#endif /* USE_OMP_TASK */

  end subroutine lb_step




  ! initialize machinery for tiles implementation
  subroutine init_tile(lb_dom, s_par, prc)
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    type(sim_parameter), intent(inout)  :: s_par
    type(mpl_var),intent(inout)     :: prc

    integer                         :: alloc_stat
    integer                         :: iTile

    ! calculate number of tiles for given size tileSize
    prc%tileSize = TILESIZE
    prc%nTiles = lb_dom%lx(1)/prc%tileSize * lb_dom%lx(2)/prc%tileSize * lb_dom%lx(3)/prc%tileSize

    ! assert compatibility of tile size and domain size
    if ((modulo(lb_dom%lx(1),prc%tileSize)/=0) &
         &.or. (modulo(lb_dom%lx(2),prc%tileSize)/=0) &
         &.or. (modulo(lb_dom%lx(3),prc%tileSize)/=0)) then
       write (*,*) "Fatal error: domain sizes must be and integer multiple of TILESIZE=", prc%tileSize
       write (*,*) "Exiting ..."
       s_par%goend = .true.        ! why is this not sufficient?
       call abort()
    end if

    ! compute tile boundaries within parent block
    call tile_bounds(prc%tileSize, prc%nTiles, lb_dom%lx, prc%n_ghost_nodes, prc%bounds_tile, prc%nOuterTiles)

    ! allocate space for a container of tile lb_blocks
    allocate(prc%lb_dom_tiles(prc%nTiles),stat=alloc_stat)
    if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tiles"

    ! localize data fields of all tile lb_blocks
#ifdef USE_OMP_COARSE
    !$omp parallel
    !$omp do schedule(runtime)
#endif /* USE_OMP_COARSE */
    do iTile = 1, prc%nTiles
       call localize_tile(lb_dom, prc%bounds_tile(:,iTile), prc)
       !write (*,"(Z20 Z20)") loc(lb_dom), loc(prc%bounds_tile(:,iTile))
    end do
#ifdef USE_OMP_COARSE
    !$omp end do
    !$omp end parallel
#endif /* USE_OMP_COARSE */

    ! CSS: do not merge these two loops (above and below) or will serialize
    ! allocate data fields of all tile lb_blocks
    do iTile = 1, prc%nTiles
       call allocate_tile(lb_dom, prc%lb_dom_tiles(iTile), prc%bounds_tile(:,iTile), prc)
    end do
  end subroutine init_tile


  ! gather data from tiles into block for output, etc.
  subroutine gather_tiles(prc, lb_dom)
    implicit none
    type(mpl_var),intent(inout)     :: prc
    type(lb_block),intent(inout)    :: lb_dom

    integer                         :: alloc_stat
    integer                         :: iTile
    integer,dimension(6)            :: tb
    type(lb_block),pointer          :: lb_dom_tile
    integer                         :: lx, ly, lz
    integer                         :: nghost
    integer,dimension(4)            :: lb        ! lower index of physical cells in lb_dom_tile%pnt_cur/nxt
    integer,dimension(4)            :: tS        ! tileSize excluding ghost cells

    nghost = prc%n_ghost_nodes

    do iTile = 1, prc%nTiles
       lb_dom_tile => prc%lb_dom_tiles(iTile)

       ! set tile size from tile boundaries
       tb = prc%bounds_tile(:,iTile)
       lb = lbound(lb_dom_tile%pnt_cur_tile)              ! get lower boundaries from buffer
       tS = shape(lb_dom_tile%pnt_cur_tile)-nghost*2
#ifdef LAYOUT_LIJK
       lb(1:3) = lb(2:4)                         ! spatial indices should be in lb(1:3) regardless of LAYOUT_IJKL/LIJK
       tS(1:3) = tS(2:4)
#endif /* LAYOUT_LIJK */

       ! copy data from tile to block, exclude ghost cells and first physical cell (aka boundary cell) as are invalid
       lb_dom%pnt_nxt(NDX(:,tb(1)+nghost+1:tb(2)-nghost-1,tb(3)+nghost+1:tb(4)-nghost-1,tb(5)+nghost+1:tb(6)-nghost-1)) = &
           & lb_dom_tile%pnt_cur_tile(NDX(:,lb(1)+nghost+1:lb(1)+tS(1)-1,lb(2)+nghost+1:lb(2)+tS(2)-1,lb(3)+nghost+1:lb(3)+tS(3)-1))
    end do

  end subroutine gather_tiles


  ! calculates boundaries of tile with given number i within parent block
  subroutine tile_bounds(ts, nT, lx, ng, tb, nOT)
    implicit none
    integer,intent(in)           :: ts, nT, lx(3), ng        ! tile size, tile number, domain size and # of ghostcells
    integer,intent(out),dimension(:,:),pointer          :: tb               ! tile boundaries
    integer,intent(out)          :: nOT                      ! number of outer tiles

    integer                      :: i,nx,ny,nz,xl,xu,yl,yu,zl,zu,tx,ty,tz
    integer                      :: alloc_stat



    ! allocate memory for tile boundaries tb(:,1..nT) and "all tiles" tb(:,0)
    allocate(tb(6,0:nT), stat=alloc_stat)

    nx = lx(1)/ts
    ny = lx(2)/ts
    nz = lx(3)/ts

    ! set boundaries of whole block
    tb(:,0) = (/ 0, ts*nx+2*ng, 0, ts*ny+2*ng, 0, ts*nz+2*ng /)
    !write (*,*) tb(:,0)

    i = 1
    ! outer tiles first, upper and lower tx
    do tx=1,nx,max(1,nx-1)
       do ty = 1, ny
          do tz = 1, nz
             ! tile boundaries including ghostcells
             xl = (tx-1)*ts; xu = xl + ts+2*ng - 1     ! fields are zero based in spatial coords?
             yl = (ty-1)*ts; yu = yl + ts+2*ng - 1
             zl = (tz-1)*ts; zu = zl + ts+2*ng - 1

             tb(1,i)=xl;tb(2,i)=xu; tb(3,i)=yl;tb(4,i)=yu; tb(5,i)=zl;tb(6,i)=zu
             i = i+1
          end do
       end do
    end do

    ! next upper and lower ty
    do ty=1,ny,max(1,ny-1)
       do tx = 2,nx-1
          do tz = 1, nz
             ! tile boundaries including ghostcells
             xl = (tx-1)*ts; xu = xl + ts+2*ng - 1     ! fields are zero based in spatial coords?
             yl = (ty-1)*ts; yu = yl + ts+2*ng - 1
             zl = (tz-1)*ts; zu = zl + ts+2*ng - 1

             tb(1,i)=xl;tb(2,i)=xu; tb(3,i)=yl;tb(4,i)=yu; tb(5,i)=zl;tb(6,i)=zu
             i = i+1
          end do
       end do
    end do

    ! finally upper and lower tz
    do tz=1,nz,max(1,nz-1)
       do tx = 2,nx-1
          do ty = 2,ny-1
             ! tile boundaries including ghostcells
             xl = (tx-1)*ts; xu = xl + ts+2*ng - 1     ! fields are zero based in spatial coords?
             yl = (ty-1)*ts; yu = yl + ts+2*ng - 1
             zl = (tz-1)*ts; zu = zl + ts+2*ng - 1

             tb(1,i)=xl;tb(2,i)=xu; tb(3,i)=yl;tb(4,i)=yu; tb(5,i)=zl;tb(6,i)=zu
             i = i+1
          end do
       end do
    end do
    nOT = i-1
    !write (*,*) "# outer tiles:", nOT

    ! now inner tiles
    do tx = 2, nx-1
       do ty = 2, ny-1
          do tz = 2, nz-1
             ! tile boundaries including ghostcells
             xl = (tx-1)*ts; xu = xl + ts+2*ng - 1     ! fields are zero based in spatial coords?
             yl = (ty-1)*ts; yu = yl + ts+2*ng - 1
             zl = (tz-1)*ts; zu = zl + ts+2*ng - 1

             tb(1,i)=xl;tb(2,i)=xu; tb(3,i)=yl;tb(4,i)=yu; tb(5,i)=zl;tb(6,i)=zu
             i = i+1
          end do
       end do
    end do
    !write (*,*) tb

  end subroutine tile_bounds

  subroutine lb_step_tile(lb_dom, pnt_nxt, pnt_cur, s_par, prc, meas, tStep)
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_nxt
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_cur
    type(sim_parameter),intent(in)  :: s_par
    type(mpl_var),intent(in)     :: prc
    type(measure),intent(inout)     :: meas
    integer, intent(in)             :: tStep

    real(R8B),pointer,dimension(:,:,:,:)                :: swapper

      !---------------------------
      ! Combined Stream-collide Routine
#ifndef COMBINE_STREAM_COLLIDE
#error "The tiled version of lbc necessarily requires to set -DCOMBINE_STREAM_COLLIDE in the Makefile."
#endif

    call stream_collide_bgk( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
  end subroutine lb_step_tile


  subroutine lb_step_tile_task(lb_dom, s_par, prc, meas, tStep)
    implicit none
    type(lb_block), intent(inout)   :: lb_dom
    type(sim_parameter),intent(in)  :: s_par
    type(mpl_var),intent(in)        :: prc
    type(measure),intent(inout)     :: meas
    integer, intent(in)             :: tStep

    real(R8B),pointer,dimension(:,:,:,:)                 :: dummy_pnt_nxt
    real(R8B),pointer,dimension(:,:,:,:)                 :: dummy_pnt_cur
    real(R8B),pointer,dimension(:,:,:,:)                 :: swapper

    ! unpack pointers to cur and nxt field
    dummy_pnt_cur => lb_dom%pnt_cur_tile
    dummy_pnt_nxt => lb_dom%pnt_nxt_tile

    call lb_step_tile(lb_dom, dummy_pnt_nxt, dummy_pnt_cur, s_par, prc, meas, tStep)

    swapper=>lb_dom%pnt_nxt_tile; lb_dom%pnt_nxt_tile=>lb_dom%pnt_cur_tile; lb_dom%pnt_cur_tile=>swapper

  end subroutine lb_step_tile_task

  subroutine set_bnds_task(lb_dom, s_par, prc, meas, tStep)
    use lbm_functions
    implicit none
    type(lb_block), intent(inout)          :: lb_dom
    type(sim_parameter),intent(in)         :: s_par
    type(mpl_var),intent(in)               :: prc
    type(measure),intent(inout)            :: meas
    integer, intent(in)                    :: tStep

    ! do boundary conditions after stream and collide steps
    call set_nrbc_antibb(lb_dom,lb_dom%pnt_nxt,lb_dom%pnt_cur,s_par,prc,meas)
    call set_bnd(lb_dom,lb_dom%state,lb_dom%pnt_nxt,lb_dom%rho,lb_dom%u,s_par,prc,meas,tStep)

  end subroutine set_bnds_task

  ! allocate memory for tiles and init them
  subroutine allocate_tile(lb_dom, lb_dom_tile, tb, prc)
    use lbm_functions
    implicit none
    type(lb_block),intent(in)       :: lb_dom
    type(lb_block),intent(inout)    :: lb_dom_tile
    integer,dimension(6),intent(inout) :: tb
    type(mpl_var),intent(in)        :: prc

    integer                         :: alloc_stat
    integer                         :: lx, ly, lz
    integer                         :: nghost


    ! set tile size from tile boundaries
    lb_dom_tile%lx = (/ tb(2)-tb(1)-1, tb(4)-tb(3)-1, tb(6)-tb(5)-1 /)
    nghost = prc%n_ghost_nodes
    lx = lb_dom_tile%lx(1); ly = lb_dom_tile%lx(2); lz = lb_dom_tile%lx(3);

    ! map lb_dom subarrays to tiles
    ! %u0 could be pointer if remapping of indices worked. Works as of GNU 4.6
    lb_dom_tile%u0(NDX(1:,0:,0:,0:)) => lb_dom%u0(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%state => lb_dom%state(tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)
    lb_dom_tile%fIn(NDX(1:,0:,0:,0:)) => lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%fOut(NDX(1:,0:,0:,0:)) => lb_dom%fOut(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))

    lb_dom_tile%pnt_cur_tile => lb_dom_tile%fIn
    lb_dom_tile%pnt_nxt_tile => lb_dom_tile%fOut

  end subroutine allocate_tile

  subroutine localize_tile(lb_dom, tb, prc)
    use lbm_functions
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    integer,dimension(6),intent(inout) :: tb    ! inout: serves also as sentinel to track deps, as lb_dom_tile not in args
    type(mpl_var),intent(in)        :: prc

    integer                         :: nghost


    nghost = prc%n_ghost_nodes
    ! do first touch in tasks to localize arrays pages
    lb_dom%u0(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6))) = 0.0
    lb_dom%state(tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost) = 0.0
    lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6))) = 0.0
    lb_dom%fOut(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6))) = 0.0

  end subroutine localize_tile

end module lbm_step_tiled

#endif /* USE_TILES (all the file) */
