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
module lbm_step_tiled
  use lbm_functions

#include "include/replace.h"
contains

  ! TODO: 
  !   
  ! lb_step() is called from the main time loop
  ! does one lbm step, ie. propagate & collide
  ! and possibly boundary conditions.
  !
  ! The block lb_dom owned by every MPI rank is first decomposed
  ! into smaller "tiles" lb_dom_tile, which are then passed on to
  ! lb_step_tile() do to th actual calculations as SMPSs tasks.
  subroutine lb_step(lb_dom, pnt_nxt, pnt_cur, s_par, prc, meas, tStep)
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_nxt  ! intent(out) might be enough
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_cur  ! intent(out) might be enough
    type(sim_parameter),intent(in)  :: s_par
    type(mpl_var),intent(inout)     :: prc
    type(measure)                   :: meas
    integer, intent(in)             :: tStep
    
#ifdef USE_TILES
    integer                         :: nTiles, iTile      ! total # of tiles, tile index
    integer,dimension(:,:),pointer  :: bounds_tile        ! boundaries of given tile in parent block
    type(lb_block),dimension(:),pointer                 :: lb_dom_tiles    ! lb_block's for tiles
    type(lb_block),pointer          :: lb_dom_tile         
    integer                         :: alloc_stat
    integer                         :: tileSize = TILESIZE, nghost

#ifdef USENOT_SMPSS_MODULES
    ! need to declare explitic interfaces for StarSs tasks
    interface
       !$CSS TASK REDUCTION(lb_dom, tb0)
       subroutine decompose_lb_block(lb_dom, tb0, nghost, lb_dom_tile, tb)
#ifdef USENOT_SMPSS_MODULES
         use lbm_functions
#endif /* USENOT_SMPSS_MODULES */ 
         implicit none
         type(lb_block),intent(inout)           :: lb_dom
         integer,dimension(6),intent(inout)     :: tb0              ! sentinel to track dependencies
         integer,intent(in)                  :: nghost           ! # of ghost cells
         type(lb_block),intent(inout)        :: lb_dom_tile
         integer,dimension(6),intent(inout)  :: tb               ! boundaries of given tile in parent block
       end subroutine decompose_lb_block
       !$CSS TASK REDUCTION(lb_dom, tb0)
       subroutine recompose_lb_block_sentinel(lb_dom_tile,tb, nghost, lb_dom,tb0)
         use lbm_functions
         implicit none
         type(lb_block),intent(inout)   :: lb_dom_tile           ! why is intent(in) not sufficient?
         integer,dimension(6),intent(inout)      :: tb           ! boundaries of given tile in parent block
         integer,intent(in)             :: nghost                ! # of ghost cells
         type(lb_block),intent(inout)   :: lb_dom
         integer,dimension(6),intent(inout)      :: tb0          ! sentinel for dependency tracking
       end subroutine recompose_lb_block_sentinel
       !$CSS TASK REDUCTION(meas)
       subroutine lb_step_tile_sentinel(lb_dom, tb, s_par, prc, meas, tStep)
         use lbm_functions
         implicit none
         type(lb_block),intent(inout)        :: lb_dom
         integer,dimension(6),intent(inout)  :: tb        ! sentinel for dependency tracking
         type(sim_parameter),intent(in)      :: s_par
         type(mpl_var),intent(in)            :: prc
         type(measure),intent(inout)         :: meas
         integer, intent(in)                 :: tStep
       end subroutine lb_step_tile_sentinel
       !$CSS TASK
       subroutine swap_and_set_bnd_sentinel(lb_dom, tb0, s_par, prc, meas, tStep)
         use lbm_functions
         implicit none
         type(lb_block), intent(inout)          :: lb_dom
         integer,dimension(6),intent(inout)     :: tb0       ! sentinel for dependency tracking
         type(sim_parameter),intent(in)         :: s_par
         type(mpl_var),intent(in)               :: prc
         type(measure),intent(inout)            :: meas
         integer, intent(in)                    :: tStep
       end subroutine swap_and_set_bnd_sentinel
    end interface
#endif  /* USENOT_SMPSS_MODULES */

    !if (prc%rk==0) then
    !   write (*,*) ' .... Doing tiled lb_step!'
    !endif

    !---------------------------
    ! Make small tasks for SMPSs by decomposing MPI lb_block into smaller tiles

    !  number of ghost cells and tiles
    nghost = prc%n_ghost_nodes
    nTiles = prc%nTiles

    ! tile boundaries within parent block 
    bounds_tile => prc%bounds_tile 
    ! bounds_tile needs to be pointer, otherwise will be overwritten -> race-condition

    ! container of tile lb_blocks 
    lb_dom_tiles => prc%lb_dom_tiles

    do iTile = 1, nTiles           ! decompose block into tiles
       call decompose_lb_block(lb_dom, bounds_tile(:,0), nghost, lb_dom_tiles(iTile), bounds_tile(:,iTile))
    end do

    do iTile = 1, nTiles
       ! select tile
       lb_dom_tile => lb_dom_tiles(iTile)

       ! do a single tiled Lattice Boltzmann step
       call lb_step_tile_sentinel(lb_dom_tile, bounds_tile(:,iTile), s_par, prc, meas, tStep)

       ! reassemble block from tiles
       call recompose_lb_block_sentinel(lb_dom_tile,bounds_tile(:,iTile), nghost,&
            & lb_dom, bounds_tile(:,0))
    end do

    call swap_and_set_bnd_sentinel(lb_dom, bounds_tile(:,0), s_par,prc,meas, tStep)
    !$CSS WAIT ON(bounds_tile(:,0))

#endif  /* USE_TILES*/

  end subroutine lb_step


#ifdef USE_TILES


  ! initialize machinery for tiles implementation
  subroutine init_tile(lb_dom, s_par, prc)
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    type(sim_parameter), intent(inout)  :: s_par 
    type(mpl_var),intent(inout)     :: prc

#ifdef USENOT_SMPSS_MODULES
    ! need to declare explitic interfaces for StarSs tasks
    interface
       !$CSS TASK 
       subroutine allocate_tile(lb_dom, lb_dom_tile, tb, prc)
         use lbm_functions
         implicit none
         type(lb_block),intent(in)       :: lb_dom
         type(lb_block),intent(inout)    :: lb_dom_tile
         integer,dimension(6),intent(in) :: tb
         type(mpl_var),intent(in)        :: prc
       end subroutine allocate_tile
    end interface
#endif /* USENOT_SMPSS_MODULES */

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
    !call tile_bounds(prc%tileSize, prc%nTiles, prc%n_ghost_nodes, prc%bounds_tile)
    call tile_bounds(prc%tileSize, prc%nTiles, lb_dom%lx, prc%n_ghost_nodes, prc%bounds_tile)
    
    ! allocate space for a container of tile lb_blocks
    allocate(prc%lb_dom_tiles(prc%nTiles),stat=alloc_stat)
    if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tiles"
    
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
    integer,dimension(4)            :: lb        ! lower index of physical cells in lb_dom_tile%fIn / %fOut
    integer,dimension(4)            :: tS        ! tileSize excluding ghost cells

    nghost = prc%n_ghost_nodes

    do iTile = 1, prc%nTiles
       lb_dom_tile => prc%lb_dom_tiles(iTile)

       ! set tile size from tile boundaries
       tb = prc%bounds_tile(:,iTile)
       lb = lbound(lb_dom_tile%fIn)              ! get lower boundaries from buffer
       tS = shape(lb_dom_tile%fIn)-nghost*2
#ifdef LAYOUT_LIJK
       lb(1:3) = lb(2:4)                         ! spatial indices should be in lb(1:3) regardless of LAYOUT_IJKL/LIJK
       tS(1:3) = tS(2:4)
#endif

       ! copy data from tile to block, exclude ghost cells and first physical cell (aka boundary cell) as are invalid
       lb_dom%fIn(NDX(:,tb(1)+nghost+1:tb(2)-nghost-1,tb(3)+nghost+1:tb(4)-nghost-1,tb(5)+nghost+1:tb(6)-nghost-1)) = &
           & lb_dom_tile%fIn(NDX(:,lb(1)+nghost+1:lb(1)+tS(1)-1,lb(2)+nghost+1:lb(2)+tS(2)-1,lb(3)+nghost+1:lb(3)+tS(3)-1))
    end do

  end subroutine gather_tiles

 
  ! calculates boundaries of tile with given number i within parent block    
  subroutine tile_bounds(ts, nT, lx, ng, tb)
    implicit none
    integer,intent(in)           :: ts, nT, lx(3), ng        ! tile size, tile number, domain size and # of ghostcells  
    integer,intent(out),dimension(:,:),pointer          :: tb               ! tile boundaries   
    
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
    do tx = 1, nx
       do ty = 1, ny
          do tz = 1, nz
             ! tile boundaries including ghostcells
             xl = (tx-1)*ts; xu = xl + ts+2*ng - 1     ! fields are zero based in spatial coords?
             yl = (ty-1)*ts; yu = yl + ts+2*ng - 1
             zl = (tz-1)*ts; zu = zl + ts+2*ng - 1
             !write (*,*) '...tilebounds tx,ty,tz,xl,yl,zl: ', tx,ty,tz, xl,yl,zl

             tb(1,i)=xl;tb(2,i)=xu; tb(3,i)=yl;tb(4,i)=yu; tb(5,i)=zl;tb(6,i)=zu
             i = i+1
          end do
       end do
    end do
    !write (*,*) tb

  end subroutine tile_bounds

  subroutine recompose_lb_block(lb_dom_tile, tb, nghost, lb_dom)
    implicit none
    type(lb_block),intent(in)    :: lb_dom_tile
    integer,dimension(6),intent(in)           :: tb        ! boundaries of given tile in parent block
    integer,intent(in)           :: nghost    ! # of ghost cells
    type(lb_block),intent(inout) :: lb_dom

    integer,dimension(4)         :: lb        ! lower index of physical cells in lb_dom_tile%fIn / %fOut
    integer,dimension(4)         :: tS        ! tileSize excluding ghost cells

    lb = lbound(lb_dom_tile%fIn)              ! get lower boundaries from buffer
    tS = shape(lb_dom_tile%fIn)-nghost*2
#ifdef LAYOUT_LIJK
    lb(1:3) = lb(2:4)                         ! spatial indices should be in lb(1:3) regardless of LAYOUT_IJKL/LIJK
    tS(1:3) = tS(2:4)
#endif


    ! copy from tile to block; ignore ghost cell
    !lb_dom%fOut(NDX(:,tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)) = &
    !     & lb_dom_tile%fIn(NDX(:,lb(1)+1:lb(1)+tS(1),lb(2)+1:lb(2)+tS(2),lb(3)+1:lb(3)+tS(3))) ! assumes one ghost cell

    ! copy first physical cell from tile to block, assumes one ghost cell
    lb_dom%fOut(NDX(:,tb(1)+nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)) = &
         & lb_dom_tile%fIn(NDX(:,lb(1)+1,lb(2)+1:lb(2)+tS(2),lb(3)+1:lb(3)+tS(3)))
    lb_dom%fOut(NDX(:,tb(1)+nghost:tb(2)-nghost,tb(3)+nghost,tb(5)+nghost:tb(6)-nghost)) = &
         & lb_dom_tile%fIn(NDX(:,lb(1)+1:lb(1)+tS(1),lb(2)+1,lb(3)+1:lb(3)+tS(3)))
    lb_dom%fOut(NDX(:,tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost)) = &
         & lb_dom_tile%fIn(NDX(:,lb(1)+1:lb(1)+tS(1),lb(2)+1:lb(2)+tS(2),lb(3)+1))
    lb_dom%fOut(NDX(:,tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)) = &
         & lb_dom_tile%fIn(NDX(:,lb(1)+tS(1),lb(2)+1:lb(2)+tS(2),lb(3)+1:lb(3)+tS(3)))
    lb_dom%fOut(NDX(:,tb(1)+nghost:tb(2)-nghost,tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)) = &
         & lb_dom_tile%fIn(NDX(:,lb(1)+1:lb(1)+tS(1),lb(2)+tS(2),lb(3)+1:lb(3)+tS(3)))
    lb_dom%fOut(NDX(:,tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(6)-nghost)) = &
         & lb_dom_tile%fIn(NDX(:,lb(1)+1:lb(1)+tS(1),lb(2)+1:lb(2)+tS(2),lb(3)+tS(3)))

    !write (*,*) '.... Done recomposing block!'

  end subroutine recompose_lb_block
#endif  /* USE_TILES*/

  subroutine lb_step_tile(lb_dom, pnt_nxt, pnt_cur, s_par, prc, meas, tStep)
    implicit none
    type(lb_block),intent(inout)    :: lb_dom
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_nxt
    real(R8B),pointer,dimension(:,:,:,:), intent(inout) :: pnt_cur
    type(sim_parameter),intent(in)  :: s_par
    type(mpl_var),intent(in)     :: prc
    type(measure),intent(inout)     :: meas
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
  end subroutine lb_step_tile

#ifdef USENOT_SMPSS_MODULES
end module lbm_step_tiled
#endif /* USENOT_SMPSS_MODULES */

  !$CSS TASK REDUCTION(lb_dom, tb0)
  subroutine decompose_lb_block(lb_dom,tb0, nghost, lb_dom_tile, tb)
#ifdef USENOT_SMPSS_MODULES
    use lbm_functions
#endif /* USENOT_SMPSS_MODULES */ 
    implicit none
    type(lb_block),intent(inout)           :: lb_dom
    integer,dimension(6),intent(inout)     :: tb0            ! sentinel to track dependencies
    integer,intent(in)                  :: nghost         ! # of ghost cells
    type(lb_block),intent(inout)        :: lb_dom_tile
    integer,dimension(6),intent(inout)  :: tb             ! boundaries of given tile in parent block
    
    integer                             :: lx,ly,lz       ! block sizes
    integer                             :: alloc_stat

    lx = lb_dom_tile%lx(1); ly = lb_dom_tile%lx(2); lz = lb_dom_tile%lx(3); 

    ! copy from block to tile
    !lb_dom_tile%fIn(NDX(:,0:lx+1,0:ly+1,0:lz+1)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))

    ! copy ghost cells from block to tile, assumes 1 ghost cell only
    lb_dom_tile%fIn(NDX(:,0,:,:)) = lb_dom%fIn(NDX(:,tb(1),tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,0,:)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3),tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,:,0)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5)))
    lb_dom_tile%fIn(NDX(:,lx+1,:,:)) = lb_dom%fIn(NDX(:,tb(2),tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,ly+1,:)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(4),tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,:,lz+1)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(6)))

    ! copy first physical cell too
    lb_dom_tile%fIn(NDX(:,1,:,:)) = lb_dom%fIn(NDX(:,tb(1)+1,tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,1,:)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3)+1,tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,:,1)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5)+1))
    lb_dom_tile%fIn(NDX(:,lx,:,:)) = lb_dom%fIn(NDX(:,tb(2)-1,tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,ly,:)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(4)-1,tb(5):tb(6)))
    lb_dom_tile%fIn(NDX(:,:,:,lz)) = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(6)-1))
 
    !write (*,*) '.... Done decomposing block!'

  end subroutine decompose_lb_block

  !$CSS TASK REDUCTION(lb_dom, tb0)
  subroutine recompose_lb_block_sentinel(lb_dom_tile,tb, nghost, lb_dom,tb0)
#ifdef USENOT_SMPSS_MODULES
    use lbm_functions
    use lbm_step_tiled
#endif /* USENOT_SMPSS_MODULES */
    implicit none
    type(lb_block),intent(inout)   :: lb_dom_tile          ! why is intent(in) not sufficient?
    integer,dimension(6),intent(inout)      :: tb          ! boundaries of given tile in parent block
    integer,intent(in)             :: nghost               ! # of ghost cells
    type(lb_block),intent(inout)   :: lb_dom
    integer,dimension(6),intent(inout)      :: tb0         ! sentinel for dependency tracking

    call recompose_lb_block(lb_dom_tile,tb, nghost, lb_dom)

  end subroutine recompose_lb_block_sentinel

  !$CSS TASK REDUCTION(meas)
  subroutine lb_step_tile_sentinel(lb_dom, tb, s_par, prc, meas, tStep)
#ifdef USENOT_SMPSS_MODULES
    use lbm_functions
    use lbm_step_tiled
#endif /* USENOT_SMPSS_MODULES */
    implicit none
    type(lb_block), intent(inout)   :: lb_dom
    integer,dimension(6),intent(inout)  :: tb        ! sentinel for dependency tracking
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

#ifdef COMBINE_STREAM_COLLIDE
       swapper=>lb_dom%pnt_nxt_tile; lb_dom%pnt_nxt_tile=>lb_dom%pnt_cur_tile; lb_dom%pnt_cur_tile=>swapper
#endif

  end subroutine lb_step_tile_sentinel

!$CSS TASK
  subroutine swap_and_set_bnd_sentinel(lb_dom, tb0, s_par, prc, meas, tStep)
    use lbm_functions
    implicit none
    type(lb_block), intent(inout)          :: lb_dom
    integer,dimension(6),intent(inout)     :: tb0       ! sentinel for dependency tracking
    type(sim_parameter),intent(in)         :: s_par
    type(mpl_var),intent(in)               :: prc
    type(measure),intent(inout)            :: meas
    integer, intent(in)                    :: tStep
    
    real(R8B),pointer,dimension(:,:,:,:)                :: swapper, pnt_nxt, pnt_cur
    
    ! copy boundaries from buffer to block and swap pointers
    lb_dom%fOut(NDX(:,0,:,:)) = lb_dom%fIn(NDX(:,0,:,:))
    lb_dom%fOut(NDX(:,:,0,:)) = lb_dom%fIn(NDX(:,:,0,:))
    lb_dom%fOut(NDX(:,:,:,0)) = lb_dom%fIn(NDX(:,:,:,0)) 
    lb_dom%fOut(NDX(:,lb_dom%lx(1)+1,:,:)) = lb_dom%fIn(NDX(:,lb_dom%lx(1)+1,:,:))
    lb_dom%fOut(NDX(:,:,lb_dom%lx(2)+1,:)) = lb_dom%fIn(NDX(:,:,lb_dom%lx(2)+1,:))
    lb_dom%fOut(NDX(:,:,:,lb_dom%lx(3)+1)) = lb_dom%fIn(NDX(:,:,:,lb_dom%lx(3)+1))
    swapper => lb_dom%fIn; lb_dom%fIn => lb_dom%fOut; lb_dom%fOut => swapper 
    
    ! set cur/nxt pointers on block 
    pnt_cur => lb_dom%fOut
    pnt_nxt => lb_dom%fIn
    
    ! do boundary conditions after stream and collide steps
    call set_nrbc_antibb(lb_dom,pnt_nxt,pnt_cur,s_par,prc,meas)
    call set_bnd(lb_dom,lb_dom%state,pnt_nxt,lb_dom%rho,lb_dom%u,s_par,prc,meas,tStep)

  end subroutine swap_and_set_bnd_sentinel

  ! allocate memory for tiles and init them
  !$CSS TASK 
  subroutine allocate_tile(lb_dom, lb_dom_tile, tb, prc)
    use lbm_functions
    implicit none
    type(lb_block),intent(in)       :: lb_dom
    type(lb_block),intent(inout)    :: lb_dom_tile
    integer,dimension(6),intent(in) :: tb
    type(mpl_var),intent(in)        :: prc
    
    integer                         :: alloc_stat
    integer                         :: lx, ly, lz
    integer                         :: nghost
    
    
    ! set tile size from tile boundaries
    lb_dom_tile%lx = (/ tb(2)-tb(1)-1, tb(4)-tb(3)-1, tb(6)-tb(5)-1 /)
    nghost = prc%n_ghost_nodes
    
    ! allocate memory for tile arrays
    lx = lb_dom_tile%lx(1); ly = lb_dom_tile%lx(2); lz = lb_dom_tile%lx(3); 
    allocate(lb_dom_tile%fOut(NDX(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
        if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fOut"
    allocate(lb_dom_tile%fIn(NDX(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
        if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fIn"
    allocate(lb_dom_tile%state(1:lx,1:ly,1:lz),stat=alloc_stat)
        if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fOut"
    allocate(lb_dom_tile%u0(NDX(NDIM,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
        if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fOut"
    
    ! %u0 could be pointer if remapping of indices worked
    lb_dom_tile%state = lb_dom%state(tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)
    lb_dom_tile%u0 = lb_dom%u0(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%fIn = lb_dom%fIn(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))
    lb_dom_tile%fOut = 0.  ! first touch

    lb_dom_tile%pnt_cur_tile => lb_dom_tile%fIn
    lb_dom_tile%pnt_nxt_tile => lb_dom_tile%fOut
    
  end subroutine allocate_tile


#ifndef USENOT_SMPSS_MODULES
end module lbm_step_tiled
#endif /* !USENOT_SMPSS_MODULES */

#endif /* USE_TILES (all the file) */ 
