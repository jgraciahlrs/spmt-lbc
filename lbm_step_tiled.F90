!!
!! Copyright (c) 2010-2011   HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


module lbm_step_tiled
  use lbm_functions

#include "include/replace.h"
contains

  ! TODO: 
  !   * assert TILESIZE and block dimensions are compatible

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
    real(R8B),pointer,dimension(:,:,:,:)                :: pnt_nxt_tile
    real(R8B),pointer,dimension(:,:,:,:)                :: pnt_cur_tile
    real(R8B),pointer,dimension(:,:,:,:)                :: swapper

#ifdef USENOT_SMPSS_MODULES
    ! need to declare explitic interfaces for StarSs tasks
    interface
       !$CSS TASK
       subroutine decompose_lb_block(lb_dom, tb, nghost, lb_dom_tile)
#ifdef USENOT_SMPSS_MODULES
         use lbm_functions
#endif /* USENOT_SMPSS_MODULES */ 
         implicit none
         type(lb_block),intent(in)    :: lb_dom
         integer,dimension(6),intent(in)  :: tb ! boundaries of given tile in parent block
         integer,intent(in)           :: nghost    ! # of ghost cells
         type(lb_block),intent(inout) :: lb_dom_tile
       end subroutine decompose_lb_block
       !$CSS TASK REDUCTION(lb_dom_fIn)
       subroutine recompose_lb_block_array(lb_dom_tile_fIn, tS1,tS2,tS3, tb, nghost, lx1,lx2,lx3, lb_dom_fIn)
         use lbm_functions
         implicit none
         real(R8B),dimension(NDX(nnod,0:tS1+1,0:tS2+1,0:tS3+1)), intent(in) :: lb_dom_tile_fIn
         integer, intent(in)          :: tS1,tS2,tS3             ! tileSize excluding ghost cells
         integer,dimension(6),intent(in)            :: tb        ! boundaries of given tile in parent block
         integer,intent(in)           :: nghost    ! # of ghost cells
         !integer,dimension(3),intent(in) :: lx
         integer,intent(in)               :: lx1,lx2,lx3
         !real(R8B),dimension(NDX(nnod,0:lx(1)+1,0:lx(2)+1,0:lx(3)+1)),intent(inout) :: lb_dom_fIn
         !real(R8B),dimension(NDX(nnod,0:65,0:65,0:65)),intent(inout) :: lb_dom_fIn
         real(R8B),dimension(NDX(nnod,0:lx1+1,0:lx2+1,0:lx3+1)),intent(inout) :: lb_dom_fIn
       end subroutine recompose_lb_block_array
       !$CSS TASK REDUCTION(meas)
       subroutine lb_step_tile_array(lx1,lx2,lx3, state, u0, nxt, cur, s_par, prc, meas, tStep)
         use lbm_functions
         implicit none
         !integer,dimension(3),intent(in) :: lx
         integer,intent(in)               :: lx1,lx2,lx3
         integer,target,dimension(1:lx1,1:lx2,1:lx3),intent(in)        :: state
         real(R8B),target,dimension(NDX(NDIM,0:lx1+1,0:lx2+1,0:lx3+1)),intent(in)     :: u0
         real(R8B),target,dimension(NDX(nnod,0:lx1+1,0:lx2+1,0:lx3+1)), intent(inout) :: nxt
         real(R8B),target,dimension(NDX(nnod,0:lx1+1,0:lx2+1,0:lx3+1)), intent(inout) :: cur
         type(sim_parameter),intent(in)  :: s_par
         type(mpl_var),intent(in)        :: prc
         type(measure),intent(inout)         :: meas   ! this is cheating as it really is inout, wait for reduction
         integer, intent(in)             :: tStep
       end subroutine lb_step_tile_array
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
       call decompose_lb_block(lb_dom, bounds_tile(:,iTile), nghost, lb_dom_tiles(iTile))
    end do
    ! $CSS BARRIER

    do iTile = 1, nTiles
       ! select tile
       lb_dom_tile => lb_dom_tiles(iTile)
       ! $CSS WAIT ON(lb_dom_tile)
       !$CSS WAIT ON(bounds_tile(:,iTile))

       pnt_cur_tile => lb_dom_tile%fIn
       pnt_nxt_tile => lb_dom_tile%fOut

       ! do a single tiled Lattice Boltzmann step
#ifndef SMPSS_NO_ARRAYS /* version with arrays */
       call lb_step_tile_array(lb_dom_tile%lx(1), lb_dom_tile%lx(2), lb_dom_tile%lx(3), &
            &lb_dom_tile%state, lb_dom_tile%u0, pnt_nxt_tile, pnt_cur_tile, s_par, prc, meas, tStep)
#ifdef COMBINE_STREAM_COLLIDE
       swapper => pnt_nxt_tile; pnt_nxt_tile => pnt_cur_tile; pnt_cur_tile => swapper
#endif
       lb_dom_tile%fIn => pnt_cur_tile; lb_dom_tile%fOut => pnt_nxt_tile    ! fix pointers, as above works with arrays

       ! reassemble block from tiles
       call recompose_lb_block_array(lb_dom_tile%fIn, tileSize,tileSize,tileSize, &
            & bounds_tile(:,iTile), nghost,&
            & lb_dom%lx(1),lb_dom%lx(2),lb_dom%lx(3), lb_dom%fBuffer)
#else     /* version without arrays */
       call lb_step_tile(lb_dom_tile, pnt_nxt_tile, pnt_cur_tile, s_par, prc, meas, tStep)
#ifdef COMBINE_STREAM_COLLIDE
       swapper => lb_dom_tile%fIn; lb_dom_tile%fIn => lb_dom_tile%fOut; lb_dom_tile%fOut => swapper
#endif
       call recompose_lb_block(lb_dom_tile, bounds_tile(:,iTile), nghost, lb_dom)
#endif
    end do

    !$CSS BARRIER

    ! copy boundaries from buffer to block and swap pointers
    lb_dom%fBuffer(NDX(:,0,:,:)) = lb_dom%fIn(NDX(:,0,:,:))
    lb_dom%fBuffer(NDX(:,:,0,:)) = lb_dom%fIn(NDX(:,:,0,:))
    lb_dom%fBuffer(NDX(:,:,:,0)) = lb_dom%fIn(NDX(:,:,:,0)) 
    lb_dom%fBuffer(NDX(:,lb_dom%lx(1)+1,:,:)) = lb_dom%fIn(NDX(:,lb_dom%lx(1)+1,:,:))
    lb_dom%fBuffer(NDX(:,:,lb_dom%lx(2)+1,:)) = lb_dom%fIn(NDX(:,:,lb_dom%lx(2)+1,:))
    lb_dom%fBuffer(NDX(:,:,:,lb_dom%lx(3)+1)) = lb_dom%fIn(NDX(:,:,:,lb_dom%lx(3)+1))
    swapper => lb_dom%fIn; lb_dom%fIn => lb_dom%fBuffer; lb_dom%fBuffer => swapper 


    ! set cur/nxt pointers on block 
    pnt_cur => lb_dom%fOut
    pnt_nxt => lb_dom%fIn

    ! do boundary conditions after stream and collide steps
    call set_nrbc_antibb(lb_dom,pnt_nxt,pnt_cur,s_par,prc,meas)
    call set_bnd(lb_dom,lb_dom%state,pnt_nxt,lb_dom%rho,lb_dom%u,s_par,prc,meas,tStep)
#endif  /* USE_TILES*/

  end subroutine lb_step


#ifdef USE_TILES

  ! initialize machinery for tiles implementation
  subroutine init_tile(prc, lb_dom)
    implicit none
    type(mpl_var),intent(inout)     :: prc
    type(lb_block),intent(inout)    :: lb_dom

    integer                         :: alloc_stat
    integer                         :: iTile
    integer,dimension(6)            :: tb
    type(lb_block),pointer          :: lb_dom_tile
    integer                         :: lx, ly, lz
    integer                         :: nghost

    ! calculate number of tiles for given size tileSize
    prc%tileSize = TILESIZE
    prc%nTiles = lb_dom%lx(1)/prc%tileSize * lb_dom%lx(2)/prc%tileSize * lb_dom%lx(3)/prc%tileSize
    
    ! compute tile boundaries within parent block 
    call tile_bounds(prc%tileSize, prc%nTiles, prc%n_ghost_nodes, prc%bounds_tile)
    
    ! allocate space for a container of tile lb_blocks
    allocate(prc%lb_dom_tiles(prc%nTiles),stat=alloc_stat)
    if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tiles"
    
    ! allocate data fields of all tile lb_blocks
    do iTile = 1, prc%nTiles
       lb_dom_tile => prc%lb_dom_tiles(iTile)

       ! set tile size from tile boundaries
       tb = prc%bounds_tile(:,iTile)
       lb_dom_tile%lx = (/ tb(2)-tb(1)-1, tb(4)-tb(3)-1, tb(6)-tb(5)-1 /)
       nghost = prc%n_ghost_nodes

       ! allocate memory for tile arrays
       lx = lb_dom_tile%lx(1); ly = lb_dom_tile%lx(2); lz = lb_dom_tile%lx(3); 
       allocate(lb_dom_tile%fOut(NDX(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
           if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fOut"
       allocate(lb_dom_tile%fIn(NDX(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
           if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fOut"
       allocate(lb_dom_tile%state(1:lx,1:ly,1:lz),stat=alloc_stat)
           if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fOut"
       allocate(lb_dom_tile%u0(NDX(NDIM,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
           if(alloc_stat /=0) write(*,*) "Error allocating lb_dom_tile%fOut"

       lb_dom_tile%state = lb_dom%state(tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)
       lb_dom_tile%u0 = lb_dom%u0(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))
    end do

  end subroutine init_tile
 
  ! calculates boundaries of tile with given number i within parent block    
  subroutine tile_bounds(ts, nT, ng, tb)
    implicit none
    integer,intent(in)           :: ts, nT, ng       ! tile size, and # of ghostcells  
    integer,intent(out),dimension(:,:),pointer          :: tb               ! tile boundaries   
    
    integer                      :: i,n,xl,xu,yl,yu,zl,zu,tx,ty,tz
    integer                      :: alloc_stat

    ! allocate memory for tile boundaries
    allocate(tb(6,nT), stat=alloc_stat)

    n = nT**(1./3.)                              ! number of tiles in each direction

    do i = 1, nT                         ! loop over all tiles
       tx = mod( i-1, n) + 1                     ! (i-1) % ntiles + 1
       ty = mod( (i-1)/n, n) + 1                 ! ( (i-1)//ntiles ) % ntiles + 1
       tz = (i-1) / (n*n) + 1                    ! (i-1) // (ntiles*ntiles) + 1
       !write (*,*) '...tilebounds i,tx,ty,tz: ', i, tx,ty,tz

       ! tile boundaries including ghostcells
       xl = (tx-1)*ts; xu = xl + ts+2*ng - 1     ! fields are zero based in spatial coords?
       yl = (ty-1)*ts; yu = yl + ts+2*ng - 1
       zl = (tz-1)*ts; zu = zl + ts+2*ng - 1
       !write (*,*) '...tilebounds i,xl,yl,zl: ', i, xl,yl,zl

       tb(1,i)=xl;tb(2,i)=xu; tb(3,i)=yl;tb(4,i)=yu; tb(5,i)=zl;tb(6,i)=zu
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

    ! copy from block to tile; ignore ghost cell
    lb_dom%fBuffer(NDX(:,tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)) = &
         & lb_dom_tile%fIn(NDX(:,lb(1)+1:lb(1)+tS(1),lb(2)+1:lb(2)+tS(2),lb(3)+1:lb(3)+tS(3))) ! assumes 1 ghost cell

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

  !$CSS TASK
  subroutine decompose_lb_block(lb_dom, tb, nghost, lb_dom_tile)
#ifdef USENOT_SMPSS_MODULES
    use lbm_functions
#endif /* USENOT_SMPSS_MODULES */ 
    implicit none
    type(lb_block),intent(in)    :: lb_dom
    integer,dimension(6),intent(in)  :: tb ! boundaries of given tile in parent block
    integer,intent(in)           :: nghost    ! # of ghost cells
    type(lb_block),intent(inout) :: lb_dom_tile

    integer                      :: lx,ly,lz  ! block sizes
    integer                      :: alloc_stat

    ! copy from block to tile
    lb_dom_tile%fIN   = lb_dom%fIN(NDX(:,tb(1):tb(2),tb(3):tb(4),tb(5):tb(6)))

    !write (*,*) '.... Done decomposing block!'

  end subroutine decompose_lb_block


  !$CSS TASK REDUCTION(lb_dom_fIn)
  subroutine recompose_lb_block_array(lb_dom_tile_fIn, tS1,tS2,tS3, tb, nghost, lx1,lx2,lx3, lb_dom_fIn)
#ifdef USENOT_SMPSS_MODULES
    use lbm_functions
#endif /* USENOT_SMPSS_MODULES */
    implicit none
    real(R8B),dimension(NDX(nnod,0:tS1+1,0:tS2+1,0:tS3+1)), intent(in) :: lb_dom_tile_fIn
    integer, intent(in)          :: tS1,tS2,tS3            ! tileSize excluding ghost cells
    integer,dimension(6),intent(in) :: tb  ! boundaries of given tile in parent block
    integer,intent(in)           :: nghost    ! # of ghost cells
    !integer,dimension(3),intent(in) :: lx
    integer,intent(in)               :: lx1,lx2,lx3
    !real(R8B),dimension(NDX(nnod,0:lx(1)+1,0:lx(2)+1,0:lx(3)+1)),intent(inout) :: lb_dom_fIn
    !real(R8B),dimension(NDX(nnod,0:65,0:65,0:65)),intent(inout) :: lb_dom_fIn
    real(R8B),dimension(NDX(nnod,0:lx1+1,0:lx2+1,0:lx3+1)),intent(inout) :: lb_dom_fIn

    integer,dimension(4)         :: lb        ! lower index of physical cells in lb_dom_tile%fIn / %fOut
    integer,dimension(4)         :: tS        ! tileSize excluding ghost cells

    lb = lbound(lb_dom_tile_fIn)              ! get lower boundaries from buffer
    tS = shape(lb_dom_tile_fIn)-nghost*2
#ifdef LAYOUT_LIJK
    lb(1:3) = lb(2:4)                         ! spatial indices should be in lb(1:3) regardless of LAYOUT_IJKL/LIJK
    tS(1:3) = tS(2:4)
#endif
    

    ! copy from block to tile; ignore ghost cell
    lb_dom_fIn(NDX(:,tb(1)+nghost:tb(2)-nghost,tb(3)+nghost:tb(4)-nghost,tb(5)+nghost:tb(6)-nghost)) = &
         & lb_dom_tile_fIn(NDX(:,lb(1)+1:lb(1)+tS(1),lb(2)+1:lb(2)+tS(2),lb(3)+1:lb(3)+tS(3))) ! assumes one ghost cell


    !write (*,*) '.... Done recomposing block!'

  end subroutine recompose_lb_block_array

  !$CSS TASK REDUCTION(meas)
  subroutine lb_step_tile_array(lx1,lx2,lx3, state, u0, nxt, cur, s_par, prc, meas, tStep)
#ifdef USENOT_SMPSS_MODULES
    use lbm_functions
    use lbm_step_tiled
#endif /* USENOT_SMPSS_MODULES */
    implicit none
    !integer,dimension(3),intent(in) :: lx
    integer,intent(in) :: lx1,lx2,lx3
    integer,target,dimension(1:lx1,1:lx2,1:lx3),intent(in)        :: state
    real(R8B),target,dimension(NDX(NDIM,0:lx1+1,0:lx2+1,0:lx3+1)),intent(in)     :: u0
    real(R8B),target,dimension(NDX(nnod,0:lx1+1,0:lx2+1,0:lx3+1)), intent(inout) :: nxt
    real(R8B),target,dimension(NDX(nnod,0:lx1+1,0:lx2+1,0:lx3+1)), intent(inout) :: cur
    type(sim_parameter),intent(in)  :: s_par
    type(mpl_var),intent(in)        :: prc
    type(measure),intent(inout)     :: meas
    integer, intent(in)             :: tStep

    type(lb_block)                  :: dummy_lb_dom
    real(R8B),pointer,dimension(:,:,:,:)                 :: dummy_pnt_nxt
    real(R8B),pointer,dimension(:,:,:,:)                 :: dummy_pnt_cur
    
    integer, dimension(3)           :: lx
    lx = (/ lx1,lx2,lx3 /)


    ! convert from array to type inteface
    dummy_pnt_nxt => nxt
    dummy_pnt_cur => cur
    dummy_lb_dom%lx = lx
    dummy_lb_dom%fIn => dummy_pnt_cur
    dummy_lb_dom%fOut => dummy_pnt_nxt
    dummy_lb_dom%state => state
    dummy_lb_dom%u0 => u0

    call lb_step_tile(dummy_lb_dom, dummy_pnt_nxt, dummy_pnt_cur, s_par, prc, meas, tStep)

  end subroutine lb_step_tile_array

#ifndef USENOT_SMPSS_MODULES
end module lbm_step_tiled
#endif /* !USENOT_SMPSS_MODULES */

