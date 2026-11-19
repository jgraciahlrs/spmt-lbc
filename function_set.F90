!!
!! Copyright (c) 2010        HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


module function_set
   use lbmodel
   use tools
   use mpl_set
   implicit none
#include "include/replace.h"
contains

#ifdef DEBUG
subroutine dump_matrix(matrix,lb_dom,prc,s_par,description,tstep)
!------------------------------------------------------------------------
!
! Master thread gathers the output from all threads in order to write it out.
!


   type(lb_block),intent(in)       :: lb_dom
   real(R8B) :: matrix(NDX(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))
   type(sim_parameter),intent(in)  :: s_par
   type(mpl_var)      :: prc
   integer,intent(in) :: tstep
   integer            :: i,j,k,l,lx,ly,lz
   character(*)       :: description
 
   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)

   open(71,file='dump.mtr',position='append')
   write(71,'(a9,i4)') "Timestep ",tstep 
   write(71,*) description 
   write(71,*) " " 
   do k=0,lb_dom%lx(3)
   do j=0,lb_dom%lx(2)
   do i=0,lb_dom%lx(1)
   do l=1,nnod
      write(71,'(5i3,f11.6)') tstep,i,j,k,l,matrix(NDX(l,i,j,k))
   enddo
   enddo
   enddo
   enddo

   write(71,*) " " 
   write(71,*) "-----------------------------------------------------------------------------" 
   write(71,*) " "
   close(71) 

end subroutine dump_matrix 
!------------------------------------------------------------------------
#endif




subroutine gather_output(lb_dom,countOut,tStep,s_par,prc)
!------------------------------------------------------------------------
!
! Master thread gathers the output from all threads in order to write it out.
!


   type(lb_block),intent(inout)    :: lb_dom
   type(sim_parameter),intent(inout)  :: s_par
   type(mpl_var)      :: prc
   integer :: countOut,tStep
   integer :: lx,ly,lz
   integer,parameter :: tagrho=101,tagu=102,tagomega=103,tagstate=104 
   integer,parameter :: tagrho_inc=105
   integer :: b_l(3),b_u(3)
   integer :: ab(prc%np(1),prc%np(2),prc%np(3),3),ae(prc%np(1),prc%np(2),prc%np(3),3)

   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)


   ab = prc%bnd(:,:,:,:,1)
   ae = prc%bnd(:,:,:,:,2)


end subroutine gather_output

!------------------------------------------------------------------------










   subroutine read_params(s_par,filename,prc)


   !------------------------------------------------------------------------
   !
   ! Simulation parameters are read from the file lbc.params
   ! These include
   ! umax         the max velocity, here: of the lid
   ! Re           Reynolds number 
   ! tMax         max. timesteps), tOut (interval for writing out results
   ! problem      1=cylinder, 2=cavity
   ! s_par%gx(1),s_par%gx(2),s_par%gx(3)  global domain size
   ! s_par%setObstacles set small obstacles in equal distance in the whole domain
   ! prc%np(1:3)   number of processes/threads in each direction (cartesian topolos_par%gx(2))
   !

      type(sim_parameter),intent(inout)   :: s_par
      type(mpl_var)   :: prc
      character(10), intent(in)           :: filename
      character(10)           :: temp_name 
      logical :: file_exists     
      integer :: tRamp           ! read from params, must be -1 (no ramp)

      !------------------------------------
      ! Reset some variables and set Name variables for output

      !   Set rho0 
      s_par%rho0 = 1.0d0 

      s_par%goend        = .false.
      s_par%setObstacles = .false.
      s_par%read_obs     = .false.
      s_par%periodic     = .false.
      s_par%initial      = .false.
      s_par%calc_rho_inc = .false.
      s_par%comm_method  =  ""//&
#if defined SENDRECV_BUF
     "sendrecv_buf" &
#elif defined ISEND_IRECV_BUF
 "isend_irecv_buf" &
#elif defined ISEND_IRECV
 "isend_irecv" &
#elif defined SENDRECV
 "sendrecv" &
#elif defined MPI_SUBARRAY 
 "subarray" &
#elif defined PERSISTENT 
 "persistent" &
#elif defined PERSISTENT2
 "persistent" &
#else
 "serial" &
#endif
#if defined MPI_OLD
//"_old"&
#endif
#if defined BUFFER_EACH  
//"B" &
#endif
#if defined COMM_REDUCED 
//"R" &
#endif
//""

      s_par%modelname = "d3q19"
      temp_name       = s_par%modelname
#ifdef INCOMPRESSIBLE
      s_par%modelname = "inc_"//s_par%modelname
#endif
#if defined LAYOUT_LIJK
      s_par%layoutname = "lijk"
#elif defined LAYOUT_IJKL
      s_par%layoutname = "ijkl"
#endif

      !------------------------------------
      ! Read in parameter file

#ifdef USE_MPI
      if(prc%rk == 0) then
#endif
        ! Root thread reads in simulation parameters
         inquire(file=filename,exist=file_exists)
         if(file_exists .eqv. .false.) then
            write(*,*) "Error:"
            write(*,*) "Parameter file '",filename,"' does not exist. "
            write(*,*) "Exiting ... "  
            stop
         endif
         open(1,file=filename)
         read(1,*) s_par%umax                               ! boundary condition, lid speed (0.01<u<0.1)
         read(1,*) s_par%omega                              ! Reynolds number Re<5000?
         read(1,*) s_par%tMax,s_par%tOut,tRamp,s_par%tOutBeg,s_par%tInitial   ! max timesteps, output interval, tRamp (must be -1), tOutBegin, tInitial
         read(1,*) s_par%problem                            
! problem: 1 cylinder in channel, 2 cavity, 3 shock,gaussian=4, cavity_same=5 
         read(1,*) s_par%gx(1),s_par%gx(2),s_par%gx(3)      ! global domain size x, y, z
         read(1,*) s_par%periodic(1),s_par%periodic(2),s_par%periodic(3) ! periodic settings
         read(1,*) prc%np(1),prc%np(2),prc%np(3)            ! # processes in each direction
         read(1,*) s_par%save_output                        ! save output t / f
         read(1,*) s_par%obst_l                             ! Additional obstacle Parameter. For CRVP: Distance of vortices 
         read(1,*) s_par%read_obs                           ! Read obstacles from lbc.obs file??
         close(1)

      !------------------------------------
      ! Enforced invariants: this codebase only supports the lid-driven cavity
      ! (problem=2) without a velocity ramp (tRamp=-1). All other problem
      ! cases and the ramp mechanism have been removed.

         if(s_par%problem /= cavity .or. tRamp /= -1) then
            write(*,*) "Error: unsupported configuration in '",filename,"'."
            write(*,*) "  This build only supports the lid-driven cavity:"
            write(*,*) "    problem = 2 (cavity), tRamp = -1 (no ramp)"
            write(*,*) "  Got: problem = ",s_par%problem,", tRamp = ",tRamp
            write(*,*) "Exiting ... "
            stop
         endif

      !------------------------------------
      ! set characteristic lengths and problem lengths (cylinder radius etc.) 

         if(NDIM==2) s_par%gx(3) = 1
     ! Set position for obstacle (only problem cylinder)
         if(s_par%problem==cylinder) then
         s_par%obst_x=s_par%gx(1)/5+1
         s_par%obst_y=int(real(s_par%gx(2))*(7.d0/12.d0))
         s_par%obst_z=int(real(s_par%gx(3))*(1.d0/2.d0))
         s_par%obst_r=int(real(s_par%gx(2))*(1.d0/15.d0))+1
         elseif(s_par%problem==flute) then
         s_par%obst_x=s_par%gx(1)/5
         s_par%obst_r=int(real(s_par%gx(2))/90.)
         s_par%obst_y=s_par%gx(2)/4-s_par%obst_r
         s_par%obst_z=0
         elseif(s_par%problem==corotating_vortex) then
         s_par%obst_r = s_par%gx(1)/100*s_par%obst_l+1
         elseif(s_par%problem==gaussian .or. s_par%problem==gauss_convect) then
            if(s_par%obst_l < 1.0) then
               s_par%obst_r = s_par%gx(1)/100 +1
            else 
               s_par%obst_r = int(s_par%obst_l) 
            endif
         endif    

     ! Set sponge layer thickness
         s_par%sponge_size = 20
 
     ! set characteristic length lchar
         if(s_par%problem==cylinder .or. s_par%problem==flute) then
            s_par%lchar = real(s_par%obst_r)*2.
         else
         if(NDIM==2) then
            s_par%lchar = real(MAX(s_par%gx(1),s_par%gx(2)))
         else
            s_par%lchar = real(MAX(s_par%gx(1),s_par%gx(2),s_par%gx(3)))
         endif
      endif


      !------------------------------------
      ! set problem names for output 

      if(s_par%problem == cylinder) then
         s_par%problem_name = "cylinder"
      elseif(s_par%problem == cavity) then
         s_par%problem_name = "cavity"
      elseif(s_par%problem == shock) then
         s_par%problem_name = "shock"
      elseif(s_par%problem == gaussian) then
         s_par%problem_name = "gaussian"
      elseif(s_par%problem == gauss_convect) then
         s_par%problem_name = "gauss_conv"
      elseif(s_par%problem == cavity_same) then
         s_par%problem_name = "cavity_same"
      elseif(s_par%problem == flute ) then
         s_par%problem_name = "flute"
      elseif(s_par%problem == corotating_vortex ) then
         s_par%problem_name = "crvp"
      elseif(s_par%problem == cylinder_impulse ) then
         s_par%problem_name = "cyl_imp"
      elseif(s_par%problem == plate_trail ) then
         s_par%problem_name = "plate_trail"
      elseif(s_par%problem == taylor_vortex) then
         s_par%problem_name = "taylor_vortex"
      elseif(s_par%problem == planar_standing_wave) then
         s_par%problem_name = "plane_wave"
      else
         s_par%problem_name = ""
      endif

      ! Calculate Viscosity
       s_par%nu = (2./s_par%omega-1)/6.0
      ! Calculate Re
       s_par%Re = s_par%lchar*s_par%umax/s_par%nu 



      !------------------------------------
      ! Output information about parameters 

      if(prc%rk == 0) then
         ! print out some info 
         write(*,*) 'PARAMETERS    '
         write(*,*) '-----------------------------'
         write(*,'(a12,a)') 'Problem:     ',  s_par%problem_name
         write(*,'(a12,i3,i3,i3)') 'prc%np       ',  prc%np(1),prc%np(2),prc%np(3)
         write(*,'(a12,f6.4)')     'umax        ',   s_par%umax
         write(*,'(a12,f10.4)')     'Re          ',   s_par%Re
         write(*,'(a12,f15.4)')     'Re_acc      ',   s_par%lchar*cs/s_par%nu
         write(*,'(a12,f6.4)')     'nu          ',   s_par%nu
         write(*,'(a12,f8.5)')     'omega       ',   s_par%omega
         write(*,'(a12,i4,i4,i4)') 'Problem size',   s_par%gx(1),s_par%gx(2),s_par%gx(3)
         write(*,'(a12,i8)')     'Max Tstep   ',   s_par%tMax
         write(*,'(a12,i8)')     'Out Tstep   ',   s_par%tOut
         write(*,'(a12,i8)')     'Beg Out at  ',   s_par%tOutBeg
         write(*,'(a12,i8)')     'Init until  ',   s_par%tInitial
         if(s_par%setObstacles .EQV. .TRUE.) write(*,*) 'Insert regular obstacles ', s_par%setObstacles 
         write(*,*) "Model:    D3Q19"
         write(*,*) "Relaxation is: ", &
"bgk"

         write(*,*) "----------------------------------------------"
         write(*,*) 
         if(s_par%tInitial>2) then
            write(*,*) "omega for initial relaxation (u fixed): ",1.1," ...  ",s_par%omega
            write(*,*) "----------------------------------------------"
            write(*,*) 
         endif
     end if
#ifdef USE_MPI
  end if

      !------------------------------------
      ! Broadcast input variables 

   call MPI_BCAST(prc%np,3,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%omega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%bulkViscosity,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%Re,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%umax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_x,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_y,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_z,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%gx,3,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%setObstacles,1,MPI_LOGICAL,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%save_output ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tMax,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tOut,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tOutBeg,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tInitial,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%problem,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%problem_name,16,MPI_CHARACTER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%modelname,16,MPI_CHARACTER,0,MPI_COMM_WORLD,prc%ierr)
#endif /* USE_MPI */
   end subroutine read_params
!------------------------------------------------------------------------








subroutine write_performance_results(result_state,s_par,prc,meas)
!
!  Performance results (MLUPs) are written to a file for statistical review
!  The file-format is simple can be evaluated with gnuplot
!
   implicit none
   type(sim_parameter)  :: s_par
   type(mpl_var)  :: prc  
   type(measure)  :: meas
   integer  :: result_state,lxyztot
   logical       :: file_exists
   character(15) :: filename

  filename='performance.res'
  inquire(file=filename,exist=file_exists)
  open(11,file=filename,position='append') !TRIM(ADJUSTL(filename)))
  if(file_exists .EQV. .false.) then
     write(11,*) '#loglx*ly*lz lx*ly*lz MLUPs  lx   ly   lz   nprc  tTot   tComm  %Comm mpi_meth    layo  resok om bulk umax'
  end if
   lxyztot = s_par%gx(1)*s_par%gx(2)*s_par%gx(3)
  write(11,'(f6.2,i11,f8.2,4i5,3f8.2,5a,i3,3f7.3)') log10(dble(lxyztot)),       &
            & int(s_par%gx(1)*s_par%gx(2)*s_par%gx(3),4),real(meas%mlups,4),    &
            & int(s_par%gx(1),2),int(s_par%gx(2),2), int(s_par%gx(3),2),                                             &
            & int(prc%size,2),meas%mean_total,meas%mean_comm,meas%mean_comm/meas%mean_total*100.,' ',            &
            & s_par%comm_method,s_par%layoutname," ",s_par%modelname,result_state,s_par%omega,s_par%bulkViscosity,s_par%umax
  close(11)
end subroutine write_performance_results
!------------------------------------------------------------------------






subroutine write_temp_data(tStep,nr_unit,size_x,size_y,size_z,size_nnod,fIn,limit_x,limit_y,limit_z,sx,sy,sz)

   integer, intent(in) :: nr_unit,size_x,size_y,size_z,size_nnod,tStep,limit_x,limit_y,limit_z,sx,sy,sz
   real(R8B), intent(in)    :: fIn(NDX(size_nnod,0:size_x+1,0:size_y+1,0:size_z+1))
   integer             :: i,j,k,l
   write(nr_unit,*) "Timestep",tStep
   do k=sz,limit_z
      do j=sy,limit_y
         do i=sx,limit_x
            do l=1,size_nnod
               write(nr_unit,'(3i3,f23.10)') i+1-sx,j+1-sy,k+1-sz,fIn(NDX(l,i,j,k))
            end do 
         end do 
      end do
   end do
end subroutine write_temp_data







subroutine  dealloc_mem(lb_dom,s_par)

!------------------------------------------------------------------------
! Deallocate allocated arrays
!


   type(lb_block)             :: lb_dom
   type(sim_parameter)        :: s_par 
   integer                    :: alloc_stat
   deallocate(lb_dom%obs,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating obs       "
   deallocate(lb_dom%obs_nrwall,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating obs_nrwall"
   deallocate(lb_dom%nrwall_0val,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating nrwall_0val"
   deallocate(lb_dom%nrwall_prev,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating nrwall_prev"
!   deallocate(lb_dom%nrfplus,stat=alloc_stat)
!      if(alloc_stat /=0) write(*,*) "Error deallocating nrfplus    "

   deallocate(lb_dom%fIn,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%fin "
   deallocate(lb_dom%fOut,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%fOut"
   deallocate(lb_dom%u,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%u  "
if(s_par%tInitial > 0) then
   deallocate(lb_dom%u0,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%u0  "
endif
   deallocate(lb_dom%rho,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%rho"
   deallocate(lb_dom%state,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%state"
   deallocate(lb_dom%x,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%x  "
   deallocate(lb_dom%y,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%y  "
   deallocate(lb_dom%z,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%z  "
#ifdef SPARSE
   deallocate(lb_dom%fsIn,stat=alloc_stat)
#endif

end subroutine dealloc_mem






   subroutine write_status(meas,result_state,cur_tstep,s_par,prc)

      implicit none
      type(measure),intent(inout) :: meas
      type(sim_parameter),intent(inout) :: s_par
      type(mpl_var)   :: prc
      integer,intent(in)  :: result_state,cur_tstep
      real(R8B) :: prop_min,prop_max,prop_sum,comm_min,comm_max,comm_sum
      real(R8B) :: sbnd_min,sbnd_max,sbnd_sum,coll_min,coll_max,coll_sum
      real(R8B) :: mpl_cp_buf_max,mpl_cp_buf_min,mpl_cp_buf_sum,mpl_cp_buf_back_max,mpl_cp_buf_back_sum
      real(R8B) :: mpl_cp_buf_back_min,mpl_exch_max,mpl_exch_min,mpl_exch_sum
      real(R8B) :: mpl_sync_max,mpl_sync_min,mpl_sync_sum
      real(R8B) :: totl_min,totl_max,totl_sum
      real(R8B) :: mpl_comm_size_min,mpl_comm_size_max,mpl_comm_size_sum 

      comm_min             = 1.e37 
      mpl_cp_buf_min       = 1.e37 
      mpl_exch_min         = 1.e37 
      mpl_cp_buf_back_min  = 1.e37 
      mpl_exch_min         = 1.e37 
      mpl_sync_min         = 1.e37 
      mpl_comm_size_min    = 1.e37 
      prop_min             = 1.e37 
      sbnd_min             = 1.e37 
      coll_min             = 1.e37 

      comm_max             = 0.
      mpl_cp_buf_max       = 0.
      mpl_exch_max         = 0.
      mpl_cp_buf_back_max  = 0.
      mpl_exch_max         = 0.
      mpl_sync_max         = 0.
      mpl_comm_size_max    = 0.
      prop_max             = 0.
      sbnd_max             = 0.
      coll_max             = 0.

      !-------------------------------
      ! get the min and max values of measured communication times

#ifdef USE_MPI
      call mpi_reduce(meas%duration,totl_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%duration,totl_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%duration,totl_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%sbnd_duration,sbnd_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%sbnd_duration,sbnd_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%sbnd_duration,sbnd_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%coll_duration,coll_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%coll_duration,coll_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%coll_duration,coll_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%prop_duration,prop_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%prop_duration,prop_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%prop_duration,prop_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%comm_duration,comm_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%comm_duration,comm_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%comm_duration,comm_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf,mpl_cp_buf_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf,mpl_cp_buf_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf,mpl_cp_buf_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf_back,mpl_cp_buf_back_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf_back,mpl_cp_buf_back_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf_back,mpl_cp_buf_back_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_exch,mpl_exch_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_exch,mpl_exch_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_exch,mpl_exch_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_sync,mpl_sync_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_sync,mpl_sync_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_sync,mpl_sync_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_comm_size,mpl_comm_size_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_comm_size,mpl_comm_size_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_comm_size,mpl_comm_size_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)


#endif

! Mean values
   if(prc%rk == 0) then
      meas%mean_total = totl_sum/real(prc%size)
      meas%mean_sbnd  = sbnd_sum/real(prc%size)
      meas%mean_coll  = coll_sum/real(prc%size)
      meas%mean_prop  = prop_sum/real(prc%size)
      meas%mean_comm  = comm_sum/real(prc%size)
      meas%mean_mpl_cp_buf  = mpl_cp_buf_sum/real(prc%size)
      meas%mean_mpl_cp_buf_back  = mpl_cp_buf_back_sum/real(prc%size)
      meas%mean_mpl_exch  = mpl_exch_sum/real(prc%size)
      meas%mean_mpl_sync  = mpl_sync_sum/real(prc%size)
      meas%mean_mpl_comm_size  = mpl_comm_size_sum/real(prc%size)
      meas%max_total = totl_max

      meas%mlups = real(cur_tstep,8)*real(s_par%gx(1),8)*real(s_par%gx(2),8)*real(s_par%gx(3),8)/ meas%mean_total/10**6
      write(*,*) 
      write(*,*) "  Statistics" 
      write(*,*) 
      write(*,'(a35,a32)')          "Communication method:           ",s_par%comm_method
      write(*,'(a35,a32)')          "Memory layout of densities:     ",s_par%layoutname
      write(*,*) 
      write(*,'(a35,e10.3)')        "Max. Total wall time:           ", meas%max_total
      write(*,*) 
      write(*,'(a35,f8.3,a12,f8.3,a2)') "Mean Total duration:            ", &
      & meas%mean_total," s     (th0: " , meas%duration, " )" !meas%comm_duration+meas%sbnd_duration+meas%cceq_duration+ &
    !  & meas%coll_duration+meas%prop_duration+meas%ccmv_duration+meas%bncb_duration
      write(*,'(a35,e8.3,a2)') "initialization duration         ", meas%init_duration," s       min     max"
#if defined USE_MPI
      if(prc%np(1)*prc%np(2)*prc%np(3) > 1) then
      write(*,'(a35,f8.3,a2)') "Communication duration          ", &
      & meas%mean_comm," s"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   buffer copy                  ", &
      & meas%mean_mpl_cp_buf   ," s", meas%mean_mpl_cp_buf/meas%mean_comm*100.,"%    (", & 
      & mpl_cp_buf_min,mpl_cp_buf_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   buffer exchange              ", &
      & meas%mean_mpl_exch     ," s", meas%mean_mpl_exch/meas%mean_comm*100.,"%    (",mpl_exch_min,mpl_exch_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   synchronization              ", &
      & meas%mean_mpl_sync     ," s", meas%mean_mpl_sync/meas%mean_comm*100.,"%    (",mpl_sync_min,mpl_sync_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   buffer copy back             ", &
      & meas%mean_mpl_cp_buf_back  ," s", meas%mean_mpl_cp_buf_back/meas%mean_comm*100.,"%    (", & 
      & mpl_cp_buf_back_min,mpl_cp_buf_back_max,")"

      write(*,*) 
      write(*,'(a35,f12.6,a15,f6.3,a4,f6.3)') "   Data size                    ", &
      & meas%mean_mpl_comm_size/1000000.d0,"MB      min,",mpl_comm_size_min/1000000.," max",mpl_comm_size_max/1000000.
      write(*,'(a35,f12.6,a14,f12.6,a5)') "   Data rate                    ", &
      & meas%mean_mpl_comm_size/mpl_exch_min/1000000.,"MB/s     mean ",meas%mean_mpl_comm_size/meas%mean_mpl_exch/1000000.,"MB/s"
      end if
#endif
      write(*,*)
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Communication duration          ", &
      & meas%mean_comm," s", meas%mean_comm/meas%mean_total*100.,"%    (", & 
      & comm_min,comm_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Set boundary  duration          ",  &
& meas%mean_sbnd," s", meas%mean_sbnd/meas%mean_total*100.,"%    (",sbnd_min,sbnd_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a2)') "Calc fEq      duration          ",  &
& meas%mean_cceq," s", meas%mean_cceq/meas%mean_total*100.,"%"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Collision     duration          ",  &
& meas%mean_coll," s", meas%mean_coll/meas%mean_total*100.,"%    (",coll_min,coll_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Propagation   duration          ",  &
& meas%mean_prop," s", meas%mean_prop/meas%mean_total*100.,"%    (",prop_min,prop_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a2)') "Bounceback    duration          ",  &
& meas%mean_bncb," s", meas%mean_bncb/meas%mean_total*100.,"%"
      write(*,'(a35,f8.3,a2,f8.3,a2)') "Calcmacr vals duration          ",  &
& meas%mean_ccmv," s", meas%mean_ccmv/meas%mean_total*100.,"%"
         write(*,*)
      write(*,'(a35,f8.3,a6)') "Performance [MLUPs]             ",meas%mlups," MLUPs"

      call write_performance_results(result_state,s_par,prc,meas)

#ifdef VERBOSE
      write(*,*)
      write(*,*) "Timings of all processes"
      write(*,*)
      write(*,*) "Rk,     set_bnd,     coll,       comm,       prop,       bncb"
#endif /* VERBOSE */
   end if

#ifdef VERBOSE
   call mpl_barrier()
   write(*,'(i4,5f12.8)') prc%rk, meas%sbnd_duration, meas%coll_duration, meas%comm_duration, meas%prop_duration, meas%bncb_duration
#endif /* VERBOSE */

   end subroutine write_status 






subroutine alloc_mem(lb_dom,s_par,prc)

   type(lb_block)        :: lb_dom
   type(mpl_var)         :: prc
   type(sim_parameter)   :: s_par
   integer               :: alloc_stat,lx,ly,lz

   lx = lb_dom%lx(1) 
   ly = lb_dom%lx(2) 
   lz = lb_dom%lx(3) 

   if(prc%rk==0) call calc_mem(prc,lb_dom,s_par)
   allocate(lb_dom%fOut(NDX(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%fOut"
!   if(associated(lb_dom%fOut)) write(*,*) "fOut OK"

if(s_par%tInitial > 0) then
   allocate(lb_dom%u0(  NDX(NDIM,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%u0 "
endif

   allocate(lb_dom%u(  NDX(NDIM,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%u  "
!   if(allocated(lb_dom%u)) write(*,*) "u OK"
   allocate(lb_dom%rho(0:(lx+1),0:(ly+1),0:(lz+1)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%rho"
   allocate(lb_dom%state(lx,ly,lz),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%state"
   allocate(lb_dom%x(lx),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%x  "
   allocate(lb_dom%y(ly),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%y  "
   allocate(lb_dom%z(lz),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%z  "
   allocate(lb_dom%fIn(NDX(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%fIn"
end subroutine alloc_mem 

end module function_set
!------------------------------------------------------------------------


