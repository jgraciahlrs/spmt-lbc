!!
!! Copyright (c) 2010        HLRS, University of Stuttgart. All rights reserved.
!!
!! $COPYRIGHT$
!!
!! Additional copyrights may follow
!!
!! $HEADER$
!!


module timing 
   use nrt_lib
   implicit none
#include "include/replace.h"
contains


!------------------------------------------------------------------------
  subroutine cpu_time_measure(curr_time)
   !
   ! The Cpu Time is measured 
   !
   implicit none
   real(R8B)     curr_time
   curr_time = mpi_wtime()
 end subroutine cpu_time_measure
!------------------------------------------------------------------------


end module timing 
!------------------------------------------------------------------------


