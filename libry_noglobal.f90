
! INITIALISE MPI
!---------------------------------------------------------------

subroutine init_mpi ()
use libry_par
implicit none

rank = 0

para = .false.

return
end subroutine init_mpi

! Dummy LIBRY_GLOBAL*
!---------------------------------------------------------------

module libry_global2

end module libry_global2

!---------------------------------------------------------------

module libry_global
contains

! other subroutines
!---------------------------------------------------------------
subroutine read_land ()
end subroutine read_land

!---------------------------------------------------------------
subroutine alloc_global ()
end subroutine alloc_global

!---------------------------------------------------------------
subroutine bc_pp ()
end subroutine bc_pp

!---------------------------------------------------------------
subroutine bc_namelist ()
end subroutine bc_namelist

!---------------------------------------------------------------
subroutine bc_specpar ()
end subroutine bc_specpar

!---------------------------------------------------------------
subroutine libry_readBCg ()
end subroutine libry_readBCg

!---------------------------------------------------------------
subroutine sc_BC ()
end subroutine sc_BC

!---------------------------------------------------------------
subroutine libry_openG ()
end subroutine libry_openG 

!---------------------------------------------------------------
subroutine libry_readHourG (dummy)
  integer :: dummy
end subroutine libry_readHourG

!---------------------------------------------------------------
subroutine sc_clim ()
end subroutine sc_clim 

!---------------------------------------------------------------
subroutine libry_outputG ()
end subroutine libry_outputG 

!---------------------------------------------------------------
subroutine def_varG ()
end subroutine def_varG 

!---------------------------------------------------------------
subroutine write_varG ()
end subroutine write_varG 

!---------------------------------------------------------------
subroutine libry_outspecG ()
end subroutine libry_outspecG 

!---------------------------------------------------------------
subroutine libry_closeG ()
end subroutine libry_closeG 

!---------------------------------------------------------------
subroutine dealloc_global ()
end subroutine dealloc_global 

!---------------------------------------------------------------
subroutine libry_read_restart (dummy)
  integer :: dummy
end subroutine libry_read_restart 

!---------------------------------------------------------------
subroutine sc_restart ()
end subroutine sc_restart

!---------------------------------------------------------------
subroutine libry_write_restart ()
end subroutine libry_write_restart 

end module libry_global

