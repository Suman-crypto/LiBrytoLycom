
module libry_opt
contains

! DEFINE BSC TYPES
!---------------------------------------------------------------

subroutine libry_defBSC (k)
integer :: k
end subroutine libry_defBSC

! SUM UP BSC COVER
!---------------------------------------------------------------

subroutine libry_covsBSC (k,k2)
integer :: k,k2
end subroutine libry_covsBSC

! READ NOHONO FLUXES
!---------------------------------------------------------------
subroutine libry_readNOHONO ()
end subroutine libry_readNOHONO

! BROADCAST NOHONO DATA
!---------------------------------------------------------------
subroutine bc_NOHONO ()
end subroutine bc_NOHONO

! ALLOCATE BSC AND NO/HONO VARIABLES
!---------------------------------------------------------------

subroutine libry_allocBSC ()
end subroutine libry_allocBSC

! DEALLOCATE BSC AND NO/HONO VARIABLES
!---------------------------------------------------------------

subroutine libry_deallocBSC ()
end subroutine libry_deallocBSC

! CALCULATE NO AND HONO EMISSIONS
!---------------------------------------------------------------

subroutine libry_fNOHONO (k,st,ac)
integer :: k
real    :: st,ac
end subroutine libry_fNOHONO

! ACCUMULATE BSC PROPERTIES
!---------------------------------------------------------------

subroutine libry_accBSC (i,i2,j)
integer :: i,i2,j
end subroutine libry_accBSC

! AVERAGE BSC PROPERTIES
!---------------------------------------------------------------

subroutine libry_avBSC (i)
integer :: i
end subroutine libry_avBSC

! AVERAGE GRID CELL BSC PROPERTIES
!---------------------------------------------------------------

subroutine libry_avgcBSC (i,i2,j)
integer :: i,i2,j
end subroutine libry_avgcBSC

end module libry_opt

