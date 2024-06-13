
!###############################################################
! LIBRY_NC
!###############################################################

module libry_nc
contains

! CHECK STATUS
!---------------------------------------------------------------

subroutine check(status)
use netcdf

  integer, intent(in) :: status
                     
  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop 2
  end if

end subroutine check

end module libry_nc

!###############################################################
! LIBRY_GLOBAL2
!###############################################################

module libry_global2
contains

! DEFINE GLOBAL OUTPUT VARIABLES
!---------------------------------------------------------------

subroutine def_varG (outIDx, code0, ovID0, nDimsX)
use libry_par
use netcdf
use libry_nc
implicit none

integer, intent(in)                     :: code0, outIDx, nDimsX
integer, intent(out)                    :: ovID0
character (len=8)                       :: vd0

! Define output variables
write(vd0,'(A4,I4)') "code",code0

select case( nDimsX )
  case( 0 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs0, ovID0)) ! returns outvarID
  case( 1 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs1, ovID0))
  case( 2 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs2, ovID0))
  case( 3 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs3, ovID0))
  case( 4 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs4, ovID0))
  case( 5 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs5, ovID0))
  case( 6 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs6, ovID0))
  case( 7 )
    call check(nf90_def_var(outIDx, vd0, NF90_REAL, dimIDs7, ovID0))
  case default
    write(*,*) "error in define out var"
    stop
end select

return
end subroutine def_varG


! WRITE GLOBAL OUTPUT VARIABLES
!---------------------------------------------------------------

subroutine write_varG (outIDx, var0, ovID1)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,k,l
integer                 :: mperr

integer, intent(in)     :: outIDx, ovID1
integer, dimension(1)   :: tcode

real, dimension(pp) :: var0

real, allocatable, dimension(:) :: var1      !(ppnp)  
real, allocatable, dimension(:) :: outvarL   !(nland) 
real, allocatable, dimension(:,:) :: outvar    !(nx,ny) 

allocate(var1(ppnp))
allocate(outvarL(nland))
allocate(outvar(nx,ny))

!DEBUG
write(*,*) "grid size"
write(*,*) "nx",nx
write(*,*) "ny",ny
write(*,*) " "

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

! Gather output variable
call MPI_GATHER( var0, pp, MPI_DOUBLE_PRECISION, var1, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

if (rank .eq. 0) then

  ! assemble list of output points from all processors
  k=1
  do i = 1,numproc

    outvarL(k:k+ppvec(i)-1) = var1((i-1)*pp+1:(i-1)*pp+ppvec(i))

    k = k + ppvec(i)
  enddo

  ! fillvalue
  outvar(:,:) = -9999.0

  ! distribute list of output land points to output map
  do l=1,nland

    outvar(indvec(l,1),indvec(l,2)) = outvarL(l)
  enddo

  ! Write output variable
   tcode(1) = (year*100 + month)*100 + day

   call check(nf90_put_var(outIDx, t_varID, tcode, start = (/tpos/), count = (/1/) ))

   call check(nf90_put_var(outIDx, ovID1, outvar(:,:), start = (/ 1, 1, tpos /), count = (/ nx, ny, 1 /) ))
endif

deallocate(var1)
deallocate(outvarL)
deallocate(outvar)

return
end subroutine write_varG

! WRITE GLOBAL RESTART VARIABLES - MULTIPLE DIMENSIONS
!---------------------------------------------------------------

subroutine write_varG0 (rsfIDx, varX, rvIDx)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,v2,v3,j,m,n,l,k
integer                 :: mperr

integer                 :: rvIDx, rsfIDx
integer, dimension(1)   :: tcode

real, dimension(pp) :: varX

real, allocatable, dimension(:) :: varXb      !(ppnp)  
real, allocatable, dimension(:) :: rsvarL   !(nland) 
real, allocatable, dimension(:,:) :: rsvar    !(nx,ny) 

!DEBUG
! write(*,*) "RANK", rank
! write(*,*) "write var G0 START"

allocate(varXb(ppnp))
allocate(rsvarL(nland))
allocate(rsvar(nx,ny))

!DEBUG
! write(*,*) "RANK", rank
! write(*,*) "MPI barrier start"

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

!DEBUG
! write(*,*) "RANK", rank
! write(*,*) "MPI barrier end"
! write(*,*) "size varx", size(varX)
! write(*,*) "sum varx", sum(varX(:))

! Gather restart variable
call MPI_GATHER( varX, pp, MPI_DOUBLE_PRECISION, varXb, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

if (rank .eq. 0) then

!DEBUG
! write(*,*) "assemble start"

  ! assemble list of restart points from all processors
  k=1
  do i = 1,numproc

    rsvarL(k:k+ppvec(i)-1) = varXb((i-1)*pp+1:(i-1)*pp+ppvec(i))

    k = k + ppvec(i)
  enddo

  ! fillvalue
  rsvar(:,:) = -9999.0

  ! distribute list of restart land points to restart map
  do l=1,nland

    rsvar(indvec(l,1),indvec(l,2)) = rsvarL(l)
  enddo

!DEBUG
! write(*,*) "assemble finished, start write"

!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "SNOW -- DIRECTLY BEFORE WRITE"
!write(*,*) " "
!write(*,*) "cell   80", rsvar(indvec(80,1),indvec(80,2))
!write(*,*) "cell    1", rsvar(indvec(01,1),indvec(01,2))
!write(*,*) "cell  345", rsvar(indvec(345,1),indvec(345,2))
!write(*,*) "cell  871", rsvar(indvec(871,1),indvec(871,2))
!write(*,*) "cell 1450", rsvar(indvec(1450,1),indvec(1450,2))
!write(*,*) " "
!write(*,*) " "
!write(*,*) " "



  ! Write restart variable
   call check(nf90_put_var(rsfIDx, rvIDx, rsvar(:,:), start = (/ 1, 1 /), &
                                                      count = (/ nx, ny /) ))
endif

deallocate(varXb)
deallocate(rsvarL)
deallocate(rsvar)

!DEBUG
! write(*,*) "RANK", rank
! write(*,*) "write var G0 END"

return
end subroutine write_varG0

!---------------------------------------------------------------

subroutine write_varG1 (rsfIDx, varX, rvIDx)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,v2,v3,j,m,n,l,k
integer                 :: mperr

integer                 :: rvIDx, rsfIDx
integer, dimension(1)   :: tcode

real, dimension(pp,p_ntiles) :: varX
real, dimension(pp) :: varX1d

real, allocatable, dimension(:) :: varXb        !(ppnp)  
real, allocatable, dimension(:,:) :: rsvarL     !(nland) 
real, allocatable, dimension(:,:,:) :: rsvar    !(nx,ny) 

allocate(varXb(ppnp))
allocate(rsvarL(nland,p_ntiles))
allocate(rsvar(nx,ny,p_ntiles))

! assemble list of restart points from all processors
do n = 1,p_ntiles

  varX1d(:) = varX(:, n)

  call MPI_BARRIER(MPI_COMM_WORLD, mperr)

  ! Gather restart variable
  call MPI_GATHER( varX1d, pp, MPI_DOUBLE_PRECISION, varXb, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

  call MPI_BARRIER(MPI_COMM_WORLD, mperr)

  if (rank .eq. 0) then
    k=1
    do i = 1,numproc

      rsvarL(k:k+ppvec(i)-1, n) = varXb((i-1)*pp+1:(i-1)*pp+ppvec(i))

      k = k + ppvec(i)
    enddo
  endif
enddo

if (rank .eq. 0) then

  ! fillvalue
  rsvar(:,:,:) = -9999.0

  ! distribute list of restart land points to restart map
  do l=1,nland

    rsvar(indvec(l,1),indvec(l,2),:) = rsvarL(l,:)
  enddo

  ! Write restart variable
   call check(nf90_put_var(rsfIDx, rvIDx, rsvar(:,:,:), start = (/ 1, 1, 1 /), &
                                                      count = (/ nx, ny, p_ntiles /) ))
endif

deallocate(varXb)
deallocate(rsvarL)
deallocate(rsvar)

return
end subroutine write_varG1

!---------------------------------------------------------------

subroutine write_varG2 (rsfIDx, varX, rvIDx)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,v2,v3,j,m,n,l,k
integer                 :: mperr

integer                 :: rvIDx, rsfIDx
integer, dimension(1)   :: tcode

real, dimension(pp,p_ntiles,p_nhabM) :: varX
real, dimension(pp) :: varX1d

real, allocatable, dimension(:) :: varXb        !(ppnp)  
real, allocatable, dimension(:,:,:) :: rsvarL   !(nland) 
real, allocatable, dimension(:,:,:,:) :: rsvar  !(nx,ny) 

allocate(varXb(ppnp))
allocate(rsvarL(nland,p_ntiles,p_nhabM))
allocate(rsvar(nx,ny,p_ntiles,p_nhabM))

! assemble list of restart points from all processors
do m = 1,p_nhabM
  do n = 1,p_ntiles

    varX1d(:) = varX(:, n,m)

    call MPI_BARRIER(MPI_COMM_WORLD, mperr)

    ! Gather restart variable
    call MPI_GATHER( varX1d, pp, MPI_DOUBLE_PRECISION, varXb, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

    call MPI_BARRIER(MPI_COMM_WORLD, mperr)

    if (rank .eq. 0) then
      k=1
      do i = 1,numproc

        rsvarL(k:k+ppvec(i)-1, n,m) = varXb((i-1)*pp+1:(i-1)*pp+ppvec(i))

        k = k + ppvec(i)
      enddo
    endif
  enddo
enddo

if (rank .eq. 0) then

  ! fillvalue
  rsvar(:,:,:,:) = -9999.0

  ! distribute list of restart land points to restart map
  do l=1,nland

    rsvar(indvec(l,1),indvec(l,2),:,:) = rsvarL(l,:,:)
  enddo

  ! Write restart variable
   call check(nf90_put_var(rsfIDx, rvIDx, rsvar(:,:,:,:), start = (/ 1, 1, 1, 1 /), &
                                                      count = (/ nx, ny, p_ntiles, p_nhabM /) ))
endif

deallocate(varXb)
deallocate(rsvarL)
deallocate(rsvar)

return
end subroutine write_varG2

!---------------------------------------------------------------

subroutine write_varG4 (rsfIDx, varX, rvIDx)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,v2,v3,j,m,n,l,k
integer                 :: mperr

integer                 :: rvIDx, rsfIDx
integer, dimension(1)   :: tcode

real, dimension(pp,p_ntiles,p_nvertM) :: varX
real, dimension(pp) :: varX1d

real, allocatable, dimension(:) :: varXb        !(ppnp)  
real, allocatable, dimension(:,:,:) :: rsvarL   !(nland) 
real, allocatable, dimension(:,:,:,:) :: rsvar  !(nx,ny) 

allocate(varXb(ppnp))
allocate(rsvarL(nland,p_ntiles,p_nvertM))
allocate(rsvar(nx,ny,p_ntiles,p_nvertM))

! assemble list of restart points from all processors
do v2 = 1,p_nvertM
  do n = 1,p_ntiles

    varX1d(:) = varX(:, n,v2)

    call MPI_BARRIER(MPI_COMM_WORLD, mperr)

    ! Gather restart variable
    call MPI_GATHER( varX1d, pp, MPI_DOUBLE_PRECISION, varXb, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

    call MPI_BARRIER(MPI_COMM_WORLD, mperr)

    if (rank .eq. 0) then
      k=1
      do i = 1,numproc

        rsvarL(k:k+ppvec(i)-1, n,v2) = varXb((i-1)*pp+1:(i-1)*pp+ppvec(i))

        k = k + ppvec(i)
      enddo
    endif
  enddo
enddo

if (rank .eq. 0) then

  ! fillvalue
  rsvar(:,:,:,:) = -9999.0

  ! distribute list of restart land points to restart map
  do l=1,nland

    rsvar(indvec(l,1),indvec(l,2),:,:) = rsvarL(l,:,:)
  enddo

  ! Write restart variable
   call check(nf90_put_var(rsfIDx, rvIDx, rsvar(:,:,:,:), start = (/ 1, 1, 1, 1 /), &
                                                      count = (/ nx, ny, p_ntiles, p_nvertM /) ))
endif

deallocate(varXb)
deallocate(rsvarL)
deallocate(rsvar)

return
end subroutine write_varG4

!---------------------------------------------------------------

subroutine write_varG5 (rsfIDx, varX, rvIDx)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,v2,v3,j,m,n,l,k
integer                 :: mperr

integer                 :: rvIDx, rsfIDx
integer, dimension(1)   :: tcode

real, dimension(pp,p_ntiles,p_nvertM,p_nhabM) :: varX
real, dimension(pp) :: varX1d

real, allocatable, dimension(:) :: varXb        !(ppnp)  
real, allocatable, dimension(:,:,:,:) :: rsvarL !(nland) 
real, allocatable, dimension(:,:,:,:,:) :: rsvar !(nx,ny) 

allocate(varXb(ppnp))
allocate(rsvarL(nland,p_ntiles,p_nvertM,p_nhabM))
allocate(rsvar(nx,ny,p_ntiles,p_nvertM,p_nhabM))

! assemble list of restart points from all processors
do m = 1,p_nhabM
  do v3 = 1,p_nvertM
    do n = 1,p_ntiles

      varX1d(:) = varX(:, n,v3,m)

      call MPI_BARRIER(MPI_COMM_WORLD, mperr)

      ! Gather restart variable
      call MPI_GATHER( varX1d, pp, MPI_DOUBLE_PRECISION, varXb, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

      call MPI_BARRIER(MPI_COMM_WORLD, mperr)

      if (rank .eq. 0) then

        k=1
        do i = 1,numproc

          rsvarL(k:k+ppvec(i)-1, n,v3,m) = varXb((i-1)*pp+1:(i-1)*pp+ppvec(i))

          k = k + ppvec(i)
        enddo
      endif
    enddo
  enddo
enddo

if (rank .eq. 0) then

  ! fillvalue
  rsvar(:,:,:,:,:) = -9999.0

  ! distribute list of restart land points to restart map
  do l=1,nland

    rsvar(indvec(l,1),indvec(l,2),:,:,:) = rsvarL(l,:,:,:)
  enddo

  ! Write restart variable
   call check(nf90_put_var(rsfIDx, rvIDx, rsvar(:,:,:,:,:), start = (/ 1, 1, 1, 1, 1 /), &
                                                      count = (/ nx, ny, p_ntiles,p_nvertM,p_nhabM /) ))
endif

deallocate(varXb)
deallocate(rsvarL)
deallocate(rsvar)

return
end subroutine write_varG5

!---------------------------------------------------------------

subroutine write_varG6 (rsfIDx, varX, rvIDx)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,v2,v3,j,m,n,l,k
integer                 :: mperr

integer                 :: rvIDx, rsfIDx
integer, dimension(1)   :: tcode
!DEBUG
!integer :: x1,x2,x3

real, dimension(pp,p_ntiles,p_nvertM,p_nhabM,p_nspec) :: varX
real, dimension(pp) :: varX1d

real, allocatable, dimension(:) :: varXb                !(ppnp)  
real, allocatable, dimension(:,:,:,:,:) :: rsvarL       !(nland) 
real, allocatable, dimension(:,:,:,:,:,:) :: rsvar      !(nx,ny) 

allocate(varXb(ppnp))
allocate(rsvarL(nland,p_ntiles,p_nvertM,p_nhabM,p_nspec))
allocate(rsvar(nx,ny,p_ntiles,p_nvertM,p_nhabM,p_nspec))

!DEBUG
!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "AREA fractions in cell 80 -- VAR G6 START"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", varX(80,x1,1,1,x2)
!    write(*,*) "hab 2", varX(80,x1,1,2,x2)
!    write(*,*) "level 2", varX(80,x1,2,1,x2)
!    write(*,*) "hab 2", varX(80,x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo

! assemble list of restart points from all processors
do j = 1,p_nspec
  do m = 1,p_nhabM
    do v3 = 1,p_nvertM
      do n = 1,p_ntiles

        varX1d(:) = varX(:, n,v3,m,j)

        call MPI_BARRIER(MPI_COMM_WORLD, mperr)

        ! Gather restart variable
        call MPI_GATHER( varX1d, pp, MPI_DOUBLE_PRECISION, varXb, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )
        
        call MPI_BARRIER(MPI_COMM_WORLD, mperr)

        if (rank .eq. 0) then

          k=1
          do i = 1,numproc

            rsvarL(k:k+ppvec(i)-1, n,v3,m,j) = varXb((i-1)*pp+1:(i-1)*pp+ppvec(i))

            k = k + ppvec(i)
          enddo
        endif

      enddo
    enddo
  enddo
enddo

!DEBUG
!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "AREA fractions in cell 80 -- VAR G6 AFTER GATHER & ASSEMBLE"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", rsvarL(80,x1,1,1,x2)
!    write(*,*) "hab 2", rsvarL(80,x1,1,2,x2)
!    write(*,*) "level 2", rsvarL(80,x1,2,1,x2)
!    write(*,*) "hab 2", rsvarL(80,x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo


if (rank .eq. 0) then

  ! fillvalue
  rsvar(:,:,:,:,:,:) = -9999.0

  ! distribute list of restart land points to restart map
  do l=1,nland

    rsvar(indvec(l,1),indvec(l,2),:,:,:,:) = rsvarL(l,:,:,:,:)
  enddo


!DEBUG
!!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "AREA fractions in cell 80 -- VAR G6 DIRECTLY BEFORE WRITE"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", rsvar(indvec(80,1),indvec(80,2),x1,1,1,x2)
!    write(*,*) "hab 2", rsvar(indvec(80,1),indvec(80,2),x1,1,2,x2)
!    write(*,*) "level 2", rsvar(indvec(80,1),indvec(80,2),x1,2,1,x2)
!    write(*,*) "hab 2", rsvar(indvec(80,1),indvec(80,2),x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo


  ! Write restart variable
   call check(nf90_put_var(rsfIDx, rvIDx, rsvar(:,:,:,:,:,:), start = (/ 1, 1, 1, 1, 1, 1 /), &
                                                      count = (/ nx, ny, p_ntiles,p_nvertM,p_nhabM,p_nspec /) ))
endif

deallocate(varXb)
deallocate(rsvarL)
deallocate(rsvar)

return
end subroutine write_varG6

!---------------------------------------------------------------

subroutine write_varG7 (rsfIDx, varX, rvIDx)
use libry_par
use netcdf
use libry_nc
use mpi
implicit none

integer                 :: i,v2,v3,j,m,n,l,k
integer                 :: mperr

integer                 :: rvIDx, rsfIDx
integer, dimension(1)   :: tcode

real, dimension(pp,p_ntiles,p_nhabM,p_nspec) :: varX
real, dimension(pp) :: varX1d

real, allocatable, dimension(:) :: varXb        !(ppnp)  
real, allocatable, dimension(:,:,:,:) :: rsvarL !(nland) 
real, allocatable, dimension(:,:,:,:,:) :: rsvar !(nx,ny) 

allocate(varXb(ppnp))
allocate(rsvarL(nland,p_ntiles,p_nhabM,p_nspec))
allocate(rsvar(nx,ny,p_ntiles,p_nhabM,p_nspec))

! assemble list of restart points from all processors
do j = 1,p_nspec
  do m = 1,p_nhabM
    do n = 1,p_ntiles

      varX1d(:) = varX(:, n,m,j)

      call MPI_BARRIER(MPI_COMM_WORLD, mperr)

      ! Gather restart variable
      call MPI_GATHER( varX1d, pp, MPI_DOUBLE_PRECISION, varXb, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )

      call MPI_BARRIER(MPI_COMM_WORLD, mperr)

      if (rank .eq. 0) then

        k=1
        do i = 1,numproc

          rsvarL(k:k+ppvec(i)-1, n,m,j) = varXb((i-1)*pp+1:(i-1)*pp+ppvec(i))

          k = k + ppvec(i)
        enddo
      endif
    enddo
  enddo
enddo

if (rank .eq. 0) then

  ! fillvalue
  rsvar(:,:,:,:,:) = -9999.0

  ! distribute list of restart land points to restart map
  do l=1,nland

    rsvar(indvec(l,1),indvec(l,2),:,:,:) = rsvarL(l,:,:,:)
  enddo

  ! Write restart variable
   call check(nf90_put_var(rsfIDx, rvIDx, rsvar(:,:,:,:,:), start = (/ 1, 1, 1, 1, 1 /), &
                                                      count = (/ nx, ny, p_ntiles,p_nhabM,p_nspec /) ))
endif

deallocate(varXb)
deallocate(rsvarL)
deallocate(rsvar)

return
end subroutine write_varG7

end module libry_global2


!###############################################################
! LIBRY_GLOBAL
!###############################################################

module libry_global
contains

! INITIALISE MPI
!---------------------------------------------------------------

subroutine init_mpi ()
use libry_par
use mpi
implicit none

integer :: mperr

call MPI_INIT(mperr)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mperr)

call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, mperr)

para = .true.

allocate(ppvec(numproc))

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

return
end subroutine init_mpi


! READ LAND FILE
!---------------------------------------------------------------

subroutine read_land ()
use libry_par
use netcdf
use libry_nc
implicit none

integer                 :: i, j, l
integer                 :: ncID1, ncID2, varID
integer                 :: ndims, xtype
integer                 :: pp0, addp, cntp
integer, dimension(2)   :: dimIDs2n
character (len=20)      :: dimName, varName0

! Open land mask file and determine resolution
call check(nf90_open(landmask, nf90_nowrite, ncID1) )
!!!DEBUG -- TEMPORARY FIX FOR PALEO !!!
!!call check(nf90_inquire_dimension(ncID1, 1, dimName, nx) )
!!call check(nf90_inquire_dimension(ncID1, 2, dimName, ny) )

call check(nf90_inquire_dimension(ncID1, 2, dimName, nx) )      !! TEMPORARY
call check(nf90_inquire_dimension(ncID1, 3, dimName, ny) )      !! TEMPORARY

!!!DEBUG!!!
write(*,*) "read_land"
write(*,*) "nx=",nx
write(*,*) "ny=",ny
write(*,*) "  "


allocate(lsdata(nx,ny))
allocate(xpos(nx))
allocate(ypos(ny))

! Get the values of the coordinates and put them in xpos & ypos
!!call check(nf90_inquire_variable(ncID1, 1, varName0, xtype, ndims, dimIDs2n))
call check(nf90_inquire_variable(ncID1, 2, varName0, xtype, ndims, dimIDs2n))     !! TEMPORARY
call check(nf90_inq_varID(ncID1, varName0, varID))
call check(nf90_get_var(ncID1, varID, xpos))

!!call check(nf90_inquire_variable(ncID1, 2, varName0, xtype, ndims, dimIDs2n))
call check(nf90_inquire_variable(ncID1, 3, varName0, xtype, ndims, dimIDs2n))     !! TEMPORARY
call check(nf90_inq_varID(ncID1, varName0, varID))
call check(nf90_get_var(ncID1, varID, ypos))

! Get the values of the land data and put them in lsdata
call check(nf90_inq_varID(ncID1, varName, varID))
call check(nf90_get_var(ncID1, varID, lsdata))
call check(nf90_close(ncID1))

! Open one climate data file (tair) and read length
call check(nf90_open(tairfileG, nf90_nowrite, ncID2) )
!!call check(nf90_inquire_dimension(ncID2, 3, dimName, nt) )
call check(nf90_inquire_dimension(ncID2, 1, dimName, nt) )   !! TEMPORARY
call check(nf90_close(ncID2))


!!!DEBUG!!!
write(*,*) "dimName=",dimName
write(*,*) "nt=",nt
write(*,*) "  "



nland = sum(lsdata(:,:))

allocate(indvec(nland,2))

! Create list of land points
l = 0
do i = 1,nx
  do j = 1,ny

    if (lsdata(i,j) .gt. 0.0) then
      l = l + 1
      indvec(l,1) = i
      indvec(l,2) = j
    endif
  enddo
enddo

! Compute land points per processor (pp)

pp0 = int(floor(real(nland / numproc)))

ppvec(:) = pp0

addp = modulo(nland, numproc)

if (addp .gt. 0) then
  pp = pp0 + 1

  cntp = 0
  do i = 1, numproc
    cntp = cntp + 1
    if (cntp .le. addp) ppvec(i) = ppvec(i) + 1
  enddo
else
  pp = pp0
endif

! total land points + extra points
ppnp = pp * numproc

deallocate(lsdata)

return
end subroutine read_land


! BROADCAST POINTS PER PROCESSOR
!---------------------------------------------------------------

subroutine bc_pp ()
use libry_par
use mpi
implicit none

integer :: mperr

call MPI_BCAST( pp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr ) ! pp: dim=1, bcast process is root (0)

call MPI_BCAST( ppnp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( nland, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

call MPI_BCAST( ppvec, numproc, MPI_INTEGER, 0, MPI_COMM_WORLD, mperr )

return
end subroutine bc_pp


! ALLOCATE GLOBAL FIELDS
!---------------------------------------------------------------

subroutine alloc_global ()
use libry_par
implicit none

if (rank .eq. 0) then

  allocate( indata(nx,ny), &
  indataL(nland), &
  bcdata(nx,ny,12), &
  bcdata0(nx,ny), &    !! TEMPORARY
  !!bcdata(12,nx,ny), &    !! TEMPORARY
  bcdataL(nland,12), &
  bcdata2(nx,ny), &
  bcdata2L(nland), &
  rsdata1(nx,ny), &
  rsdata2(nx,ny,p_ntiles,p_nvertM,p_nhabM), &
  rsdata3(nx,ny,p_ntiles), &
  rsdata4(nx,ny,p_ntiles), &
  rsdata5(nx,ny,p_ntiles,p_nvertM), &
  rsdata6(nx,ny,p_ntiles), &
  rsdata7(nx,ny,p_ntiles,p_nhabM), &
  rsdata1L(nland), &
  rsdata2L(nland,p_ntiles,p_nvertM,p_nhabM), &
  rsdata3L(nland,p_ntiles), &
  rsdata4L(nland,p_ntiles), &
  rsdata5L(nland,p_ntiles,2), &
  rsdata6L(nland,p_ntiles), &
  rsdata7L(nland,p_ntiles,p_nhabM), &
  rsdataV1(nx,ny,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
  rsdataV2(nx,ny,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
  rsdataV3(nx,ny,p_ntiles,p_nhabM,p_nspec), &
  rsdataV1L(nland,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
  rsdataV2L(nland,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
  rsdataV3L(nland,p_ntiles,p_nhabM,p_nspec) )
endif

allocate( laidata(ppnp,12), &
saidata(ppnp,12), &
ssadata(ppnp), &
biomedata(ppnp), &
indata7(ppnp,7), &
restartdata1(ppnp), &
restartdata2(ppnp,p_ntiles,p_nvertM,p_nhabM), &
restartdata3(ppnp,p_ntiles), &
restartdata4(ppnp,p_ntiles), &
restartdata5(ppnp,p_ntiles,2), &
restartdata6(ppnp,p_ntiles), &
restartdata7(ppnp,p_ntiles,p_nhabM), &
restartdataV1(ppnp,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
restartdataV2(ppnp,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
restartdataV3(ppnp,p_ntiles,p_nhabM,p_nspec) )

! Initialize with fill value
laidata(:,:) = 1.0E12
saidata(:,:) = 1.0E12
ssadata(:)   = 1.0E12
biomedata(:) = 1.0E12

indata7(:,:) = 1.0E12

restartdata1(:)       = 1.0E12
restartdata2(:,:,:,:) = 1.0E12
restartdata3(:,:)     = 1.0E12
restartdata4(:,:)     = 1.0E12
restartdata5(:,:,:)   = 1.0E12
restartdata6(:,:)     = 1.0E12
restartdata7(:,:,:)   = 1.0E12

restartdataV1(:,:,:,:,:) = 1.0E12
restartdataV2(:,:,:,:,:) = 1.0E12
restartdataV3(:,:,:,:)   = 1.0E12

return
end subroutine alloc_global


! BROADCAST NAMELIST
!---------------------------------------------------------------

subroutine bc_namelist ()
use libry_par
use mpi
implicit none

integer :: mperr

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

call MPI_BCAST( year0,          1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( cyear0,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( tsindata0,      1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( accts0,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( tpos0,          1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( lastyear,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( runperiod,      1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( tsl,            1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( yearout1,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( yearoutX,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( outint,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( nSites,         1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( p_nspec,        1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( fracTH_s_init,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( lrestart,       1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( NOHONO,         1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( BSCtypes,       1, MPI_LOGICAL,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( llevels,        1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )
call MPI_BCAST( ltraits,        1, MPI_INTEGER,          0, MPI_COMM_WORLD, mperr )

return
end subroutine bc_namelist


! BROADCAST SPECIES PARAMETERS
!---------------------------------------------------------------

subroutine bc_specpar ()
use libry_par
use mpi
implicit none

integer :: mperr, k

do k = 1,p_nspecpar

  call MPI_BCAST( vec_o(:,k),  p_nspec, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )
enddo

return
end subroutine bc_specpar


! READ BOUNDARY CONDITIONS (GLOBAL)
!---------------------------------------------------------------

subroutine libry_readBCg ()
use libry_par
use netcdf
use libry_nc
implicit none

integer :: ncID, varID, i,k,l,m
!!DEBUG
integer :: nx2,ny2,nt2,lv
character (len=20)      :: dimName



!!DEBUG
write(*,*) "read BC started"

! LAI
call check(nf90_open(laifileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))

!!DEBUG
!!call check(nf90_inquire_dimension(ncID, 1, dimName, nx2) )      !! TEMPORARY
!!call check(nf90_inquire_dimension(ncID, 2, dimName, ny2) )      !! TEMPORARY
!!call check(nf90_inquire_dimension(ncID, 3, dimName, nt2) )      !! TEMPORARY
!!write(*,*) "dim1  :",nx2  !12
!!write(*,*) "dim2  :",ny2  !144
!!write(*,*) "dim3  :",nt2  !143
!!write(*,*) "  "
!!write(*,*) "calling get_var with nx :",11, " and ny :", 11, "time 1"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count = (/ 1, 11, 11 /) ))    !!  TEMPORARY
!!
!!write(*,*) "calling get_var with nx :",12, " and ny :", 12, "time 1"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count = (/ 1, 12, 12 /) ))    !!  TEMPORARY
!!
!!!!write(*,*) "calling get_var with nx :",nx-1, " and ny :", ny-1, "time 1"
!!!!
!!!! call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count = (/ 1, nx-1, ny-1 /) ))    !!  TEMPORARY
!!
!!write(*,*) "calling get_var with nx :",nx-1, " and ny :", ny-1, "time 1 at POS 3"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count = (/ nx-1, ny-1, 1 /) ))    !!  TEMPORARY
!!
!!write(*,*) "calling get_var with nx :",nx, " and ny :", ny, "time 1 at POS 3"
!!
 !!call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
!!
!!!!write(*,*) "calling get_var with nx :",nx, " and ny :", ny, "time 1"
!!!!
!!!! call check(nf90_get_var(ncID, varID, bcdata0, start = (/ 1, 1, 1 /), count = (/ 1, nx, ny /) ))    !!  TEMPORARY
!!
!!write(*,*) "bcdata0 finished"
!!


!!call check(nf90_get_var(ncID, varID, bcdata))
do lv = 1,12
 call check(nf90_get_var(ncID, varID, bcdata(:,:,lv), start = (/ 1, 1, lv /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
enddo

call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdataL(l,:) = bcdata(indvec(l,1),indvec(l,2),:)
enddo

! distribute boundary condition data to processors
do m = 1,12
  k=1
  do i = 1,numproc
    laidata((i-1)*pp+1:(i-1)*pp+ppvec(i), m) = bcdataL(k:k+ppvec(i)-1, m)

    k = k + ppvec(i)
  enddo
enddo


!!DEBUG
write(*,*) "read SAI started "


! SAI
call check(nf90_open(saifileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
!!call check(nf90_get_var(ncID, varID, bcdata))
do lv = 1,12
 call check(nf90_get_var(ncID, varID, bcdata(:,:,lv), start = (/ 1, 1, lv /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
enddo
call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdataL(l,:) = bcdata(indvec(l,1),indvec(l,2),:)
enddo

! distribute boundary condition data to processors
do m = 1,12
  k=1
  do i = 1,numproc
    saidata((i-1)*pp+1:(i-1)*pp+ppvec(i), m) = bcdataL(k:k+ppvec(i)-1, m)

    k = k + ppvec(i)
  enddo
enddo

! SSA
call check(nf90_open(ssafileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
!!call check(nf90_get_var(ncID, varID, bcdata2))
write(*,*) "start SSA"
 call check(nf90_get_var(ncID, varID, bcdata2, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY

call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdata2L(l) = bcdata2(indvec(l,1),indvec(l,2))
enddo

! distribute boundary condition data to processors
k=1
do i = 1,numproc
  ssadata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = bcdata2L(k:k+ppvec(i)-1)

  k = k + ppvec(i)
enddo

! biomes
call check(nf90_open(biomefileG, nf90_nowrite, ncID) )
call check(nf90_inq_varID(ncID, varName, varID))
!!call check(nf90_get_var(ncID, varID, bcdata2))
write(*,*) "start biome"
 call check(nf90_get_var(ncID, varID, bcdata2, start = (/ 1, 1, 1 /), count = (/ nx, ny, 1 /) ))    !!  TEMPORARY
call check(nf90_close(ncID))

! assign boundary condition data to list of land points
do l=1,nland
  bcdata2L(l) = bcdata2(indvec(l,1),indvec(l,2))
enddo

! distribute boundary condition data to processors
k=1
do i = 1,numproc
  biomedata((i-1)*pp+1:(i-1)*pp+ppvec(i)) = bcdata2L(k:k+ppvec(i)-1)

  k = k + ppvec(i)
enddo


!!DEBUG
write(*,*) "read BC finished"

return
end subroutine libry_readBCg


! SCATTER BOUNDARY CONDITIONS
!---------------------------------------------------------------

subroutine sc_BC ()
use libry_par
use mpi
implicit none

integer :: m, mperr

real, allocatable :: laimon(:)
real, allocatable :: saimon(:)

allocate(laimon(pp))
allocate(saimon(pp))

do m = 1,12

  call MPI_BARRIER(MPI_COMM_WORLD, mperr)

! parts of v0 have length pp*8760, root=0
  call MPI_SCATTER( laidata(:,m), pp, MPI_DOUBLE_PRECISION, &
                    laimon, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

  call MPI_SCATTER( saidata(:,m), pp, MPI_DOUBLE_PRECISION, &
                    saimon, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

  areaLEAF(:,m) = laimon(:)
  areaSTEM(:,m) = saimon(:)
enddo

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

call MPI_SCATTER( ssadata, pp, MPI_DOUBLE_PRECISION, areaSOIL, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( biomedata, pp, MPI_DOUBLE_PRECISION, biome, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

deallocate(laimon)
deallocate(saimon)

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

return
end subroutine sc_BC

! OPEN AND READ RESTART FILES
!---------------------------------------------------------------

subroutine libry_read_restart (nCPts3)
use libry_par
use netcdf
use libry_nc

implicit none

integer :: i,v2,v3,j,m,n,l,k,x
integer :: ncID1,ncID2,varID
integer :: nCPts3

character (len=8), dimension(7) :: vd1
character (len=8), dimension(3) :: vd2

!DEBUG
character (len=8)       :: name0, name1
integer                 :: len0, xtype0, ndims0, varV1
!integer, dimension(4)   :: dimIDs2B                                     ! dimension IDs for restart file
!integer, dimension(6)   :: dimIDs6B                                     ! dimension IDs for restart file
!integer :: x1,x2,x3

 write(*,*) "  "
 write(*,*) "LIBRY READ RESTART STARTED"
 write(*,*) "  "


do x = 1,7
  write(vd1(x),'(A4,I4)') "code",kfile_restartL+x
enddo

do x = 1,3
  write(vd2(x),'(A4,I4)') "code",kfile_restartV+x
enddo

! open restart files

call check(nf90_open(restartfileL, nf90_nowrite, ncID1) )
call check(nf90_open(restartfileV, nf90_nowrite, ncID2) )

!DEBUG
! write(*,*) "restart files opened"

! land surface
call check(nf90_inq_varID(ncID1, vd1(1), varID))
call check(nf90_get_var(ncID1, varID, rsdata1))
call check(nf90_inq_varID(ncID1, vd1(2), varID))
call check(nf90_get_var(ncID1, varID, rsdata2))
call check(nf90_inq_varID(ncID1, vd1(3), varID))
call check(nf90_get_var(ncID1, varID, rsdata3))
call check(nf90_inq_varID(ncID1, vd1(4), varID))
call check(nf90_get_var(ncID1, varID, rsdata4))
call check(nf90_inq_varID(ncID1, vd1(5), varID))
call check(nf90_get_var(ncID1, varID, rsdata5))
call check(nf90_inq_varID(ncID1, vd1(6), varID))
call check(nf90_get_var(ncID1, varID, rsdata6))
call check(nf90_inq_varID(ncID1, vd1(7), varID))
call check(nf90_get_var(ncID1, varID, rsdata7))

! write(*,*) "vd1(1):", vd1(1)
! write(*,*) "vd1(2):", vd1(2)
! write(*,*) "vd1(3):", vd1(3)
! write(*,*) "vd1(4):", vd1(4)
! write(*,*) "vd1(5):", vd1(5)
! write(*,*) "vd1(6):", vd1(6)
! write(*,*) "vd1(7):", vd1(7)
! write(*,*) "size of rsdata7:", size(rsdata7)
! write(*,*) "  "
! write(*,*) "rsdata7(40,25,1,1):", rsdata7(40,25,1,1)
! write(*,*) "rsdata7(40,25,4,1):", rsdata7(40,25,4,1)
! write(*,*) "rsdata5(80,50,1,2):", rsdata5(80,50,1,2)


call check(nf90_close(ncID1))

!DEBUG
! write(*,*) "get land vars - finished"


 !  call check(nf90_inquire_dimension(ncID2, 1 , name0, len0))
 !  write(*,*) "name", name0
 !  write(*,*) "length", len0
 !  call check(nf90_inquire_dimension(ncID2, 2 , name0, len0))
 !  write(*,*) "name", name0
 !  write(*,*) "length", len0
 !  call check(nf90_inquire_dimension(ncID2, 3 , name0, len0))
 !  write(*,*) "name", name0
 !  write(*,*) "length", len0
 !  call check(nf90_inquire_dimension(ncID2, 4 , name0, len0))
 !  write(*,*) "name", name0
 !  write(*,*) "length", len0
 !  call check(nf90_inquire_dimension(ncID2, 5 , name0, len0))
 !  write(*,*) "name", name0
 !  write(*,*) "length", len0
 !  call check(nf90_inquire_dimension(ncID2, 6 , name0, len0))
 !  write(*,*) "name", name0
 !  write(*,*) "length", len0

 !  write(*,*) "  "
 !  write(*,*) "rsdataV1:d1",size(rsdataV1,1)
 !  write(*,*) "rsdataV1:d2",size(rsdataV1,2)
 !  write(*,*) "rsdataV1:d3",size(rsdataV1,3)
 !  write(*,*) "rsdataV1:d4",size(rsdataV1,4)
 !  write(*,*) "rsdataV1:d5",size(rsdataV1,5)
 !  write(*,*) "rsdataV1:d6",size(rsdataV1,6)
 !  write(*,*) "  "

 !  call check(nf90_inq_varID(ncID2, "code6201", varV1))
 !  call check(nf90_inquire_variable(ncID2, varV1, name1, xtype0, ndims0, dimIDs6B))
 !  write(*,*) "name: ", name1
 !  write(*,*) "xtype: ", xtype0
 !  write(*,*) "ndims: ", ndims0
 !  write(*,*) "dimids: ", dimIDs6B
 !  write(*,*) "  "
 !  call check(nf90_get_var(ncID2, varV1, rsdataV1))
 !  write(*,*) "size of rsdataV1:", size(rsdataV1)

 !  call check(nf90_inq_varID(ncID2, "code6202", varV1))
 !  call check(nf90_inquire_variable(ncID2, varV1, name1, xtype0, ndims0, dimIDs6B))
 !  write(*,*) "name: ", name1
 !  write(*,*) "xtype: ", xtype0
 !  write(*,*) "ndims: ", ndims0
 !  write(*,*) "dimids: ", dimIDs6B
 !  write(*,*) "  "
 !  call check(nf90_get_var(ncID2, varV1, rsdataV2))
 !  write(*,*) "size of rsdataV2:", size(rsdataV2)



! vegetation

call check(nf90_inq_varID(ncID2, vd2(1), varID))
call check(nf90_get_var(ncID2, varID, rsdataV1))
call check(nf90_inq_varID(ncID2, vd2(2), varID))
call check(nf90_get_var(ncID2, varID, rsdataV2))
call check(nf90_inq_varID(ncID2, vd2(3), varID))
call check(nf90_get_var(ncID2, varID, rsdataV3))

call check(nf90_close(ncID2))

!DEBUG
! write(*,*) "get veg vars - finished"
! write(*,*) "   "
! write(*,*) "   "
! write(*,*) " CHECK READ-IN VALUES -- BEFORE ASSIGN"
! write(*,*) "   "
!!! ONLY WORKS WITH ONE PROCESSOR!!!
!
!write(*,*) "SNOW"
!write(*,*) " "
!write(*,*) "cell   80", rsdata1(indvec(80,1),indvec(80,2))
!write(*,*) "cell    1", rsdata1(indvec(01,1),indvec(01,2))
!write(*,*) "cell  345", rsdata1(indvec(345,1),indvec(345,2))
!write(*,*) "cell  871", rsdata1(indvec(871,1),indvec(871,2))
!write(*,*) "cell 1450", rsdata1(indvec(1450,1),indvec(1450,2))
!write(*,*) " "
!write(*,*) " "
!write(*,*) " "
!write(*,*) "AREA fractions in cell 80:"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", rsdataV2(indvec(80,1),indvec(80,2),x1,1,1,x2)
!    write(*,*) "hab 2", rsdataV2(indvec(80,1),indvec(80,2),x1,1,2,x2)
!    write(*,*) "level 2", rsdataV2(indvec(80,1),indvec(80,2),x1,2,1,x2)
!    write(*,*) "hab 2", rsdataV2(indvec(80,1),indvec(80,2),x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo



! assign land restart data to list of land points
do l=1,nland
  rsdata1L(l)       = rsdata1(indvec(l,1),indvec(l,2))
  rsdata2L(l,:,:,:) = rsdata2(indvec(l,1),indvec(l,2),:,:,:)
  rsdata3L(l,:)     = rsdata3(indvec(l,1),indvec(l,2),:)
  rsdata4L(l,:)     = rsdata4(indvec(l,1),indvec(l,2),:)
  rsdata5L(l,:,:)   = rsdata5(indvec(l,1),indvec(l,2),:,:)
  rsdata6L(l,:)     = rsdata6(indvec(l,1),indvec(l,2),:)
  rsdata7L(l,:,:)   = rsdata7(indvec(l,1),indvec(l,2),:,:)
enddo

! assign ground/canopy restart data to list of land points
do l=1,nland
  rsdataV1L(l,:,:,:,:) = rsdataV1(indvec(l,1),indvec(l,2),:,:,:,:)
  rsdataV2L(l,:,:,:,:) = rsdataV2(indvec(l,1),indvec(l,2),:,:,:,:)
  rsdataV3L(l,:,:,:)   = rsdataV3(indvec(l,1),indvec(l,2),:,:,:)
enddo



!DEBUG
! write(*,*) "   "
! write(*,*) "   "
! write(*,*) " CHECK READ-IN VALUES -- BEFORE DISTRIB"
! write(*,*) "   "

!write(*,*) sum(rsdata1L(:)      )
!write(*,*) sum(rsdata2L(:,:,:,:))
!write(*,*) sum(rsdata3L(:,:)    )
!write(*,*) sum(rsdata4L(:,:)    )
!write(*,*) sum(rsdata5L(:,:,:)  )
!write(*,*) sum(rsdata6L(:,:)    )
!write(*,*) sum(rsdata7L(:,:,:)  )
!write(*,*) "   "
!write(*,*) sum(rsdataV1L(:,:,:,:,:) )
!write(*,*) sum(rsdataV2L(:,:,:,:,:) )
!write(*,*) sum(rsdataV3L(:,:,:,:)   )
!write(*,*) "   "
!write(*,*) "rsdata7(40,25,1,1):", rsdata7(40,25,1,1)


!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "area fractions in cell 80:"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", rsdataV2L(80,x1,1,1,x2)
!    write(*,*) "hab 2", rsdataV2L(80,x1,1,2,x2)
!    write(*,*) "level 2", rsdataV2L(80,x1,2,1,x2)
!    write(*,*) "hab 2", rsdataV2L(80,x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo



! distribute land restart data to processors
k=1
do i = 1,numproc
  restartdata1((i-1)*pp+1:(i-1)*pp+ppvec(i)) = rsdata1L(k:k+ppvec(i)-1)
  k = k + ppvec(i)
enddo

do n = 1,p_ntiles
  k=1
  do i = 1,numproc
    restartdata3((i-1)*pp+1:(i-1)*pp+ppvec(i), n) = rsdata3L(k:k+ppvec(i)-1, n)
    restartdata4((i-1)*pp+1:(i-1)*pp+ppvec(i), n) = rsdata4L(k:k+ppvec(i)-1, n)
    restartdata6((i-1)*pp+1:(i-1)*pp+ppvec(i), n) = rsdata6L(k:k+ppvec(i)-1, n)
    k = k + ppvec(i)
  enddo
enddo

do m = 1,p_nhabM
  do n = 1,p_ntiles
    k=1
    do i = 1,numproc
      restartdata7((i-1)*pp+1:(i-1)*pp+ppvec(i), n,m) = rsdata7L(k:k+ppvec(i)-1, n,m)
      k = k + ppvec(i)
    enddo
  enddo
enddo

do v2 = 1,p_nvertM
  do n = 1,p_ntiles
    k=1
    do i = 1,numproc
      restartdata5((i-1)*pp+1:(i-1)*pp+ppvec(i), n,v2) = rsdata5L(k:k+ppvec(i)-1, n,v2)
      k = k + ppvec(i)
    enddo
  enddo
enddo

do m = 1,p_nhabM
  do v3 = 1,p_nvertM
    do n = 1,p_ntiles
      k=1
      do i = 1,numproc
        restartdata2((i-1)*pp+1:(i-1)*pp+ppvec(i), n,v3,m) = rsdata2L(k:k+ppvec(i)-1, n,v3,m)
        k = k + ppvec(i)
      enddo
    enddo
  enddo
enddo

! distribute ground/canopy restart data to processors

do j = 1,p_nspec
  do m = 1,p_nhabM
    do v3 = 1,p_nvertM
      do n = 1,p_ntiles
        k=1
        do i = 1,numproc
          restartdataV1((i-1)*pp+1:(i-1)*pp+ppvec(i), n,v3,m,j) = rsdataV1L(k:k+ppvec(i)-1, n,v3,m,j)
          restartdataV2((i-1)*pp+1:(i-1)*pp+ppvec(i), n,v3,m,j) = rsdataV2L(k:k+ppvec(i)-1, n,v3,m,j)
          if(v3 .eq. 1) restartdataV3((i-1)*pp+1:(i-1)*pp+ppvec(i), n,m,j) = rsdataV3L(k:k+ppvec(i)-1, n,m,j)
          k = k + ppvec(i)
        enddo
      enddo
    enddo
  enddo
enddo


!DEBUG
! write(*,*) "  "
! write(*,*) "LIBRY READ RESTART FINISHED -- AFTER DISTRIB"
! write(*,*) "  "
! write(*,*) "  "


!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "area fractions in cell 80:"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", restartdataV2(80,x1,1,1,x2)
!    write(*,*) "hab 2", restartdataV2(80,x1,1,2,x2)
!    write(*,*) "level 2", restartdataV2(80,x1,2,1,x2)
!    write(*,*) "hab 2", restartdataV2(80,x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo


return
end subroutine libry_read_restart

! SCATTER RESTART DATA
!---------------------------------------------------------------

subroutine sc_restart ()
use libry_par
use mpi
implicit none

integer :: j,m,n,v2,v3, mperr

real, allocatable :: vartemp(:)


!DEBUG
! write(*,*) "RANK", rank
! write(*,*) "LIBRY scatter restart -- STARTED"
! write(*,*) "SIZE restartdata1", size(restartdata1)
! write(*,*) "SIZE rH2Os_g", size(rH2Os_g)
! write(*,*) "SUM rH2Os_g", sum(rH2Os_g)
! write(*,*) "pp", pp
! write(*,*) "  "
!integer :: x1,x2,x3


allocate(vartemp(pp))

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

! write(*,*) "  "
! write(*,*) "SCATTER RESTART: BARRIER FINISHED"
! write(*,*) "  "

call MPI_SCATTER( restartdata1, pp, MPI_DOUBLE_PRECISION, rH2Os_g, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

do n = 1,p_ntiles

  call MPI_BARRIER(MPI_COMM_WORLD, mperr)

  call MPI_SCATTER( restartdata3(:,n), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
  rH2Ol_g1(:,n) = vartemp(:)

  call MPI_SCATTER( restartdata4(:,n), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
  rH2Ol_g2(:,n) = vartemp(:)

  call MPI_SCATTER( restartdata6(:,n), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
  fH2Ol_ug2(:,n) = vartemp(:)
enddo

do m = 1,p_nhabM
  do n = 1,p_ntiles

    call MPI_BARRIER(MPI_COMM_WORLD, mperr)

    call MPI_SCATTER( restartdata7(:,n,m), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
    xT_g0(:,n,m) = vartemp(:)
  enddo
enddo

do v2 = 1,p_nvertM
  do n = 1,p_ntiles

    call MPI_BARRIER(MPI_COMM_WORLD, mperr)

    call MPI_SCATTER( restartdata5(:,n,v2), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
    atN_fH2Ol_xd(:,n,v2) = vartemp(:)
  enddo
enddo

do m = 1,p_nhabM
  do v3 = 1,p_nvertM
    do n = 1,p_ntiles
      call MPI_BARRIER(MPI_COMM_WORLD, mperr)

      call MPI_SCATTER( restartdata2(:,n,v3,m), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
      rH2Ol_0(:,n,v3,m) = vartemp(:)
    enddo
  enddo
enddo

do j = 1,p_nspec
  do m = 1,p_nhabM
    do v3 = 1,p_nvertM
      do n = 1,p_ntiles

        call MPI_BARRIER(MPI_COMM_WORLD, mperr)

        call MPI_SCATTER( restartdataV1(:,n,v3,m,j), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
        rH2Ol_t(:,n,v3,m,j) = vartemp(:)

        call MPI_SCATTER( restartdataV2(:,n,v3,m,j), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
        areaTH_s(:,n,v3,m,j) = vartemp(:)
      enddo
    enddo
  enddo
enddo

do j = 1,p_nspec
  do m = 1,p_nhabM
    do n = 1,p_ntiles

      call MPI_BARRIER(MPI_COMM_WORLD, mperr)

      call MPI_SCATTER( restartdataV3(:,n,m,j), pp, MPI_DOUBLE_PRECISION, vartemp, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)
        xT_g(:,n,m,j) = vartemp(:)
    enddo
  enddo
enddo

deallocate(vartemp)

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

!DEBUG

!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "area fractions in cell 80:"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", areaTH_s(80,x1,1,1,x2)
!    write(*,*) "hab 2", areaTH_s(80,x1,1,2,x2)
!    write(*,*) "level 2", areaTH_s(80,x1,2,1,x2)
!    write(*,*) "hab 2", areaTH_s(80,x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo



return
end subroutine sc_restart

! OPEN GLOBAL FILES
!---------------------------------------------------------------

subroutine libry_openG ()
use libry_par
use netcdf
use libry_nc
implicit none

! Open climate data files
call check(nf90_open(tairfileG, nf90_nowrite, ncID0(1)) )
call check(nf90_open(rhumfileG, nf90_nowrite, ncID0(2)) )
call check(nf90_open(windfileG, nf90_nowrite, ncID0(3)) )
call check(nf90_open(rainfileG, nf90_nowrite, ncID0(4)) )
call check(nf90_open(snowfileG, nf90_nowrite, ncID0(5)) )
call check(nf90_open(sradfileG, nf90_nowrite, ncID0(6)) )
call check(nf90_open(lradfileG, nf90_nowrite, ncID0(7)) )

return
end subroutine libry_openG


! READ GLOBAL INPUT DATA FILES
!---------------------------------------------------------------

subroutine libry_readHourG (h)
use libry_par
use netcdf
use libry_nc
implicit none

integer :: varID, c,i,h,k,l

! loop over all 7 climate forcing files


!!!DEBUG!!!
write(*,*) "starting read hour G "
write(*,*) "hour =",h
write(*,*) " "

do c = 1,7
  call check(nf90_inq_varID(ncID0(c), varName, varID))
  call check(nf90_get_var(ncID0(c),varID,indata(:,:), start = (/ 1, 1, h /), count = (/ nx, ny, 1 /) ))

  !! fails !!!call check(nf90_get_var(ncID0(c),varID,indata(:,:), start = (/ h, 1, 1 /), count = (/ 1, nx, ny /) ))    !!  TEMPORARY




  ! assign boundary condition data to list of land points
  do l = 1,nland
    indataL(l) = indata(indvec(l,1),indvec(l,2))
  enddo

  ! distribute boundary condition data to processors
  k=1
  do i = 1,numproc
    indata7((i-1)*pp+1:(i-1)*pp+ppvec(i), c) = indataL(k:k+ppvec(i)-1)

    k = k + ppvec(i)
  enddo
enddo


!!!DEBUG!!!
write(*,*) "read hour G finished "

return
end subroutine libry_readHourG


! SCATTER GLOBAL INPUT DATA
!---------------------------------------------------------------

subroutine sc_clim ()
use libry_par
use mpi
implicit none

integer :: mperr

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,1), pp, MPI_DOUBLE_PRECISION, &
                  xT_a, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,2), pp, MPI_DOUBLE_PRECISION, &
                  rH2Og_RH, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,3), pp, MPI_DOUBLE_PRECISION, &
                  fAIR_s, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,4), pp, MPI_DOUBLE_PRECISION, &
                  fH2Ol_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,5), pp, MPI_DOUBLE_PRECISION, &
                  fH2Os_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,6), pp, MPI_DOUBLE_PRECISION, &
                  fRADs_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_SCATTER( indata7(:,7), pp, MPI_DOUBLE_PRECISION, &
                  fRADl_ad, pp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr)

call MPI_BARRIER(MPI_COMM_WORLD, mperr)

return
end subroutine sc_clim


! GLOBAL OUTPUT
!---------------------------------------------------------------

subroutine libry_outputG ()
use libry_par
use libry_global2
use netcdf
use libry_nc
implicit none

integer :: x_dimID, y_dimID, t_dimID, k
integer :: x_varID, y_varID

if (rank .eq. 0) then

  if (out1) then

    ! Create global output file
    call check(nf90_create(sfile_outputG, NF90_CLOBBER, outID))
   
    ! Define dimensions.
    call check(nf90_def_dim(outID, "lon", nx, x_dimID))    ! returns dimID
    call check(nf90_def_dim(outID, "lat", ny, y_dimID))
    call check(nf90_def_dim(outID, "time", NF90_UNLIMITED, t_dimID))
    dimIDs3 = (/ x_dimID, y_dimID, t_dimID /)
   
    ! Define coordinate variables
    call check(nf90_def_var(outID, "lon", NF90_REAL, x_dimID, x_varID))    ! returns varID
    call check(nf90_put_att(outID, x_varID, "long_name", "longitude"))
    call check(nf90_put_att(outID, x_varID, "units", "degrees_east"))
    call check(nf90_put_att(outID, x_varID, "standard_name", "longitude"))
    call check(nf90_def_var(outID, "lat", NF90_REAL, y_dimID, y_varID))
    call check(nf90_put_att(outID, y_varID, "long_name", "latitude"))
    call check(nf90_put_att(outID, y_varID, "units", "degrees_north"))
    call check(nf90_put_att(outID, y_varID, "standard_name", "latitude"))
   
    call check(nf90_def_var(outID, "time", NF90_INT, t_dimID, t_varID))
    call check(nf90_put_att(outID, t_varID, "units", "day as %Y%m%d.%f"))
    call check(nf90_put_att(outID, t_varID, "calendar", "proleptic_gregorian"))

    ! Define output variables

    do k = 1,2
      call def_varG(outID, k*1000+kfile_rCO2d          ,  outvarID(1+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_rCb            ,  outvarID(2+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_rH2Ol          ,  outvarID(3+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_rH2Os          ,  outvarID(4+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_rmaxH2Ol       ,  outvarID(5+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_area           ,  outvarID(6+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_act            ,  outvarID(7+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fCO2gc         ,  outvarID(8+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fCcg           ,  outvarID(9+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fCcb           , outvarID(10+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fCbo           , outvarID(11+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fH2Ol_ux       , outvarID(12+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fH2Ol_xd       , outvarID(13+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fH2Ogl_ux      , outvarID(14+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_fH2Olg_xu      , outvarID(15+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_Ts             , outvarID(16+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_H              , outvarID(17+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_E              , outvarID(18+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_C              , outvarID(19+20*(k-1)), 3)
      call def_varG(outID, k*1000+kfile_EB             , outvarID(20+20*(k-1)), 3)
    enddo

    call def_varG(outID, 3000+kfile_rH2Ol_g1           , outvarID(100), 3) ! must be larger than highest index of outvarID 3 lines above
    call def_varG(outID, 3000+kfile_rH2Ol_g2           , outvarID(101), 3)
    call def_varG(outID, 3000+kfile_rH2Os_g            , outvarID(102), 3)
    call def_varG(outID, 3000+kfile_fH2Ol_ug           , outvarID(103), 3)
    call def_varG(outID, 3000+kfile_fH2Ol_ug2          , outvarID(104), 3)
    call def_varG(outID, 3000+kfile_fH2Ol_go           , outvarID(105), 3)
    call def_varG(outID, 3000+kfile_fH2Ol_gb           , outvarID(106), 3)
    call def_varG(outID, 3000+kfile_fH2Olg_ga          , outvarID(107), 3)
    call def_varG(outID, 3000+kfile_fH2Osl_g           , outvarID(108), 3)
    call def_varG(outID, 3000+kfile_fH2Os_ad           , outvarID(109), 3)
    call def_varG(outID, 3000+kfile_fH2Ol_ad           , outvarID(110), 3)
    call def_varG(outID, 3000+kfile_xT_a               , outvarID(111), 3)
    call def_varG(outID, 3000+kfile_Tg                 , outvarID(112), 3) 
    call def_varG(outID, 3000+kfile_G                  , outvarID(113), 3)
    call def_varG(outID, 3000+kfile_fRADs              , outvarID(114), 3)

    ! Optional variables
    if (BSCtypes) then
      call def_varG(outID, kfile_MareaLCgS     , outvarID2(1), 3)                
      call def_varG(outID, kfile_MareaDCgS     , outvarID2(2), 3)                
      call def_varG(outID, kfile_MareaCCgS     , outvarID2(3), 3)
      call def_varG(outID, kfile_MareaMCgS     , outvarID2(4), 3)
      call def_varG(outID, kfile_MrH2OlLC_S    , outvarID2(5), 3)
      call def_varG(outID, kfile_MrH2OlDC_S    , outvarID2(6), 3)
      call def_varG(outID, kfile_MrH2OlCC_S    , outvarID2(7), 3)
      call def_varG(outID, kfile_MrH2OlMC_S    , outvarID2(8), 3)
      call def_varG(outID, kfile_MTsLC_S       , outvarID2(9), 3)
      call def_varG(outID, kfile_MTsDC_S       , outvarID2(10), 3)
      call def_varG(outID, kfile_MTsCC_S       , outvarID2(11), 3)
      call def_varG(outID, kfile_MTsMC_S       , outvarID2(12), 3)
      call def_varG(outID, kfile_MfCcbLC_S     , outvarID2(13), 3)
      call def_varG(outID, kfile_MfCcbDC_S     , outvarID2(14), 3)
      call def_varG(outID, kfile_MfCcbCC_S     , outvarID2(15), 3)
      call def_varG(outID, kfile_MfCcbMC_S     , outvarID2(16), 3)

      if (NOHONO) then

        call def_varG(outID, kfile_MfNO_N      , outvarID2(20), 3)
        call def_varG(outID, kfile_MfHONO_N    , outvarID2(21), 3)
      endif
    endif

    call check(nf90_enddef(outID)) !End Definitions  

    ! Write grid variables
    call check(nf90_put_var(outID, x_varID, xpos))
    call check(nf90_put_var(outID, y_varID, ypos))
  else

    ! Open output file for next write
    call check(nf90_open(sfile_outputG, NF90_WRITE, outID))
  endif
endif

! Write output (variable, code, outfileID, dimensions)

do k = 1,2
  call write_varG(outID, ag_rCO2d(:,k)     ,  outvarID(1+20*(k-1)))
  call write_varG(outID, ag_rCb(:,k)       ,  outvarID(2+20*(k-1)))
  call write_varG(outID, ag_rH2Ol(:,k)     ,  outvarID(3+20*(k-1)))
  call write_varG(outID, ag_rH2Os(:,k)     ,  outvarID(4+20*(k-1)))
  call write_varG(outID, ag_rmaxH2Ol(:,k)  ,  outvarID(5+20*(k-1)))
  call write_varG(outID, ag_areaTH_s(:,k)  ,  outvarID(6+20*(k-1)))
  call write_varG(outID, ag_act(:,k)       ,  outvarID(7+20*(k-1)))
  call write_varG(outID, ag_fCO2gc(:,k)    ,  outvarID(8+20*(k-1)))
  call write_varG(outID, ag_fCcg(:,k)      ,  outvarID(9+20*(k-1)))
  call write_varG(outID, ag_fCcb(:,k)      , outvarID(10+20*(k-1)))
  call write_varG(outID, ag_fCbo(:,k)      , outvarID(11+20*(k-1)))
  call write_varG(outID, ag_fH2Ol_ux(:,k)  , outvarID(12+20*(k-1)))
  call write_varG(outID, ag_fH2Ol_xd(:,k)  , outvarID(13+20*(k-1)))
  call write_varG(outID, ag_fH2Ogl_ux(:,k) , outvarID(14+20*(k-1)))
  call write_varG(outID, ag_fH2Olg_xu(:,k) , outvarID(15+20*(k-1)))
  call write_varG(outID, ag_Ts(:,k)        , outvarID(16+20*(k-1)))
  call write_varG(outID, ag_H(:,k)         , outvarID(17+20*(k-1)))
  call write_varG(outID, ag_E(:,k)         , outvarID(18+20*(k-1)))
  call write_varG(outID, ag_C(:,k)         , outvarID(19+20*(k-1)))
  call write_varG(outID, ag_EB(:,k)        , outvarID(20+20*(k-1)))
enddo

call write_varG(outID, ag_rH2Ol_g1(:)        , outvarID(100))
call write_varG(outID, ag_rH2Ol_g2(:)        , outvarID(101))
call write_varG(outID, ag_rH2Os_g(:)         , outvarID(102))
call write_varG(outID, ag_fH2Ol_ug(:)        , outvarID(103))
call write_varG(outID, ag_fH2Ol_ug2(:)       , outvarID(104))
call write_varG(outID, ag_fH2Ol_go(:)        , outvarID(105))
call write_varG(outID, ag_fH2Ol_gb(:)        , outvarID(106))
call write_varG(outID, ag_fH2Olg_ga(:)       , outvarID(107))
call write_varG(outID, ag_fH2Osl_g(:)        , outvarID(108))
call write_varG(outID, ag_fH2Os_ad(:)        , outvarID(109))
call write_varG(outID, ag_fH2Ol_ad(:)        , outvarID(110))
call write_varG(outID, ag_xT_a(:)            , outvarID(111))
call write_varG(outID, ag_Tg(:)              , outvarID(112)) 
call write_varG(outID, ag_G(:)               , outvarID(113))
call write_varG(outID, ag_fRADs(:)           , outvarID(114))

! Optional output variables
if (BSCtypes) then
  call write_varG(outID, a_MareaTHLC_g_S(:), outvarID2(1))                 
  call write_varG(outID, a_MareaTHDC_g_S(:), outvarID2(2))                 
  call write_varG(outID, a_MareaTHCC_g_S(:), outvarID2(3))
  call write_varG(outID, a_MareaTHMC_g_S(:), outvarID2(4))
  call write_varG(outID, a_MrH2OlLC_S(:) , outvarID2(5))
  call write_varG(outID, a_MrH2OlDC_S(:) , outvarID2(6))
  call write_varG(outID, a_MrH2OlCC_S(:) , outvarID2(7))
  call write_varG(outID, a_MrH2OlMC_S(:) , outvarID2(8))
  call write_varG(outID, a_MTsLC_S(:)    , outvarID2(9))
  call write_varG(outID, a_MTsDC_S(:)    , outvarID2(10))
  call write_varG(outID, a_MTsCC_S(:)    , outvarID2(11))
  call write_varG(outID, a_MTsMC_S(:)    , outvarID2(12))
  call write_varG(outID, a_MfCcbLC_S(:)  , outvarID2(13))
  call write_varG(outID, a_MfCcbDC_S(:)  , outvarID2(14))
  call write_varG(outID, a_MfCcbCC_S(:)  , outvarID2(15))
  call write_varG(outID, a_MfCcbMC_S(:)  , outvarID2(16))

  if (NOHONO) then

    call write_varG(outID, a_MfNO_N(:)   , outvarID2(20))
    call write_varG(outID, a_MfHONO_N(:) , outvarID2(21))
  endif
endif

if (rank .eq. 0) call check(nf90_close(outID))

return
end subroutine libry_outputG




! WRITE RESTART FILES
!---------------------------------------------------------------

subroutine libry_write_restart ()
use libry_par
use libry_global2
use netcdf
use mpi
use libry_nc
implicit none

integer :: k
integer :: mperr

integer :: x_dimID, y_dimID, tl_dimID, vr_dimID, nh_dimID, sp_dimID
!integer :: x_dimID2, y_dimID2, tl_dimID2, vr_dimID2, nh_dimID2
integer :: x_varID, y_varID, tl_varID, vr_varID, nh_varID, sp_varID
!integer :: x_varID2, y_varID2, tl_varID2, vr_varID2, nh_varID2, sp_varID2


!DEBUG
! write(*,*) "RANK", rank
! write(*,*) "LIBRY write restart -- STARTED"
!integer :: x1,x2,x3


if (rank .eq. 0) then

  ! Create global restart files
  call check(nf90_create(restartfileL, NF90_CLOBBER, rsfID))
  call check(nf90_create(restartfileV, NF90_CLOBBER, rsfID2))
   
  ! Define dimensions.
  call check(nf90_def_dim(rsfID, "lon", nx, x_dimID))    ! returns dimID
  call check(nf90_def_dim(rsfID, "lat", ny, y_dimID))
  call check(nf90_def_dim(rsfID, "tile", p_ntiles, tl_dimID))
  call check(nf90_def_dim(rsfID, "vert", p_nvertM, vr_dimID))
  call check(nf90_def_dim(rsfID, "nhab", p_nhabM, nh_dimID))

  call check(nf90_def_dim(rsfID2, "lon", nx, x_dimID))    ! returns dimID
  call check(nf90_def_dim(rsfID2, "lat", ny, y_dimID))
  call check(nf90_def_dim(rsfID2, "tile", p_ntiles, tl_dimID))
  call check(nf90_def_dim(rsfID2, "vert", p_nvertM, vr_dimID))
  call check(nf90_def_dim(rsfID2, "nhab", p_nhabM, nh_dimID))
  call check(nf90_def_dim(rsfID2, "spec", p_nspec, sp_dimID))

  dimIDs0 = (/ x_dimID, y_dimID /)
  dimIDs1 = (/ x_dimID, y_dimID, tl_dimID /)
  dimIDs2 = (/ x_dimID, y_dimID, tl_dimID, nh_dimID /)

  dimIDs4 = (/ x_dimID, y_dimID, tl_dimID, vr_dimID /)
  dimIDs5 = (/ x_dimID, y_dimID, tl_dimID, vr_dimID, nh_dimID /)
  dimIDs6 = (/ x_dimID, y_dimID, tl_dimID, vr_dimID, nh_dimID, sp_dimID /)
  dimIDs7 = (/ x_dimID, y_dimID, tl_dimID, nh_dimID, sp_dimID /)

!DEBUG
! write(*,*) "DIM IDs defined"

!  dimIDs6 = (/ x_dimID2, y_dimID2, tl_dimID2, vr_dimID2, nh_dimID2, sp_dimID2 /)

  ! Define coordinate variables
  call check(nf90_def_var(rsfID, "lon", NF90_REAL, x_dimID, x_varID))    ! returns varID
  call check(nf90_put_att(rsfID, x_varID, "long_name", "longitude"))
  call check(nf90_put_att(rsfID, x_varID, "units", "degrees_east"))
  call check(nf90_put_att(rsfID, x_varID, "standard_name", "longitude"))
  call check(nf90_def_var(rsfID, "lat", NF90_REAL, y_dimID, y_varID))
  call check(nf90_put_att(rsfID, y_varID, "long_name", "latitude"))
  call check(nf90_put_att(rsfID, y_varID, "units", "degrees_north"))
  call check(nf90_put_att(rsfID, y_varID, "standard_name", "latitude"))
  call check(nf90_def_var(rsfID, "tile", NF90_REAL, tl_dimID, tl_varID))
  call check(nf90_put_att(rsfID, tl_varID, "long_name", "tiles"))
  call check(nf90_put_att(rsfID, tl_varID, "units", "none"))
  call check(nf90_put_att(rsfID, tl_varID, "standard_name", "tiles"))
  call check(nf90_def_var(rsfID, "vert", NF90_REAL, vr_dimID, vr_varID))
  call check(nf90_put_att(rsfID, vr_varID, "long_name", "verticalLev"))
  call check(nf90_put_att(rsfID, vr_varID, "units", "none"))
  call check(nf90_put_att(rsfID, vr_varID, "standard_name", "verticalLev"))
  call check(nf90_def_var(rsfID, "nhab", NF90_REAL, nh_dimID, nh_varID))
  call check(nf90_put_att(rsfID, nh_varID, "long_name", "nrHabitats"))
  call check(nf90_put_att(rsfID, nh_varID, "units", "none"))
  call check(nf90_put_att(rsfID, nh_varID, "standard_name", "nrHabitats"))
   
  call check(nf90_def_var(rsfID2, "lon", NF90_REAL, x_dimID, x_varID))    ! returns varID
  call check(nf90_put_att(rsfID2, x_varID, "long_name", "longitude"))
  call check(nf90_put_att(rsfID2, x_varID, "units", "degrees_east"))
  call check(nf90_put_att(rsfID2, x_varID, "standard_name", "longitude"))
  call check(nf90_def_var(rsfID2, "lat", NF90_REAL, y_dimID, y_varID))
  call check(nf90_put_att(rsfID2, y_varID, "long_name", "latitude"))
  call check(nf90_put_att(rsfID2, y_varID, "units", "degrees_north"))
  call check(nf90_put_att(rsfID2, y_varID, "standard_name", "latitude"))
  call check(nf90_def_var(rsfID2, "tile", NF90_REAL, tl_dimID, tl_varID))
  call check(nf90_put_att(rsfID2, tl_varID, "long_name", "tiles"))
  call check(nf90_put_att(rsfID2, tl_varID, "units", "none"))
  call check(nf90_put_att(rsfID2, tl_varID, "standard_name", "tiles"))
  call check(nf90_def_var(rsfID2, "vert", NF90_REAL, vr_dimID, vr_varID))
  call check(nf90_put_att(rsfID2, vr_varID, "long_name", "verticalLev"))
  call check(nf90_put_att(rsfID2, vr_varID, "units", "none"))
  call check(nf90_put_att(rsfID2, vr_varID, "standard_name", "verticalLev"))
  call check(nf90_def_var(rsfID2, "nhab", NF90_REAL, nh_dimID, nh_varID))
  call check(nf90_put_att(rsfID2, nh_varID, "long_name", "nrHabitats"))
  call check(nf90_put_att(rsfID2, nh_varID, "units", "none"))
  call check(nf90_put_att(rsfID2, nh_varID, "standard_name", "nrHabitats"))
  call check(nf90_def_var(rsfID2, "spec", NF90_REAL, sp_dimID, sp_varID))
  call check(nf90_put_att(rsfID2, sp_varID, "long_name", "nrSpecies"))
  call check(nf90_put_att(rsfID2, sp_varID, "units", "none"))
  call check(nf90_put_att(rsfID2, sp_varID, "standard_name", "nrSpecies"))
 
!DEBUG
! write(*,*) "Coordinate vars defined"

  ! Define restart variables

  call def_varG(rsfID, kfile_restartL+1, rsvarIDR(1), 0 )
  call def_varG(rsfID, kfile_restartL+2, rsvarIDR(2), 5 )
  call def_varG(rsfID, kfile_restartL+3, rsvarIDR(3), 1 )
  call def_varG(rsfID, kfile_restartL+4, rsvarIDR(4), 1 )
  call def_varG(rsfID, kfile_restartL+5, rsvarIDR(5), 4 )
  call def_varG(rsfID, kfile_restartL+6, rsvarIDR(6), 1 )
  call def_varG(rsfID, kfile_restartL+7, rsvarIDR(7), 2 )

  call def_varG(rsfID2, kfile_restartV+1, rsvarIDR2(1), 6 )
!DEBUG
!  call def_varG(rsfID2, kfile_restartV+1, rsvarIDR2(1), 2 )
  call def_varG(rsfID2, kfile_restartV+2, rsvarIDR2(2), 6 )
  call def_varG(rsfID2, kfile_restartV+3, rsvarIDR2(3), 7 )

  call check(nf90_enddef(rsfID)) !End Definitions  
  call check(nf90_enddef(rsfID2))

!DEBUG
! write(*,*) "End definitions"

  ! Write grid variables
  call check(nf90_put_var(rsfID, x_varID, xpos))
  call check(nf90_put_var(rsfID, y_varID, ypos))
  call check(nf90_put_var(rsfID2, x_varID, xpos))
  call check(nf90_put_var(rsfID2, y_varID, ypos))

!DEBUG
! write(*,*) "Grid vars written"

endif

! Write restart (variable, code, rsfileID, dimensions)

call write_varG0(rsfID,  rH2Os_g(:),            rsvarIDR(1)) 
call write_varG5(rsfID,  rH2Ol_0(:,:,:,:),      rsvarIDR(2)) 
call write_varG1(rsfID,  rH2Ol_g1(:,:),         rsvarIDR(3)) 
call write_varG1(rsfID,  rH2Ol_g2(:,:),         rsvarIDR(4)) 
call write_varG4(rsfID,  atN_fH2Ol_xd(:,:,:),   rsvarIDR(5)) 
call write_varG1(rsfID,  fH2Ol_ug2(:,:),        rsvarIDR(6)) 
call write_varG2(rsfID,  xT_g0(:,:,:),          rsvarIDR(7)) 
!DEBUG
!call write_varG2(rsfID,  xT_g0B(:,:,:),          rsvarIDR(7)) 

!DEBUG
! if (rank .eq. 0) write(*,*) "land vars written"


call write_varG6(rsfID2,  rH2Ol_t(:,:,:,:,:),   rsvarIDR2(1)) 
!DEBUG
!call write_varG2(rsfID2,  xT_g0(:,:,:),         rsvarIDR2(1)) 
call write_varG6(rsfID2,  areaTH_s(:,:,:,:,:),  rsvarIDR2(2)) 
call write_varG7(rsfID2,  xT_g(:,:,:,:),        rsvarIDR2(3)) 

!DEBUG
! if (rank .eq. 0) write(*,*) "NVP vars written"


if (rank .eq. 0) call check(nf90_close(rsfID))
if (rank .eq. 0) call check(nf90_close(rsfID2))


!DEBUG
! write(*,*) "RANK", rank
! write(*,*) "LIBRY write restart -- FINISHED"

!!! ONLY WORKS WITH ONE PROCESSOR!!!
!write(*,*) "SNOW  -- WRITE RESTART BEFORE DISTRIB"
!write(*,*) " "
!write(*,*) "cell   80", rH2Os_g(80)
!write(*,*) "cell    1", rH2Os_g(1)
!write(*,*) "cell  345", rH2Os_g(345)
!write(*,*) "cell  871", rH2Os_g(871)
!write(*,*) "cell 1450", rH2Os_g(1450)
!write(*,*) " "
!write(*,*) " "
!write(*,*) " "
!write(*,*) "AREA fractions in cell 80:"
!write(*,*) " "
!do x1 = 1, p_ntiles
!  write(*,*) "TILE ",x1
!  do x2 = 1, p_nspec
!    write(*,*) "hab 1", areaTH_s(80,x1,1,1,x2)
!    write(*,*) "hab 2", areaTH_s(80,x1,1,2,x2)
!    write(*,*) "level 2", areaTH_s(80,x1,2,1,x2)
!    write(*,*) "hab 2", areaTH_s(80,x1,2,2,x2)
!    write(*,*) " "
!  enddo
!  write(*,*) " "
!enddo


return
end subroutine libry_write_restart




! WRITE SPECIES PROPERTIES
!---------------------------------------------------------------

subroutine libry_outspecG ()
use libry_par
use libry_global2
use netcdf
use libry_nc
implicit none

integer :: x_dimID, y_dimID, t_dimID
integer :: x_varID, y_varID

integer :: k,l,i0,s,t,x
integer :: kfile_o_spec, kfile_habitat, kfilex

! Open output file
if (rank .eq. 0) then

!DEBUG
write(*,*) "LIBRY OUTSPEC G -- STARTED"

!  if (out1) then

    ! Create global output file
    call check(nf90_create(sfile_outputGS, NF90_CLOBBER, outIDS))

    ! Define dimensions.
    call check(nf90_def_dim(outIDS, "lon", nx, x_dimID))    ! returns dimID
    call check(nf90_def_dim(outIDS, "lat", ny, y_dimID))
    call check(nf90_def_dim(outIDS, "time", NF90_UNLIMITED, t_dimID))
    dimIDs3 = (/ x_dimID, y_dimID, t_dimID /)
   
    ! Define coordinate variables
    call check(nf90_def_var(outIDS, "lon", NF90_REAL, x_dimID, x_varID))    ! returns varID
    call check(nf90_put_att(outIDS, x_varID, "long_name", "longitude"))
    call check(nf90_put_att(outIDS, x_varID, "units", "degrees_east"))
    call check(nf90_put_att(outIDS, x_varID, "standard_name", "longitude"))
    call check(nf90_def_var(outIDS, "lat", NF90_REAL, y_dimID, y_varID))
    call check(nf90_put_att(outIDS, y_varID, "long_name", "latitude"))
    call check(nf90_put_att(outIDS, y_varID, "units", "degrees_north"))
    call check(nf90_put_att(outIDS, y_varID, "standard_name", "latitude"))
   
    call check(nf90_def_var(outIDS, "time", NF90_INT, t_dimID, t_varID))
    call check(nf90_put_att(outIDS, t_varID, "units", "day as %Y%m%d.%f"))
    call check(nf90_put_att(outIDS, t_varID, "calendar", "proleptic_gregorian"))

    ! Define output variables
    call def_varG(outIDS, kfile_count_spec, outvarIDS(1), 3)

    i0 = 1

    do l = 1,p_nhabA

      i0 = i0 + 1

      kfile_o_spec = kfile_count_spec2 * l

      call def_varG(outIDS, kfile_o_spec, outvarIDS(i0), 3)

      do s = 1,6
        do k = 1,p_nspecpar

          kfile_o_spec = kfile_count_spec2 * l + s * 100 + k

          i0 = i0 + 1

          call def_varG(outIDS, kfile_o_spec, outvarIDS(i0), 3)
        enddo
      enddo
    enddo

    do x = 1,3
      do t = 1,p_ntiles

        i0 = i0 + 1

        kfilex = kfile_areas + x * 100 + t

        call def_varG(outIDS, kfilex, outvarIDS(i0), 3)
      enddo
    enddo

    call check(nf90_enddef(outIDS)) !End Definitions  

    ! Write grid variables
    call check(nf90_put_var(outIDS, x_varID, xpos))
    call check(nf90_put_var(outIDS, y_varID, ypos))

!  else
!
!    ! Open output file for next write
!    call check(nf90_open(sfile_outputGS, NF90_WRITE, outIDS))
!  endif
endif

! write species properties
call write_varG(outIDS, count_spec(:), outvarIDS(1))

i0 = 1

do l = 1,p_nhabA

  i0 = i0 + 1

  call write_varG(outIDS, count_spec_h(l,:), outvarIDS(i0))

  do s = 1,6
    do k = 1,p_nspecpar

      i0 = i0 + 1

      call write_varG(outIDS, vec_o_avg(s,k,l,:), outvarIDS(i0))
    enddo
  enddo
enddo

do x = 1,3
  do t = 1,p_ntiles

    i0 = i0 + 1

    if (x .eq. 1) then
      call write_varG(outIDS, frac_tile(:,t), outvarIDS(i0))
    elseif (x .eq. 2) then
      call write_varG(outIDS, Acano(:,t), outvarIDS(i0))
    else
      if (t .eq. 1) then
        call write_varG(outIDS, LAIforest(:), outvarIDS(i0))
      elseif (t .eq. 2) then
        call write_varG(outIDS, LAIgrass(:), outvarIDS(i0))
      else
        call write_varG(outIDS, LAIgrass(:)-LAIgrass(:), outvarIDS(i0)) ! quick fix
      endif
    endif
  enddo
enddo

if (rank .eq. 0) call check(nf90_close(outIDS))


!DEBUG
if (rank .eq. 0) write(*,*) "LIBRY OUTSPEC G -- FINISHED"

return
end subroutine libry_outspecG


! CLOSE GLOBAL FILES
!---------------------------------------------------------------

subroutine libry_closeG ()
use libry_par
use netcdf
use libry_nc
implicit none

! Close climate data files
call check(nf90_close(ncID0(1)))
call check(nf90_close(ncID0(2)))
call check(nf90_close(ncID0(3)))
call check(nf90_close(ncID0(4)))
call check(nf90_close(ncID0(5)))
call check(nf90_close(ncID0(6)))
call check(nf90_close(ncID0(7)))

return
end subroutine libry_closeG


! DEALLOCATE GLOBAL FIELDS
!---------------------------------------------------------------

subroutine dealloc_global ()
use libry_par
use mpi
implicit none

integer                 :: mperr

if (rank .eq. 0) then
  deallocate( indata, &
  indataL, &
  bcdata, &
  bcdataL, &
  bcdata2, &
  bcdata2L, &
  rsdata1, &
  rsdata2, &
  rsdata3, &
  rsdata4, &
  rsdata5, &
  rsdata6, &
  rsdata7, &
  rsdata1L, &
  rsdata2L, &
  rsdata3L, &
  rsdata4L, &
  rsdata5L, &
  rsdata6L, &
  rsdata7L, &
  rsdataV1, &
  rsdataV2, &
  rsdataV3, &
  rsdataV1L, &
  rsdataV2L, &
  rsdataV3L, &
  xpos, &
  ypos, &
  indvec )
endif

!DEBUG
!write(*,*) "Proc:  ",rank, "  LAI",size(laidata)
!write(*,*) "Proc:  ",rank, "  SAI",size(saidata)
!write(*,*) "Proc:  ",rank, "  SSA",size(ssadata)
!write(*,*) "Proc:  ",rank, "biome",size(biomedata)
!write(*,*) "Proc:  ",rank, "  IN7",size(indata7)
!write(*,*) "Proc:  ",rank, "ppvec",size(ppvec)
!write(*,*) "Proc:  ",rank, " xpos",size(xpos)
!write(*,*) "Proc:  ",rank, " ypos",size(ypos)
!write(*,*) "Proc:  ",rank, " indv",size(indvec)
!write(*,*) "       "

deallocate( laidata, &
saidata, &
ssadata, &
biomedata, &
indata7, &
restartdata1, &
restartdata2, &
restartdata3, &
restartdata4, &
restartdata5, &
restartdata6, &
restartdata7, &
restartdataV1, &
restartdataV2, &
restartdataV3, &
ppvec )

! Stop MPI
call MPI_FINALIZE(mperr)

return
end subroutine dealloc_global

end module libry_global

