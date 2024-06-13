
!***********************************************************************
! MAIN PROGRAM
!***********************************************************************

program libry
use libry_par
use libry_common
use libry_interface
use libry_local
use libry_nova
use libry_land
use libry_global
use libry_global2
use libry_opt
implicit none

integer :: month0, day0

! initialisation

call init_mpi ! initialise MPI

if (para) then

  if (rank .eq. 0) call read_land ! read land mask file

  call bc_pp ! broadcast points per processor and other parameters

  nCPts = pp

  if (rank .eq. 0) call libry_namelist ! read namelist
    
  call bc_namelist

  call alloc_global ! allocate global fields     # if (rank .eq. 0) 
else

  call libry_namelist ! read namelist

  nCPts = nSites
endif

call libry_alloc ! allocate variables

if (BSCtypes) call libry_allocBSC

if (rank .eq. 0) call libry_readSpec ! read strategy parameters

if (para) then
  call bc_specpar ! broadcast strategies

  if (rank .eq. 0) call libry_readBCg ! read boundary conditions (global)

  call sc_BC ! scatter boundary conditions
else
  call inq_ntLoc

  call libry_readBCl ! local
endif

if (BSCtypes .and. NOHONO) then
  if (rank .eq. 0) call libry_readNOHONO ! read NO/HONO fluxes (optional)
  
  if (para) call bc_NOHONO
endif

if (para) then

  call libry_init(ppvec(rank+1)) ! initialise variables (global)

  if (lrestart) then

    if (rank .eq. 0) call libry_read_restart(ppvec(rank+1))

    call sc_restart

    call nova_update_init(ppvec(rank+1))
  endif
else
  call libry_init(nCPts) ! local

!  call libry_init_speciesL
endif

call libry_reset(nCPts) ! set accumulation variables to zero

if (para) then
  if (rank .eq. 0) call libry_openG ! open forcing files (global)
else
  call libry_openL ! local
endif

if (.not. para) open(kfile_survspec, file=sfile_survspec, status='replace', action='write') ! open survivors file

year = year0
cyear = cyear0
accts = accts0
tspd = 86400 / tsl

months30(:) = (/ 4, 6, 9, 11 /)

tsindata = tsindata0
tpos = tpos0
out1 = .true.

if (para .and. rank .eq. 0) then
  if (lrestart) then
    open( kstatus, file=sstatus, status='old', action='write', position='append' )
  else
    open( kstatus, file=sstatus, status='replace', action='write' )
  endif
endif

if (.not. para) runperiod = lastyear

! main yearly loop

do while (( year .le. year0+runperiod-1 ) .and. ( year .le. lastyear ))

  if ( mod( cyear, 4 ) .eq. 0 ) then ! determine leap year
    leapyear = .true.
    if ( mod( cyear, 100 ) .eq. 0 .and. mod( cyear, 400 ) .ne. 0 ) leapyear = .false.
  else
    leapyear = .false.
  endif

  if ( year .ge. yearout1 .and. year .le. yearoutX ) then ! check for writing output
    writeout = .true.
  else
    writeout = .false.
  endif

! monthly loop

  do month0 = 1, 12

    month = month0

    if ( any( months30 .eq. month )) then ! determine days per month
      dpm = 30
    elseif ( month .eq. 2 ) then
      dpm = 28
      if ( leapyear ) dpm = 29
    else
      dpm = 31
    endif

    call libry_readMon(month) ! read monthly forcing (from variable)

! daily loop

    do day0 = 1, dpm

      day = day0

! timestep loop

      do ts = 1, tspd

        if (para) then

          if (rank .eq. 0) call libry_readHourG(tsindata) ! read hourly climate forcing (global)

          call sc_clim ! scatter climate forcing

          fH2Ol_ad(:) = fH2Ol_ad(:) * 0.001 ! adapt units for ERA5 or WATCH data (mm/s -> m/s)
          fH2Os_ad(:) = fH2Os_ad(:) * 0.001 
        else
          call libry_readHourL ! local

          fH2Ol_ad(:) = fH2Ol_ad(:) * 0.001
          fH2Os_ad(:) = fH2Os_ad(:) * 0.001 
        endif

        if (para) then

          call libry_step(ppvec(rank+1)) ! perform calculations
        else
          call libry_step(nCPts)
        endif

        if ( writeout .and. mod( accts, outint ) .eq. 0 ) then ! output

          if (para) then

            call libry_av_output(ppvec(rank+1)) ! average accumulated variables

            call libry_outputG ! write global output

            tpos = tpos + 1
          else

            call libry_av_output(nCPts) ! local

            if (out1) call libry_openOutL

            call libry_outputL
          endif

          if (out1) out1 = .false.

          call libry_reset(nCPts) ! reset accumulated variables
        endif

!        if ( year .eq. lastyear .and. outdiurnal .and. .not. para) call libry_outdnl ! write diurnal output

        accts = accts + 1

        tsindata = tsindata + 1

        if (tsindata .gt. nt) then  ! shift to end of year ?
          
          tsindata = 1

          cyear = cyear0 -1
        endif

      enddo ! time step loop

      if (.not. para) call libry_outsurv ! write number of surviving species

    enddo ! daily loop

  enddo ! monthly loop

  if (para) then
    if (rank .eq. 0) write( kstatus,* ) "year ",year," finished"
  else
    write(*,*) "year ",year," finished"
  endif

  year = year + 1 ! update simulation year

  cyear = cyear + 1 ! update calendar year

enddo ! yearly loop

! output at end of simulation

if (para) then

  call libry_write_restart

  if ( year-1 .eq. lastyear ) then

    call libry_av_species(ppvec(rank+1)) ! average species properties

    tpos = 1 ! rewind

    call libry_outspecG ! write averaged properties
  endif

  if (rank .eq. 0) call libry_closeG ! close forcing
else

  call libry_av_species(nCPts) ! local

!  call libry_av_speciesL

  call libry_outspecL ! write avg properties and weight of each species per habitat

  call libry_closeL ! close forcing and output files

  close( kfile_survspec )
endif

call libry_dealloc

if (BSCtypes) call libry_deallocBSC

if (para) then

  if (rank .eq. 0) then

    lrestart = .true. 

    if ( year-1 .eq. lastyear ) then

      open( krestartX, file=skrestartX, status='replace', action='write' )
      write( krestartX,* ) "Last year: no restart"
      close( krestartX )
    endif

    open( krestart, file=skrestart, status='replace', action='write' )

    write( krestart, '(I9)') year
    write( krestart, '(I9)') cyear
    write( krestart, '(I9)') tsindata
    write( krestart, '(I9)') accts
    write( krestart, '(I9)') tpos
    write( krestart, '(L1)') lrestart

    close( krestart )

  endif

  call dealloc_global ! and stop MPI

  if (rank .eq. 0) then
    write( kstatus,* ) "simulation finished"
    close( kstatus )
  endif
endif

end

