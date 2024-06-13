
!#######################################################################
! LIBRY_LOCAL
!#######################################################################

module libry_local
contains

! INQ_NTLOC -- INQUIRE LENGTH OF INPUT DATA
!-----------------------------------------------------------------------

subroutine inq_ntLoc ()
use libry_par
implicit none

integer                 :: stat
real, allocatable       :: dummy(:)

allocate( dummy(nCPts) )

stat = 0
nt = 0

open( ktairfile, file=stairfile, status='old', action='read', iostat=stat )

do while ( stat .ge. 0 )

  read( ktairfile, *, iostat=stat ) dummy(:)

  if ( stat .ge. 0 ) nt = nt + 1

enddo

close ( ktairfile )

deallocate ( dummy )

end subroutine inq_ntLoc


! LIBRY_READBCL -- READ BOUNDARY CONDITIONS
!-----------------------------------------------------------------------

subroutine libry_readBCl ()
use libry_par
implicit none

integer :: stat,j,m

integer, dimension(4) :: stat2

! read canopy and soil areas

open( klaifile, file=slaifile, status='old', action='read', iostat=stat2(1) )
open( ksaifile, file=ssaifile, status='old', action='read', iostat=stat2(2) )
open( kssafile, file=sssafile, status='old', action='read', iostat=stat2(3) )
open( kbiomefile, file=sbiomefile, status='old', action='read', iostat=stat2(4) )

if ( any( stat2 .ne. 0 ) ) then
  write(*,*) "ERROR opening ecosystem files"
  stop
endif

do m = 1, 12
  read( klaifile, * ) areaLEAF(:,m) ! read monthly leaf area index
  read( ksaifile, * ) areaSTEM(:,m) ! read monthly stem area index
enddo

read( kssafile, * ) areaSOIL(:) ! read available soil area
read( kbiomefile, * ) biome(:)  ! read biome indices

close( klaifile )
close( ksaifile )
close( kssafile )
close( kbiomefile )

return
end subroutine libry_readBCl


! LIBRY_OPENL -- OPEN FORCING FILES
!-----------------------------------------------------------------------

subroutine libry_openL ()
use libry_par
implicit none

integer, dimension(7) :: stat3

! open climate forcing files

open( ksradfile, file=ssradfile, status='old', action='read', iostat=stat3(1) )
open( klradfile, file=slradfile, status='old', action='read', iostat=stat3(2) )
open( krainfile, file=srainfile, status='old', action='read', iostat=stat3(3) )
open( ksnowfile, file=ssnowfile, status='old', action='read', iostat=stat3(4) )
open( krhumfile, file=srhumfile, status='old', action='read', iostat=stat3(5) )
open( ktairfile, file=stairfile, status='old', action='read', iostat=stat3(6) )
open( kwindfile, file=swindfile, status='old', action='read', iostat=stat3(7) )

if ( any( stat3 .ne. 0 ) ) then
  write(*,*) "ERROR opening climate forcing files"
  stop
endif

return
end subroutine libry_openL


! LIBRY_READHOURL -- READ HOURLY CLIMATE INPUT
!-----------------------------------------------------------------------

subroutine libry_readHourL ()
use libry_par
implicit none

integer :: stat

! read( ksradfile, '('//fd1//'E9.3)', iostat=stat ) fRADs_as(:)
  read( ksradfile, *, iostat=stat ) fRADs_ad(:)
  if ( stat .lt. 0 ) then
    rewind( ksradfile )
    read( ksradfile, *, iostat=stat ) fRADs_ad(:)
  endif
  read( klradfile, *, iostat=stat ) fRADl_ad(:)
  if ( stat .lt. 0 ) then
    rewind( klradfile )
    read( klradfile, *, iostat=stat ) fRADl_ad(:)
  endif
  read( krainfile, *, iostat=stat ) fH2Ol_ad(:)
  if ( stat .lt. 0 ) then
    rewind( krainfile )
    read( krainfile, *, iostat=stat ) fH2Ol_ad(:)
  endif
  read( ksnowfile, *, iostat=stat ) fH2Os_ad(:)
  if ( stat .lt. 0 ) then
    rewind( ksnowfile )
    read( ksnowfile, *, iostat=stat ) fH2Os_ad(:)
  endif
  read( krhumfile, *, iostat=stat ) rH2Og_RH(:)
  if ( stat .lt. 0 ) then
    rewind( krhumfile )
    read( krhumfile, *, iostat=stat ) rH2Og_RH(:)
  endif
  read( ktairfile, *, iostat=stat ) xT_a(:)
  if ( stat .lt. 0 ) then
    rewind( ktairfile )
    read( ktairfile, *, iostat=stat ) xT_a(:)
  endif
  read( kwindfile, *, iostat=stat ) fAIR_s(:)
  if ( stat .lt. 0 ) then
    rewind( kwindfile )
    read( kwindfile, *, iostat=stat ) fAIR_s(:)
  endif

return
end subroutine libry_readHourL

! LIBRY_OPENOUTL -- OPEN OUTPUT FILES
!-----------------------------------------------------------------------

subroutine libry_openOutL ()
use libry_par
implicit none

integer :: i,k,l
integer :: kfile_diurnalgc

character(9) :: ix

character(30) :: sfile_habitat

! open output files

call system('mkdir -p '//fluxdir//'')
call system('mkdir -p '//specdir//'')

do k = 1,2

  select case( k )
    case( 1 )
      sfile_habitat = trim(sfile_ground)
    case( 2 )
      sfile_habitat = trim(sfile_canopy)
    case default
      sfile_habitat = "Dummy0"
  end select

  open(k*1000+kfile_rCO2d          , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_rCO2d      ,  status='replace',    action='write'  )
  open(k*1000+kfile_rCb            , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_rCb        ,  status='replace',    action='write'  )
  open(k*1000+kfile_rH2Ol          , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_rH2Ol      ,  status='replace',    action='write'  )
  open(k*1000+kfile_rH2Os          , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_rH2Os      ,  status='replace',    action='write'  )
  open(k*1000+kfile_rmaxH2Ol       , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_rmaxH2Ol   ,  status='replace',    action='write'  )
  open(k*1000+kfile_area           , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_area       ,  status='replace',    action='write'  )
  open(k*1000+kfile_act            , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_act        ,  status='replace',    action='write'  )
  open(k*1000+kfile_actB           , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_actB       ,  status='replace',    action='write'  )
  open(k*1000+kfile_fCO2gc         , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fCO2gc     ,  status='replace',    action='write'  )
  open(k*1000+kfile_fCcg           , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fCcg       ,  status='replace',    action='write'  )
  open(k*1000+kfile_fCcb           , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fCcb       ,  status='replace',    action='write'  )
  open(k*1000+kfile_fCcb_l         , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fCcb_l     ,  status='replace',    action='write'  )
  open(k*1000+kfile_fCcb_c         , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fCcb_c     ,  status='replace',    action='write'  )

  open(k*1000+kfile_fCbo           , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fCbo       ,  status='replace',    action='write'  )
  open(k*1000+kfile_frgr           , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_frgr       ,  status='replace',    action='write'  )
  open(k*1000+kfile_fH2Ol_ux       , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fH2Ol_ux   ,  status='replace',    action='write'  )
  open(k*1000+kfile_fH2Ol_xd       , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fH2Ol_xd   ,  status='replace',    action='write'  )
  open(k*1000+kfile_fH2Ogl_ux      , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fH2Ogl_ux  ,  status='replace',    action='write'  )
  open(k*1000+kfile_fH2Olg_xu      , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fH2Olg_xu  ,  status='replace',    action='write'  )
  open(k*1000+kfile_fH2Ol_bx       , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fH2Ol_bx   ,  status='replace',    action='write'  )
  open(k*1000+kfile_fH2Ol_bd       , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_fH2Ol_bd   ,  status='replace',    action='write'  )
  open(k*1000+kfile_Ts             , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_Ts         ,  status='replace',    action='write'  )
  open(k*1000+kfile_H              , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_H          ,  status='replace',    action='write'  )
  open(k*1000+kfile_E              , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_E          ,  status='replace',    action='write'  )
  open(k*1000+kfile_C              , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_C          ,  status='replace',    action='write'  )
  open(k*1000+kfile_EB             , file=fluxdir//trim(adjustl(sfile_habitat))//sfile_EB         ,  status='replace',    action='write'  )
enddo

open( kfile_fH2Ol_ug2   , file=fluxdir//sfile_fH2Ol_ug2   ,  status='replace',    action='write'  )
open( kfile_rH2Ol_g1    , file=fluxdir//sfile_rH2Ol_g1    ,  status='replace',    action='write'  )
open( kfile_rH2Ol_g2    , file=fluxdir//sfile_rH2Ol_g2    ,  status='replace',    action='write'  )
open( kfile_rH2Os_g     , file=fluxdir//sfile_rH2Os_g     ,  status='replace',    action='write'  )
open( kfile_fH2Ol_ug    , file=fluxdir//sfile_fH2Ol_ug    ,  status='replace',    action='write'  )
open( kfile_fH2Ol_go    , file=fluxdir//sfile_fH2Ol_go    ,  status='replace',    action='write'  )
open( kfile_fH2Ol_gb    , file=fluxdir//sfile_fH2Ol_gb    ,  status='replace',    action='write'  )
open( kfile_fH2Olg_ga   , file=fluxdir//sfile_fH2Olg_ga   ,  status='replace',    action='write'  )
open( kfile_fH2Osl_g    , file=fluxdir//sfile_fH2Osl_g    ,  status='replace',    action='write'  )
open( kfile_fH2Os_ad    , file=fluxdir//sfile_fH2Os_ad    ,  status='replace',    action='write'  )
open( kfile_fH2Ol_ad    , file=fluxdir//sfile_fH2Ol_ad    ,  status='replace',    action='write'  )
open( kfile_xT_a        , file=fluxdir//sfile_xT_a        ,  status='replace',    action='write'  )
open( kfile_Tg          , file=fluxdir//sfile_Tg          ,  status='replace',    action='write'  )
open( kfile_G           , file=fluxdir//sfile_G           ,  status='replace',    action='write'  )
open( kfile_fRADs       , file=fluxdir//sfile_fRADs       ,  status='replace',    action='write'  )
!
if (BSCtypes) then
  open( kfile_MareaLCgS     , file=fluxdir//sfile_MareaLCgS   ,  status='replace',    action='write'  )
  open( kfile_MareaDCgS     , file=fluxdir//sfile_MareaDCgS   ,  status='replace',    action='write'  )
  open( kfile_MareaCCgS     , file=fluxdir//sfile_MareaCCgS   ,  status='replace',    action='write'  )
  open( kfile_MareaMCgS     , file=fluxdir//sfile_MareaMCgS   ,  status='replace',    action='write'  )
  open( kfile_MrH2OlLC_S    , file=fluxdir//sfile_MrH2OlLC_S  ,  status='replace',    action='write'  )
  open( kfile_MrH2OlDC_S    , file=fluxdir//sfile_MrH2OlDC_S  ,  status='replace',    action='write'  )
  open( kfile_MrH2OlCC_S    , file=fluxdir//sfile_MrH2OlCC_S  ,  status='replace',    action='write'  )
  open( kfile_MrH2OlMC_S    , file=fluxdir//sfile_MrH2OlMC_S  ,  status='replace',    action='write'  )
  open( kfile_MTsLC_S       , file=fluxdir//sfile_MTsLC_S     ,  status='replace',    action='write'  )
  open( kfile_MTsDC_S       , file=fluxdir//sfile_MTsDC_S     ,  status='replace',    action='write'  )
  open( kfile_MTsCC_S       , file=fluxdir//sfile_MTsCC_S     ,  status='replace',    action='write'  )
  open( kfile_MTsMC_S       , file=fluxdir//sfile_MTsMC_S     ,  status='replace',    action='write'  )
  open( kfile_MfCcbLC_S     , file=fluxdir//sfile_MfCcbLC_S   ,  status='replace',    action='write'  )
  open( kfile_MfCcbDC_S     , file=fluxdir//sfile_MfCcbDC_S   ,  status='replace',    action='write'  )
  open( kfile_MfCcbCC_S     , file=fluxdir//sfile_MfCcbCC_S   ,  status='replace',    action='write'  )
  open( kfile_MfCcbMC_S     , file=fluxdir//sfile_MfCcbMC_S   ,  status='replace',    action='write'  )

  if (NOHONO) then
    open( kfile_MfNO_N      , file=fluxdir//sfile_MfNO_N      ,  status='replace',    action='write'  )
    open( kfile_MfHONO_N    , file=fluxdir//sfile_MfHONO_N    ,  status='replace',    action='write'  )
  endif
endif

return
end subroutine libry_openOutL


! LIBRY_OUTPUTL -- WRITE OUTPUT
!-----------------------------------------------------------------------

subroutine libry_outputL ()
use libry_par
implicit none

integer :: k

character(9) :: fd1

write( fd1 ,'(I9)' ) nCPts

! write averaged fluxes and reservoirs

do k = 1,2

  write(k*1000+kfile_rCO2d          , '('//fd1//'(E9.3,2X))') ag_rCO2d(:,k)
  write(k*1000+kfile_rCb            , '('//fd1//'(E9.3,2X))') ag_rCb(:,k)
  write(k*1000+kfile_rH2Ol          , '('//fd1//'(E9.3,2X))') ag_rH2Ol(:,k)
  write(k*1000+kfile_rH2Os          , '('//fd1//'(E9.3,2X))') ag_rH2Os(:,k)
  write(k*1000+kfile_rmaxH2Ol       , '('//fd1//'(E9.3,2X))') ag_rmaxH2Ol(:,k)
  write(k*1000+kfile_area           , '('//fd1//'(E9.3,2X))') ag_areaTH_s(:,k)
  write(k*1000+kfile_act            , '('//fd1//'(E9.3,2X))') ag_act(:,k)
  write(k*1000+kfile_actB           , '('//fd1//'(E9.3,2X))') ag_actB(:,k)
  write(k*1000+kfile_fCO2gc         , '('//fd1//'(E9.3,2X))') ag_fCO2gc(:,k)
  write(k*1000+kfile_fCcg           , '('//fd1//'(E9.3,2X))') ag_fCcg(:,k)
  write(k*1000+kfile_fCcb           , '('//fd1//'(E9.3,2X))') ag_fCcb(:,k)
  write(k*1000+kfile_fCcb_l         , '('//fd1//'(E9.3,2X))') ag_fCcb_l(:,k)
  write(k*1000+kfile_fCcb_c         , '('//fd1//'(E9.3,2X))') ag_fCcb_c(:,k)

  write(k*1000+kfile_fCbo           , '('//fd1//'(E9.3,2X))') ag_fCbo(:,k)
  write(k*1000+kfile_frgr           , '('//fd1//'(E9.3,2X))') ag_frgr(:,k)
  write(k*1000+kfile_fH2Ol_ux       , '('//fd1//'(E9.3,2X))') ag_fH2Ol_ux(:,k)
  write(k*1000+kfile_fH2Ol_xd       , '('//fd1//'(E9.3,2X))') ag_fH2Ol_xd(:,k)
  write(k*1000+kfile_fH2Ogl_ux      , '('//fd1//'(E9.3,2X))') ag_fH2Ogl_ux(:,k)
  write(k*1000+kfile_fH2Olg_xu      , '('//fd1//'(E9.3,2X))') ag_fH2Olg_xu(:,k)
  write(k*1000+kfile_fH2Ol_bx       , '('//fd1//'(E9.3,2X))') ag_fH2Ol_bx(:,k)
  write(k*1000+kfile_fH2Ol_bd       , '('//fd1//'(E9.3,2X))') ag_fH2Ol_bd(:,k)
  write(k*1000+kfile_Ts             , '('//fd1//'(E9.3,2X))') ag_Ts(:,k)
  write(k*1000+kfile_H              , '('//fd1//'(E9.3,2X))') ag_H(:,k)
  write(k*1000+kfile_E              , '('//fd1//'(E9.3,2X))') ag_E(:,k)
  write(k*1000+kfile_C              , '('//fd1//'(E9.3,2X))') ag_C(:,k)
  write(k*1000+kfile_EB             , '('//fd1//'(E9.3,2X))') ag_EB(:,k)
enddo

write( kfile_fH2Ol_ug2    , '('//fd1//'(E9.3,2X))') ag_fH2Ol_ug2(:)
write( kfile_rH2Ol_g1     , '('//fd1//'(E9.3,2X))') ag_rH2Ol_g1(:)
write( kfile_rH2Ol_g2     , '('//fd1//'(E9.3,2X))') ag_rH2Ol_g2(:)
write( kfile_rH2Os_g      , '('//fd1//'(E9.3,2X))') ag_rH2Os_g(:)
write( kfile_fH2Ol_ug     , '('//fd1//'(E9.3,2X))') ag_fH2Ol_ug(:)
write( kfile_fH2Ol_go     , '('//fd1//'(E9.3,2X))') ag_fH2Ol_go(:)
write( kfile_fH2Ol_gb     , '('//fd1//'(E9.3,2X))') ag_fH2Ol_gb(:)
write( kfile_fH2Olg_ga    , '('//fd1//'(E9.3,2X))') ag_fH2Olg_ga(:)
write( kfile_fH2Osl_g     , '('//fd1//'(E9.3,2X))') ag_fH2Osl_g(:)
write( kfile_fH2Os_ad     , '('//fd1//'(E9.3,2X))') ag_fH2Os_ad(:)
write( kfile_fH2Ol_ad     , '('//fd1//'(E9.3,2X))') ag_fH2Ol_ad(:)
write( kfile_xT_a         , '('//fd1//'(E9.3,2X))') ag_xT_a(:)
write( kfile_Tg           , '('//fd1//'(E9.3,2X))') ag_Tg(:)
write( kfile_G            , '('//fd1//'(E9.3,2X))') ag_G(:)
write( kfile_fRADs        , '('//fd1//'(E9.3,2X))') ag_fRADs(:)
!
if (BSCtypes) then
  write( kfile_MareaLCgS    , '('//fd1//'(E9.3,2X))') a_MareaTHLC_g_S(:)
  write( kfile_MareaDCgS    , '('//fd1//'(E9.3,2X))') a_MareaTHDC_g_S(:)
  write( kfile_MareaCCgS    , '('//fd1//'(E9.3,2X))') a_MareaTHCC_g_S(:)
  write( kfile_MareaMCgS    , '('//fd1//'(E9.3,2X))') a_MareaTHMC_g_S(:)
  write( kfile_MrH2OlLC_S   , '('//fd1//'(E9.3,2X))') a_MrH2OlLC_S(:)
  write( kfile_MrH2OlDC_S   , '('//fd1//'(E9.3,2X))') a_MrH2OlDC_S(:)
  write( kfile_MrH2OlCC_S   , '('//fd1//'(E9.3,2X))') a_MrH2OlCC_S(:)
  write( kfile_MrH2OlMC_S   , '('//fd1//'(E9.3,2X))') a_MrH2OlMC_S(:)
  write( kfile_MTsLC_S      , '('//fd1//'(E9.3,2X))') a_MTsLC_S(:)
  write( kfile_MTsDC_S      , '('//fd1//'(E9.3,2X))') a_MTsDC_S(:)
  write( kfile_MTsCC_S      , '('//fd1//'(E9.3,2X))') a_MTsCC_S(:)
  write( kfile_MTsMC_S      , '('//fd1//'(E9.3,2X))') a_MTsMC_S(:)
  write( kfile_MfCcbLC_S    , '('//fd1//'(E9.3,2X))') a_MfCcbLC_S(:)
  write( kfile_MfCcbDC_S    , '('//fd1//'(E9.3,2X))') a_MfCcbDC_S(:)
  write( kfile_MfCcbCC_S    , '('//fd1//'(E9.3,2X))') a_MfCcbCC_S(:)
  write( kfile_MfCcbMC_S    , '('//fd1//'(E9.3,2X))') a_MfCcbMC_S(:)

  if (NOHONO) then
    write( kfile_MfNO_N     , '('//fd1//'(E9.3,2X))') a_MfNO_N(:)
    write( kfile_MfHONO_N   , '('//fd1//'(E9.3,2X))') a_MfHONO_N(:)
  endif
endif

return
end subroutine libry_outputL


! LIBRY_OUTSURV -- SURVIVING SPECIES OUTPUT
!-----------------------------------------------------------------------

subroutine libry_outsurv ()
use libry_par
implicit none

integer                                 :: i,t,v,h,g,k
integer, dimension(nCPts,p_nhabA)       :: klifesum
character(9)                            :: fd1

write( fd1 ,'(I9)' ) nCPts * p_nhabA

do i=1,nCPts
  k = 0
  do t=1,p_ntiles-2 ! not all tiles used
    do v=1,p_nvert(t)
      do h=1,p_nhab(t,v)
        k = k + 1
        klifesum(i,k) = sum(klife(i,t,v,h,:))
      end do
    end do
  end do
end do

write( kfile_survspec, '(I9,2X,'//fd1//'(I6,2X))') accts/tspd, ( klifesum(g,:), g=1,nCPts )

return
end subroutine libry_outsurv


! LIBRY_OUTSPECL -- SPECIES PROPERTIES OUTPUT
!-----------------------------------------------------------------------

subroutine libry_outspecL ()
use libry_par
implicit none

integer                                 :: l,k,j,g,s
integer                                 :: i,t,v,h
integer                                 :: kfile_o_spec, kfile_habitat
real                                    :: Asum
real, dimension(nCPts,p_nspec,p_nhabA)  :: weight

character(30)                           :: sfile_habitat
character(30)                           :: sfile_property

character(9)                            :: fd1, fd2

write( fd1 ,'(I9)' ) nCPts
write( fd2 ,'(I9)' ) p_nhabA * nCPts

! species numbers 

open( kfile_count_spec, file=specdir//sfile_count_spec, status='replace', action='write' )

write( kfile_count_spec, '('//fd1//'(E9.3,2X))') count_spec(:)

close( kfile_count_spec )

do l = 1,p_nhabA

  ! species numbers per habitat
  kfile_habitat = kfile_count_spec + l * 100

  select case( l )
    case( 1 )
      sfile_habitat = trim(sfile_forFloor)
    case( 2 )
      sfile_habitat = trim(sfile_forStems)
    case( 3 )
      sfile_habitat = trim(sfile_forLeaves)
    case( 4 )
      sfile_habitat = trim(sfile_grass)
    case default
      sfile_habitat = "Dummy0"
  end select

  open( kfile_habitat, file=specdir//trim(adjustl(sfile_habitat))//sfile_count_spec, &
  status='replace', action='write' )

  write( kfile_habitat, '('//fd1//'(E9.3,2X))') count_spec_h(l,:)

  close( kfile_habitat )

  ! physiological parameters

  do k = 1,p_nspecpar

    select case( k )
      case( 1 )
        sfile_property = trim(sfile_Albedo)
      case( 2 )
        sfile_property = trim(sfile_Height)
      case( 3 )
        sfile_property = trim(sfile_dSublayer)
      case( 4 )
        sfile_property = trim(sfile_FracCells)
      case( 5 )
        sfile_property = trim(sfile_PoreWater)
      case( 6 )
        sfile_property = trim(sfile_MolVcmax)
      case( 7 )
        sfile_property = trim(sfile_MolVomax)
      case( 8 )
        sfile_property = trim(sfile_PhsynCap)
      case( 9 )
        sfile_property = trim(sfile_rCellWall)
      case( 10 )
        sfile_property = trim(sfile_Topt)
      case( 11 )
        sfile_property = trim(sfile_Q10value)
      case default
        sfile_property = "Dummy0"
    end select

    kfile_o_spec = kfile_count_spec + l * 100 + k

    open( kfile_o_spec, file=specdir//trim(adjustl(sfile_habitat))//trim(adjustl(sfile_property)), &
    status='replace', action='write' )

    do s = 1,6
      write( kfile_o_spec, '('//fd1//'(E9.3,2X))') vec_o_avg(s,k,l,:)
    enddo

    close( kfile_o_spec )
  enddo
enddo

! weight of each strategy

do i=1,nCPts
  k = 0
  do t=1,p_ntiles-2 ! not all tiles used
    do v=1,p_nvert(t)
      do h=1,p_nhab(t,v)
        
        k = k + 1

        Asum = sum(areaTH_s(i,t,v,h,:))

        do j = 1,p_nspec

          weight(i,j,k) = areaTH_s(i,t,v,h,j) / Asum

        enddo
      enddo
    enddo
  enddo
enddo

open( kfile_weight, file=specdir//sfile_weightS, status='replace', action='write' )

do j = 1,p_nspec

  write( kfile_weight, '('//fd2//'(E9.3,2X))') ( weight(g,j,:), g=1,nCPts )

enddo

close( kfile_weight )

return
end subroutine libry_outspecL

! LIBRY_CLOSEL -- CLOSE OUTPUT FILES
!-----------------------------------------------------------------------

subroutine libry_closeL ()
use libry_par
implicit none

integer :: i,l
integer :: kfile_diurnalgc

close( ksradfile )
close( klradfile )
close( krainfile )
close( ksnowfile )
close( krhumfile )
close( ktairfile )
close( kwindfile )


close( kfile_rCO2d )
close( kfile_rCb )
close( kfile_rH2Ol )
close( kfile_rH2Os )
close( kfile_rH2Ol_g1 )
close( kfile_rH2Ol_g2 )
close( kfile_rmaxH2Ol )
close( kfile_rH2Os_g )
close( kfile_area )
close( kfile_act )
close( kfile_actB )
close( kfile_fCO2gc )
close( kfile_fCcg )
close( kfile_fCcb )
close( kfile_fCcb_l )
close( kfile_fCcb_c )
close( kfile_fCbo )
close( kfile_frgr )
close( kfile_fH2Ol_ux )
close( kfile_fH2Ol_xd )
close( kfile_fH2Ogl_ux )
close( kfile_fH2Olg_xu )
close( kfile_fH2Ol_bx )
close( kfile_fH2Ol_bd )
close( kfile_fH2Ol_ug )
close( kfile_fH2Ol_go )
close( kfile_fH2Ol_gb )
close( kfile_fH2Olg_ga )
close( kfile_fH2Ol_ug2 )
close( kfile_fH2Osl_g )
close( kfile_fH2Os_ad )
close( kfile_fH2Ol_ad )
close( kfile_Ts )
close( kfile_xT_a )
close( kfile_Tg )
close( kfile_H )
close( kfile_G )
close( kfile_E )
close( kfile_C )
close( kfile_EB )
close( kfile_fRADs )

!
if (BSCtypes) then
  close( kfile_MareaLCgS )
  close( kfile_MareaDCgS )
  close( kfile_MareaCCgS )
  close( kfile_MareaMCgS )
  close( kfile_MrH2OlLC_S )
  close( kfile_MrH2OlDC_S )
  close( kfile_MrH2OlCC_S )
  close( kfile_MrH2OlMC_S )
  close( kfile_MTsLC_S )
  close( kfile_MTsDC_S )
  close( kfile_MTsCC_S )
  close( kfile_MTsMC_S )
  close( kfile_MfCcbLC_S )
  close( kfile_MfCcbDC_S )
  close( kfile_MfCcbCC_S )
  close( kfile_MfCcbMC_S )

  if (NOHONO) then
    close( kfile_MfNO_N )
    close( kfile_MfHONO_N )
  endif
endif

return
end subroutine libry_closeL

end module libry_local

