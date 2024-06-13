
module libry_opt
contains

! DEFINE BSC TYPES
!---------------------------------------------------------------

subroutine libry_defBSC (k)
use libry_par
implicit none

integer :: k

! distinguish cyanobacteria, lichens & mosses

! height > 2 mm : moss / lichen                 otherwise : cyanobacterium / cyanobacterium0
! DCO2 > X      : moss                          otherwise : cyanobacterium / lichen
! PScap > Y     : moss / cyanobacterium         otherwise : lichen / moss /cyanobacterium0

o_LC(k)                         = 0.0
o_DC(k)                         = 0.0
o_CC(k)                         = 0.0
o_MC(k)                         = 0.0

if (vec_o(k,2) .lt. 0.15) then        ! 0.15 ~ 2 mm -> height limit for cyanobacteria
  if (vec_o(k,8) .lt. 0.25) then       ! no detailed knowledge about PSCAP -> 50:50
    o_LC(k)                     = 1.0
  else
    o_DC(k)                     = 1.0
  endif
else
  if (vec_o(k,9) .lt. 0.75) then       ! no detailed knowledge about DCO2 -> 50:50
    o_CC(k)                     = 1.0
  else
    o_MC(k)                     = 1.0
  endif
endif

return
end subroutine libry_defBSC


! SUM UP BSC COVER
!---------------------------------------------------------------

subroutine libry_covsBSC (k,k2)
use libry_par
implicit none

integer :: k,k2,j

csumL                       = 0.0
csumD                       = 0.0
csumC                       = 0.0
csumM                       = 0.0

do j = 1,p_nspec

! Check if lichen is alive

  if (klife(k,k2,j) .eq. 1) then

! Sum up cover fractions and heights
    csumL                   = csumL + areaTH_s(k,k2,j)*o_LC(j)
    csumD                   = csumD + areaTH_s(k,k2,j)*o_DC(j)
    csumC                   = csumC + areaTH_s(k,k2,j)*o_CC(j)
    csumM                   = csumM + areaTH_s(k,k2,j)*o_MC(j)

  endif ! check for survival

enddo ! End of loop over all species - II -

return
end subroutine libry_covsBSC


! READ NOHONO FLUXES
!---------------------------------------------------------------

subroutine libry_readNOHONO ()
use libry_par
implicit none

integer :: stat,j

! read NO / HONO fluxes as function of water saturation

open( kNOHONO, file=sNOHONO, status='old', action='read', iostat=stat )

if ( stat .ne. 0 ) then
  write(*,*) "ERROR opening NO / HONO emissions file"
  stop
endif

do j = 1,100
  read( kNOHONO, * ) NO_HONO_fSat(j,:)
enddo

close( kNOHONO )

return
end subroutine libry_readNOHONO


! BROADCAST NOHONO DATA
!---------------------------------------------------------------

subroutine bc_NOHONO ()
use libry_par
use mpi
implicit none

integer :: mperr, k

do k = 1,8

  call MPI_BCAST( NO_HONO_fSat(:,k), 100, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mperr )
enddo

return
end subroutine bc_NOHONO


! ALLOCATE BSC AND NO/HONO VARIABLES
!---------------------------------------------------------------

subroutine libry_allocBSC ()
use libry_par
implicit none

allocate( o_LC(p_nspec), &
o_DC(p_nspec), &
o_CC(p_nspec), &
o_MC(p_nspec), &
!
a_areaTHLC_g(nCPts,p_nhab,p_nspec), &
a_areaTHDC_g(nCPts,p_nhab,p_nspec), &
a_areaTHCC_g(nCPts,p_nhab,p_nspec), &
a_areaTHMC_g(nCPts,p_nhab,p_nspec), &
a_rH2OlLC(nCPts,p_nhab,p_nspec), &
a_rH2OlDC(nCPts,p_nhab,p_nspec), &
a_rH2OlCC(nCPts,p_nhab,p_nspec), &
a_rH2OlMC(nCPts,p_nhab,p_nspec), &
a_TsLC(nCPts,p_nhab,p_nspec), &
a_TsDC(nCPts,p_nhab,p_nspec), &
a_TsCC(nCPts,p_nhab,p_nspec), &
a_TsMC(nCPts,p_nhab,p_nspec), &
a_fCcbLC(nCPts,p_nhab,p_nspec), &
a_fCcbDC(nCPts,p_nhab,p_nspec), &
a_fCcbCC(nCPts,p_nhab,p_nspec), &
a_fCcbMC(nCPts,p_nhab,p_nspec), &
!
a_MareaTHLC_g_S(nCPts), &
a_MareaTHDC_g_S(nCPts), &
a_MareaTHCC_g_S(nCPts), &
a_MareaTHMC_g_S(nCPts), &
a_MrH2OlLC_S(nCPts), &
a_MrH2OlDC_S(nCPts), &
a_MrH2OlCC_S(nCPts), &
a_MrH2OlMC_S(nCPts), &
a_MTsLC_S(nCPts), &
a_MTsDC_S(nCPts), &
a_MTsCC_S(nCPts), &
a_MTsMC_S(nCPts), &
a_MfCcbLC_S(nCPts), &
a_MfCcbDC_S(nCPts), &
a_MfCcbCC_S(nCPts), &
a_MfCcbMC_S(nCPts))

if (NOHONO) then
  allocate( NO_HONO_fSat(100,8), &
  fNO_N(p_nspec), &
  fHONO_N(p_nspec), &
  a_fNO_N(nCPts,p_nhab,p_nspec), &
  a_fHONO_N(nCPts,p_nhab,p_nspec), &
  a_MfNO_N(nCPts), &
  a_MfHONO_N(nCPts))
endif

return
end subroutine libry_allocBSC


! DEALLOCATE BSC AND NO/HONO VARIABLES
!---------------------------------------------------------------

subroutine libry_deallocBSC ()
use libry_par
implicit none

deallocate( o_LC, &
o_DC, &
o_CC, &
o_MC, &
!
a_areaTHLC_g, &
a_areaTHDC_g, &
a_areaTHCC_g, &
a_areaTHMC_g, &
a_rH2OlLC, &
a_rH2OlDC, &
a_rH2OlCC, &
a_rH2OlMC, &
a_TsLC, &
a_TsDC, &
a_TsCC, &
a_TsMC, &
a_fCcbLC, &
a_fCcbDC, &
a_fCcbCC, &
a_fCcbMC, &
!
a_MareaTHLC_g_S, &
a_MareaTHDC_g_S, &
a_MareaTHCC_g_S, &
a_MareaTHMC_g_S, &
a_MrH2OlLC_S, &
a_MrH2OlDC_S, &
a_MrH2OlCC_S, &
a_MrH2OlMC_S, &
a_MTsLC_S, &
a_MTsDC_S, &
a_MTsCC_S, &
a_MTsMC_S, &
a_MfCcbLC_S, &
a_MfCcbDC_S, &
a_MfCcbCC_S, &
a_MfCcbMC_S)

if (NOHONO) then
  deallocate( NO_HONO_fSat, &
  fNO_N, &
  fHONO_N, &
  a_fNO_N, &
  a_fHONO_N, &
  a_MfNO_N, &
  a_MfHONO_N)
endif

return
end subroutine libry_deallocBSC


! CALCULATE NO AND HONO EMISSIONS
!---------------------------------------------------------------

subroutine libry_fNOHONO (k,st,ac)
use libry_par
implicit none

integer :: k
real    :: st,ac

! NO / HONO Emissions [ ng / m2 / s ]

if (nint(st*100.0) .gt. 0 .and. ac .gt. 0.0) then

  if (o_LC(k) .eq. 1.0) then

    fHONO_N(k)          = NO_HONO_fSat(nint(st*100.0),1)
    fNO_N(k)            = NO_HONO_fSat(nint(st*100.0),5)

  elseif (o_DC(k) .eq. 1.0) then

    fHONO_N(k)          = NO_HONO_fSat(nint(st*100.0),2)
    fNO_N(k)            = NO_HONO_fSat(nint(st*100.0),6)

  elseif (o_CC(k) .eq. 1.0) then

    fHONO_N(k)          = NO_HONO_fSat(nint(st*100.0),3)
    fNO_N(k)            = NO_HONO_fSat(nint(st*100.0),7)
  else

    fHONO_N(k)          = NO_HONO_fSat(nint(st*100.0),4)
    fNO_N(k)            = NO_HONO_fSat(nint(st*100.0),8)
  endif

  fHONO_N(k)            = fHONO_N(k) * 2.0**((xT_s(k)-298.15)/10.0)
  fNO_N(k)              = fNO_N(k) * 2.0**((xT_s(k)-298.15)/10.0)
else

  fHONO_N(k)            = 0.0
  fNO_N(k)              = 0.0

endif ! NO / HONO

return
end subroutine libry_fNOHONO


! ACCUMULATE BSC PROPERTIES
!---------------------------------------------------------------

subroutine libry_accBSC (i,i2,j)
use libry_par
implicit none

integer :: i,i2,j

a_areaTHLC_g(i,i2,j)= a_areaTHLC_g(i,i2,j) +Agrid*o_LC(j)

a_areaTHDC_g(i,i2,j)= a_areaTHDC_g(i,i2,j) +Agrid*o_DC(j)

a_areaTHCC_g(i,i2,j)= a_areaTHCC_g(i,i2,j) +Agrid*o_CC(j)

a_areaTHMC_g(i,i2,j)= a_areaTHMC_g(i,i2,j) +Agrid*o_MC(j)

a_act(i,i2,j)       = a_act(i,i2,j) + real(ceiling(act(j))) * cweight

if (NOHONO) then
  a_fNO_N(i,i2,j)   = a_fNO_N(i,i2,j) + fNO_N(j) * Agrid

  a_fHONO_N(i,i2,j) = a_fHONO_N(i,i2,j) + fHONO_N(j) * Agrid
endif

a_rH2OlLC(i,i2,j)   = a_rH2OlLC(i,i2,j) + o_prs(j)*o_LC(j) &
                    * areaTH_s(i,i2,j) / max(p_critD,csumL)

a_rH2OlDC(i,i2,j)   = a_rH2OlDC(i,i2,j) + o_prs(j)*o_DC(j) &
                    * areaTH_s(i,i2,j) / max(p_critD,csumD)

a_rH2OlCC(i,i2,j)   = a_rH2OlCC(i,i2,j) + o_prs(j)*o_CC(j) &
                    * areaTH_s(i,i2,j) / max(p_critD,csumC)

a_rH2OlMC(i,i2,j)   = a_rH2OlMC(i,i2,j) + o_prs(j)*o_MC(j) &
                    * areaTH_s(i,i2,j) / max(p_critD,csumM)

a_TsLC(i,i2,j)      = a_TsLC(i,i2,j) +o_zt(j)*o_LC(j) &     ! xT_s
                    * areaTH_s(i,i2,j) / max(p_critD,csumL)
                                        
a_TsDC(i,i2,j)      = a_TsDC(i,i2,j) +o_zt(j)*o_DC(j) &
                    * areaTH_s(i,i2,j) / max(p_critD,csumD)
                                        
a_TsCC(i,i2,j)      = a_TsCC(i,i2,j) +o_zt(j)*o_CC(j) &
                    * areaTH_s(i,i2,j) / max(p_critD,csumC)
                                        
a_TsMC(i,i2,j)      = a_TsMC(i,i2,j) +o_zt(j)*o_MC(j) &
                    * areaTH_s(i,i2,j) / max(p_critD,csumM)

a_fCcbLC(i,i2,j)    = a_fCcbLC(i,i2,j) + fCcb(j)*o_LC(j) &
                    * c_MC * c_year * 1000.0 &
                    * areaTH_s(i,i2,j) / max(p_critD,csumL)

a_fCcbDC(i,i2,j)    = a_fCcbDC(i,i2,j) + fCcb(j)*o_DC(j) &
                    * c_MC * c_year * 1000.0 &
                    * areaTH_s(i,i2,j) / max(p_critD,csumD)

a_fCcbCC(i,i2,j)    = a_fCcbCC(i,i2,j) + fCcb(j)*o_CC(j) &
                    * c_MC * c_year * 1000.0 &
                    * areaTH_s(i,i2,j) / max(p_critD,csumC)

a_fCcbMC(i,i2,j)    = a_fCcbMC(i,i2,j) + fCcb(j)*o_MC(j) &
                    * c_MC * c_year * 1000.0 &
                    * areaTH_s(i,i2,j) / max(p_critD,csumM)

return
end subroutine libry_accBSC


! AVERAGE BSC PROPERTIES
!---------------------------------------------------------------

subroutine libry_avBSC (i)
use libry_par
implicit none

integer :: i

a_areaTHLC_g(i,:,:)             = a_areaTHLC_g(i,:,:) / real(naccu(i))
a_areaTHDC_g(i,:,:)             = a_areaTHDC_g(i,:,:) / real(naccu(i))
a_areaTHCC_g(i,:,:)             = a_areaTHCC_g(i,:,:) / real(naccu(i))
a_areaTHMC_g(i,:,:)             = a_areaTHMC_g(i,:,:) / real(naccu(i))

if (NOHONO) then
  a_fNO_N(i,:,:)                = a_fNO_N(i,:,:) / real(naccu(i))
  a_fHONO_N(i,:,:)              = a_fHONO_N(i,:,:) / real(naccu(i))
endif
a_rH2OlLC(i,:,:)                = a_rH2OlLC(i,:,:) / real(naccu(i))
a_rH2OlDC(i,:,:)                = a_rH2OlDC(i,:,:) / real(naccu(i))
a_rH2OlCC(i,:,:)                = a_rH2OlCC(i,:,:) / real(naccu(i))
a_rH2OlMC(i,:,:)                = a_rH2OlMC(i,:,:) / real(naccu(i))
a_TsLC(i,:,:)                   = a_TsLC(i,:,:) / real(naccu(i))
a_TsDC(i,:,:)                   = a_TsDC(i,:,:) / real(naccu(i))
a_TsCC(i,:,:)                   = a_TsCC(i,:,:) / real(naccu(i))
a_TsMC(i,:,:)                   = a_TsMC(i,:,:) / real(naccu(i))
a_fCcbLC(i,:,:)                 = a_fCcbLC(i,:,:) / real(naccu(i))
a_fCcbDC(i,:,:)                 = a_fCcbDC(i,:,:) / real(naccu(i))
a_fCcbCC(i,:,:)                 = a_fCcbCC(i,:,:) / real(naccu(i))
a_fCcbMC(i,:,:)                 = a_fCcbMC(i,:,:) / real(naccu(i))

return
end subroutine libry_avBSC


! AVERAGE GRID CELL BSC PROPERTIES
!---------------------------------------------------------------

subroutine libry_avgcBSC (i,i2,j)
use libry_par
implicit none

integer :: i,i2,j

a_MareaTHLC_g_S(i)              = a_MareaTHLC_g_S(i) + a_areaTHLC_g(i,i2,j)
a_MareaTHDC_g_S(i)              = a_MareaTHDC_g_S(i) + a_areaTHDC_g(i,i2,j)
a_MareaTHCC_g_S(i)              = a_MareaTHCC_g_S(i) + a_areaTHCC_g(i,i2,j)
a_MareaTHMC_g_S(i)              = a_MareaTHMC_g_S(i) + a_areaTHMC_g(i,i2,j)
if (NOHONO) then
  a_MfNO_N(i)                   = a_MfNO_N(i) + a_fNO_N(i,i2,j)
  a_MfHONO_N(i)                 = a_MfHONO_N(i) + a_fHONO_N(i,i2,j)
endif
a_MrH2OlLC_S(i)                 = a_MrH2OlLC_S(i) + a_rH2OlLC(i,i2,j)
a_MrH2OlDC_S(i)                 = a_MrH2OlDC_S(i) + a_rH2OlDC(i,i2,j)
a_MrH2OlCC_S(i)                 = a_MrH2OlCC_S(i) + a_rH2OlCC(i,i2,j)
a_MrH2OlMC_S(i)                 = a_MrH2OlMC_S(i) + a_rH2OlMC(i,i2,j)
a_MTsLC_S(i)                    = a_MTsLC_S(i) + a_TsLC(i,i2,j)
a_MTsDC_S(i)                    = a_MTsDC_S(i) + a_TsDC(i,i2,j)
a_MTsCC_S(i)                    = a_MTsCC_S(i) + a_TsCC(i,i2,j)
a_MTsMC_S(i)                    = a_MTsMC_S(i) + a_TsMC(i,i2,j)
a_MfCcbLC_S(i)                  = a_MfCcbLC_S(i) + a_fCcbLC(i,i2,j)
a_MfCcbDC_S(i)                  = a_MfCcbDC_S(i) + a_fCcbDC(i,i2,j)
a_MfCcbCC_S(i)                  = a_MfCcbCC_S(i) + a_fCcbCC(i,i2,j)
a_MfCcbMC_S(i)                  = a_MfCcbMC_S(i) + a_fCcbMC(i,i2,j)

return
end subroutine libry_avgcBSC

end module libry_opt
