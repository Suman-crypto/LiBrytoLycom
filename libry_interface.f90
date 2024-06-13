
module libry_interface
contains

!***********************************************************************
! LIBRY_INIT
!***********************************************************************

subroutine libry_init (nCPts2)
use libry_par
use libry_nova
use libry_land
use libry_opt
implicit none

integer :: i,i2,j,m
integer :: nCPts2

! set time step length

p_dt = real(tsl)

! initialise tiles!!!

call surf_init(nCPts2)

! initialise land surface

call land_init(nCPts2)

! initialise non-vascular vegetation

call nova_init(nCPts2)

return
end subroutine libry_init

!***********************************************************************
! LIBRY_STEP
!***********************************************************************

subroutine libry_step (nCPts2)
use libry_par
use libry_nova
use libry_land
use libry_opt
implicit none

integer :: i,t,v,h
integer :: nCPts2

do i = 1,nCPts2

  call land_stepG(i) ! same input for all tiles, levels and habitats

  do t = 1,p_ntiles

    if (frac_tile(i,t) .gt. p_critD) then

      do v = 1,p_nvert(t)

        call land_stepV(i,t,v) ! same input for all habitats

        do h = 1,p_nhab(t,v)

          ! land surface

          call land_step(i,t,v,h)

          ! non-vascular vegetation

          call nova_step(i,t,v,h)

          ! average properties per habitat

          call avg_land_nova(i,t,v,h)

        enddo

        ! average properties per level

        call avg_habitats(i,t,v)

      enddo

      ! average properties per tile

      call avg_levels(i,t)

    endif

  enddo

  ! average properties per tile

  call avg_tiles(i)

enddo

return
end subroutine libry_step

!***********************************************************************
! SURF_INIT
!***********************************************************************

subroutine surf_init (nCPts3)
use libry_par
use libry_opt
implicit none

integer :: i,t,v,h,j,m
integer :: nCPts3

!real :: MareaSTEM, MINareaSTEM, MAXareaSTEM
!real :: MareaLEAF, MINareaLEAF, MAXareaLEAF
real :: MSLAI, LAIcorr, fracForest

! general settings

p_nvert(1)                      = 2                                     ! 2 levels in forest
p_nvert(2)                      = 1                                     ! 1 level on bare/grass
p_nvert(3)                      = 1                                     ! 1 level in wetlands
p_nvert(4)                      = 1                                     ! 1 level on rocks

p_nhab(1,2)                     = 2                                     ! 2 habitats in forest canopy
p_nhab(1,1)                     = 1                                     ! 1 habitat on forest floor
p_nhab(2,:)                     = 1                                     ! 1 habitat on bare/grass
p_nhab(3,:)                     = 1                                     ! 1 habitat in wetlands
p_nhab(4,:)                     = 1                                     ! 1 habitat on rocks

! grid cell -specific surface properties
do i = 1,nCPts3

  ! cover fraction without any vegetation

  frac_noland(i)                = 1.0 - areaSOIL(i)

  if (frac_noland(i) .le. p_critD) frac_noland(i) = 0.0

  MareaSTEM(i)                  = sum(areaSTEM(i,:)) / 12.0
  MareaLEAF(i)                  = sum(areaLEAF(i,:)) / 12.0
  
  if (MareaSTEM(i) .gt. p_critD .and. MareaLEAF(i) .gt. p_critD) then
    fracSTEM(i)                 = MareaSTEM(i) / (MareaSTEM(i)+MareaLEAF(i))
    fracLEAF(i)                 = MareaLEAF(i) / (MareaSTEM(i)+MareaLEAF(i))
  else
    fracSTEM(i)                 = 0.5
    fracLEAF(i)                 = 0.5
  endif

  ! biome-specific surface properties per grid cell

  if (biome(i) .le. 6.0 .or. biome(i) .eq. 12.0) then                   ! forest biomes
  
    fracForest                  = min(1.0, (MareaLEAF(i)/p_LAImax)**p_LAIcorrF)

    LAIforest(i)                = fracForest &
                                * (p_LAImax - p_LAImin) + p_LAImin

    LAIgrass(i)                 = max(p_LAIGmin, &
                                  (MareaLEAF(i)-LAIforest(i)*fracForest)) &
                                / max(1.0 - fracForest, p_critD)
  else                                                                  ! no forest

    fracForest                  = 0.0

    LAIforest(i)                = 0.0

    LAIgrass(i)                 = MareaLEAF(i)
  endif                                                                  

! does not work as intended
!  LAIscale(i)                   = MareaLEAF(i) &
!                                / (p_LAIforest*fracForest +p_LAIgrass*(1.0-fracForest))


  ! tile-specific surface properties
  do t = 1,p_ntiles

    ! cover fraction of each tile

    if (t .eq. 1) then ! forest

      frac_tile(i,t)            = fracForest

      albLsf(i,t)               = p_albVEG

    elseif (t .eq. 2) then ! grassland/bare

      frac_tile(i,t)            = 1.0 - fracForest

      if (biome(i) .eq. 13.0) then
        albLsf(i,t)             = p_albDES
      else
        albLsf(i,t)             = p_albVEG
      endif
    else ! wetlands / rocks not yet considered

      frac_tile(i,t)            = 0.0

      albLsf(i,t)               = p_albVEG
    endif

    if (frac_tile(i,t) .le. p_critD) frac_tile(i,t) = 0.0

    ! Radiation

    MSLAI                       = LAIforest(i) + MareaSTEM(i)           ! averaged SLAI [m2 SL / m2 G]
                                
    !MSLAI                       = (MareaLEAF + MareaSTEM)

    if (t .eq. 1 .and. MSLAI .gt. p_critD &
          .and. frac_tile(i,t) .gt. p_critD) then ! forest

      Acano(i,t)                = ( 1.0 - exp(-p_beer_s * MSLAI)) &     ! vegetation area in canopy [m2]
                                / ( 1.0 *(1.0 - p_xcano) &
                                + exp(-p_beer_s * MSLAI) *p_xcano )

      Acano(i,t)                = min( Acano(i,t), MSLAI ) ! should not be necessary

    else ! no forest

      Acano(i,t)                = 0.0                                   ! no vegetation area in canopy
    endif

    ! level -specific surface properties
    do v = 1,p_nvert(t)

      ! ...

      ! habitat -specific surface properties
      do h = 1,p_nhab(t,v)

        if (biome(i) .eq. 1.0) then                                           ! tropical evergreen
       
          if (t .eq. 1) then ! forest

            if (v .eq. 1) then ! ground

              tauD(i,t,v,h)     = p_tauLEAF_trop !disturbance by leaves, turnover time in months
            else
              if (h .eq. 1) then ! stems

                tauD(i,t,v,h)   = p_tauG_tropwet*12.0

              else ! leaves
                tauD(i,t,v,h)   = p_tauLEAF_trop
              endif
            endif
          else ! no forest
            tauD(i,t,v,h)       = p_tauG_tropwet*12.0
          endif
       
        elseif (biome(i) .eq. 2.0) then                                       ! tropical deciduous

          tauD(i,t,v,h)         = p_tauG_tropdry*12.0

        elseif (biome(i) .eq. 3.0) then                                       ! tropical coniferous

          tauD(i,t,v,h)         = p_tauG_tropwet*12.0

        elseif (biome(i) .eq. 4.0) then                                       ! temperate deciduous

          tauD(i,t,v,h)         = p_tauG_temp*12.0

        elseif (biome(i) .eq. 5.0) then                                       ! temperate evergreen

          if (t .eq. 1) then ! forest

            if (v .eq. 1) then ! ground

              tauD(i,t,v,h)     = p_tauLEAF_trop !disturbance by leaves
            else
              if (h .eq. 1) then ! stems

                tauD(i,t,v,h)   = p_tauG_temp*12.0

              else ! leaves
                tauD(i,t,v,h)   = p_tauLEAF_trop
              endif
            endif
          else ! no forest
            tauD(i,t,v,h)       = p_tauG_temp*12.0
          endif

        elseif (biome(i) .eq. 6.0) then                                       ! boreal evergreen

          tauD(i,t,v,h)         = p_tauG_bor*12.0

        elseif (biome(i) .eq. 7.0) then                                       ! savanna

          tauD(i,t,v,h)         = p_tauG_sav*12.0

        elseif ((biome(i) .eq. 8.0).or.(biome(i) .eq. 10.0)) then             ! grassland

          tauD(i,t,v,h)         = p_tauG_gra*12.0

        elseif (biome(i) .eq. 11.0) then                                      ! tundra

          tauD(i,t,v,h)         = p_tauG_tun*12.0

        elseif (biome(i) .eq. 13.0) then                                      ! desert

          tauD(i,t,v,h)         = p_tauG_des*12.0

        elseif (biome(i) .eq. 12.0) then                                      ! mediterranean vegetation

          if (t .eq. 1) then ! forest

            if (v .eq. 1) then ! ground

              tauD(i,t,v,h)     = p_tauG_med*12.0
            else
              if (h .eq. 1) then ! stems

                tauD(i,t,v,h)   = p_tauG_med*12.0

              else ! leaves
                tauD(i,t,v,h)   = p_tauLEAF_med
              endif
            endif
          else ! no forest
            tauD(i,t,v,h)       = p_tauG_med*12.0
          endif

        elseif (biome(i) .eq. 20.0) then                                      ! ordovician landscape

          tauD(i,t,v,h)         = p_tauG_pal*12.0

        else                                                                  ! no vegetation (e.g. flooded or ice)

          tauD(i,t,v,h)         = 1.0                                           

        endif                                                                  
   
        ! boundary layer conductance

        if (t .eq. 1) then ! forest

          if (v .eq. 1) then ! ground

            roughlen(i,t,v,h)   = p_roughlen_floor
          else
            if (h .eq. 1) then ! stems

              roughlen(i,t,v,h) = p_roughlen_stems
            else ! leaves

              roughlen(i,t,v,h) = p_roughlen_leaves
            endif
          endif
        else ! no forest
          roughlen(i,t,v,h)     = p_roughlen_grd
        endif

!        dheight                 = p_ratio_dr * roughlen
!        roughlen_h              = p_ratio_mh * roughlen
!   
!        dn_kH2Og(i,t,v,h)       = log((p_mheight-dheight)/roughlen) &          
!                                * log((p_mheight-dheight)/roughlen_h)          

        ! water reservoir below thallus
   
        rmaxH2Ol_b(i,t,v,h)     = 0.0

        if (t .eq. 1) then ! forest
          if (v .eq. 2) then ! canopy
            if (h .eq. 1) then ! stems

              rmaxH2Ol_b(i,t,v,h) = p_rmaxH2Ol_bs
            else
              rmaxH2Ol_b(i,t,v,h) = 0.0
            endif
          else
            rmaxH2Ol_b(i,t,v,h) = p_rmaxH2Ol_bg
          endif
        else
          rmaxH2Ol_b(i,t,v,h) = p_rmaxH2Ol_bg
        endif

      enddo
    enddo

  enddo
enddo

return
end subroutine surf_init

!***********************************************************************
! AVG_LAND_NOVA
!***********************************************************************

subroutine avg_land_nova (i,t,v,h)
use libry_par
use libry_opt
implicit none

integer :: i,t,v,h

! properties which are needed as input for other habitats / levels / tiles
! -> needed every time step !!

ah_fH2Ol_xd(h)                  = fH2Ol_xd0 *1000.0*c_year * (1.0-as_areaTH_s) &
                                + as_fH2Ol_td

ah_fH2Ol_bd(h)                  = fH2Ol_bd0 *1000.0*c_year * (1.0-as_areaTH_s) &
                                + as_fH2Ol_bd

! -> needed only for output !!
if (writeout) then

  ! properties which are specific for lichens and mosses
  ah_areaTH_s(h)                = as_areaTH_s

  ah_rH2Ol(h)                   = as_rH2Ol_t
                                  !rH2Ol_0(i,t,v,h) /rmaxH2Ol_0 * (1.0-as_areaTH_s) &
  
  ah_rH2Os(h)                   = as_rH2Os_t

  ah_rmaxH2Ol(h)                = as_rmaxH2Ol_t
                                  !rmaxH2Ol_0 *1000.0 * (1.0-as_areaTH_s)

  ah_act(h)                     = as_act
  ah_actB(h)                    = as_actB
  
  ah_rCO2d(h)                   = as_rCO2d
  
  ah_rCb(h)                     = as_rCb
  
  ah_fCO2gc(h)                  = as_fCO2gc
  
  ah_fCcg(h)                    = as_fCcg
  
  ah_fCcb(h)                    = as_fCcb
  ah_fCcb_l(h)                  = as_fCcb_l
  ah_fCcb_c(h)                  = as_fCcb_c
  
  ah_fCbo(h)                    = as_fCbo
  ah_frgr(h)                    = as_frgr

  ! properties which are averaged over the surface and the lichen and moss cover

  ah_fH2Ogl_ux(h)               = min(0.0,fH2Olg_xu0) *1000.0*c_year * (1.0-as_areaTH_s) &
                                + as_fH2Ogl_ut
  
  ah_fH2Olg_xu(h)               = max(0.0,fH2Olg_xu0) *1000.0*c_year * (1.0-as_areaTH_s) &
                                + as_fH2Olg_tu
  
  ah_fH2Ol_bx(h)                = as_fH2Ol_bx

    ah_Ts(h)                    = a0_Ts * (1.0-as_areaTH_s) &
                                + as_Ts * as_areaTH_s
                        
  ah_H(h)                       = fRAD_H0 * (1.0-as_areaTH_s) &
                                + as_H
                        
  ah_E(h)                       = fQ_ta_L0 * (1.0-as_areaTH_s) &
                                + as_E
                        
  ah_C(h)                       = fQ_ta_S0 * (1.0-as_areaTH_s) &
                                + as_C
                        
  ah_EB(h)                      = (fRAD_H0-fQ_ta_L0-fQ_ta_S0) * (1.0-as_areaTH_s) &
                                + as_EB

  if (v .eq. 1) then ! ground

    ah_Tg(h)                    = xT_g0(i,t,h) * (1.0-as_areaTH_s) &
                                + as_Tg * as_areaTH_s
                      
    ah_G(h)                     = fQ_tg0 * (1.0-as_areaTH_s) &
                                + as_G
  endif
endif
                      
return
end subroutine avg_land_nova

!***********************************************************************
! AVG_HABITATS
!***********************************************************************

subroutine avg_habitats (i,t,v)
use libry_par
use libry_opt
implicit none

integer :: i,t,v
integer :: k

real    :: wgthabC,wgthabB,wgthab

! initialise accumulated variables *per level* with zero
av_fH2Ol_xd(v)                  = 0.0
av_fH2Ol_bd(v)                  = 0.0

! loop over all habitats per level
do k=1,p_nhab(t,v)

  ! determine weight of each habitat

  if (v .eq. 2) then ! canopy
    if (k .eq. 1) then ! stems
      wgthabC                   = fracSTEM(i)
      wgthabB                   = 1.0 ! all stemflow from stems, flux per BARK area
    else
      wgthabC                   = fracLEAF(i)
      wgthabB                   = 0.0 ! zero stemflow from leaves
    endif
  else ! forest floor or other tiles
    wgthabC                     = 1.0 ! only one habitat
    wgthabB                     = 1.0
  endif

  ! average habitats per level (v)
  av_fH2Ol_xd(v)                = av_fH2Ol_xd(v) + ah_fH2Ol_xd(k) * wgthabC

  av_fH2Ol_bd(v)                = av_fH2Ol_bd(v) + ah_fH2Ol_bd(k) * wgthabB
enddo

! only if output is written
if (writeout) then

  ! same flux for all habitats
  av_fH2Ol_ux(v)                = fH2Ol_ux * 1000.0 *c_year

  ! flux/reservoir depends on habitat
  av_areaTH_s(v)                = 0.0

  av_rH2Ol(v)                   = 0.0
  av_rH2Os(v)                   = 0.0
  av_rmaxH2Ol(v)                = 0.0
  av_act(v)                     = 0.0
  av_actB(v)                    = 0.0
  av_rCO2d(v)                   = 0.0
  av_rCb(v)                     = 0.0
  av_fCO2gc(v)                  = 0.0
  av_fCcg(v)                    = 0.0
  av_fCcb(v)                    = 0.0
  av_fCcb_l(v)                  = 0.0
  av_fCcb_c(v)                  = 0.0
  av_fCbo(v)                    = 0.0
  av_frgr(v)                    = 0.0
  
  av_fH2Ogl_ux(v)               = 0.0
  av_fH2Olg_xu(v)               = 0.0
  av_fH2Ol_bx(v)                = 0.0
  av_Ts(v)                      = 0.0
  av_H(v)                       = 0.0
  av_E(v)                       = 0.0
  av_C(v)                       = 0.0
  av_EB(v)                      = 0.0
  
  if (v .eq. 1) then ! ground
    av_Tg                       = 0.0
    av_G                        = 0.0
  endif

  ! loop over all habitats per level
  do k=1,p_nhab(t,v)

    ! determine weight of each habitat
   
    if (v .eq. 2) then ! canopy
      if (k .eq. 1) then ! stems
        wgthab                  = fracSTEM(i)
      else
        wgthab                  = fracLEAF(i)
      endif
    else ! forest floor or other tiles
      wgthab                    = 1.0 ! only one habitat
    endif
   
    ! average habitats per level (v)
    av_areaTH_s(v)              = av_areaTH_s(v) + ah_areaTH_s(k) * wgthab

    av_rH2Ol(v)                 = av_rH2Ol(v) + ah_rH2Ol(k) * wgthab
    av_rH2Os(v)                 = av_rH2Os(v) + ah_rH2Os(k) * wgthab
    av_rmaxH2Ol(v)              = av_rmaxH2Ol(v) + ah_rmaxH2Ol(k) * wgthab
    av_act(v)                   = av_act(v) + ah_act(k) * wgthab
    av_actB(v)                  = av_actB(v) + ah_actB(k) * wgthab
    av_rCO2d(v)                 = av_rCO2d(v) + ah_rCO2d(k) * wgthab
    av_rCb(v)                   = av_rCb(v) + ah_rCb(k) * wgthab
    av_fCO2gc(v)                = av_fCO2gc(v) + ah_fCO2gc(k) * wgthab
    av_fCcg(v)                  = av_fCcg(v) + ah_fCcg(k) * wgthab
    av_fCcb(v)                  = av_fCcb(v) + ah_fCcb(k) * wgthab
    av_fCcb_l(v)                = av_fCcb_l(v) + ah_fCcb_l(k) * wgthab
    av_fCcb_c(v)                = av_fCcb_c(v) + ah_fCcb_c(k) * wgthab
    av_fCbo(v)                  = av_fCbo(v) + ah_fCbo(k) * wgthab
    av_frgr(v)                  = av_frgr(v) + ah_frgr(k) * wgthab
   
    av_fH2Ogl_ux(v)             = av_fH2Ogl_ux(v) + ah_fH2Ogl_ux(k) * wgthab
    av_fH2Olg_xu(v)             = av_fH2Olg_xu(v) + ah_fH2Olg_xu(k) * wgthab

    av_fH2Ol_bx(v)              = av_fH2Ol_bx(v) + ah_fH2Ol_bx(k) * wgthab

    av_Ts(v)                    = av_Ts(v) + ah_Ts(k) * wgthab
    av_H(v)                     = av_H(v) + ah_H(k) * wgthab
    av_E(v)                     = av_E(v) + ah_E(k) * wgthab
    av_C(v)                     = av_C(v) + ah_C(k) * wgthab
    av_EB(v)                    = av_EB(v) + ah_EB(k) * wgthab
   
    if (v .eq. 1) then ! ground
   
      av_Tg                     = av_Tg + ah_Tg(k) * wgthab
      av_G                      = av_G + ah_G(k) * wgthab
    endif
  enddo

  if (v .eq. 1) then ! ground
    av_rH2Ol_g1                 = rH2Ol_g1(i,t) /p_rmaxH2Ol_g1 !*1000.0           !av_rH2Ol_g1 + 
    av_rH2Ol_g2                 = rH2Ol_g2(i,t) /p_rmaxH2Ol_g2 !*1000.0           !av_rH2Ol_g2 + 

    av_fH2Ol_ug                 = fH2Ol_ug *1000.0*c_year       !av_fH2Ol_ug + 
    av_fH2Ol_go                 = fH2Ol_go *1000.0*c_year       !av_fH2Ol_go + 
    av_fH2Ol_gb                 = fH2Ol_gb *1000.0*c_year       !av_fH2Ol_gb + 
    av_fH2Olg_ga                = fH2Olg_ga *1000.0*c_year     !av_fH2Olg_ga +

    av_fH2Ol_ug2                = 0.0

  else ! canopy

    av_fH2Ol_ug2                = fH2Ol_ug2(i,t) *1000.0*c_year

  endif
endif

return
end subroutine avg_habitats

!***********************************************************************
! AVG_LEVELS
!***********************************************************************

subroutine avg_levels (i,t)
use libry_par
use libry_opt
implicit none

integer :: i,t,k

real    :: ratioCG, ratioBG

do k=1,p_nvert(t)

  ! convert flux per canopy area back to grid cell area
  ratioCG = 1.0
  ratioBG = 1.0
  if (k .eq. 2) then
    ratioCG = max(1.0,Acano(i,t))
    ratioBG = max(1.0,Acano(i,t)*fracSTEM(i))
  endif

  at_fH2Ol_xd(t,k)                 = av_fH2Ol_xd(k)*ratioCG

  at_fH2Ol_bd(t,k)                 = av_fH2Ol_bd(k)*ratioBG
enddo

! only if output is written
if (writeout) then

  do k=1,p_nvert(t)

    ratioCG = 1.0
    if (k .eq. 2) ratioCG = max(1.0,Acano(i,t))

    ! average levels per tile (t) - separate output (C,G) for several properties
    at_fH2Ol_ux(t,k)               = av_fH2Ol_ux(k) *ratioCG

    at_areaTH_s(t,k)               = av_areaTH_s(k) ! this always refers to cover per available area

    at_rH2Ol(t,k)                  = av_rH2Ol(k) *ratioCG
    at_rH2Os(t,k)                  = av_rH2Os(k) *ratioCG
    at_rmaxH2Ol(t,k)               = av_rmaxH2Ol(k) *ratioCG
    at_act(t,k)                    = av_act(k)
    at_actB(t,k)                   = av_actB(k)
    at_rCO2d(t,k)                  = av_rCO2d(k)
    at_rCb(t,k)                    = av_rCb(k) *ratioCG
    at_fCO2gc(t,k)                 = av_fCO2gc(k) *ratioCG
    at_fCcg(t,k)                   = av_fCcg(k) *ratioCG
    at_fCcb(t,k)                   = av_fCcb(k) *ratioCG
    at_fCcb_l(t,k)                 = av_fCcb_l(k) *ratioCG
    at_fCcb_c(t,k)                 = av_fCcb_c(k) *ratioCG
    at_fCbo(t,k)                   = av_fCbo(k) *ratioCG
    at_frgr(t,k)                   = av_frgr(k) *ratioCG

    at_fH2Ogl_ux(t,k)              = av_fH2Ogl_ux(k) *ratioCG
    at_fH2Olg_xu(t,k)              = av_fH2Olg_xu(k) *ratioCG
    at_fH2Ol_bx(t,k)               = av_fH2Ol_bx(k) *ratioCG
    at_Ts(t,k)                     = av_Ts(k)
    at_H(t,k)                      = av_H(k)
    at_E(t,k)                      = av_E(k)
    at_C(t,k)                      = av_C(k)
    at_EB(t,k)                     = av_EB(k)
  enddo

  at_Tg(t)                      = av_Tg
  at_G(t)                       = av_G

  at_rH2Ol_g1(t)                = av_rH2Ol_g1
  at_rH2Ol_g2(t)                = av_rH2Ol_g2
                        
  at_fH2Ol_ug(t)                = av_fH2Ol_ug
  at_fH2Ol_go(t)                = av_fH2Ol_go
  at_fH2Ol_gb(t)                = av_fH2Ol_gb
  at_fH2Olg_ga(t)               = av_fH2Olg_ga

  at_fH2Ol_ug2(t)               = av_fH2Ol_ug2
endif

return
end subroutine avg_levels

!***********************************************************************
! AVG_TILES
!***********************************************************************

subroutine avg_tiles (i)
use libry_par
use libry_opt
implicit none

integer :: i
integer :: k,l
real :: ft2

! save averaged flux-state variables for next time step
do k=1,p_ntiles

  atN_fH2Ol_xd(i,k,:)           = at_fH2Ol_xd(k,:)

  atN_fH2Ol_bd(i,k,:)           = at_fH2Ol_bd(k,:)
enddo

! only if output is written
if (writeout) then

! accumulate at the grid cell level -> set only to zero when output is written and reset is called!!

  ! loop over all tiles per grid cell
  do k=1,p_ntiles

    if (frac_tile(i,k) .gt. p_critD) then

      do l=1,p_nvert(k)

        if (l .eq. 1) then
          ft2 = frac_tile(i,k)
        else
          ft2 = 1.0 !only works if there is max one tile with a canopy
        endif

        ! average tiles per grid cell (i) - separate output (C,G) for several properties
        ag_fH2Ol_ux(i,l)        = ag_fH2Ol_ux(i,l)  + at_fH2Ol_ux(k,l) * frac_tile(i,k)
        ag_fH2Ol_xd(i,l)        = ag_fH2Ol_xd(i,l)  + at_fH2Ol_xd(k,l) * frac_tile(i,k)

        ag_areaTH_s(i,l)        = ag_areaTH_s(i,l)  + at_areaTH_s(k,l) * frac_tile(i,k)

        ag_rH2Ol(i,l)           = ag_rH2Ol(i,l)     + at_rH2Ol(k,l) * frac_tile(i,k)
        ag_rH2Os(i,l)           = ag_rH2Os(i,l)     + at_rH2Os(k,l) * frac_tile(i,k)
        ag_rmaxH2Ol(i,l)        = ag_rmaxH2Ol(i,l)  + at_rmaxH2Ol(k,l) * frac_tile(i,k)
        ag_act(i,l)             = ag_act(i,l)       + at_act(k,l) * ft2
        ag_actB(i,l)            = ag_actB(i,l)      + at_actB(k,l) * ft2
        ag_rCO2d(i,l)           = ag_rCO2d(i,l)     + at_rCO2d(k,l) * ft2
        ag_rCb(i,l)             = ag_rCb(i,l)       + at_rCb(k,l) * frac_tile(i,k)
        ag_fCO2gc(i,l)          = ag_fCO2gc(i,l)    + at_fCO2gc(k,l) * frac_tile(i,k)
        ag_fCcg(i,l)            = ag_fCcg(i,l)      + at_fCcg(k,l) * frac_tile(i,k)
        ag_fCcb(i,l)            = ag_fCcb(i,l)      + at_fCcb(k,l) * frac_tile(i,k)
        ag_fCcb_l(i,l)          = ag_fCcb_l(i,l)    + at_fCcb_l(k,l) * frac_tile(i,k)
        ag_fCcb_c(i,l)          = ag_fCcb_c(i,l)    + at_fCcb_c(k,l) * frac_tile(i,k)
        ag_fCbo(i,l)            = ag_fCbo(i,l)      + at_fCbo(k,l) * frac_tile(i,k)
        ag_frgr(i,l)            = ag_frgr(i,l)      + at_frgr(k,l) * ft2

        ag_fH2Ogl_ux(i,l)       = ag_fH2Ogl_ux(i,l) + at_fH2Ogl_ux(k,l) * frac_tile(i,k)
        ag_fH2Olg_xu(i,l)       = ag_fH2Olg_xu(i,l) + at_fH2Olg_xu(k,l) * frac_tile(i,k)
        ag_fH2Ol_bx(i,l)        = ag_fH2Ol_bx(i,l)  + at_fH2Ol_bx(k,l) * frac_tile(i,k)
        ag_fH2Ol_bd(i,l)        = ag_fH2Ol_bd(i,l)  + at_fH2Ol_bd(k,l) * frac_tile(i,k)

        ag_Ts(i,l)              = ag_Ts(i,l)        + at_Ts(k,l) * ft2
        ag_H(i,l)               = ag_H(i,l)         + at_H(k,l) * ft2
        ag_E(i,l)               = ag_E(i,l)         + at_E(k,l) * ft2
        ag_C(i,l)               = ag_C(i,l)         + at_C(k,l) * ft2
        ag_EB(i,l)              = ag_EB(i,l)        + at_EB(k,l) * ft2
      enddo
     
      ag_Tg(i)                  = ag_Tg(i)         + at_Tg(k) * frac_tile(i,k)
      ag_G(i)                   = ag_G(i)          + at_G(k) * frac_tile(i,k)

      ag_fH2Ol_ug2(i)           = ag_fH2Ol_ug2(i)  + at_fH2Ol_ug2(k) * frac_tile(i,k)
      ag_rH2Ol_g1(i)            = ag_rH2Ol_g1(i)   + at_rH2Ol_g1(k) * frac_tile(i,k)
      ag_rH2Ol_g2(i)            = ag_rH2Ol_g2(i)   + at_rH2Ol_g2(k) * frac_tile(i,k)
     
      ag_fH2Ol_ug(i)            = ag_fH2Ol_ug(i)   + at_fH2Ol_ug(k) * frac_tile(i,k)
      ag_fH2Ol_go(i)            = ag_fH2Ol_go(i)   + at_fH2Ol_go(k) * frac_tile(i,k)
      ag_fH2Ol_gb(i)            = ag_fH2Ol_gb(i)   + at_fH2Ol_gb(k) * frac_tile(i,k)
      ag_fH2Olg_ga(i)           = ag_fH2Olg_ga(i)  + at_fH2Olg_ga(k) * frac_tile(i,k)
    endif
  enddo
  
  ag_rH2Os_g(i)                 = ag_rH2Os_g(i)    + dsnow *1000.0
  ag_fH2Osl_g(i)                = ag_fH2Osl_g(i)   + fH2Osl_g *1000.0*c_year
  
  ag_fH2Os_ad(i)                = ag_fH2Os_ad(i)   + fH2Os_ad(i) *1000.0*c_year
  ag_xT_a(i)                    = ag_xT_a(i)       + xT_a(i)
  ag_fH2Ol_ad(i)                = ag_fH2Ol_ad(i)   + fH2Ol_ad(i) *1000.0*c_year
  ag_fRADs(i)                   = ag_fRADs(i)      + fRADs_ad(i)
endif

return
end subroutine avg_tiles

end module libry_interface

