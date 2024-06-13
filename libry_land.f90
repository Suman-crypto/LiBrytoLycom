
module libry_land
contains

!***********************************************************************
! LAND_INIT
!***********************************************************************

subroutine land_init (nCPts3)
use libry_par
use libry_opt
implicit none

integer :: i,i2,j,m
integer :: nCPts3

! initialise atmospheric conditions

rCO2g_a                         = 360.0                                 ! [ppm]
rO2g_a                          = 210000.0                              ! [ppm]

! initialise fields

rH2Os_g(:)                      = 0.0                                   ! snow layer [m3 H2O / m2 G]
                               
rH2Ol_0(:,:,:,:)                = 0.0                                   ! water on leaves [m3 H2O / m2 C]

rH2Ol_g1(:,:)                   = 0.0                                   ! water in the top soil [m3 H2O / m2 C]
                               
rH2Ol_g2(:,:)                   = 0.0                                   ! water in the bulk soil [m3 H2O / m2 C]

atN_fH2Ol_xd(:,:,:)             = 0.0                                   ! water overflow from layers [m/s]

atN_fH2Ol_bd(:,:,:)             = 0.0                                   ! overflow from bark/substrate [m/s]

fH2Ol_ug2(:,:)                  = 0.0                                   ! throughfall [m/s]

! initialise Soil & air properties

do i=1,nCPts3

  ! initial ground temperature
  xT_g(i,:,:,:)                 = 298.0 !xT_a(i)                               ! ground temperature [K]
  xT_g0(i,:,:)                  = 298.0 !xT_a(i)                               ! bare ground temperature [K]

  ! heat conduction
  if ( (biome(i) .le. 6.0) .or. (biome(i) .eq. 12.0) ) then

    CSOIL(i)                    = p_CSOIL_for
    kSOIL(i)                    = p_kSOIL_for
  
    dewmax(i)                   = -p_dewmax_for &                       ! transform value mm/yr to m/s, assuming 12h night
                                / 1000.0 / 365.25 / 12.0 / 3600.0
  else
    CSOIL(i)                    = p_CSOIL_des
    kSOIL(i)                    = p_kSOIL_des
  
    dewmax(i)                   = -p_dewmax_des & 
                                / 1000.0 / 365.25 / 12.0 / 3600.0
  endif
enddo

return
end subroutine land_init

!***********************************************************************
! LAND_STEP G
!***********************************************************************

subroutine land_stepG (i)
use libry_par
use libry_opt
implicit none

integer :: i

! day and time step counter
if (fRADs_ad(i) .gt. 0.0) then
  kday                          = 1
  naccu_d(i)                    = naccu_d(i) + 1
else
  kday                          = 0
  naccu_n(i)                    = naccu_n(i) + 1
endif

naccu(i)                        = naccu(i) + 1

! saturation vapor pressure
Ta3                             = xT_a(i) * xT_a(i) * xT_a(i)

Ta4                             = Ta3 * xT_a(i)

zT_a                            = xT_a(i) - c_TH2Osl                    ! air temperature [deg C]
                 
esatAIR                         = p_esatAIR3 *exp(p_esatAIR1*zT_a &     ! saturation vapour pressure [kg * m / (s2 * m2) = Pa]                  REF: Allen,1998
                                / (p_esatAIR2 + zT_a))
                 
desatdT                         = exp(p_esatAIR1*zT_a/(p_esatAIR2 &     ! slope of saturation vapour pressure vs temperature relationship []
                                + zT_a)) &
                                * (p_esatAIR1 *p_esatAIR2 *p_esatAIR3 &
                                / ((p_esatAIR2 + zT_a) &
                                * (p_esatAIR2 + zT_a)))

! Snow balance
fH2Osl_g                        = min(3.22 * max(0.0,xT_a(i)-c_TH2Osl)& ! snow melt [m/s]                                                       REF: Bergstrom,1992
                                / c_day / 1000.0, &
                                  rH2Os_g(i) / p_dt + fH2Os_ad(i))

rH2Os_g(i)                      = max(0.0, rH2Os_g(i) &
                                + fH2Os_ad(i) * p_dt &
                                - fH2Osl_g * p_dt &
                                - rH2Os_g(i) * p_H2Os_loss*p_dt)

dsnow                           = rH2Os_g(i) * c_rhoH2Ol / p_rhoH2Os    ! snow depth [m]

return
end subroutine land_stepG

!***********************************************************************
! LAND_STEP V
!***********************************************************************

subroutine land_stepV (i,t,v)
use libry_par
use libry_opt
implicit none

integer :: i,t,v

real :: fracRADs_0

! Switches

if (v .eq. 2) then ! canopy
  lground                       = 0.0
else
  lground                       = 1.0
endif

! Radiation

Aslai                           = LAIforest(i) &
                                * areaLEAF_month(i) &
                                / max(p_critD,MareaLEAF(i)) &
                                + areaSTEM_month(i)                     ! SLAI [m2 SL / m2 G]

AlaiG                           = LAIgrass(i) &
                                * areaLEAF_month(i) &
                                / max(p_critD,MareaLEAF(i))             ! LAI grassland [m2 SL / m2 G]

if (Aslai .le. p_critD) Aslai = 0.0
if (AlaiG .le. p_critD) AlaiG = 0.0

if (t .eq. 1) then ! forest

  if (v .eq. 2) then ! canopy
    fracRADs_0                  = 1.0 - exp(-p_beer_s * Aslai)          ! fraction of total absorbed short-wave radiation in canopy             REF: Bonan,2008

    fracRADs                    = fracRADs_0 / max(1.0,Acano(i,t))      ! average intensity of absorbed short-wave radiation in canopy

    fracRADl                    = 1.0 - exp(-p_beer_l * Aslai)          ! fraction of total absorbed long-wave radiation in canopy              REF: Bonan,2008
  else ! ground

    fracRADs                    = exp(-p_beer_s * Aslai)                ! fraction of total absorbed short-wave radiation on ground             REF: Bonan,2008

    fracRADl                    = exp(-p_beer_l * Aslai)                ! fraction of total absorbed long-wave radiation on ground              REF: Bonan,2008
  endif
else ! no forest

  fracRADs                      = exp(-p_beer_s * AlaiG)                ! fraction of total absorbed short-wave radiation on ground             REF: Bonan,2008

  fracRADl                      = exp(-p_beer_l * AlaiG)                ! fraction of total absorbed long-wave radiation on ground              REF: Bonan,2008
endif

! Water fluxes

if (v .eq. 2) then ! canopy -> only called for forest tile! 

  fracrain                      = min(1.0, (areaLEAF_month(i) &         ! [ ]
                                / p_LAImax )**p_fracIcpt)

  fH2Ol_ux                      = fH2Ol_ad(i)*fracrain &                ! water input into canopy [m3 H2O / (m2 C * s)]
                                / max(1.0,Acano(i,t))

  fH2Ol_ug2(i,t)                = fH2Ol_ad(i)*(1.0 - fracrain) &        ! throughfall (direct+drip) [ m3 H2O / m2 G / s ]
                                + atN_fH2Ol_xd(i,t,2) /1000.0 /c_year &
                                * (1.0 - p_fracBark)

  fH2Ol_ub                      = atN_fH2Ol_xd(i,t,2) /1000.0 /c_year & ! water input into bark [m3 H2O / (m2 Bark * s)]
                                * p_fracBark /max(1.0, Acano(i,t)*fracSTEM(i))
else

!  fracrain                      = 1.0                                   ! [ ] not used yet

  if (t .eq. 1) then ! forest

    fH2Ol_ux                    = atN_fH2Ol_bd(i,t,2) /1000.0 /c_year & ! water input onto ground (throughfall+stemflow) [m3 H2O / (m2 G * s)]
                                + fH2Ol_ug2(i,t)

    fH2Ol_ug2(i,t)              = 0.0                                   ! no throughfall from ground layer
  else

    fH2Ol_ux                    = fH2Ol_ad(i)                           ! water input onto ground [m3 H2O / (m2 C * s)]

    fH2Ol_ug2(i,t)              = 0.0                                   ! no throughfall without canopy
  endif

  fH2Ol_ux                      = fH2Ol_ux + fH2Osl_g                   ! add snow melt

  fH2Ol_ub                      = atN_fH2Ol_xd(i,t,1) /1000.0 /c_year   ! water input into substrate [m3 H2O / (m2 G * s)]

  fH2Ol_ug                      = atN_fH2Ol_bd(i,t,1) /1000.0 /c_year   ! water input into soil [m3 H2O / (m2 G * s)]
endif

return
end subroutine land_stepV

!***********************************************************************
! LAND_STEP
!***********************************************************************

subroutine land_step (i,t,v,h)
use libry_par
use libry_opt
implicit none

integer :: i,t,v,h,j,m

real    :: rmaxH2Ol_0, fracWup_0
real    :: kH2Og0, ETpot0
real    :: grh0, cdg0, xT_s_wet0, xT_s_dry0, crh0
real    :: wetfrac_0
real    :: dRAD, fRAD_Hw0, fRAD_Hd0
real    :: Rnet, ETpot_v
real    :: percolation, soilevap, rootuptk, transpiration

! Interception by leaves -----------------------------------------------

if (t .le. 2) then ! only forests and grassland

  if (v .eq. 2) then ! canopy

    if (h .eq. 1 ) then

      rmaxH2Ol_0                = p_rmaxH2Ol_cs0                        ! [m]                                                                   REF: Miralles,2010

      fracWup_0                 = p_fracWup0S                           ! [ ]
    else
      rmaxH2Ol_0                = p_rmaxH2Ol_cl0                        ! [m]                                                                   REF: Miralles,2010

      fracWup_0                 = p_fracWup0L                           ! [ ]
    endif

  else ! ground

    rmaxH2Ol_0                  = p_rmaxH2Ol_s0                         ! [m]                                                                   REF: Miralles,2010

    fracWup_0                   = p_fracWup0L                           ! [ ]
  endif
                          
  kH2Og0                        = fAIR_s(i) / (75.0 &                   ! [m / s]
                                * 8.0**roughlen(i,t,v,h))

!  kH2Og0                        = p_vonKarman * p_vonKarman &           ! [m / s]                                                               REF: Allen,1998
!                                * max(p_critD,fAIR_s(i)) &
!                                / dn_kH2Og(i,t,v,h)
                          
  if (v .eq. 1 .and. dsnow .ge. p_H2Os_crit) then ! snow layer at ground too thick

    ETpot0                      = 0.0
    fH2Olg_xu0                  = 0.0
    fH2Ol_xd0                   = 0.0
    fH2Ol_bd0                   = 0.0
    xT_s_dry0                   = min(c_TH2Osl - 0.1, xT_a(i))
    xT_s0                       = xT_s_dry0
    fRAD_H0                     = 0.0
    fQ_ta_L0                    = 0.0
    fQ_ta_S0                    = 0.0
  else

    grh0                        = c_gamma / kH2Og0
    cdg0                        = c_CAIR*(desatdT+c_gamma)
                               
    xT_s_wet0                   = ( grh0 * (fRADs_ad(i)*fracRADs &      ! surface temperature [K]                                               REF: Monteith,1981
                                * (1.0-albLsf(i,t)) &
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) &
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL*xT_g0(i,t,h) * lground ) &
                                + xT_a(i)*cdg0 &
                                - c_CAIR*(esatAIR-rH2Og_RH(i)*esatAIR) ) &
                                / ( grh0*(4.0*p_eps*c_sigma*Ta3 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL * lground) + cdg0 )

    crh0                        = 1.0/(c_CAIR*kH2Og0)
                               
    xT_s_dry0                   = ( xT_a(i) &
                                + crh0*(fracRADs*(1.0-albLsf(i,t))*fRADs_ad(i) &
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) &
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL*xT_g0(i,t,h) * lground) ) &
                                / ( 1.0 + crh0*(4.0*p_eps*c_sigma*Ta3 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL * lground) )

    dRAD                        = fRADs_ad(i)*fracRADs*(1.0-albLsf(i,t))&  ! downwelling radiation [W / m2 T]
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) & ! lrad from below
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground)

    fRAD_Hd0                    = dRAD &                                ! net radiation (dry) [W / m2 T]
                                - 4.0*p_eps*c_sigma*Ta3*xT_s_dry0 * (2.0-lground) &
                                - kSOIL(i) *(xT_s_dry0 - xT_g0(i,t,h)) /p_dz_SOIL * lground

    if (xT_s_dry0 .lt. c_TH2Osl) then ! frozen surface

      fRAD_Hw0                  = 0.0
      ETpot0                    = 0.0
      fH2Olg_xu0                = 0.0
      fH2Ol_xd0                 = fH2Ol_ux
      wetfrac_0                 = 0.0
    else

      fRAD_Hw0                  = dRAD &                                ! net radiation (wet) [W / m2 T]
                                - 4.0*p_eps*c_sigma*Ta3*xT_s_wet0 * (2.0-lground) &
                                - kSOIL(i) *(xT_s_wet0 - xT_g0(i,t,h)) /p_dz_SOIL * lground

      ETpot0                    = ( fRAD_Hw0 * desatdT &                ! potential evaporation [m3 H2O / (m2 C * s)]                           REF: Monteith,1965
                                + c_CAIR*(esatAIR-rH2Og_RH(i)*esatAIR) *kH2Og0 ) &
                                / (desatdT+c_gamma) / c_HH2Olg / c_rhoH2Ol

      ! evaporation
      if (rH2Ol_0(i,t,v,h) .gt. p_critD) then

        if (rH2Ol_0(i,t,v,h) / p_dt .lt. ETpot0) then

          fH2Olg_xu0            = rH2Ol_0(i,t,v,h) / p_dt               ! correction for mass balance
          wetfrac_0             = fH2Olg_xu0/max(p_critD,ETpot0)
        else
          fH2Olg_xu0            = ETpot0
          wetfrac_0             = 1.0
        endif
      else
        if (ETpot0 .lt. 0.0) then

          fH2Olg_xu0            = ETpot0
          wetfrac_0             = 1.0
        else
          fH2Olg_xu0            = 0.0
          wetfrac_0             = 0.0
        endif
      endif

      if (ETpot0 .lt. 0.0) then

        if (dewmax(i) .lt. ETpot0) then
    
          fH2Olg_xu0            = ETpot0
          wetfrac_0             = 1.0
        else
          fH2Olg_xu0            = dewmax(i)                             ! correction for mass balance
          wetfrac_0             = fH2Olg_xu0/min(-p_critD,ETpot0)
        endif
      endif

      ! water balance
      rH2Ol_0(i,t,v,h)          = max(0.0, rH2Ol_0(i,t,v,h) &           ! [m3 H2O / m2 C]
                                + fH2Ol_ux * fracWup_0 * p_dt &
                                - fH2Olg_xu0 * p_dt)

      fH2Ol_xd0                 = max(0.0, rH2Ol_0(i,t,v,h) &           ! [m3 H2O / (m2 C * s)]
                                - rmaxH2Ol_0) / p_dt &
                                + fH2Ol_ux * (1.0 - fracWup_0)

      rH2Ol_0(i,t,v,h)          = min(rmaxH2Ol_0, rH2Ol_0(i,t,v,h))

    endif ! surface not frozen

    if (v .eq. 2) then ! canopy

      if (h .eq. 1) then ! stems

        fH2Ol_bd0               = fH2Ol_ub ! no explicit below reservoir, flux per BARK area
      else

        fH2Ol_bd0               = 0.0 ! no input below leaves  
      endif
    else

      fH2Ol_bd0                 = fH2Ol_ub ! no explicit below reservoir
    endif

    xT_s0                       = wetfrac_0*xT_s_wet0 + (1.0-wetfrac_0)*xT_s_dry0

    fRAD_H0                     = wetfrac_0*fRAD_Hw0 + (1.0-wetfrac_0)*fRAD_Hd0

   ! Latent heat
    fQ_ta_L0                    = fH2Olg_xu0 * c_HH2Olg * c_rhoH2Ol     ! Latent heat flux [W / m2 T]

   ! Sensible heat
    fQ_ta_S0                    = c_CAIR * (xT_s_wet0 - xT_a(i)) &      ! Sensible heat flux [W / m2 T]
                                * kH2Og0 * wetfrac_0 &
                                + c_CAIR * (xT_s_dry0 - xT_a(i)) &
                                * kH2Og0 * (1.0 - wetfrac_0)

  endif ! no snow layer at ground

else ! wetland or rocks

  ETpot0                        = 0.0
  fH2Olg_xu0                    = 0.0
  fH2Ol_xd0                     = 0.0
  fH2Ol_bd0                     = 0.0
  xT_s0                         = xT_a(i)
  fRAD_H0                       = 0.0
endif

! soil water balance ---------------------------------------------------

if (v .eq. 1) then ! only for ground layer

  if (xT_s0 .lt. c_TH2Osl) then ! inactivity below zero degrees

    fH2Ol_go                    = fH2Ol_ug

    fH2Ol_gb                    = 0.0

    fH2Olg_ga                   = 0.0

  else

! Potential evaporation

    Rnet                        = fRADs_ad(i)*0.85 +p_eps*fRADl_ad(i) &
                                - p_eps*c_sigma*Ta4
               
    ETpot_v                     = 1.4 *Rnet*desatdT/(desatdT+c_gamma) & ! [m3 H2O / (m2 G * s)]
                                / c_HH2Olg/c_rhoH2Ol

! top soil

    rH2Ol_g1(i,t)               = max(0.0, rH2Ol_g1(i,t) &              ! [m3 H2O / m2 G]
                                + fH2Ol_ug * p_dt ) !&
!                                - fH2Ol_xu )    future work

    percolation                 = max(0.0, rH2Ol_g1(i,t) &              ! [m3 H2O / (m2 G * s)]
                                - p_rmaxH2Ol_g1 ) / p_dt

    rH2Ol_g1(i,t)               = min(p_rmaxH2Ol_g1, rH2Ol_g1(i,t))     ! [m3 H2O / m2 G]

    soilevap                    = min(rH2Ol_g1(i,t) / p_dt, &           ! [m3 H2O / (m2 G * s)]
                                  max(0.0, ETpot_v) )

    rH2Ol_g1(i,t)               = rH2Ol_g1(i,t) - soilevap*p_dt         ! [m3 H2O / m2 G]

! bulk soil

    rH2Ol_g2(i,t)               = rH2Ol_g2(i,t) + percolation*p_dt      ! [m3 H2O / m2 G]

    fH2Ol_go                    = max(0.0, rH2Ol_g2(i,t) &              ! [m3 H2O / (m2 G * s)]
                                - p_rmaxH2Ol_g2 ) / p_dt

    rH2Ol_g2(i,t)               = min(p_rmaxH2Ol_g2, rH2Ol_g2(i,t))     ! [m3 H2O / m2 G]

    rootuptk                    = min( rH2Ol_g2(i,t)/p_dt, p_kH2Ol_sv & ! [m3 H2O / (m2 G * s)]         p_kH2Ol_sv = 5.0e-8
                                * (rH2Ol_g2(i,t)/p_rmaxH2Ol_g2)**2 )

    transpiration               = min( max(0.0,ETpot_v), rootuptk )     ! [m3 H2O / (m2 G * s)]

    rH2Ol_g2(i,t)               = rH2Ol_g2(i,t) - transpiration*p_dt    ! [m3 H2O / m2 G]

    fH2Ol_gb                    = min( rH2Ol_g2(i,t)/p_dt, p_kH2Ol_sb & ! [m3 H2O / (m2 G * s)]         p_kH2Ol_sb = 3.0e-8
                                * (rH2Ol_g2(i,t)/p_rmaxH2Ol_g2)**2 )

    rH2Ol_g2(i,t)               = max(0.0,rH2Ol_g2(i,t)-fH2Ol_gb*p_dt) ! [m3 H2O / m2 G]

    fH2Olg_ga                   = transpiration + soilevap

  endif ! below zero degrees ?

  ! Ground heat
  fQ_tg0                        = kSOIL(i) * (xT_s0 - xT_g0(i,t,h)) &   ! Ground heat flux [W / m2 T]
                                / p_dz_SOIL * lground

  ! Heat balance
  xT_g0(i,t,h)                  = xT_g0(i,t,h) &                        ! balance of ground heat reservoir [K]
                                + fQ_tg0 /CSOIL(i) /p_dz_SOIL *p_dt

endif

! variables for averaging

a0_Ts                           = xT_s0

return
end subroutine land_step

end module libry_land

