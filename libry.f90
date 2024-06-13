
!!!   jvResp25,Rresat/dynAct?


module libry_nova
contains

!***********************************************************************
! NOVA_INIT
!***********************************************************************

subroutine nova_init (nCPts3)
use libry_par
use libry_opt
implicit none

integer :: i,t,v,h,j
integer :: nCPts3

real :: frac_A, frac_B, frac_Bt, frac_W, frac_WB, frac_A2, frac_Wi
real :: Rtest, NCratio

! initialise fields

klife(:,:,:,:,:)                = 0                                     ! switch between dead and alive

rH2Ol_t(:,:,:,:,:)              = 0.0                                   ! thallus water content [m3 H2O / m2 T]

rH2Os_t(:,:,:,:,:)              = 0.0                                   ! thallus ice content [m3 H2O / m2 T]

rH2Ol_b(:,:,:,:,:)              = 0.0                                   ! water content of substrate [m3 H2O / m2 T]

act_state(:,:,:,:,:)            = 0.0                                   ! initial active state [s]

netgrowth(:,:,:,:,:)            = 0.0                                   ! net growth [1 / ts]

areaTH_s(:,:,:,:,:)             = 0.0                                   ! surface cover [ m2 T / m2 V ]

gpp0(:,:,:,:,:)                 = 0.0                                   ! initial GPP


rNd_t(:,:,:,:,:)                = 0.0                                   ! initial N pool [ kg N / m2 T ]

Nsupply(:,:,:,:,:)              = 0.0                                   ! N supply [ kg N / m2 T / month ]

Ndemand(:,:,:,:,:)              = 0.0                                   ! N demand [ kg N / m2 T / month ]

! initialise species parameters

do j = 1,p_nspec ! loop over all species

  o_albedo2(j)                  = vec_o(j,1)*(p_alb_h-p_alb_l)+p_alb_l  ! pure NVV albedo []
                                                                        
  o_zt(j)                       = p_zt_l * exp(vec_o(j,2) &             ! thallus height (including fracAir) [m]
                                * log(p_zt_h / p_zt_l))                 
                                                                        
  o_dzt(j)                      = min(o_zt(j), p_dzt_l*exp(vec_o(j,3) & ! sublayer thickness [m]                                        --> EXP/LITERATURE
                                * log(p_dzt_h / p_dzt_l)) )
                                                                        
  frac_A2                       = 1.0-(1.0-p_fracA2_h) * exp((1.0-vec_o(j,2)) & ! air space increases with height []                    --> EXP/LITERATURE
                                * log((1.0-p_fracA2_l)/(1.0-p_fracA2_h)))
                                                                        
  frac_Bt                       = vec_o(j,4) &                          ! fraction of turgid cells in leaves/lobes []
                                * (p_fracBmax-p_fracBmin)+p_fracBmin
                                                                        
  frac_B                        = p_ratioBW * frac_Bt                   ! fraction of dry matter in fresh biomass []                    --> EXP/LITERATURE
                                                                        
  frac_WB                       = (1.0-p_ratioBW) * frac_Bt             ! cell water []
                                                                        
  frac_W                        = vec_o(j,5) * (1.0 - frac_Bt)          ! water-filled pore space fraction in leaves/lobes []
                                                                        
  frac_A                        = max(0.0, 1.0 -frac_Bt -frac_W)        ! fraction of air in pore space at saturation []
                                                                        
  o_spec_area(j)                = 1.0/(frac_B *c_rhoOrg *0.4 *o_zt(j))  ! specific area (fresh) [m2 T / kg C] -> volume &height must be measured in same state --> CHECK

  o_theta_max(j)                = (frac_Bt*(1.0-p_ratioBW) +frac_W) &   ! water storage capacity [kg H2O / kg C]                        --> CHECK
                                * (1.0-frac_A2) &
                                * o_zt(j) *c_rhoH2Ol *o_spec_area(j)
                                                                        
  LAInvv(j)                     = max(1.0, min(p_LAInvvMax, &           ! number of sublayers = LAI; limited to 12 for light extinction! [] --> CHECK
                                  real( floor( o_zt(j)*(1.0-frac_A2) &
                                / o_dzt(j) ))))
                                                                        
  o_redkB(j)                    = (1.0-(vec_o(j,2)+vec_o(j,3))/2.0) &   ! reduction of boundary layer conductance []                    --> EXP/LITERATURE
                                * (p_redkB_h-p_redkB_l) +p_redkB_l

  o_redkB(j)                    = max(p_critD, o_redkB(j))

  frac_Wi                       = frac_WB / (frac_W+frac_WB)            ! fraction of internal water pool []                            --> EXP/LITERATURE

  o_rs(j)                       = ((1.0-vec_o(j,9)) &                   ! cell wall resistance to H2O [s / m]                           --> EXP/LITERATURE
                                * (p_rs_h-p_rs_l)+p_rs_l) &
                                * frac_Wi                               ! weighted by internal water pool

  o_sat_X(j)                    = frac_Wi*(1.0-0.2) +0.2                ! fraction of internal water pool [] -> water potential         --> EXP/LITERATURE

  o_sat_A(j)                    = (1.0-frac_Wi)*(0.5-0.1) +0.1          ! fraction of internal water pool [] -> full activity           --> EXP/LITERATURE !!

  o_DCO2(j)                     = p_kCO2g_satl &                        ! DCO2 sat [mol / (m2 T * s)]                                   --> CHECK
                                * exp((2.0-vec_o(j,4)-vec_o(j,5))/2.0 &
                                * log(p_kCO2g_sath / p_kCO2g_satl))     

  o_DCO2M(j)                    = p_kCO2gMax_l*exp(vec_o(j,9) &         ! max. CO2 diffusivity [mol/m2/s]
                                * log(p_kCO2gMax_h /p_kCO2gMax_l))

  o_satHph(j)                   = 0.0 !vec_o(j,9)*(p_satHph_h-p_satHph_l) &  ! saturation below which hydrophobicity occurs []
                                !+ p_satHph_l
                                                                        
  o_vcmax_M(j)                  = p_vcmaxM_l * exp(vec_o(j,6) &         ! carboxylation rate of Rubisco (molar vcmax) [1 / s]
                                * log(p_vcmaxM_h / p_vcmaxM_l))         
                                                                        
  o_vomax_M(j)                  = vec_o(j,7) * (p_vomaxM_h &            ! oxygenation rate of Rubisco (molar vomax) [1 / s]
                                - p_vomaxM_l) + p_vomaxM_l              
                                                                        
  o_resp_main(j)                = p_Rref_l * exp(vec_o(j,8) &           ! reference maintenance Respiration  [mol / (m2 T * s)]
                                * log(p_Rref_h / p_Rref_l)) &           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                * LAInvv(j) /p_LAInvvMax &              ! scaling by LAI                                                --> CHECK
                                / o_spec_area(j)                        

  o_spec_Rubisco(j)             = p_RR * o_resp_main(j)                 ! specific Rubisco content [mol / m2 T]                         --> adapt pRR to account for new LAI scaling


!  NCratio                       = vec_o(j,8) * ( p_NC_h - p_NC_l ) &    ! N:C ratio [ (kg N) / (kg C)]
!                                + p_NC_l
!
!  o_Nspec(j)                    = NCratio / o_spec_area(j) !&            ! specific nitrogen content [kg N / m2 T]
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                * LAInvv(j) / p_LAInvvMax               ! scaling by LAI   
!
!  o_spec_Rubisco(j)             = p_RubpN * o_Nspec(j)                  ! Rubisco per area [mol Rubisco / (m2 T)]
!                                                                        
!  o_resp_main(j)                = 4.14e-03 * o_vcmax_M(j) * o_spec_Rubisco(j) & ! empirical relationship
!                                + 6.98e-05 * o_Nspec(j) + 4.69e-07


  o_turnover(j)                 = p_turnover_l * exp(vec_o(j,8) &       ! turnover [1 / yr]
                                * log(p_turnover_h / p_turnover_l))                      
                                                                        
  o_ToptP(j)                    = vec_o(j,10) * (p_ToptP_max &          ! optimum temperature for photosynthesis [K]                    --> make J:V a function of temperature?
                                - p_ToptP_min) + p_ToptP_min            
                                                                        
  o_Q10_resp(j)                 = vec_o(j,11) * (p_Q10R_h &             ! Q10 value of respiration []
                                - p_Q10R_l) + p_Q10R_l                  
                                                                        
  o_Eact_Kc(j)                  = vec_o(j,11) * (p_EaKc_h - p_EaKc_l) & ! enzyme activation energy of Kc [J / mol] ! vec_o(j,12)
                                + p_EaKc_l                              
                                                                        
  o_Eact_Ko(j)                  = vec_o(j,11) * (p_EaKo_h - p_EaKo_l) & ! enzyme activation energy of Ko [J / mol] ! vec_o(j,13)
                                + p_EaKo_l                              

  o_Eact_Vm(j)                  = vec_o(j,11) * (p_EaVm_h - p_EaVm_l) & ! enzyme activation energy of Vcmax [J / mol] ! vec_o(j,12)
                                + p_EaVm_l                              
                                                                        
  o_Eact_Jm(j)                  = vec_o(j,11) * (p_EaJm_h - p_EaJm_l) & ! enzyme activation energy of Vcmax [J / mol] ! vec_o(j,12)
                                + p_EaJm_l                              
                                                                        
  o_CCM(j)                      = 0.9 ! NO CCM   vec_o(j,12)             ! switch for CCM []

  ! define BSC types

!  if (BSCtypes) call libry_defBSC(j)

enddo

!-----------------------------------------------------------------------
! Additional output to evaluate initial trait distributions

if (ltraits .gt. 0 .and. (.not. para .or. (para .and. rank .eq. 0))) then

  open(666, file='test_initTraits', status='replace', action='write' )
   
  do j = 1,p_nspec ! loop over all species
  
    ! constrain ranges (optional)
    !if (o_spec_area(j) .lt. 6.25 .or.  o_spec_area(j) .gt. 62.5 ) klife(:,:,:,:,j) = 0
  
  
    frac_Bt                     = vec_o(j,4) &                          ! fraction of turgid cells in leaves/lobes []
                                * (p_fracBmax-p_fracBmin)+p_fracBmin
                                                                        
    frac_WB                     = (1.0-p_ratioBW) * frac_Bt             ! cell water []
                                                                        
    frac_W                      = vec_o(j,5) * (1.0 - frac_Bt)          ! water-filled pore space fraction in leaves/lobes []
                                                                        
    frac_Wi                     = frac_WB / (frac_W+frac_WB)            ! fraction of internal water pool []                            --> EXP/LITERATURE
  
  
        
    write(666, '(8(E10.3,5X))') 1.0/o_spec_area(j)/0.4 *100.0, &                ! Aspec [m2/kg C] -> STM [mg B/cm2]
                                o_theta_max(j) / o_spec_area(j) *100.0, &       ! WHC [kg H2O/kg C] -> WHC [mg H2O/cm2]
                                frac_Wi, &                                      ! fraction of WHCint 
                                frac_W, &                                       ! fraction of WHCext (ext+int < 1.0)
                                o_zt(j) *100.0, &                               ! Thallus height [cm]
                                1.0/o_spec_area(j)/0.4/o_zt(j), &               ! Bulk density [mg B/cm3]
                                LAInvv(j), &                                    ! LAI []
                                o_resp_main(j)                                  ! base respiration [mol/m2/s]
    !endif
  enddo

  close(666)

  if (ltraits .eq. 2)   stop

endif
!-----------------------------------------------------------------------

! initialise canopy properties

do i = 1,nCPts3

! initialise cover & gpp

  do t = 1,p_ntiles

    if (frac_tile(i,t) .gt. p_critD) then

      do v = 1,p_nvert(t)
        do h = 1,p_nhab(t,v)
          do j = 1,p_nspec
 
            if (t .eq. 1) then
 
              !if ((h .ne. 2) ) klife(i,t,v,h,j) = 1  !no leaves 

              !if ((v .eq. 2 .and. h .eq. 1) ) klife(i,t,v,h,j) = 1  !only stems 

              if ( v .eq. 1 .and. h .eq. 1 ) klife(i,t,v,h,j) = 1  !only forest floor

              !klife(i,t,v,h,j) = 0  !no forest

              !klife(i,t,v,h,j) = 1
 
            elseif (t .eq. 2) then
             
              if (v .eq. 1 .and. h .eq. 1) klife(i,t,v,h,j) = 1

              !klife(i,t,v,h,j) = 0  !no open ground
            endif
 
            if (klife(i,t,v,h,j) .eq. 1) then
 
              areaTH_s(i,t,v,h,j) = fracTH_s_init / real(p_nspec)       ! thallus area[m2] per AVAILABLE AREA[m2]
 
              gpp0(i,t,v,h,j)     = 0.001 / 1.0 /c_MC /p_dt             ! GPP in [mol C / (m2 T * s)]
            endif
          enddo
        enddo
      enddo

    endif

  enddo

enddo ! loop over all grid cells

fracTH_s_crit = fracTH_s_init/real(p_nspec) *fracratiocrit

return
end subroutine nova_init

!***********************************************************************

subroutine nova_update_init (nCPts3)
use libry_par
use libry_opt
implicit none

integer :: i,t,v,h,j
integer :: nCPts3

do i = 1,nCPts3
  do t = 1,p_ntiles

    if (frac_tile(i,t) .gt. p_critD) then

      do v = 1,p_nvert(t)
        do h = 1,p_nhab(t,v)
          do j = 1,p_nspec

            if (areaTH_s(i,t,v,h,j) .le. fracTH_s_crit) then

              klife(i,t,v,h,j)        = 0
              areaTH_s(i,t,v,h,j)     = 0.0

            endif ! check for sufficient cover
          enddo
        enddo
      enddo
    endif
  enddo
enddo

return
end subroutine nova_update_init

!***********************************************************************
! LIBRY_STEP
!***********************************************************************

subroutine nova_step (i,t,v,h)
use libry_par
use libry_opt
implicit none

! temporary variables

integer :: i,t,v,h,j

real    :: kH2Og
real    :: csum

real    :: RH_red, satb
real    :: kH2Og_sM
real    :: kCO2g_t
real    :: ETpot
real    :: gamma2
real    :: grh, cdg, crh
real    :: xT_s_wet, xT_s_dry
real    :: dRAD, fRAD_Hw, fRAD_Hd
real    :: dew, wetfrac
real    :: waterUp0, waterUp, waterUp_b
real    :: fH2Ol_ub2

real    :: vcmaxTo, vcmax25
real    :: HcpO,Hcp
real    :: Resp, IR, JeT, JeB
real    :: D1l, D2l, al, bl, cl, discl, xl, convLiq
real    :: D1c, D2c, ac, bc, cc, discc, xc, K0

real    :: ngsum, wgtsum, expsum, hsum
real    :: retreat, disturbance

real    :: resatResp
real    :: frac_A, frac_B, frac_W

!!!!!!!!!!!!!!!!!! CHECK
real    :: dnW1, v_wat, v_ice, ksatORG, KE
real    :: dWmelt, dWfreeze, CORG, kORG

real    :: waterPrev, Nleaching, Ninput, nitro_lim

integer :: kawfile, kawfileb

character(12) :: fdza,fdzb

! surface coverage
csum                            = 0.0

!-----------------------------------------------------------------------
! Start loop over all species in a grid cell - I -
!-----------------------------------------------------------------------

do j = 1,p_nspec

! Check if lichen is alive

  if (klife(i,t,v,h,j) .eq. 1) then

! *** WATER and ENERGY BALANCE ***

! Boundary layer conductance

    kH2Og                       = fAIR_s(i) / (75.0 &                   ! [m / s]
                                * 8.0**roughlen(i,t,v,h)) &
                                * o_redkB(j)

! Maximum lichen water storage capacity:

    rmaxH2Ol_t(j)               = max(0.0, o_theta_max(j) / c_rhoH2Ol & ! [m3 H2O / m2 T]
                                / o_spec_area(j) )!&
!                                - rH2Os_t(i,t,v,h,j) )


! Lichen water saturation:

    sat(j)                      = rH2Ol_t(i,t,v,h,j) / rmaxH2Ol_t(j)    ! []

    if (o_CCM(j) .lt. p_FCCM) then                                      ! lichen with CCM
      sat_act                   = p_sat_act0CCM
    else
      sat_act                   = p_sat_act0
    endif

! Lichen water potential:

    xH2Ol                       = min(0.0, max( -50.0, &                ! [MPa = MJ / m3]
                                  15.0 * (1.0 - o_sat_X(j) &
                                / max(p_critD, sat(j))))) 

    RH_red                      = xH2Ol * 1000.0 * c_MH2O / c_Rgas &    ! correction for capillary forces                                       REF: Nikolov,1995
                                / xT_a(i)

! Below thallus water potential:

    if (v .eq. 2 .and. h .eq. 2) then

      satb                      = -1.0 ! test for consistency
      xH2Ol_b                   = -1000.0 ! ensure no uptake for leaves
    else
      satb                      = rH2Ol_b(i,t,v,h,j)/rmaxH2Ol_b(i,t,v,h) ! []

      xH2Ol_b                   = max( -45.0, 5.0 &                    ! [MPa = MJ / m3]
                                * (1.0 - 1.0/max(p_critD,satb)) )
    endif

! Saturation water vapour pressure

    esatAIR2                    = exp(RH_red) * esatAIR

    desatdT2                    = exp(p_esatAIR1 * zT_a / (p_esatAIR2 & ! slope of saturation vapour pressure vs temperature relationship []
                                + zT_a) + RH_red) &
                                * (p_esatAIR1 * p_esatAIR2 *p_esatAIR3 &
                                / ((p_esatAIR2 + zT_a) &
                                * (p_esatAIR2 + zT_a)) - RH_red/xT_a(i))

! Lichen CO2 diffusion:

    kCO2g_t                     = (o_DCO2M(j) - o_DCO2(j)) &            ! [mol / (m2 T * s)]                                                   REF: Cowan,1992; Williams,1998
                                * (1.0-sat(j))**7 + o_DCO2(j)

!        kH2Og_sM                = kH2Og * rH2Og_RH(i) * esatAIR2 &      ! transform kH2Og [m/s] into [mol/m2/s]
!                                / c_Rgas / xT_a(i)

    kCO2g                       = kCO2g_t !* o_int_area(j) !1.0/(1.0/kCO2g_t + 1.0/(kH2Og_sM/1.6))

! Lichen H2O diffusion:

    gamma2                      = c_gamma * (1.0+o_rs(j)*kH2Og)


!!!!!!!! CHECK
!
!    dnW1                        = rH2Ol_t(i,t,v,h,j) +rH2Os_t(i,t,v,h,j)
!     
!    if (dnW1 .gt. p_critD) then
!      v_wat                     = rH2Ol_t(i,t,v,h,j) / o_zt(j)          ! volumetric water content []
!      v_ice                     = rH2Os_t(i,t,v,h,j) / o_zt(j)          ! volumetric ice content []
!    else
!      v_wat                     = 0.0
!      v_ice                     = 0.0
!    endif
!
!    ksatORG                     = (p_kORG**(1.0-(v_wat+v_ice)))*(p_kH2O**v_wat)*(p_kICE**v_ice)  ! [W / m / K]
!
!    KE                          = min(1.0, (v_wat+v_ice)/o_prs(j))      ! Kersten number []
!
!! Thermal conductivity
!
!    kORG                        = ksatORG *KE + (1.0-KE)*p_kdry         ! [W / m / K]
!
!! NVV heat capacity
!
!    CORG                        = (1.0-o_prs(j))*p_CORG +v_wat*c_CH2O +v_ice*c_CICE  ! [J / m3 / K]


! Lichen surface temperature:

    if (dsnow .ge. p_H2Os_crit .and. v .eq. 1) then             ! snow layer too thick

      ETpot                     = 0.0
      fH2Olg_xu(j)              = 0.0
      fH2Ol_bx(j)               = 0.0
      fH2Ol_bd(j)               = 0.0
      fH2Ol_xd(j)               = 0.0
      xT_s_dry                  = min(c_TH2Osl - 0.1, xT_a(i))
      xT_s(j)                   = xT_s_dry
      fRAD_H(j)                 = 0.0
      fQ_ta_L(j)                = 0.0
      fQ_ta_S(j)                = 0.0
      Ninput                    = 0.0
      Nleaching                 = 0.0
    else

      grh                       = gamma2 / kH2Og
      cdg                       = c_CAIR*(desatdT2+gamma2)
 
      xT_s_wet                  = ( grh * (fRADs_ad(i)*fracRADs*(1.0-o_albedo2(j)) & ! surface temperature [K]                                  REF: Monteith,1981
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) &
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground)&
                                + kSOIL(i)/p_dz_SOIL*xT_g(i,t,h,j) * lground ) &
                                + xT_a(i)*cdg &
                                - c_CAIR*(esatAIR2-rH2Og_RH(i)*esatAIR) ) &
                                / ( grh*(4.0*p_eps*c_sigma*Ta3 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL * lground) + cdg )
 
      crh                       = 1.0/(c_CAIR*kH2Og)
 
      xT_s_dry                  = ( xT_a(i) &
                                + crh*(fracRADs*(1.0-o_albedo2(j))*fRADs_ad(i) &
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) &
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL*xT_g(i,t,h,j) * lground) ) &
                                / ( 1.0 + crh*(4.0*p_eps*c_sigma*Ta3 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL * lground) )
 
! Net radiation

      dRAD                      = fRADs_ad(i)*fracRADs*(1.0-o_albedo2(j)) & ! downwelling radiation [W / m2 T]
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) &
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground) 
 
      fRAD_Hd                   = dRAD &                                ! net radiation [W / m2 T]
                                - 4.0*p_eps*c_sigma*Ta3*xT_s_dry * (2.0-lground) &
                                - kSOIL(i) *(xT_s_dry - xT_g(i,t,h,j)) /p_dz_SOIL *lground

      if (xT_s_dry .lt. c_TH2Osl) then ! frozen surface

        fRAD_Hw                 = 0.0
        ETpot                   = 0.0
        fH2Olg_xu(j)            = 0.0
        fH2Ol_bx(j)             = 0.0

        if (v .eq. 2 .and. h .eq. 2) then
          fH2Ol_bd(j)           = 0.0
        else
          fH2Ol_bd(j)           = fH2Ol_ub
        endif

        fH2Ol_xd(j)             = fH2Ol_ux
        wetfrac                 = 0.0
        Ninput                  = 0.0
        Nleaching               = 0.0
      else

        fRAD_Hw                 = dRAD &                                ! net radiation [W / m2 T]
                                - 4.0*p_eps*c_sigma*Ta3*xT_s_wet * (2.0-lground) &
                                - kSOIL(i) *(xT_s_dry - xT_g(i,t,h,j)) /p_dz_SOIL *lground

        ETpot                   = ( fRAD_Hw * desatdT2 &                ! total evaporation [m3 H2O / (m2 T * s)]                               REF: Monteith,1965
                                + c_CAIR*(esatAIR2-rH2Og_RH(i)*esatAIR) *kH2Og ) &
                                / (desatdT2+gamma2) /c_HH2Olg /c_rhoH2Ol

        ! dewfall
        if (ETpot .lt. 0.0) then

          if (dewmax(i) .lt. ETpot) then
 
            dew                 = ETpot
            wetfrac             = 1.0
          else
            dew                 = dewmax(i)                             ! correction for mass balance
            wetfrac             = dew /min(-p_critD,ETpot)              ! sensible heat flow, air has dried out 
          endif
        else
          dew                   = 0.0
        endif

        ! water uptake from below
        waterUp_b               = max( 0.0, (xH2Ol_b-xH2Ol)/50.0 &      ! [ m/s ]
                                * p_ksatBark)

        waterUp_b               = min(waterUp_b, rH2Ol_b(i,t,v,h,j)/p_dt)

        fH2Ol_bx(j)             = min(waterUp_b, &
                                  (rmaxH2Ol_t(j)-rH2Ol_t(i,t,v,h,j))/p_dt)

        !if (v .eq. 2 .and. h .eq. 2) fH2Ol_bx(j) = 0                    ! no uptake on leaves
        fH2Ol_bx(j) = 0                                                 ! no uptake from below

        ! rainfall uptake efficiency
        waterUp0                = fH2Ol_ux * p_fracWup - dew

        ! hydrophobicity
        if (sat(j) .lt. o_satHph(j)) then

          waterUp               = min(waterUp0, &
                                  rmaxH2Ol_t(j)/c_hour *p_redHph)
        else

          waterUp               = waterUp0
        endif

        ! actual evaporation
        if (ETpot .ge. 0.0) then

          if (rH2Ol_t(i,t,v,h,j) / p_dt + waterUp .le. ETpot) then

            fH2Olg_xu(j)        = rH2Ol_t(i,t,v,h,j) / p_dt +waterUp    ! correction for mass balance
            wetfrac             = fH2Olg_xu(j)/max(p_critD,ETpot)
          else

            fH2Olg_xu(j)        = ETpot
            wetfrac             = 1.0
          endif
        else

          fH2Olg_xu(j)          = dew
        endif

        ! Water balance    
        waterPrev               = max(p_critD, rH2Ol_t(i,t,v,h,j))      ! save old water reservoir

        rH2Ol_t(i,t,v,h,j)      = max(0.0, rH2Ol_t(i,t,v,h,j) &         ! water reservoir [m3 H2O / m2 T]
                                + waterUp * p_dt &                      ! rainfall/throughfall + dew
                                + fH2Ol_bx(j) * p_dt &                  ! from below
                                - max(0.0, fH2Olg_xu(j)) * p_dt)
      
        fH2Ol_xd(j)             = max(0.0, rH2Ol_t(i,t,v,h,j) &         ! runoff / throughfall [m3 H2O / (m2 T * s)]
                                - rmaxH2Ol_t(j)) / p_dt

        fH2Ol_xd(j)             = fH2Ol_xd(j) + waterUp0 - waterUp &
                                + fH2Ol_ux * (1.0 - p_fracWup)

        rH2Ol_t(i,t,v,h,j)      = min(rmaxH2Ol_t(j),rH2Ol_t(i,t,v,h,j)) ! water reservoir II

        if (v .eq. 2) then ! canopy

          if (h .eq. 1) then ! stems

            fH2Ol_ub2           = fH2Ol_ub
          else

            fH2Ol_ub2           = 0.0 ! no input below leaves
          endif
        else

          fH2Ol_ub2             = fH2Ol_ub
        endif

        rH2Ol_b(i,t,v,h,j)      = max(0.0, rH2Ol_b(i,t,v,h,j) &         ! below reservoir (bark/substrate) [m3 H2O / m2 T]
                                + fH2Ol_ub2 * p_dt &
                                - fH2Ol_bx(j) * p_dt)

        fH2Ol_bd(j)             = max(0.0, rH2Ol_b(i,t,v,h,j) &
                                - rmaxH2Ol_b(i,t,v,h)) / p_dt

        rH2Ol_b(i,t,v,h,j)      = min(rmaxH2Ol_b(i,t,v,h),rH2Ol_b(i,t,v,h,j))


        Ninput                  = waterUp * p_dt * Nconc                ! uptake of N deposition

        Nleaching               = min( waterPrev, fH2Ol_xd(j)*p_dt ) &  ! leaching of N from pool
                                * rNd_t(i,t,v,h,j) / waterPrev 

      endif ! surface not frozen

      xT_s(j)                   = wetfrac*xT_s_wet + (1.0-wetfrac)*xT_s_dry
 
      fRAD_H(j)                 = wetfrac*fRAD_Hw + (1.0-wetfrac)*fRAD_Hd
 
      ! Latent heat
      fQ_ta_L(j)                = fH2Olg_xu(j) * c_HH2Olg * c_rhoH2Ol   ! Latent heat flux [W / m2 T]

      ! Sensible heat
      fQ_ta_S(j)                = c_CAIR * (xT_s_wet - xT_a(i)) &       ! Sensible heat flux [W / m2 T]
                                * kH2Og * wetfrac &
                                + c_CAIR * (xT_s_dry - xT_a(i)) &
                                * kH2Og * (1.0 - wetfrac)

    endif ! no snow layer at ground


!!!!!!!!!!!!! CHECK
!
!    ! Freezing and thawing
! 
!    if (xT_s(j) .lt. c_TH2Osl .and. rH2Ol_t(i,t,v,h,j) .gt. 0.0) then 
!    
!      dWmelt                    = 0.0
!      dWfreeze                  = min(rH2Ol_t(i,t,v,h,j), &
!                                  rmaxH2Ol_t(j)*CORG*(c_TH2Osl-xT_s(j))/(c_HH2Ols*c_rhoH2Ol))
! 
!      xT_s(j)                   = min(xT_s(j) + dWfreeze*c_HH2Ols*c_rhoH2Ol/(CORG*rmaxH2Ol_t(j)), c_TH2Osl)
! 
!      rH2Ol_t(i,t,v,h,j)        = max(0.0, rH2Ol_t(i,t,v,h,j) - dWfreeze)
!      rmaxH2Ol_t(j)             = max(0.0, rmaxH2Ol_t(j) - dWfreeze)            ! only needed for output below
!      rH2Os_t(i,t,v,h,j)        = min(o_theta_max(j) /c_rhoH2Ol /o_spec_area(j), rH2Os_t(i,t,v,h,j) + dWfreeze)
! 
!    else if (xT_s(j) .gt. c_TH2Osl .and. rH2Os_t(i,t,v,h,j) .gt. 0.0) then
! 
!      dWfreeze                  = 0.0
!      dWmelt                    = min(rH2Os_t(i,t,v,h,j), &
!                                  rmaxH2Ol_t(j)*CORG*(xT_s(j)-c_TH2Osl)/(c_HH2Ols*c_rhoH2Os))
! 
!      xT_s(j)                   = max(xT_s(j) - dWmelt*c_HH2Ols*c_rhoH2Os/(CORG*rmaxH2Ol_t(j)), c_TH2Osl)
! 
!      rmaxH2Ol_t(j)             = min(o_theta_max(j) /c_rhoH2Ol /o_spec_area(j), rmaxH2Ol_t(j) + dWmelt)
!      rH2Ol_t(i,t,v,h,j)        = min(rmaxH2Ol_t(j), rH2Ol_t(i,t,v,h,j) + dWmelt)
!      rH2Os_t(i,t,v,h,j)        = max(0.0, rH2Os_t(i,t,v,h,j) - dWmelt)
!    else
!      dWfreeze                  = 0.0
!      dWmelt                    = 0.0
!    endif
! 
!    ! Ground heat
!    fQ_tg(j)                    = 1.0/ (1.0/kORG *o_zt(j)/(o_zt(j)+p_dz_SOIL) & ! Ground heat flux [W / m2 T]
!                                + 1.0/kSOIL(i) *p_dz_SOIL/(o_zt(j)+p_dz_SOIL) ) & 
!                                * (xT_s(j) - xT_g(i,t,h,j)) &
!                                / p_dz_SOIL * lground
!
!    ! Heat balance
!    xT_g(i,t,h,j)               = xT_g(i,t,h,j) + fQ_tg(j) &            ! balance of ground heat reservoir [K]
!                                / ( CORG *o_zt(j)/(o_zt(j)+p_dz_SOIL) &
!                                + CSOIL(i) *p_dz_SOIL/(o_zt(j)+p_dz_SOIL) ) &
!                                /p_dz_SOIL *p_dt

! OLD HEAT DIFFUSION SCHEME

    ! Ground heat
    fQ_tg(j)                    = kSOIL(i) * (xT_s(j) - xT_g(i,t,h,j))& ! Ground heat flux [W / m2 T]
                                / p_dz_SOIL * lground

    ! Heat balance
    xT_g(i,t,h,j)               = xT_g(i,t,h,j) &
                                + fQ_tg(j) /CSOIL(i) /p_dz_SOIL *p_dt   ! balance of ground heat reservoir [K]


! *** NITROGEN BALANCE ***

    rNd_t(i,t,v,h,j)            = max(0.0, rNd_t(i,t,v,h,j) + Ninput - Nleaching) ! update N pool [ kg N / m2 T ]

    Nsupply(i,t,v,h,j)          = Nsupply(i,t,v,h,j) + rNd_t(i,t,v,h,j) ! accumulate N supply [ kg N / m2 T / month ]

! *** CARBON BALANCE ***

    ! Photosynthesis model
    KcfT                        = exp((xT_s(j) - o_ToptP(j)) &          ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
                                * o_Eact_Kc(j) &
                                / (o_ToptP(j) * c_Rgas * xT_s(j)))
                                
    KofT                        = exp((xT_s(j) - o_ToptP(j)) &          ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
                                * o_Eact_Ko(j) &
                                / (o_ToptP(j) * c_Rgas * xT_s(j)))
                                
    VmfT                        = exp((xT_s(j) - o_ToptP(j)) &          ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
                                * o_Eact_Vm(j) &
                                / (o_ToptP(j) * c_Rgas * xT_s(j)))
                                
    JmfT                        = exp((xT_s(j) - 298.15) &              ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
                                * o_Eact_Jm(j) &
                                / (xT_s(j) * c_Rgas * 298.15))
                                
    vcmaxTo                     = o_vcmax_M(j) * o_spec_Rubisco(j)      ! vcmax [mol / (m2 T * s)] *** at T opt ***
                                
    vcmax(j)                    = vcmaxTo * VmfT                        ! vcmax at current Ts
                                
    vcmax25                     = vcmaxTo*exp((298.15 - o_ToptP(j)) &   ! vcmax at standard temperature
                                * o_Eact_Vm(j) &
                                / (o_ToptP(j) * c_Rgas * 298.15))
                                
    jmax(j)                     = vcmax25 * 2.1 * JmfT                  ! jmax [mol / (m2 T * s)] at current Ts                                 REF: Wullschleger,1993
                                
    KcM                         = p_KcM1 * o_vcmax_M(j)**p_KcM2         ! Michaelis-Menten-Constant for CO2 [muM]                               REF: Savir,2009
                                
    KoM                         = o_vomax_M(j) / (p_KoM1 &              ! Michaelis-Menten-Constant for O2 [muM]                                REF: Savir,2009
                                * (o_vcmax_M(j) / KcM)**p_KoM2)
                                
    Kc                          = KcM * 0.001 * KcfT                    ! temperature correction [mol / m3]
                                
    Ko                          = KoM * 0.001 * KofT                    ! temperature correction [mol / m3]
                                
    HcpO                        = p_sO2gl*exp(1700.0*(1.0/xT_s(j) - 1.0/298.15))
                                
    Gs                          = 0.5*rO2g_a/1.0E6*HcpO*1000.0 &        ! gammastar [mol / m3]                                                  REF: Farquhar,1981
                                * o_vomax_M(j)/ o_vcmax_M(j) * Kc/Ko
 
    ! Respiration:
    Rspec                       = o_resp_main(j) &                      ! [mol / (m2 T * s)]                                                    REF: Kruse,2010
                                * o_Q10_resp(j) &
                                **((xT_s(j) - o_ToptP(j)) / 10.0)

    ! Activity
    act(j)                      = max( 0.0, min(1.0, &                  ! activity as a function of water content []
                                  (sat(j) - sat_act) &
                                / max(p_critD, o_sat_A(j) -sat_act) ) )

    if (act(j) .le. p_critD) then

      act(j) = 0.0

      act_state(i,t,v,h,j)      = 0.0
    endif

    if (act_state(i,t,v,h,j) .lt. 0.9 .and. act(j) .gt. p_critD) then

      act_state(i,t,v,h,j)      = 1.0

      !resatResp                 = Rspec*0.5 *7200.0/p_dt        ! 50% increased respiration rate for 2 hours, converted to 1 time step
      resatResp                 = 0.0                           ! no resaturation respiration
    else

      resatResp                 = 0.0
    endif

    ! Respiration, NPP and litterfall

    if (xT_s(j) .lt. c_TH2Osl) then ! inactivity below zero degrees

      act(j)                    = 0.0
      fCcg_M(j)                 = 0.0
      fCcb(j)                   = 0.0
      fCcg_G(j)                 = 0.0
!      fCbo(j)                   = 0.0
      fCbo(j)                   = o_turnover(j) /c_year /c_MC &         ! litterfall of photobiont [mol / (m2 T * s)]
                                / o_spec_area(j)
      fCO2gc_L(j)               = 0.0
      fCO2gc_W(j)               = 0.0
      fCO2gc(j)                 = 0.0
      gpp0(i,t,v,h,j)           = 0.0
      CO2_p(j)                  = rCO2g_a
    else

      fCcg_M(j)                 = Rspec * act(j) +resatResp             ! maintenance resp. [mol / (m2 T * s)]

      fCcb(j)                   = gpp0(i,t,v,h,j) - fCcg_M(j)           ! NPP of photobiont [mol / (m2 T * s)]

      if (fCcb(j) .gt. 0.0) fCcb(j) = fCcb(j) * p_NCcb ! growth efficiency

      fCcg_G(j)                 = (1.0-p_NCcb)*max(fCcb(j),0.0)/p_NCcb  ! growth resp. [mol / (m2 T * s)]

      fCbo(j)                   = o_turnover(j) /c_year /c_MC &         ! litterfall of photobiont [mol / (m2 T * s)]
                                / o_spec_area(j) ! * act(j)

      Resp                      = fCcg_M(j) + fCcg_G(j)                 ! total respiration

      IR                        = fRADs_ad(i) &                         ! potential electron flow per m2 T
                                * (1.0-o_albedo2(j)) * p_conv_PAR 

      if (o_CCM(j) .lt. p_FCCM) IR = IR * p_conv_CCM                    ! lichen with CCM

      Je                        = (1.0-exp(-p_beer_s*LAInvv(j))) &      ! electron flow per m2 T
                                * IR * jmax(j) &       
                                / (IR + p_cvxL*(LAInvv(j)/p_LAInvvMax)*2.1*jmax(j)) &
                                * act(j)

!      Je                        = IR * jmax(j) &                        ! electron flow per m2 T
!                                / (IR + 3.5*2.1*jmax(j)) &
!                                * act(j)

!      JeT                       = min(IR/2.1, jmax(j))
!
!      JeB                       = min(IR*exp(-p_beer_s *LAInvv(j)) /2.1, jmax(j))
!
!      Je                        = (JeT*(1.0-p_xcano) + JeB*p_xcano) &   ! electron flow per m2 T
!                                * act(j)

!      Je                        = IR * jmax(j) &                        ! electron flow per m2 T
!                                / (IR + o_cvxL(j)*2.1*jmax(j)) &
!                                * act(j)
 
      VA                        = vcmax(j) * act(j)                     ! vcmax per m2 T
 
      ! gpp(Ci[mol/l]) - R = D * (Ca[atm] - Cp[atm])
      !                                     Ci[mol/l] / H [mol/l/atm]   (Henry's law)
      !     Ci[mol/m3]                      Ci[mol/m3] / 1000 [l/m3] / H [mol/l/atm]
      !
      ! CCM:
      !     Ci[mol/m3]                      Ci[mol/m3] / 1000 [l/m3] * 22.7 [l/mol*atm] / p_NCCM

      Hcp                       = p_sCO2gl*exp(2400.0*(1.0/xT_s(j) - 1.0/298.15))

      if (o_CCM(j) .lt. p_FCCM) then ! with CCM
        convLiq                 =  p_NCCM / 22.7
      else
        convLiq                 =  Hcp
      endif

      ! light-limited rate                                                                                                                      ! REF: Farquhar,1981
      D1l                       = 4.0 * kCO2g /convLiq /1000.0          ! Ci[mol/m3] to Cp[atm]
      D2l                       = 4.0 * kCO2g * rCO2g_a / 1.0E6         ! [ppm] to [atm]
      al                        = -D1l
      bl                        = D2l - 2.0*D1l*Gs - Je + 4.0*Resp
      cl                        = (2.0*D2l + Je  + 8.0*Resp)*Gs

      discl                     = sqrt(bl*bl  -  4.0 * al  * cl)
      xl                        = (-bl -  discl)/ (2.0  *  al)  
      fCO2gc_L(j)               = Je * (xl - Gs) / (4.0 * xl + 8.0 * Gs) 
                                                                                                              
      ! CO2-limited rate
      D1c                       = kCO2g /convLiq /1000.0
      D2c                       = kCO2g * rCO2g_a / 1.0E6
      K0                        = Kc * (1.0 + rO2g_a/1.0E6 *HcpO*1000.0 / Ko) 
      ac                        = -D1c
      bc                        = D2c - D1c*K0 - VA + Resp
      cc                        = D2c*K0 + VA*Gs + Resp*K0

      discc                     = sqrt(bc*bc - 4.0 * ac * cc) 
      xc                        = (-bc - discc) / (2.0 * ac) 

      fCO2gc_W(j)               = VA * (xc - Gs) / (xc + K0) 

      if (fCO2gc_L(j) .lt. fCO2gc_W(j)) then
        fCO2gc(j)               = fCO2gc_L(j)
        CO2_p(j)                = xl /convLiq /1000.0 *1.0E6            ! [atm] to [ppm]
      else
        fCO2gc(j)               = fCO2gc_W(j)
        CO2_p(j)                = xc /convLiq /1000.0 *1.0E6
      endif

      gpp0(i,t,v,h,j)           = fCO2gc(j)

    endif ! above zero degrees

    ! nitrogen limitation of growth

    Ndemand(i,t,v,h,j)          = Ndemand(i,t,v,h,j) &                  ! N demand per growth [kg N /m2 T /month]
                                + fCcb(j) * p_dt * c_MC &
                                * o_spec_area(j) * o_Nspec(j)

    ! net growth (NPP - litterfall) per m2 T in one timestep

    netgrowth(i,t,v,h,j)        = netgrowth(i,t,v,h,j) &                ! [1 / month]
                                + (fCcb(j)*p_dt - fCbo(j)*p_dt) &
                                * o_spec_area(j) * c_MC

    ! accumulate cover
    csum                        = csum + areaTH_s(i,t,v,h,j)

!        if (BSCtypes .and. NOHONO) call libry_fNOHONO(j,sat(j),act(j))

  endif ! check for survival

enddo ! End of loop over all species - I -

!    if (BSCtypes) call libry_covsBSC(i,i2)

!-----------------------------------------------------------------------
! Execute once per month
!-----------------------------------------------------------------------

if (day .eq. dpm .and. ts .eq. tspd) then

!-----------------------------------------------------------------------
! Start loop over all species - II -
!-----------------------------------------------------------------------

  hsum                          = 0.0
!  ngsum                         = 0.0

  do j = 1,p_nspec

! Determine nitrogen limitation

    if (Ndemand(i,t,v,h,j) .gt. p_critD) then 
      if (Ndemand(i,t,v,h,j) .le. Nsupply(i,t,v,h,j)) then                  ! when N demand is smaller than supply
        nitro_lim               = 1.0
      else
        nitro_lim               = Nsupply(i,t,v,h,j) / Ndemand(i,t,v,h,j)
      endif
    else
      nitro_lim                 = 1.0
    endif

    Ndemand(i,t,v,h,j)          = 0.0
    Nsupply(i,t,v,h,j)          = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SWITCHED OFF    netgrowth(i,t,v,h,j)        = netgrowth(i,t,v,h,j) * nitro_lim 

! Sum up net growth of last month (only positive)

    if (klife(i,t,v,h,j) .eq. 1) hsum = hsum + o_zt(j)
!    if (klife(i,t,v,h,j) .eq. 1) ngsum = ngsum + max(0.0,netgrowth(i,t,v,h,j))

  enddo

! Start loop over all species - IIb -

  if (hsum .gt. p_critD) then
!  if (ngsum .gt. p_critD) then

    wgtsum                      = 0.0

    do j = 1,p_nspec

      if (klife(i,t,v,h,j) .eq. 1) then

        wgtspec(j)              = 1.0 &                                 ! equal weights at low cover (no competition)
                                * ( (o_zt(j)/hsum)**csum )**2           ! competition under low disturbance -> weighting by growth height

!        wgtspec(j)              = 1.0 &                                 ! equal weights at low cover (no competition)
!                                * ( ( max(0.0,netgrowth(i,t,v,h,j)) &                                 
!                                / ngsum )**csum )**2                    ! competition under low disturbance -> weighting by net growth

!        wgtspec(j)              = 1.0                                   ! equal weights (neutral model)

        wgtsum                  = wgtsum + wgtspec(j)
      endif
    enddo ! End of loop over all species - IIb -

! Start loop over all species - IIc -

    do j = 1,p_nspec

      if (klife(i,t,v,h,j) .eq. 1) then

! *** Cover change I ***

        expansion(j)            = max(0.0, min( &                       ! [m2 T / m2 V / month]
                                  netgrowth(i,t,v,h,j) *areaTH_s(i,t,v,h,j) & ! mass balance constraint
                                * p_NCbt, &                             ! establishment
                                  (1.0 - csum) &                        ! available area
                                * wgtspec(j) /wgtsum) )                 ! competition 

      else

        expansion(j)            = 0.0

      endif ! check for survival
    enddo ! End of loop over all species - IIc -

  else 
    expansion(:)                = 0.0
  endif

  expsum                        = sum(expansion(:))

  if (expsum .gt. p_critD) then

    expansion(:)                = expansion(:) * min(1.0, &
                                  (1.0-csum) / expsum)
  else
    expansion(:)                = 0.0
  endif

!-----------------------------------------------------------------------
! Start loop over all species - III -
!-----------------------------------------------------------------------

  do j = 1,p_nspec

! Check if lichen is alive

    if (klife(i,t,v,h,j) .eq. 1) then

! *** Cover change II ***

      retreat                   = min( max(0.0, -netgrowth(i,t,v,h,j)) &   ! [m2 T / m2 V / month]
                                * areaTH_s(i,t,v,h,j), &
                                  areaTH_s(i,t,v,h,j) )
                            
      disturbance               = 1.0/tauD(i,t,v,h) * areaTH_s(i,t,v,h,j)  ! [m2 T / m2 V / month]
                              
      areaTH_s(i,t,v,h,j)       = max(0.0, areaTH_s(i,t,v,h,j) &           ! [m2 T / m2 V]
                                + expansion(j) &
                                - disturbance &
                                - retreat )
                            
      netgrowth(i,t,v,h,j)      = 0.0

! Set lichen to dead if cover is too low

      if (areaTH_s(i,t,v,h,j) .le. fracTH_s_crit) then

        klife(i,t,v,h,j)        = 0
        areaTH_s(i,t,v,h,j)     = 0.0

      endif ! check for sufficient cover

    endif ! check for survival

  enddo ! End of loop over all species - III -

endif ! End of execute once per month


! *** Average strategies ***

! properties which affect other habitats/levels/tiles

as_fH2Ol_td                     = 0.0

as_fH2Ol_bd                     = 0.0

as_areaTH_s                     = csum ! total area of all strategies

! properties which are only needed for output

if (writeout) then

  ! initialise accumulated variables with zero
  as_rCO2d                      = 0.0
  as_rCb                        = 0.0
  as_rH2Ol_t                    = 0.0
  as_rH2Os_t                    = 0.0
  as_rmaxH2Ol_t                 = 0.0
  as_act                        = 0.0
  as_actB                       = 0.0
  as_fCO2gc                     = 0.0
  as_fCcg                       = 0.0
  as_fCcb                       = 0.0
  as_fCcb_l                     = 0.0
  as_fCcb_c                     = 0.0
  as_fCbo                       = 0.0
  as_frgr                       = 0.0
  as_fH2Ogl_ut                  = 0.0
  as_fH2Olg_tu                  = 0.0
  as_fH2Ol_bx                   = 0.0
  as_Ts                         = 0.0
  as_Tg                         = 0.0
  as_H                          = 0.0
  as_G                          = 0.0
  as_E                          = 0.0
  as_C                          = 0.0
  as_EB                         = 0.0
endif

!-----------------------------------------------------------------------
! Start loop over all species - IV -
!-----------------------------------------------------------------------

! loop over all strategies
do j = 1,p_nspec

  ! Check if lichen is alive

  if (klife(i,t,v,h,j) .eq. 1) then

    as_fH2Ol_td                 = as_fH2Ol_td + fH2Ol_xd(j) *1000.0*c_year *areaTH_s(i,t,v,h,j)

    as_fH2Ol_bd                 = as_fH2Ol_bd + fH2Ol_bd(j) *1000.0*c_year *areaTH_s(i,t,v,h,j)
  endif
enddo

!-----------------------------------------------------------------------
! Start loop over all species - V -
!-----------------------------------------------------------------------

if (writeout) then

  ! loop over all strategies
  do j = 1,p_nspec
   
    ! Check if lichen is alive
   
    if (klife(i,t,v,h,j) .eq. 1) then
   
      ! intensive variables are weighted by cover
   
      cweight                   = areaTH_s(i,t,v,h,j) / csum
   
      as_rCO2d                  = as_rCO2d + CO2_p(j) *cweight
                                    
      as_rCb                    = as_rCb + 1.0/o_spec_area(j)*1000.0 *areaTH_s(i,t,v,h,j)
                                    
      as_rH2Ol_t                = as_rH2Ol_t + sat(j) *cweight
                                    
      as_rH2Os_t                = as_rH2Os_t + rH2Os_t(i,t,v,h,j)/(rmaxH2Ol_t(j)+rH2Os_t(i,t,v,h,j)) *cweight

      as_rmaxH2Ol_t             = as_rmaxH2Ol_t + rmaxH2Ol_t(j) *1000.0 *areaTH_s(i,t,v,h,j)

      as_act                    = as_act + act(j) *cweight
                                    
      as_actB                   = as_actB + ceiling(act(j)) *cweight
                                    
      as_fCO2gc                 = as_fCO2gc + fCO2gc(j) *c_MC*c_year*1000.0 *areaTH_s(i,t,v,h,j)
                                
      as_fCcg                   = as_fCcg + (fCcg_M(j) + fCcg_G(j)) *c_MC*c_year*1000.0 *areaTH_s(i,t,v,h,j)
                                
      as_fCcb                   = as_fCcb + fCcb(j) *c_MC*c_year*1000.0 *areaTH_s(i,t,v,h,j)
                                
      as_frgr                   = as_frgr +(fCcb(j)*p_dt-fCbo(j)*p_dt)*o_spec_area(j)*c_MC  /p_dt*c_year *cweight


      as_fCcb_l                 = as_fCcb_l + (fCO2gc_L(j) -fCcg_M(j)-fCcg_G(j)) *c_MC*c_year*1000.0 *areaTH_s(i,t,v,h,j)
                                
      as_fCcb_c                 = as_fCcb_c + (fCO2gc_W(j) -fCcg_M(j)-fCcg_G(j)) *c_MC*c_year*1000.0 *areaTH_s(i,t,v,h,j)

      as_fCbo                   = as_fCbo + fCbo(j) *c_MC*c_year*1000.0 *areaTH_s(i,t,v,h,j)
   
      as_fH2Ogl_ut              = as_fH2Ogl_ut + min(0.0,fH2Olg_xu(j)) *1000.0*c_year *areaTH_s(i,t,v,h,j)
                                
      as_fH2Olg_tu              = as_fH2Olg_tu + max(0.0,fH2Olg_xu(j)) *1000.0*c_year *areaTH_s(i,t,v,h,j)
                                
      as_fH2Ol_bx               = as_fH2Ol_bx + fH2Ol_bx(j) *1000.0*c_year *areaTH_s(i,t,v,h,j)

        as_Ts                   = as_Ts + xT_s(j) *cweight
   
      as_Tg                     = as_Tg + xT_g(i,t,h,j) *lground *cweight
   
      as_H                      = as_H + fRAD_H(j) *areaTH_s(i,t,v,h,j)
   
      as_G                      = as_G + fQ_tg(j) *areaTH_s(i,t,v,h,j)
   
      as_E                      = as_E + fQ_ta_L(j) *areaTH_s(i,t,v,h,j)
   
      as_C                      = as_C + fQ_ta_S(j) *areaTH_s(i,t,v,h,j)
   
      as_EB                     = as_EB + (fRAD_H(j)-fQ_ta_L(j)-fQ_ta_S(j)) *areaTH_s(i,t,v,h,j)
   
      ! average BSC-related properties
   
      !     if (i2 .eq. 1 .and. BSCtypes) call libry_accBSC(i,i2,j)
   
      endif ! check for survival
   
  enddo ! End of loop over all species - VI -

endif ! write output?

return
end subroutine nova_step

end module libry_nova

