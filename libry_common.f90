
!#######################################################################
! LIBRY_COMMON
!#######################################################################

module libry_common
contains


! LIBRY_NAMELIST
!-----------------------------------------------------------------------

subroutine libry_namelist ()
use libry_par
implicit none

integer :: stat

! define namelist parameters

namelist /librypar/ year0,        &
                    cyear0,       &
                    tsindata0,    &
                    accts0,       &
                    tpos0,        &
                    lastyear,     &
                    runperiod,    &
                    tsl,          & ! [s]
                    yearout1,     &
                    yearoutX,     &
                    outint,       &
                    nSites,       &
                    p_nspec,      &
                    fracTH_s_init,&
                    lrestart,     &
                    NOHONO,       &
                    BSCtypes,     &
                    llevels,      &
                    ltraits

! read namelist

open( knamelist, file=snamelist, status='old', iostat=stat )

if ( stat .eq. 0 ) then
  read( knamelist, librypar )
  close( knamelist )
else
  write(*,*) "ERROR reading namelist, using default values"
endif

! check if tsl is le 1 h

if ( tsl .gt. 3600 ) then
  write(*,*) "ERROR: time step may not be longer than 1 hour"
  stop
endif

return
end subroutine libry_namelist


! LIBRY_ALLOC
!-----------------------------------------------------------------------

subroutine libry_alloc ()
use libry_par
implicit none

! allocate variables

allocate( naccu(nCPts), &
naccu_d(nCPts), &
naccu_n(nCPts), &
naccu_m(nCPts), &
!
vec_o(p_nspec,p_nspecpar), &
o_albedo2(p_nspec), &
o_theta_max(p_nspec), &
o_spec_area(p_nspec), &
o_sat_A(p_nspec), &
o_sat_X(p_nspec), &
o_vcmax_M(p_nspec), &
o_spec_Rubisco(p_nspec), &
o_resp_main(p_nspec), &
o_turnover(p_nspec), &
o_vomax_M(p_nspec), &
o_ToptP(p_nspec), &
o_Q10_resp(p_nspec), &
o_Eact_Kc(p_nspec), &
o_Eact_Ko(p_nspec), &
o_Eact_Vm(p_nspec), &
o_Eact_Jm(p_nspec), &
o_CCM(p_nspec), &
o_DCO2(p_nspec), &
o_DCO2M(p_nspec), &
o_redkB(p_nspec), &
o_rs(p_nspec), &
o_zt(p_nspec), &
o_dzt(p_nspec), &
o_Nspec(p_nspec), &
!o_prs(p_nspec), &
o_satHph(p_nspec), &
LAInvv(p_nspec), &
!o_extL(p_nspec), &
!o_cvxL(p_nspec), &
!
xT_a(nCPts), &
fRADs_ad(nCPts), &
fRADl_ad(nCPts), &
fH2Ol_ad(nCPts), &
fH2Os_ad(nCPts), &
rH2Og_RH(nCPts), &
fAIR_s(nCPts), &
areaLEAF(nCPts,p_stepLAI), &
areaSTEM(nCPts,p_stepLAI), &
areaLEAF_month(nCPts), &
areaSTEM_month(nCPts), &
MareaLEAF(nCPts), &
MareaSTEM(nCPts), &
fracLEAF(nCPts), &
fracSTEM(nCPts), &
areaSOIL(nCPts), &
LAIforest(nCPts), &
LAIgrass(nCPts), &
Acano(nCPts,p_ntiles), &
biome(nCPts), &
frac_noland(nCPts), &
frac_tile(nCPts,p_ntiles), &
albLsf(nCPts,p_ntiles), &
tauD(nCPts,p_ntiles,p_nvertM,p_nhabM), &
roughlen(nCPts,p_ntiles,p_nvertM,p_nhabM), &
!dn_kH2Og(nCPts,p_ntiles,p_nvertM,p_nhabM), &
rmaxH2Ol_b(nCPts,p_ntiles,p_nvertM,p_nhabM), &
rH2Os_g(nCPts), &
rH2Ol_0(nCPts,p_ntiles,p_nvertM,p_nhabM), &
rH2Ol_g1(nCPts,p_ntiles), &
rH2Ol_g2(nCPts,p_ntiles), &
xT_g0(nCPts,p_ntiles,p_nhabM), &
!
CSOIL(nCPts), &
kSOIL(nCPts), &
dewmax(nCPts), &
!
!kH2Og(nCPts,p_nhabM), &
rmaxH2Ol_t(p_nspec), &
sat(p_nspec), &
xT_s(p_nspec), &
CO2_p(p_nspec), &
jmax(p_nspec), &
vcmax(p_nspec), &
act(p_nspec), &
expansion(p_nspec), &
wgtspec(p_nspec), &
!
klife(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
rH2Ol_t(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
rH2Ol_b(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
rH2Os_t(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
act_state(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
gpp0(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
areaTH_s(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
xT_g(nCPts,p_ntiles,p_nhabM,p_nspec), &
netgrowth(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
rNd_t(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
Ndemand(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
Nsupply(nCPts,p_ntiles,p_nvertM,p_nhabM,p_nspec), &
!
fCO2gc_L(p_nspec), &
fCO2gc_W(p_nspec), &
fCO2gc(p_nspec), &
fCcg_M(p_nspec), &
fCcg_G(p_nspec), &
fCcb(p_nspec), &
fCbo(p_nspec), &
fH2Ol_xd(p_nspec), &
fH2Olg_xu(p_nspec), &
fH2Ol_bx(p_nspec), &
fH2Ol_bd(p_nspec), &
fQ_tg(p_nspec), &
fQ_ta_L(p_nspec), &
fQ_ta_S(p_nspec), &
fRAD_H(p_nspec), &
!
fH2Ol_ug2(nCPts,p_ntiles), &
!
ah_areaTH_s(p_nhabM), &
ah_rH2Ol(p_nhabM), &
ah_rH2Os(p_nhabM), &
ah_rmaxH2Ol(p_nhabM), &
ah_act(p_nhabM), &
ah_actB(p_nhabM), &
ah_fH2Ol_xd(p_nhabM), &
ah_fH2Ogl_ux(p_nhabM), &
ah_fH2Olg_xu(p_nhabM), &
ah_fH2Ol_bd(p_nhabM), &
ah_fH2Ol_bx(p_nhabM), &
ah_Ts(p_nhabM), &
ah_Tg(p_nhabM), &
ah_H(p_nhabM), &
ah_G(p_nhabM), &
ah_E(p_nhabM), &
ah_C(p_nhabM), &
ah_EB(p_nhabM), &
ah_rCO2d(p_nhabM), &
ah_rCb(p_nhabM), &
ah_fCO2gc(p_nhabM), &
ah_fCcg(p_nhabM), &
ah_fCcb(p_nhabM), &
ah_fCcb_l(p_nhabM), &
ah_fCcb_c(p_nhabM), &
ah_fCbo(p_nhabM), &
ah_frgr(p_nhabM), &
!
av_areaTH_s(p_nvertM), &
av_rH2Ol(p_nvertM), &
av_rH2Os(p_nvertM), &
av_rmaxH2Ol(p_nvertM), &
av_act(p_nvertM), &
av_actB(p_nvertM), &
av_fH2Ol_ux(p_nvertM), &
av_fH2Ol_xd(p_nvertM), &
av_fH2Ogl_ux(p_nvertM), &
av_fH2Olg_xu(p_nvertM), &
av_fH2Ol_bx(p_nvertM), &
av_fH2Ol_bd(p_nvertM), &
av_Ts(p_nvertM), &
av_H(p_nvertM), &
av_E(p_nvertM), &
av_C(p_nvertM), &
av_EB(p_nvertM), &
av_rCO2d(p_nvertM), &
av_rCb(p_nvertM), &
av_fCO2gc(p_nvertM), &
av_fCcg(p_nvertM), &
av_fCcb(p_nvertM), &
av_fCcb_l(p_nvertM), &
av_fCcb_c(p_nvertM), &
av_fCbo(p_nvertM), &
av_frgr(p_nvertM), &
!
at_fH2Ol_ux(p_ntiles,2), &
at_fH2Ol_xd(p_ntiles,2), &
at_areaTH_s(p_ntiles,2), &
at_rH2Ol(p_ntiles,2), &
at_rH2Os(p_ntiles,2), &
at_rmaxH2Ol(p_ntiles,2), &
at_act(p_ntiles,2), &
at_actB(p_ntiles,2), &
at_rCO2d(p_ntiles,2), &
at_rCb(p_ntiles,2), &
at_fCO2gc(p_ntiles,2), &
at_fCcg(p_ntiles,2), &
at_fCcb(p_ntiles,2), &
at_fCcb_l(p_ntiles,2), &
at_fCcb_c(p_ntiles,2), &
at_fCbo(p_ntiles,2), &
at_frgr(p_ntiles,2), &
at_fH2Ogl_ux(p_ntiles,2), &
at_fH2Olg_xu(p_ntiles,2), &
at_fH2Ol_bx(p_ntiles,2), &
at_fH2Ol_bd(p_ntiles,2), &
at_fH2Ol_ug2(p_ntiles), &
at_Ts(p_ntiles,2), &
at_H(p_ntiles,2), &
at_E(p_ntiles,2), &
at_C(p_ntiles,2), &
at_EB(p_ntiles,2), &
at_Tg(p_ntiles), &
at_G(p_ntiles), &
at_rH2Ol_g1(p_ntiles), &
at_rH2Ol_g2(p_ntiles), &
at_fH2Ol_ug(p_ntiles), &
at_fH2Ol_go(p_ntiles), &
at_fH2Ol_gb(p_ntiles), &
at_fH2Olg_ga(p_ntiles), &
!
atN_fH2Ol_xd(nCPts,p_ntiles,2), &
atN_fH2Ol_bd(nCPts,p_ntiles,2), &
!
ag_fH2Ol_ux(nCPts,2), &
ag_fH2Ol_xd(nCPts,2), &
ag_areaTH_s(nCPts,2), &
ag_rH2Ol(nCPts,2), &
ag_rH2Os(nCPts,2), &
ag_rmaxH2Ol(nCPts,2), &
ag_act(nCPts,2), &
ag_actB(nCPts,2), &
ag_rCO2d(nCPts,2), &
ag_rCb(nCPts,2), &
ag_fCO2gc(nCPts,2), &
ag_fCcg(nCPts,2), &
ag_fCcb(nCPts,2), &
ag_fCcb_l(nCPts,2), &
ag_fCcb_c(nCPts,2), &
ag_fCbo(nCPts,2), &
ag_frgr(nCPts,2), &
ag_fH2Ogl_ux(nCPts,2), &
ag_fH2Olg_xu(nCPts,2), &
ag_fH2Ol_bx(nCPts,2), &
ag_fH2Ol_bd(nCPts,2), &
ag_fH2Ol_ug2(nCPts), &
ag_Ts(nCPts,2), &
ag_H(nCPts,2), &
ag_E(nCPts,2), &
ag_C(nCPts,2), &
ag_EB(nCPts,2), &
ag_Tg(nCPts), &
ag_G(nCPts), &
ag_rH2Ol_g1(nCPts), &
ag_rH2Ol_g2(nCPts), &
ag_fH2Ol_ug(nCPts), &
ag_fH2Ol_go(nCPts), &
ag_fH2Ol_gb(nCPts), &
ag_fH2Olg_ga(nCPts), &
ag_rH2Os_g(nCPts), &
ag_fH2Osl_g(nCPts), &
ag_fH2Os_ad(nCPts), &
ag_xT_a(nCPts), &
ag_fH2Ol_ad(nCPts), &
ag_fRADs(nCPts), &
!
count_spec(nCPts), &
count_spec_h(p_nhabA,nCPts), &
vec_o_avg(6,p_nspecpar,p_nhabA,nCPts) )

return
end subroutine libry_alloc


! LIBRY_DEALLOC
!-----------------------------------------------------------------------

subroutine libry_dealloc ()
use libry_par
implicit none

! deallocate variables

deallocate( naccu, &
naccu_d, &
naccu_n, &
naccu_m, &
!
vec_o, &
o_albedo2, &
o_theta_max, &
o_spec_area, &
o_sat_A, &
o_sat_X, &
o_vcmax_M, &
o_spec_Rubisco, &
o_resp_main, &
o_turnover, &
o_vomax_M, &
o_ToptP, &
o_Q10_resp, &
o_Eact_Kc, &
o_Eact_Ko, &
o_Eact_Vm, &
o_Eact_Jm, &
o_CCM, &
o_DCO2, &
o_DCO2M, &
o_redkB, &
o_rs, &
o_zt, &
o_dzt, &
o_Nspec, &
!o_prs, &
o_satHph, &
LAInvv, &
!o_extL, &
!o_cvxL, &
!
xT_a, &
fRADs_ad, &
fRADl_ad, &
fH2Ol_ad, &
fH2Os_ad, &
rH2Og_RH, &
fAIR_s, &
areaLEAF, &
areaSTEM, &
areaLEAF_month, &
areaSTEM_month, &
MareaLEAF, &
MareaSTEM, &
fracLEAF, &
fracSTEM, &
areaSOIL, &
LAIforest, &
LAIgrass, &
Acano, &
biome, &
frac_noland, &
frac_tile, &
albLsf, &
tauD, &
roughlen, &
!dn_kH2Og, &
rmaxH2Ol_b, &
rH2Os_g, &
rH2Ol_0, &
rH2Ol_g1, &
rH2Ol_g2, &
xT_g0, &
!
CSOIL, &
kSOIL, &
dewmax, &
!
rmaxH2Ol_t, &
sat, &
xT_s, &
CO2_p, &
jmax, &
vcmax, &
act, &
netgrowth, &
rNd_t, &
Ndemand, &
Nsupply, &
expansion, &
wgtspec, &
!
klife, &
rH2Ol_t, &
rH2Ol_b, &
rH2Os_t, &
act_state, &
gpp0, &
areaTH_s, &
xT_g, &
!
fCO2gc_L, &
fCO2gc_W, &
fCO2gc, &
fCcg_M, &
fCcg_G, &
fCcb, &
fCbo, &
fH2Ol_xd, &
fH2Olg_xu, &
fH2Ol_bx, &
fH2Ol_bd, &
fQ_tg, &
fQ_ta_L, &
fQ_ta_S, &
fRAD_H, &
!
fH2Ol_ug2, &
!
ah_areaTH_s, &
ah_rH2Ol, &
ah_rH2Os, &
ah_rmaxH2Ol, &
ah_act, &
ah_actB, &
ah_fH2Ol_xd, &
ah_fH2Ogl_ux, &
ah_fH2Olg_xu, &
ah_fH2Ol_bd, &
ah_fH2Ol_bx, &
ah_Ts, &
ah_Tg, &
ah_H, &
ah_G, &
ah_E, &
ah_C, &
ah_EB, &
ah_rCO2d, &
ah_rCb, &
ah_fCO2gc, &
ah_fCcg, &
ah_fCcb, &
ah_fCcb_l, &
ah_fCcb_c, &
ah_fCbo, &
ah_frgr, &
!
av_areaTH_s, &
av_rH2Ol, &
av_rH2Os, &
av_rmaxH2Ol, &
av_act, &
av_actB, &
av_fH2Ol_ux, &
av_fH2Ol_xd, &
av_fH2Ogl_ux, &
av_fH2Olg_xu, &
av_fH2Ol_bx, &
av_fH2Ol_bd, &
av_Ts, &
av_H, &
av_E, &
av_C, &
av_EB, &
av_rCO2d, &
av_rCb, &
av_fCO2gc, &
av_fCcg, &
av_fCcb, &
av_fCcb_l, &
av_fCcb_c, &
av_fCbo, &
av_frgr, &
!
at_fH2Ol_ux, &
at_fH2Ol_xd, &
at_areaTH_s, &
at_rH2Ol, &
at_rH2Os, &
at_rmaxH2Ol, &
at_act, &
at_actB, &
at_rCO2d, &
at_rCb, &
at_fCO2gc, &
at_fCcg, &
at_fCcb, &
at_fCcb_l, &
at_fCcb_c, &
at_fCbo, &
at_frgr, &
at_fH2Ogl_ux, &
at_fH2Olg_xu, &
at_fH2Ol_bx, &
at_fH2Ol_bd, &
at_Ts, &
at_H, &
at_E, &
at_C, &
at_EB, &
at_Tg, &
at_G, &
at_rH2Ol_g1, &
at_rH2Ol_g2, &
at_fH2Ol_ug, &
at_fH2Ol_go, &
at_fH2Ol_gb, &
at_fH2Olg_ga, &
at_fH2Ol_ug2, &
!
atN_fH2Ol_xd, &
atN_fH2Ol_bd, &
!
ag_fH2Ol_ux, &
ag_fH2Ol_xd, &
ag_areaTH_s, &
ag_rH2Ol, &
ag_rH2Os, &
ag_rmaxH2Ol, &
ag_act, &
ag_actB, &
ag_rCO2d, &
ag_rCb, &
ag_fCO2gc, &
ag_fCcg, &
ag_fCcb, &
ag_fCcb_l, &
ag_fCcb_c, &
ag_fCbo, &
ag_frgr, &
ag_fH2Ogl_ux, &
ag_fH2Olg_xu, &
ag_fH2Ol_bx, &
ag_fH2Ol_bd, &
ag_Ts, &
ag_H, &
ag_E, &
ag_C, &
ag_EB, &
ag_Tg, &
ag_G, &
ag_rH2Ol_g1, &
ag_rH2Ol_g2, &
ag_fH2Ol_ug, &
ag_fH2Ol_go, &
ag_fH2Ol_gb, &
ag_fH2Olg_ga, &
ag_fH2Ol_ug2, &
ag_rH2Os_g, &
ag_fH2Osl_g, &
ag_fH2Os_ad, &
ag_xT_a, &
ag_fH2Ol_ad, &
ag_fRADs, &
!
count_spec, &
count_spec_h, &
vec_o_avg )

return
end subroutine libry_dealloc


! LIBRY_READSPEC -- READ STRATEGY PARAMETERS
!-----------------------------------------------------------------------

subroutine libry_readSpec ()
use libry_par
implicit none

integer :: stat,j

character(9) :: fd0

write( fd0 ,'(I9)' ) p_nspecpar

! read species parameter file with random numbers [0, 1]

open( kspecpar, file=sspecpar, status='old', action='read', iostat=stat )

if ( stat .ne. 0 ) then
  write(*,*) "ERROR opening species parameter file"
  stop
endif

do j = 1,p_nspec
  read( kspecpar, '('//fd0//'F6.3)' ) vec_o(j,:)
enddo

close( kspecpar )

return
end subroutine libry_readSpec


! LIBRY_READMON -- READ MONTHLY INPUT
!-----------------------------------------------------------------------

subroutine libry_readMon(m)
use libry_par
implicit none

integer :: m

! assign monthly canopy properties
areaLEAF_month(:) = areaLEAF(:,m)
areaSTEM_month(:) = areaSTEM(:,m)

return
end subroutine libry_readMon


! LIBRY_RESET
!-----------------------------------------------------------------------

subroutine libry_reset (nCPts2)
use libry_par
implicit none

integer :: i,c
integer :: nCPts2
real :: fillv

naccu(:) = 0

naccu_d(:) = 0
naccu_n(:) = 0

do i = 1,nCPts2
  do c = 1,2

    if (c .eq. 2 .and. frac_tile(i,1) .le. p_critD) then
      fillv = -9999.0
    else
      fillv = 0.0
    endif

    ag_fH2Ol_ux(i,c)            = fillv
    ag_fH2Ol_xd(i,c)            = fillv
   
    ag_areaTH_s(i,c)            = fillv
    ag_rH2Ol(i,c)               = fillv
    ag_rH2Os(i,c)               = fillv
    ag_rmaxH2Ol(i,c)            = fillv
    ag_act(i,c)                 = fillv
    ag_actB(i,c)                = fillv
    ag_rCO2d(i,c)               = fillv
    ag_rCb(i,c)                 = fillv
    ag_fCO2gc(i,c)              = fillv
    ag_fCcg(i,c)                = fillv
    ag_fCcb(i,c)                = fillv
    ag_fCcb_l(i,c)              = fillv
    ag_fCcb_c(i,c)              = fillv
    ag_fCbo(i,c)                = fillv
    ag_frgr(i,c)                = fillv
   
    ag_fH2Ogl_ux(i,c)           = fillv
    ag_fH2Olg_xu(i,c)           = fillv
    ag_fH2Ol_bx(i,c)            = fillv
   
    ag_fH2Ol_bd(i,c)            = fillv
   
    ag_Ts(i,c)                  = fillv
    ag_H(i,c)                   = fillv
    ag_E(i,c)                   = fillv
    ag_C(i,c)                   = fillv
    ag_EB(i,c)                  = fillv
  enddo
enddo

ag_fH2Ol_ug2(:)                 = 0.0

ag_Tg(:)                        = 0.0
ag_G(:)                         = 0.0
                      
ag_rH2Ol_g1(:)                  = 0.0
ag_rH2Ol_g2(:)                  = 0.0
                      
ag_fH2Ol_ug(:)                  = 0.0
ag_fH2Ol_go(:)                  = 0.0
ag_fH2Ol_gb(:)                  = 0.0
ag_fH2Olg_ga(:)                 = 0.0
                      
ag_rH2Os_g(:)                   = 0.0
ag_fH2Osl_g(:)                  = 0.0
                      
ag_fH2Os_ad(:)                  = 0.0
ag_xT_a(:)                      = 0.0
ag_fH2Ol_ad(:)                  = 0.0
ag_fRADs(:)                     = 0.0

return
end subroutine libry_reset

! LIBRY_AV_OUTPUT
!-----------------------------------------------------------------------

subroutine libry_av_output(nCPts2)
use libry_par
! use libry_opt
implicit none

integer :: i
integer :: nCPts2

do i = 1,nCPts2

  ag_fH2Ol_ux(i,:)             = ag_fH2Ol_ux(i,:) / real(naccu(i))
  ag_fH2Ol_xd(i,:)             = ag_fH2Ol_xd(i,:) / real(naccu(i))

  ag_areaTH_s(i,:)             = ag_areaTH_s(i,:) / real(naccu(i))
  ag_rH2Ol(i,:)                = ag_rH2Ol(i,:)    / real(naccu(i))
  ag_rH2Os(i,:)                = ag_rH2Os(i,:)    / real(naccu(i))
  ag_rmaxH2Ol(i,:)             = ag_rmaxH2Ol(i,:) / real(naccu(i))
  ag_act(i,:)                  = ag_act(i,:)      / real(naccu(i))
  ag_actB(i,:)                  = ag_actB(i,:)     / real(naccu(i))
  ag_rCO2d(i,:)                = ag_rCO2d(i,:)    / real(naccu(i))
  ag_rCb(i,:)                  = ag_rCb(i,:)      / real(naccu(i)) 
  ag_fCO2gc(i,:)               = ag_fCO2gc(i,:)   / real(naccu(i))
  ag_fCcg(i,:)                 = ag_fCcg(i,:)     / real(naccu(i))
  ag_fCcb(i,:)                 = ag_fCcb(i,:)     / real(naccu(i))
  ag_fCcb_l(i,:)                = ag_fCcb_l(i,:)   / real(naccu(i))
  ag_fCcb_c(i,:)                = ag_fCcb_c(i,:)   / real(naccu(i))
  ag_fCbo(i,:)                 = ag_fCbo(i,:)     / real(naccu(i))
  ag_frgr(i,:)                  = ag_frgr(i,:)     / real(naccu(i))

  ag_fH2Ogl_ux(i,:)            = ag_fH2Ogl_ux(i,:)/ real(naccu(i))
  ag_fH2Olg_xu(i,:)            = ag_fH2Olg_xu(i,:)/ real(naccu(i))
  ag_fH2Ol_bx(i,:)              = ag_fH2Ol_bx(i,:) / real(naccu(i))
  ag_fH2Ol_bd(i,:)              = ag_fH2Ol_bd(i,:) / real(naccu(i))

  ag_fH2Ol_ug2(i)               = ag_fH2Ol_ug2(i) / real(naccu(i))

  ag_Ts(i,:)                    = ag_Ts(i,:)       / real(naccu(i))      ! _d
  ag_H(i,:)                    = ag_H(i,:)        / real(naccu(i)) 
  ag_E(i,:)                    = ag_E(i,:)        / real(naccu(i)) 
  ag_C(i,:)                    = ag_C(i,:)        / real(naccu(i)) 
  ag_EB(i,:)                   = ag_EB(i,:)       / real(naccu(i)) 

  ag_Tg(i)                    = ag_Tg(i)        / real(naccu(i)) 
  ag_G(i)                     = ag_G(i)         / real(naccu(i)) 

  ag_rH2Ol_g1(i)              = ag_rH2Ol_g1(i)  / real(naccu(i))
  ag_rH2Ol_g2(i)              = ag_rH2Ol_g2(i)  / real(naccu(i))

  ag_fH2Ol_ug(i)              = ag_fH2Ol_ug(i)  / real(naccu(i))
  ag_fH2Ol_go(i)              = ag_fH2Ol_go(i)  / real(naccu(i))
  ag_fH2Ol_gb(i)              = ag_fH2Ol_gb(i)  / real(naccu(i))
  ag_fH2Olg_ga(i)             = ag_fH2Olg_ga(i) / real(naccu(i))

  ag_rH2Os_g(i)               = ag_rH2Os_g(i)   / real(naccu(i)) 
  ag_fH2Osl_g(i)              = ag_fH2Osl_g(i)  / real(naccu(i)) 

  ag_fH2Os_ad(i)              = ag_fH2Os_ad(i)  / real(naccu(i)) 
  ag_xT_a(i)                  = ag_xT_a(i)      / real(naccu(i)) 
  ag_fH2Ol_ad(i)              = ag_fH2Ol_ad(i)  / real(naccu(i)) 
  ag_fRADs(i)                 = ag_fRADs(i)     / real(naccu(i)) 
enddo

!  if (BSCtypes) call libry_avBSC(i)
!  if (BSCtypes) call libry_avgcBSC(i,i2,j)

return
end subroutine libry_av_output

!***********************************************************************
! LIBRY_AV_SPECIES
!***********************************************************************

subroutine libry_av_species(nCPts2)
use libry_par
implicit none

integer :: i,j,k,l
integer :: nCPts2, lMd,lQ1,lQ3,lMin,lMax
real    :: Asum

! distribution of species parameters

vec_o_avg(:,:,:,:) = 0.0
count_spec_h(:,:) = 0.0
count_spec(:) = 0.0

do i = 1,nCPts2

  count_spec(i) = real(sum(klife(i,:,:,:,:)))

  do k = 1,p_nhabA ! total number of tile/level/habitat combinations

    select case( k )
      case( 1 ) ! forest floor
        Asum                    = sum(areaTH_s(i,1,1,1,:))
        count_spec_h(k,i)       = real(sum(klife(i,1,1,1,:)))

        if (Asum .gt. p_critD .and. count_spec_h(k,i) .gt. 0.1) then
          do l = 1,p_nspecpar
            vec_o_avg(1,l,k,i)  = sum( vec_o(:,l) * areaTH_s(i,1,1,1,:) / Asum )

            call sortvec(vec_o(:,l), lMd,lQ1,lQ3, lMin,lMax, areaTH_s(i,1,1,1,:),klife(i,1,1,1,:),Asum ) 

            vec_o_avg(2,l,k,i)  = vec_o(lMd,l) 
            vec_o_avg(3,l,k,i)  = vec_o(lQ1,l) 
            vec_o_avg(4,l,k,i)  = vec_o(lQ3,l) 
            vec_o_avg(5,l,k,i)  = vec_o(lMin,l) 
            vec_o_avg(6,l,k,i)  = vec_o(lMax,l) 
          end do
        else
          vec_o_avg(:,:,k,i)    = 0.5
        endif
      case( 2 ) ! forest canopy stems
        Asum                    = sum(areaTH_s(i,1,2,1,:))
        count_spec_h(k,i)       = real(sum(klife(i,1,2,1,:)))

        if (Asum .gt. p_critD .and. count_spec_h(k,i) .gt. 0.1) then
          do l = 1,p_nspecpar
            vec_o_avg(1,l,k,i)  = sum( vec_o(:,l) * areaTH_s(i,1,2,1,:) / Asum )

            call sortvec(vec_o(:,l), lMd,lQ1,lQ3, lMin,lMax, areaTH_s(i,1,2,1,:),klife(i,1,2,1,:),Asum ) 

            vec_o_avg(2,l,k,i)  = vec_o(lMd,l) 
            vec_o_avg(3,l,k,i)  = vec_o(lQ1,l) 
            vec_o_avg(4,l,k,i)  = vec_o(lQ3,l) 
            vec_o_avg(5,l,k,i)  = vec_o(lMin,l) 
            vec_o_avg(6,l,k,i)  = vec_o(lMax,l) 
          end do
        else
          vec_o_avg(:,:,k,i)    = 0.5
        endif
      case( 3 ) ! forest canopy leaves
        Asum                    = sum(areaTH_s(i,1,2,2,:))
        count_spec_h(k,i)       = real(sum(klife(i,1,2,2,:)))

        if (Asum .gt. p_critD .and. count_spec_h(k,i) .gt. 0.1) then
          do l = 1,p_nspecpar
            vec_o_avg(1,l,k,i)  = sum( vec_o(:,l) * areaTH_s(i,1,2,2,:) / Asum )

            call sortvec(vec_o(:,l), lMd,lQ1,lQ3, lMin,lMax, areaTH_s(i,1,2,2,:),klife(i,1,2,2,:),Asum ) 

            vec_o_avg(2,l,k,i)  = vec_o(lMd,l) 
            vec_o_avg(3,l,k,i)  = vec_o(lQ1,l) 
            vec_o_avg(4,l,k,i)  = vec_o(lQ3,l) 
            vec_o_avg(5,l,k,i)  = vec_o(lMin,l) 
            vec_o_avg(6,l,k,i)  = vec_o(lMax,l) 
          end do
        else
          vec_o_avg(:,:,k,i)    = 0.5
        endif
      case( 4 ) ! bare / grass
        Asum                    = sum(areaTH_s(i,2,1,1,:))
        count_spec_h(k,i)       = real(sum(klife(i,2,1,1,:)))

        if (Asum .gt. p_critD .and. count_spec_h(k,i) .gt. 0.1) then
          do l = 1,p_nspecpar
            vec_o_avg(1,l,k,i)  = sum( vec_o(:,l) * areaTH_s(i,2,1,1,:) / Asum )

            call sortvec(vec_o(:,l), lMd,lQ1,lQ3, lMin,lMax, areaTH_s(i,2,1,1,:),klife(i,2,1,1,:),Asum ) 

            vec_o_avg(2,l,k,i)  = vec_o(lMd,l) 
            vec_o_avg(3,l,k,i)  = vec_o(lQ1,l) 
            vec_o_avg(4,l,k,i)  = vec_o(lQ3,l) 
            vec_o_avg(5,l,k,i)  = vec_o(lMin,l) 
            vec_o_avg(6,l,k,i)  = vec_o(lMax,l) 
          end do
        else
          vec_o_avg(:,:,k,i)    = 0.5
        endif
      case default
        Asum                    = 0.0
        count_spec_h(k,i)       = 0.0
        vec_o_avg(:,:,k,i)      = 0.5
    end select
  end do
enddo ! loop over all grid cells

return
end subroutine libry_av_species

!***********************************************************************
! LIBRY_AV_SPECIES
!***********************************************************************

subroutine sortvec(x,lm,l1,l3,lmin,lmax,A,K,Az)
use libry_par
implicit none

integer, dimension(p_nspec), intent(in) :: K
real, dimension(p_nspec), intent(in)    :: x,A
integer, allocatable                    :: snb(:)
real, allocatable                       :: y(:)

integer, intent(out)    :: lm,l1,l3, lmin, lmax
integer                 :: zm,z1,z3
integer                 :: i,s2,j,sK,k1,k2

real, intent(in)        :: Az
real                    :: y2,accA

logical                 :: order

sK = sum(K)

allocate(y(sK), snB(sK))

k2 = 0

do k1=1,p_nspec

  if (K(k1) .gt. 0) then

    k2 = k2 + 1

    snb(k2)     = k1
    y(k2)       = x(k1)
  endif
enddo

!y   = pack(x, K .ne. 0)
!snb = pack(sn, K .ne. 0)

!if (size(y) .gt. 1) then
if (k2 .gt. 1) then

  order = .false.

  do while (.not. order)

    do i=1,k2-1

      if (y(i) .gt. y(i+1)) then

        y2 = y(i)
        y(i) = y(i+1)
        y(i+1) = y2

        s2 = snb(i)
        snb(i) = snb(i+1)
        snb(i+1) = s2

        order = .false.
      else 
        if(i .eq. 1) order = .true.
      endif
    enddo
  enddo

  accA=0.0
  zm = 0
  z1 = 0
  z3 = 0

  do j=1,k2

    accA = accA + A(snb(j))/Az

    if (accA .ge. 0.25 .and. z1 .eq. 0) then 
      l1 = snb(j)
      z1 = 1
    endif
    if (accA .ge. 0.5 .and. zm .eq. 0) then 
      lm = snb(j)
      zm = 1
    endif
    if (accA .ge. 0.75 .and. z3 .eq. 0) then 
      l3 = snb(j)
      z3 = 1
    endif
  enddo

  lmin = snb(1)
  lmax = snb(k2)
else

  l1 = snb(1)
  lm = snb(1)
  l3 = snb(1)
  lmin = snb(1)
  lmax = snb(1)
endif

deallocate(y, snB)

return
end subroutine sortvec

end module libry_common

