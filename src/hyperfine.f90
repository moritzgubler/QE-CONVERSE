!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Ported to QE-CONVERSE / QE 7.5 from qe-gipaw-6.2 (Knight-shift fork).
!
! Computes the isotropic Fermi-contact hyperfine coupling, i.e. the
! electron spin density at the nuclei, which is the microscopic
! ingredient of the metallic Knight shift (PRB 76, 165122 (2007);
! J. Phys. Chem. C 120, 27059 (2016)).  The contact spin density is the
! sum of three terms:
!    bare       : smooth (pseudo) spin density at the nucleus
!    GIPAW      : all-electron reconstruction near the nucleus
!    core-relax : spin polarization of the core s shells (core_relax.f90)
!
! Only the Fermi-contact (isotropic) part is implemented; the anisotropic
! dipolar hyperfine tensor is not computed here.
!
! Spin convention: channel 1 = up, channel 2 = down; the contact density
! is the magnetization n_up - n_dn at the nucleus.
!
!-----------------------------------------------------------------------
SUBROUTINE hyperfine
  !-----------------------------------------------------------------------
  USE kinds,                  ONLY : dp
  USE io_global,              ONLY : stdout
  USE parameters,             ONLY : ntypx
  USE constants,              ONLY : pi, fpi, c_si, bohr_radius_si, rytoev, electronvolt_si
  USE fft_base,               ONLY : dffts
  USE ions_base,              ONLY : nat, atm, ityp
  USE mp,                     ONLY : mp_sum
  USE lsda_mod,               ONLY : nspin
  USE klist,                  ONLY : two_fermi_energies
  USE ener,                   ONLY : ef_up, ef_dw
  USE gipaw_module,           ONLY : hfi_nuclear_g_factor, hfi_output_unit, &
                                     core_relax_method, g_e

  !-- constants ----------------------------------------------------------
  IMPLICIT NONE
  real(dp), parameter :: mu0_by_fpi = 1e-7
  real(dp), parameter :: mu_n = 5.05078324e-27_dp
  real(dp), parameter :: bohr_radius = bohr_radius_si
  real(dp), parameter :: gamma_e = 28024953.64_dp
  real(dp), parameter :: lambda = C_SI / 1.0e+8_dp
  real(dp), parameter :: common_factor = mu0_by_fpi * mu_n / Bohr_radius ** 3

  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: spin_den(:), rho_s(:,:)
  real(dp), allocatable :: hfi_fc_bare(:)
  real(dp), allocatable :: hfi_fc_gipaw(:)
  real(dp), allocatable :: hfi_fc_core(:), hfi_fc_tot(:)
  integer :: na
  real(dp) :: output_factor, fact
  ! Knight shift (spin shielding) sigma_s
  real(dp), parameter :: mu_b_si = 9.2740100783e-24_dp   ! Bohr magneton, J/T
  real(dp) :: delta_ef, mu0, ry_to_j, sigma_factor, b_ext, sigma_s

  call start_clock('hyperfine')

  if (nspin /= 2) call errore('hyperfine', 'Fermi-contact requires a spin-polarized (nspin=2) calculation', 1)

  allocate( hfi_fc_bare(nat) )
  allocate( hfi_fc_gipaw(nat) )
  allocate( hfi_fc_core(nat), hfi_fc_tot(nat) )

  !--------------------------------------------------------------------
  ! Fermi-contact contribution to hyperfine
  !--------------------------------------------------------------------
  ! bare contribution: smooth spin density at the nucleus
  allocate( rho_s(dffts%nnr,2), spin_den(dffts%nnr) )
  call get_smooth_density(rho_s)        ! per-spin (up,down) on the smooth grid
  spin_den(:) = rho_s(:,1) - rho_s(:,2)
  call hfi_fc_bare_el(spin_den, hfi_fc_bare)
  deallocate( rho_s, spin_den )

  ! GIPAW reconstruction contribution
  call hfi_fc_gipaw_correction(hfi_fc_gipaw)

  ! core-relaxation contribution
  call hfi_fc_core_relax(core_relax_method, hfi_fc_core)

  !--------------------------------------------------------------------
  ! Print results
  !--------------------------------------------------------------------
  select case (trim(hfi_output_unit))
      case('MHz')
          output_factor = 1d-3 * common_factor * gamma_e
      case('mT')
          output_factor = 1d3  * common_factor
      case('G', 'Gauss')
          output_factor = 1d4  * common_factor
      case('10e-4cm^-1')
          output_factor = 1d-3 * common_factor * gamma_e / lambda
      case default
          call errore( "hyperfine", "unknown units for output", 1)
  end select

  write(stdout,*)
  write(stdout,'(5X,''KNIGHT SHIFT / FERMI-CONTACT HYPERFINE COUPLINGS'')')
  write(stdout,'(5X,''USING CORE-RELAXATION METHOD: PRB 76, 035124 (2007)'')')
  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  write(stdout,*)

  write(stdout,'(5X,''NUCLEAR G-TENSORS FROM INPUT:'')')
  do na = 1, nat
      write(stdout,'(5X,A,I3,2X,F12.6)') atm(ityp(na)), na, hfi_nuclear_g_factor(ityp(na))
  enddo
  write(stdout,*)

  write(stdout,'(5X,''----- contact spin-densities in bohrradius^-3 -----'')')
  write(stdout,'(5X,8X,''  bare            GIPAW           core-relax      total'')')
  do na = 1, nat
      hfi_fc_tot(na) = hfi_fc_bare(na) + hfi_fc_gipaw(na) + hfi_fc_core(na)
      write(stdout,1002) atm(ityp(na)), na, hfi_fc_bare(na), &
          hfi_fc_gipaw(na), hfi_fc_core(na), hfi_fc_tot(na)
  enddo
  write(stdout,*)

  write(stdout,'(5X,''----- Fermi contact in '',A,'' -----'')') trim(hfi_output_unit)
  write(stdout,'(5X,8X,''  bare            GIPAW           core-relax      total'')')
  do na = 1, nat
      hfi_fc_tot(na) = hfi_fc_bare(na) + hfi_fc_gipaw(na) + hfi_fc_core(na)
      fact = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
      write(stdout,1002) atm(ityp(na)), na, fact * hfi_fc_bare(na), &
          fact * hfi_fc_gipaw(na), fact * hfi_fc_core(na), fact * hfi_fc_tot(na)
  enddo
1002 FORMAT(5X,A,I3,2X,4(F14.6,2X))

  !--------------------------------------------------------------------
  ! Knight shift as a spin shielding sigma_s (ppm), addable to the
  ! orbital (converse) shielding sigma_o.  Following Ferreira, Reuter,
  ! Scheurer, J. Phys. Chem. C 120, 25530 (2016), eqs. (1) and (5):
  !     B_FC(R) = (2/3) mu0 g_e mu_B rho_s(R)            [contact field]
  !     B_ext   = (ef_up - ef_dw) / (g_e mu_B)           [effective field]
  !     sigma_s = - B_FC(R) / B_ext      (shielding convention, eq.1)
  ! so that  sigma_tot = sigma_o + sigma_s.
  ! This needs a fixed-moment metallic SCF (two Fermi energies): the
  ! self-consistent spin splitting plays the role of the (inverse) spin
  ! susceptibility, so no external chi is required.
  !--------------------------------------------------------------------
  write(stdout,*)
  if (.not. two_fermi_energies) then
     write(stdout,'(5X,''Knight shift sigma_s (ppm) not computed:'')')
     write(stdout,'(5X,''  it needs a fixed-moment metallic SCF (two Fermi energies):'')')
     write(stdout,'(5X,''  occupations=''''smearing'''' together with a (small) tot_magnetization.'')')
  else
     delta_ef = ef_up - ef_dw                       ! Ry
     if (abs(delta_ef) < 1.d-12) then
        write(stdout,'(5X,''Knight shift sigma_s (ppm) not computed: spin splitting ~ 0'')')
     else
        mu0     = fpi * 1.0e-7_dp                   ! vacuum permeability (SI)
        ry_to_j = rytoev * electronvolt_si          ! Ry -> Joule
        ! sigma_s = -(2 mu0/3)(g_e mu_B)^2 rho_s[m^-3] / (delta_ef[J])
        sigma_factor = (2.d0*mu0/3.d0) * (g_e*mu_b_si)**2 / (bohr_radius_si**3 * ry_to_j)
        b_ext = delta_ef * ry_to_j / (g_e * mu_b_si)               ! Tesla
        write(stdout,'(5X,''----- Knight shift (spin shielding) sigma_s -----'')')
        write(stdout,'(5X,''spin splitting ef_up-ef_dw = '',F12.6,'' Ry'',4X, &
             &''effective B_ext = '',ES12.4,'' T'')') delta_ef, b_ext
        write(stdout,'(5X,''(add sigma_s to the orbital/converse shielding sigma_o)'')')
        write(stdout,'(5X,8X,''   sigma_s (ppm)'')')
        do na = 1, nat
           sigma_s = -1.d6 * sigma_factor * hfi_fc_tot(na) / delta_ef
           write(stdout,'(5X,A,I3,2X,F16.4)') atm(ityp(na)), na, sigma_s
        enddo
     endif
  endif

  deallocate( hfi_fc_bare, hfi_fc_gipaw, hfi_fc_core, hfi_fc_tot )

  call stop_clock('hyperfine')

END SUBROUTINE hyperfine



!-----------------------------------------------------------------------
SUBROUTINE get_smooth_density(rho)
  !-----------------------------------------------------------------------
  !
  ! ... Get the (per-spin) charge density on the smooth grid
  !
  USE kinds,                  ONLY : dp
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : inter_pool_comm
  USE lsda_mod,               ONLY : current_spin, isk, nspin
  USE wvfct,                  ONLY : nbnd, wg, g2kin, current_k
  USE gvecw,                  ONLY : gcutw
  USE klist,                  ONLY : nks, xk, igk_k, ngk
  USE gvect,                  ONLY : ngm, g
  USE wavefunctions,          ONLY : evc
  USE cell_base,              ONLY : tpiba2, omega
  USE io_files,               ONLY : nwordwfc, iunwfc
  USE buffers,                ONLY : get_buffer
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : invfft
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: rho(dffts%nnr,nspin)
  !-- local variables ----------------------------------------------------
  complex(dp), allocatable :: psic(:)
  integer :: ibnd, ik
  integer :: npw

  allocate( psic(dffts%nnr) )
  rho = 0.d0

  ! loop over k-points
  do ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     npw = ngk(ik)

     ! initialize at k-point k and read wfcs from (in-memory) buffer
     call gk_sort(xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
     call get_buffer (evc, nwordwfc, iunwfc, ik)

     ! loop over bands
     do ibnd = 1, nbnd
       psic(:) = (0.d0,0.d0)
       psic(dffts%nl(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
       call invfft ('Wave', psic, dffts)
       rho(:,current_spin) = rho(:,current_spin) + wg(ibnd,ik) * &
                             (dble(psic(:))**2 + aimag(psic(:))**2) / omega
     enddo
  enddo
#ifdef __MPI
  ! reduce over k-points
  call mp_sum( rho, inter_pool_comm )
#endif

  deallocate( psic )

END SUBROUTINE get_smooth_density



!-----------------------------------------------------------------------
SUBROUTINE hfi_fc_bare_el(rho_s, hfi_bare)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the bare contribution to the Fermi-contact hyperfine
  ! ... (smooth spin density evaluated at the nuclear positions)
  !
  USE kinds,                  ONLY : dp
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : intra_pool_comm
  USE constants,              ONLY : tpi, fpi
  USE gvecs,                  ONLY : ngms
  USE gvect,                  ONLY : g, gstart
  USE ions_base,              ONLY : nat, tau
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : fwfft
  USE gipaw_module,           ONLY : hfi_via_reconstruction_only

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(in) :: rho_s(dffts%nnr)
  real(dp), intent(out) :: hfi_bare(nat)
  !-- local variables ----------------------------------------------------
  complex(dp), allocatable :: rhoaux(:)
  integer :: ig, na
  real(dp) :: arg
  complex(dp) :: phase

  hfi_bare(:) = 0.d0
  IF ( hfi_via_reconstruction_only ) RETURN

  ! transform to reciprocal space (smooth grid)
  allocate(rhoaux(dffts%nnr))
  rhoaux(:) = cmplx(rho_s(:), kind=dp)
  CALL fwfft('Rho', rhoaux, dffts)

  ! fourier transform on the atomic position
  do na = 1, nat
      do ig = gstart, ngms
          arg = sum(tau(1:3,na) * g(1:3,ig)) * tpi
          phase = cmplx(cos(arg),sin(arg), kind=dp)
          hfi_bare(na) = hfi_bare(na) + real(rhoaux(dffts%nl(ig)) * phase, kind=dp)
      enddo
  enddo
  call mp_sum(hfi_bare, intra_pool_comm)

  deallocate(rhoaux)
  return
END SUBROUTINE hfi_fc_bare_el



!-----------------------------------------------------------------------
SUBROUTINE hfi_fc_gipaw_correction(fc_gipaw)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate the GIPAW (all-electron reconstruction) contribution
  ! ... to the Fermi-contact hyperfine coupling.
  !
  USE kinds,                 ONLY : dp
  USE parameters,            ONLY : ntypx
  USE atom,                  ONLY : rgrid, msh
  USE io_global,             ONLY : stdout
  USE gvect,                 ONLY : g, ngm
  USE klist,                 ONLY : nks, xk, igk_k, ngk
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp, atm
  USE wvfct,                 ONLY : g2kin, nbnd
  USE gvecw,                 ONLY : gcutw
  USE wavefunctions,         ONLY : evc
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp,                    ONLY : mp_sum
  USE buffers,               ONLY : get_buffer
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE lsda_mod,              ONLY : current_spin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE control_flags,         ONLY : iverbosity
  USE gipaw_module,          ONLY : alpha, hfi_via_reconstruction_only, use_rt_avg

  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: fc_gipaw(nat)
  !-- local variables ----------------------------------------------------
  real(dp), allocatable :: at_hfi(:,:,:)
  real(dp), allocatable :: work(:)
  integer :: j, nt, ibnd, il1, il2, ik, nbs1, nbs2, kkpsi
  integer :: nrt, inrt
  real(dp) :: mav, mavac
  integer :: m1, m2, lm1, lm2, l1, l2, nrc
  integer :: ijkb0, ih, jh, na, ikb, jkb
  integer :: s_weight, r_first
  real(dp) :: r_Thomson
  complex(dp) :: bec_product
  integer, external :: atomic_number
  integer :: npw

  allocate( at_hfi(paw_nkb,paw_nkb,ntyp) )
  at_hfi = 0.0_dp

  ! calculate radial integrals: <aephi|aephi> - <psphi|psphi> at the origin
  ! (or averaged over the Thomson sphere)
  do nt = 1, ntyp
     kkpsi = paw_recon(nt)%aephi(1)%kkpsi
     allocate( work(kkpsi) )

     r_first = 1
     if ( abs ( rgrid(nt)%r(1) ) < 1d-8 ) r_first = 2
     r_thomson = atomic_number(atm(nt)) * alpha**2
     nrt = COUNT ( rgrid(nt)%r(1:msh(nt)) <= r_thomson )
     nrt = max(nrt, r_first)

     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc
        l1 = paw_recon(nt)%psphi(il1)%label%l
        if (l1 /= 0) cycle

        do il2 = 1, paw_recon(nt)%paw_nbeta
           l2 = paw_recon(nt)%psphi(il2)%label%l
           if (l2 /= 0) cycle

           work = 0.0_dp
           IF ( hfi_via_reconstruction_only ) THEN
              do j = r_first, nrc
                 work(j) = &
                      ( paw_recon(nt)%aephi(il1)%psi(j) &
                      * paw_recon(nt)%aephi(il2)%psi(j) ) &
                      / rgrid(nt)%r(j) ** 2 / fpi
              end do
           ELSE
              do j = r_first, nrc
                 work(j) = &
                      ( paw_recon(nt)%aephi(il1)%psi(j) &
                      * paw_recon(nt)%aephi(il2)%psi(j) &
                      - paw_recon(nt)%psphi(il1)%psi(j) &
                      * paw_recon(nt)%psphi(il2)%psi(j) ) &
                      / rgrid(nt)%r(j) ** 2 / fpi
              enddo
           END IF

           if (use_rt_avg) then
               ! average density over the Thomson sphere
               mavac = 0.0_dp
               do inrt = r_first, nrt
                   mavac = mavac + work(inrt)
               enddo
               mav = mavac / max(nrt - r_first + 1, 1)
               at_hfi(il1,il2,nt) = mav
           else
               ! density at the "origin"
               at_hfi(il1,il2,nt) = work(r_first)
           end if
        enddo
     enddo

     deallocate ( work )
  enddo

  ! calculate the reconstruction part
  fc_gipaw = 0.d0

  do ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     npw = ngk(ik)

     ! magnetization: + for spin-up channel, - for spin-down channel
     if (current_spin == 2) then
        s_weight = -1
     else
        s_weight = +1
     endif

     call gk_sort ( xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin )
     call get_buffer ( evc, nwordwfc, iunwfc, ik)

     call init_gipaw_2 ( npw, igk_k(1,ik), xk(1,ik), paw_vkb )
     call calbec ( npw, paw_vkb, evc, paw_becp, nbnd )

     do ibnd = 1, nbnd
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat

              if ( ityp(na) == nt ) then
                 do ih = 1, paw_recon(nt)%paw_nh
                    ikb = ijkb0 + ih
                    nbs1 = paw_recon(nt)%paw_indv(ih)
                    l1 = paw_recon(nt)%paw_nhtol(ih)
                    m1 = paw_recon(nt)%paw_nhtom(ih)
                    lm1 = m1 + l1**2
                    if (l1 /= 0) cycle

                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh)
                       l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh)
                       lm2 = m2 + l2**2
                       if (l2 /= 0) cycle

                       bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))
                       fc_gipaw(na) = fc_gipaw(na) + s_weight * at_hfi(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
                    enddo
                 enddo

                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              endif

           enddo
        enddo
     enddo

  enddo

  call mp_sum( fc_gipaw, inter_pool_comm )

  deallocate( at_hfi )

END SUBROUTINE hfi_fc_gipaw_correction
