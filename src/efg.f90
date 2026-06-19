!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! EFG tensor and NMR quadrupolar coupling constants.
! Adapted from qe-gipaw (D. Ceresoli et al.) for QE-CONVERSE.
!
! The EFG is a ground-state property; call calc_efg() after gipaw_setup
! and before newscf, using the charge density and wavefunctions from read_file.
!

!-----------------------------------------------------------------------
SUBROUTINE calc_efg()
  !-----------------------------------------------------------------------
  !
  ! Compute the electric field gradient (EFG) tensor at each atomic site.
  ! Three contributions:
  !   efg_bare  - electronic valence contribution (G-space formula)
  !   efg_ion   - ionic contribution (Ewald sum)
  !   efg_gipaw - PAW reconstruction correction
  !
  USE kinds,         ONLY : dp
  USE io_global,     ONLY : stdout
  USE constants,     ONLY : pi, tpi, fpi, angstrom_au, rytoev, electronvolt_si
  USE fft_base,      ONLY : dffts, dfftp
  USE scf,           ONLY : rho
  USE lsda_mod,      ONLY : nspin
  USE ions_base,     ONLY : nat, atm, ityp, zv
  USE symme,         ONLY : symtensor
  USE mp,            ONLY : mp_sum
  USE gipaw_module,  ONLY : q_efg, i_efg
  USE io_files,      ONLY : iunwfc, nwordwfc
  USE control_flags, ONLY : io_level
  USE buffers,       ONLY : open_buffer, close_buffer
  USE efg_mod,       ONLY : efg_tensor

  IMPLICIT NONE

  real(dp), allocatable :: rho_valence(:)
  real(dp), allocatable :: efg_bare(:,:,:), efg_ion(:,:,:)
  real(dp), allocatable :: efg_gipaw(:,:,:), efg_tot(:,:,:)
  real(dp), allocatable :: zion(:)
  complex(dp), allocatable :: tmp(:,:,:)
  integer :: alpha, beta, na
  real(dp) :: v(3), axis(3,3), eta, Cq, spinI, denom, nu_Q
  logical :: exst

  write(stdout,*)
  write(stdout,'(5X,A)') 'Computing EFG tensors ...'
  write(stdout,*)

  allocate( efg_bare(3,3,nat), efg_ion(3,3,nat) )
  allocate( efg_gipaw(3,3,nat), efg_tot(3,3,nat) )

  ! --- electronic (valence) contribution ---
  ! Use total charge density from read_file. For nspin=2, rho%of_r(:,1)
  ! is the total charge (rho_up + rho_down) in (rho,zeta) representation.
  allocate( rho_valence(dfftp%nnr) )
  rho_valence(:) = rho%of_r(:,1)
  call efg_bare_el(rho_valence, efg_bare)
  deallocate( rho_valence )

  ! --- ionic contribution (Ewald sum) ---
  allocate( zion(nat), tmp(nat,3,3) )
  do na = 1, nat
    zion(na) = zv(ityp(na))
  enddo
  call ewald_dipole(tmp, zion)
  do na = 1, nat
    efg_ion(:,:,na) = real(tmp(na,:,:), kind=dp)
  enddo
  deallocate( tmp, zion )

  ! --- GIPAW reconstruction correction ---
  ! read_file (io_level=1) opens then closes iunwfc; reopen it so efg_correction
  ! can call get_buffer for the ground-state wavefunctions, then close again.
  call open_buffer(iunwfc, 'wfc', nwordwfc, io_level, exst)
  call efg_correction(efg_gipaw)
  call close_buffer(iunwfc, 'KEEP')

  ! --- sum contributions and symmetrize ---
  do na = 1, nat
    efg_tot(:,:,na) = efg_bare(:,:,na) + efg_ion(:,:,na) + efg_gipaw(:,:,na)
  enddo

  write(stdout,'(5X,A)') '----- bare electronic EFG (Ha/bohr^2) -----'
  do na = 1, nat
    do beta = 1, 3
      write(stdout,'(5X,A,I3,2X,3(F14.6,2X))') atm(ityp(na)), na, &
            (efg_bare(alpha,beta,na), alpha=1,3)
    enddo
    write(stdout,*)
  enddo

  write(stdout,'(5X,A)') '----- ionic EFG (Ha/bohr^2) -----'
  do na = 1, nat
    do beta = 1, 3
      write(stdout,'(5X,A,I3,2X,3(F14.6,2X))') atm(ityp(na)), na, &
            (efg_ion(alpha,beta,na), alpha=1,3)
    enddo
    write(stdout,*)
  enddo

  write(stdout,'(5X,A)') '----- GIPAW correction (Ha/bohr^2) -----'
  do na = 1, nat
    do beta = 1, 3
      write(stdout,'(5X,A,I3,2X,3(F14.6,2X))') atm(ityp(na)), na, &
            (efg_gipaw(alpha,beta,na), alpha=1,3)
    enddo
    write(stdout,*)
  enddo

  write(stdout,'(5X,A)') '----- total EFG (Ha/bohr^2) -----'
  do na = 1, nat
    do beta = 1, 3
      write(stdout,'(5X,A,I3,2X,3(F14.6,2X))') atm(ityp(na)), na, &
            (efg_tot(alpha,beta,na), alpha=1,3)
    enddo
    write(stdout,*)
  enddo

  call symtensor(nat, efg_tot)

  write(stdout,'(5X,A)') '----- total EFG symmetrized (Ha/bohr^2) -----'
  do na = 1, nat
    do beta = 1, 3
      write(stdout,'(5X,A,I3,2X,3(F14.6,2X))') atm(ityp(na)), na, &
            (efg_tot(alpha,beta,na), alpha=1,3)
    enddo
    write(stdout,*)
  enddo

  ! --- spectroscopic parameters ---
  write(stdout,*)
  write(stdout,'(5X,A)') 'NMR/NQR QUADRUPOLAR PARAMETERS:'
  write(stdout,'(5X,A)') 'Vxx, Vyy, Vzz: EFG principal values (Ha/bohr^2), ordered |Vzz|>=|Vyy|>=|Vxx|'
  write(stdout,'(5X,A)') 'axis: corresponding eigenvectors in Cartesian crystal coordinates (a,b,c)'
  write(stdout,'(5X,A)') 'Q: nuclear quadrupole moment (input, 1e-30 m^2 = 10 mbarn),  I: nuclear spin (input)'
  write(stdout,'(5X,A)') 'Phi_zz = Vzz (largest |eigenvalue| of EFG)'
  write(stdout,'(5X,A)') 'Cq = e*Q*Phi_zz/h (MHz),  eta = (Vxx-Vyy)/Vzz'
  write(stdout,'(5X,A)') 'nu_Q = 3*e*Q*Phi_zz / (2I(2I-1)h) = 3*Cq / (2I(2I-1))  (MHz)'
  write(stdout,*)

  do na = 1, nat
    call principal_axis(efg_tot(:,:,na), v, axis)

    write(stdout,'(5X,A,I3,4X,A,F10.4,4X,A,3F10.6,A)') &
          atm(ityp(na)), na, 'Vxx=', v(1), 'axis=(', axis(1:3,1), ')'
    write(stdout,'(5X,A,I3,4X,A,F10.4,4X,A,3F10.6,A)') &
          atm(ityp(na)), na, 'Vyy=', v(2), 'axis=(', axis(1:3,2), ')'
    write(stdout,'(5X,A,I3,4X,A,F10.4,4X,A,3F10.6,A)') &
          atm(ityp(na)), na, 'Vzz=', v(3), 'axis=(', axis(1:3,3), ')'

    eta = 0.d0
    if (abs(v(3)) > 1d-5) eta = (v(1) - v(2)) / v(3)

    if (abs(q_efg(ityp(na))) > 1d-10) then
      ! Cq[MHz] = e * Q * Phi_zz / h, Phi_zz = Vzz (largest EFG eigenvalue).
      ! angstrom_au = 1/bohr_radius_angs = 1.88973, so with Phi_zz in Ha/bohr^2 and
      ! Q in 1e-30 m^2 this evaluates to  Cq[MHz] = 2.3499 * Q[1e-30 m^2] * Phi_zz[a.u.]
      ! (identical to the reference qe-gipaw efg.f90).
      Cq = v(3) * q_efg(ityp(na)) * rytoev * 2.d0 * angstrom_au**2 &
           * electronvolt_si * 1.d18 / 6.62620d0
      write(stdout,'(5X,A,I3,2X,A,F10.4,A,2X,A,F12.4,A,2X,A,F8.5)') &
            atm(ityp(na)), na, &
            'Q=', q_efg(ityp(na)), ' 1e-30 m^2', &
            'Cq=', Cq, ' MHz', &
            'eta=', eta

      ! quadrupolar frequency nu_Q = 3*e*Q*Phi_zz/(2I(2I-1)h) = 3*Cq/(2I(2I-1))
      ! (Phi_zz = Vzz = v(3) is carried inside Cq); defined only for I >= 1
      spinI = i_efg(ityp(na))
      denom = 2.d0 * spinI * (2.d0 * spinI - 1.d0)
      if (denom > 1d-10) then
        nu_Q = 3.d0 * Cq / denom
        write(stdout,'(5X,A,I3,2X,A,F6.1,4X,A,F12.4,A)') &
              atm(ityp(na)), na, 'I=', spinI, 'nu_Q=', nu_Q, ' MHz'
      endif
    else
      write(stdout,'(5X,A,I3,2X,A,F8.5)') atm(ityp(na)), na, 'eta=', eta
    endif
    write(stdout,*)
  enddo

  ! Save symmetrized tensor for summary printing at end of run
  allocate( efg_tensor(3,3,nat) )
  efg_tensor(:,:,:) = efg_tot(:,:,:)

  deallocate( efg_bare, efg_ion, efg_gipaw, efg_tot )

END SUBROUTINE calc_efg


!-----------------------------------------------------------------------
SUBROUTINE print_efg_summary()
  !-----------------------------------------------------------------------
  !
  ! Print Cq and eta for all atoms at end of run, alongside chemical shifts.
  !
  USE kinds,        ONLY : dp
  USE io_global,    ONLY : stdout
  USE constants,    ONLY : angstrom_au, rytoev, electronvolt_si
  USE ions_base,    ONLY : nat, atm, ityp
  USE gipaw_module, ONLY : q_efg, i_efg
  USE efg_mod,      ONLY : efg_tensor

  IMPLICIT NONE
  integer :: na
  real(dp) :: v(3), axis(3,3), eta, Cq, spinI, denom, nu_Q

  if (.not. allocated(efg_tensor)) return

  write(stdout,*)
  write(stdout,'(5X,A)') '=========== NMR/NQR QUADRUPOLAR PARAMETERS ==========='
  write(stdout,'(5X,A)') 'Phi_zz = Vzz (largest |eigenvalue| of EFG);  Q: input quadrupole moment (1e-30 m^2)'
  write(stdout,'(5X,A)') 'Cq = e*Q*Phi_zz/h (MHz),  eta = (Vxx-Vyy)/Vzz'
  write(stdout,'(5X,A)') 'nu_Q = 3*e*Q*Phi_zz / (2I(2I-1)h) = 3*Cq / (2I(2I-1))  (MHz),  I: nuclear spin (nu_Q printed when I >= 1)'
  write(stdout,*)

  do na = 1, nat
    call principal_axis(efg_tensor(:,:,na), v, axis)
    eta = 0.d0
    if (abs(v(3)) > 1d-5) eta = (v(1) - v(2)) / v(3)

    if (abs(q_efg(ityp(na))) > 1d-10) then
      ! Cq[MHz] = 2.3499 * Q[1e-30 m^2] * Phi_zz[a.u.]  (see calc_efg for derivation)
      Cq = v(3) * q_efg(ityp(na)) * rytoev * 2.d0 * angstrom_au**2 &
           * electronvolt_si * 1.d18 / 6.62620d0
      spinI = i_efg(ityp(na))
      denom = 2.d0 * spinI * (2.d0 * spinI - 1.d0)
      if (denom > 1d-10) then
        nu_Q = 3.d0 * Cq / denom
        write(stdout,'(5X,A,I3,2X,A,F10.4,A,2X,A,F12.4,A,2X,A,F8.5,2X,A,F6.1,2X,A,F12.4,A)') &
              atm(ityp(na)), na, &
              'Q=', q_efg(ityp(na)), ' 1e-30 m^2', &
              'Cq=', Cq, ' MHz', 'eta=', eta, &
              'I=', spinI, 'nu_Q=', nu_Q, ' MHz'
      else
        write(stdout,'(5X,A,I3,2X,A,F10.4,A,2X,A,F12.4,A,2X,A,F8.5)') &
              atm(ityp(na)), na, &
              'Q=', q_efg(ityp(na)), ' 1e-30 m^2', &
              'Cq=', Cq, ' MHz', 'eta=', eta
      endif
    else
      write(stdout,'(5X,A,I3,2X,A,F10.4,4X,A,F8.5)') &
            atm(ityp(na)), na, 'Vzz=', v(3), 'eta=', eta
    endif
  enddo
  write(stdout,*)

END SUBROUTINE print_efg_summary


!-----------------------------------------------------------------------
SUBROUTINE efg_bare_el(rho_in, efg_bare)
  !-----------------------------------------------------------------------
  !
  ! Electronic contribution to the EFG in G-space.
  ! V_{alpha,beta}(R) = (4pi/Omega) Sum_{G/=0} rho(G) *
  !     (G_alpha G_beta / G^2 - delta_{alpha,beta}/3) * exp(iG.R) / G^2
  ! computed in Hartree units (e2 = 1).
  !
  USE kinds,       ONLY : dp
  USE mp,          ONLY : mp_sum
  USE mp_pools,    ONLY : intra_pool_comm
  USE constants,   ONLY : tpi, fpi
  USE gvect,       ONLY : g, gg, gstart, ngm
  USE gvecs,       ONLY : ngms
  USE fft_base,    ONLY : dffts
  USE fft_interfaces, ONLY : fwfft
  USE ions_base,   ONLY : nat, tau

  IMPLICIT NONE
  real(dp), intent(in)  :: rho_in(dffts%nnr)
  real(dp), intent(out) :: efg_bare(3,3,nat)

  complex(dp), allocatable :: efg_g(:,:,:), rhoaux(:)
  integer :: alpha, beta, ig, na
  real(dp) :: arg, trace
  real(dp), parameter :: e2 = 1.0_dp   ! Hartree units
  complex(dp) :: phase

  allocate( efg_g(ngms,3,3), rhoaux(dffts%nnr) )
  efg_g(:,:,:) = (0.d0, 0.d0)
  rhoaux(:) = cmplx(rho_in(:), 0.d0, kind=dp)

  call fwfft('Rho', rhoaux, dffts)

  do ig = gstart, ngms
    trace = gg(ig) / 3.d0
    do alpha = 1, 3
      efg_g(ig,alpha,alpha) = -trace
      do beta = 1, 3
        efg_g(ig,alpha,beta) = ( efg_g(ig,alpha,beta) + &
            g(alpha,ig) * g(beta,ig) ) * fpi * e2 * rhoaux(dffts%nl(ig)) / gg(ig)
      enddo
    enddo
  enddo

  efg_bare(:,:,:) = 0.d0
  do alpha = 1, 3
    do beta = 1, 3
      do na = 1, nat
        do ig = gstart, ngms
          arg = tpi * ( tau(1,na)*g(1,ig) + tau(2,na)*g(2,ig) + tau(3,na)*g(3,ig) )
          phase = cmplx(cos(arg), sin(arg), kind=dp)
          efg_bare(alpha,beta,na) = efg_bare(alpha,beta,na) + &
               real(efg_g(ig,alpha,beta) * phase, kind=dp)
        enddo
      enddo
    enddo
  enddo

#ifdef __MPI
  call mp_sum(efg_bare, intra_pool_comm)
#endif

  deallocate(efg_g, rhoaux)

END SUBROUTINE efg_bare_el


!-----------------------------------------------------------------------
SUBROUTINE efg_correction(efg_corr)
  !-----------------------------------------------------------------------
  !
  ! GIPAW reconstruction correction to the EFG.
  ! Radial integrals: <aephi|1/r^3|aephi> - <psphi|1/r^3|psphi>
  ! contracted with the l=2 Gaunt coefficients ap(lm=5..9, lm1, lm2).
  !
  USE io_files,    ONLY : nwordwfc, iunwfc
  USE kinds,       ONLY : dp
  USE uspp,        ONLY : ap
  USE parameters,  ONLY : ntypx
  USE atom,        ONLY : rgrid
  USE gvect,       ONLY : g, ngm
  USE klist,       ONLY : nks, xk, igk_k, ngk
  USE cell_base,   ONLY : tpiba2
  USE ions_base,   ONLY : nat, ityp, ntyp => nsp
  USE wvfct,       ONLY : g2kin, current_k, wg
  USE gvecw,       ONLY : gcutw
  USE lsda_mod,    ONLY : current_spin, isk
  USE wavefunctions, ONLY : evc
  USE paw_gipaw,   ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE becmod,      ONLY : calbec
  USE constants,   ONLY : pi, fpi
  USE buffers,     ONLY : get_buffer
  USE gipaw_module, ONLY : nbnd_occ
  USE mp_pools,    ONLY : inter_pool_comm
  USE mp,          ONLY : mp_sum

  IMPLICIT NONE
  real(dp), intent(out) :: efg_corr(3,3,nat)

  integer :: j, nt, ibnd, il1, il2, ik, nbs1, nbs2, kkpsi
  integer :: lm, m1, m2, lm1, lm2, l1, l2, nrc
  integer :: ijkb0, ih, jh, na, ikb, jkb, r_first
  integer :: npw
  complex(dp) :: bec_product
  real(dp), allocatable :: at_efg(:,:,:), work(:)
  complex(dp), allocatable :: efg_c(:,:)

  allocate( efg_c(9,nat) )
  efg_c = (0.d0, 0.d0)

  allocate( at_efg(paw_nkb, paw_nkb, ntypx) )
  at_efg = 0.d0

  ! radial integrals: <aephi|1/r^3|aephi> - <psphi|1/r^3|psphi>
  do nt = 1, ntyp
    kkpsi = paw_recon(nt)%aephi(1)%kkpsi
    allocate( work(kkpsi) )

    r_first = 1
    if ( abs(rgrid(nt)%r(1)) < 1d-8 ) r_first = 2

    do il1 = 1, paw_recon(nt)%paw_nbeta
      nrc = paw_recon(nt)%psphi(il1)%label%nrc

      do il2 = 1, paw_recon(nt)%paw_nbeta
        work = 0.d0
        do j = r_first, nrc
          work(j) = ( paw_recon(nt)%aephi(il1)%psi(j) &
                    * paw_recon(nt)%aephi(il2)%psi(j) &
                    - paw_recon(nt)%psphi(il1)%psi(j) &
                    * paw_recon(nt)%psphi(il2)%psi(j) ) &
                    / rgrid(nt)%r(j)**3
        enddo
        call simpson(nrc, work, rgrid(nt)%rab, at_efg(il1,il2,nt))
      enddo
    enddo

    deallocate(work)
  enddo

  ! projection onto PAW basis and accumulation
  do ik = 1, nks
    current_k = ik
    current_spin = isk(ik)
    npw = ngk(ik)

    call gk_sort(xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
    call get_buffer(evc, nwordwfc, iunwfc, ik)

    call init_gipaw_2(npw, igk_k(1,ik), xk(1,ik), paw_vkb)
    call calbec(npw, paw_vkb, evc, paw_becp)

    do ibnd = 1, nbnd_occ(ik)
      ijkb0 = 0
      do nt = 1, ntyp
        do na = 1, nat
          if (ityp(na) == nt) then
            do ih = 1, paw_recon(nt)%paw_nh
              ikb = ijkb0 + ih
              nbs1 = paw_recon(nt)%paw_indv(ih)
              l1   = paw_recon(nt)%paw_nhtol(ih)
              m1   = paw_recon(nt)%paw_nhtom(ih)
              lm1  = m1 + l1**2

              do jh = 1, paw_recon(nt)%paw_nh
                jkb = ijkb0 + jh
                nbs2 = paw_recon(nt)%paw_indv(jh)
                l2   = paw_recon(nt)%paw_nhtol(jh)
                m2   = paw_recon(nt)%paw_nhtom(jh)
                lm2  = m2 + l2**2

                bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))

                do lm = 5, 9
                  efg_c(lm,na) = efg_c(lm,na) + &
                       bec_product * at_efg(nbs1,nbs2,nt) * ap(lm,lm1,lm2) * wg(ibnd,ik)
                enddo
              enddo
            enddo
            ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
          endif
        enddo
      enddo
    enddo
  enddo

#ifdef __MPI
  call mp_sum(efg_c, inter_pool_comm)
#endif

  ! convert from spherical harmonics to Cartesian tensor components
  efg_corr(1,1,:) =  sqrt(3.d0) * real(efg_c(8,:)) - real(efg_c(5,:))
  efg_corr(2,2,:) = -sqrt(3.d0) * real(efg_c(8,:)) - real(efg_c(5,:))
  efg_corr(3,3,:) =  2.d0 * real(efg_c(5,:))
  efg_corr(1,2,:) =  sqrt(3.d0) * real(efg_c(9,:))
  efg_corr(2,1,:) =  efg_corr(1,2,:)
  efg_corr(1,3,:) = -sqrt(3.d0) * real(efg_c(6,:))
  efg_corr(3,1,:) =  efg_corr(1,3,:)
  efg_corr(2,3,:) = -sqrt(3.d0) * real(efg_c(7,:))
  efg_corr(3,2,:) =  efg_corr(2,3,:)

  efg_corr = -sqrt(4.d0 * pi / 5.d0) * efg_corr

  deallocate(efg_c, at_efg)

END SUBROUTINE efg_correction
