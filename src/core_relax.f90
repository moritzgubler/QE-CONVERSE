!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Ported to QE-CONVERSE / QE 7.5 from qe-gipaw-6.2 (Knight-shift fork).
! Computes the core-relaxation contribution to the Fermi-contact
! hyperfine coupling (the spin polarization of the core s shells induced
! by the valence spin density), following Bahramy, Sluiter, Kawazoe,
! PRB 76, 035124 (2007), eqs. (15)-(19).
!
! Knight-shift / metal extensions kept from the fork:
!   * spin-resolved band occupations are handled through the wg weights
!     (loop over all nbnd), so partial occupations in metals are correct;
!   * optional averaging of the contact density over the Thomson sphere
!     r_T = Z*alpha^2 (use_rt_avg), more robust than the singular r->0
!     value for metals and heavy nuclei.
!
!-----------------------------------------------------------------------
SUBROUTINE hfi_fc_core_relax(method, fc_core)
  !-----------------------------------------------------------------------
  USE kinds,                 ONLY : dp
  USE constants,             ONLY : pi, tpi, fpi, e2
  USE parameters,            ONLY : ntypx
  USE ions_base,             ONLY : ntyp => nsp, atm, nat, tau, ityp
  USE atom,                  ONLY : rgrid, msh
  USE radial_grids,          ONLY : ndmx
  USE scf,                   ONLY : rho
  USE gvect,                 ONLY : g, ngm
  USE fft_base,              ONLY : dfftp
  USE fft_interfaces,        ONLY : fwfft
  USE lsda_mod,              ONLY : nspin, isk, current_spin
  USE buffers,               ONLY : get_buffer
  USE control_flags,         ONLY : iverbosity
  USE klist,                 ONLY : nks, xk, igk_k, ngk
  USE wvfct,                 ONLY : nbnd, g2kin, wg, current_k
  USE gvecw,                 ONLY : gcutw
  USE becmod,                ONLY : calbec
  USE wavefunctions,         ONLY : evc
  USE io_global,             ONLY : stdout
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE mp_pools,              ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                    ONLY : mp_sum
  USE paw_gipaw,             ONLY : paw_recon, paw_vkb, paw_becp
  USE xc_lib,                ONLY : xc
  USE gipaw_module,          ONLY : alpha, core_relax_r_max, use_rt_avg, nbrx

  !-- parameters --------------------------------------------------------
  IMPLICIT NONE
  integer, intent(in) :: method
  real(dp), intent(out) :: fc_core(nat)

  ! -- constants ---------------------------------------------------------
  real(dp) :: r_max                         ! max core radius (from input)
  integer, parameter :: n_max = 10          ! max number of s orbitals

  !-- local variables ----------------------------------------------------
  real(dp) :: eigenvalue(n_max,ntypx)
  real(dp) :: ae_orb(ndmx,n_max,ntypx)
  real(dp), allocatable :: vpot(:)
  integer :: nt, nn, zz, nstop, nin
  integer, external :: atomic_number

  integer :: nrt, inrt, nr_max
  real(dp) :: mav, mavac, cav, r_thomson

  integer :: s_maj, s_min, na, ispin, j
  complex(dp), allocatable :: aux(:), rho_g(:)
  real(dp), allocatable :: sph_rho_bare(:,:), sph_rho_gipaw(:,:)
  real(dp) :: sph_rho_core(ndmx)

  integer :: il1, il2, nrc, l1, l2
  real(dp), allocatable :: rho_recon(:,:,:,:)

  integer :: ik, ibnd, ih, jh, ikb, jkb, m1, m2, lm1, lm2, ijkb0, nbs1, nbs2
  complex(dp) :: bec_product

  real(dp), allocatable :: delta_v(:,:), work(:)
  integer :: n1, n2, ncore, r_first
  real(dp) :: b(2), coeff, norm, contrib
  integer :: mesh
  real(dp) :: arho, zeta
  ! XClib (array) interface for methods 2 and 3
  real(dp) :: xc_rho(1,2), xc_ex(1), xc_ec(1), xc_vx(1,2), xc_vc(1,2)
  integer :: npw

  if (method < 1 .or. method > 3) call errore('core-relax', 'unknown method', abs(method))

  call start_clock('core_relax')
  fc_core = 0.d0

  !====================================================================
  ! recalculate AE orbitals (s only)
  !====================================================================
  write(stdout,'(5X,''Warning: core-relaxation is an experimental feature'')')
  if (iverbosity > 1) write(stdout,'(5X,''core-relax: calculating AE orbitals'')')
  if (iverbosity > 1) write(stdout,'(5X,''core-relax: method = '',I1)') method
  if (method == 3) write(stdout,'(5X,''core-relax: method 3 uses LDA XC only '', &
       &''(radial GGA gradient term omitted)'')')
  eigenvalue(:,:) = 0.d0
  ae_orb(:,:,:) = 0.d0
  do nt = 1, ntyp
      if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle
      allocate( vpot(rgrid(nt)%mesh), work(rgrid(nt)%mesh) )
      do nn = 1, n_max
          ! setup grid
          zz = atomic_number(atm(nt))
          rgrid(nt)%dx = log( rgrid(nt)%r(2)/rgrid(nt)%r(1) )
          rgrid(nt)%r2 = rgrid(nt)%r**2
          rgrid(nt)%sqr = sqrt(rgrid(nt)%r)

          ! setup potential
          vpot(:) = paw_recon(nt)%gipaw_ae_vloc(:) / rgrid(nt)%r(:)

          ! solve radial schroedinger equation (s only, l=0), non-relativistic
          eigenvalue(nn,nt) = -(dble(zz)/dble(nn))**2.0
          nin = 0
          nstop = 0
          call lschps(2, 2.d0*zz, 1d-12, rgrid(nt), nin, nn, 0, &
                      eigenvalue(nn,nt), vpot, ae_orb(1,nn,nt), nstop)
          if (nstop /= 0 .and. nstop /= 50) then
              eigenvalue(nn,nt) = 0.d0
              exit
          endif

          if (iverbosity > 0) then
              do j = 1, rgrid(nt)%mesh
                  work(j) = ae_orb(j,nn,nt)**2.d0
              enddo
              call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, norm)
              write(stdout,'(5X,A,4X,I1,''S   eig='',F12.4,'' Ry   norm='',F12.4)') &
                    atm(nt), nn, eigenvalue(nn,nt), norm
          endif

      enddo
      deallocate( vpot, work )
      if (iverbosity > 1) write(stdout,*)
  enddo

  ! Spin channels: 1 = up, 2 = down (QE convention isk/current_spin).
  ! The Fermi-contact coupling is driven by the magnetization n_up - n_dn,
  ! so we work with per-spin (up,down) densities throughout and form the
  ! difference explicitly; no majority/minority relabelling is needed.
  s_maj = 1
  s_min = 2

  ! prepare reconstruction terms: (aephi*aephi - psphi*psphi)/(r^2 4pi)
  allocate( rho_recon(ndmx,nbrx,nbrx,ntyp) )
  rho_recon = 0.d0
  do nt = 1, ntyp
    do il1 = 1, paw_recon(nt)%paw_nbeta
      nrc = paw_recon(nt)%psphi(il1)%label%nrc
      l1 = paw_recon(nt)%psphi(il1)%label%l
      if (l1 /= 0) cycle
      do il2 = 1, paw_recon(nt)%paw_nbeta
        l2 = paw_recon(nt)%psphi(il2)%label%l
        if (l2 /= 0) cycle
        do j = 1, nrc
          rho_recon(j,il1,il2,nt) = &
                   ( paw_recon(nt)%aephi(il1)%psi(j) &
                   * paw_recon(nt)%aephi(il2)%psi(j) &
                   - paw_recon(nt)%psphi(il1)%psi(j) &
                   * paw_recon(nt)%psphi(il2)%psi(j) ) &
                   / rgrid(nt)%r(j) ** 2 / fpi
        enddo
      enddo ! il2
    enddo ! il1
   enddo ! nt

  !====================================================================
  ! loop over atoms with core electrons
  !====================================================================
  allocate( sph_rho_bare(ndmx,nspin) )
  allocate( sph_rho_gipaw(ndmx,nspin) )
  allocate( aux(dfftp%nnr), rho_g(ngm) )
  allocate( delta_v(ndmx,nat) )
  delta_v = 0.d0

  do na = 1, nat
    nt = ityp(na)
    if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle
    if (iverbosity > 1) write(stdout,'(5X,''core-relax: projecting around atom '',I3)') na

    !====================================================================
    ! setting the max core radius
    !====================================================================
    r_max = core_relax_r_max
    nr_max = COUNT ( rgrid(nt)%r(1:msh(nt)) <= r_max )
    if (iverbosity > 1) then
       write(stdout,*)
       write(stdout,'(5X,''core-relax: max core radius (r_max) =  '',F12.4)') r_max
       write(stdout,'(5X,''core-relax: number of points inside (nr_max) =  '',I4)') nr_max
       write(stdout,*)
    endif

    !====================================================================
    ! project the density around each atom
    !====================================================================
    sph_rho_bare = 0.d0
    do ispin = 1, nspin
      aux(1:dfftp%nnr) = rho%of_r(1:dfftp%nnr,ispin)
      call fwfft('Rho',aux,dfftp)
      rho_g(1:ngm) = aux(dfftp%nl(1:ngm))
      call spherical_average(rgrid(nt)%mesh, rgrid(nt)%r, tau(1,na), r_max, rho_g, sph_rho_bare(1,ispin))
    enddo
    call mp_sum(sph_rho_bare, intra_pool_comm)
    ! In QE LSDA, rho%of_r(:,1)=rho_up+rho_dn (total) and rho%of_r(:,2)=rho_up-rho_dn
    ! (magnetization).  The reconstruction (sph_rho_gipaw) and core (sph_rho_core)
    ! densities below are per-spin (up,down), so convert the bare density to (up,down)
    ! too, to keep a single consistent convention in the b(1:2) arithmetic.
    do j = 1, rgrid(nt)%mesh
       arho = sph_rho_bare(j,1)            ! total
       zeta = sph_rho_bare(j,2)            ! magnetization (reuse as temporary)
       sph_rho_bare(j,1) = 0.5d0*(arho + zeta)   ! up
       sph_rho_bare(j,2) = 0.5d0*(arho - zeta)   ! down
    enddo

    !====================================================================
    ! do GIPAW reconstruction
    !====================================================================
    sph_rho_gipaw = 0.d0
    do ik = 1, nks
      current_k = ik
      current_spin = isk(ik)
      npw = ngk(ik)

      call gk_sort (xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
      call get_buffer (evc, nwordwfc, iunwfc, ik)
      call init_gipaw_2 (npw, igk_k(1,ik), xk(1,ik), paw_vkb)
      call calbec (npw, paw_vkb, evc, paw_becp, nbnd)

      do ibnd = 1, nbnd
        ijkb0 = 0
        do nt = 1, ntyp
              if ( ityp(na) == nt ) then
                 do ih = 1, paw_recon(nt)%paw_nh
                    ikb = ijkb0 + ih
                    nbs1 = paw_recon(nt)%paw_indv(ih)
                    l1 = paw_recon(nt)%paw_nhtol(ih)
                    m1 = paw_recon(nt)%paw_nhtom(ih)
                    lm1 = m1 + l1**2
                    nrc = paw_recon(nt)%psphi(nbs1)%label%nrc
                    if (l1 /= 0) cycle

                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh)
                       l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh)
                       lm2 = m2 + l2**2
                       if (l2 /= 0) cycle

                       bec_product = paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))

                       sph_rho_gipaw(1:nrc,current_spin) = &
                            sph_rho_gipaw(1:nrc,current_spin) + &
                            rho_recon(1:nrc,nbs1,nbs2,nt) * &
                            bec_product * wg(ibnd,ik)
                    end do  ! jh
                 end do  ! ih
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              end if
        end do  ! nt
      end do  ! ibnd
    end do  ! ik
    call mp_sum(sph_rho_gipaw, inter_pool_comm)

    !====================================================================
    ! compute perturbing potential
    !====================================================================
    nt = ityp(na)
    if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle

    mesh = rgrid(nt)%mesh

    ! contribution of the core orbitals, eq.(19) of PRB 76, 035124
    sph_rho_core = 0.d0
    do n1 = 1, paw_recon(nt)%gipaw_ncore_orbital
      if (paw_recon(nt)%gipaw_core_orbital_l(n1) /= 0) cycle
      do j = 1, mesh
        sph_rho_core(j) = sph_rho_core(j) + 2.d0*paw_recon(nt)%gipaw_core_orbital(j,n1)**2.d0 / rgrid(nt)%r2(j) / fpi
      enddo
    enddo

    allocate(work(mesh))
    do j = 1, mesh
      if (rgrid(nt)%r(j) > r_max) cycle

      ! calculate density and magnetization at radial point j (eq.19)
      b(1:2) = sph_rho_bare(j,1:2) + sph_rho_gipaw(j,1:2) + sph_rho_core(j)
      arho = abs(b(1)+b(2))
      if (arho < 1.d-30) cycle
      zeta = (b(s_maj)-b(s_min))/arho

      ! compute the perturbing potential, three possibilities
      select case (method)
         case (1) ! simple local exchange: eq.(20) of PRB 76, 035124
         delta_v(j,na) = -(e2*2.d0/pi) * (b(s_maj)-b(s_min)) / arho**(2.d0/3.d0)

         case (2) ! Exchange only (LDA), via XClib
         xc_rho(1,1) = arho
         xc_rho(1,2) = b(s_maj) - b(s_min)
         call xc(1, 2, 2, xc_rho, xc_ex, xc_ec, xc_vx, xc_vc)
         delta_v(j,na) = e2*(xc_vx(1,1) - xc_vx(1,2))

         case (3) ! Full XC (LDA exchange + correlation), via XClib
         xc_rho(1,1) = arho
         xc_rho(1,2) = b(s_maj) - b(s_min)
         call xc(1, 2, 2, xc_rho, xc_ex, xc_ec, xc_vx, xc_vc)
         delta_v(j,na) = e2*((xc_vx(1,1)+xc_vc(1,1)) - (xc_vx(1,2)+xc_vc(1,2)))
      end select
    enddo
    deallocate(work)

  !====================================================================
  ! end of the loop over atoms
  !====================================================================
  enddo
  deallocate(aux, rho_g, sph_rho_bare, sph_rho_gipaw)
  deallocate(rho_recon)

  !====================================================================
  ! now, do the core relaxation via perturbation theory (PRB 76, 035124)
  !====================================================================
  allocate( work(ndmx) )
  fc_core(1:nat) = 0.d0

  do na = 1, nat
    nt = ityp(na)
    if (paw_recon(nt)%gipaw_ncore_orbital == 0) cycle

    do j = 1, rgrid(nt)%mesh
      if (rgrid(nt)%r(j) > 1d-5) then
        r_first = j
        exit
      endif
    enddo

    r_thomson = atomic_number(atm(nt)) * alpha**2
    nrt = COUNT ( rgrid(nt)%r(1:msh(nt)) <= r_thomson )
    if (iverbosity > 1) then
       write(stdout,*)
       write(stdout,'(5X,'' the computed rT is '',F12.6)') r_thomson
       write(stdout,'(5X,'' the number of points inside rT is '',I6)') nrt
       write(stdout,*)
    endif

    ! count number of s core orbitals
    ncore = 0
    do n1 = 1, paw_recon(nt)%gipaw_ncore_orbital
      if (paw_recon(nt)%gipaw_core_orbital_l(n1) /= 0) cycle
      ncore = ncore + 1
    enddo

    if (iverbosity > 1) then
      do n1 = 1, ncore
        work = 0.d0
        do j = 1, rgrid(nt)%mesh
          work(j) = ae_orb(j,n1,nt) * delta_v(j,na) * ae_orb(j,n1,nt)
        end do
        call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, coeff)
        write(stdout,'(5X,A,I3,2X,I1,''S splitting:'',F12.6)') atm(ityp(na)), na, n1, coeff
      enddo
    endif

    do n1 = 1, ncore
      do n2 = n1+1, n_max
        if (eigenvalue(n2, nt) == 0.d0) cycle  ! unbound

        work = 0.d0
        do j = 1, rgrid(nt)%mesh
          work(j) = ae_orb(j,n1,nt) * delta_v(j,na) * ae_orb(j,n2,nt)
        end do

        call simpson(rgrid(nt)%mesh, work, rgrid(nt)%rab, coeff)
        contrib = 2.d0 * 4.d0 * ae_orb(r_first,n1,nt) * ae_orb(r_first,n2,nt) * &
                  coeff / (eigenvalue(n1,nt) - eigenvalue(n2,nt)) / &
                  rgrid(nt)%r(r_first)**2 / fpi

        if (use_rt_avg) then
           ! average the contact density over the Thomson sphere
           mavac = 0.0_dp
           do inrt = r_first, nrt
               cav = 2.d0 * 4.d0 * ae_orb(inrt,n1,nt) * ae_orb(inrt,n2,nt) * &
                         coeff / (eigenvalue(n1,nt) - eigenvalue(n2,nt)) / &
                         rgrid(nt)%r(inrt)**2 / fpi
               mavac = mavac + cav
           enddo
           mav = mavac / max(nrt,1)
           fc_core(na) = fc_core(na) + mav
        else
           ! density at the "origin"
           fc_core(na) = fc_core(na) + contrib
        end if

        if (iverbosity > 0) &
            write(stdout,'(5X,A,I3,2X,I1,''S -> '',I1,''S :'',F12.6)') &
                 atm(ityp(na)), na, n1, n2, contrib
      enddo
    enddo

  enddo  ! na

  write(stdout,*)
  deallocate( work, delta_v )

  call stop_clock('core_relax')

END SUBROUTINE hfi_fc_core_relax
