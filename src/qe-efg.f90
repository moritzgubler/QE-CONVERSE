!
! Copyright (C) 2001-2009 GIPAW and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM qe_efg
  !-----------------------------------------------------------------------
  !
  ! Standalone driver to compute EFG tensors and NMR/NQR quadrupolar
  ! coupling constants from a ground-state QE charge density.
  !
  ! Workflow: read_file -> gipaw_setup -> calc_efg -> print_efg_summary
  !
  ! Input namelist: &input_qeefg
  !   prefix    : must match the pw.x SCF prefix
  !   outdir    : must match the pw.x SCF outdir
  !   q_efg(n)  : nuclear quadrupole moment (barn) for atom type n
  !               (indexed as in ATOMIC_SPECIES block of pw.x input)
  !               omit or set to 0 to skip Cq for that type
  !
  USE gipaw_module,  ONLY : q_efg
  USE constants,     ONLY : eps12
  USE environment,   ONLY : environment_start, environment_end
  USE io_files,      ONLY : prefix, tmp_dir
  USE io_global,     ONLY : ionode, ionode_id
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : world_comm
  USE mp_global,     ONLY : mp_startup
  USE mp_images,     ONLY : my_image_id
  USE control_flags, ONLY : io_level
  USE cellmd,        ONLY : cell_factor

  IMPLICIT NONE

  character(len=9)  :: codename = 'QE-EFG'
  character(len=256) :: outdir
  character(len=256), external :: trimcheck
  integer :: ios

  NAMELIST / input_qeefg / prefix, outdir, q_efg

#if defined(__MPI)
  call mp_startup(start_images=.true., images_only=.false.)
#else
  call mp_startup(start_images=.false.)
#endif
  call environment_start(codename)

  if (.not. ionode .or. my_image_id > 0) goto 100

  call input_from_file()

  ! defaults
  prefix = 'prefix'
  CALL get_environment_variable('ESPRESSO_TMPDIR', outdir)
  IF (TRIM(outdir) == ' ') outdir = './'
  q_efg(:) = 0.d0

  read(5, input_qeefg, iostat=ios)
  tmp_dir = trimcheck(outdir)

100 continue

#ifdef __MPI
  CALL mp_bcast(tmp_dir, 0, world_comm)
  CALL mp_bcast(prefix,  0, world_comm)
  CALL mp_bcast(q_efg,   0, world_comm)
#endif

  io_level   = 1
  cell_factor = 1.1d0

  call read_file()
  call gipaw_setup()
  call calc_efg()
  call print_efg_summary()

  call environment_end(codename)
  call stop_code(.true.)

  STOP

END PROGRAM qe_efg
