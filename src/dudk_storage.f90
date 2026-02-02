!-----------------------------------------------------------------------
MODULE dudk_storage
  !-----------------------------------------------------------------------
  !
  ! Module for in-memory storage of wavefunction derivatives (dudk)
  ! This provides an alternative to disk I/O for better performance
  ! when memory is available.
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  SAVE

  ! Control flag: if .true., store dudk in memory instead of disk
  LOGICAL :: dudk_in_memory = .false.

  ! Storage array: dudk_mem(npwx, nbnd, 3_directions, nks)
  ! Only allocated if dudk_in_memory = .true.
  COMPLEX(DP), ALLOCATABLE :: dudk_mem(:,:,:,:)

  ! Flag to track if memory has been allocated
  LOGICAL, PRIVATE :: dudk_mem_allocated = .false.

!-----------------------------------------------------------------------
  CONTAINS
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE allocate_dudk_storage(npwx, nbnd, nks)
  !-----------------------------------------------------------------------
    !
    ! Allocate memory for in-memory dudk storage
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npwx  ! max number of plane waves
    INTEGER, INTENT(IN) :: nbnd  ! number of bands
    INTEGER, INTENT(IN) :: nks   ! number of k-points

    IF (.NOT. dudk_in_memory) RETURN

    IF (dudk_mem_allocated) THEN
      CALL errore('allocate_dudk_storage', 'dudk_mem already allocated', 1)
    END IF

    ALLOCATE(dudk_mem(npwx, nbnd, 3, nks))
    ! Note: No initialization - array will be filled by store_dudk
    dudk_mem_allocated = .true.

  END SUBROUTINE allocate_dudk_storage

  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_dudk_storage()
  !-----------------------------------------------------------------------
    !
    ! Deallocate in-memory dudk storage
    !
    IMPLICIT NONE

    IF (dudk_mem_allocated) THEN
      DEALLOCATE(dudk_mem)
      dudk_mem_allocated = .false.
    END IF

  END SUBROUTINE deallocate_dudk_storage

  !-----------------------------------------------------------------------
  SUBROUTINE store_dudk(dudk, ik, ipol, npw, nbnd_loc)
  !-----------------------------------------------------------------------
    !
    ! Store dudk for a given k-point and polarization direction
    !
    IMPLICIT NONE
    COMPLEX(DP), INTENT(IN) :: dudk(:,:)  ! dudk array (npwx, nbnd)
    INTEGER, INTENT(IN) :: ik             ! k-point index
    INTEGER, INTENT(IN) :: ipol           ! polarization direction (1-3)
    INTEGER, INTENT(IN) :: npw            ! number of plane waves for this k
    INTEGER, INTENT(IN) :: nbnd_loc       ! number of bands

    IF (.NOT. dudk_in_memory) RETURN

    IF (.NOT. dudk_mem_allocated) THEN
      CALL errore('store_dudk', 'dudk_mem not allocated', 1)
    END IF

    IF (ipol < 1 .OR. ipol > 3) THEN
      CALL errore('store_dudk', 'invalid ipol value', ipol)
    END IF

    ! Store the data
    dudk_mem(1:npw, 1:nbnd_loc, ipol, ik) = dudk(1:npw, 1:nbnd_loc)

  END SUBROUTINE store_dudk

  !-----------------------------------------------------------------------
  SUBROUTINE retrieve_dudk(dudk, ik, ipol, npw, nbnd_loc)
  !-----------------------------------------------------------------------
    !
    ! Retrieve dudk for a given k-point and polarization direction
    !
    IMPLICIT NONE
    COMPLEX(DP), INTENT(OUT) :: dudk(:,:) ! dudk array (npwx, nbnd)
    INTEGER, INTENT(IN) :: ik              ! k-point index
    INTEGER, INTENT(IN) :: ipol            ! polarization direction (1-3)
    INTEGER, INTENT(IN) :: npw             ! number of plane waves for this k
    INTEGER, INTENT(IN) :: nbnd_loc        ! number of bands

    IF (.NOT. dudk_in_memory) RETURN

    IF (.NOT. dudk_mem_allocated) THEN
      CALL errore('retrieve_dudk', 'dudk_mem not allocated', 1)
    END IF

    IF (ipol < 1 .OR. ipol > 3) THEN
      CALL errore('retrieve_dudk', 'invalid ipol value', ipol)
    END IF

    ! Retrieve the data (no initialization needed - we copy what's needed)
    dudk(1:npw, 1:nbnd_loc) = dudk_mem(1:npw, 1:nbnd_loc, ipol, ik)

  END SUBROUTINE retrieve_dudk

  !-----------------------------------------------------------------------
  FUNCTION get_dudk_memory_mb(npwx, nbnd, nks) RESULT(mem_mb)
  !-----------------------------------------------------------------------
    !
    ! Calculate memory requirement in MB for dudk storage
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npwx, nbnd, nks
    REAL(DP) :: mem_mb
    REAL(DP) :: bytes_per_complex

    ! COMPLEX(DP) = 2 * 8 bytes = 16 bytes
    bytes_per_complex = 16.0_DP

    ! Calculate: npwx * nbnd * 3 * nks * bytes_per_complex / (1024^2)
    mem_mb = DBLE(npwx) * DBLE(nbnd) * 3.0_DP * DBLE(nks) * &
             bytes_per_complex / (1024.0_DP * 1024.0_DP)

  END FUNCTION get_dudk_memory_mb

!-----------------------------------------------------------------------
END MODULE dudk_storage
!-----------------------------------------------------------------------
