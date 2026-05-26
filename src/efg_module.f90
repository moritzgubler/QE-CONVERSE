!-----------------------------------------------------------------------
MODULE efg_mod
  !-----------------------------------------------------------------------
  ! Stores the symmetrized EFG tensor for later summary printing.
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  SAVE

  real(dp), allocatable :: efg_tensor(:,:,:)   ! (3,3,nat)

END MODULE efg_mod
