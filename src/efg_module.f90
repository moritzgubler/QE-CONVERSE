!-----------------------------------------------------------------------
MODULE efg_mod
  !-----------------------------------------------------------------------
  ! EFG parameters and results.
  !
  USE kinds,      ONLY : dp
  USE parameters, ONLY : ntypx
  IMPLICIT NONE
  SAVE

  ! nuclear quadrupole moment Q for the EFG (input, per atom type)
  real(dp) :: q_efg ( ntypx )

  ! nuclear spin I for the EFG (quadrupolar frequency nu_Q)
  real(dp) :: i_efg ( ntypx )

  ! symmetrized EFG tensor, stored for later summary printing
  real(dp), allocatable :: efg_tensor(:,:,:)   ! (3,3,nat)

END MODULE efg_mod
