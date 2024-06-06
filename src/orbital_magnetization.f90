!-----------------------------------------------------------------------
MODULE orbital_magnetization
  !-----------------------------------------------------------------------
  !
  ! ... the title is self explaining
  !
  USE kinds, ONLY : DP
  USE gipaw_module, ONLY: dudk_method
  USE lsda_mod,   ONLY: nspin, current_spin
  USE fft_base,      ONLY : dfftp, dffts 
  IMPLICIT NONE
  SAVE
  REAL(DP) :: orb_magn_LC(3)
  REAL(DP) :: orb_magn_IC(3)
  REAL(DP) :: berry_curvature(3)
  REAL(DP) :: orb_magn_tot(3)
  REAL(DP) :: delta_M_bare(3)
  REAL(DP) :: delta_M_para(3)
  REAL(DP) :: delta_M_dia(3)
  COMPLEX(dp), ALLOCATABLE :: dbecp(:,:,:)
  COMPLEX(dp), ALLOCATABLE :: paw_dbecp(:,:,:)
  INTEGER, ALLOCATABLE :: igk(:,:)
  real(dp) :: delta_k

  REAL(DP), ALLOCATABLE :: &
       dvrs(:,:,:)      ! gradient of (the total pot. in real space...)
   ! ... local variables
  REAL(DP), PARAMETER :: alpha = 1.0_dp / 137.03599911_dp
  ! g_e
  REAL(DP), PARAMETER :: g_e = 2.0023192778_dp
 ! g'
  REAL(DP), PARAMETER :: gprime = 2.0046385556_dp
  ! (alpha^2 g')/4
  REAL(DP), PARAMETER :: a2gp4 = alpha*alpha*gprime/4_dp
  ! (alpha^2 g')/8
  REAL(DP), PARAMETER :: a2gp8 = alpha*alpha*gprime/8_dp

 !-----------------------------------------------------------------------
END MODULE orbital_magnetization
!-----------------------------------------------------------------------


