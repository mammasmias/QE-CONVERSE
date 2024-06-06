!--------------------------------------------------------------------
subroutine set_dvrs (dvrs, vrs, nrxxs, nspin)
  !--------------------------------------------------------------------
  ! gradient of the total local potential vrs on the smooth mesh
  !
  USE kinds
  USE gipaw_module,      ONLY : lambda_so
  USE cell_base,         ONLY : tpiba
  USE gvect,             ONLY : ngm, g
  USE gvecs,             ONLY : ngms
  USE fft_base,          ONLY : dffts, dfftp, dfftb
  implicit none

  integer :: nspin, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  real(dp), intent(in) ::  vrs (dffts%nnr, nspin)
  real(dp), intent(inout) :: dvrs ( dffts%nnr, nspin, 3)
  real(dp), allocatable :: aux(:,:)
  integer:: is, ig, ipol


  nr1s=dffts%nr1
  nr2s=dffts%nr2
  nr3s=dffts%nr3
  nrx1s=dffts%nr1x
  nrx2s=dffts%nr2x
  nrx3s=dffts%nr3x
  

  if (all(lambda_so == 0.d0)) return

  allocate( aux(3,nrxxs) )
  do is = 1, nspin
    call gradient( nrx1s, nrx2s, nrx3s, nr1s, nr2s, nr3s, nrxxs, &
                   vrs(1,is), ngm, g, aux )
    do ipol = 1, 3
     dvrs(1:nrxxs,is,ipol) = aux(ipol,1:nrxxs)
    enddo
  enddo
deallocate (aux)
!
end subroutine set_dvrs

!----------------------------------------------------------------------------
SUBROUTINE gradient( nrx1, nrx2, nrx3, &
                     nr1, nr2, nr3, nrxx, a, ngm, g, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE fft_base,          ONLY : dffts, dfftp, dfftb
!  USE gvect,     ONLY : nlm
!  USE wvfct,     ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx
  INTEGER,  INTENT(IN)  :: ngm
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga(3,nrxx)
  !
  INTEGER                  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
!  print*, 'entra in gradient'
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 )
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Rho',aux, dfftp)
!  CALL cft3( aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  ga(:,:) = 0.D0
  !
  DO ipol = 1, 3
     !
     gaux(:) = 0.D0
     !

    gaux(dfftp%nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(dfftp%nl(:)) ), REAL( aux(dfftp%nl(:)) ) )
     !
!     IF ( gamma_only ) THEN
        !
!        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) )
        !
!     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
!     CALL cft3( gaux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  CALL invfft ('Rho',gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = ga(ipol,:) + tpiba * DBLE( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE gradient
!
!------------------------------------------
