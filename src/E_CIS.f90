! subroutine E_CIS(A_mat, E_CIS)
! 
!         implicit none
! 
!         ! input variables 
! 
!         double precision,intent(in)                   :: A_mat(nBas, nBas, nBas, nBas) 
!         integer, intent(in)                           :: nBas
!         double precision, intent(in)                  :: ERI_MO(nBas, nBas, nBas, nBas)
! 
!         ! Local variables 
! 
!         ! Output variables
! 
!         double precision, intent(out)                 :: E_CIS(nBas)
! 
! 
! call diagonalize_matrix(nBas, A_mat, E_CIS)
! 
! 
! 
! end subroutine
