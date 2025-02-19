subroutine TD_HF(nBas, E, ERI_MO, A_mat, E_TDHF, no, nv)
        
        implicit none

        ! Input variables

         integer, intent(in)           :: nBas
         integer, intent(in)           :: no, nv
         double precision, intent(in)  :: E(nBas)
         double precision, intent(in)  :: ERI_MO(nBas, nBas, nBas, nBas)
        
        ! Output variables
         double precision, intent(out) :: A_mat(no*nv, no*nv)
         double precision              :: B_mat(no*nv, no*nv)
         double precision              :: C_mat(no*nv, no*nv)
         double precision, intent(out) :: E_TDHF(no*nv)

         ! Local variables 
         integer                       :: i,b,j, a, k
         integer                       :: ia, jb
         integer                       :: scr2(no*nv, 2)
         double precision              :: ab_sqrt(no*nv, no*nv)
         double precision              :: a_plus_b(no*nv, no*nv)
         double precision              :: a_minus_b(no*nv, no*nv)
         double precision              :: omega_square(no*nv) 
         double precision              :: eigen_sqrt(no*nv, no*nv)
         double precision              :: eigen(no*nv)
         double precision              :: x_mat(no*nv, no*nv)
         double precision              :: norm_val, normalized_value

          double precision              :: U(no*nv, no*nv)
    double precision              :: Vt(no*nv, no*nv)
    double precision              :: D(no*nv)
    double precision              :: D_sqrt(no*nv, no*nv)

         interface
        function Kronecker_delta(i,j) result(delta)
            integer, intent(in) :: i,j
            double precision    :: delta
        end function Kronecker_delta
    end interface


! = 0.d0

! allocate variables
!allocate(scr(2, no*nv))

k = 1

! First, we construct a matrix filled with the pairs of excitations
 
   do i = no, 1, -1
     do j = no + 1, nv + no

      scr2(k, 1) = i
      scr2(k, 2) = j
        
      k = k + 1

      end do
   end do


A_mat = 0.d0        

do ia = 1, no*nv
   do jb = 1, no*nv

      i = scr2(ia,1)
      j = scr2(jb,1)
      a = scr2(ia,2)
      b = scr2(jb,2)


      A_mat(ia,jb) = ( E(a) - E(i) )*kronecker_delta(i,j) * &
      kronecker_delta(a,b) + 2.d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a) 

   end do
end do



B_mat = 0.d0

do ia = 1, no*nv
   do jb = 1, no*nv

      i = scr2(ia,1)
      j = scr2(jb,1)
      a = scr2(ia,2)
      b = scr2(jb,2)

      B_mat(ia,jb) = 2.d0 * ERI_MO(i,j,a,b) - ERI_MO(i,j,b,a)

   end do
end do

! To compute the square root of the A-B matrix, we need to perform a symmetric orthogonalization

 write(*,*) "A-B matrix. This is a non-diagonal matrix" 
 
 do i = 1, no*nv
    write(*,'(*(f6.3))') A_mat(i,:) - B_mat(i,:)
 end do
 
 ! First, we diagonalize the matrix A-B
 
a_minus_b = A_mat - B_mat
 
! call diagonalize_matrix(no*nv, a_minus_b, eigen)
 
 ! We compute the square root of the eigenvalues
 
! eigen_sqrt = 0.d0
 
! do i = 1, no*nv
!    eigen_sqrt(i,i) = sqrt(eigen(i))
! end do
 
!  write(*,*) "Diagonalized A-B matrix. This matrix should be diagonal"
!  
!  do i = 1, no*nv
!          write(*, '(*(f6.3))') a_minus_b(i,:)
!  end do
!  
!  
!  write(*,*) "Eigenvalues of A-B matrix. This matrix should be diagonal"
!  
!  do i = 1, no*nv
!     write(*, '(*(f6.3))') eigen_sqrt(i,:)
!  end do
!  
 ! We compute the multiplication (A-B) @ eigen_sqrt @ (A-B)†
 
! ab_sqrt = matmul( a_minus_b ,matmul(eigen_sqrt, transpose(a_minus_b)))
 
 
!  write(*,*)
!  write(*,*) "sqrt of the A-B matrix: "
!  write(*,*)
!  do i = 1, no*nv
!     write(*, '(*(f6.3))') ab_sqrt(i,:)
!  end do
 

call svd(no*nv, a_minus_b, U, D, Vt)

D_sqrt = 0.d0
    do i = 1, no*nv
        D_sqrt(i,i) = sqrt(D(i))
    end do

 ! Construct matrix square root: ab_sqrt = U * sqrt(D) * Vt
 ab_sqrt = matmul(U, matmul(D_sqrt, Vt))

a_plus_b = A_mat + B_mat

C_mat = matmul(ab_sqrt, matmul(a_plus_b,ab_sqrt))


!write(*,*) "The C matrix is"

!do i = 1, no*nv
!   write(*, '(*(f6.3))') C_mat(i,:)
!end do

! C_mat diagonalization

 call diagonalize_matrix(no*nv, C_mat, omega_square)

 E_TDHF(:) = 0.d0

 do i = 1, no*nv 
    E_TDHF(i) = sqrt(omega_square(i))
 end do
 
!write(*,*) "The first excitation energy is:"
!write(*, '(*(f6.3, a))') E_TDHF(1)*27.2114, "eV"


write(*,*)
write(*,*) " TD HF Excitation energies"
write(*,*) "--------------------------"
write(*,*)
 do i = 1, no*nv
 write(*,'(*(a,i3,a,f10.4,a,f10.4,a))') "Excited State", i,": ", E_TDHF(i)*27.2114, " eV", E_TDHF(i)*45.5640, " nm"

 write(*,*) "Printout of CI-coefficients larger than 0.01"
 write(*,*)
 write(*,'(T24, a, T20, a)') "Weight Coeff"
    do j = 1, no*nv
        norm_val = sum(C_mat(i, :)**2)

        if ( norm_val  > 1.0e-10 ) then
           normalized_value = (C_mat(i,j)**2) / norm_val

            if (normalized_value > 0.01 ) then

                b = scr2(j, 1)
                a = scr2(j, 2)
                write(*, '(I4," ->",I4,T20,F9.5, F9.5)') &
                b, a, normalized_value, normalized_value ** 2

                write(*,*)

            end if
        end if
    end do
  end do






end subroutine
