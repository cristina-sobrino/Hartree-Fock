subroutine orthogonalization_matrix(number_basis,s_mat,x_mat)

        implicit none
        
        ! input variables
        integer, intent(in)  :: number_basis
        real*8, intent(in)   :: s_mat(number_basis,number_basis)
        real*8, intent(out)  :: x_mat(number_basis,number_basis)

        ! local variable
        real*8 :: s_eigen(number_basis) !son los autovalores 
        real*8 :: s_diagonal(number_basis,number_basis) !la matriz U
        real*8 :: s_miniscula(number_basis, number_basis) !matriz diagonal con eigenvalues de S
        real*8 :: s_sqrt(number_basis,number_basis)
        real*8 :: test(number_basis, number_basis)
        integer :: i

        !/code/!
        
        ! x = s^1/2


        ! first we have to diagonalize S
        
        s_diagonal(:,:) = s_mat(:,:)


        call diagonalize_matrix(number_basis,s_diagonal,s_eigen) 
        
        s_miniscula(:,:) = 0.d0
        
        !copy of eigenvalues

        do i = 1, number_basis
           s_miniscula(i,i) = s_eigen(i)
        end do

        
        do i = 1,number_basis
           if (s_eigen(i) > 1e-12) then
                s_sqrt(i,i) =  1d0/sqrt(s_eigen(i))
           endif
        end do
        open(1, access='append', file = "matrix_print.out")
        

        write(1,*) "" 
        write(1,*) "The s_sqrt matrix is"
        write(1,*) ""

        do i = 1, number_basis
           write(1,'(*(f16.6))') s_sqrt(i,:)
        end do

        !close(1)

        ! Symmetric orthogonalization

        x_mat = matmul(s_diagonal, matmul(s_sqrt, transpose(s_diagonal)))
      
        ! Canonical orthogonalization
!        x_mat = matmul(s_diagonal,s_sqrt)

        ! TEST: XSX^T = I
        test = matmul(x_mat, matmul(s_mat, transpose(x_mat)))

        write(1,*) ""        
        write(1,*) "This should be the identity matrix"
         do i = 1, number_basis
           write(1,'(*(f16.6))') test(i,:)
        end do

end subroutine
        

