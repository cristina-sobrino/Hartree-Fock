subroutine CIS_matrix(nBas, E, ERI_MO, A_mat, E_CIS, no, nv)
        
        implicit none

        ! Input variables

         integer, intent(in)           :: nBas
         integer, intent(in)           :: no, nv
         double precision, intent(in)  :: E(nBas)
         double precision, intent(in)  :: ERI_MO(nBas, nBas, nBas, nBas)
        
        ! Output variables
         double precision, intent(out) :: A_mat(no*nv, no*nv)
         double precision, intent(out) :: E_CIS(no*nv)

         ! Local variables 
         integer                       :: i,b,j, a, k
         integer                       :: ia, jb
         integer                       :: scr2(no*nv, 2)
         double precision              :: contributions(no*nv)
         double precision              :: norm_val, normalized_value

         interface
        function Kronecker_delta(i,j) result(delta)
            integer, intent(in) :: i,j
            double precision    :: delta
        end function Kronecker_delta
    end interface


E_CIS = 0.d0

! allocate variables
!allocate(scr(2, no*nv))

k = 1

! First, we construct a matrix filled with the pairs of excitations
 
write(*,*) "The number of occupied orbitals is: ", no
write(*,*) "The number of virtual orbitals is: ", nv
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

   ! Extract the corresponding orbitals from the excitation matrix (scr2)

      i = scr2(ia,1) !occupied orbital
      j = scr2(jb,1) !occupied orbital
      a = scr2(ia,2) !virtual orbital
      b = scr2(jb,2) !virtual orbital

   ! Fill the A matrix
 
      A_mat(ia,jb) = ( E(a) - E(i) )*kronecker_delta(i,j) * &
      kronecker_delta(a,b) + 2.d0 * ERI_MO(i,b,a,j) - ERI_MO(i,b,j,a) 

   end do
end do

!


 call diagonalize_matrix(no*nv, A_mat, E_CIS)


 !write(*,*) "CIS eigenvectors"

 !do i = 1, no*nv

! write(*,'(*(f10.4))') A_mat(i,:)
! write(*,*) 
! write(*, '(*(a,f10.4))' )  "coeff", A_mat(i,:)**2
! write(*, '(*(a,f10.4))') "weight", A_mat(i,:)**2 / norm2(A_mat(i,:)**2)

 !end do


write(*,*)
write(*,*) "Excitation energies"
write(*,*) "-------------------"
write(*,*) 
 do i = 1, no*nv
 write(*,'(*(a,i3,a,f10.4,a,f10.4,a))') "Excited State", i,": ", E_CIS(i)*27.2114, " eV", E_CIS(i)*45.5640, " nm"

 write(*,*) "Printout of CI-coefficients larger than 0.01"
 write(*,*)
 write(*,'(T24, a, T20, a)') "Weight Coeff"
    do j = 1, no*nv
        norm_val = sum(A_mat(i, :)**2)
        
        if ( norm_val  > 1.0e-10 ) then  
           normalized_value = (A_mat(i,j)**2) / norm_val

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
