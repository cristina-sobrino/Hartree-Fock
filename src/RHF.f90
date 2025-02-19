subroutine RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)

! Perform a restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

! Local variables

  integer,parameter             :: maxSCF = 64
  double precision,parameter    :: thresh = 1d-5
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: Gap
  double precision              :: ET,EV,EJ
  double precision              :: EK
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: F(:,:),Fp(:,:)
  double precision,allocatable  :: error(:,:)
  double precision,external     :: trace_matrix
  integer                       :: mu
  integer                       :: nu 
  integer                       :: lambda 
  integer                       :: sigma
  double precision              :: G(nBas,nBas)
  integer                       :: i

! Output variables

  double precision,intent(out)  :: EHF
  double precision,intent(out)  :: e(nBas)
  double precision,intent(out)  :: c(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Restricted Hartree-Fock calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(cp(nBas,nBas),P(nBas,nBas),      &
           J(nBas,nBas),K(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           error(nBas,nBas))

! Guess coefficients and eigenvalues

  F(:,:) = Hc(:,:)

! Initialization

  nSCF = 0
  Conv = 1d0
  P(:,:) = 0.d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| RHF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','HF energy','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

  ! 4.A Calculate G 

   G(:,:) = 0.d0

   ! do mu = 1, nBas
   !    do nu = 1, nBas
   !       do lambda = 1, nBas
   !           do sigma = 1, nBas 
                
   !           G(mu,nu) = G(mu,nu) + P(lambda,sigma) * (eri(mu,nu,sigma,lambda) - 0.5d0 * eri(mu,lambda,sigma,nu))

    !          end do
    !      end do
    !   end do
    ! end do

    ! F(:,:) = Hc(:,:) + G(:,:)
    
    Fp = matmul(transpose(X),matmul(F,X))

   ! Diagonalize Fp to obtain Cp and e
    call diagonalize_matrix(nbas,Fp,e)

    ! Cp is composed of the eigenvectors of the diagonalized Fp
    Cp(:,:) = Fp(:,:)
     
    ! Compute C = X @ Cp
    C = matmul(X,Cp)


    ! Gap = E(LUMO) - E(HOMO)

     Gap = abs(E(no) - E(no+1))

    ! Conv = FPS - SPF and then take the maximum absolute value of this matrix
    ! then recalculate J and K and 

!   Increment 
    nSCF = nSCF + 1


      ! Reinitialize P
    P = 0.d0
    G = 0.d0
    ! New P matrix

        
     P(:,:) = 2.d0 * matmul( C(:,1:no), transpose( C(:, 1:no)))


     K(:,:) = 0.d0
     J(:,:) = 0.d0
     do mu = 1, nBas
        do nu = 1, nBas
            do lambda = 1, nBas
                do sigma = 1, nBas

                G(mu,nu) = G(mu,nu) + P(lambda,sigma) * (eri(mu,lambda,nu,sigma) - 0.5d0 * eri(mu,lambda,sigma,nu))
                K(mu,nu) = -0.5d0 * P(lambda,sigma) * eri(mu,lambda,sigma,nu) + K(mu,nu)
                J(mu,nu) =  P(lambda,sigma) * eri(mu,lambda,nu,sigma) + j(mu,nu)
             
              end do
          end do
       end do
     end do

    

    F = Hc + G

   !  EHF = 1/2 * tr(P * (H + F))
      
    EHF = 0.5d0 * trace_matrix(nbas,(matmul(P, Hc + F)))

    

    ! Commutator to evaluate the convergency

    Conv = maxval(abs(matmul(F, matmul(P,S)) - matmul(S, matmul(P,F))))


!   Dump results
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',EHF+ENuc,'|',Conv,'|',Gap,'|'
 
  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif


    ET = trace_matrix(nbas,matmul(P,T))

    EV = trace_matrix(nBas, matmul(p,V))

    EJ = trace_matrix(nBas, matmul(P,J)) * 0.5d0

    EK = trace_matrix(nBas,matmul(p, K)) * 0.5d0


! Compute final HF energy

  call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)

end subroutine RHF
