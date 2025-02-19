subroutine AO_to_MO(nBas,c,ERI_AO,ERI_MO)

! Expression of bi-electronic integrals in the MO basis set

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  integer                       :: p,q,r,s
  double precision              :: scr(nBas,nBas,nBas,nBas)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

! Memory allocation
  
!--------------------------------------
! AO to MO transformation starts here !
!--------------------------------------

! To be consistent with the previous code this integral is computed in physicist's notation (mu,la,nu,si) = (mu la | nu si ) = < mu nu | la si >
! and (p,q,r,s) = ( p q | r s ) = < p r | q s >

ERI_MO(:,:,:,:) = 0.d0 
 
 !do p = 1, nBas
 !   do q = 1, nBas
 !      do r = 1, nBas
 !         do s = 1, nBas
 !            do mu = 1, nBas
 !               do nu = 1, nBas
 !                  do la = 1, nBas
 !                     do si = 1, nBas
 !               
 !                      ERI_MO(p,r,q,s) = ERI_MO(p,r,q,s) + c(mu,p) * c(nu,q) * c(la,r) * c(si, s) * ERI_AO(mu,la,nu,si)
 !
 !                     end do
 !                  end do
 !               end do
 !            end do
 !         end do
 !      end do
 !   end do
 !end do
 !
!write(*,*) "Before transformation", ERI_MO(1,1,1,1)

! This code is very inefficient as having eight nested loops results in O(N)^8. To avoid this, we can do the transformation variable
! by variable so that we only have 5 nested loops -> O(N)^5. 


  scr(:,:,:,:) = 0.0d0
  ERI_MO(:,:,:,:) = 0.0d0
  
  ! Copy AO integrals to scratch space
  scr(:,:,:,:) = ERI_AO(:,:,:,:)
 
  ! First transformation (s -> si) 

  do mu = 1, nBas
     do nu = 1, nBas
        do la = 1, nBas 
           do s = 1, nBas
              do si = 1, nBas 

              ERI_MO(mu,la,nu,s) =  ERI_MO(mu,la,nu,s) + c(si,s) * scr(mu,la,nu,si)

              end do
           end do
        end do
     end do
   end do

   scr(:,:,:,:) = ERI_MO(:,:,:,:)

   ! Second transformation (nu -> q) 

   ERI_MO(:,:,:,:) = 0.d0

   do mu = 1, nBas
      do s = 1, nBas
         do la = 1, nBas
            do q = 1, nBas
               do nu = 1, nBas

                ERI_MO(mu,la,q,s) =  ERI_MO(mu,la,q,s) + c(nu,q) * scr(mu,la,nu,s)

               end do
            end do
         end do
      end do
   end do

  scr(:,:,:,:) = ERI_MO(:,:,:,:)

   ! Third transformation (la -> r) 

   ERI_MO(:,:,:,:) = 0.d0

   do mu = 1, nBas
      do q = 1, nBas
         do s = 1, nBas 
            do r = 1, nBas
               do la = 1, nBas

                ERI_MO(mu,r,q,s) =  ERI_MO(mu,r,q,s) + c(la,r) * scr(mu,la,q,s)

               end do
            end do
         end do
      end do
   end do

  scr(:,:,:,:) = ERI_MO(:,:,:,:)


   ! Fourth transformation (mu -> p) 

   ERI_MO(:,:,:,:) = 0.d0

   do s = 1, nBas
      do r = 1, nBas
         do q = 1, nBas
            do p = 1, nBas
               do mu = 1, nBas

                ERI_MO(p,r,q,s) =  ERI_MO(p,r,q,s) + c(mu,p) * scr(mu,r,q,s)

               end do
            end do
         end do
      end do
   end do


   write(*,*) "My first efficient code", ERI_MO(1,1,1,1)

  
end subroutine AO_to_MO
