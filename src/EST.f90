program EST

! Electronic Structure Theory program

  include 'parameters.h'

  integer                       :: nAt,nBas,nEl,nO,nV
  double precision              :: ENuc,EHF

  double precision,allocatable  :: ZNuc(:)
  double precision,allocatable  :: rAt(:,:)

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:)
  integer,allocatable           :: KShell(:)
  double precision,allocatable  :: CenterShell(:,:)
  double precision,allocatable  :: DShell(:,:)
  double precision,allocatable  :: ExpShell(:,:)

  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: T(:,:)
  double precision,allocatable  :: V(:,:)
  double precision,allocatable  :: Hc(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: ERI(:,:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: A_mat(:,:), E_TDHF(:)
  double precision,allocatable  :: e(:), E_CIS(:)
  double precision,allocatable  :: c(:,:)

  double precision              :: start_HF,end_HF,t_HF, start_TD_HF, end_TD_HF, t_TD_HF
  double precision              :: start_AOtoMO,end_AOtoMO,t_AOtoMO
  double precision              :: start_MP2,end_MP2,t_MP2
  double precision              :: start_CCD,end_CCD,t_CCD
  double precision              :: start_CIS_E, end_CIS_E, t_CIS_E

! Hello World

  write(*,*)
  write(*,*) '***************************'
  write(*,*) '* TCCM winter school 2023 *'
  write(*,*) '***************************'
  write(*,*)

!------------------------------------------------------------------------
! Read input information
!------------------------------------------------------------------------

! Read number of atoms, number of electrons of the system
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nBas = number of basis functions (see below)
!      = nO + nV

  call read_molecule(nAt,nEl,nO)
  allocate(ZNuc(nAt),rAt(nAt,3))

! Read geometry

  call read_geometry(nAt,ZNuc,rAt,ENuc)

  allocate(CenterShell(maxShell,3),TotAngMomShell(maxShell),KShell(maxShell), &
           DShell(maxShell,maxK),ExpShell(maxShell,maxK))

!------------------------------------------------------------------------
! Read basis set information
!------------------------------------------------------------------------

  call read_basis(nAt,rAt,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)

!------------------------------------------------------------------------
! Read one- and two-electron integrals
!------------------------------------------------------------------------

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas), &
           ERI(nBas,nBas,nBas,nBas),e(nBas),c(nBas,nBas),A_mat(no*nv, no*nv), &
           E_CIS(no*nv), E_TDHF(no*nv))

! Read integrals

  call read_integrals(nBas,S,T,V,Hc,ERI)  

!------------------------------------------------------------------------
! Orthogonalization X = S^(-1/2)
!------------------------------------------------------------------------
 
 open(1, file = "matrix_print.out")
 
 write(1,*) "the matrix S before diagonalization"
 write(1,*) ""

 do i = 1,nBas
    write(1,'(*(f16.6))') S(i,:)
 end do

 call orthogonalization_matrix(nBas,S,X)

 write(1,*) ""
 write(1,*) ""
 write(1,*) "The X matrix is"

 do i = 1, nBas 
         write(1,'(*(f16.6))') X(i,:)
 end do

 close(1)
 ! get the matrices X

!------------------------------------------------------------------------
! Compute restricted HF energy
!------------------------------------------------------------------------

    call cpu_time(start_HF)
    call RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)
    call cpu_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! AO to MO transformation
!------------------------------------------------------------------------

    allocate(ERI_MO(nBas,nBas,nBas,nBas))

    call cpu_time(start_AOtoMO)

!    write(*,*) "start of the transformation"
    call AO_to_MO(nBas,c,ERI,ERI_MO)
    write(*,*) "the new integral is", ERI_MO(1,1,1,1)
    open(25, file="MOmatrix") 
    write(25, *) ERI_MO 
    close(25)
    call cpu_time(end_AOtoMO)

    t_AOtoMO = end_AOtoMO - start_AOtoMO
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! Compute CIS energy
!------------------------------------------------------------------------
    call cpu_time(start_CIS_E)

    write(*,*) "------------------------"
    write(*,*) "Start of CIS calculation"
    write(*,*) "------------------------"

    call CIS_matrix(nBas, E, ERI_MO, A_mat, E_CIS, no, nv) 
!    write(*,'(A, F12.6, A)') "Energy diff S1/S0", E_CIS(1) * 27.2114, "eV"

    t_CIS_E = end_CIS_E - start_CIS_E 
    call cpu_time(end_CIS_E)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CIS energy calculation = "', t_CIS_E, 'seconds'
    write(*,*)
   
!------------------------------------------------------------------------
! Compute TDHF energy
!------------------------------------------------------------------------
    call cpu_time(start_TD_HF)

    write(*,*) "-------------------------"
    write(*,*) "Start of TDHF calculation"
    write(*,*) "-------------------------"

    call TD_HF(nBas, E, ERI_MO, A_mat, E_TDHF, no, nv)
!    write(*,'(A, F12.6, A)') "Energy diff S1/S0", E_CIS(1) * 27.2114, "eV"

    t_TD_HF = end_TD_HF - start_TD_HF
    call cpu_time(end_TD_HF)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CIS energy calculation = "', t_CIS_E, 'seconds'
    write(*,*)

!------------------------------------------------------------------------
! Compute MP2 energy
!------------------------------------------------------------------------

    call cpu_time(start_MP2)
!   call MP2(nBas,nO,nV,e,ERI_MO,ENuc,EHF)
    call cpu_time(end_MP2)

    t_MP2 = end_MP2 - start_MP2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 = ',t_MP2,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! Compute CCD energy
!------------------------------------------------------------------------

    call cpu_time(start_CCD)
!   call CCD(nBas,nO,nV,ENuc,EHF,e,ERI_MO)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCD = ',t_CCD,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! End of EST
!------------------------------------------------------------------------
end program EST
