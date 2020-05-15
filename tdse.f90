!!!!!!!!!!!!!!!!!!      One Gaussian Basis Function         !!!!!!!!!!!!!!
  ! Hartree-Fock SCF code: Restricted Hartree Fock Code using 6-31++g** basis set
  !

  ! gfortran -Wall -fbounds-check /work/klopata/repo/fortran/matutils.f90 /usr/lib/lapack/liblapack.so tdhf_01.f90

  ! xxx some change

  ! gfortran -I /work/adonay/repo/klopata/library -c -fbounds-check mag_prop_tdhf.f90
  ! gfortran mag_prop_tdhf.o /work/adonay/repo/klopata/library/matutils.o /usr/lib/lapack/liblapack.so
  !
module hf_module
  use matutils

  implicit none

  ! Constants (dont change through SCF or time-propagation)
  integer, parameter :: n = 30                      ! number of basis functions
  double precision   :: Hcore(n,n)                 ! core (1e) part of Fock matrix
  double precision   :: S(n,n)                 ! overlap (n x n)
  double precision   :: Dx(n,n), Dy(n,n), Dz(n,n)  ! transition dipole matrices (3D)
  double precision   :: X(n,n), Y(n,n), svals(n), svecs(n,n), stemp(n,n)             ! transform mats
  double precision   :: twoe(n,n,n,n)              ! 4center 2e integrals
  double complex     :: zi = (0d0, 1d0)            ! imaginary "i"
  !  double complex     :: Tr, k1(n,n), k2(n,n),  k3(n,n), k4(n,n), P1mo(n,n), P2mo(n,n), P3mo(n,n)
  !  double precision   :: t                         ! current time
  double precision   :: dt = 0.001d0                ! time step
  integer, parameter :: nt = 50000                 ! number of time steps
  double precision, Parameter :: Pi = 3.1415927d0
  double complex     :: n_occ(n,n), p_gamma(n,n)
  !  integer :: it

contains

  ! Read in integrals from file
  subroutine init_integrals ()
    implicit none

    integer, parameter :: fnum = 19823

    double precision :: T(n,n), V(n,n)  ! kinetic and potential (will be stored in Hcore eventually)
    integer :: i, j, k, l
    integer :: i_in, j_in, k_in, l_in
    double precision :: val
    character(len=20) :: tag
    character(len=10) :: e1, g1, e2, g2
    !
    ! Read in various n x n matrices
    !
    open (unit=fnum, file="overlap_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) tag, i_in, e1, g1, j_in, e2, g2, val

       if (tag .ne. "1eov") then
          write (6,*) "Did not find 1eov tag in input file"
          stop
       endif

       S(i_in, j_in) = val
    enddo
    close (fnum)

    open (unit=fnum, file="kinetic_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) tag, i_in, e1, g1, j_in, e2, g2, val

       if (tag .ne. "1eke") then
          write (6,*) "Did not find 1eke tag in input file"
          stop
       endif

       T(i_in, j_in) = val
    enddo
    close (fnum)

    open (unit=fnum, file="potential_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) tag, i_in, e1, g1, j_in, e2, g2, val

       if (tag .ne. "1epe") then
          write (6,*) "Did not find 1epe tag in input file"
          stop
       endif

       V(i_in, j_in) = val
    enddo
    close (fnum)

    !
    ! Read in two electron integrals
    !
    open (unit=fnum, file="2E_6_31++g_h2o.dat", status="old")
    do i = 1, n*n*n*n
       read (fnum, *) tag, i_in, j_in, k_in, l_in, val

       if (tag .ne. "ao") then
          write (6,*) "Did not find ao tag in input file"
          stop
       endif

       twoe(i_in, j_in, k_in, l_in) = val
    end do
    close (fnum)



    Hcore = T + V


    !
    ! Read occ. matrix
    !
    open (unit=fnum, file="occu_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) i_in, e1, g1, j_in, e2, g2, val
       n_occ(i_in, j_in) = val
    enddo
    close (fnum)

    !
    ! Read imaginary potential matrix
    !
    open (unit=fnum, file="p_gamma.dat", status="old")
    do i = 1, n*n
       read (fnum, *) i_in, j_in, val

       p_gamma(i_in, j_in) = val

       ! if (p_gamma(i_in, j_in) .ne. p_gamma(j_in, i_in)) then
       !    write (6,*) "ABC potential is not symmetric"
       !    stop
       ! endif


       ! p_gamma(i_in, j_in) = val
    enddo
    close (fnum)


    !
    ! Print input n x n matrices to screen
    !
    write (6, *) " Input N x N matrices: "
    write (6,*) "  i   j              S                   T                   V                   Hcore         "
    write (6,*) "-----------------------------------------------------------------------------------------------"
    do i = 1, n
       do j = 1, n
          write (6,"(i4, i4, 5f20.4)") i, j, S(i,j), T(i,j), V(i,j), Hcore(i,j)
       enddo
    enddo
    write (6, *) ""
    write (6, *) ""


    !
    ! Print 2e integrals
    !
    ! write (6, *) " Input twoe integrals: "
    ! write (6, *) "   i             j             k            l          val"
    ! write (6, *) "------------------------------------------------------------------------"
    ! do i = 1, n
    !    do j = 1, n
    !       do k = 1, n
    !          do l = 1, n
    !             write (6,*) i, j, k, l, twoe(i, j, k, l)
    !          enddo
    !       enddo
    !    end do
    ! end do
    ! write (6, *) ""
    ! write (6, *) ""


    !
    ! Transition dipole matrix
    !
    ! XXX HARDCODED FOR NOWXXXXX NONONONONONONO
    !


    ! Dipole matrix directed in the 'X' direction is read in from NWChem

    ! open (unit=fnum, file="Dx.dat", status="old")
    ! do i = 1, n*n
    !    read (fnum, *)  i_in, j_in,val
    !    Dx(i_in, j_in) = val
    ! end do
    ! close (fnum)

    ! Dipole matrix directed in the 'Y' direction is read in from NWChem

    ! open (unit=fnum, file="Dy.dat", status="old")
    ! do i = 1, n*n
    !    read (fnum, *) i_in, j_in,val
    !    Dy(i_in, j_in) = val
    ! end do
    ! close (fnum)

    ! Dipole matrix directed in the 'Z' direction is read in from NWChem

    open (unit=fnum, file="Dz_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) i_in, j_in, val
       Dz(i_in, j_in) = val
    end do
    close (fnum)


  end subroutine init_integrals


  !
  ! Initialize the n x n density matrix (in AO basis)
  !
  subroutine init_densmat (Pao)
    implicit none
    integer :: i, j
    integer :: i_in, j_in
    double precision :: val
    double complex, intent(out)  :: Pao(n,n)   ! density matrix
    integer, parameter :: fnum = 19823

    !  Dummy density matrix
    !  XXX ADONAY INPUT FROM NWCHEM
    !    Pao(1,1) = 0.5445200870d0
    !   Pao(1,2) = 0.5445200870d0
    !   Pao(2,1) = 0.5445200870d0
    !   Pao(2,2) = 0.5445200870d0

    ! Density matrix is read in from NWChem
    open (unit=fnum, file="density_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) i_in, j_in, val
       Pao(i_in, j_in) = val
    end do
    close (fnum)


  end subroutine init_densmat

  subroutine init_transform_new ()
    implicit none
    integer :: i, j
    integer :: i_in, j_in
    double precision :: val
    integer, parameter :: fnum = 19823

    ! Canonical transform matrix is read in from NWChem

    open (unit=fnum, file="transform_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) i_in, j_in, val
       X(i_in, j_in) = val
    end do
    close (fnum)

    !   X(1,1)= 0.1748634146d+01
    !   X(1,2)= 0.5217854382d+00
    !   X(2,1)= -0.1748634146d+01
    !   X(2,2)= 0.5217854382d+00

    ! X(1,1)= 0.5217854382d+00
    ! X(1,2)= 0.1748634146d+01
    ! X(2,1)= 0.5217854382d+00
    ! X(2,2)= 0.1748634146d+01

    ! svals (1) = 0.1635204581d0
    ! svals (2) = 1.836479542d0

    ! Eigenvalues of the overlap  matrix after diagonalizing are read in from NWChem

    open (unit=fnum, file="eigen_6_31++g_h2o.dat", status="old")
    do i = 1, n
       read (fnum, *) i_in, val
       svals(i_in) = val
    end do
    close (fnum)


    stemp = 0d0
    do i = 1, n
       stemp(i,i) = svals(i)
    enddo

    ! Y = X s = U s^(1/2)
    Y = matmul(X, stemp)


    write (6, *) " X Transform matrix"
    write (6,*) "  i   j              X     "
    write (6,*) "---------------------------"
    do i = 1, n
       do j = 1, n
          write (6,"(i4, i4, 5f20.4)") i, j, dble(X(i,j))
       enddo
    enddo
    write (6, *) ""
    write (6, *) ""

    write (6, *) " Y Transform matrix"
    write (6,*) "  i   j              Y     "
    write (6,*) "---------------------------"
    do i = 1, n
       do j = 1, n
          write (6,"(i4, i4, 5f20.4)") i, j, dble(Y(i,j))
       enddo
    enddo
    write (6, *) ""
    write (6, *) ""


  end subroutine init_transform_new

  subroutine init_transform ()
    implicit none
    !XXX ADONAY:  X = U s^-1/2 or X = U * S^-1/2 * U --------> But we get two different answers! The first method doesn't give a diagonally
    !symmetrical X matrix. With the second method, the minus signs do not match up when comparing eigenvector matrix obtained from octave to svecs obtained
    ! straight from nwchem.

    integer :: i , j
    integer :: i_in, j_in
    double precision :: val
    integer, parameter :: fnum = 19823

    !OLD WAY
    ! svecs (1,1) = 0.7071067812E+00
    ! svecs (1,2) = 0.7071067812E+00
    ! svecs (2,1) = -0.7071067812E+00
    ! svecs (2,2) =  0.7071067812E+00

    !OLD WAY
    ! svecs (1,1) = -0.7071067812E+00
    ! svecs (1,2) = 0.7071067812E+00
    ! svecs (2,1) = 0.7071067812E+00
    ! svecs (2,2) =  0.7071067812E+00

    ! Eigenvectors of the overlap  matrix after diagonalizing are read in from NWChem

    open (unit=fnum, file="evecs_6_31++g_h2o.dat", status="old")
    do i = 1, n*n
       read (fnum, *) i_in, j_in, val
       svecs(i_in, j_in) = val
    end do
    close (fnum)


    ! svals (1) = 0.1635204581d0
    ! svals (2) = 1.836479542d0

    !
    ! generate matrix: s^{-1/2}
    !
    stemp = 0d0
    do i = 1, n
       stemp(i,i) = 1d0/sqrt(svals(i))
    enddo

    ! X = U s^(-1/2)
    X = matmul(svecs, stemp)


    !
    ! generate s matrix: s
    !
    stemp = 0d0
    do i = 1, n
       stemp(i,i) = svals(i)
    enddo

    ! Y = X s = U s^(1/2)
    Y = matmul(X, stemp)

    do i = 1, n
       do j = 1, n
          write (6, *) "svecs: ", i, j, svecs(i,j)
       enddo
    enddo
    write (6,*) svecs

    do i = 1, n
       do j = 1, n
          write (6, *) "Y: ", i, j, Y(i,j)
       enddo
    enddo
    write (6,*) Y

    !    write(6,*) svecs
    !    write(6,*) stemp
    !    write(6,*) X
    !    write(6,*) Y

  end subroutine init_transform




  !
  ! Construct the Fock matrix
  !
  ! From S&O notation:
  !
  ! mu -> i
  ! nu -> j
  ! lambda -> k
  ! sigma -> l
  !
  !
  ! F(i,j) = T(i,j) + V(i,j) + G(i,j)
  !
  ! where G(i,j) = sum_kl P_kl * [ (ij|lk) - 0.5*(ik|lj) ]

  subroutine build_fock (Pao, efield, Fao)
    implicit none

    double complex, intent(in)   :: Pao(n,n)   ! density matrix
    double precision, intent(in) :: efield     ! applied electric field
    double complex, intent(out)  :: Fao(n,n)   ! Fock matrix

    double complex :: G(n,n), trans_Fao(n,n)
    double complex :: val
    integer :: i, j, k, l
    double precision :: differ
    !
    ! Compute G matrix
    !
    do i = 1, n
       do j = 1, n

          ! sum over kl
          val = 0d0
          do k = 1, n
             do l = 1, n
                val = val + Pao(k,l) * ( twoe(i,j,l,k) - 0.5d0 * twoe(i,k,l,j) )
             end do
          end do
          G(i, j) = val

       end do
    end do

    Fao = Hcore + G


    ! add applied electric field and imaginary potential
    ! XXX in z-direction for now XXX
    Fao = Fao - Dz * efield  - zi * p_gamma

    !
    ! Check to see whether the effective Hamiltonian is a Hermitian
    !
    ! trans_Fao = transpose(conjg(Fao))
    ! ! call zmat_print (n, n, trans_Fao , "Transpose of Fock matrix")

    ! differ = sum(abs(Fao - trans_Fao))

    ! if (differ > 1e-5) then
    !      write (6,*) "Error: Effective Hamiltonian (Fao) is not a Hermitian that implies non-sym"
    !    write (6,*) "differ = ", differ
    !    stop
    ! else
    !    write (6,*) "Check passed: X transform orthogonalizes S"
    ! endif

    !
    ! Print output Fock matrix
    !
    !   write (6, *) " Calculated N x N Fock matrix"
    !    write (6,*) "  i   j              Re[F]                Im[F] "
    !    write (6,*) "-------------------------------------------------"
    !     do i = 1, n
    !        do j = 1, n
    !           write (6,"(i4, i4, 5f20.4)") i, j,dble(Fao(i,j)), aimag(Fao(i,j))
    !        enddo
    !     enddo
    !     write (6, *) ""
    !     write (6, *) ""
  end subroutine build_fock


  ! see /home/klopata/nwchem/nwchem-dev-stock/src/nwdft/rt_tddft/canorg/
  !
  ! Transform Fock matrix from AO to MO basis
  !
  ! F' = X^+ F X
  !
  subroutine transform_fock_ao2mo (Fao, Fmo)
    implicit none
    double complex, intent(in) ::  Fao(n,n)
    double complex, intent(out) :: Fmo(n,n)

    Fmo = matmul(transpose(X), matmul(Fao, X))
  end subroutine transform_fock_ao2mo


  !
  ! Transform density matrix from AO to MO basis
  !
  ! P' = Y^+ P Y
  !
  subroutine transform_dens_ao2mo (Pao, Pmo)
    implicit none
    double complex, intent(in) ::  Pao(n,n)
    double complex, intent(out) :: Pmo(n,n)

    double complex :: tmp(n,n)
    integer :: i, j

    !    tmp = matmul(Pao, Y)
    Pmo = matmul(transpose(Y), matmul(Pao, Y))

    !    do i = 1, n
    !       do j = 1, n
    !          write (6, *) "tmp: ", i, j, tmp(i,j)
    !       enddo
    !   enddo
  end subroutine transform_dens_ao2mo


  !
  ! Transform density matrix from MO to AO basis
  !
  ! P = X P' X^+
  !
  subroutine transform_dens_mo2ao (Pmo, Pao)
    implicit none
    double complex, intent(in) ::  Pmo(n,n)
    double complex, intent(out) :: Pao(n,n)

    Pao = matmul(X, matmul(Pmo, transpose(X)))
  end subroutine transform_dens_mo2ao

  !
  ! Transform fock matrix from MO to AO basis
  ! F = Y^+ F' Y
  !
  !  subroutine transform_fock_mo2ao (Fmo, Fao)
  !    implicit none
  !    double complex, intent(in) ::  Fmo(n,n)
  !    double complex, intent(out) :: Fao(n,n)
  !
  !    Fao = matmul(transpose(Y), matmul(Fmo, Y))
  !  end subroutine transform_fock_mo2ao


  !
  ! Check that transforms work ok:
  !
  ! Orthogonalization of overlap:
  ! X^+ S X = 1
  !
  ! Apply forward and reverse transform to dens mat recovers original matrix:
  ! P = X P' X^+ ,
  ! where P' = Y P Y^+
  !
  subroutine check_transforms (Pmo)
    implicit none

    double complex, intent(in)  :: Pmo(n,n)
    double complex              :: tmp(n,n), Imat(n,n), Pmo_new(n,n)
    integer                     :: i
    double precision            :: diff

    ! complex identity matrix
    Imat = 0d0
    do i = 1, n
       Imat(i,i) = (1d0, 0d0)   !on-diagonal is 1
    enddo

    ! 1 = X^+ S X
    tmp = matmul(transpose(X), matmul(S, X))

    diff = sum(abs(tmp - Imat))

    if (diff > 1e-5) then
       write (6,*) "ERROR: X transform does not orthogonalize S"
       write (6,*) "diff = ", diff
       stop
    else
       write (6,*) "Check passed: X transform orthogonalizes S"
    endif


    !
    ! Convert from Pmo -> Pao, then convert back to Pmo and check that
    ! the answer is the same.
    !
    tmp = 0d0
    call transform_dens_mo2ao (Pmo, tmp)
    call transform_dens_ao2mo (tmp, Pmo_new)

    diff = sum(abs(Pmo - Pmo_new))  !diff between old and new Pmo mats

    if (diff > 1e-6) then
       write (6,*) "ERROR: Bad density matrix transforms"
       write (6,*) "diff = ", diff
       stop
    else
       write (6,*) "Check passed: Density matrix transforms"
    endif


    !
    ! Check that Y is the inverse of X
    ! XXX IS THIS TRUE?
    !
    ! tmp = 0d0
    ! tmp = matmul(Y, X)

    ! diff = sum(abs(tmp - Imat))  !diff from ident mat

    ! if (diff > 1e-6) then
    !    write (6,*) "ERROR: X is not the inverse of Y"
    !    write (6,*) "diff = ", diff
    !    stop
    ! else
    !    write (6,*) "Check passed: X is the inverse of Y"
    ! endif

  end subroutine check_transforms



  !
  ! Compute dP'/dt = -i [F', P'] (in MO basis)
  !
  subroutine calc_dPdt (Pmo, Fmo, dPdt)
    implicit none
    double complex, intent(in)  :: Pmo(n,n), Fmo(n,n)
    double complex, intent(out) :: dPdt(n,n)

    dPdt = -zi * (matmul(Fmo, Pmo) - matmul(Pmo, Fmo))
  end subroutine calc_dPdt


  !
  ! Compute electronic energy:
  !
  ! E = 0.5 * sum_{mu nu} P(mu,nu) * [ Hcore(mu,nu) + F(mu,nu) ]
  !
  subroutine calc_energy (Pao, Fao, energy)
    implicit none
    double complex, intent(in)   :: Pao(n,n), Fao(n,n)
    double precision, intent(out)  :: energy

    integer :: i, j

    energy = 0d0

    do i = 1, n
       do j = 1, n
          energy = energy + 0.5d0 * dble(Pao(i,j)) * dble(Hcore(i,j) + dble(Fao(i,j)))
       end do
    end do

    !    write (6, *) "Energy = ", energy

  end subroutine calc_energy

  subroutine calc_efield (freq, t, efield)
    implicit none
    double precision :: enve
    double precision, intent(in)   :: t, freq
    double precision, intent(out)  :: efield



    if (t<1d-6) then
       efield=0.0001d0
    else
       efield = 0d0
    endif


  end subroutine calc_efield
end module hf_module



!Occupation matrix: normalizing the denisty of electrons
!n_occ(1,1)=2;
!n_occ(1,2)=0;
!n_occ(2,1)=0;
!n_occ(2,2)=0;


! Main time propagation routine
!
! i dP'/dt = [ F'(t), P'(t) ]
!
! where prime means MO basis.  Note, F'(t) is built from P(t) (in the
! AO basis) then transformed to MO basis.
!
! XXX First let's do that assuming F'(t) = F'0 (time-independent)
!
! 1) Start with P (ground state dens mat in AO basis)
! 2) Build F[P] in AO basis (includes E-field)
! 3) Convert P and F into MO basis
! 4) Step P' forward in time using: i dP'/dt = [ F', P' ]
! 5) Convert P'(t+dt) back into AO basis.  You are now at time t+dt
! 6) go back to step (2)
!
program main
  use hf_module
  use matutils
  implicit none


  double complex :: Pao(n,n), Pmo(n,n)   ! density matrix in AO/MO basis
  double complex :: Fao(n,n), Fmo(n,n)   ! Fock matrix in AO/MO basis
  double complex :: dPdt(n,n)
  double precision :: energy, efield
  double complex :: Tr, Fmo_old(n,n),Fmo_new(n,n),Pmo_crude(n,n), DP(n,n), Pmid(n,n)
  double precision   :: t                ! current time
  double precision :: a, b, c, freq, eStrength
  integer :: it, i, j
  double precision :: dipmom, charge, diffe, differe
  double complex :: zvals(n) , zvecs(n,n), zvecs_ao(n,n), inv(n,n), omega(n,n), evals(n,n), U(n,n), inv_U(n,n),expA(n,n)
  double complex :: k1, k2
  double complex :: ident(n,n), ident_const(n,n)
  logical :: do_interpolate


  !  xxxxxx define: it, nt, dt, etc etc etc
  do_interpolate = .false.

  call init_integrals ()
  call init_densmat (Pao)
  call init_transform_new ()
  call check_transforms (Pao)


  ! compute ground state Fock and energies
  efield = 0d0
  call build_fock (Pao, efield, Fao)
  call calc_energy (Pao, Fao, energy)
  write (6, *) "Ground State Electronic Energy = ", energy


  freq = 0.1d0
  !  eStrength = 0.001d0
  !  eStrength = 0d0

  !
  ! Main time loop
  !
  ! Fmo_old=Fmo ! Starting matrix (current matrix, F(t))
  do it = 1, nt
     t  = (it-1) * dt
     ! efield = eStrength * Sin(freq * t)

     !     call calc_efield (freq, t, efield)
     efield = 0d0  ! no e-field if we are doing imag time prop
     call build_fock (Pao, efield, Fao)
     call transform_dens_ao2mo (Pao, Pmo)
     call transform_fock_ao2mo (Fao, Fmo)


     call calc_dPdt(Pmo, Fmo, dPdt)
     k1 = dt * dpdt(n,n)


     Pmid = Pmo + (k1/2d0)

     call calc_dPdt(Pmid, Fmo, dPdt)
     k2 = dt * dPdt(n,n)

     Pmo = Pmo + k2

     ! Diagonalization is performed of the Fock matrix is performed here

     ! double complex :: zvals(n), zvecs(n,n), Fmo(n,n)
     !integer :: i, j

     !    call zmat_print (n, n, Fmo, "Input Matrix")
     call zmat_diag (n, Fmo, zvals, zvecs)
     !   call zmat_print (n, n, zvecs, "Eigenvector Matrix")



     !    call zmat_print (n, n, Fmo, "Coeffiecent matrix")

     ! Y= matmul(zvecs,s12)
     ! X= matmul (Y, transpose(zvecs))

     ! Sanity check
     ! q= matmul(transpose(X), S)
     ! a= matmul(q,X)

     ! write(6,*) "Transform Matrix, X", X
     ! write(6,*) "Identity Matrix, a",a

     ! End of diagonalization



     !Check to see the coefficient matrix in AO basis (eigenvectors of Fock matrix) match NWchem output.
     zvecs_ao = matmul(X, zvecs);
     !    call zmat_print (n, n, zvecs_ao, "AO Eigenvectors")
     !    call zmat_print (n, n, zvecs, "MO Eigenvectors")
     !   end subroutine transform_zvecs_mo2ao
     !

     !
     ! Normalization of the density matrix is performed here
     !


     call zmat_inv (n, transpose(conjg(zvecs)), inv_U)
     Pmo = matmul(matmul(inv_U,  n_occ), inv)

     !
     !call zmat_inv (n, transpose(conjg(zvecs)), inv_U)
     !Pmo = matmul(matmul(zvecs,  n_occ), inv)

     !    subroutine transform_zvecs_mo2ao (zvecs, zvecs_ao)
     !      implicit none
     !      double complex, intent(in) ::  zvecs(n,n)
     !      double complex, intent(out) :: zvecs_ao(n,n)


     ! Check to see coefficient matrix in MO basis is unitary

     call zmat_inv (n, zvecs, inv)

     ! if (.not. zmat_inv_check (n, zvecs, inv)) then
     !    write (6,*) "*** ERROR: Inversion check failed ***";
     !    stop
     ! else
     !    write (6,*) "Inversion check passed";
     ! endif

     ident = matmul(zvecs, inv)
     !ident = matmul(zvecs, transpose(conjg(zvecs)))
     !    call zmat_print (n, n, ident, "Unitary check")



     !
     !Check to see the transpose of coefficient matrix is the same as its inverse
     !
     differe = abs(sum(inv_U -  transpose(conjg(zvecs))))
     if (diffe > 1e-6) then
        write (6,*) "ERROR: coefficient matrix inverse NOT the same as its complex cojugate transpose"
        write (6,*) "differe = ", differe
        stop
     else
        write (6,*) "Check passed: coefficient matrix inverse is the same as its complex cojugate transpose"
     endif



     do i = 1, n
        do j = 1, n

           if (i == j) then  !on-diagonal
              ident_const(i,i)= 1d0
           else  !off-diagonal
              ident_const(i,j) = 0d0
           endif

        enddo
     enddo

     diffe = abs(sum(ident - ident_const))
     if (diffe > 1e-6) then
        write (6,*) "ERROR: coefficient matrix in MO basis is not unitary"
        write (6,*) "diffe = ", diffe
        stop
     else
        write (6,*) "Check passed: coefficient matrix in MO basis is unitary"
     endif

     !
     ! Damping is performed to get convergence of the energy
     !

     Pmo = 0.01d0 * Pmo + 0.99d0 * Pmid
     ! Compute time-dependent properties (dipole moment and energy)
     !
     call transform_dens_mo2ao (Pmo, Pao)

     ! compute dipole moment
     DP = matmul(Dz,Pao)
     Tr = 0d0
     do i = 1, n
        Tr = Tr +  DP(i,i)
     enddo
     dipmom = dble (Tr)

     ! compute energy
     call calc_energy (Pao, Fao, energy)

     ! compute total electronic charge

     Tr = 0d0
     do i = 1, n
        Tr = Tr +  Pmo(i,i)
     enddo
     charge = dble (Tr)


     write(700,*) t, efield, dipmom, energy
     write (750, *) t,  charge
     write (800, *) t, dble(zvecs_ao)
     write(1010, *) t, dble(zvals)
     write(1020, *) t, aimag(zvals)
     write (900, *) ident
     write(6,*) "TDHF: ", t, efield, dipmom, energy
     call flush(700)




  enddo



  !  efield = 0d0
  !  freq = 0.23d0
  !  eStrength = 0.1
  !  a = dble(Pao(1,1))
  !  b = dble (Pao(2,2))
  !  c=a+b
  !  efield = eStrength * Sin(freq* t)

  !  Fao = Fao - (Dz * efield)


  ! call build_fock (Pao, efield, Fao)
  ! call calc_energy (Pao, Fao, energy)




  ! call transform_fock_ao2mo (Fao, Fmo)
  ! call transform_dens_ao2mo (Pao, Pmo)

  ! call calc_dPdt (Pmo, Fmo, dPdt)
  ! write (6, *) dPdt





end program main




!   t = -dt
!   do  it = 1, nt
!       t = t + dt
! !!! k1 = dt * dP/dt(t) = -i * [Fmo(t), P(t)]

!      call calc_dPdt(n, Fmo, Pmo, dPdt)
!      k1 = dt * dPdt

! !!! k2 = dt * dP/dt(t+dt/2)
!      P1mo = Pmo + (k1/2d0)
!      call calc_dPdt(n, Fmo, P1mo, dPdt)
!      k2 = dt * dPdt


! !!! K3 = dt * dP1/dt(t+dt/2)
!      P2mo = Pmo + (k2/2d0)
!      call calc_dPdt(n, Fmo, P2mo, dPdt)
!      k3 = dt * dPdt


! !!! K4 = dt * dP2/dt(t+dt)
!      P3mo = Pmo + k3
!      call calc_dPdt(n, Fmo, P3mo, dPdt)
!      k4 = dt * dPdt


!      Pmo = Pmo + (k1/6d0) + (k2/3d0) + (k3/3d0) + (k4/6d0)

!   enddo



!
!    double complex :: zvals(n), zvecs(n,n), s12(n,n)
!    integer :: i, j
!
!    call print_zmat (n, n, S, "Input Matrix")
!    call mat_diag (n, S, zvals, zvecs)
!    call print_zm                                                                                                                         at (n, 1, zvals, "Eigenvalues")
!    call print_zmat (n, n, zvecs, "Eigenvector Matrix")
!
!   do i = 1, n
!       do j = 1, n
!
!          if (i == j) then  !on-diagonal = 1 / s^1/2
!             s12(i,i)= 1d0/(sqrt(zvals(i)))
!          else  !off-diagonal
!             s12(i,j) = 0d0
!          endif
!
!       enddo
!    enddo
!
!    call print_zmat (n, n, s12, "s12 Matrix")
!
!    Y= matmul(zvecs,s12)
!    X= matmul (Y, transpose(zvecs))
!
! Sanity check
!    q= matmul(transpose(X), S)
!    a= matmul(q,X)

!    write(6,*) "Transform Matrix, X", X
!    write(6,*) "Identity Matrix, a",a


!   stop



! ! XXX TMP: Test diagonalization
! call zmat_diag (n, omega, zvals, zvecs)
!  if (.not. zmat_diag_check (n, omega, zvals, zvecs)) then
!    write (6,*) "*** Diag check failed ***";
!    stop
! else
!    write (6,*) "Diag check passed";
! endif

! call zmat_inv (n, zvecs, inv)

! if (.not. zmat_inv_check (n, zvecs, inv)) then
!    write (6,*) "*** ERROR: Inversion check failed ***";
!    stop
! else
!    write (6,*) "Inversion check passed";
     ! endif
