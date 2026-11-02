!***********************************************************************
!
!  HOPO_PEF.f90 — Standalone evaluation of the HOPO potential energy
!                 function (PEF) from ab initio data.
!
!  Molecule:  HOPO  (metaphosphorous acid, H–O–P=O)
!  Method:    6-dimensional analytic PEF in Morse–cosine coordinates
!
!
!  Usage:
!    Compile:  gfortran -O2 -o HOPO_PEF HOPO_PEF.f90
!    Run:      ./HOPO_PEF < HOPO_PEF.inp
!
!  Input file format (HOPO_PEF.inp):
!    Lines beginning with '!' are comments and are ignored.
!    First data block:
!      Line 1:  npar  (number of non-zero linear parameters)
!    Equilibrium geometry and Morse parameters (8 values):
!      re_OP    (Angstrom)
!      re_PO    (Angstrom)
!      re_OH    (Angstrom)
!      ae_OPO   (degrees)
!      ae_POH   (degrees)
!      a1       (Morse parameter for O=P stretch, Angstrom^-1)
!      a2       (Morse parameter for P-O stretch, Angstrom^-1)
!      a3       (Morse parameter for O-H stretch, Angstrom^-1)
!    Linear expansion coefficients (npar lines):
!      n1 n2 n3 n4 n5 n6   coeff   (integers + real in cm^-1)
!    Geometry block (one geometry per line):
!      r_OP(A)  r_PO(A)  r_OH(A)  ang_OPO(deg)  ang_POH(deg)  tau(deg)
!      terminated by a line with r_OP < 0.
!
!  Output:
!    For each geometry: coordinates and V (cm^-1), relative to equilibrium.
!
!  Internal coordinate definitions:
!    xi1 = 1 - exp(-a1*(r1 - r1e))              [Morse, O=P stretch]
!    xi2 = 1 - exp(-a2*(r2 - r2e))              [Morse, P-O stretch]
!    xi3 = 1 - exp(-a3*(r3 - r3e))              [Morse, O-H stretch]
!    xi4 = cos(beta1) - cos(beta1e)             [OPO bend, beta1 = ang_OPO]
!    xi5 = cos(rho2e) - cos(rho2)               [POH umbrella,
!                                                 rho = pi - ang_POH]
!    xi6 = tau                                  [torsion, stored but
!                                                 enters via cos(n6*tau)]
!
!  The potential is evaluated as:
!    V = sum_k  c_k * xi1^n1 * xi2^n2 * xi3^n3 * xi4^n4 * xi5^n5
!                   * F_tors(n6, beta1, rho2, tau)
!  where
!    F_tors = 1                                  for n6 = 0
!    F_tors = sin(beta1) * sin(rho2) * cos(n6*tau)  for n6 > 0
!
!***********************************************************************

program HOPO_PEF

  implicit none

  ! ---- parameters -------------------------------------------------------
  integer,  parameter :: rk      = selected_real_kind(15, 307)
  integer,  parameter :: nparmax = 4000   ! max number of linear parameters
  integer,  parameter :: ndeg    = 6      ! dimensionality of the PEF

  ! ---- local variables --------------------------------------------------
  real(rk), parameter :: pi = 4.0_rk * atan2(1.0_rk, 1.0_rk)
  real(rk), parameter :: deg2rad = pi / 180.0_rk

  ! non-linear parameters
  real(rk) :: r1e, r2e, r3e        ! equilibrium bond lengths (Angstrom)
  real(rk) :: beta1e, beta2e        ! equilibrium bond angles (radians)
  real(rk) :: a1, a2, a3           ! Morse parameters (Angstrom^-1)

  ! linear expansion
  integer  :: npar
  integer  :: ipow(ndeg, nparmax)  ! power indices
  real(rk) :: coeff(nparmax)       ! expansion coefficients (cm^-1)

  ! geometry input
  real(rk) :: r1, r2, r3           ! bond lengths (Angstrom)
  real(rk) :: beta1_deg, beta2_deg ! angles (degrees)
  real(rk) :: tau_deg              ! torsion (degrees)

  ! working variables
  real(rk) :: beta1, beta2, tau    ! angles in radians
  real(rk) :: rho2, rho2e         ! umbrella coordinates
  real(rk) :: xi(ndeg)             ! internal coordinates
  real(rk) :: term, V, Ve

  integer  :: i, k, igeom
  character(len=200) :: linebuf

  ! ---- read input -------------------------------------------------------
  ! skip comment lines
  do
    read(*,'(a)') linebuf
    if (linebuf(1:1) /= '!') exit
  end do

  ! number of linear parameters
  read(linebuf, *) npar
  if (npar > nparmax) then
    write(*,'(a,i0,a,i0)') 'ERROR: npar = ', npar, ' exceeds nparmax = ', nparmax
    stop
  end if

  ! equilibrium geometry and Morse parameters (skipping comment lines)
  r1e     = read_real()
  r2e     = read_real()
  r3e     = read_real()
  beta1e  = read_real() * deg2rad   ! convert degrees -> radians
  beta2e  = read_real() * deg2rad
  a1      = read_real()
  a2      = read_real()
  a3      = read_real()

  rho2e = pi - beta2e

  ! linear expansion coefficients
  do k = 1, npar
    call read_coeff_line(ipow(:, k), coeff(k))
  end do

  ! ---- compute equilibrium energy Ve (needed to zero the PEF) ----------
  call eval_PEF(r1e, r2e, r3e, beta1e, beta2e*180.0_rk/pi, 0.0_rk, Ve)

  ! ---- print header -----------------------------------------------------
  write(*,'(a)') repeat('*', 72)
  write(*,'(a)') '  HOPO Potential Energy Function — standalone evaluation'
  write(*,'(a)') repeat('*', 72)
  write(*,'(a)')
  write(*,'(a)') '  Equilibrium geometry:'
  write(*,'(a,f12.8,a)') '    r(O=P)    = ', r1e,         ' Angstrom'
  write(*,'(a,f12.8,a)') '    r(P-O)    = ', r2e,         ' Angstrom'
  write(*,'(a,f12.8,a)') '    r(O-H)    = ', r3e,         ' Angstrom'
  write(*,'(a,f12.8,a)') '    ang(OPO)  = ', beta1e/deg2rad, ' degrees'
  write(*,'(a,f12.8,a)') '    ang(POH)  = ', beta2e/deg2rad, ' degrees'
  write(*,'(a,f12.8,a)') '    a1        = ', a1,           ' Angstrom^-1'
  write(*,'(a,f12.8,a)') '    a2        = ', a2,           ' Angstrom^-1'
  write(*,'(a,f12.8,a)') '    a3        = ', a3,           ' Angstrom^-1'
  write(*,'(a)')
  write(*,'(a,i0,a)') '  Number of linear parameters: ', npar
  write(*,'(a,f12.4,a)') '  V(eq) = ', Ve, ' cm^-1 (absolute; output is relative)'
  write(*,'(a)')
  write(*,'(a)') repeat('-', 72)
  write(*,'(a)') '   r(O=P)/A   r(P-O)/A   r(O-H)/A  a(OPO)/deg  a(POH)/deg  tau/deg       V/cm-1'
  write(*,'(a)') repeat('-', 72)

  ! ---- geometry loop ----------------------------------------------------
  igeom = 0
  do
    ! read next geometry (skip comment lines)
    do
      read(*,'(a)', end=999) linebuf
      if (linebuf(1:1) /= '!') exit
    end do

    read(linebuf, *) r1, r2, r3, beta1_deg, beta2_deg, tau_deg

    if (r1 < 0.0_rk) exit   ! terminator

    igeom = igeom + 1
    call eval_PEF(r1, r2, r3, beta1_deg * deg2rad, beta2_deg, tau_deg, V)

    write(*,'(3f11.6, 2f12.5, f10.3, f14.4)') &
          r1, r2, r3, beta1_deg, beta2_deg, tau_deg, V - Ve
  end do

999 continue
  write(*,'(a)') repeat('-', 72)
  write(*,'(a,i0,a)') '  ', igeom, ' geometries evaluated.'
  write(*,'(a)') repeat('*', 72)

contains

  !-------------------------------------------------------------------
  ! eval_PEF: evaluate potential at one geometry
  !
  ! Input:
  !   r1_      bond length r(O=P)  in Angstrom
  !   r2_      bond length r(P-O)  in Angstrom
  !   r3_      bond length r(O-H)  in Angstrom
  !   beta1_   angle OPO           in RADIANS
  !   beta2_deg angle POH          in DEGREES  (converted internally)
  !   tau_deg  dihedral tau        in DEGREES  (converted internally)
  !
  ! Output:
  !   Vout     potential energy in cm^-1 (absolute, before Ve subtraction)
  !-------------------------------------------------------------------
  subroutine eval_PEF(r1_, r2_, r3_, beta1_, beta2_deg_, tau_deg_, Vout)

    real(rk), intent(in)  :: r1_, r2_, r3_, beta1_, beta2_deg_, tau_deg_
    real(rk), intent(out) :: Vout

    real(rk) :: xi_(ndeg)
    real(rk) :: beta2_, tau_, rho2_, term_
    real(rk) :: ftors, sb1, sr2
    integer  :: k_, n6_

    beta2_ = beta2_deg_ * deg2rad
    tau_   = tau_deg_   * deg2rad
    rho2_  = pi - beta2_

    ! Morse-type stretching coordinates
    xi_(1) = 1.0_rk - exp(-a1 * (r1_ - r1e))
    xi_(2) = 1.0_rk - exp(-a2 * (r2_ - r2e))
    xi_(3) = 1.0_rk - exp(-a3 * (r3_ - r3e))

    ! Bending coordinate for OPO (cos-type)
    xi_(4) = cos(beta1_) - cos(beta1e)

    ! Umbrella bending coordinate for POH (via rho = pi - angle)
    xi_(5) = cos(rho2e) - cos(rho2_)

    ! Torsion angle (radians; enters only through cos(n6*tau) * sin*sin)
    xi_(6) = tau_

    ! Precompute sin factors needed for the torsional coupling
    sb1 = sin(beta1_)
    sr2 = sin(rho2_)

    Vout = 0.0_rk

    do k_ = 1, npar

      n6_ = ipow(6, k_)

      ! torsional factor
      if (n6_ == 0) then
        ftors = 1.0_rk
      else
        ftors = sb1 * sr2 * cos(real(n6_, rk) * tau_)
      end if

      ! polynomial term in stretching/bending coordinates
      term_ = coeff(k_)                       &
              * xi_(1)**ipow(1, k_)           &
              * xi_(2)**ipow(2, k_)           &
              * xi_(3)**ipow(3, k_)           &
              * xi_(4)**ipow(4, k_)           &
              * xi_(5)**ipow(5, k_)           &
              * ftors

      Vout = Vout + term_

    end do

  end subroutine eval_PEF

  !-------------------------------------------------------------------
  ! read_real: read one real value from stdin, skipping comment lines
  !-------------------------------------------------------------------
  function read_real() result(val)
    real(rk) :: val
    character(len=200) :: buf
    do
      read(*,'(a)') buf
      if (buf(1:1) /= '!') then
        read(buf, *) val
        return
      end if
    end do
  end function read_real

  !-------------------------------------------------------------------
  ! read_coeff_line: read one "n1 n2 n3 n4 n5 n6  coeff" line,
  !                  skipping comment lines
  !-------------------------------------------------------------------
  subroutine read_coeff_line(ipow_, c_)
    integer,  intent(out) :: ipow_(ndeg)
    real(rk), intent(out) :: c_
    character(len=200) :: buf
    do
      read(*,'(a)') buf
      if (buf(1:1) /= '!') then
        read(buf, *) ipow_(1), ipow_(2), ipow_(3), &
                     ipow_(4), ipow_(5), ipow_(6), c_
        return
      end if
    end do
  end subroutine read_coeff_line

end program HOPO_PEF
