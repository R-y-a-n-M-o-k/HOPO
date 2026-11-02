!***********************************************************************
!
!  HOPO_DMS.f90 — Standalone evaluation of the HOPO dipole moment
!                 surface (DMS) from Cartesian coordinates.
!
!  Molecule:  HOPO  (metaphosphorous acid, H–O–P=O)
!  Method:    6-dimensional analytic DMS: polynomial in the bond-length
!             and bond-angle displacements times a Fourier series in
!             the torsion angle.  Units: Angstrom, degrees, Debye.
!
!  Usage:
!    Compile:  gfortran -O2 -o HOPO_DMS HOPO_DMS.f90
!    Run:      ./HOPO_DMS < HOPO_DMS.inp
!
!  Atom numbering (TROVE), with the molecule read as  H - O - P = O :
!
!        H4 - O2 - P1 = O3
!
!    atom 1 = P
!    atom 2 = O   (hydroxyl oxygen, bonded to H)
!    atom 3 = O   (terminal oxygen)
!    atom 4 = H
!
!  Input file format (HOPO_DMS.inp):
!    Lines beginning with '!' and blank lines are comments and are ignored.
!    First data line:
!      npar_x  npar_y  npar_z   (number of non-zero linear parameters
!                                for mu_x, mu_y and mu_z)
!    Expansion-centre geometry (5 values):
!      re_PO1   r(P1-O2)       (Angstrom)
!      re_PO2   r(P1=O3)       (Angstrom)
!      re_OH    r(O2-H4)       (Angstrom)
!      ae_OPO   ang(O2-P1-O3)  (degrees)
!      ae_POH   ang(P1-O2-H4)  (degrees)
!    Linearity threshold (1 value):
!      lin_tol  if sin(ang(OPO)) or sin(ang(POH)) is below lin_tol the
!               angle is treated as linear and the torsion is set to 0.
!    Linear expansion coefficients, three consecutive blocks
!    (mu_x: npar_x lines, mu_y: npar_y lines, mu_z: npar_z lines):
!      n1 n2 n3 n4 n5 n6   coeff   (integers + real in Debye)
!    Geometry block (one geometry per line, 12 numbers in Angstrom):
!      x1 y1 z1   x2 y2 z2   x3 y3 z3   x4 y4 z4
!      atoms in TROVE order (P1, O2, O3, H4), in ANY Cartesian frame;
!      terminated by a line beginning with END (or by end of file).
!
!  Output:
!    For each geometry: the z-matrix coordinates, the dipole in the DMS
!    (fit) frame and the dipole transformed back to the Cartesian frame
!    of the input geometry (the final result), all in Debye.
!
!  Internal (z-matrix) coordinates computed from the Cartesians:
!    r1     = |x2 - x1|                    r(P-O)   P-OH single bond
!    r2     = |x3 - x1|                    r(P=O)   terminal bond
!    r3     = |x4 - x2|                    r(O-H)
!    alpha1 = ang(O2-P1-O3)
!    alpha2 = ang(P1-O2-H4)
!    tau    = dihedral H4-O2-P1-O3, in (-pi, pi]
!             (tau = 0: H and terminal O cis; tau > 0 when H lies on
!             the +y side of the DMS frame).  If alpha1 or alpha2 is
!             linear the dihedral is undefined and tau is set to 0.
!
!  Expansion variables:
!    xi1 = r1 - r1e
!    xi2 = r2 - r2e
!    xi3 = r3 - r3e
!    xi4 = alpha1 - alpha1e                        (radians)
!    xi5 = rho2 - rho2e = alpha2e - alpha2         (radians,
!                                                   rho2 = pi - alpha2)
!
!  The DMS (components in the DMS frame) is evaluated as:
!    mu_x = sum_k cx_k * xi1^n1 ... xi5^n5 *            cos(n6*tau)
!    mu_y = sum_k cy_k * xi1^n1 ... xi5^n5 * sin(tau) * cos(n6*tau)
!    mu_z = sum_k cz_k * xi1^n1 ... xi5^n5 *            cos(n6*tau)
!
!  DMS (fit) frame, built from the input Cartesians:
!    e_z = b1 / |b1|                 b1 = x2 - x1   (along P -> O(H))
!    e_y = (b1 x b2) / |b1 x b2|     b2 = x3 - x1   (normal to OPO plane)
!    e_x = e_y x e_z                 in the OPO plane, O3 at x > 0
!    If O-P-O is linear, the x-z plane is taken as the P-O-H plane with
!    H at x > 0 (consistent with tau = 0); if the whole molecule is
!    linear, e_y is an arbitrary unit vector perpendicular to e_z.
!
!  Back-transformation to the input frame (final result):
!    mu = mu_x * e_x + mu_y * e_y + mu_z * e_z
!
!***********************************************************************

program HOPO_DMS

  implicit none

  ! ---- parameters -------------------------------------------------------
  integer,  parameter :: rk      = selected_real_kind(15, 307)
  integer,  parameter :: nparmax = 2000   ! max number of linear parameters per component
  integer,  parameter :: ndeg    = 6      ! dimensionality of the DMS
  integer,  parameter :: natoms  = 4      ! number of atoms
  integer,  parameter :: maxpow  = 12     ! max power / torsional multiplicity allowed

  ! ---- local variables --------------------------------------------------
  real(rk), parameter :: pi = 4.0_rk * atan2(1.0_rk, 1.0_rk)
  real(rk), parameter :: deg2rad = pi / 180.0_rk
  real(rk), parameter :: rmin = 1.0e-8_rk  ! minimum allowed bond length (Angstrom)

  ! non-linear parameters
  real(rk) :: r1e, r2e, r3e        ! expansion-centre bond lengths (Angstrom)
  real(rk) :: alpha1e, alpha2e     ! expansion-centre bond angles (radians)
  real(rk) :: lin_tol              ! linearity threshold on sin(angle)

  ! linear expansion (third index: 1 = mu_x, 2 = mu_y, 3 = mu_z)
  integer  :: npar(3)
  integer  :: ipow(ndeg, nparmax, 3)  ! power indices
  real(rk) :: coeff(nparmax, 3)       ! expansion coefficients (Debye)

  ! geometry input / results
  real(rk) :: xyz(3, natoms)       ! Cartesian coordinates (Angstrom), TROVE order
  real(rk) :: zmat(ndeg)           ! r1, r2, r3 (A), alpha1, alpha2, tau (rad)
  real(rk) :: mu_mol(3)            ! dipole in DMS frame (Debye)
  real(rk) :: mu(3)                ! dipole in input frame (Debye)

  integer  :: ic, k, igeom, ios
  logical  :: eof
  character(len=512) :: linebuf
  character(len=1), parameter :: cname(3) = ['x', 'y', 'z']

  ! ---- read input -------------------------------------------------------
  ! number of linear parameters for each component
  call next_data_line(linebuf, eof)
  if (eof) call input_error('number of parameters')
  read(linebuf, *) npar(1), npar(2), npar(3)
  do ic = 1, 3
    if (npar(ic) > nparmax) then
      write(*,'(a,a,a,i0,a,i0)') 'ERROR: npar_', cname(ic), ' = ', npar(ic), &
                                ' exceeds nparmax = ', nparmax
      stop
    end if
  end do

  ! expansion-centre geometry and linearity threshold (skipping comment lines)
  r1e     = read_real()
  r2e     = read_real()
  r3e     = read_real()
  alpha1e = read_real() * deg2rad   ! convert degrees -> radians
  alpha2e = read_real() * deg2rad
  lin_tol = read_real()

  ! linear expansion coefficients: mu_x block, mu_y block, mu_z block
  do ic = 1, 3
    do k = 1, npar(ic)
      call read_coeff_line(ipow(:, k, ic), coeff(k, ic))
      if (any(ipow(:, k, ic) < 0) .or. any(ipow(:, k, ic) > maxpow)) then
        write(*,'(a,a,a,i0,a,i0)') 'ERROR: power out of range in mu_', cname(ic), &
                                  ' parameter ', k, '; maxpow = ', maxpow
        stop
      end if
    end do
  end do

  ! ---- print header -----------------------------------------------------
  write(*,'(a)') repeat('*', 72)
  write(*,'(a)') '  HOPO Dipole Moment Surface — standalone evaluation from Cartesians'
  write(*,'(a)') repeat('*', 72)
  write(*,'(a)')
  write(*,'(a)') '  Atom order (TROVE):  1 = P,  2 = O(H),  3 = O,  4 = H'
  write(*,'(a)')
  write(*,'(a)') '  Expansion-centre geometry:'
  write(*,'(a,f12.8,a)') '    r(P-O)    = ', r1e,             ' Angstrom'
  write(*,'(a,f12.8,a)') '    r(P=O)    = ', r2e,             ' Angstrom'
  write(*,'(a,f12.8,a)') '    r(O-H)    = ', r3e,             ' Angstrom'
  write(*,'(a,f12.8,a)') '    ang(OPO)  = ', alpha1e/deg2rad, ' degrees'
  write(*,'(a,f12.8,a)') '    ang(POH)  = ', alpha2e/deg2rad, ' degrees'
  write(*,'(a,es12.4)')  '    lin_tol   = ', lin_tol
  write(*,'(a)')
  write(*,'(a,i0,a,i0,a,i0)') '  Number of linear parameters: mu_x ', npar(1), &
                              ',  mu_y ', npar(2), ',  mu_z ', npar(3)
  write(*,'(a)')
  write(*,'(a)') '  DMS frame: z along P->O(H), y normal to OPO plane, x = y cross z'
  write(*,'(a)') '  mu_X, mu_Y, mu_Z: dipole in the Cartesian frame of the input (Debye)'
  write(*,'(a)')
  write(*,'(a)') repeat('-', 139)
  write(*,'(a)') '  geom  r(P-O)/A  r(P=O)/A  r(O-H)/A  a(OPO)/deg a(POH)/deg   tau/deg' // &
                 '   mu_x(DMS)   mu_y(DMS)   mu_z(DMS)        mu_X        mu_Y        mu_Z'
  write(*,'(a)') repeat('-', 139)

  ! ---- geometry loop ----------------------------------------------------
  igeom = 0
  do
    ! read next geometry (skip comment lines), stop at END or end of file
    call next_data_line(linebuf, eof)
    if (eof) exit
    if (is_end(linebuf)) exit

    read(linebuf, *, iostat=ios) xyz
    if (ios /= 0) then
      write(*,'(a)') 'ERROR: could not read 12 Cartesian coordinates from line:'
      write(*,'(a)') trim(linebuf)
      stop
    end if

    igeom = igeom + 1
    call dipole_HOPO(xyz, mu, zmat, mu_mol)

    write(*,'(i6, 3f10.6, 2f11.5, f10.4, 6f12.8)') &
          igeom, zmat(1:3), zmat(4:6)/deg2rad, mu_mol, mu
  end do

  write(*,'(a)') repeat('-', 139)
  write(*,'(a,i0,a)') '  ', igeom, ' geometries evaluated.'
  write(*,'(a)') repeat('*', 72)

contains

  !-------------------------------------------------------------------
  ! dipole_HOPO: dipole moment of HOPO at one Cartesian geometry
  !
  ! Input:
  !   xyz_(3,4)   Cartesian coordinates in Angstrom, TROVE order:
  !               xyz_(:,1) = P, xyz_(:,2) = O(H), xyz_(:,3) = O,
  !               xyz_(:,4) = H.  Any orientation / origin.
  !
  ! Output:
  !   mu_(3)      dipole moment (Debye) in the SAME Cartesian frame as
  !               the input coordinates  <-- final result
  !   zmat_(6)    (optional) r1, r2, r3 (Angstrom),
  !               alpha1, alpha2, tau (radians)
  !   mu_mol_(3)  (optional) dipole components in the DMS (fit) frame
  !-------------------------------------------------------------------
  subroutine dipole_HOPO(xyz_, mu_, zmat_, mu_mol_)

    real(rk), intent(in)            :: xyz_(3, natoms)
    real(rk), intent(out)           :: mu_(3)
    real(rk), intent(out), optional :: zmat_(ndeg)
    real(rk), intent(out), optional :: mu_mol_(3)

    real(rk) :: q_(ndeg)     ! z-matrix coordinates
    real(rk) :: e_(3, 3)     ! DMS frame unit vectors, e_(:,1)=e_x, (:,2)=e_y, (:,3)=e_z
    real(rk) :: mm_(3)       ! dipole in DMS frame

    ! Cartesians -> z-matrix coordinates and DMS frame
    call cart_to_zmat(xyz_, q_, e_)

    ! dipole components in the DMS frame
    call eval_DMS(q_, mm_)

    ! rotate back to the input Cartesian frame
    mu_(:) = mm_(1) * e_(:, 1) + mm_(2) * e_(:, 2) + mm_(3) * e_(:, 3)

    if (present(zmat_))   zmat_   = q_
    if (present(mu_mol_)) mu_mol_ = mm_

  end subroutine dipole_HOPO

  !-------------------------------------------------------------------
  ! cart_to_zmat: z-matrix coordinates and DMS frame from Cartesians
  !
  ! Input:
  !   xyz_(3,4)  Cartesian coordinates (Angstrom), TROVE order
  !
  ! Output:
  !   q_(6)      r1, r2, r3 (Angstrom), alpha1, alpha2 (radians, 0..pi),
  !              tau (radians, -pi..pi; 0 if OPO or POH is linear)
  !   e_(3,3)    unit vectors of the DMS frame as columns (e_x, e_y, e_z)
  !-------------------------------------------------------------------
  subroutine cart_to_zmat(xyz_, q_, e_)

    real(rk), intent(in)  :: xyz_(3, natoms)
    real(rk), intent(out) :: q_(ndeg)
    real(rk), intent(out) :: e_(3, 3)

    real(rk) :: b1(3), b2(3), b3(3)   ! P1->O2, P1->O3, O2->H4
    real(rk) :: n1(3), n2(3)          ! normals to the OPO and POH planes
    real(rk) :: r1_, r2_, r3_
    real(rk) :: s1, c1, s2, c2        ! sin/cos of alpha1, alpha2
    logical  :: lin1, lin2

    b1 = xyz_(:, 2) - xyz_(:, 1)
    b2 = xyz_(:, 3) - xyz_(:, 1)
    b3 = xyz_(:, 4) - xyz_(:, 2)

    ! bond lengths
    r1_ = norm2(b1)
    r2_ = norm2(b2)
    r3_ = norm2(b3)
    if (min(r1_, r2_, r3_) < rmin) then
      write(*,'(a)') 'ERROR: cart_to_zmat: coincident atoms in input geometry'
      stop
    end if

    ! bond angles via atan2 (accurate also close to 0 and 180 degrees)
    n1 = cross(b1, b2)                    ! |n1| = r1*r2*sin(alpha1)
    s1 = norm2(n1) / (r1_ * r2_)
    c1 = dot_product(b1, b2) / (r1_ * r2_)

    n2 = cross(b1, b3)                    ! |n2| = r1*r3*sin(alpha2)
    s2 = norm2(n2) / (r1_ * r3_)
    c2 = -dot_product(b1, b3) / (r1_ * r3_)

    lin1 = (s1 < lin_tol)                 ! O-P-O linear
    lin2 = (s2 < lin_tol)                 ! P-O-H linear

    q_(1) = r1_
    q_(2) = r2_
    q_(3) = r3_
    q_(4) = atan2(s1, c1)
    q_(5) = atan2(s2, c2)

    ! dihedral H4-O2-P1-O3 in (-pi, pi]; undefined (-> 0) if any angle is linear
    if (lin1 .or. lin2) then
      q_(6) = 0.0_rk
    else
      q_(6) = atan2(dot_product(cross(n1, n2), b1) / r1_, dot_product(n1, n2))
    end if

    ! DMS frame
    e_(:, 3) = b1 / r1_                               ! z along P -> O(H)
    if (.not. lin1) then
      e_(:, 2) = n1 / norm2(n1)                       ! y normal to OPO plane
    else if (.not. lin2) then
      e_(:, 2) = n2 / norm2(n2)                       ! OPO linear: use POH plane (H at x > 0)
    else
      e_(:, 2) = perpendicular(e_(:, 3))              ! fully linear: arbitrary
    end if
    e_(:, 1) = cross(e_(:, 2), e_(:, 3))              ! x = y cross z

  end subroutine cart_to_zmat

  !-------------------------------------------------------------------
  ! eval_DMS: evaluate the dipole components in the DMS frame
  !
  ! Input:
  !   q_(6)      r1, r2, r3 (Angstrom), alpha1, alpha2, tau (radians)
  !
  ! Output:
  !   mu_mol_(3) mu_x, mu_y, mu_z in the DMS frame (Debye)
  !-------------------------------------------------------------------
  subroutine eval_DMS(q_, mu_mol_)

    real(rk), intent(in)  :: q_(ndeg)
    real(rk), intent(out) :: mu_mol_(3)

    real(rk) :: xi_(5)
    real(rk) :: xpow_(0:maxpow, 5)     ! xi_i ** n
    real(rk) :: cosn_(0:maxpow)        ! cos(n*tau)
    real(rk) :: term_
    integer  :: i_, n_, k_, ic_

    ! displacement coordinates
    xi_(1) = q_(1) - r1e
    xi_(2) = q_(2) - r2e
    xi_(3) = q_(3) - r3e
    xi_(4) = q_(4) - alpha1e
    xi_(5) = alpha2e - q_(5)           ! rho2 - rho2e, with rho2 = pi - alpha2

    ! tables of powers and torsional cosines
    do i_ = 1, 5
      xpow_(0, i_) = 1.0_rk
      do n_ = 1, maxpow
        xpow_(n_, i_) = xpow_(n_ - 1, i_) * xi_(i_)
      end do
    end do
    do n_ = 0, maxpow
      cosn_(n_) = cos(real(n_, rk) * q_(6))
    end do

    do ic_ = 1, 3

      mu_mol_(ic_) = 0.0_rk

      do k_ = 1, npar(ic_)
        term_ = coeff(k_, ic_)                      &
                * xpow_(ipow(1, k_, ic_), 1)        &
                * xpow_(ipow(2, k_, ic_), 2)        &
                * xpow_(ipow(3, k_, ic_), 3)        &
                * xpow_(ipow(4, k_, ic_), 4)        &
                * xpow_(ipow(5, k_, ic_), 5)        &
                * cosn_(ipow(6, k_, ic_))
        mu_mol_(ic_) = mu_mol_(ic_) + term_
      end do

    end do

    ! mu_y is antisymmetric in tau: overall factor sin(tau)
    mu_mol_(2) = mu_mol_(2) * sin(q_(6))

  end subroutine eval_DMS

  !-------------------------------------------------------------------
  ! cross: vector (cross) product a x b
  !-------------------------------------------------------------------
  function cross(a, b) result(c)
    real(rk), intent(in) :: a(3), b(3)
    real(rk) :: c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

  !-------------------------------------------------------------------
  ! perpendicular: a unit vector perpendicular to the unit vector u
  !-------------------------------------------------------------------
  function perpendicular(u) result(p)
    real(rk), intent(in) :: u(3)
    real(rk) :: p(3), t(3)
    integer  :: imin
    imin = minloc(abs(u), dim=1)       ! Cartesian axis least aligned with u
    t = 0.0_rk
    t(imin) = 1.0_rk
    p = cross(u, t)
    p = p / norm2(p)
  end function perpendicular

  !-------------------------------------------------------------------
  ! next_data_line: next non-blank, non-comment line from stdin
  !                 (eof = .true. at end of file)
  !-------------------------------------------------------------------
  subroutine next_data_line(buf, eof)
    character(len=*), intent(out) :: buf
    logical,          intent(out) :: eof
    integer :: ios, j
    eof = .false.
    do
      read(*,'(a)', iostat=ios) buf
      if (ios /= 0) then
        eof = .true.
        return
      end if
      ! tolerate tabs and Windows (CRLF) line endings
      do j = 1, len_trim(buf)
        if (buf(j:j) == char(9) .or. buf(j:j) == char(13)) buf(j:j) = ' '
      end do
      buf = adjustl(buf)
      if (len_trim(buf) == 0) cycle
      if (buf(1:1) == '!')    cycle
      return
    end do
  end subroutine next_data_line

  !-------------------------------------------------------------------
  ! is_end: .true. if the line starts with the keyword END
  !-------------------------------------------------------------------
  logical function is_end(buf)
    character(len=*), intent(in) :: buf
    character(len=3) :: w
    integer :: j
    w = buf(1:3)
    do j = 1, 3
      if (w(j:j) >= 'a' .and. w(j:j) <= 'z') w(j:j) = achar(iachar(w(j:j)) - 32)
    end do
    is_end = (w == 'END')
  end function is_end

  !-------------------------------------------------------------------
  ! read_real: read one real value from stdin, skipping comment lines
  !-------------------------------------------------------------------
  function read_real() result(val)
    real(rk) :: val
    character(len=512) :: buf
    logical :: eof_
    call next_data_line(buf, eof_)
    if (eof_) call input_error('non-linear parameters')
    read(buf, *) val
  end function read_real

  !-------------------------------------------------------------------
  ! read_coeff_line: read one "n1 n2 n3 n4 n5 n6  coeff" line,
  !                  skipping comment lines
  !-------------------------------------------------------------------
  subroutine read_coeff_line(ipow_, c_)
    integer,  intent(out) :: ipow_(ndeg)
    real(rk), intent(out) :: c_
    character(len=512) :: buf
    logical :: eof_
    call next_data_line(buf, eof_)
    if (eof_) call input_error('linear parameters')
    read(buf, *) ipow_(1), ipow_(2), ipow_(3), &
                 ipow_(4), ipow_(5), ipow_(6), c_
  end subroutine read_coeff_line

  !-------------------------------------------------------------------
  ! input_error: premature end of input
  !-------------------------------------------------------------------
  subroutine input_error(what)
    character(len=*), intent(in) :: what
    write(*,'(a,a)') 'ERROR: unexpected end of input while reading ', what
    stop
  end subroutine input_error

end program HOPO_DMS
