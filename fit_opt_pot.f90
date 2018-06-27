!======================================================================!
!     Title  : fit_opt_pot.f90                                         !
!     Author : Yusa Shusaku                                            !
!     date   : 2008-7-3-Thu ~ 2008-7-4-Fri                             !
!     Last modified : 2008-9-30-Tue                                    !
!                                                                      !
!     A program which computes a set of good parameters of the         !
!     potential by kai square fitting.                                 !
!     For the fitting, I adopted the Oak Ridge and Oxford method.      !
!     In this method, we have to solve coupled linear equation to find !
!     the increment of parameters. We use conjugate gradient method    !
!     for this purpose.                                                !
!     If the difference of kai squares with new and old parameters     !
!     becomes sufficiently small, we terminate the procedure. With     !
!     the obtainded parameters, we compute differential cross sections !
!     and compare them with experiment.                                !
!                                                                      !
!     *** Structure of this program ***                                !
!                                                                      !
!     module      Com_var                                              !
!     module      fitting_parameter                                    !
!     module      exp_data                                             !
!     module      Disp_pm                                              !
!       subroutine   Disp_parameters                                   !
!     program     main                                                 !
!     subroutine  Read_Data                                            !
!     subroutine  OROF_method                                          !
!     subroutine  Solve_dp                                             !
!     subroutine  Deriv                                                !
!     subrouitne  Parab_Approx                                         !
!     subroutine  Conj_Grad                                            !
!     subroutine  Guess_x                                              !
!     function    kai2                                                 !
!     subroutine  Phase_Shift                                          !
!     subroutine  Inn_Sol                                              !
!     subroutine  Numerov                                              !
!     function    kk                                                   !
!     subroutine  Disp_diff_CS                                         !
!     function    dsigma                                               !
!     function    dsc                                                  !
!     function    Vc                                                   !
!     function    V                                                    !
!     function    W                                                    !
!     function    Vopt                                                 !
!     function    FD                                                   !
!     function    PL                                                   !
!     subroutine  Gnuplot                                              !
!     subroutine  Show_Matrix                                          !
! *** Coulomb wave function ***                                        !
!     module Coulomb                                                   !
!       subroutine  jflgam                                             !
!       subroutine  yfclen                                             !
!       subroutine  yfasym                                             !
!       subroutine  dfcoul                                             !
!       subroutine  yfireg                                             !
!       subroutine  yfrica                                             !
!       subroutine  dfcz0                                              !
!       subroutine  jfdelg                                             !
!======================================================================!
      module fitting_parameter
!----------------------------------------------------------------------!
!     Definition of a derived data type.                               !
!     This is used for fitting parameters.                             !
!----------------------------------------------------------------------!
      implicit none
      type fitp
         real*8 :: V0, W0, a1, a2, R1, R2
      end type
      end module
!======================================================================!
      module exp_data
!----------------------------------------------------------------------!
!     Definition of a derived data dype.                               !
!     This is used for stocking experimental data and Rutherford cross !
!     sections of measured angles.                                     !
!     th   ... Measured angle in radian.                               !
!     ds   ... Experimental differential cross section.                !
!     er   ... Absolute error of the cross section.                    !
!     Ruth ... Rutherford cross section of an angle 'th'.              !
!----------------------------------------------------------------------!
      implicit none
      type xdata
         real*8 :: th, ds, er, Ruth
      end type
      end module
!======================================================================!
      module Disp_pm
      contains
      subroutine  Disp_parameters(pms)
      use Com_var, only : pn
      implicit none
      integer :: i
      real*8, intent(in) :: pms(pn)
      character(len=30), parameter :: FM='(1x,a,f10.4,a)'

      write(6,*)
      write(6,FM) 'V0 =', pms(1), ' MeV'
      write(6,FM) 'W0 =', pms(2), ' MeV'
      write(6,FM) 'a1 =', pms(3), ' fm'
      write(6,FM) 'a2 =', pms(4), ' fm'
      write(6,FM) 'R1 =', pms(5), ' fm'
      write(6,FM) 'R2 =', pms(6), ' fm'
      write(6,*)

      return
      end subroutine
      end module
!======================================================================!
      program main
!----------------------------------------------------------------------!
!     Main program.                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, Apro, Atar, pn, Maxdata
      use fitting_parameter
      use exp_data
      use Disp_pm
      implicit none
      integer :: n
      real*8, dimension(Maxdata) :: th, ds, er, Ruth
      real*8, external :: kai2
      real*8 :: sigma(0:Lc), pms(pn), dkai, z, tmp
      complex*16 :: Sn(0:Lc)
      type(fitp) :: pm
      type(xdata), allocatable :: xd(:)
      character(len=30), parameter :: FM='(1x,a,f12.3)'

      pm = fitp(-10.0d0, -10.0d0, 1.0d0, 1.0d0, 6.0d0, 6.0d0)
      pms(1) = pm%V0  ;  pms(2) = pm%W0
      pms(3) = pm%a1  ;  pms(4) = pm%a2
      pms(5) = pm%R1  ;  pms(6) = pm%R2

      call Read_Data(n, th, ds, er, Ruth, 'experiment.dat')
      allocate(xd(n))
      xd(1:n)%th = th(1:n)  ;  xd(1:n)%ds   = ds(1:n)
      xd(1:n)%er = er(1:n)  ;  xd(1:n)%Ruth = Ruth(1:n)

      write(6,*)
      write(6,*) '** Before fitting **'
      write(6,FM) 'kai square =', kai2(pms, n, xd)
      call  Disp_parameters(pms)
      call  Phase_Shift(Sn, sigma, pms)
      call  Gnuplot(sigma, Sn, 7, 'gnubefore', '12C_on_28Si_bf.eps')
      call  OROF_method(n, xd, pms)
      write(6,*)
      write(6,*) '** After  fitting **'
      write(6,FM) 'kai square =', kai2(pms, n, xd)
      call  Disp_parameters(pms)
      call  Phase_Shift(Sn, sigma, pms)
      call  Disp_diff_CS(sigma, Sn)
      call  Gnuplot(sigma, Sn, 8, 'gnuopt', '12C_on_28Si.eps')

      stop
      end program
!======================================================================!
      subroutine  Read_Data(n, th, ds, er, Ruth, Fname)
!----------------------------------------------------------------------!
!     A program which imports the experimental data from a file.       !
!     We import angles, scattering differential cross sections and     !
!     its errors.  We convert unit of angles from degree to radian     !
!     and errors from relative errors to absolute errors.              !
!     Not only these experimental data, we also compute Rutherford     !
!     cross sections of the measured angles.                           !
!----------------------------------------------------------------------!
      use Com_var, only : PI, Lc, Maxdata
      implicit none
      integer :: k, io, n
      real*8, dimension(Maxdata), intent(out) :: th, ds, er, Ruth
      real*8, external :: dsc
      character(len=*), intent(in) :: Fname

      open(unit=7, action='read', file=Fname)
      do k=1, Maxdata
         read(7,*,iostat=io) th(k), ds(k), er(k)
         if (io < 0) then
            n = k - 1
            exit
         end if
      end do
      close(7)
      if (k == Maxdata) then
         stop 'Too small size of array for the experimental data.'
      end if
      th(1:n) = PI / 180.0d0 * th(1:n)
      er(1:n) = er(1:n) * ds(1:n) * 0.01d0
!$omp parallel do
      do k=1, n
         Ruth(k) = dsc(th(k))
      end do
!$omp end parallel do

      return
      end subroutine
!======================================================================!
      subroutine OROF_method(n, xd, pms)
!----------------------------------------------------------------------!
!     A subroutine which executes Oak Ridge and Oxfort method.         !
!     We compute increment of parameters by subroutine 'Solve_dp' and  !
!     compute difference of kai squares with new and old parameters.   !
!     If the difference becomes sufficiently small, we terminate the   !
!     iteration.                                                       !
!----------------------------------------------------------------------!
      use Com_var, only : pn, MaxIter, epsr
      use fitting_parameter
      use exp_data
      implicit none
      integer, intent(in) :: n
      integer :: i, k
      real*8, intent(inout) :: pms(pn)
      real*8, external :: kai2
      real*8 :: dp(pn), z, pms0(pn), kai0, dkai
      type(xdata), intent(in) :: xd(n)
      character(len=11), parameter :: FM='(1x,a,i5)'

      do i=1, MaxIter
         pms0 = pms
         call Solve_dp(dp, pms, n, xd)
         kai0 = kai2(pms0, n, xd)
         dkai = (abs(kai2(pms, n, xd) - kai0)) / kai0
         if (dkai < epsr) then
            write(6,FM) 'Iteration :',i
            exit
         end if
      end do
      if (i==MaxIter .or. kai0>10.0d0) print *,  'We failed to fitting.'

      return
      end subroutine
!======================================================================!
      subroutine  Solve_dp(dp, pms, n, xd)
!----------------------------------------------------------------------!
!     This subroutine returns the increment of fitting parameters      !
!     'dp(pn)'.                                                        !
!----------------------------------------------------------------------!
      use Com_var, only : pn
      use exp_data
      implicit none
      integer, intent(in) :: n
      integer :: i
      real*8, intent(in)  :: pms(pn)
      real*8, intent(out) :: dp(pn)
      real*8 :: dF(pn), ddF(pn,pn), phi, z
      type(xdata), intent(in) :: xd(n)

      call  Deriv(n, xd, pms, dF, ddF)
      call  Conj_Grad(pn, ddF, -dF, dp)
      phi = maxval(dp / pms)
      if (phi < 0.1d0) then
         z = 1.0d0
      else if (phi < 0.3d0) then
         z = 0.5d0
      else if (phi < 0.5d0) then
         z = 0.2d0
      else if (phi < 1.0d0) then
         z = 0.1d0
      else if (phi >= 1.0d0) then
         z = 0.05d0 / phi
      end if
      call  Parab_Approx(dp, z, pms, n, xd)

      return
      end subroutine
!======================================================================!
      subroutine  Deriv(n, xd, pms0, dF, ddF)
!----------------------------------------------------------------------!
!     This subroutine computes the first and second derivatives of     !
!     kai square(F) as a function of fitting parameters.               !
!     We increment the parameters by one percent and compute the       !
!     derivatives of the cross secion.                                 !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, PI, pn
      use exp_data
      implicit none
      integer, intent(in) :: n
      integer :: i, j, k
      real*8, intent(in) :: pms0(pn)
      real*8, intent(out) :: dF(pn), ddF(pn,pn)
      real*8, external :: dsigma, dsc
      real*8 :: pms(pn), kai(n), S, sigma(0:Lc), rf(pn,n), tds(n)
      complex*16 :: Sn(0:Lc)
      type(xdata) :: xd(n)

      ddF = 0.0d0

      call Phase_Shift(Sn, sigma, pms0)
!$omp parallel do
      do k=1, n
         tds(k) = dsigma(Sn, sigma, xd(k)%th) / xd(k)%Ruth
      end do
!$omp end parallel do 

      pms = pms0
      do i=1, pn
         if (i /= 1) pms(i-1) = pms0(i-1)
         pms(i) = 1.01d0 * pms0(i)
         call  Phase_Shift(Sn, sigma, pms)
!$omp parallel do 
         do k=1, n
            kai(k) = (tds(k) - xd(k)%ds) / xd(k)%er ** 2
            rf(i,k) = dsigma(Sn, sigma, xd(k)%th) / xd(k)%Ruth - tds(k)
         end do
!$omp end parallel do 
         dF(i) = 2.0d2 * Dot_product(kai(1:n), rf(i,1:n)) / pms0(i)
         do j=1, i
            S = 0.0d0
!$omp parallel do reduction(+:S)
            do k=1, n
               S = S + rf(i,k) * rf(j, k) / xd(k)%er ** 2
            end do
!$omp end parallel do 
            ddF(i,j) = 2.0d4 * S / (pms0(i) * pms0(j))
         end do
      end do
      ddF = ddF + transpose(ddF)

      return
      end subroutine
!====================================================================!
      subroutine  Parab_Approx(dp, z, p0, n, xd)
!--------------------------------------------------------------------!
!     This subroutine computes a factor 'z' which is multipled to    !
!     the increment of fitting parameters 'dp(pn)' and returns new   !
!     fitting parameters. The factor 'z' is computed according to    !
!     the parabolic approximation.                                   !
!--------------------------------------------------------------------!
      use Com_var, only : pn
      use exp_data
      implicit none
      integer, intent(in) :: n
      real*8, intent(inout) :: p0(pn), z
      real*8, intent(in) :: dp(pn)
      real*8, dimension(pn) :: pa, pb
      real*8, external :: kai2
      real*8 :: y0, ya, yb, z0, za, zb, zm, ym, w, concave
      type(xdata), intent(in) :: xd(n)

      y0 = kai2(p0, n, xd)
      pa = p0 + z * dp
      ya = kai2(pa, n, xd)
      if (ya <= y0) then
         za = z
         zb = 2.0d0 * z
      else
         za = 0.5d0 * z
         zb = z
      end if
      pa = p0 + za * dp     ;  pb = p0 + zb * dp   
      ya = kai2(pa, n, xd)  ;  yb = kai2(pb, n, xd) 
      concave = (yb - y0) * za + (y0 - ya) * zb
      if (concave > 0.0d0) then
         w  = 3.0d0 * y0 - 4.0d0 * ya + yb
         zm = 0.5d0 * za * w / (y0 - 2.0d0 * ya + yb)
         if (zm <= - zb) then
            z = - zb
         else if (zm <= zb) then
            z = zm
         else if (zm > zb) then
            z = 2.0d0 * zb
         end if
      else
         if (yb < y0) then
            z = 3.0d0 * za
         else if (yb >= y0) then
            z = - za
         end if
      end if
      p0 = p0 + z * dp

      return
      end subroutine
!======================================================================!
      subroutine  Conj_Grad(N, a, b, x)
!----------------------------------------------------------------------!
!     A subroutine which solves coupled linear equation by conjugate-  !
!     gradient method.                                                 !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: N
      integer :: i, Step
      real*8, intent(in) :: a(N,N), b(N)
      real*8, intent(inout) :: x(N)
      real*8, dimension(N) :: r, p, q
      real*8 :: C1, C2, alpha, beta

      call  Guess_x(N, x)
      r = b - Matmul(a, x)
      p = r
      C2 = Dot_Product(x, b + r)

      Step = 0
      do 
          Step = Step + 1
          C1 = C2
          q = Matmul(a, p)
          alpha = Dot_product(p, r) / Dot_Product(p, q)
          x = x + alpha * p
          if (Mod(Step, 3) == 0) then
              r = b - Matmul(a, x)
          else
              r = r - alpha * q
          end if
          C2 = Dot_Product(x, b + r)
          if (C1 >= C2) exit
          beta = - Dot_Product(r, q) / Dot_Product(p, q)
          p = r + beta * p
      end do
          
      return
      end subroutine
!======================================================================!
      subroutine  Guess_x(N, x)
!----------------------------------------------------------------------!
!     A subroutine which gives the initial guess of the solutions of   !
!     coupled linear equation.                                         !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: N
      integer :: i
      real*8, intent(out) :: x(N)

      do i=1 ,N
          x(i) = 1.0d0
      end do

      return
      end subroutine
!====================================================================!
      function  kai2(pms, n, xd)  result(kai)
!--------------------------------------------------------------------!
!     A subroutine which computes kai square from the experimental   !
!     data and given parameters.                                     !
!--------------------------------------------------------------------!
      use Com_var, only : Lc, PI, pn
      use exp_data
      implicit none
      integer, intent(in) :: n
      integer :: k
      real*8, intent(in) :: pms(pn)
      real*8, external :: dsigma, dsc
      real*8 :: tds(n), t(n), error(n), kai, sigma(0:Lc)
      complex*16 :: Sn(0:Lc)
      type(xdata), intent(in) :: xd(n)

      call  Phase_Shift(Sn, sigma, pms)

      kai = 0.0d0
!$omp parallel do
      do k=1, n
         tds(k) = dsigma(Sn, sigma, xd(k)%th) / dsc(xd(k)%th)
      end do
!$omp end parallel do
!$omp parallel do reduction(+:kai)
      do k=1, n
         kai = kai + ((tds(k) - xd(k)%ds) / xd(k)%er) ** 2
      end do
!$omp end parallel do
      kai = kai / dble(n)

      return
      end function
!======================================================================!
      subroutine  Phase_Shift(Sn, sigma, pms)
!----------------------------------------------------------------------!
!     A subroutine which computes Coulomb phase shifts and nuclear     !
!     S-matrix. Since nuclear phase shifts are complex value, it is    !
!     difficult to compute them (the intrinsic function 'atan' is not  !
!     allowed to use complex values). Therefore, we compute the        !
!     nuclear S-matrix instead of the phase shifts.                    !
!     To compute the nuclear S-matrix, we need to know the values of   !
!     Coulomb wave functions at two points which are larger than the   !
!     range of nuclear potential. We use a subroutine 'dfcoul' to do   !
!     it. The Coulomb phase shifts are also obtained from the          !
!     subroutine.                                                      !
!                                                                      !
!     Fcw  --- Regular Coulomb wave function.                          !
!     Fpcw --- Derivative of Regular Coulomb wave function.            !
!     Gcw  --- Irreglar Coulomb wave function.                         !
!     Gpcw --- Derivative of Irregular Coulomb wave function.          !
!     Hp   --- Outgoing Coulomb wave function.                         !
!     Hm   --- Ingoing Coulomb wave function.                          !
!     sigma --- Coulomb phase shift.                                   !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, k, eta, pn
      use Coulomb, only : dfcoul
!$    use omp_lib
      implicit none
      integer :: L, iexp(0:Lc)
      real*8, intent(in) :: pms(pn)
      real*8, intent(out) :: sigma(0:Lc)
      real*8, dimension(0:Lc) :: Fcwa, Fpcwa, Gcwa, Gpcwa
      real*8, dimension(0:Lc) :: Fcwb, Fpcwb, Gcwb, Gpcwb
      real*8 :: rho, a, b
      complex*16, intent(out) :: Sn(0:Lc)
      complex*16, dimension(0:Lc) :: Hpa, Hma, Hpb, Hmb
      complex*16, parameter :: i=(0.0d0, 1.0d0)
      complex*16 :: s

      call Inn_Sol(0, a, b, s, pms)
!--Coulomb wave function at 'ka'.
      rho = k * a
      call dfcoul(eta,rho,Fcwa,Fpcwa,Gcwa,Gpcwa,sigma,Lc,iexp)
!--Coulomb wave function at 'kb'.
      rho = k * b
      call dfcoul(eta,rho,Fcwb,Fpcwb,Gcwb,Gpcwb,sigma,Lc,iexp)
      Hpa = Gcwa + i * Fcwa
      Hma = Gcwa - i * Fcwa
      Hpb = Gcwb + i * Fcwb
      Hmb = Gcwb - i * Fcwb
      Sn(0) = (s * Hmb(0) - Hma(0)) / (s * Hpb(0) - Hma(0))

!$omp parallel do
      do L=1, Lc
         call Inn_Sol(L, a, b, s, pms)
         Sn(L) = (s * Hmb(L) - Hma(L)) / (s * Hpb(L) - Hpa(L))
      end do
!$omp end parallel do
      return
      end subroutine
!======================================================================!
      subroutine  Inn_Sol(L, a, b, s, pms)
!----------------------------------------------------------------------!
!     A subroutine which integrates schroedinger equation from         !
!     origin to rmax using Numerov method.                             !
!     After that we compute 's' which is necessary for calculating     !
!     phase shifts.                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : rmax, dx, pn
      implicit none
      integer, intent(in) :: L
      real*8, intent(in) :: pms(pn)
      real*8, intent(out) :: a, b
      real*8 :: x
      complex*16, intent(out) :: s
      complex*16 :: y(3), ya

      x = 1.0d-3
      y(3) = x ** (L + 1)
      y(2) = y(3) + dble(L + 1) * x ** L * dx

      x = x + 2.0d0 * dx
      do
         call Numerov(L, x, y, pms)
         if (x > rmax) exit
         x = x + dx
         y(3) = y(2)
         y(2) = y(1)
      end do

      a = x - 2.0d0 * dx
      ya = y(3)
      b = x
      s = ya / y(1)

      return
      end subroutine
!======================================================================!
      subroutine  Numerov(L, x, y, pms)
!----------------------------------------------------------------------!
!     A subroutine which computes 'y(x+dx)' from the values 'y(x)'     !
!     and 'y(x-dx).                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : dx, pn
      use fitting_parameter
      implicit none
      integer, intent(in) :: L
      real*8, intent(in) :: x, pms(pn)
      real*8, parameter :: wh=dx*dx/12.0d0
      complex*16, intent(inout) :: y(3)
      complex*16, external :: kk
      complex*16 :: ynm2, ynm1

      ynm1 = (2.0d0 - 10.0d0 * wh * kk(L,x-dx,pms)) * y(2)
      ynm2 = (1.0d0 + wh * kk(L,x-2.0d0*dx,pms)) * y(3)
      y(1) = (ynm1 - ynm2) / (1.0d0 + wh * kk(L,x,pms))

      return
      end subroutine
!======================================================================!
      function kk(L, x, pms)  result(k)
!----------------------------------------------------------------------!
!     Definition of function 'kk' which appears in the Schroedinger    !
!     equation when we write it in the following form :                !
!      d^2y                                                            !
!     ------ + kk * y = 0 .                                            !
!      dx^2                                                            !
!----------------------------------------------------------------------!
      use Com_var, only : rmass, hbarc, E, pn
      implicit none
      integer, intent(in) :: L
      real*8, intent(in) :: x, pms(pn)
      real*8, external :: Vc
      complex*16, external :: Vopt
      complex*16 :: k
      
      k = 2.0d0 * rmass * (E - Vopt(x,pms) - Vc(x)) &
     &     / (hbarc * hbarc) - dble((L * (L + 1))) / (x * x)

      return
      end function
!======================================================================!
      subroutine  Disp_Diff_CS(sigma, Sn)
!----------------------------------------------------------------------!
!     This subroutine displays elastic differential cross section on   !
!     the screen.                                                      !
!----------------------------------------------------------------------!
      use Com_var, only : PI, Lc
      implicit none
      integer :: i
      real*8, intent(in) :: sigma(0:Lc)
      real*8, external :: dsc, dsigma
      real*8 :: t
      complex*16, intent(in) :: Sn(0:Lc)
      character(len=30) :: FM1='(2x,i4,3x,a,5x,1pd11.4)'

      write(6,*) '  Angle | Diff. Cross Section'
      write(6,*) '---------------------------------'
      do i=10, 180, 10
         t = PI / 180.0d0 * dble(i)
         write(6,FM1) i, '|',dsigma(Sn,sigma,t)/dsc(t)
      end do
      write(6,*)

      return
      end subroutine
!======================================================================!
      function  dsigma(Sn, sigma, t)  result(ds)
!----------------------------------------------------------------------!
!     A function which returns differential cross section for the      !
!     elastic scattering at angle 't'.                                 !
!     To compute the Coulomb scattering amplitude, we use an exact     !
!     formula and we use partial wave expansion for the nuclear        !
!     scattering amplitude.                                            !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, k, eta
      implicit none
      integer :: L
      real*8, intent(in) :: t, sigma(0:Lc)
      real*8, external :: PL, dsc
      real*8 :: ds
      complex*16, intent(in) :: Sn(0:Lc)
      complex*16, parameter :: i=(0.0d0, 1.0d0)
      complex*16 :: fn, fc, Sc(0:Lc)

      fn = (0.0d0, 0.0d0)
      do L=Lc, 0, -1
         Sc(L) = exp(2.0d0*i*sigma(L))
      end do
      do L=Lc, 0, -1
         fn = fn + dble(2*L + 1) * (Sn(L) - 1.0d0) * PL(L,cos(t))*Sc(L)
      end do
      fn = - 0.5d0 * i / k * fn
      fc = - 0.5d0 * eta / (k * sin(0.5d0*t) ** 2) &
     &      * exp(- 2.0d0*i*eta*log(sin(0.5d0*t)) + 2.0d0*i*sigma(0))
      ds = abs(fc + fn) ** 2

      return
      end function
!======================================================================!
      function  dsc(t) result(Ruthfd)
!----------------------------------------------------------------------!
!     Definition of Rutherford differential cross section.             !
!----------------------------------------------------------------------!
      use Com_var, only : eta, k
      implicit none
      real*8, intent(in) :: t
      real*8 :: Ruthfd

      Ruthfd = (0.5d0 * eta / (k * (sin(0.5d0 * t) ** 2))) ** 2

      return
      end function
!======================================================================!
      function  Vc(x)  result(V)
!----------------------------------------------------------------------!
!     Definition of Coulomb potential. We are assuming uniform charge  !
!     distribution.                                                    !
!----------------------------------------------------------------------!
      use Com_var, only : Zpro, Ztar, hbarc, Rc
      implicit none
      real*8, intent(in) :: x
      real*8 :: V

      if (x <= Rc) then
         V = 0.5d0 * Zpro * Ztar * hbarc / (137.0d0 * Rc) &
     &       * (3.0d0 - (x / Rc) ** 2)
      else
         V = Zpro * Ztar * hbarc / (137.0d0 * x)
      end if

      return
      end function
!======================================================================!
      function  V(x, pms)  result(ReOpt) 
!----------------------------------------------------------------------!
!     Definition of the real part of the optical potential.            !
!----------------------------------------------------------------------!
      use Com_var, only : pn
      implicit none
      real*8, intent(in) :: x, pms(pn)
      real*8, external :: FD
      real*8 :: ReOpt

      ReOpt = pms(1) * FD(x, pms(3), pms(5))

      return
      end function
!======================================================================!
      function  W(x, pms) result(ImOpt)
!----------------------------------------------------------------------!
!     Definition of the imaginary part of the Optical potential.       !
!----------------------------------------------------------------------!
      use Com_var, only : pn
      implicit none
      real*8, intent(in) :: x, pms(pn)
      real*8, external :: FD
      real*8 :: ImOpt

      ImOpt = pms(2) * FD(x, pms(4), pms(6))

      return
      end function
!======================================================================!
      function  Vopt(x, pms)  result(Opt)
!----------------------------------------------------------------------!
!     Definition of optical potential.                                 !
!----------------------------------------------------------------------!
      use Com_var, only : pn
      implicit none
      real*8, intent(in) :: x, pms(pn)
      real*8, external :: V, W
      complex*16, parameter :: i=(0.0d0, 1.0d0)
      complex*16 :: Opt

      Opt = V(x,pms) + i * W(x,pms)

      return
      end function
!======================================================================!
      function  FD(x, a, R)  result(f)
!----------------------------------------------------------------------!
!     Definition of Fermi-Dirac distribution function.                 !
!----------------------------------------------------------------------!
      implicit none
      real*8, intent(in) :: x, R, a
      real*8 :: f

      f = 1.0d0 / (1.0d0 + exp((x - R) / a))

      return
      end function
!======================================================================!
      function  PL(L, x)  result(P)
!----------------------------------------------------------------------!
!     Definition of Legendre polynomials.                              !
!     We use recursion relation in computing it.                       !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real*8, intent(in) :: x
      real*8 :: P0, P1, P

      P0 = 1.0d0
      P1 = x

      if (L == 0) then
         P = P0
      else if (L == 1) then
         P = P1
      else
         do i=2, L
            P = (dble(2*i - 1) * x * P1 - dble(i - 1) * P0) / dble(i)
            P0 = P1
            P1 = P
         end do
      end if

      return
      end function
!======================================================================!
      subroutine  Gnuplot(sigma, Sn, Fnum, Fname, eps_name)
!----------------------------------------------------------------------!
!     This subroutine makes a gnuplot file.                            !
!     If we redirect it into the gnuplot, we can obtain eps file.      !
!----------------------------------------------------------------------!
      use Com_var, only : Lc, PI
      implicit none
      integer, intent(in) :: Fnum
      integer :: i
      real*8, intent(in) :: sigma(0:Lc)
      real*8, external :: dsc, dsigma
      real*8 :: t
      complex*16, intent(in) :: Sn(0:Lc)
      character(len=*), intent(in) :: Fname, eps_name
      character(len=17) :: FM='(1x,f8.3,f13.8)'
      character(len=35) :: CS1='(d{/Symbol s} / d{/Symbol W})'
      character(len=35) :: CS2='(d{/Symbol s}_c / d{/Symbol W})'

      open(Fnum, file=Fname, action='write')
      write(Fnum,*) 'set term postscript eps enhanced color'
      write(Fnum,*) "set output '" // eps_name // "'"
      write(Fnum,*) 'set title ''Differential cross section'''
      write(Fnum,*) 'set xlabel ''{/Symbol q}_{cm} (degree)'''
      write(Fnum,*) "set ylabel '" //Trim(CS1)//" / "//Trim(CS2)//"'"
      write(Fnum,*) 'set label ''^{12}C + ^{28}Si'' at 20, 0.7'
      write(Fnum,*) 'set label ''E_{lab} = 131.5 MeV'' at 20,0.2'
      write(Fnum,*) 'unset key'
      write(Fnum,*) 'set size 0.8, 0.8'
      write(Fnum,*) 'set logscale y'
      write(Fnum,*) 'set xrange[0:50]'
      write(Fnum,*) 'set yrange[0.000001:2]'
      write(Fnum,*) 'set format y "10^{%L}"'
      write(Fnum,*) 'set size ratio 1.5'
      write(Fnum,*) 'plot ''experiment.dat''using 1:2:($2*$3/100.0) \'
      write(Fnum,*) ' with yerrorbars 12,''-'' with line 1'

      do i=1, 900
         t = dble(i) * PI / 900.0d0
         write(Fnum,FM) 0.2d0 * dble(i), dsigma(Sn, sigma, t) / dsc(t)
      end do
      close(Fnum)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real*8, intent(in) :: a(N,N)
      character(len=30) :: FM

      FM = '(x,a,' // CHAR(48+N) // 'es9.1,a)'
      do i=1, N
          write(6,FM) ' |', (a(i,j), j=1, N), ' |'
      end do
      write(6,*)

      return
      end subroutine
