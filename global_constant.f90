!======================================================================!
      module Com_var
!----------------------------------------------------------------------!
! *** Definitions of various global constants ***                      !
!     dx    ... Step width of Numerov method.                          !
!     rmax  ... Range of nuclear potential.                            !
!     Apro  ... Mass number of projectile nucleus.                     !
!     Atar  ... Mass number of target nucleus.                         !
!     Zpro  ... Proton number of projectile nucleus.                   !
!     Ztar  ... Proton number of target nucleus.                       !
!     Elab  ... Bombering energy in laboratory system.                 !
!     E     ... Bombering energy in center of mass system.             !
!     rmass ... Reduced mass.                                          !
!     eta   ... Sommerfeld parameter.                                  !
!     k     ... Wave number of free particle.                          !
!     Lc    ... Cutoff angular momentum.                               !
!     Rc    ... Coulomb radius.                                        !
!                                                                      !
!     pn    ... The number of fitting parameters.                      !
!     MaxIter ... Maximum number of iteration.                         !
!     Maxdata ... Maximum number of importable experimental data.      !
!     epsr  ... Condition of convergence of fitting parameters.        !
!----------------------------------------------------------------------!
      implicit none
      real*8, parameter :: PI=3.1415926535897932d0, dx=0.05d0
      real*8, parameter :: rmax=15.0d0
      real*8, parameter :: Mn=938.0d0, hbarc=197.329d0
      real*8, parameter :: Apro=12.0d0, Atar=28.0d0
      real*8, parameter :: Zpro=6.0d0,  Ztar=14.0d0
      real*8, parameter, private :: Elab=131.5d0
      real*8, parameter :: E=Atar*Elab/(Apro+Atar)
      real*8, parameter :: rmass=Mn*Apro*Atar/(Apro+Atar) 
      real*8, parameter :: eta=Zpro*Ztar*sqrt(0.5d0*rmass/E)/137.0d0
      real*8, parameter :: k=sqrt(2.0d0*E*rmass)/hbarc
      integer :: Lc
      parameter(Lc=k*rmax*(sqrt(1.0d0-2.0d0*eta/(k*rmax))) + 1.0d0)
      real*8, private :: Rsum
      parameter(Rsum=Apro**(1.0d0/3.0d0)+Atar**(1.0d0/3.0d0))
      real*8, parameter :: Rc=1.2d0*Rsum
!-----
      integer, parameter :: pn=6, MaxIter=1000, Maxdata=200
      real*8, parameter :: epsr=1.0d-8
      end module
