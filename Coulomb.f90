!======================================================================!
      module  Coulomb
      contains
!**********************************************************************!
!                                                                      !
!     *** Subroutine for Coulomb wave functions ***                    !
!                                                                      !
!     call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,L,iexp)             !
!                                                                      !
!     eta  --- Sommerfeld parameter. real*8, intent(in)                !
!     rho  --- Independent variable. real*8, intent(in)                !
!     fcw(L)  --- Regular Coulomb wave function of angular momentum L. !
!                 real*8, intent(out)                                  !
!     fpcw(L) --- Derivative of regular Coulomb wave function.         !
!                 real*8, intent(out)                                  !
!     gcw(L)  --- Irregular Coulomb wave function. real*8, intent(out) !
!     gpcw(L) --- Derivative of irregular Coulomb wave function.       !
!                 real*8, intent(out)                                  !
!     sigmad(L) --- Coulomb phase shift. Maybe... real*8, intent(out)  !
!     L --- Angular momentum. We can compute above quantities up to    !
!           this angular momentum. integer, intent(in)                 !
!     iexp(L) --- ??? integer                                          !
!----------------------------------------------------------------------!
      subroutine jflgam(xd,yd,par,pai,nbchif)
      real*8 :: xd, yd, par, pai, test, c, hlo2pi, pi, pis2, pis4
      real*8 :: x, y, u, v, tra, tra1, alo2pi, rac2, cosi, sini
      real*8 :: cos2i, sin2i, zmod, depi
      real*8 :: xx
      real*8 :: alopi
      real*8 :: supint
      dimension test(7),c(6)
      data test/2.9152d+7,2.2958d+3,1.4124d+2,3.9522d+1,19.6611d0,  &
     &12.791d0,-10.d0/
      data c/8.333333333333333d-2,-2.777777777777777d-3, &
     &7.936507936507937d-4,-5.952380952380952d-4,8.417508417508418d-4, &
     &-1.917526917526918d-3/
      data hlo2pi/0.918938533204672/,pi/3.141592653589793/
      data pis2/1.570796326794897/,pis4/0.785398163397448/
      data alo2pi/1.837877066409345/,rac2/0.3465735902799726/
      data depi/6.283185307179586/,alopi/1.1447298858494002/
      data supint/2147483647.d0/
      
      nbchif = 15
      x = dabs(xd)
      xx = x
      if (yd) 1,2,1
    1 y = dabs(yd)
      kr = 1
      i = dmod(10.99d0-x,supint)
!     translation
      if (i) 3,3,4
    4 tra = i
      x = x + tra
!     logarithme(x+iy) (x,y positifs)
    3 if (x - y) 5,6,7
    5 tra1 = x / y
      if (tra1) 8,8,9
    8 u = dlog(y)
      v = pis2
      go to 10
    6 u = rac2 + dlog(x)
      v = pis4
      go to 10
    9 tra = y * dsqrt(1.d0 + tra1 * tra1)
      tra1 = y / x
   11 u = dlog(tra)
      v = datan(tra1)
   10 go to (12,19,23),kr
    7 tra1 = y / x
      tra=x*dsqrt(1.d0+tra1*tra1)
      go to 11
   12 kr=2
!     developpement asymptotique ( x superieur a 10 )
      tra=x-0.5d0
      pai=v*tra+y*(u-1.d0)
      par=-x+hlo2pi+u*tra-v*y
      zmod=x+y
      if(zmod-test(1))13,13,14
   13 tra=x*x+y*y
      cosi=x/tra
      sini=y/tra
      sin2i=(sini*cosi)+(sini*cosi)
      cos2i=(cosi+sini)*(cosi-sini)
      k=1
      go to 15
   16 tra=cosi*cos2i-sini*sin2i
      sini=sini*cos2i+cosi*sin2i
      cosi=tra
   15 par=par+c(k)*cosi
      pai=pai-c(k)*sini
      k=k+1
      if(zmod-test(k))16,16,14
!     translation inverse
   17 i=i-1
      x=i
      x=xx+x
      go to 3
   19 par=par-u
      pai=pai-v
   14 if(i-1)18,60,17
   60 if(xd)17,61,17
!     controle du quadrant
   18 if(xd)20,61,21
   61 tra=pi*y
      if(tra-1.d-2)300,300,301
  300 par= tra*(2.d0+tra*(-2.d0+tra*(1.333333333333333d0+tra*(&
     &-0.6666666666666666d0+tra*(0.2666666666666666d0+tra*(   &
     &-0.08888888888888888d0+tra*0.02539682539682540d0))))))
      tra1=-dlog(y)-dlog(par)
      go to 302
  301 par=1.d0-dexp(-tra-tra)
      tra1=-dlog(y*par)
  302 par=0.5d0*(alo2pi-tra+tra1)
      pai=pai-pis2
   21 if(yd)28,100,100
!     x+iy change en -x-iy
   20 tra=pi*y
      par=alo2pi-u-par-tra
      pai=pi-v-pai
      tra=dexp(-tra-tra)
      x=pi*dmod(x,2.d0)
      sini=(1.d0-tra)*dcos(x)
      cosi=(1.d0+tra)*dsin(x)
      kr=3
      x=dabs(cosi)
      y=dabs(sini)
      go to 3
   23 if(cosi)24,25,25
   24 v=pi-dsign(v,sini)
      go to 26
   25 if(sini)27,26,26
   27 v=-v
   26 par=par-u
      pai=pai-v
      if(yd)100,100,28
   28 pai=-pai
!     argument dans -pi,pi
  100 tra=dabs(pai/depi)
      if(tra-1.d+15)203,204,204
  204 nbchif=0
      pai=0.d0
      go to 29
  203 if(tra-1.d0)205,205,206
  206 nbchif=dlog10(tra)
      nbchif=14-nbchif
      tra=dmod(tra,supint)
      pai=dmod(tra,1.d0)*dsign(depi,pai)
  205 if(dabs(pai)-pi)29,29,207
  207 pai=pai-dsign(depi,pai)
      go to 29
!     jflgam reel
    2 pai=0.d0
      if(xd)31,32,33
!     conditions d existence
   32 write (6,1000)
 1000 format (21h jflgam(0) est infini)
      go to 50
   31 if(x-4503599627370496.d0)34,35,35
   35 write (6,1001)
 1001 format (30h argument de jflgam trop grand)
      go to 50
   34 y=dmod(x,supint)
      if(y)400,36,400
  400 if(y-0.99d0)33,33,405
  405 tra=idint(y+0.1d0)
      if(dabs(y-tra)-5.d-15)36,36,33
   36 write (6,1002)
 1002 format (28h jflgam (-entier) est infini)
   50 par=1.d+74
      nbchif=0
      go to 29
!     translation
   33 i=dmod(10.99d0-x,supint)
      if(i)37,37,38
   38 tra=i
      x=x+tra
!     developpement asymptotique
   37 y=dlog(x)
      par=-x+hlo2pi+y*(x-0.5d0)
      if(x-test(1))39,39,43
   39 cosi=1.d0/x
      cos2i=cosi*cosi
      k=1
      go to 41
   42 cosi=cosi*cos2i
   41 par=par+c(k)*cosi
      k=k+1
      if(x-test(k))42,42,40
!     translation inverse
   44 x=x-1.d0
   48 y=dlog(x)
      par=par-y
      i=i-1
   40 if(i-1)43,49,44
   49 x=dabs(xd)
      go to 48
!     x negatif
   43 if(xd)45,29,29
   45 par=alopi-par-y
      y=pi*dmod(x,2.d0)
      y=-dsin(y)
      if(y)46,36,47
   46 y=-y
      pai=pi
   47 par=par-dlog(y)
      entry jflgv1
   29 return
      end subroutine

!====================================================================!
      subroutine yfclen(eta,ro,u,up,v,vp,sigma0,idiv,nn)
      implicit complex*16(a-d,f-h),real*8(e,p-z)
!
      if(nn.eq.1)go to 20
!
      eta2=eta*eta
      fa=dcmplx(1.d0,eta)
      m=0.25*eta+4.d1
!
!          polynomes de tchebichev jusqu'au rang m
!
      k=m+2
      x=(eta+eta)/ro
      xx=x+x-1.d0
      t0=1.d0
      t1=xx
      xx=xx+xx
      do 6 j=2,k
      tj=xx*t1-t0
      t0=t1
      t1=tj
    6 continue
      tm=t1
      tl=t0
!
!          initialisation
!
      am=(0.d0,0.d0)
      al=(1.d-40,1.d-40)
      bn=(0.d0,0.d0)
      bm=(1.d-40,1.d-40)
      bl=(0.d0,0.d0)
      bk=dcmplx(4.d0*dfloat(m+1),0.d0)*al+bm
      f=(0.d0,0.d0)
      g =(0.d0,0.d0)
      gp=(0.d0,0.d0)
!
!          recurrence descendante
!
      k=m
   10 r=k
      tk=xx*tl-tm
      tm=tl
      tl=tk
      hk=dcmplx(tk,0.d0)
      c1=dcmplx(r*(r+1.d0)-eta2,eta*(r+r+1.d0))
      c2=(4.d0,0.d0)*dcmplx(r+1.d0,0.d0)
      c2=c2*dcmplx(-r-1.d0,eta*3.d0)
      c3=fa*dcmplx(-r-r-4.d0,eta)
      c4=dcmplx((7.d0*r+5.d0)/4.d0,0.d0)
      c5=dcmplx(r+r+2.d0,0.d0)
      c6=dcmplx((r+3.d0)/4.d0,0.d0)
      ak=(c2*al+c3*am-c4*bl-c5*bm-c6*bn)/c1
      j=k/2
      j=k-j-j
      if(j)1,2,1
    1 f=f-ak
      go to 3
    2 f=f+ak
    3 z=abs(ak)
      g=g+hk*ak
      gp=gp+hk*bk
!
!          f=a0/2-a1+a2-a3+a4-a5+...
!
!          congruence modulo 10**60
!
      if(z-1.d60)4,5,5
    5 d=(1.d60,0.d0)
      ak=ak/d
      al=al/d
      am=am/d
      bk=bk/d
      bl=bl/d
      bm=bm/d
      bn=bn/d
      f=f/d
      g=g/d
      gp=gp/d
    4 if(k)8,8,9
    9 d=dcmplx(4.d0*r,0.d0)
      bj=d*ak+bl
      am=al
      al=ak
      bn=bm
      bm=bl
      bl=bk
      bk=bj
      k=k-1
      go to 10
!
!          normalisation et calcul de z(ro)
!
    8 d=(0.5d0,0.d0)*ak
      f=f-d
      g=g-d
      gp=gp-(0.5d0,0.d0)*bk
      d=dcmplx(0.d0,-eta*dlog(2.d0)+sigma0)
      axpo=exp(d)
      f=f/axpo
      g=g/f
      gp=gp/f
!
!          calcul de f et g
!
      d=dcmplx(0.d0,ro-eta*dlog(ro))
      axpo=exp(d)
      d=g*axpo
      gp=axpo*(dcmplx(0.d0,1.d0-eta/ro)*g-dcmplx(x/ro,0.d0)*gp)
      v=d
      d=(0.d0,-1.d0)*d
      u=d
      vp=gp
      gp=(0.d0,-1.d0)*gp
      up=gp
      idiv=0
      return
!
!          serie origine
!
   20 pi=3.141592653589793d0
      xa=0.577215664901533d0
      ro2=ro*ro
      etap=eta+eta
      pieta=pi*eta
      z=138.15510557964276d0
      idiv=0
      if(dabs(pieta)-z)21,21,22
   22 indg=idint(pieta/z)
      idiv=60*indg
      if(eta.lt.0) idiv=0
      pieta=pieta-z*dfloat(indg)
   21 pieta2=0.5d0*pieta
      p=dexp(pieta2)*dsqrt(dsinh(pieta)/pieta)
      call jfdelg(1.d0,eta,par,pai,nb)
      z1=etap*(xa+xa+dlog(2.d0)-1.d0+par)
      u0=0.d0
      u1=ro
      v0=1.d0
      v1=z1*ro
      u=u0+u1
      v=v0+v1
      up=1.d0
      vp=z1
      xn=2.d0
      do 104 n=2,10000
      xn1=xn*(xn-1.d0)
      u2=(etap*ro*u1-ro2*u0)/xn1
      u=u+u2
      v2=(etap*ro*v1-ro2*v0-etap*(xn+xn-1.d0)*u2)/xn1
      v=v+v2
      up=up+xn*u2/ro
      vp=vp+xn*v2/ro
      if(dabs(u2/u).gt.1.d-14)goto 18
      if(dabs(v2/v).le.1.d-14)goto 19
   18 u0=u1
      u1=u2
      v0=v1
      v1=v2
      xn=xn+1.d0
  104 continue
   19 pp=v+etap*u*dlog(ro)
      w=u/p
      wp=up/p
      v=p*pp
      vp=p*(vp+etap*(up*dlog(ro)+u/ro))
      u=w
      up=wp
      return
      end subroutine

!====================================================================!
      subroutine yfasym(eta,rau,fo,fpo,go,gpo,sigo,iexp)
      implicit real*8 (a-h,o-z)
      iexp=0
      trb=0.d0
      rau2=rau+rau
      etac=eta*eta
      call jflgam(1.d0,eta,tra,sigo,ntruc)
   40 n=0
      ps=1.d0
      gs=0.d0
      pt=0.d0
      gt=1.d0-eta/rau
      sf=ps
      sg=gs
      spf=pt
      spg=gt
   45 denom= dfloat(n+1)*rau2
      an= dfloat(n+n+1)*eta/denom
      bn= (etac-dfloat(n*(n+1)))/denom
      ps1=an*ps-bn*pt
      gs1=an*gs-bn*gt-ps1/rau
      pt1=an*pt+bn*ps
      gt1=an*gt+bn*gs-pt1/rau
   42 sf=sf+ps1
      sg=sg+gs1
      spf=spf+pt1
      spg=spg+gt1
      n=n+1
      if(n-17)46,48,44
   48 tra=ps*ps+pt*pt
   44 trb=ps1*ps1+pt1*pt1
      test=tra-trb
      if(test)47,47,46
   46 ps=ps1
      gs=gs1
      pt=pt1
      gt=gt1
      tra=trb
      goto  45
   47 tetao= rau-eta*dlog (rau2)+sigo
      tra= dsin(tetao)
      trb=dcos(tetao)
      go=sf*trb-spf*tra
      gpo=sg*trb-spg*tra
      fo=spf*trb+sf*tra
      fpo=spg*trb+sg*tra
      return
      end subroutine

!====================================================================!
      subroutine dfcoul(eta,ro,f,fp,g,gp,sigma,l,iexp)
      implicit real*8(a-h,o-z)
      dimension f(*),fp(*),g(*),gp(*),iexp(*),sigma(*)
      data depi/6.283185307179586d0/
      etac=eta*eta
      l1=l+1
      call dfcz0(eta,ro,f0,fp0,g0,gp0,s,i)
      f(1)=f0
      fp(1)=fp0
      g(1)=g0
      gp(1)=gp0
      iexp(1)=i
      sigma(1)=s
      if(l)1,1,2
    1 return
    2 linf=0
      ind=0
      if((eta.gt.0).and.(ro.lt.(eta+eta)))go to 21
      z=eta+dsqrt(etac+dfloat(l*(l+1)))
      if(ro.ge.z)go to 20
    7 roinf=eta+dsqrt(etac+dfloat(linf*(linf+1)))
      if(ro-roinf)3,4,4
    4 if(linf-l)5,6,6
    5 linf=linf+1
      go to 7
    3 ind=1
    6 lin=linf+1
   20 xm=1.d0
      if(ind.eq.0)lin=l1
      do 8 j=2,lin
      zig=(dsqrt(etac+xm*xm))/xm
      zag=eta/xm+xm/ro
      f(j)=(zag*f(j-1)-fp(j-1))/zig
      fp(j)=zig*f(j-1)-zag*f(j)
      g(j)=(zag*g(j-1)-gp(j-1))/zig
      gp(j)=zig*g(j-1)-zag*g(j)
      iexp(j)=i
      sig=sigma(j-1)+datan(eta/(j-1))
      ipi=sig/depi
      sig=sig-ipi*depi
      if(sig)60,50,70
   60 if(sig.lt.(-depi/2.d0))sig=sig+depi
      go to 50
   70 if(sig.gt.(depi/2.d0))sig=sig-depi
   50 sigma(j)=sig
    8 xm=xm+1.d0
      if(ind.eq.0)return
      go to 22
   21 lin=1
   22 ftest=f(lin)
      fptest=fp(lin)
      lmax=linf+25+idint(5.d0*dabs(eta))
      if(lmax-l)9,10,10
    9 lmax=l
   10 fi=1.d0
      fpi=1.d0
   13 xm=lmax+1
      zig=(dsqrt(etac+xm*xm))/xm
      zag=eta/xm+xm/ro
      fl=(zag*fi+fpi)/zig
      fpl=zag*fl-zig*fi
      if(dabs(fl)-1.d15)26,27,27
   26 if(dabs(fpl)-1.d15)28,27,27
   27 fl=fl*1.d-15
      fpl=fpl*1.d-15
   28 fi=fl
      fpi=fpl
      if(lmax-l)11,11,12
   12 lmax=lmax-1
      go to 13
   11 f(lmax+1)=fl
      fp(lmax+1)=fpl
      if(lmax-linf)15,15,14
   14 go to 12
   15 fact=ftest/f(lin)
      factp=fptest/fp(lin)
      indice=i/60
      xm=linf
      do 16 j=lin,l1
      f(j)=f(j)*fact
      fp(j)=fp(j)*factp
   25 if(j.eq.1)go to 16
      zig=(dsqrt(etac+xm*xm))/xm
      zag=eta/xm+xm/ro
      g(j)=(zag*g(j-1)-gp(j-1))/zig
      gp(j)=zig*g(j-1)-zag*g(j)
      if(dabs(g(j))-1.d60)17,18,18
   17 if(dabs(gp(j))-1.d60)19,18,18
   18 g(j)=g(j)/1.d60
      gp(j)=gp(j)/1.d60
      indice=indice+1
   19 iexp(j)=indice*60
      a=fp(j)*g(j)
      b=-f(j)*gp(j)
      if(a-1.d0)29,30,30
   29 i1=idint(dlog10(a))
      i2=idint(dlog10(b))
      go to 31
   30 i1=idint(dlog10(a))+1
      i2=idint(dlog10(b))+1
   31 f(j)=f(j)*1.d1**(-i2)
      fp(j)=fp(j)*1.d1**(-i1)
      sig=sigma(j-1)+datan(eta/(j-1))
      ipi=sig/depi
      sig=sig-ipi*depi
      if(sig)61,51,71
   61 if(sig.lt.(-depi/2.d0))sig=sig+depi
      go to 51
   71 if(sig.gt.(depi/2.d0))sig=sig-depi
   51 sigma(j)=sig
   16 xm=xm+1.d0
      return
      end subroutine

!====================================================================!
      subroutine yfireg(eta,ro,g0,gp0)
      implicit real*8(a-h,o-z)
      if(eta.le.0.d0)goto 250
      if(eta.le.3.d0)goto 251
      if(eta.le.1.d1)goto 252
      if(eta.le.18.d0)goto 253
      if(eta.le.22.d0)goto 254
      if(ro.le.0.3d0+(3.d1-eta)/8.d1)goto 200
!   serie de taylor depart rau0
  300 continue
      rau0=1.666666666666667d0*dabs(eta)+7.5d0
      call yfasym(eta,rau0,f0,fp0,g0,gp0,sigma0,iexp)
      x=rau0-ro
      x2=x*x
      x3=x*x2
      unr=1.d0/rau0
      etr0=1.d0-2.d0*eta*unr
      u0=g0
      u1=-x*gp0
      u2=-0.5d0*etr0*x2*u0
      s=u0+u1+u2
      v1=u1/x
      v2=2.d0*u2/x
      t=v1+v2
      xn=3.d0
      do 10 n=3,10000
!     n=n
      xn1=xn-1.d0
      xn1=xn*xn1
      u3=x*u2*unr*(1.d0-2.d0/xn)-etr0*u1*x2/xn1+x3*u0*unr/xn1
      s=s+u3
      v3=xn*u3/x
      t=t+v3
   16 if(dabs(u3/s).gt.1.d-11)go to 11
      if(dabs(v3/t).le.1.d-11)go to 20
   11 u0=u1
      u1=u2
      u2=u3
      xn=xn+1.d0
   10 continue
   20 g0=s
      gp0=-t
      return
!   serie  origine
  200 continue
      pi=3.141592653589793d0
      ga=0.577215664901533d0
      eta2=eta*eta
      ro2=ro*ro
      etap=eta+eta
      pieta=pi*eta
      pieta2=0.5d0*pieta
      b=dexp(pieta2)*dsqrt(dsinh(pieta)/pieta)
      call jfdelg(1.d0,eta,par,pai,nb)
      c1=etap*(ga+ga+dlog(2.d0)-1.d0+par)
      u0=0.d0
      u1=ro
      v0=1.d0
      v1=c1*ro
      u=u0+u1
      v=v0+v1
      up=1.d0
      vp=c1
      xn=2.d0
      do 104 n=2,10000
      xn1=xn*(xn-1.d0)
      u2=(etap*ro*u1-ro2*u0)/xn1
      u=u+u2
      v2=(etap*ro*v1-ro2*v0-etap*(xn+xn-1.d0)*u2)/xn1
      v=v+v2
      up=up+xn*u2/ro
      vp=vp+xn*v2/ro
   17 if(dabs(u2/u).gt.1.d-14)goto 18
      if(dabs(v2/v).le.1.d-14)goto 19
   18 u0=u1
      u1=u2
      v0=v1
      v1=v2
      xn=xn+1.d0
  104 continue
   19 gp=v+etap*u*dlog(ro)
      g0=b*gp
      gp0=b*(vp+etap*(up*dlog(ro)+u/ro))
      return
  250 if(ro.le.0.5d0*eta+9.d0)goto 200
      goto  300
  251 if(ro.le.2.25d0+7.35d0*(3.d0-eta))goto 200
      goto  300
  252 if(ro.le.1.2d0+1.5d-1*(1.d1-eta))goto 200
      goto  300
  253 if(ro.le.0.6d0+0.75d-1*(18.d0-eta))goto 200
      goto  300
  254 if(ro.le.0.4d0+0.5d-1*(22.d0-eta))goto 200
      goto  300
      end subroutine

!====================================================================!
      subroutine yfrica(eta,rau,fo,fpo,go,gpo,sigma0,idiv)
      implicit real*8 (a-h,o-z)
      dimension q(5),qp(5)
!
!        coefficients riccati
!
      data g61,g62,g63,g64,g65,g66,g67,g68,g69,g610,  &
     &g611/ 0.1159057617187498d-01,0.3863525390624998d-01,  &
     &0.4660034179687498d-01,0.4858398437499998d-01,  &
     &0.1156514485677080d1,0.5687475585937496d1,  &
     &0.1323888288225445d2,0.1713083224826384d2,  &
     &0.1269003295898436d2,0.5055236816406248d1,  &
     &0.8425394694010415d0/
      data g81,g82,g83,g84,g85,g86,g87,g88,g89,g810,g811,g812,g813,g814,  &
     &g815/ 0.1851092066083633d-01,0.8638429641723630d-01,  &
     &0.1564757823944092d0,0.1430139541625977d0,  &
     &0.1924622058868408d0,0.8500803152720129d1,  &
     &0.7265429720878595d2,0.3057942376817972d3,  &
     &0.7699689544836672d3,0.1254157054424285d4,  &
     &0.1361719536066055d4,0.9831831171035763d3,  &
     &0.4547869927883148d3,0.1222640538215636d3,  &
     &0.1455524450256709d2/
      data gp61,gp62,gp63,gp64,gp65,gp66/ 0.2897644042968748d-01,  &
     &0.2318115234375000d0,0.8056640625000000d0, &
     &0.1601562499999998d1,0.3046875000000000d0,&
     &0.5624999999999998d1/
      data gp81,gp82,gp83,gp84,gp85,gp86,gp87, &
     &gp88/ 0.6478822231292720d-01,0.6910743713378906d0, &
     &0.3322952270507811d1,0.9483032226562498d1, &
     &0.1769653320312499d2,0.3478710937499998d2, &
     &0.5020312499999999d2,0.7874999999999999d2/
      data q /0.4959570165d-1,0.8888888889d-2,0.2455199181d-2, &
     &0.9108958061d-3,0.2534684115d-3/
      data qp /0.1728260369d0,0.3174603174d-3,0.3581214850d-2, &
     &0.3117824680d-3,0.9073966427d-3/
      call jflgam(1.d0,eta,tra,sigma0,ind)
      rau2=rau+rau
      rauc=rau*rau
      etac=eta*eta
      eta2=eta+eta
      etaro=eta*rau
      etaro2=etaro+etaro
      pieta=3.141592653589793*eta
      ind=0
      jnd=0
      ig=0
      if(eta)20,20,21
   20 if(-etaro-14.0625d0)32,22,22
   22 indice=1
 
!             riccati 3
      idiv=0
      go to 2
   21 if(dabs(rau-eta2).le.1.d-9)go to 18
      if(rau-eta2)30,18,31
   31 if(rau-eta2-2.d1*(eta**0.25d0))34,33,33
   33 indice=0
!             riccati  2
      idiv=0
      go to 2
   32 nn=1
      go to 35
   34 nn=0
   35 call yfclen(eta,rau,fo,fpo,go,gpo,sigma0,idiv,nn)
      return
   30 if(etaro-12.d0)32,32,23
   23 tra=eta2-6.75d0*(eta**0.4d0)
      if(rau-tra)6,6,24
   24 ind=1
      jnd=1
      ro=rau
      rau=tra
      rau0=tra
!             riccati  1
 
    6 x=rau/eta2
      u= (1.d0-x)/x
      x2=x*x
      ru= dsqrt(u)
      rx= dsqrt(x)
      tre= 1.d0/(u*ru*eta2)
      trb=tre*tre
      fi= (dsqrt((1.d0-x)*x)+dasin(rx)-1.570796326794897d0)*eta2
      tr1= -0.25d0*dlog(u)
  602 tr2= -((9.d0*u+6.d0)*u+5.d0)/48.d0
  603 tr3= ((((-3.d0*u-4.d0)*u+6.d0)*u+12.d0)*u+5.d0)/64.d0
  604 tr4=- ((((((u+2.d0)*945.d0*u+1395.d0)*u+12300.d0)*u+25191.d0) &
     &*u+19890.d0)*u+5525.d0)/46080.d0
  605 tr5= ((((((((-27.d0*u-72.d0)*u-68.d0)*u+360.d0)*u+2190.d0) &
     &*u+4808.d0)*u+5148.d0)*u+2712.d0)*u+565.d0)/2048.d0
  606 tr6=- (((((((((g61*u+g62)*u+g63)*u+g64)*u+g65)*u+g66)*u+g67) &
     &*u+g68)*u+g69)*u+g610)*u+g611
  607 tr7= ((((((((((((-81.*u-324.)*u-486.)*u-404.)*u+4509.)*u+52344.) &
     &*u+233436.)*u+567864.)*u+838521.)*u+775884.)*u+441450.) &
     &*u+141660.)*u+19675.) /6144.
  608 tr8= (((((((((((((g81*u+g82)*u+g83)*u+g84)*u+g85)*u+g86)*u+g87) &
     &*u+g88)*u+g89)*u+g810)*u+g811)*u+g812)*u+g813)*u+g814)*u+g815
      psip=psip+tra
      xxx=138.1551055796428d0
      fi= fi+tre*(tr2+trb*(tr4+trb*(tr6+trb*tr8)))
      psi=-fi
      indg=idint(psi/xxx)
      idiv=60*indg
      tra= tr1+trb*(tr3+trb*(tr5+trb*tr7))
      fi=fi+tra
      psi=psi+tra
 
      fip=ru*eta2
      tra=1.d0/(x2*u)
      tr1=0.25d0
      tre= tre/(x2*x2*u)
      trb=trb/(x2*x2)
  702 tr2= -(8.d0*x-3.d0)/32.d0
  703 tr3= ((24.d0*x-12.d0)*x+3.d0)/64.d0
  704 tr4= (((-1536.d0*x+704.d0)*x-336.d0)*x+63.d0)/2048.d0
  705 tr5= ((((1920.d0*x-576.d0)*x+504.d0)*x-180.d0)*x+27.d0)/1024.d0
  706 tr6 = ((((-gp66*x+gp65)*x-gp64)*x+gp63)*x-gp62)*x+gp61
  707 tr7= - ((((((-40320.d0*x-10560.d0)*x-13248.d0)*x+7560.d0) &
     &*x-3132.d0)*x+756.d0)*x-81.d0) / 2048.d0
  708 tr8 =- ((((((gp88*x+gp87)*x+gp86)*x-gp85)*x+gp84)*x-gp83)*x+gp82) &
     &*x-gp81
      fip=fip +tre*(tr2+trb*(tr4+trb*(tr6+trb*tr8)))
      tra= tra*(tr1+trb*(tr3+trb*(tr5+trb*tr7)))
      fip=fip+tra
      psip=-fip
      if(indg.eq.0)go to 8
      psi=psi-xxx*dfloat(indg)
      fi =fi +xxx*dfloat(indg)
    8 fo=0.5d0*dexp(fi)
      go= dexp(psi)
      fpo= fo*fip/eta2
      gpo=go*psip/eta2
      if(jnd.eq.0)return
      rau=ro
      go=fo
      gpo=fpo
   27 x=rau0-ro
      x2=x*x
      x3=x*x2
      unr=1.d0/rau0
      etr0=1.d0-2.d0*eta*unr
      u0=go
      u1=-x*gpo
      u2=-0.5d0*etr0*x2*u0
      s=u0+u1+u2
      v1=u1/x
      v2=2.d0*u2/x
      t=v1+v2
      xn=3.d0
 
      do 10 n=3,10000
!     n=n
      xn1=xn-1.d0
      xn1=xn*xn1
      u3=x*u2*unr*(1.d0-2.d0/xn)-etr0*u1*x2/xn1+x3*u0*unr/xn1
      s=s+u3
      v3=xn*u3/x
      t=t+v3
   16 if(dabs(u3/s).gt.1.d-10)go to 11
      if(dabs(v3/t).le.1.d-10)go to 17
   11 u0=u1
      u1=u2
      u2=u3
      xn=xn+1.d0
   10 continue
   17 if(ig)25,26,25
   25 go=s
      gpo=-t
      fo=ho
      fpo=hpo
      return
   26 ho=s
      hpo=-t
   18 et0=eta**(0.166666666666667)
      etad=etac*etac
      et=eta**(0.6666666666666667)
      et1=et*et
      et2=et1*et1
      et3=et2*et
      et4=etad*et
      et5=et4*et
      fo=1.d0-q(1)/et1-q(2)/etac-q(3)/et3-q(4)/etad-q(5)/et5
      go=1.d0+q(1)/et1-q(2)/etac+q(3)/et3-q(4)/etad+q(5)/et5
      fpo=1.d0+qp(1)/et+qp(2)/etac+qp(3)/et2+qp(4)/etad+qp(5)/et4
      gpo=1.d0-qp(1)/et+qp(2)/etac-qp(3)/et2+qp(4)/etad-qp(5)/et4
      fo=0.7063326373d0*et0*fo
      go=1.223404016d0*et0*go
      fpo=0.4086957323d0*fpo/et0
      gpo=-0.7078817734d0*gpo/et0
      idiv=0
      if(ind.eq.0)return
      ig=1
      rau0=eta2
      go to 27
    2 x=eta2/rau
      x2=x*x
      u=1.d0-x
      ru= dsqrt(u)
      u3=u*u*u
      trd= 1.d0/(u3*eta2*eta2)
      trc=x2*trd
      tre=1.d0/(u*ru*eta2)
      fi= -0.25d0*dlog(u)
      trb=trd/64.d0
      tr3= (((3.d0*u-4.d0)*u-6.d0)*u+12.d0)*u-5.d0
  501 tr5= ((((((((-27.d0*u+72.d0)*u-68.d0)*u-360.d0)*u+2190.d0) &
     &*u-4808.d0)*u+5148.d0)*u-2712.d0)*u+565.d0)/32.d0
  502 tr7= ((((((((((((81.d0*u-324.d0)*u+486.d0)*u-404.d0)*u-4509.d0)  &
     &*u+52344.d0)*u-233436.d0)*u+567864.d0)*u-838521.d0)*u+775884.d0) &
     &*u-441450.d0)*u+141660.d0)*u-19675.d0)/96.d0
      fi= fi+trb*(tr3+trd*(tr5+trd*tr7))
 
      fip=0.25d0/u
      trb=3.d0*trc/(64.d0*u)
      tr3= (x-4.d0)*x+8.d0
  511 tr5= ((((9.d0*x-60.d0)*x+168.d0)*x-192.d0)*x+640.d0)/16.d0
  512 tr7= ((((((-27.d0*x+252.d0)*x-1044.d0)*x+2520.d0)*x-4416.d0) &
     &*x-3520.d0)*x-13440.d0)/32.d0
      fip =fip+trb*(tr3+trc*(tr5+trc*tr7))
      tra= dabs((ru-1.d0)/(ru+1.d0))
      psi= (0.5d0*dlog(tra)+ru/x)*eta2+0.785398163397448d0
      tr2= -((9.d0*u-6.d0)*u+5.d0)/48.d0
  521 tr4= ((((((u-2.d0)*945.d0*u+1395.d0)*u-12300.d0)*u+25191.d0) &
     &*u-19890.d0)*u+5525.d0)/46080.d0
  522 tr6 = (((((((((-g61*u+g62)*u-g63)*u+g64)*u-g65)*u+g66)*u-g67) &
     &*u+g68)*u-g69)*u+g610)*u-g611
  523 tr8= (((((((((((((g81*u-g82)*u+g83)*u-g84)*u+g85)*u-g86)*u+g87) &
     &*u-g88)*u+g89)*u-g810)*u+g811)*u-g812)*u+g813)*u-g814)*u+g815
      psi= psi+tre*(tr2+trd*(tr4+trd*(tr6+trd*tr8)))
      psip = -ru*eta2/x2
      trb=tre*x/u
      tr2= (3.d0*x-8.d0)/32.d0
  531 tr4= - (((63.d0*x-336.d0)*x+704.d0)*x-1536.d0)/2048.d0
  532 tr6 = ((((gp61*x-gp62)*x+gp63)*x-gp64)*x+gp65)*x-gp66
  533 tr8 = ((((((-gp81*x+gp82)*x-gp83)*x+gp84)*x-gp85)*x+gp86)*x+gp87) &
     &*x+gp88
      psip =psip+ trb*(tr2+trc*(tr4+trc*(tr6+trc*tr8)))
      tra= dexp(fi)
      fo= tra*dsin(psi)
      go= tra*dcos(psi)
      if(indice)535,536,535
  535 tra=fo
      fo=go
      go=-tra
  536 tra=-eta2/rauc
      fpo=(fip*fo+psip*go)*tra
      gpo=(fip*go-psip*fo)*tra
      return
      end subroutine
!====================================================================!
      subroutine dfcz0(eta,ro,f0,fp0,g0,gp0,sigma0,iexp)
      implicit real*8 (a-h,o-z)
      dimension a1(110),a2(110),b1(110),b2(110)
      if(ro)2,2,1
    2 write (6,1000)
 1000 format (21h ro negatif ou nul **)
      return
    1 if(eta-30.d0)3,3,4
    4 if(dabs(eta)-5.d2)28,28,29
   28 call yfrica(eta,ro,f0,fp0,g0,gp0,sigma0,iexp)
      return
   29 write (6,1001)
 1001 format (42h valeur absolue de eta supe-&eu-e a 500 **)
      return
    3 if(eta+8.d0)4,5,5
    5 if(eta)6,7,6
    7 f0=dsin(ro)
      g0=dcos(ro)
      fp0=g0
      gp0=-f0
      iexp=0
      sigma0=0.d0
      return
    6 borne=1.666666666666667d0*dabs(eta)+7.5d0
      if(ro-borne)8,9,9
    9 call yfasym(eta,ro,f0,fp0,g0,gp0,sigma0,iexp)
      return
    8 if(eta-10.d0)10,11,11
   10 if(eta)12,12,13
   13 if(ro-2.d0)14,12,12
   11 if(eta-(5.d0*ro+6.d1)/7.d0)12,12,14
   12 call yfasym(eta,borne,f0,fp0,g0,gp0,sigma0,iexp)
      h=borne
      dh=f0/h
      if(eta)20,21,21
   20 n=-0.5d0*eta+5.d0
      go to 22
   21 n=eta/5.d0+5.d0
   22 n=5*(n+1)
      z=4.d0/h
      y=1.d0-(eta+eta)*z
      a1(n+2)=1.d-55
      a1(n+3)=0.d0
      a1(n+4)=1.d-64
      b1(n+3)=1.d-50
      b1(n+4)=1.d-68
      a2(n+2)=0.d0
      a2(n+3)=1.d-74
      a2(n+4)=1.d-53
      b2(n+3)=0.d0
      b2(n+4)=1.d-66
      m=n+2
      di=n
      do 23 ii=2,m
      i=m-ii+2
      b1(i)=b1(i+2)+z*(di+1.d0)*a1(i+1)
      s=a1(i+2)+y*(a1(i+1)-a1(i))
      q=(di+2.d0)*b1(i)+(di-1.d0)*b1(i+1)
      a1(i-1)=s-z*q
      b2(i)=b2(i+2)+z*(di+1.d0)*a2(i+1)
      s=a2(i+2)+y*(a2(i+1)-a2(i))
      q=(di+2.d0)*b2(i)+(di-1.d0)*b2(i+1)
      a2(i-1)=s-z*q
      if(i.ge.n)go to 23
      d=-(b2(i+2)+b2(i))/(b1(i+2)+b1(i))
      do 24 j=i,m
      a2(j)=a2(j)+d*a1(j)
      b2(j)=b2(j)+d*b1(j)
   24 continue
      a2(i-1)=a2(i-1)+d*a1(i-1)
   23 di=di-1.d0
      q=a1(3)-a1(1)
      c=a2(3)-a2(1)
      c=q/c
      x1=0.5d0*(a1(2)-c*a2(2))
      do 25 i=3,m
      x1=x1+a1(i)-c*a2(i)
   25 continue
      x1=dh/x1
      x2=-c*x1
!$omp parallel do
      do 26 i=2,m
      b1(i)=x1*b1(i)+x2*b2(i)
      a1(i)=x1*a1(i)+x2*a2(i)
   26 continue
      a1(1)=x1*a1(1)+x2*a2(1)
      b1(1)=0.d0
      x=ro/h
      y=2.d0*x-1.d0
      t1=1.d0
      t2=y
      result=0.5d0*a1(2)+y*a1(3)
      derive=0.5d0*b1(2)+y*b1(3)
      do 27 i=2,n
      ti=2.d0*y*t2-t1
      t1=t2
      t2=ti
      result=result+ti*a1(i+2)
      derive=derive+ti*b1(i+2)
   27 continue
      f0=result*ro
      fp0=derive*ro+result
      go to 30
!   serie origine reguliere
   14 pi=3.141592653589793d0
      call jflgam(1.d0,eta,tra,sigma0,ntruc)
      iexp=0
      ro2=ro*ro
      etap=eta+eta
      pieta=pi*eta
      pieta2=0.5d0*pieta
      b=dexp(pieta2)*dsqrt(dsinh(pieta)/pieta)
      u0=0.d0
      u1=ro
      u=u0+u1
      up=1.d0
      xn=2.d0
      do 15 n=2,10000
      xn1=xn*(xn-1.d0)
      u2=(etap*ro*u1-ro2*u0)/xn1
      u=u+u2
      up=up+xn*u2/ro
   17 if(dabs(u2/u).lt.1.d-10)go to 19
   18 u0=u1
      u1=u2
      xn=xn+1.d0
   15 continue
   19 f0=u/b
      fp0=up/b
   30 call yfireg(eta,ro,g0,gp0)
      return
      end subroutine

!====================================================================!
      subroutine jfdelg (xd,yd,par,pai,nbchif)
      real*8 xd,yd,par,pai,test,c,pi
      real*8 x,y,u,v,tra,tra1,cosi,sini
      real*8 cos2i,sin2i,zmod,depi
      real*8 trb,xx
      real*8 rac2,pis4
      real*8 supint
      dimension test(7),c(6)
      data test/2.9152d+7,2.2958d+3,1.4124d+2,3.9522d+1,19.6611d0, &
     &12.791d0,-10.d0/
      data rac2/0.3465735902799726/,pis4/0.785398163397448/
      data c/8.333333333333333d-2,-8.33333333333333d-3, &
     &3.968253968253968d-3,-4.166666666666667d-3,7.575757575757576d-3, &
     &-2.109279609279609d-2/  
      data pi/3.141592653589793/
      data depi/6.283185307179586/
      data supint/2147483647.d0/
      x=dabs(xd)
      xx=x
      nbchif=15
      if(yd)1,2,1
    1 y=dabs(yd)
      kr=1
      i=dmod(10.99d0-x,supint)
!     translation
      if(i)3,3,4
    4 tra=i
      x=x+tra
!     logarithme(x+iy) (x,y positifs)
    3 if(x-y)5,6,7
    5 tra1=x/y
      trb=1.d0+tra1*tra1
      tra=y*dsqrt(trb)
      sini=1./(trb*y)
      cosi=sini*tra1
      tra1=y/x
      go to 11
    6 u=rac2+dlog(x)
      v=pis4
      sini=0.5d0/x
      cosi=sini
      go to 10
    7 tra1=y/x
      trb=1.d0+tra1*tra1
      tra=x*dsqrt(trb)
      cosi=1./(trb*x)
      sini=cosi*tra1
   11 u=dlog(tra)
      v=datan(tra1)
!     developpement asymptotique ( x superieur a 10 )
   10 par=u-0.5*cosi
      pai=v+0.5*sini
      zmod=x+y
      if(zmod-test(1))13,13,14
   13 sin2i=(sini*cosi)+(sini*cosi)
      cos2i=(cosi+sini)*(cosi-sini)
      sini=sin2i
      cosi=cos2i
      k=1
      go to 15
   16 tra=cosi*cos2i-sini*sin2i
      sini=sini*cos2i+cosi*sin2i
      cosi=tra
   15 par=par-c(k)*cosi
      pai=pai+c(k)*sini
      k=k+1
      if(zmod-test(k))16,16,14
!     translation inverse
   17 i=i-1
      x=i
      x=xx+x
   56 if(x-y)55,55,57
   55 tra1=x/y
      trb=x*tra1+y
      sini=1.d0/trb
      cosi=tra1/trb
      go to 19
   57 tra1=y/x
      trb=x+y*tra1
      cosi=1.d0/trb
      sini=tra1/trb
   19 par=par-cosi
      pai=pai+sini
   14 if(i)18,18,17
 
!     controle du quadrant
   18 if(xd)20,61,21
   61 tra=pi*y
      if(tra-1.d-2)300,300,301
  300 trb= tra*(2.d0+tra*(-2.d0+tra*(1.333333333333333d0+tra*( &
     &-0.6666666666666666d0+tra*(0.2666666666666666d0+tra*( &
     &-0.08888888888888888d0+tra*0.02539682539682540d0))))))
      trb=(2.d0-trb)/trb
      go to 302
  301 trb= dexp(-tra-tra)
      trb=(1.d0+trb)/(1.d0-trb)
  302 pai=0.5d0*(1.d0/y+pi*trb)
   21 if(yd)28,100,100
!     x+iy change en -x-iy
   20 tra=dexp(-depi*y)
      trb=tra*tra
      cos2i=depi*dmod(x,1.d0)
      sin2i=-2.d0*tra*dcos(cos2i)+1.d0+trb
      par=par+cosi+depi*tra*dsin(cos2i)/sin2i
      pai=pai-sini+pi*(trb-1.d0)/sin2i
      if(yd)100,100,28
   28 pai=-pai
   35 write (6,1001)
 1001 format (30h argument de jfdelg trop grand)
      go to 50
   34 y=dmod(x,supint)
      if(y) 400,36,400
  400 if(y-0.99d0) 33,33,405
  405 tra= idint(y+0.1d0)
      if(dabs(y-tra)-5.d-15)36,36,33
   31 if(x-4503599627370496.d0)34,35,35
 
!     argument dans -pi,pi
  100 tra=dabs(pai/depi)
      if(tra-1.d+15)203,204,204
  204 nbchif=0
      pai=0.d0
      go to 29
  203 if(tra-1.d0)205,205,206
  206 nbchif=dlog10(tra)
      nbchif=14-nbchif
      tra=dmod(tra,supint)
      pai=dmod(tra,1.d0)*dsign(depi,pai)
  205 if(dabs(pai)-pi)29,29,207
  207 pai=pai-dsign(depi,pai)
      go to 29
!        delgamma reel
    2 pai=0.d0
      if(xd)31,32,33
!     conditions d existence
   32 write (6,1000)
 1000 format (21h jfdelg(0) est infini)
      go to 50
   36 write (6,1002)
 1002 format (28h jfdelg (-entier) est infini)
   50 par=1.d+74
      nbchif=0
      go to 29
!     translation
   33 i=dmod(10.99d0-x,supint)
      if(i)37,37,38
   38 tra=i
      x=x+tra
!     developpement asymptotique
   37 y=dlog(x)
      par=y-0.5d0/x
      if(x-test(1))39,39,43
   39 cos2i=1.d0/(x*x)
      cosi=cos2i
      k=1
      go to 41
   42 cosi=cosi*cos2i
   41 par=par-c(k)*cosi
      k=k+1
      if(x-test(k))42,42,40
!     translation inverse
   44 i=i-1
      x=i
      x=xx+x
      par=par-1.d0/x
   40 if(i)43,43,44
!     x negatif
   43 if(xd)45,29,29
   45 par=par+1.d0/x
      y=pi*dmod(x,2.d0)
      par=par+pi*dcos(y)/dsin(y)
      entry jfdev1
   29 return
      end subroutine
!**********************************************************************!
      end module
