
      subroutine wzsub(x,y,u,v)
!  *********************************************************************
!
!  This subroutine sets u=real(w(z)) and v=imag(w(z)), where z=x+i*y and
!  where w(z) is the complex error function defined by formula 7.1.3 in
!  "Handbook of Mathematical functions [eds. M.Abramowitz & I.A.Stegun,
!  Washington, 1966].  The absolute error of the computed value is less
!  than 1E-8.
!
!  *** Note.  Subroutine WZSET must have been called before this sub-
!  routine can be used.
!
!  For (x,y) inside the rectangle with opposite corners (xcut,0) and
!  (0,ycut), where xcut and ycut have been set by WZSET, an interpo-
!  lation formula is used.  For (x,y) outside this rectangle, a two-
!  term rational approximation is used.
!
!  (G.A.Erskine, 29.09.1997)
!
!
!  Third-order divided-difference interpolation over the corners of a
!  square [e.g. formula (2.5.1) in "Introduction to Numerical Analysis"
!  (F.B.Hildebrand New York, 1957), but with complex nodes and
!  function values].
!
!  In the interpolation formula the corners of the grid square contain-
!  ing (x,y) are numbered (0,0)=3, (h,0)=4, (h,h)=1, (0,h)=2.
!  Identifiers d, dd and ddd denote divided-differences of orders 1, 2
!  and 3 respectively, and a preceding 't' indicates twice the value.
!
!  *********************************************************************
      implicit none
      integer k,mu,nu
      double precision a1,a2,b1,b2,d12i,d12r,d23i,d23r,d34i,d34r,p,     &
     &q,qsq,r,simag,sreal,t,tdd13i,tdd13r,tdd24i,tdd24r,tdddi,tdddr,ti, &
     &tr,u,usum,usum3,v,vsum,vsum3,w1i,w1r,w2i,w2r,w3i,w3r,w4i,w4r,x,xh,&
     &xhrel,y,yh,yhrel
      integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblo,nblz,ncom,ncor1,    &
     &nelb,nele,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran,    &
     &nrco,ntr,nzfz
      parameter(npart = 64,nmac = 1)
      parameter(nele=1200,nblo=600,nper=16,nelb=140,nblz=20000,         &
     &nzfz = 300000,mmul = 20)
      parameter(nran = 2000000,ncom = 100,mran = 500,mpa = 6,nrco = 5,  &
     &nema = 15)
      parameter(mcor = 10,mcop = mcor+6, mbea = 99)
      parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
      parameter(nmon1 = 600,ncor1 = 600)
      parameter(ntr = 20,nbb = 350)
      integer idim,kstep,nx,ny
      double precision h,half,hrecip,one,wtimag,wtreal,xcut,ycut
      parameter ( xcut = 7.77d0, ycut = 7.46d0 )
      parameter ( h = 1.d0/63.d0 )
      parameter ( nx = 490, ny = 470 )
      !parameter ( idim = (nx+2)*(ny+2) )
      parameter ( idim = 232224 )
      parameter ( half = 0.5d0, one = 1.d0 )
      common /wzcom1/ hrecip, kstep
      common /wzcom2/ wtreal(idim), wtimag(idim)
      parameter ( a1 = 0.5124242248d0, a2 = 0.0517653588d0 )
      parameter ( b1 = 0.2752551286d0, b2 = 2.7247448714d0 )

!-----------------------------------------------------------------------
      if ( x.ge.xcut .or. y.ge.ycut ) goto 1000
      xh = hrecip*x
      yh = hrecip*y
      mu = int(xh)
      nu = int(yh)
!  Compute divided differences.
      k = 2 + mu + nu*kstep
      w4r = wtreal(k)
      w4i = wtimag(k)
      k = k - 1
      w3r = wtreal(k)
      w3i = wtimag(k)
      d34r = w4r - w3r
      d34i = w4i - w3i
      k = k + kstep
      w2r = wtreal(k)
      w2i = wtimag(k)
      d23r = w2i - w3i
      d23i = w3r - w2r
      tr = d23r - d34r
      ti = d23i - d34i
      tdd24r = ti - tr
!hr05 tdd24i = - ( tr + ti )
      tdd24i = -1d0* ( tr + ti )                                         !hr05
      k = k + 1
      w1r = wtreal(k)
      w1i = wtimag(k)
      d12r = w1r - w2r
      d12i = w1i - w2i
      tr = d12r - d23r
      ti = d12i - d23i
      tdd13r = tr + ti
      tdd13i = ti - tr
      tdddr = tdd13i - tdd24i
      tdddi = tdd24r - tdd13r
!  Evaluate polynomial.
      xhrel = xh - dble(mu)
      yhrel = yh - dble(nu)
      usum3 = half*( tdd13r + ( xhrel*tdddr - yhrel*tdddi ) )
      vsum3 = half*( tdd13i + ( xhrel*tdddi + yhrel*tdddr ) )
      yhrel = yhrel - one
      usum = d12r + ( xhrel*usum3 - yhrel*vsum3 )
      vsum = d12i + ( xhrel*vsum3 + yhrel*usum3 )
      xhrel = xhrel - one
      u = w1r + ( xhrel*usum - yhrel*vsum )
      v = w1i + ( xhrel*vsum + yhrel*usum )
      return
!
!  Two-term rational approximation to w(z) [Footnote to Table 7.9
!  in "Handbook of Mathematical Functions (eds. M.Abramowitz &
!  I.A.Stegun, Washington, 1966), but with additional digits in
!  the constants]:
!              u+i*v = i*z*( a1/(z**2-b1) + a2/(z**2-b2) ).
!  Maximum absolute error:
!        <1.E-6  for  x>=4.9  or  y>=4.4
!        <1.E-7  for  x>=6.1  or  y>=5.7
!        <1.E-8  for  x>=7.8  or  y>=7.5
!
 1000 p=x**2-y**2
!hr05 q=2.d0*x*y
      q=(2.d0*x)*y                                                       !hr05
      qsq=q**2
!  First term.
      t=p-b1
      r=a1/(t**2+qsq)
      sreal=r*t
!hr05 simag=-r*q
      simag=(-1d0*r)*q                                                   !hr05
!  Second term
      t=p-b2
      r=a2/(t**2+qsq)
      sreal=sreal+r*t
      simag=simag-r*q
!  Multiply by i*z.
!hr05 u=-(y*sreal+x*simag)
      u=-1d0*(y*sreal+x*simag)                                           !hr05
      v=x*sreal-y*simag
      return
!
      end
