c=======================================================================
c FFTLog
c    computes the discrete Fast Fourier Transform
c    or Fast Hankel Transform (of arbitrary real index)
c    of a periodic logarithmic sequence.
c
c---Start of revision history-------------------------------------------
c
c Original version 27 March 1999.
c
c 5 July 1999:
c Bug fix in singular transform.
c Store amp, arg of last element of wsave separately, for even n,
c enlarging wsave by 1, for even n and q != 0.

c 17 July 1999:
c Reorganized wsave more logically.
c drffti now gets first, not last, 2*n+15 elements of wsave.
c fhti now treats last element of wsave with others, not by itself;
c wsave is enlarged from original, for even n,
c by 1 for q = 0, and by 2 for q != 0.
c
c 18 July 1999:
c A backward transform (dir = -1) is the same as
c a forward transform with q -> -q [and rk -> 1/rk]
c for any kr if n is odd,
c but only for low-ringing kr if n is even.
c Original falsely stated that transforms were same irrespective of kr.
c
c 19 Dec 1999:
c Slightly spiffed up interactive option to change kr to low-ringing kr.
c
c 13 Mar 2000:
c Made g77-safe:
c 1. Checked for any automatic variables that needed explicitly
c    `save'ing - actually there weren't any.
c    g77 is unusual in that it does NOT save variables not in
c    common or data statements.
c 2. All double precision constants explicitly suffixed d0.
c    g77 is unusual in that it does NOT automatically promote
c    real constants such as 1.1 to double precision 1.1d0.
c 
c---End of revision history---------------------------------------------
c
c For more information about FFTLog, see
c
c http://casa.colorado.edu/~ajsh/FFTLog/
c
c Andrew J S Hamilton March 1999
c email: Andrew.Hamilton@Colorado.EDU
c
c Refs:	Talman J. D., 1978, J. Comp. Phys., 29, 35
c	Hamilton A. J. S., 2000, MNRAS, 312, 257
c	( http://xxx.lanl.gov/abs/astro-ph/9905191 )
c
c FFTLog uses the NCAR suite of FFT routines,
c and a modified version of the complex Gamma function
c from the gamerf package at
c http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html .
c The original gamerf copyright statement states:
c   Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
c   You may use, copy, modify this code for any purpose and
c   without fee. You may distribute this ORIGINAL package.
c
c Permission to distribute the modified gamma function code
c with the FFTLog package has been granted
c (email from Takuya Ooura to Andrew Hamilton dated 16 March 1999).
c
c-----------------------------------------------------------------------
c Note: email programs may zap the 8-bit characters in the program
c       comments by converting them into 7-bit characters.
c       If this happens, get an unzapped copy from
c       http://casa.colorado.edu/~ajsh/FFTLog/ .
c       The fortran code itself contains only 7-bit characters;
c       the 8-bit characters occur only in comments.
c
c       If you are on a UNIX machine and 8-bit characters are not
c       showing up properly (they may appear as \nnn where n are
c       integers), then try putting
c          setenv LC_CTYPE iso_8859_1
c       in your .login file.
c
c-----------------------------------------------------------------------
c FFTLog computes a discrete version of
c the Hankel Transform (= Fourier-Bessel Transform)
c with a power law bias (k r)^q
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J  (k r) k dr
c           /              mu
c          0
c
c          infinity
c           /           -q
c    a(r) = | ã(k) (k r)   J  (k r) r dk
c           /               mu
c          0
c
c where J_mu is the Bessel function of order mu.
c The index mu may be any real number, positive or negative.
c
c The input array a_j is a periodic sequence of length n,
c uniformly logarithmically spaced with spacing dlnr
c    a_j = a(r_j)   at   r_j = r_c exp[(j-j_c) dlnr]
c centred about the point r_c.  The central index j_c = (n+1)/2
c is 1/2 integral if n is even.
c Similarly, the output array ã_j is a periodic sequence of length n,
c also uniformly logarithmically spaced with spacing dlnr
c    ã_j = ã(k_j)   at   k_j = k_c exp[(j-j_c) dlnr]
c centred about the point k_c.
c
c The centre points r_c and k_c of the periodic intervals may be
c chosen arbitrarily; but it would be normal to choose the product
c kr = k_c r_c = k_j r_(n+1-j) = k_(n+1-j) r_j
c to be about 1 (or 2, or pi, to taste).
c
c The FFTLog algorithm is (see Hamilton 2000):
c 1. FFT the input array a_j to obtain the Fourier coefficients c_m ;
c 2. Multiply c_m by 
c       u_m = (kr)^[- i 2 m pi/(n dlnr)] U_mu[q + i 2 m pi/(n dlnr)]
c    where
c       U_mu(x) = 2^x Gamma[(mu+1+x)/2] / Gamma[(mu+1-x)/2]
c    to obtain c_m u_m ;
c 3. FFT c_m u_m back to obtain the discrete Hankel transform ã_j .
c
c-----------------------------------------------------------------------
c The Fourier sine and cosine transforms
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) sin(k r) dr
c                      /
c                     0
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) cos(k r) dr
c                      /
c                     0
c
c may be regarded as special cases of the Hankel transform
c with mu = 1/2 and -1/2 since
c
c    sqrt(2/pi) sin(x) = sqrt(x) J   (x)
c                                 1/2
c
c    sqrt(2/pi) cos(x) = sqrt(x) J    (x)
c                                 -1/2
c
c The Fourier transforms may be done by making the substitutions
c                 q-(1/2)                      -q-(1/2)
c    A(r) = a(r) r          and   Ã(k) = ã(k) k    
c
c and Hankel transforming a(r) with a power law bias (k r)^q
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J    (k r) k dr
c           /              ±1/2
c          0
c
c Different choices of power law bias q lead to different discrete
c Fourier transforms of A(r), because the assumption of periodicity of
c a(r) = A(r) r^[-q+(1/2)] is different for different q.
c
c If A(r) is a power law, A(r) proportional to r^[q-(1/2)],
c then applying a bias q yields a discrete Fourier transform Ã(k)
c that is exactly equal to the continuous Fourier transform,
c because then a(r) is a constant, which is a periodic function.
c
c-----------------------------------------------------------------------
c The Hankel transform
c
c          infinity
c           /
c    Ã(k) = | A(r) J  (k r) k dr
c           /       mu
c          0
c
c may be done by making the substitutions
c                 q                      -q
c    A(r) = a(r) r    and   Ã(k) = ã(k) k    
c
c and Hankel transforming a(r) with a power law bias (k r)^q
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J  (k r) k dr
c           /              mu
c          0
c
c Different choices of power law bias q lead to different discrete
c Hankel transforms of A(r), because the assumption of periodicity of
c a(r) = A(r) r^-q is different for different q.
c
c If A(r) is a power law, A(r) proportional to r^q,
c then applying a bias q yields a discrete Hankel transform Ã(k)
c that is exactly equal to the continuous Hankel transform,
c because then a(r) is a constant, which is a periodic function.
c
c-----------------------------------------------------------------------
c ------------------------
c There are five routines:
c ------------------------
c Comments in the subroutines contain further details.
c
c (1) subroutine fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
c     is an initialization routine.
c
c (2) subroutine fftl(n,a,rk,dir,wsave)
c     computes the discrete Fourier sine or cosine transform
c     of a logarithmically spaced periodic sequence.
c     This is a driver routine that calls fhtq.
c
c (3) subroutine fht(n,a,dir,wsave)
c     computes the discrete Hankel transform
c     of a logarithmically spaced periodic sequence.
c     This is a driver routine that calls fhtq.
c
c (4) subroutine fhtq(n,a,dir,wsave)
c     computes the biased discrete Hankel transform
c     of a logarithmically spaced periodic sequence.
c     This is the basic FFTLog routine.
c
c (5) real*8 function krgood(mu,q,dlnr,kr)
c     takes an input kr and returns the nearest low-ringing kr.
c     This is an optional routine called by fhti.
c
c=======================================================================
c THE FFTLog CODE
c=======================================================================
      subroutine fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
      integer n,kropt
      logical ok
      real*8 mu,q,dlnr,kr,wsave(*)

c
c fhti initializes the working array wsave
c used by fftl, fht, and fhtq, and by NCAR routines drfftf and drfftb.
c fhti need be called once, whereafter fftl, fht, or fhtq may be
c called many times, as long as n, mu, q, dlnr, and kr remain unchanged.
c fhti should be called each time n, mu, q, dlnr, or kr is changed.
c The work array wsave should not be changed between calls to fftl,
c fht, or fhtq.
c
c If you are using the g77 fortran compiler, which by default does not
c save locally declared variables (unlike most other fortran compilers),
c then you should explicitly
c       save wsave
c in the calling program.
c
c  Input: n = number of points in the array to be transformed;
c             n may be any positive integer, but the NCAR FFT routines
c             run fastest if n is a product of small primes 2, 3, 5.
c	  mu = index of J_mu in Hankel transform;
c	       mu may be any real number, positive or negative.
c	  q = exponent of power law bias;
c	      q may be any real number, positive or negative.
c             If in doubt, use q = 0, for which case the Hankel
c             transform is orthogonal, i.e. self-inverse,
c             provided also that, for n even, kr is low-ringing.
c             Non-zero q may yield better approximations to the
c             continuous Hankel transform for some functions.
c         dlnr = separation between natural log of points;
c		 dlnr may be positive or negative.
c	  kr = k_c r_c where c is central point of array
c              = k_j r_(n+1-j) = k_(n+1-j) r_j .
c	       Normally one would choose kr to be about 1
c	       (or 2, or pi, to taste).
c         kropt = 0 to use input kr as is;
c                 1 to change kr to nearest low-ringing kr, quietly;
c                 2 to change kr to nearest low-ringing kr, verbosely;
c                 3 for option to change kr interactively.
c Output: wsave = working array used by fftl, fht, and fhtq,
c                 and by NCAR FFT routines drfftf and drfftb.
c                 wsave should be dimensioned at least:
c                 for q = 0 (unbiased transform):
c                    2*n+2*(n/2)+18
c                       [= 3*n+18 if n even, 3*n+17 if n odd]
c                 for q != 0 (biased transform):
c                    2*n+3*(n/2)+19
c                       [= (7*n)/2+19 if n even, (7*n)/2+18 if n odd]
c                 The first 2*n+15 elements of wsave are used by
c                 the NCAR FFT routines.
c         ok = .true. if all went ok;
c              .false. if not; currently only occurs if interactive
c                      response causes exit.
c
c        parameters
      real*8 PI
      parameter (PI=3.141592653589793238462643383279502884197d0)
      real*8 ZERO,ONE,TWO
      parameter (ZERO=0.d0,ONE=1.d0,TWO=2.d0)
      real*8 ROUND
      parameter (ROUND=1.d-15)
c        externals
      real*8 krgood
      complex*16 cdgamma
c        local (automatic) variables
      character*1 go
      character*64 temp
      integer l,m
      real*8 amp,arg,d,ln2,ln2kr,xm,xp,y
      complex*16 zm,zp
c
c--------adjust kr
c        keep kr as is
      if (kropt.eq.0) then
        continue
c        change kr to low-ringing kr quietly
      elseif (kropt.eq.1) then
        kr=krgood(mu,q,dlnr,kr)
c        change kr to low-ringing kr verbosely
      elseif (kropt.eq.2) then
        d=krgood(mu,q,dlnr,kr)
        if (abs(kr/d-ONE).gt.ROUND) then
          kr=d
c          write (*,'(" kr changed to",g24.16)') kr
        endif
c        option to change kr to low-ringing kr interactively
      else
        d=krgood(mu,q,dlnr,kr)
        if (abs(kr/d-ONE).gt.ROUND) then
c        fortran demonstrates its inferiority to C
c          write (*,'(" change kr = ",$)')
c          write (temp,*) kr
c          write (*,'(" to low-ringing kr = ",$)')
c          write (temp,*) d
c          write (*,'("? [CR,y=yes, n=no, x=exit]: ",$)')
c          read (*,'(a1)') go
           go='y'
          if (go.eq.' '.or.go.eq.'y'.or.go.eq.'Y') then
            kr=d
c            write (*,'(" kr changed to",g24.16)') kr
          elseif (go.eq.'n'.or.go.eq.'N') then
c            print *,'kr left unchanged at',kr
          else
c            print *,'exit'
            goto 300
          endif
        endif
      endif
c--------return if n is <= 0
      if (n.le.0) goto 200
c--------initialize normal FFT
      call drffti(n,wsave)
c        drfft uses first 2*n+15 elements of wsave
      l=2*n+15
c--------put q, dlnr, kr in next 3 elements of wsave
      wsave(l+1)=q
      wsave(l+2)=dlnr
      wsave(l+3)=kr
c        so far used 2*n+18 elements of wsave
      l=l+3
c--------rest of wsave used by fhtq - unbiased case (q = 0)
      if (q.eq.ZERO) then
        ln2kr=log(TWO/kr)
        xp=(mu+ONE)/TWO
        d=PI/(n*dlnr)
        do m=1,n/2
c        y = m pi/(n dlnr)
          y=m*d
          zp=dcmplx(xp,y)
          zp=cdgamma(zp,1)
c        Argument of kr^(-2 i y) U_mu(2 i y)
          arg=TWO*(ln2kr*y+dimag(zp))
          wsave(l+2*m-1)=cos(arg)
          wsave(l+2*m)=sin(arg)
        enddo
c Altogether 2*n+18 + 2*(n/2) = 2*((3*n)/2)+18 elements used for q = 0,
c which is 3*n+18 for even n, 3*n+17 for odd n.
c--------rest of wsave used by fhtq - biased case (q != 0)
      else
        ln2=log(TWO)
        ln2kr=log(TWO/kr)
        xp=(mu+ONE+q)/TWO
        xm=(mu+ONE-q)/TWO
c........first element of rest of wsave
        y=ZERO
c        case where xp or xm is a negative integer
        if ((anint(xp).eq.xp.and.anint(xp).le.ZERO)
     *    .or.(anint(xm).eq.xm.and.anint(xm).le.ZERO)) then
c        case where xp and xm are both negative integers
c U_mu(q) = 2^q Gamma[xp]/Gamma[xm] is finite in this case
          if ((anint(xp).eq.xp.and.anint(xp).le.ZERO)
     *      .and.(anint(xm).eq.xm.and.anint(xm).le.ZERO)) then
c        Amplitude and Argument of U_mu(q)
            amp=exp(ln2*q)
            if (xp.gt.xm) then
              do m=1,nint(xp-xm)
                amp=amp*(xm+m-ONE)
              enddo
            elseif (xp.lt.xm) then
              do m=1,nint(xm-xp)
                amp=amp/(xp+m-ONE)
              enddo
            endif
            arg=anint(xp+xm)*PI
c        one of xp or xm is a negative integer
          else
c Transformation is singular if xp is -ve integer,
c and inverse transformation is singular if xm is -ve integer,
c but transformation may be well-defined if sum_j a_j = 0,
c as may well occur in physical cases.
c Policy is to drop the potentially infinite constant in the transform.
            if (anint(xp).eq.xp.and.anint(xp).le.ZERO) then
c              print *,'fhti: (mu+1+q)/2 =',nint(xp),
c     &          ' is -ve integer, yields singular transform:'
c              print *,'   transform will omit additive constant that is 
c     &generically infinite,'
c              print *,'   but that may be finite or zero if'
c              print *,'   the sum of the elements of the input array a_j
c     & is zero.'
            else
c              print *,'fhti: (mu+1-q)/2 =',nint(xm),
c     &          ' is -ve integer, yields singular inverse transform:'
c              print *,'   inverse transform will omit additive constant 
c     &that is generically'
c              print *,'   infinite, but that may be finite or zero if'
c              print *,'   the sum of the elements of the input array a_j
c     & is zero.'
            endif
            amp=ZERO
            arg=ZERO
          endif
c        neither xp nor xm is a negative integer
        else
          zp=dcmplx(xp,y)
          zm=dcmplx(xm,y)
          zp=cdgamma(zp,1)
          zm=cdgamma(zm,1)
c        Amplitude and Argument of U_mu(q)
          amp=exp(ln2*q+dble(zp)-dble(zm))
c        note +Im(zm) to get conjugate value below real axis
          arg=dimag(zp)+dimag(zm)
        endif
c        cos(arg) = ±1, sin(arg) = 0
        wsave(l+1)=amp*cos(arg)
c........remaining elements of wsave
        d=PI/(n*dlnr)
        do m=1,n/2
c        y = m pi/(n dlnr)
          y=m*d
          zp=dcmplx(xp,y)
          zm=dcmplx(xm,y)
          zp=cdgamma(zp,1)
          zm=cdgamma(zm,1)
c        Amplitude and Argument of kr^(-2 i y) U_mu(q + 2 i y)
          amp=exp(ln2*q+dble(zp)-dble(zm))
          arg=2*ln2kr*y+dimag(zp)+dimag(zm)
          wsave(l+3*m-1)=amp
          wsave(l+3*m)=cos(arg)
          wsave(l+3*m+1)=sin(arg)
        enddo
c Altogether 2*n+18 + 3*(n/2)+1 elements used for q != 0,
c which is (7*n)/2+19 for even n, (7*n)/2+18 for odd n.
c For even n, the very last element of wsave
c [i.e. wsave(l+3*m+1)=sin(arg) for m=n/2] is not used within FFTLog;
c if a low-ringing kr is used, this element should be zero.
c The last element is computed in case somebody wants it.
      endif
  200 ok=.true.
      return
c
c--------error returns
  300 ok=.false.
      return
      end
c
c=======================================================================
      subroutine fftl(n,a,rk,dir,wsave)
      integer n,dir
      real*8 a(n),rk,wsave(*)
c
c This is a driver routine that calls fhtq.
c
c fftl computes a discrete version of the Fourier
c sine (if mu = 1/2) or cosine (if mu = -1/2) transform
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) sin(k r) dr
c                      /
c                     0
c
c                     infinity
c                      /
c    Ã(k) = sqrt(2/pi) | A(r) cos(k r) dr
c                      /
c                     0
c
c by making the substitutions
c                 q-(1/2)                      -q-(1/2)
c    A(r) = a(r) r          and   Ã(k) = ã(k) k    
c
c and applying a biased Hankel transform to a(r).
c
c The steps are:
c 1. a(r) = A(r) r^[-dir*(q-.5)]
c 2. call fhtq to transform a(r) -> ã(k)
c 3. Ã(k) = ã(k) k^[-dir*(q+.5)]
c
c fhti must be called before the first call to fftl,
c with mu = 1/2 for a sine transform,
c or mu = -1/2 for a cosine transform.
c
c A call to fftl with dir=1 followed by
c a call to fftl with dir=-1 (and rk unchanged), or vice versa,
c leaves the array a unchanged.
c
c  Input: n = length of a array.
c         rk = r_c/k_c
c              = r_j/k_j (a constant, the same constant for any j);
c              rk is not (necessarily) the same quantity as kr.
c              rk is used only to multiply the output array by
c              sqrt(rk)^dir, so if you want to do the normalization
c              later, or you don't care about the normalization,
c              you can set rk = 1.
c	  dir = 1 for forward transform,
c		-1 for backward transform.
c               A backward transform (dir = -1) is the same as
c               a forward transform with q -> -q and rk -> 1/rk,
c               for any kr if n is odd,
c               for low-ringing kr if n is even.
c	  wsave = working array set up by fhti.
c Input/Output:
c         a on  input is the array A(r) to transform:
c             a(j) is A(r_j) at r_j = r_c exp[(j-jc) dlnr]
c                where jc = (n+1)/2 = central index of array.
c         a on output is the transformed array Ã(k):
c             a(j) is Ã(k_j) at k_j = k_c exp[(j-jc) dlnr].
c
c        parameters
      real*8 ONE,TWO,HALF
      parameter (ONE=1.d0,TWO=2.d0,HALF=ONE/TWO)
c        local (automatic) variables
      integer j,l
      real*8 dlnr,jc,kr,lnkr,lnrk,q
c
      l=2*n+15
      q=wsave(l+1)
      dlnr=wsave(l+2)
      kr=wsave(l+3)
c........a(r) = A(r) (r/rc)^[-dir*(q-.5)]
c        centre point of array
      jc=dble(n+1)/TWO
      do j=1,n
        a(j)=a(j)*exp(-dir*(q-HALF)*(j-jc)*dlnr)
      enddo
c........transform a(r) -> ã(k)
      call fhtq(n,a,dir,wsave)
c........Ã(k) = ã(k) k^[-dir*(q+.5)] rc^[-dir*(q-.5)]
c        = ã(k) (k/kc)^[-dir*(q+.5)] (kc rc)^(-dir*q) (rc/kc)^(dir*.5)
      lnkr=log(kr)
      lnrk=log(rk)
      do j=1,n
        a(j)=a(j)*exp(-dir*((q+HALF)*(j-jc)*dlnr+q*lnkr-HALF*lnrk))
      enddo
      return
      end
c
c=======================================================================
      subroutine fht(n,a,dir,wsave)
      integer n,dir
      real*8 a(n),wsave(*)
c
c This is a driver routine that calls fhtq.
c
c fht computes a discrete version of the Hankel transform
c
c          infinity
c           /
c    Ã(k) = | A(r) J  (k r) k dr
c           /       mu
c          0
c
c by making the substitutions
c                 q                      -q
c    A(r) = a(r) r    and   Ã(k) = ã(k) k    
c
c and applying a biased Hankel transform to a(r).
c
c The steps are:
c 1. a(r) = A(r) r^(-dir*q)
c 2. call fhtq to transform a(r) -> ã(k)
c 3. Ã(k) = ã(k) k^(-dir*q)
c
c fhti must be called before the first call to fht.
c
c A call to fht with dir=1 followed by
c a call to fht with dir=-1, or vice versa,
c leaves the array a unchanged.
c
c  Input: n = length of a array.
c	  dir = 1 for forward transform,
c		-1 for backward transform.
c               A backward transform (dir = -1) is the same as
c               a forward transform with q -> -q,
c               for any kr if n is odd,
c               for low-ringing kr if n is even.
c	  wsave = working array set up by fhti.
c Input/Output:
c         a on  input is the array A(r) to transform:
c             a(j) is A(r_j) at r_j = r_c exp[(j-jc) dlnr]
c                where jc = (n+1)/2 = central index of array.
c         a on output is the transformed array Ã(k):
c             a(j) is Ã(k_j) at k_j = k_c exp[(j-jc) dlnr].
c
c        parameters
      real*8 ZERO,TWO
      parameter (ZERO=0.d0,TWO=2.d0)
c        local (automatic) variables
      integer j,l
      real*8 dlnr,jc,kr,lnkr,q
c
      l=2*n+15
      q=wsave(l+1)
      dlnr=wsave(l+2)
      kr=wsave(l+3)
c........a(r) = A(r) (r/rc)^(-dir*q)
      if (q.ne.ZERO) then
c        centre point of array
        jc=dble(n+1)/TWO
        do j=1,n
          a(j)=a(j)*exp(-dir*q*(j-jc)*dlnr)
        enddo
      endif
c........transform a(r) -> ã(k)
      call fhtq(n,a,dir,wsave)
c........Ã(k) = ã(k) (k rc)^(-dir*q)
c             = ã(k) (k/kc)^(-dir*q) (kc rc)^(-dir*q)
      if (q.ne.ZERO) then
        lnkr=log(kr)
        do j=1,n
          a(j)=a(j)*exp(-dir*q*((j-jc)*dlnr+lnkr))
        enddo
      endif
      return
      end
c
c=======================================================================
      subroutine fhtq(n,a,dir,wsave)
      integer n,dir
      real*8 a(n),wsave(*)
c
c This is the basic FFTLog routine.
c
c fhtq computes a discrete version of the biased Hankel transform
c
c          infinity
c           /           q
c    ã(k) = | a(r) (k r)  J  (k r) k dr
c           /              mu
c          0
c
c fhti must be called before the first call to fhtq.
c
c A call to fhtq with dir=1 followed by
c a call to fhtq with dir=-1, or vice versa,
c leaves the array a unchanged.
c
c  Input: n = length of a array.
c	  dir = 1 for forward transform,
c		-1 for backward transform.
c               A backward transform (dir = -1) is the same as
c               a forward transform with q -> -q,
c               for any kr if n is odd,
c               for low-ringing kr if n is even.
c Input/Output:
c         a on  input is the periodic array a(r) to transform:
c             a(j) is a(r_j) at r_j = r_c exp[(j-jc) dlnr]
c                 where jc = (n+1)/2 = central index of array.
c         a on output is the transformed periodic array ã(k):
c             a(j) is ã(k_j) at k_j = k_c exp[(j-jc) dlnr].
c
c        parameters
      real*8 ZERO
      parameter (ZERO=0.d0)
c        local (automatic) variables
      integer l,m
      real*8 ai,ar,q
c
      l=2*n+15
      q=wsave(l+1)
      l=l+3
c--------normal FFT
      call drfftf(n,a,wsave)
c--------unbiased (q = 0) transform
      if (q.eq.ZERO) then
c........multiply by
c        (kr)^[- i 2 m pi/(n dlnr)] U_mu[i 2 m pi/(n dlnr)]
        do m=1,(n-1)/2
          ar=a(2*m)
          ai=a(2*m+1)
          a(2*m)=ar*wsave(l+2*m-1)-ai*wsave(l+2*m)
          a(2*m+1)=ar*wsave(l+2*m)+ai*wsave(l+2*m-1)
        enddo
c        problematical last element, for even n
        if (mod(n,2).eq.0) then
          ar=wsave(l+n-1)
c        forward transform: multiply by real part
c Why?  See http://casa.colorado.edu/~ajsh/FFTLog/index.html#ure
          if (dir.eq.1) then
            a(n)=a(n)*ar
c        backward transform: divide by real part
          elseif (dir.eq.-1) then
c Real part ar can be zero for maximally bad choice of kr.
c This is unlikely to happen by chance, but if it does,
c policy is to let it happen.
c For low-ringing kr, imaginary part ai is zero by construction,
c and real part ar is guaranteed nonzero.
            a(n)=a(n)/ar
          endif
        endif
c--------biased (q != 0) transform
      else
c........multiply by
c        (kr)^[- i 2 m pi/(n dlnr)] U_mu[q + i 2 m pi/(n dlnr)]
c        phase
        do m=1,(n-1)/2
          ar=a(2*m)
          ai=a(2*m+1)
          a(2*m)=ar*wsave(l+3*m)-ai*wsave(l+3*m+1)
          a(2*m+1)=ar*wsave(l+3*m+1)+ai*wsave(l+3*m)
        enddo
c        forward transform: multiply by amplitude
        if (dir.eq.1) then
          a(1)=a(1)*wsave(l+1)
          do m=1,(n-1)/2
            a(2*m)=a(2*m)*wsave(l+3*m-1)
            a(2*m+1)=a(2*m+1)*wsave(l+3*m-1)
          enddo
c        backward transform: divide by amplitude
        elseif (dir.eq.-1) then
c        amplitude of m=0 element
          ar=wsave(l+1)
          if (ar.eq.ZERO) then
c Amplitude of m=0 element can be zero for some mu, q combinations
c (singular inverse); policy is to drop potentially infinite constant.
            a(1)=ZERO
          else
            a(1)=a(1)/ar
          endif
c        remaining amplitudes should never be zero
          do m=1,(n-1)/2
            a(2*m)=a(2*m)/wsave(l+3*m-1)
            a(2*m+1)=a(2*m+1)/wsave(l+3*m-1)
          enddo
        endif
c        problematical last element, for even n
        if (mod(n,2).eq.0) then
          m=n/2
          ar=wsave(l+3*m)*wsave(l+3*m-1)
c        forward transform: multiply by real part
          if (dir.eq.1) then
            a(n)=a(n)*ar
c        backward transform: divide by real part
          elseif (dir.eq.-1) then
c Real part ar can be zero for maximally bad choice of kr.
c This is unlikely to happen by chance, but if it does,
c policy is to let it happen.
c For low-ringing kr, imaginary part ai is zero by construction,
c and real part ar is guaranteed nonzero.
            a(n)=a(n)/ar
          endif
        endif
      endif
c--------normal FFT back
      call drfftb(n,a,wsave)
c--------reverse the array
c        and at the same time undo the FFTs' multiplication by n 
      do m=1,n/2
        ar=a(m)
        a(m)=a(n+1-m)/n
        a(n+1-m)=ar/n
      enddo
      if (mod(n,2).eq.1) then
        m=(n+1)/2
        a(m)=a(m)/n
      endif
      return
      end
c
c=======================================================================
      real*8 function krgood(mu,q,dlnr,kr)
      real*8 mu,q,dlnr,kr
c
c Use of this routine is optional.
c
c Choosing kr so that
c     (kr)^(- i pi/dlnr) U_mu(q + i pi/dlnr)
c is real may reduce ringing of the discrete Hankel transform,
c because it makes the transition of this function across the period
c boundary smoother.
c
c  Input: mu = index of J_mu in Hankel transform.
c	  q = exponent of power law bias.
c         dlnr = separation between natural log of points.
c	  kr = suggested value of kr.
c Output: krgood = low-ringing value of kr nearest to input kr.
c                  ln(krgood) is always within dlnr/2 of ln(kr).
c
c        parameters
      real*8 PI
      parameter (PI=3.141592653589793238462643383279502884197d0)
      real*8 ZERO,ONE,TWO
      parameter (ZERO=0.d0,ONE=1.d0,TWO=2.d0)
c        externals
      complex*16 cdgamma
c        local (automatic) variables
      real*8 arg,iarg,xm,xp,y
      complex*16 zm,zp
c
      krgood=kr
      if (dlnr.eq.ZERO) return
      xp=(mu+ONE+q)/TWO
      xm=(mu+ONE-q)/TWO
      y=PI/(TWO*dlnr)
      zp=dcmplx(xp,y)
      zm=dcmplx(xm,y)
      zp=cdgamma(zp,1)
      zm=cdgamma(zm,1)
c        low-ringing condition is that following should be integral
      arg=log(TWO/kr)/dlnr+(dimag(zp)+dimag(zm))/PI
      iarg=anint(arg)
c        if should ensure arg = +-Infinity dealt with correctly
      if (arg.ne.iarg) then
c        low-ringing kr
        krgood=kr*exp((arg-iarg)*dlnr)
      endif
      return
      end
c
c-----------------------------------------------------------------------
c This code was modified from a subroutine from the gamerf package at
c http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html .
c The original copyright statement states:
c   Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
c   You may use, copy, modify this code for any purpose and 
c   without fee. You may distribute this ORIGINAL package.
c
c Permission to distribute this modified gamma function code
c with the FFTLog package has been granted
c (email from Takuya Ooura to Andrew Hamilton dated 16 March 1999).
c
c Original gamerf2a.doc documentation states:
c
c   Gamma(z)=sqrt(2*pi)*(z+r)^(z-1/2)*exp(-z-r)
c            *(a_0
c            +a_1*(z-1)/z
c            +a_2*(z-1)*(z-2)/z/(z+1)
c            +a_3*(z-1)*(z-2)*(z-3)/z/(z+1)/(z+2)
c            +...)
c        a_n= f_n*(2*n)*(2*n-1)*(2*n-2)*...*(n+1)/1/2/3/.../n
c            -a_0*(2*n)*(2*n-1)*(2*n-2)*...*(n+1)/1/2/3/.../n
c            -a_1*(2*n)*(2*n-1)*(2*n-2)*...*(n+2)/1/2/3/.../(n-1)
c            -...
c            -a_(n-1)*(2*n)/1
c        f_n=1/sqrt(2*pi)*(1*2*3*...*n)*(n+1+r)^(-n-1/2)*exp(n+1+r)

c   C.Lanczos,A Precision Approximation of the Gamma Function,
c   J.SIAM Numer.Anal.Ser.B,Vol.1,1964
c
c Modified 28 Oct 98 by Andrew J S Hamilton
c http://casa.colorado.edu/~ajsh/
c (1) to return ln[Gamma(x)] with the correct phase,
c     as well as Gamma(x), and
c (2) to remain accurate for large absolute values of input x,
c     and for x near 0 or negative integers.
c-----------------------------------------------------------------------
      complex*16 function cdgamma(x,l)
      integer l
      complex*16 x
c *
c * Complex Gamma function in double precision.
c *
c * l = 0: Gamma(x)
c *     1: ln[Gamma(x)]
c *
      real*8 pi,pv,pu,pr,p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5,q6
      parameter (
     &    pi = 3.14159265358979324d+00, 
     &    pv = 7.31790632447016203d+00, 
     &    pu = 3.48064577727581257d+00, 
     &    pr = 3.27673720261526849d-02, 
     &    p1 = 1.05400280458730808d+01, 
     &    p2 = 4.73821439163096063d+01, 
     &    p3 = 9.11395751189899762d+01, 
     &    p4 = 6.62756400966213521d+01, 
     &    p5 = 1.32280130755055088d+01, 
     &    p6 = 2.93729529320536228d-01)
      parameter (
     &    q1 = 9.99999999999975753d-01, 
     &    q2 = 2.00000000000603851d+00, 
     &    q3 = 2.99999999944915534d+00, 
     &    q4 = 4.00000003016801681d+00, 
     &    q5 = 4.99999857982434025d+00, 
     &    q6 = 6.00009857740312429d+00)
      real*8 big
      parameter (big=1.d20)
      real*8 t,ui,ur,vi,vr,wr,wi,xi,xr,yi,yr,zero
      data zero /0.d0/
c
      xr = dble(x)
      xi = dimag(x)
c---x = 0, -1, -2
      if (xr.eq.dble(int(xr)).and.(xr.le.0.d0).and.(xi.eq.0.d0)) then
c...Gamma
        if (l .eq. 0) then
          wr = xr / 2.d0
c   +Infinity at even negative integers
          if (wr .eq. dble(int(wr))) then
            yr = 1.d0 / zero
c   -Infinity at odd negative integers
          else
            yr = -1.d0 / zero
          endif
          yi = 0.d0
c...lnGamma
        elseif (l .eq. 1) then
c   real part is +Infinity
          yr = 1.d0 / zero
          yi = pi * int(xr)
        endif
        goto 200
      endif
c---Re(x) < 1/2 : use reflection formula
      if (xr .lt. .5d0) then
          wr = 1.d0 - xr
          wi = -xi
      else
          wr = xr
          wi = xi
      end if
c---large |x|: keep only leading term of rational function
      t = wr * wr + wi * wi
      if (t .gt. big * big) then
c   Rational function v
        vr = wr / t + pr
        vi = - wi / t
c   ln(overall factor)
c   u = ln(x + pv) - 1
        yr = wr + pv
        ur = yr
        if (ur .lt. 0.d0) ur = - ur
        ui = wi
        if (ui .lt. 0.d0) ui = - ui
        if (ur.ge.ui) then
          t = wi / yr
          ur = log(ur) + log(1.d0 + t * t) / 2.d0 - 1.d0
        else
          t = yr / wi
          ur = log(ui) + log(1.d0 + t * t) / 2.d0 - 1.d0
        endif
        ui = atan2(wi, yr)
c---not large |x|
      else
c   u = u(x) = x + q6 = O(x)
        ur = wr + q6
c   v = v(x,u) = (x + q5) u = (x + q5)(x + q6) = O(x^2)
        vr = ur * (wr + q5) - wi * wi
        vi = wi * (wr + q5) + ur * wi
c   y = y(x,u,v) = p6 + p5 u + p4 v = p4 x^2 + ... = O(x^2)
        yr = p6 + (p5 * ur + p4 * vr)
        yi = p5 * wi + p4 * vi
c   u = u(x,v) = (x + q4) v = (x + q4)(x + q5)(x + q6) = O(x^3)
        ur = vr * (wr + q4) - vi * wi
        ui = vi * (wr + q4) + vr * wi
c   v = v(x,u) = (x + q3) u = (x + q3)(x + q4)...(x + q6) = O(x^4)
        vr = ur * (wr + q3) - ui * wi
        vi = ui * (wr + q3) + ur * wi
c   y = y(y,u,v) = y + p3 u + p2 v = p2 x^4 + ... = O(x^4)
        yr = yr + (p3 * ur + p2 * vr)
        yi = yi + (p3 * ui + p2 * vi)
c   u = u(x,v) = (x + q2) v = (x + q2)(x + q3)...(x + q6) = O(x^5)
        ur = vr * (wr + q2) - vi * wi
        ui = vi * (wr + q2) + vr * wi
c   v = v(x,u) = (x + q1) u = (x + q1)(x + q2)...(x + q6) = O(x^6)
        vr = ur * (wr + q1) - ui * wi
        vi = ui * (wr + q1) + ur * wi
c   Numerator
c   y = y(y,u,v) = y + p1 u + v = x^6 + ... = O(x^6)
c     = (x+q1)...(x+q6) + p1 (x+q2)...(x+q6) + ... + p5 (x+q6) + p6
        yr = yr + (p1 * ur + vr)
        yi = yi + (p1 * ui + vi)
c   Denominator
c   u = x v = x(x + q1)(x + q2)...(x + q6) = O(x^7)
        ur = vr * wr - vi * wi
        ui = vi * wr + vr * wi
c   t = |u|^2
        t = ur * ur + ui * ui
c   Rational function v = y u*/|u|^2 + pr = y/u + pr = pr + 1/x + ...
c     = pr + 1/x ( 1 + 1/(x+q1) ( p1 + 1/(x+q2) ( p2 + ...
        vr = (yr * ur + yi * ui) / t + pr
        vi = (yi * ur - yr * ui) / t
c   Overall factor
c   u = ln(x + pv) - 1
        yr = wr + pv
        ur = log(yr * yr + wi * wi) / 2.d0 - 1.d0
        ui = atan2(wi, yr)
      endif
c---lnGamma
c   y = u(x - .5) - pu
c     = (x - .5) [ln(x + pv) - 1] - pu
      yr = ur * (wr - 0.5d0) - ui * wi - pu
      yi = ui * (wr - 0.5d0) + ur * wi
c   y = y + ln(v)
c     = (x - .5) [ln(x + pv) - 1] - pu + ln(Rational)
c     = lnGamma(x)
      yr = yr + log(vr * vr + vi * vi) / 2.d0
      yi = yi + atan2(vi, vr)
c---Reflection formula Gamma(x) Gamma(1-x) = pi/sin(pi x)
c   sign of Gamma
      t = 1.d0
      if (xr .lt. .5d0) then
        wi = anint(xr)
        wr = xr - wi
        if (wi .gt. xr) wi = wi - 1.d0
c   case of real x
        if (xi .eq. 0.d0) then
c   w = ln[sin(pi x)]
          wr = log(sin(pi * abs(wr)))
          if (l .eq. 0) then
            if (wi .ne. 2.d0 * int(wi / 2)) t = -1.d0
            wi = 0.d0
          elseif (l .eq. 1) then
            wi = - pi * wi
          endif
c   case where imaginary part of x is < 1 in absolute value
        elseif (abs(xi) .lt. 1.d0) then
          if (l .eq. 0) then
            if (wi .ne. 2.d0 * int(wi / 2.d0)) t = -1.d0
            ui = 0.d0
          elseif (l .eq. 1) then
            ui = -pi * wi
            if (xi .lt. 0.d0) ui = -ui
          endif
          wr = pi * wr
          wi = pi * xi
          vr = sin(wr) * cosh(wi)
          vi = cos(wr) * sinh(wi)
          if (wr .lt. 0.d0) then
            vr = -vr
            vi = -vi
          endif
c   w = ln[sin(pi x)]
          wr = log(vr * vr + vi * vi) / 2.d0
          wi = ui + atan2(vi, vr)
c   case where imaginary part of x is >= 1 in absolute value
        else
          if (l .eq. 0) then
            if (wi .ne. 2.d0 * int(wi / 2)) t = -1.d0
            if (wr .ge. 0.d0) then
              ui = pi * (.5d0 - wr)
            else
              ui = pi * (- .5d0 - wr)
            endif
          elseif (l .eq. 1) then
            ui = pi * (.5d0 - xr)
          endif
          wi = exp(- 2.d0 * pi * abs(xi))
          wr = 2.d0 * pi * wr
          vr = (1.d0 - cos(wr) * wi) / 2.d0
          vi = - sin(wr) * wi / 2.d0
          ur = pi * xi
c   w = ln[sin(pi x)]
          if (xi .gt. 0.d0) then
            wr = ur + log(vr * vr + vi * vi) / 2.d0
            wi = ui + atan2(vi, vr)
          elseif (xi .lt. 0.d0) then
            wr = - ur + log(vr * vr + vi * vi) / 2.d0
            wi = - ui - atan2(vi, vr)
          endif
        endif
c   y = ln[Gamma(x)]
        yr = log(pi) - yr - wr
        yi = - yi - wi
      endif
c---Gamma
      if (l .eq. 0) then
        ur = exp(yr)
        if (xi .eq. 0.d0) then
          yr = t * ur
          yi = 0.d0
        else
          yr = t * ur * cos(yi)
          yi = ur * sin(yi)
        endif
      endif
c---finish
  200 cdgamma = dcmplx(yr, yi)
      end
c
c-----------------------------------------------------------------------
c http://www.netlib.org/bihar/
c-----------------------------------------------------------------------
      subroutine drfftb (n,r,wsave)
      integer n
c wsave is a work array which should be dimensioned at least 2*n+15
      real*8 r(n), wsave(*)
c
      if (n .eq. 1) return
c
      call drftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine dradb2 (ido,l1,cc,ch,wa1)
      integer ido, l1
      real*8 cc(ido,2,l1), ch(ido,l1,2), wa1(*)
c
      integer i, ic, idp2, k
      real*8 ti2, tr2
c
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
c
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
c
  107 return
      end
c
c-----------------------------------------------------------------------
      subroutine dradb3 (ido,l1,cc,ch,wa1,wa2)
      integer ido, l1
      real*8 cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*)
c
      real*8 TAUI, TAUR
c        -1/2, sqrt(3)/2
      parameter (TAUR = -0.5 d0,
     &  TAUI = 0.8660254037 8443864676 3723170752 93618d0)
c
      integer i, ic, idp2, k
      real*8 ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, ti2, tr2
c
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+TAUR*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = TAUI*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+TAUR*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+TAUR*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = TAUI*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = TAUI*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine dradb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      integer ido, l1
      real*8 cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*), wa3(*)
c
      real*8 SQRT2
c        sqrt(2)
      parameter (SQRT2 = 1.414213562 3730950488 0168872420 970 d0)
c
      integer i, ic, idp2, k
      real*8 ci2, ci3, ci4, cr2, cr3, cr4,
     &  ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
c
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
c
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = SQRT2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -SQRT2*(tr1+ti1)
  106 continue
c
  107 return
      end
c
c-----------------------------------------------------------------------
      subroutine dradb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      integer ido, l1
      real*8 cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*), wa3(*), wa4(*)
c
      real*8 TI11, TI12, TR11, TR12
c        sin(pi/10), sin(2 pi/5), -sin(3 pi/10), sin(pi/5)
      parameter (TR11 = 0.3090169943 7494742410 2293417182 81906d0,
     &  TI11 = 0.9510565162 9515357211 6439333379 38214d0,
     &  TR12 = -0.8090169943 7494742410 2293417182 81906d0,
     &  TI12 = 0.5877852522 9247312916 8705954639 07277d0)
c
      integer i, ic, idp2, k
      real*8 ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5,
     &  di2, di3, di4, di5, dr2, dr3, dr4, dr5,
     &  ti2, ti3, ti4, ti5, tr2, tr3, tr4, tr5
c
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+TR11*tr2+TR12*tr3
         cr3 = cc(1,1,k)+TR12*tr2+TR11*tr3
         ci5 = TI11*ti5+TI12*ti4
         ci4 = TI12*ti5-TI11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
c
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+TR11*tr2+TR12*tr3
            ci2 = cc(i,1,k)+TR11*ti2+TR12*ti3
            cr3 = cc(i-1,1,k)+TR12*tr2+TR11*tr3
            ci3 = cc(i,1,k)+TR12*ti2+TR11*ti3
            cr5 = TI11*tr5+TI12*tr4
            ci5 = TI11*ti5+TI12*ti4
            cr4 = TI12*tr5-TI11*tr4
            ci4 = TI12*ti5-TI11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine dradbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      integer ido, ip, l1
      real*8 cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     &  ch(ido,l1,ip), ch2(idl1,ip), wa(*)
c
      real*8 TPI
c        2 pi
      parameter (TPI = 6.2831853071 7958647692 5286766559 00577d0)
c
      integer i, ic, idij, idl1, idp2, ik, ipph, ipp2, is,
     &   j, jc, j2, k, l, lc, nbd
      real*8 ai1, ai2, ar1, ar1h, ar2, ar2h, arg, dc2, dcp, ds2, dsp
c
      arg = TPI/dble(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
      go to 106
c
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
c
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
c
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
c
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
c
  116 ar1 = 1.
      ai1 = 0.
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
c
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
c
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
c
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
c
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
c
      if (ido .eq. 1) return
      do 133 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  133 continue
c
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
c
      if (nbd .gt. l1) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
c
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140       continue
  141    continue
  142 continue
c
  143 return
      end
c
c-----------------------------------------------------------------------
      subroutine drftb1 (n,c,ch,wa,ifac)
      integer n, ifac(15)
      real*8 c(n), ch(n), wa(n)
c
      integer i, idl1, ido, ip, ix2, ix3, ix4, iw, k1, l1, l2, na, nf
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call dradb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call dradb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
c
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call dradb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call dradb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
c
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call dradb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call dradb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
c
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call dradb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call dradb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
c
  112    if (na .ne. 0) go to 113
         call dradbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call dradbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
c
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c http://www.netlib.org/bihar/
c-----------------------------------------------------------------------
      subroutine drfftf (n,r,wsave)
      integer n
c wsave is a work array which should be dimensioned at least 2*n+15
      real*8 r(n), wsave(*)
c
      if (n .eq. 1) return
c
      call drftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine dradf2 (ido,l1,cc,ch,wa1)
      integer ido, l1
      real*8 cc(ido,l1,2), ch(ido,2,l1), wa1(*)
c
      integer i, ic, idp2, k
      real*8 ti2, tr2
c
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
c
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
c
  107 return
      end
c
c-----------------------------------------------------------------------
      subroutine dradf3 (ido,l1,cc,ch,wa1,wa2)
      integer ido, l1
      real*8 cc(ido,l1,3), ch(ido,3,l1), wa1(*), wa2(*)
c
      real*8 TAUI, TAUR
c        -1/2, sqrt(3)/2
      parameter (TAUR = -0.5 d0,
     &  TAUI = 0.8660254037 8443864676 3723170752 93618d0)
c
      integer i, ic, idp2, k
      real*8 ci2, cr2, di2, di3, dr2, dr3, ti2, ti3, tr2, tr3
c
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = TAUI*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+TAUR*cr2
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+TAUR*cr2
            ti2 = cc(i,k,1)+TAUR*ci2
            tr3 = TAUI*(di2-di3)
            ti3 = TAUI*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine dradf4 (ido,l1,cc,ch,wa1,wa2,wa3)
      integer ido, l1
      real*8 cc(ido,l1,4), ch(ido,4,l1), wa1(*), wa2(*), wa3(*)
c
      real*8 HSQT2
      parameter (HSQT2 = .7071067811 8654752440 0844362104 85 d0)
c
      integer i, ic, idp2, k
      real*8  ci2, ci3, ci4, cr2, cr3, cr4,
     &  ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
c
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
c
      do 106 k=1,l1
         ti1 = -HSQT2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = HSQT2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
c
  107 return
      end
c
c-----------------------------------------------------------------------
      subroutine dradf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      integer ido, l1
      real*8 cc(ido,l1,5), ch(ido,5,l1), wa1(*), wa2(*), wa3(*), wa4(*)
c
      real*8 TI11, TI12, TR11, TR12
c        sin(pi/10), sin(2 pi/5), -sin(3 pi/10), sin(pi/5)
      parameter (TR11 = 0.3090169943 7494742410 2293417182 81906d0,
     &  TI11 = 0.9510565162 9515357211 6439333379 38214d0,
     &  TR12 = -0.8090169943 7494742410 2293417182 81906d0,
     &  TI12 = 0.5877852522 9247312916 8705954639 07277d0)
c
      integer i, ic, idp2, k
      real*8 ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5,
     &  di2, di3, di4, di5, dr2, dr3, dr4, dr5,
     3  ti2, ti3, ti4, ti5, tr2, tr3, tr4, tr5
c
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+TR11*cr2+TR12*cr3
         ch(1,3,k) = TI11*ci5+TI12*ci4
         ch(ido,4,k) = cc(1,k,1)+TR12*cr2+TR11*cr3
         ch(1,5,k) = TI12*ci5-TI11*ci4
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+TR11*cr2+TR12*cr3
            ti2 = cc(i,k,1)+TR11*ci2+TR12*ci3
            tr3 = cc(i-1,k,1)+TR12*cr2+TR11*cr3
            ti3 = cc(i,k,1)+TR12*ci2+TR11*ci3
            tr5 = TI11*cr5+TI12*cr4
            ti5 = TI11*ci5+TI12*ci4
            tr4 = TI12*cr5-TI11*cr4
            ti4 = TI12*ci5-TI11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine dradfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      integer ido, ip, l1
      real*8 cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     &  ch(ido,l1,ip), ch2(idl1,ip), wa(*)
c
      real*8 TPI
c        2 pi
      parameter (TPI = 6.2831853071 7958647692 5286766559 00577d0)
c
      integer i, ic, idij, idl1, idp2, ik, ipph, ipp2, is,
     &   j, jc, j2, k, l, lc, nbd
      real*8 ai1, ai2, ar1, ar1h, ar2, ar2h, arg, dc2, dcp, ds2, dsp
c
      arg = TPI/dble(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do 101 ik=1,idl1
         ch2(ik,1) = c2(ik,1)
  101 continue
      do 103 j=2,ip
         do 102 k=1,l1
            ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
c
      if (nbd .gt. l1) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
c
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108       continue
  109    continue
  110 continue
c
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
c
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
c
  119 do 120 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  120 continue
c
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
c
      ar1 = 1.d0
      ai1 = 0.d0
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125       continue
  126    continue
  127 continue
c
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
c
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
      go to 135
c
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
c
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
c
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138       continue
  139    continue
  140 continue
      return
c
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142       continue
  143    continue
  144 continue
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine drftf1 (n,c,ch,wa,ifac)
      integer n, ifac(15)
      real*8 c(n), ch(n), wa(n)
c
      integer i, idl1, ido, ip, iw, ix2, ix3, ix4,
     &  k1, kh, l1, l2, na, nf
c
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
c
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call dradf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call dradf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
c
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call dradf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call dradf2 (ido,l1,ch,c,wa(iw))
         go to 110
c
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call dradf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call dradf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
c
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call dradf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call dradf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
c
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call dradfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call dradfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
c
  110    l2 = l1
  111 continue
c
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
c
      return
      end
c
c-----------------------------------------------------------------------
c http://www.netlib.org/bihar/
c-----------------------------------------------------------------------
      subroutine drffti (n,wsave)
      integer n
c wsave is a work array which should be dimensioned at least 2*n+15
      real*8 wsave(*)
c
      if (n .eq. 1) return
c
      call drfti1 (n,wsave(n+1),wsave(2*n+1))
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine drfti1 (n,wa,ifac)
      integer n, ifac(15)
      real*8 wa(n)
c
      real*8 TPI
      parameter (TPI = 6.2831853071 7958647692 5286766559 00577d0)
c
      integer i, ib, ido, ii, ip, ipm, is, j, k1, l1, l2, ld,
     &  nf, nfm1, nl, nq, nr, ntry
      integer ntryh(4)
      real*8 arg, argh, argld, fi
c
      data ntryh /4, 2, 3, 5/
c
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
c
      argh = TPI/dble(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = dble(ld)*argh
            fi = 0.d0
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.d0
               arg = fi*argld
               wa(i-1) = dcos(arg)
               wa(i) = dsin(arg)
  108       continue
            is = is+ido
  109    continue
c
         l1 = l2
  110 continue
c
      return
      end
c
