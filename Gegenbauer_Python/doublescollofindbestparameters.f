! Last edited by Chen, Zheng on 11-08-2013

! program to use Gegenbauer expansion 
! to recover exponentially accurate point values
! based on collocation of smooth function with end-point 
! singularity like a(x)+b(x)(1+x)^(p/q)
! using specially designed variable transform, which makes 
! the 'denominator' analytic w.r.t the new varialbe,  then
! recover it with new Gegenbauer Expansion            

! Fourier(even, odd points)/Chebyshev Collocation Method
! (Gauss, Gauss Lobatto, Gauss Radau(including x=1))
! Integration performed by Chebyshev-Gauss or Chebyshev-Gauss-Lobatto
! quadratures

! fort.7 contains:  x
! fort.8 contains:  numerical solution v(x)
! fort.9 contains:  pointwise relative error 
! fort.11.... contains:  pointwise relative error for different 'N's


! This code includes MPFUN90 library with derived types:
! (mp_integer), (mp_real) and (mp_complex).
! In this code, m and lambda are of type integer, not type (mp_integer)

C
C
C      FUNCTION fs(z,p,q) result(fs1)
C! function f_s(z)
C      implicit none
C      real :: z,p,q, fs1,temp
C
C      if(abs(1.d0+z).ge.1.e-16) then
C        temp=(1.d0+z)**(p/q) 
C      else
C        temp=0.d0
C      endif
C! the original function is ff(x)/(1.+x)**\alpha
C!       fs1=1.d0                                     !ex1
C!       fs1=temp                                     !ex2
C!       fs1=1.d0+temp !(0.075,0.075)
C! N=10 m=1 lambda=5
C! N=20 m=1 lambda=10
C! N=40 m= lambda=
C! N=80 m= lambda=
C!       fs1=sin(cos(z))+temp                         !ex3
C!       fs1=cos(z)+sin(z)*temp                       !ex4
C!       fs1=sin(z)+cos(z)*temp                       !ex4
C!       fs1=sin(mpreal(1.4d0)*mppic*(z-mpreal(1.d0)))
C!       fs1=sin(z) !(0.2,0.1)
C!       fs1=sin(z)+temp                              !ex6
C!       fs1=1.+sin(z)*temp 
C      fs1=sin(z)                                !ex5
C      return 
C      end FUNCTION 


      program collo

      include 'COMMON1'
      real :: u(0:1000),er(0:1000),fg(-2000:2000)
      real :: a,alpha,beta,z,ff
      integer :: nzn,nfc,npro,nbegin,nend,npt,nr
      integer :: i,nfile,n,nh,l,lambda,m,nl,Ninput
      real :: ermax,xmax,xs,ex,temp,erav,ermin,xmin
      real :: fs
      real :: errmaxbest,xmaxb,eravb,tf(-2000:2000)
      integer :: iflag,mbest,lbest


      mpsml=1.e-16
      p=1.d0
      q=2.d0

!(0.4,0,1) good for all five cases for fs1=sin(cos(z))+temp 
      alpha=0.4d0
      beta=0.1d0

      write(6,*) 'N = ?',Ninput
      read(*,*) Ninput
      nbegin=Ninput
      nend=Ninput

!       nbegin=10
!       nend=160

      nfc=2
      npro=2
      npt=2
      nr=2
      a=-1.d0
      b=1.d0

      pi=dacos(-1.d0)
      ci = CMPLX(0.d0,1.d0) 

! nzn+1 is the number of points to compute error in [-1,1]
! (causion: those points are not uniform:  more near the boundary)
      nzn=50 

!       write(6,*) 'The singularity (1+x)**(p/q)'
!       write(6,*) 'p = ',int(p)
!      call mpread(6,p)
!       write(6,*) 'q = ',int(q)
!      call mpread(6,q)

      bet=1.d0/q  ! 1./q  or 0.5(win now)
!       write(*,*) 'bet = 1./q = '!,real(bet)
!       call mpwrite(6,bet)
!       write(6,*) ' '

!       write(6,*) '1=Fourier, 2=Chebyshev ', nfc
!      read(*,*) nfc


      if(nfc.eq.1) then
      write(6,*) '1=even(one end point), 2=odd(no boundary pts)', npt
      read(*,*) npt
      else
!       write(6,*) '1=Gauss, 2=G-Lobatto, 3=G-Radau', npt
!      read(*,*) npt
      endif

      write(6,*) '1=Fetch point values from file'
      write(6,*) '2=Calculate pt values'
      write(6,*) '3=calculate Fourier Coeffs -> pt values', nr
      read(*,*) nr


!       write(6,*) 'Choose collocation processure:', npro
!       write(6,*) '1 = Collocate (1-y**2)**(lambda-0.5)f(x)'
!       write(6,*) '2 = (1-y**2)**(lambda-0.5)f(x)(1+x)**((q-1)/2)'
!       write(6,*) '3 = f(x)  not yet edited'
!       read(*,*) npro
!  1-collocate (1-xx**2)**(\lambda-0.5)*ff(yy); 
!  2-collocate (1-xx**2)**(\lambda-0.5)*ff(yy)*(1+xx)**((q-1)/2)
!  2-collocate ff(yy)

!       write(6,*) 'input m ',m
!       read(*,*) m
!       write(6,*) 'input lambda ',lambda
!       read(*,*) lambda
!       write(*,*) ' '
    1 continue

! precomputed gamma function:
      gam(1)=1.d0
      do i=1,799
      gam(i+1)=real(i)*gam(i)
      enddo
!       write(6,*) ' '
!       write(6,*) 'given lambda (alpha) & m (beta):'!,real(alpha),real(beta)
!       call mpwrite(6,alpha,beta)
!       write(6,*) 'the total number of collocation points is:' 
!       write(6,*) 'for Fourier 2N+1, for Chebyshev N+1'
!       write(6,*) ' '
      n=nbegin/2    
      nfile=9+n/20

      open(10,file='results.dat', status='old')
      write(10,*) 'ff(z,p,q)= cos(z)+sin(z)*(mpreal(1.d0)+z)**(p/q)  '
      write(10,*) '(p/q) = ',int(p),'/',int(q)
      write(10,*) 'bet = 1./q = ',bet
      write(10,*) ' '
      write(10,*) '1=Fourier, 2=Chebyshev ', nfc
      write(10,*) 'Choose collocation processure:', npro
      write(10,*) '1 = Collocate (1-y**2)**(lambda-0.5)f(x)'
      write(10,*) '2 = (1-y**2)**(lambda-0.5)f(x)(1+x)**((q-1)/2)'
      write(10,*) '3 = f(x)  not yet edited'
      write(10,*) ' '
      write(10,*) 'given lambda (alpha) & m (beta):',alpha,beta
      write(10,*) 'the total number of collocation points is:' 
      write(10,*) 'for Fourier 2N+1, for Chebyshev N+1'
      write(10,*) ' '

    2 continue      
      n=n*2
! nt is the number of points for quadrature
      nzn=n
      nt=3*n
      nh=n/2
      nfile=nfile+1
      lambda=n
      m=n/2

      x=0.d0
      y=0.d0
      xx=0.d0
      yy=0.d0

!      write(6,*) 'setup done!'

! the collocation points 

      if(nfc.eq.1) then!Fourier collocation points
        if(npt.eq.2) then
          nl=-n+1
        else
          nl=-n
        endif
       if(npt.eq.2) then! odd case
      do i=-n+1,n
       x(i)=real(i)/real(n)
      enddo
       else! even case
      do i=-n,n
       x(i)=real(2*i)/real(2*n+1)
      enddo
       endif
      else!Chebyshev collocation points
       nl=0
       if(npt.eq.1) then! Gauss points(without endpoints)  
      do i=0,n
       x(i)=cos(pi*real(2*i+1)/real(2*n+2))
      enddo
       else if(npt.eq.2) then! Gauss-Lobatto points(both ends)
      do i=0,n
       x(i)=cos(pi*real(i)/real(n))
      enddo
       else! Gauss-Radau points(including one endpoint: x_0=1)
      do i=0,n
       x(i)=cos(pi*real(2*i)/real(2*n+1))
      enddo
       endif
      endif
!       write(6,*) 'collocation points done!'   

! transform 1:y->x  2:x->y 

      call transf(y(nl:n),x(nl:n),nl,n,bet,b,2)
!       write(6,*) 'transf done!'


      lbest=lambda
      mbest=m
      errmaxbest=100.d0
      iflag=1
      eravb=0.d0
      xmaxb=-1.d0

        if(nr.eq.1) then
          if(nfc.eq.1) then
          do i=nl,n
          read(20,*) fg(i) 
          enddo
          else 
          do i=0,n
          read(43,*) fg(n-i) 
          write(6,*) fg(n-i) 
          enddo
          endif
        else if(nr.eq.2) then
          do i=nl,n
            fg(i)=fs(x(i),p,q)
          enddo
        else
!          write(*,*) 'need to edit for calculate FC then f_i'
         call getPointValuefromFC(nfc,nl,n)
         do i=nl,n
            fg(i)=f(i)
         enddo
        endif



      do lambda=1,n/2
!       write(6,*) 'n= ',n,'m= ',m,'lambda= ',lambda
! precomputed (lambda)_m/m!

      gt=0.d0
      gt(0)=1.d0
      do i=1,400
      gt(i)=gam(lambda+i)/gam(lambda)/gam(i+1)
      enddo
!       write(6,*) 'gt s are done!'
! compute f*weight at the collocation points

      do i=nl,n
        f(i)=fg(i)
      enddo
      if(nfc.eq.1) then
      do i=nl,n
       if(npro.eq.1) then
         if(1.d0+y(i).ge.mpsml .and.1.d0-y(i).ge.mpsml) then
          f(i)=f(i)*(1.d0-y(i))**real(lambda-0.5d0)
          f(i)=f(i)*(1.d0+y(i))**real(lambda-0.5d0)
         else
          f(i)=0.d0
         end if
       else
         if(1.d0+y(i).ge.mpsml .and.1.d0-y(i).ge.mpsml) then
          f(i)=f(i)*(1.d0-y(i))**real(lambda-0.5d0)
          f(i)=f(i)*(1.d0+y(i))**(lambda-0.5d0*q)
         else
          f(i)=0.d0
         end if
       end if	 
      enddo
      else
      do i=0,n
       if(npro.eq.1) then
         if(1.d0+y(i).ge.mpsml .and.1.d0-y(i).ge.mpsml) then
          f(i)=f(i)*(1.d0-y(i))**real(lambda-0.5d0)
          f(i)=f(i)*(1.d0+y(i))**real(lambda-0.5d0)
         else
          f(i)=0.d0
         end if
       else
         if(1.d0+y(i).ge.mpsml .and.1.d0-y(i).ge.mpsml) then
          f(i)=f(i)*(1.d0-y(i))**real(lambda-0.5d0)
          f(i)=f(i)*(1.d0+y(i))**(lambda-0.5d0*q)
         else
          f(i)=0.d0
         end if
       end if	 
      enddo
      endif
!       write(6,*) 'f(i)s are done!'

! compute the first 2n+1 Fourier or the first n+1 Chebyshev 
! coefficients of f(x)*weight based on collocation: 

      call ccoe(nfc,n,npt)
!       write(6,*) 'ccoe done!'

! compute the approximate Gegenbauer coefficients

      call gcoe(nfc,npro,n,m,lambda)
!       write(6,*) 'gcoe done!'

! compute the approximation at the nzn+1 chebyshev points

      do i=0,nzn
      u(i) = 0.d0
!        xx(i)=cos(pi*i/nzn)
      enddo
      do l=0,m
!       call transf(yy(0:nzn),xx(0:nzn),0,nzn,bet,b,2)
!       call gegenP(l,nzn,lambda)
      call gegen(l,nzn,lambda)
      call transf(yy(0:nzn),xx(0:nzn),0,nzn,bet,b,1)
      do i=0,nzn
      u(i)=u(i)+hg(l)*cnl(i)
      enddo


      ermax=0.d0
      erav=0.d0
      xmax=-1.d0
      ermin=100.d0
      xmin=-1.d0

!       write(7,*) 'n= ',n,'m= ',m,'lambda= ',lambda
!       write(8,*) 'n= ',n,'m= ',m,'lambda= ',lambda
!       write(9,*) 'n= ',n,'m= ',m,'lambda= ',lambda
!       write(46,*) 'bb=['
!       write(47,*) 'ss=['
      do i=nzn,0,-1
      xs=xx(i)
      ex=fs(xs,p,q)
      er(i)=abs(u(i)-ex)
      temp=er(i)!/abs(ex)  
!         call mpwrite(46,er(i))
!         call mpwrite(47,xs)
      erav=erav+abs(er(i))
      if(temp>ermax) then
        ermax=temp
        xmax=xx(i)
      endif
      if(temp<ermin) then
        ermin=temp
        xmin=xx(i)
      endif
!       write(6,*) ' reconstructing',i,nzn
      enddo
      erav=erav/real(nzn+1)
!       write(8,*) '     '
!       write(46,*) '];'
!       write(47,*) '];'
!       write(6,*) ' '
      write(6,*) 'n= ',n,'m= ',m,'lambda= ',lambda
!       write(6,*) 'L8 error','@','L1 error'
!       call mpwrite(6,ermax)
! !       write(6,*)  'happens at '
!       call mpwrite(6,xmax)
! !       write(6,*) 'L1 error= '
!       call mpwrite(6,erav)
  110 format(5d20.5)
!       write(6,*) ' '

      if(errmaxbest.gt.ermax) then
        errmaxbest=ermax
        lbest=lambda
        mbest=l
        eravb=erav
        xmaxb=xmax
      write(6,*) 'better pair(m,lambda)',mbest,lbest
      write(6,*) errmaxbest
      write(6,*) ' '
      endif
      write(6,*) 'now best pair(m,lambda)',mbest,lbest
      write(6,*) errmaxbest
      write(6,*) ' '
      enddo
      enddo
      write(6,*) 'best pair(m,lambda)',mbest,lbest
      write(6,*) 'L8 error','@','L1 error'
      write(6,*) errmaxbest
      write(6,*)  'happens at '
      write(6,*) xmaxb
      write(6,*) 'L1 error= '
      write(6,*) eravb
      write(6,*) ' '

      if(n==nend) then
        stop
      else 
        goto 2
      endif
      stop
      end

      subroutine getPointValuefromFC(nfc,nl,n)
      include 'COMMON1'
      complex :: fort(-2000:2000),fm,fn,fp
      real :: x1,wt1,wt0,dx
      real :: xm,xn,xp
      integer :: nfc,n,nl,m,i,j
      real :: fs

      x1=sqrt(0.6d0)
      wt1=5.d0/9.d0
      wt0=8.d0/9.d0
      m=min(100*n,4000)
      dx=1.d0/m
      wt1=wt1*dx
      wt0=wt0*dx

      do i=0,n
      fort(i)=0.d0
      do j=-m+1,m-1,2
        xn=(j-x1)*dx
        xm=j*dx
        xp=(j+x1)*dx
        fn=fs(xn,p,q)*exp(-ci*i*pi*xn) 
        fm=fs(xm,p,q)*exp(-ci*i*pi*xm) 
        fp=fs(xp,p,q)*exp(-ci*i*pi*xp) 
        fort(i)=fort(i)+wt1*(fn+fp)+wt0*fm
      enddo
      fort(i)=fort(i)*0.5d0
      fort(-i)=conjg(fort(i))
!       write(6,'(i,5E13.5)') i,real(fort(i)),real(-ci*fort(i))!,real(fort(-i)),real(-ci*fort(-i))
      enddo
!       pause
      do i=nl,n
      f(i)=0.d0
      do j=-n,n
        f(i)=f(i)+fort(j)*exp(ci*j*pi*x(i))
      enddo
!        write(6,'(3E13.5)') x(i),f(i),fs(x(i),p,q)
      enddo
!       pause
      return
      end

      subroutine ccoe(nfc,n,npt)
! subroutine to compute the first 2N+1 Fourier or the first N Chebyshev 
! collocation coefficients, given the point values f(x_i)
      include 'COMMON1'
      integer nfc,n,k,l,nh,npt
      nh=n/2
      
      if(nfc.eq.1) then

      if(npt.eq.2) then!  odd case
        do k=-n,n
        hf(k)=CMPLx(0.d0,0.d0)
        do l=-n+1,n
        hf(k)=hf(k)+f(l)*exp(-ci*real(k)*pi*x(l))
        enddo
        hf(k)=hf(k)/real(2*n)
        enddo
        hf(-n)=hf(-n)*0.5d0
        hf(n) =hf(n) *0.5d0
      else!  even case
        do k=-n,n
        hf(k)=cmplx(0.d0,0.d0)
        do l=-n,n
        hf(k)=hf(k)+f(l)*exp(-ci*real(k)*pi*x(l))
        enddo
        hf(k)=hf(k)/real(2*n+1)
        enddo
      endif

      else 

      if(npt.eq.1) then!  Gauss points
       do k=0,n
       hhf(k)=0.d0
       do l=0,n
       hhf(k)=hhf(k)+f(l)*cos(pi*real((2*l+1)*k)/real(2*n+2))
       enddo
       hhf(k)=hhf(k)*real(2)/real(n+1)
       enddo
       hhf(0)=hhf(0)*0.5d0
       hhf(n)=hhf(n)*0.5d0
      else if(npt.eq.2) then!  Gauss-Lobatto points
       do k=0,n
       hhf(k)=0.d0
       do l=1,n-1
       hhf(k)=hhf(k)+f(l)*cos(pi*real(l*k)/real(n))
       enddo
       hhf(k)=hhf(k)+f(0)*0.5d0
       hhf(k)=hhf(k)+f(n)*real((-1)**k)*0.5d0
       hhf(k)=hhf(k)*real(2)/real(n)
       enddo
       hhf(0)=hhf(0)*0.5d0
       hhf(n)=hhf(n)*0.5d0
      else!  Gauss-Radau points
       do k=0,n
       hhf(k)=0.d0
       do l=1,n
       hhf(k)=hhf(k)+f(l)*cos(pi*real(2*l*k)/real(2*n+1))
       enddo
       hhf(k)=hhf(k)+f(0)*0.5d0
       hhf(k)=hhf(k)*real(4)/real(2*n+1)
       enddo
       hhf(0)=hhf(0)*0.5d0
       hhf(n)=hhf(n)*0.5d0
      end if

      endif

      return
      end



      subroutine gcoe(nfc,npro,n,m,lm)
! subroutine to compute the approximate Gegenbauer coefficients
! using the Chebyshev-Gauss quadrature.
      include 'COMMON1'
      integer nfc,npro,n,nh,lm,m,i,k,l
      real :: ff(0:nt),xr(0:nt),xi(0:nt),fff,z
      real :: temp,temp1,th,hlm,fi,hll
      complex :: ccc
      real :: fs  

! first compute and store f_N(xi_i):
      do i=1,nt-1
      ff(i)=0.d0
      th=real(i)*pi/real(nt)
      xi(i)=cos(th)
      enddo
      call transf(xi(1:nt-1),xr(1:nt-1),1,nt-1,bet,b,1)

      nh=n/2

      if(nfc.eq.1) then
        do k=-n,n
        ccc=ci*real(k)*pi
        do i=1,nt-1
        ff(i)=ff(i)+real(hf(k)*exp(ccc*xr(i)))
        enddo
        enddo
      else
       do k=0,n
       do i=1,nt-1
       ff(i)=ff(i)+hhf(k)*cos(real(k)*acos(xr(i)))
       enddo
       enddo
      end if

        do i=1,nt-1
	temp=fs(xr(i),p,q)*(real(1.d0)-xi(i)**2)**real(lm-0.5d0)
	temp1=abs(temp-ff(i))
	write(7,11110) xr(i),temp,ff(i),temp1
11110   format(5d20.5)
        enddo

      hlm=pi/gam(lm)**2*real((0.5d0)**(2*lm-1))

      do l=0,m
      hg(l)=0.d0
      call gegen(l,nt,lm)
      if(npro.eq.2) then
       do i=1,nt-1
       fi=cnl(i)*ff(i)*(1.d0-xi(i))**0.5d0
       fi=fi*(1.d0+xi(i))**(0.5d0*q)
       hg(l)=hg(l)+fi
       enddo
      else
       do i=1,nt-1
       fi=cnl(i)*ff(i)*(1.d0-xi(i))**0.5d0
       fi=fi*(1.d0+xi(i))**0.5d0
       hg(l)=hg(l)+fi
       enddo
      end if
      hll=hlm*gam(l+2*lm)/real(l+lm)/gam(l+1)
      hg(l)=hg(l)*pi/real(nt)/hll
      enddo

      return
      end


      subroutine gegenP(n,nn,lam)

      include 'COMMON1'
      integer :: i,j,n,nn,lam
      real :: theta,ta,tb,tc,tmp,dt

      do i=0,nn
       ta=1.d0
       tb=2.d0*lam*yy(i)
       if(n.eq.0) then
        tmp=ta
       else if (n.eq.1) then
        tmp=tb
       else
        do j=2,n
         tc=2.d0*(real(j-1)+lam)*yy(i)*tb
         tc=tc-(real(j-2)+2.d0*lam)*ta
         tc=tc/real(j)
         ta=tb
         tb=tc
        enddo
       tmp=tc
      endif
       cnl(i)=tmp
      enddo

      return
      end

      subroutine gegen(n,nn,lam)

      include 'COMMON1'
      integer i,j,n,nn
      real theta
      integer lam

      yy=0.d0
      do i=0,nn
      theta=pi*real(i)/real(nn) 
      yy(i)=cos(theta)
      cnl(i)=real(0.d0)
      do j=0,n
      cnl(i)=cnl(i)+gt(j)*gt(n-j)*cos(theta*real(n-2*j))
      enddo
      enddo
      return
      end


      subroutine transf(y,x,na,nb,bet,b,flag)
! transform 1:y->x  2:x->y
      implicit none
      integer na,nb,flag,i
      real :: x(na:nb),y(na:nb),bet,b,temp,mpsml

      mpsml=1.e-16
      temp=0.5d0*(1.d0+b)**bet
      if(flag==1) then
        do i=na,nb
          if(abs(1.d0+y(i)).ge.mpsml) then
            if(abs(1.d0-y(i)).ge.mpsml) then
            x(i)=(temp*(1.d0+y(i)))**(1.d0/bet)-1.d0
            else
            x(i)=1.d0
            endif
	  else
          x(i)=-1.d0
          endif
        enddo
      else
        do i=na,nb
          if(abs(1.d0+x(i)).ge.mpsml) then
            if(abs(1.d0-x(i)).ge.mpsml) then
            y(i)=-1.d0+1.d0/temp*(1.d0+x(i))**bet
            else
            y(i)=1.d0
            endif
	  else
          y(i)=-1.d0
          endif
        enddo
      endif
      return
      end


