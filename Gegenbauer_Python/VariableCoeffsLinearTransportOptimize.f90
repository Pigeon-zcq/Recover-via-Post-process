! Last edited by Chen, Zheng on 11-08-2013

! fort.7 contains:  x
! fort.8 contains:  numerical solution v(x)
! fort.9 contains:  pointwise relative error 
! fort.11.... contains:  pointwise relative error for different 'N's
!*************************************************************
! edit for calculate FC as a function in x-space
! then reconstructed in the y-space
! find best pair: m lambda
! ifort -r16 VariableCoeffsLinearTransport.f -o av.out
!*************************************************************
      FUNCTION fs(z,p,q) result(fs2)
! function f_s(z)
      implicit none
      real :: z,p,q, fs2,temp, tt, tem

      tem = 1.d0*z+1.d0 + 2.d0
      tt = tem - dble(INT((tem)/2.d0)*2)
      if(abs(tt).ge.1.e-16) then
        temp=tt**(p/q) 

!      if(abs(1.d0+z).ge.1.e-16) then
!        temp=(1.d0+z)**(p/q)        
      else
        temp=0.d0
      endif
!       fs2=1.d0 
!       fs2=cos(z) 
!       fs2=z                                    
       fs2=temp  
!       fs2=tt
!       fs2=cos(z)+sin(z)*temp                                 
      return 
      end FUNCTION 

      FUNCTION fv(z,p,q) result(fs3)
! function f_final(z)
      implicit none
      real :: z,p,q, fs3,temp, tt, tem

      tem = 2.d0*z+1.d0 + 2.d0
      tt = tem - dble(INT((tem)/2.d0)*2)

      if(abs(tt).ge.1.e-16) then
        temp=tt**(p/q) 
      else
        temp=0.d0
      endif
!       fs3=1.d0   
!	  fs3=cos(z)                                  
       fs3=temp    
!       fs3=tt                                 
      return 
      end FUNCTION 

      program collo
      include 'COMMON1'
      real :: u(0:1000),er(0:1000),fg(-2000:2000)
      real :: alpha,beta,z,ff
      integer :: nzn,nfc,npro,nbegin,nend,npt,nr
      integer :: i,nfile,n,nh,l,lambda,m,nl,Ninput
      real :: ermax,xmax,xs,ex,temp,erav,ermin,xmin
      real :: fs
      real :: fv
      real :: errmaxbest,xmaxb,eravb,tf(-2000:2000)
      integer :: iflag,mbest,lbest,nnoi
      real :: noise(-2000:2000), scale

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
      nfc=1
      npro=2
      npt=2
      nr=3
      nnoi=0
      scale = 1.e-6      
      a=-1.d0
      b=1.d0
      pi=dacos(-1.d0)
      ci = CMPLX(0.d0,1.d0) 
      nzn=50 
      bet=1.d0/q  ! 1./q  or 0.5(win now)
      if(nfc.eq.1) then
      write(6,*) '1=even(one end point), 2=odd(no boundary pts)', npt
!      read(*,*) npt
      else
!       write(6,*) '1=Gauss, 2=G-Lobatto, 3=G-Radau', npt
!      read(*,*) npt
      endif
      write(6,*) '1=Fetch point values from file'
      write(6,*) '2=Calculate pt values'
      write(6,*) '3=calculate Fourier Coeffs by definition'
      write(6,*) '4=FG for linear transport equation'
      write(6,*) '  with variable coefficients', nr
      read(*,*) nr
      if(nr.eq.4) then
        a=-0.5d0
        b=0.d0      
      write(6,*) ' To recover the solution in sub-interval [a,b] with a =', a
      read(*,*) a
      write(6,*) ' and b =', b
      read(*,*) b      
      endif

!      write(6,*) '1=add noise, 0=no noise', nnoi
!      read(*,*) nnoi
!      if (nnoi .eq. 1) then
!      write(6,*) 'noise level = ', scale
!      read(*,*) scale
!      endif
            
    1 continue
! precomputed gamma function:
      gam(1)=1.d0
      do i=1,799 !with option -r16, gam(1:35) is computable.
      gam(i+1)=real(i)*gam(i)
!      write(*,'(a,i,d20.5)') 'i,gam(i)',i,gam(i)
      enddo
!      pause
      
      n=nbegin/2    
      nfile=9+n/20

      open(10,file='results.dat', status='old')

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

! the collocation points 
      if(nfc.eq.1) then!Fourier collocation points
        if(npt.eq.2) then
          nl=-n+1
        else
          nl=-n
        endif
       if(npt.eq.2) then! odd case
      do i=-n+1,n
       y(i)=real(i)/real(n)
      enddo
       else! even case
      do i=-n,n
       y(i)=real(2*i)/real(2*n+1)
      enddo
       endif
      else!Chebyshev collocation points
       nl=0
       if(npt.eq.1) then! Gauss points(without endpoints)  
      do i=0,n
       y(i)=cos(pi*real(2*i+1)/real(2*n+2))
      enddo
       else if(npt.eq.2) then! Gauss-Lobatto points(both ends)
      do i=0,n
       y(i)=cos(pi*real(i)/real(n))
      enddo
       else! Gauss-Radau points(including one endpoint: x_0=1)
      do i=0,n
       y(i)=cos(pi*real(2*i)/real(2*n+1))
      enddo
       endif
      endif
! transform 1:y->x  2:x->y 
      call transf(y(nl:n),x(nl:n),nl,n,bet,a,b,1)
      
      lbest=lambda
      mbest=m
      errmaxbest=100.d0
      iflag=1
      eravb=0.d0
      xmaxb=-1.d0
        if(nr.eq.1) then
         write(*,*) 'nr=1'
         stop
        else if(nr.eq.2) then
         write(*,*) 'nr=2'
         stop
        else if(nr.eq.4) then
         write(*,*) 'nr=4, calculate from Fourier Galerkin methods'
         call FourierGalerkinTransport(nfc,nl,n)
        else 
         write(*,*) 'nr=3, calculate FC'
         call getPointValuefromFC(nfc,nl,n)
        endif

        noise = 0.d0
        if(nnoi .eq. 1) then
        call random_number(noise) 
        noise = (noise - 0.5d0)*scale 
        hf(-n:n) = hf(-n:n) + noise(-n:n)
        endif
        
      do lambda=1,n/2
! precomputed (lambda)_m/m!
      gt=0.d0
      gt(0)=1.d0
!      print*,'lambda',lambda
      do i=1,400
      gt(i)=gam(lambda+i)/gam(lambda)/gam(i+1)
      enddo
!      pause
! compute the approximate Gegenbauer coefficients
!      print*, 'mark 0: a,b=' ,a,b
      call gcoe(nfc,npro,n,m,lambda,a,b)
!      print*,' '
! compute the approximation at the nzn+1 chebyshev points
      do i=0,nzn
      u(i) = 0.d0
!        xx(i)=cos(pi*i/nzn)
      enddo
      do l=0,m
!       call transf(yy(0:nzn),xx(0:nzn),0,nzn,bet,a,b,2)
!       call gegenP(l,nzn,lambda)
      call gegen(l,nzn,lambda)
      call transf(yy(0:nzn),xx(0:nzn),0,nzn,bet,a,b,1)
      do i=0,nzn
      u(i)=u(i)+hg(l)*cnl(i)
      enddo
      ermax=0.d0
      erav=0.d0
      xmax=-1.d0
      ermin=100.d0
      xmin=-1.d0
      do i=nzn,0,-1
      xs=xx(i)
      if(nr.eq.4) then
      ex=fv(xs,p,q)
      else
      ex=fs(xs,p,q)
      endif
      er(i)=abs(u(i)-ex)
!      print*,ex,u(i),er(i)
      temp=er(i)!/abs(ex)  
      erav=erav+abs(er(i))
      if(temp>ermax) then
        ermax=temp
        xmax=xx(i)
      endif
      if(temp<ermin) then
        ermin=temp
        xmin=xx(i)
      endif
      enddo
      erav=erav/real(nzn+1)
      write(6,*) 'n= ',n,'m= ',l,'lambda= ',lambda
  110 format(5d20.5)
      if(errmaxbest.gt.ermax) then
        errmaxbest=ermax
        lbest=lambda
        mbest=l
        eravb=erav
        xmaxb=xmax
      write(6,*) 'better pair(m,lambda)',mbest,lbest
      write(6,'(d20.5)') errmaxbest
      write(6,*) ' '
!      pause
      endif
      write(6,*) 'now best pair(m,lambda)',mbest,lbest
      write(6,'(3d20.5)') errmaxbest, ermax, erav
      write(6,*) ' '
!      pause
      enddo
      enddo
      write(6,*) 'best pair(m,lambda)',mbest,lbest
      write(6,'(a,d20.5)') 'L8 error', errmaxbest
      write(6,'(a,d20.5)') 'happens at ', xmaxb
      write(6,'(a,d20.5)') 'L1 error= ', eravb
      write(6,*) ' '

      if(n==nend) then
        stop
      else 
        goto 2
      endif
      stop
      end


      subroutine gcoe(nfc,npro,n,m,lm,aa,bb)
! subroutine to compute the approximate Gegenbauer coefficients
! using the Chebyshev-Gauss quadrature.
      include 'COMMON1'
      integer nfc,npro,n,nh,lm,m,i,k,l
      real :: ff(0:nt),xr(0:nt),xi(0:nt),fff,z
      real :: temp,temp1,th,hlm,fi,hll,aa,bb,er4
      complex :: ccc 
      real :: fv

! first compute and store f_N(xi_i):
      do i=1,nt-1
      ff(i)=0.d0
      th=real(i)*pi/real(nt)
      xi(i)=cos(th)
      enddo
      xr = 0.d0

      call transf(xi(1:nt-1),xr(1:nt-1),1,nt-1,bet,aa,bb,1)

      nh=n/2
      if(nfc.eq.1) then
        do k=-n,n
        ccc=ci*real(k)*pi
        do i=1,nt-1
        ff(i)=ff(i)+real(hf(k)*exp(ccc*xr(i)))
        enddo
        enddo
        do i=1,nt-1            
          temp=fv(xr(i),p,q)
!          write(*,'(i3,5d20.5)') i,xi(i),xr(i),ff(i),temp,abs(ff(i)-temp)
!          write(10,'(i3,5d20.5)') i,xi(i),xr(i),ff(i),temp,abs(ff(i)-temp)
          temp=0.d0
        enddo
!        write(*,*) 'max error in FG methods = ', maxval(abs(1.d0*xr(1:nt-1)-ff(1:nt-1)))
!        write(*,*) 'L1 error in FG methods = ', sum(abs(1.d0*xr(1:nt-1)-ff(1:nt-1)))/(nt-1)        
!        pause
      else
         write(*,*) 'nfc=2'
         stop
      end if
      hlm=pi/gam(lm)**2*real((0.5d0)**(2*lm-1))
      do l=0,m
      hg(l)=0.d0
      call gegen(l,nt,lm)
      if(npro.eq.2) then
       do i=1,nt-1
       fi=cnl(i)*ff(i)*(1.d0-xi(i))**lm
       fi=fi*(1.d0+xi(i))**lm
       hg(l)=hg(l)+fi
!       if(l.eq.m) then
!         print*,'cnl(i),xi(i),ff(i)'
!         write(*,'(3d20.5)') cnl(i),xi(i),ff(i)
!         write(*,'(2i,d20.5)') i,l,hg(l)
!         print*,' '
!       endif
       enddo
!       if(l.eq.m) then
!         pause
!       endif
      else
         write(*,*) 'npro=2'
      end if
      hll=hlm*gam(l+2*lm)/real(l+lm)/gam(l+1)
      hg(l)=hg(l)*pi/real(nt)/hll
!      write(*,'(i,d20.5)') l,hg(l)
      enddo
!      print*,'The above is hg coefficients.n,m,lambda:',n,m,lm
!      pause
      
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
      integer i,j,n,nn,kk
      real theta,gtt
      integer lam
      yy=0.d0
      do i=0,nn
      theta=pi*real(i)/real(nn) 
      yy(i)=cos(theta)
      cnl(i)=real(0.d0)
      do j=0,n
      cnl(i)=cnl(i)+gt(j)*gt(n-j)*cos(theta*real(n-2*j))
!         gtt=1.d0
!         do kk=1,(lam-1)
!            gtt=dble(j+kk)/dble(kk)*dble(n-j+kk)/dble(kk)*gtt
!!            print*,'lam,n,j,kk',lam,n,j,kk
!!            print*,'dble(j+kk),dble(kk),dble(n-j+kk),dble(kk)'
!!            write(*,*) (j+kk),(kk),(n-j+kk),(kk)
!!            print*,' '
!         enddo
!        cnl(i) = cnl(i)+gtt*cos(theta*real(n-2*j))
!        write(*,'(a,i,2d20.5)') 'lam,n,j,gtt',lam,n,j,gtt
!        write(*,'(a,3i,d20.5)') 'lam,n,j,gtt',lam,n,j,gtt
      enddo
!      write(*,'(a,i,d20.5)') 'i,cnl(i) in gegen:',i,cnl(i)
!      print*,' '
      enddo
!      pause
      return
      end


      subroutine transf(y,x,na,nb,bet,a,b,flag)
! transform 1:y->x  2:x->y
      implicit none
      integer na,nb,flag,i
      real :: x(na:nb),y(na:nb),bet,a,b,temp,mpsml
      mpsml=1.e-16
      temp=0.5d0*(b-a)**bet
      if(flag==1) then
        do i=na,nb
          if(abs(1.d0+y(i)).ge.mpsml) then
            if(abs(1.d0-y(i)).ge.mpsml) then
            x(i)=(temp*(1.d0+y(i)))**(1.d0/bet)+a
            else
            x(i)=b
            endif
	      else
            x(i)=a
          endif
!          write(*,'(a,d20.5,d20.5)') 'y(i),x(i)',y(i),x(i)
!          print*,' '
        enddo
      else
        do i=na,nb
          if(abs(x(i)-a).ge.mpsml) then
            if(abs(b-x(i)).ge.mpsml) then
            y(i)=-1.d0+1.d0/temp*(x(i)-a)**bet
            else
            y(i)=1.d0
            endif
	      else
          y(i)=-1.d0
          endif
!          write(*,'(a,d20.5,d20.5)') 'x(i),y(i)',x(i),y(i)   
        enddo
      endif
!      pause
      return
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
!      write(6,'(i3,5E13.5)') i,real(fort(i)),real(-ci*fort(i))!,real(fort(-i)),real(-ci*fort(-i))
      enddo
!      pause
      do i=-n,n
        hf(i)=fort(i)
      enddo
      return
      end
      
! Using Fourier Galerkin method to solve linear transport problem
! u_t - x u_x = 0    x \in [-1,1], t>0 
! i.c.: u(x,0) = \sqrt{1+x}
! with periodic boundary condition with period 2
      subroutine FourierGalerkinTransport(nfc,nl,n)
      include 'COMMON1'
      integer :: nfc, n
      complex  ::  ak(-n:n),fm,fn,fp
      complex  ::  ak1(-n:n),ak2(-n:n),ak3(-n:n),ak4(-n:n)
      real :: AM(-n:n,-n:n), Tfinal, t, dt, cfl,filter(-n:n)
      real :: x1,wt1,wt0,dx
      real :: xm,xn,xp,a1
      integer :: nl,m,i,j,r1
      real :: fs
      
      a1=16.d0*LOG(10.d0)
      r1=int(sqrt(dble(2*n)))
      do i=-n,n
        filter(i) = exp(-a1*(dble(i)/dble(n))**r1)
      enddo
      Tfinal = LOG(2.d0)
      write(*,*) 'Final time Tfinal = ',Tfinal
!      read(*,*) Tfinal      
      cfl = 0.5d0
      write(*,*) 'cfl number for FG = ', cfl
!      read(*,*) cfl
      dt = cfl*Tfinal/dble(2*n+1)
! Initialize the Fourier coefficients ak(-n:n) at time t=0.
      x1=sqrt(0.6d0)
      wt1=5.d0/9.d0
      wt0=8.d0/9.d0
      m=min(100*n,4000)
      dx=1.d0/m
      wt1=wt1*dx
      wt0=wt0*dx

      
      do i=0,n
      ak(i)=0.d0
      do j=-m+1,m-1,2
        xn=(j-x1)*dx
        xm=j*dx
        xp=(j+x1)*dx
        fn=fs(xn,p,q)*exp(-ci*i*pi*xn) 
        fm=fs(xm,p,q)*exp(-ci*i*pi*xm) 
        fp=fs(xp,p,q)*exp(-ci*i*pi*xp) 
        ak(i)=ak(i)+wt1*(fn+fp)+wt0*fm
!        print*,' '
!        write(*,'(3d20.5)') xn, real(fn), real(-ci*fn)
!        write(*,'(3d20.5)') xm, real(fm), real(-ci*fm)
!        write(*,'(3d20.5)') xp, real(fp), real(-ci*fp)
!        write(*,'(2d20.5)') real(wt1*(fn+fp)+wt0*fm),real(-ci*(wt1*(fn+fp)+wt0*fm))
!        write(*,'(i3,a,2d20.5)') i,',',real(ak(i)),real(-ci*ak(i))
!        print*,' '        
      enddo
      ak(i)=ak(i)*0.5d0
      ak(-i)=conjg(ak(i))
!      write(*,'(i3,2d13.5)') i,real(ak(i)),real(-ci*ak(i))
      enddo
! Mass matrix AM      
      AM = 0.d0
      do i=-n,n
        do j=-n,i-1
        AM(i,j) = dble(j)/dble(j-i)*(-1.d0)**(j-i)
        enddo
        do j=i+1,n
        AM(i,j) = dble(j)/dble(j-i)*(-1.d0)**(j-i)
        enddo
      enddo
!      do i=-n,n
!        write(*,*) i,(AM(i,j), j=-n,n)
!        pause
!      enddo
!       write(*,*) 'maxval(AM) = ',maxval(AM)
!      pause
! Solve linear system ak'(t) = AM*ak(t), using RK4.
      t = 0.d0
      j=0
      do while(t<Tfinal)
      if(t+dt>Tfinal) then
        dt = Tfinal-t
      endif
      ak1(-n:n) = MATMUL(AM(-n:n,-n:n),ak(-n:n))
      ak2(-n:n) = MATMUL(AM(-n:n,-n:n),(ak(-n:n)+0.5d0*dt*ak1(-n:n)))
      ak3(-n:n) = MATMUL(AM(-n:n,-n:n),(ak(-n:n)+0.5d0*dt*ak2(-n:n)))
      ak4(-n:n) = MATMUL(AM(-n:n,-n:n),(ak(-n:n)+dt*ak3(-n:n)))
      ak1(-n:n) = ak1(-n:n)+2.d0*ak2(-n:n)+2.d0*ak3(-n:n)+ak4(-n:n)
      ak(-n:n)=ak(-n:n)+dt/6.d0*ak1(-n:n)
!      ak(-n:n)=ak(-n:n)+dt*ak1(-n:n)

      do i=-n,n
        ak(i)=ak(i)*filter(i)
      enddo
      
      t = t+dt
      j=j+1
!      print*,'current time t= ', t
      enddo
      hf = 0.d0
      hf(-n:n) = ak(-n:n)
!      print*,' '
!      print*,'computed Fourier coefficents at final time:'
!      do i=-n,n
!      write(*,'(i5,a,2d20.5)') i,',',real(hf(i)),real(-ci*hf(i))
!      enddo
!      print*,'The above are Fourier coeffs.',j
!      print*,' '
!      pause
      return
      end subroutine FourierGalerkinTransport

      
!      subroutine transf(y,x,na,nb,bet,a,b,flag)
!! transform 1:y->x  2:x->y
!      implicit none
!      integer na,nb,flag,i
!      real :: x(na:nb),y(na:nb),bet,a,b,temp,mpsml
!
!      mpsml=1.e-16
!      temp=0.5d0*(1.d0+b)**bet
!      if(flag==1) then
!        do i=na,nb
!          if(abs(1.d0+y(i)).ge.mpsml) then
!            if(abs(1.d0-y(i)).ge.mpsml) then
!            x(i)=(temp*(1.d0+y(i)))**(1.d0/bet)-1.d0
!            else
!            x(i)=1.d0
!            endif
!	      else
!          x(i)=-1.d0
!          endif
!!          print*,'x,y=',x(i),y(i)
!        enddo
!      else
!        do i=na,nb
!          if(abs(1.d0+x(i)).ge.mpsml) then
!            if(abs(1.d0-x(i)).ge.mpsml) then
!            y(i)=-1.d0+1.d0/temp*(1.d0+x(i))**bet
!            else
!            y(i)=1.d0
!            endif
!	      else
!          y(i)=-1.d0
!          endif
!!          print*,'x,y=',x(i),y(i)
!        enddo
!      endif
!      return
!      end