c..............................................................................
      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     f   rkqs,tf)
      implicit none
      integer nbad,nok,nvar,kmaxx,maxstep,nmax
      double precision eps,h1,hmin,x1,x2,ystart(nvar),tiny
      external derivs,rkqs
      parameter (maxstep=100000,nmax=50,kmaxx=200,tiny=1.d-30)
      integer i,kmax,kount,nstp
      double precision dxsav,h,hdid,hnext,x,xsav,dydx(nmax),xp(kmaxx),
     f   y(nmax),yp(nmax,kmaxx),yscal(nmax)
      logical tf
      common /path/ kmax,kount,dxsav,xp,yp
      
      x = x1
      h = sign(h1,x2-x1)
      nok = 0
      nbad = 0
      kount = 0
      do i=1,nvar
         y(i) = ystart(i)
      enddo
      if ( kmax>0 ) xsav = x-2.d0*dxsav
      do nstp = 1,maxstep
         if (nstp.eq.maxstep) tf=.false.
         call derivs(x,y,dydx)
         do i=1,nvar
            yscal(i) = abs(y(i))+abs(h*dydx(i))+tiny
         enddo
         if ( kmax> 0 ) then
            if (abs(x-xsav)>abs(dxsav)) then           
               if ( kount<kmax-1 )then
                  kount = kount+1
                  xp(kount) = x
                  do i=1,nvar
                     yp(i,kount) = y(i)
                  enddo
                  xsav = x
               endif
            endif
         endif
         if ((x+h-x2)*(x+h-x1)>0.d0) h=x2-x

         call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
         if (hdid==h) then
            nok = nok+1
         else
            nbad = nbad+1
         endif
         if ((x-x2)*(x2-x1)>=0.d0) then
            do i=1,nvar
               ystart(i) = y(i)
            enddo
            if ( kmax .ne. 0 ) then
               kount = kount+1
               xp(kount) = x
               do i=1,nvar
                  yp(i,kount) = y(i)
               enddo
            endif
            return
         endif
         if (abs(hnext)<hmin) call error(1,"stepsize smaller than min.")
         h = hnext
      enddo
      call error(1,"too many steps in odeint")
      return 
      end


c..............................................................................
      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      implicit none
      integer n,nmax
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n),
     f   htemp
      external derivs
      parameter (nmax=50)
      integer i
      double precision errmax,h,xnew,yerr(nmax),ytemp(nmax),safety,    
     f   pgrow,pshrnk,errcon
      parameter(safety=0.9d0,pgrow=-0.2d0,pshrnk=-0.25d0,errcon=1.89d-4)
      
      h = htry           ! Set step sixe to initial trial value
 1    call rkck(y,dydx,n,x,h,ytemp,yerr,derivs) !Take a step
      errmax = 0.d0      ! Evaluate accuracy
      do i=1,n
         errmax = max(errmax,abs(yerr(i)/yscal(i)))
      end do 
      errmax = errmax/eps! Scale relative to required tolerance
      if ( errmax>1.d0 ) then !Truncation error too large
         h = safety*h*(errmax**pshrnk)
         if ( h<0.1d0*h) then
            h = 0.1d0*h
         endif
         htemp = safety*h*(errmax**pshrnk)
         h = sign(max(abs(htemp),0.2d0*abs(h)),h) !No more than factor 5
         xnew = x+h
         if ( xnew==x ) call error(11,'stepsize underflow in rkqs') 
         goto 1
      else
         if ( errmax > errcon) then
            hnext = safety*h*(errmax**pgrow)
         else
            hnext = 5.d0*h
         endif
         hdid=h
         x=x+h
         do i=1,n
            y(i) = ytemp(i)
         enddo
         return
      endif
      end

c..............................................................................
      subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)
      implicit none
      integer n,nmax
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      external derivs
      parameter (nmax=6)
      integer i
      double precision ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),
     f    ak6(nmax),ytemp(nmax),a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,
     f    b43,b51,b52,b53,b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,
     f     dc1,dc3,dc4,dc5,dc6 
      parameter (a2=0.2d0,a3=0.3d0,a4=0.6d0,a5=1.d0,a6=0.875d0,
     f    b21=0.2d0,b31=3.d0/40.d0,b32=9.d0/40.d0,b41=0.3d0,
     f    b42=-0.9d0,b43=1.2d0,b51=-11.d0/54.d0,b52=2.5d0,
     f    b53=-70.d0/27.d0,b54=35.d0/27.d0,
     f    b61=1631.d0/55296.d0,b62=175.d0/512.d0,b63=575.d0/13824.d0,
     f    b64=44275.d0/110592.d0,b65=253.d0/4096.d0,c1=37.d0/378.d0,
     f    c3=250.d0/621.d0,c4=125.d0/594.d0,c6=512.d0/1771.d0,
     f    dc1=c1-2825.d0/27648.d0,dc3=c3-18575.d0/48384.d0,
     f    dc4=c4-13525.d0/55296.d0,dc5=-277.d0/14336.d0,dc6=c6-0.25d0)


      do i=1,n                           ! first step
         ytemp(i) = y(i)+b21*h*dydx(i)
      enddo
      
      call derivs(x+a2*h,ytemp,ak2)      ! second step
      do i=1,n
         ytemp(i) = y(i) + h*(b31*dydx(i)+b32*ak2(i))
      enddo
      
      call derivs(x+a3*h,ytemp,ak3)      ! third step
      do i=1,n
         ytemp(i) = y(i) + h*(b41*dydx(i)+b42*ak2(i)+b43*ak3(i))
      enddo
      
      call derivs(x+a4*h,ytemp,ak4)      ! fourth step
      do i=1,n
         ytemp(i) = y(i) + h* ( b51*dydx(i)+b52*ak2(i) +
     f                          b53*ak3(i)+b54*ak4(i)  )
      enddo
      
      call derivs(x+a5*h,ytemp,ak5)     ! fifth step
      do i=1,n
         ytemp(i) = y(i) + h* ( b61*dydx(i)+b62*ak2(i) +
     f                          b63*ak3(i)+b64*ak4(i)+b65*ak5(i) )
      enddo
      
      call derivs(x+a6*h,ytemp,ak6)     ! sixth step
      do i=1,n
         yout(i) = y(i) + h * ( c1*dydx(i)+c3*ak3(i)+c4*ak4(i) 
     f                        + c6*ak6(i) ) 
      enddo
      
      do i=1,n                          ! error estimate
         yerr(i) = h * ( dc1*dydx(i)+dc3*ak3(i)+dc4*ak4(i) 
     f                 + dc5*ak5(i)+dc6*ak6(i) ) 
      enddo
 
      return
      end

c.................................end .........................................