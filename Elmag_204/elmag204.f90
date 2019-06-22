!=============================================================================!
!=============================================================================!
!          subroutines handling cascade evolution and stack                   !
!=============================================================================!
!=============================================================================!
! (1+1)-dimensional treatment of e/m cascades in extragalactic space          !
!=============================================================================!
subroutine cascade(icq,e00,weight0,z_in)
!-----------------------------------------------------------------------------
! calls: get_particle, propagate, angle_delay, interaction, register, error;
! funcs: t_z, z_t, int_length, themf, eloss, zloss;
! input:
!        icq      - initial particle type (0 - photon, +/-1 - e+/e-);
!        e00      - initial particle energy (eV);
!        weight0  - initial particle weight;
!        z_in     - initial particle redshift position;
!-----------------------------------------------------------------------------
   use constants, only: cy,ame
   use stack, only: jcmb
   use cosmology
   use user_variables
   use internal, only : rcmb,E21
   use bfield_parameters, only : B_rms,L_c
   implicit none
   integer icq,ierr,noint,ncohl
   double precision e00,z_in,weight0
   double precision t_z,z_t,int_length,themf,eloss,zloss
   double precision e0,x,zz,t_c,t,weight,pos,z1,zz0,s
   double precision beta_red,t_in,dt,de,dethr
   double precision y(6),dthet,thex,they,theta,w2
   
   e0=e00
   weight=weight0
   y = 0.d0
   y(6)=vz0 ! y(6) is a dummy variable when using the old calc. routine
   if (old_calc.eq.0) then
      y(4)=vx0
   else
      y(1)=t_z(z_in)
   end if

   s=0.d0
   dt=0.d0

   zz0=z_in
   t_c = t_z(z_in)
   t_in = t_c*t_0
   rcmb = (1.d0-t_c)*t_0*3.d10      !light-travel time from the source (cm)
   !write(*,*) rcmb/(3.086d26)      !length in Mpc
   thex=0.d0
   they=0.d0
   theta=0.d0
   dt=0.d0
   pos=0.d0

   jcmb=0
   noint=0
   do
      ! new particle from stack
      if (jcmb.ne.0.and.noint.eq.0) call get_particle(y,e0,zz0,weight,s,dt,icq)
      if (.not.e0.gt.0.d0) call error('strange e0 in cascade',0)
      noint=0

      if (old_calc.eq.0) then
         pos = sqrt(y(1)**2+y(2)**2+y(3)**2)
      else
         pos = s
      end if

      x=min(rcmb-pos, int_length(e0,zz0,icq))      ! interaction length
      if (x.lt.0.d0) call error('interaction length x is negative',0)

      t = t_c+(pos+x)/(t_0*3.d10)
      zz = z_t(t)                                ! new redshift

      if (zz0-zz.gt..05d0) then                  ! dz_max = 0.05 
         zz=zz0-.05d0
         if (zz.lt.0.d0) zz=0.d0
         noint=1
         t=t_z(zz)
         x=(t-t_z(zz0))*t_0*3.d10
      end if

      z1 = 1.d0+zz
      beta_red = e0*H_0*sqrt(Omega_m*z1**3.d0 + Omega_v)
      e0 = e0-x/3.d10*beta_red

      if (icq.ne.0) then
         if (old_calc.eq.0) then
            E21 = icq*e0/1.d21                     ! variable used in derivs
            call propagate(y,x,e0)
         else                                      ! use routine from Elmag203
            dthet=themf(e0,B_rms)                  ! deflection angle per cm
            if (x+y(2).lt.L_c) then
               y(4)=y(4)+dthet*x                   ! accumulate deflection angle
               y(2)=y(2)+x
            else
               y(4)=y(4)+dthet*(L_c-y(2))          ! accumulate deflection angle
               y(5)=y(5)+y(4)**2                   ! accumulate deflection angle^2
               y(2)=x+y(2)-L_c
               if (y(2).ge.L_c) then               ! account for random B-field
                  ncohl=int(y(2)/L_c)
                  y(5)=y(5)+(dthet*L_c)**2*ncohl   ! accumulate deflection angle^2
                  y(2)=y(2)-L_c*ncohl
               end if
               y(4)=dthet*y(2)                     ! accumulate deflection angle
            end if
            y(3)=s+x
         end if
         dt=dt+x/cy*(ame/e0)**2/2.d0               ! kinematic time delay
         de=min(e0,eloss(e0,B_rms)*x)              ! local magnetic field
         dethr=min(e0-de,e0*zloss(e0,zz0)*x)       ! synchrotron E-loss
         e0=e0-de-dethr
      elseif (old_calc.eq.0) then
         y(1) = y(1) + y(4)*x
         y(2) = y(2) + y(5)*x
         y(3) = y(3) + y(6)*x
      end if

      s=s+x
      if (old_calc.eq.0) then
         pos = sqrt(y(1)**2+y(2)**2+y(3)**2)
      else
         pos = s
      end if

      ! if particle is deflected to negative z-vel => discard it
      if (y(6).lt.0.d0) then ! discard particle
         noint=0
      elseif (e0.le.ethr.or.pos.ge.0.999d0*rcmb) then ! end of tracking
         if (icq.eq.0) then
            if (old_calc.eq.1) then
               call angle_delay_1d(e0,y(3),y(5),thex,dt)
               they=0.d0
               w2=1.d0
            else
               call angle_delay_3d(y,e0,s,thex,they,dt,w2)
            end if
         end if
         if (e0.gt.ethr) call register(e0,thex,they,weight*w2,dt,icq)
         noint=0
      elseif (noint.eq.0) then      ! interaction with EBL
         call interaction(y,e0,zz,weight,s,dt,icq,ierr)
      else
         zz0 = zz
      end if

      if (noint.eq.0.and.jcmb.eq.0) exit                  ! empty stack
   end do

end subroutine cascade
!=============================================================================!
!=============================================================================!
subroutine propagate(y,x,e0)
!-----------------------------------------------------------------------------
! calls: odeint, normalizer;
! input:
!        y     - particle position and velocity;
!        x     - particle travel length before reaction;
! output:
!        y     - new particle position and velocity;
!        e0    - particle energy;
!-----------------------------------------------------------------------------
   use user_result, only : n_kc
   use bfield_parameters, only : B_rms,L_c
   implicit none
   integer nok,nbad
   integer, parameter :: nvar=6
   double precision y(6),x,e0
   double precision lar_rad,xcorr,eps,h1,hmin
   logical tf
   external derivs,rkqs

   eps=1.d-5

   ! correcting travel length
   lar_rad=e0*1.081d0*3.086d0/(B_rms*2.d0/3.d0)/1.d3
   xcorr=lar_rad*sin(x/lar_rad)
   
   h1=L_c/1.d3                   ! Starting step size
   hmin=L_c/1.d8                 ! Minumum step size
   tf=.true.
   call odeint(y,nvar,0.d0,xcorr,eps,h1,hmin,nok,nbad,derivs,rkqs,tf)
   call normalizer(y)
   if (tf.eqv..false.) then   ! the path could not be solved
      e0=0.d0                 ! exclude particle
      n_kc=n_kc+1             ! count bad particles
   end if
end subroutine propagate
!=============================================================================!
!           Normalizer for the propagation vector                             !
!           To make sure total velocity does not exceed c                     !
!=============================================================================!
subroutine normalizer(y)
   implicit none
   double precision y(6),n_fac

   n_fac=sqrt(y(4)**2+y(5)**2+y(6)**2)
   y(4)=y(4)/n_fac
   y(5)=y(5)/n_fac
   y(6)=y(6)/n_fac
end subroutine normalizer
!=============================================================================!
!           Derivatives to use for the numerical solver (odeint)              !
!=============================================================================!
subroutine derivs(x,y,dydx)
!-----------------------------------------------------------------------------
! calls: generate_B;
! input:
!        y     - particle position and velocity;
! output:
!        dydx  - derivatives of y;
!
! NB: this subroutine is called by rkck() from dgl_numrec.f which is the
!     numerical solver for the propagtion in the turbulent B-field
!-----------------------------------------------------------------------------
   use internal, only : E21
   implicit none
   double precision x,y(6),dydx(6)
   double precision Omega(3),s(3)
   double precision, parameter :: scale=1.d0/(1.081d0*3.086d18)

   ! x-variable should only be used for time dependent B-field
   dydx=0.d0
   s(1)=y(1)
   s(2)=y(2)
   s(3)=y(3)
   call generate_B(s,Omega)   ! retrieve magnetic field
   
   dydx(1) = y(4)
   dydx(2) = y(5)
   dydx(3) = y(6)
   dydx(4) = scale*(y(5)*Omega(3) - y(6)*Omega(2))/E21
   dydx(5) = scale*(y(6)*Omega(1) - y(4)*Omega(3))/E21
   dydx(6) = scale*(y(4)*Omega(2) - y(5)*Omega(1))/E21

end subroutine derivs
!=============================================================================!
! calculation of deflection angles/ time delays for observed gammas           !
!=============================================================================!
subroutine angle_delay_3d(y,e0,s,thetax,thetay,dt,w2)
!-----------------------------------------------------------------------------
! calls: jet_distribution;
! input:
!        y        - pos/vel variable;
!        e0       - energy of photon;
!        s        - total travel length;
! output:
!        thetax   - observation angle (degree) in x-dir for the photon;
!        thetay   - observation angle (degree) in y-dir for the photon;
!        dt       - time delay;
!        w2       - weight from jet distribution;
!-----------------------------------------------------------------------------
   use user_variables, only : th_jet,th_jetx
   use constants, only : rad_degree,degree_rad,cy
   use internal, only : rcmb
   implicit none
   double precision y(6),e0,s,thetax,thetay,dt,w2
   double precision a,b,c,f1,f2
   double precision thex,they,thevx,thevy,the
   double precision xa,xi,yi,the_s,vx,vy

   thetax=0.d0
   thetay=0.d0
   the=0.d0
   w2=1
   if (y(1).eq.0.d0.and.y(2).eq.0.d0) then
      s=rcmb
      return
   end if
   ! Correction to the geometry on the sphere
   a = y(4)**2+y(5)**2+y(6)**2
   b = 2*(y(1)*y(4)+y(2)*y(5)+y(3)*y(6))
   c = y(1)**2+y(2)**2+y(3)**2-rcmb**2
   f1 = (-b+sqrt(b**2-4*a*c))/(2*a)       ! sol. to the 2nd order eq
   f2 = (-b-sqrt(b**2-4*a*c))/(2*a)       ! sol. to the 2nd order eq
   if (abs(f2).lt.abs(f1)) f1=f2          ! choose wichever is smallest
   y(1) = y(1) + y(4)*f1                  ! new positions
   y(2) = y(2) + y(5)*f1
   y(3) = y(3) + y(6)*f1
   if (y(3).lt.0.d0) then
      e0=0.d0
      return
   end if
   s=s+f1                                 ! corrected traveled distance
   dt=dt+max(0.d0,(s-rcmb)/cy)

   the=atan(sqrt(y(1)**2+y(2)**2)/y(3))   ! small angle apprx
   if (the*rad_degree.gt.th_jet) then
      e0=0.d0
      return
   end if
   ! Calculating observed angles
   xa=th_jetx/2.d0*degree_rad          ! x-axis to mirror
   thex=atan(y(1)/y(3))
   they=atan(y(2)/y(3))      
   xi=2*xa-thex                        ! starting angle in x-dir
   yi=-they                            ! starting angle in y-dir
   the_s=sqrt((xi-2*xa)**2+yi**2)*rad_degree   !angle from jet center
   call jet_distribution(the_s,w2)     ! weight from jet distribution
   thevx=atan(y(4)/y(6))
   thevy=atan(y(5)/y(6))
   vx=thevx+(xi-2*xa)                  ! correct observed angle in x-dir
   vy=thevy+yi                         ! correct observed angle in y-dir
   thetax=-vx*rad_degree               ! observed angle in x-dir
   thetay=-vy*rad_degree               ! observed angle in y-dir

end subroutine angle_delay_3d
!=============================================================================!
! calculation of deflection angles/ time delays for observed gammas           !
!=============================================================================!
subroutine angle_delay_1d(e0,xx,the,theta,dt)
!-----------------------------------------------------------------------------
! calls: jet_distribution;
! input:
!        y        - pos/vel variable;
!        e0       - energy of photon;
!        s        - total travel length;
! output:
!        theta    - observation angle (degree) total for the photon;
!        dt       - time delay;
!-----------------------------------------------------------------------------
   use user_variables, only : th_jet
   use constants, only : rad_degree,degree_rad,pi,cy
   use internal, only : rcmb
   use cosmology, only : t_0
   implicit none
   double precision e0,xx,the,theta,dt
   double precision cbeta,psran,beta,sbeta

   if (the.gt.pi**2/4.d0) then     !beyond small-angle approx. -> take isotropic
      dt=t_0
      cbeta=2.d0*psran()-1.d0 
      beta=acos(cbeta)*rad_degree
      sbeta=sqrt(1.d0-cbeta**2)
      theta=asin(xx/rcmb*sbeta)*rad_degree               ! observation angle
      if (beta-theta.gt.th_jet) e0=0.d0                    ! dismiss photon
      return
   elseif (the.eq.0.d0) then
      theta=0.d0
   else
      sbeta=sin(sqrt(the))
      theta=asin(xx/rcmb*sbeta)*rad_degree               ! observation angle
      if (sqrt(the)*180.d0/pi-theta.gt.th_jet) e0=0.d0    ! dismiss photon
      dt=dt+2.d0*xx/cy*(1.d0-xx/rcmb)*sbeta**2         ! time delay
   end if
   if (theta<0.d0) call error('theta<0 in angle_delay',1)

end subroutine angle_delay_1d
!=============================================================================!
!=============================================================================!
! interaction with background photons                                         !
!=============================================================================!
subroutine interaction(y,e0,zz,weight,s,dt,icq,ierr)
!-----------------------------------------------------------------------------
! calls: sample_photon, sample_electron, zpair, zics, store_particle, error
! input:
!        y        - particle position and velocity;
!        e0       - particle energy;
!        zz       - particle redshift position;
!        weight   - particle weight;
!        s        - total traveled distance for the particle;
!        dt       - particle time delay;
!        icq      - particle type (0 - photon, +/-1 - electron/positron);
! output:
!        ierr     - error code (0 - o.k., 1 - error);
!-----------------------------------------------------------------------------
   use user_variables, only : old_calc
   implicit none
   integer icq,ierr
   double precision e0,y(6),zz,weight,s,dt,z,sgam
   double precision zpair,zics

   ierr=0   
   select case (icq) 
   case (0)                                   ! pair production on EBL
      call sample_photon(e0,zz,sgam,ierr)     ! sample c.m. energy for interaction
      if (ierr.eq.1) return
      z=zpair(sgam)                           ! E-share taken by electron
      if (old_calc.eq.1) y(4)=0.d0
      call store_particle(y,z*e0,zz,z,weight,s,dt,1) ! record e- 
      call store_particle(y,(1.d0-z)*e0,zz,(1.d0-z),weight,s,dt,-1)
   case (-1,1)                                ! ICS on EBL
      call sample_electron(e0,zz,sgam,ierr)
      if (ierr.eq.1) return
      z=zics(e0,sgam)                         ! E-share taken by electron/positron
      call store_particle(y,z*e0,zz,z,weight,s,dt,icq)               ! e+-
      if (old_calc.eq.1) then
         y(5)=y(5)+y(4)**2
         y(4)=0.d0
      end if
      call store_particle(y,(1.d0-z)*e0,zz,(1.d0-z),weight,s,dt,0)   ! gam
   case default
      call error('wrong icq in interaction',0)
   end select

end subroutine interaction
!=============================================================================!
!=============================================================================!
!        add a particle to stack                                              !
!=============================================================================!
subroutine store_particle(y,e0,zz,ze,weight,s,dt,icq)
!-----------------------------------------------------------------------------
! calls: psran, error
! input:
!        e0       - particle energy;
!        y        - particle position and velocity;
!        zz       - particle redshift position;
!        ze       - share of the parent particle energy;
!        weight   - particle weight;
!        s        - total traveled distance for the particle;
!        dt       - particle time delay;
!        icq      - particle type (0 - photon, +/-1 - electron/positron);
!-----------------------------------------------------------------------------
   use stack
   use user_variables, only : ethr,a_smp
   implicit none
   integer icq,i,n
   double precision y(6),e0,zz,ze,weight,s,dt
   double precision aweight,psran

   if (e0.le.ethr) return              ! disregard underthreshold particles
   aweight=ze**a_smp                   ! sampling weight for produced particle
   if (psran().gt.aweight) return      ! disregard particle with prob. (1-aweight)

   jcmb=jcmb+1                         ! enlarge stack
   if (jcmb.eq.n_max) call error('enlarge storage!',0)

   act = one_event(y,e0,zz,weight,s,dt,icq)

   if (jcmb.gt.1) then                 ! add particle to E-ordered stack
      do i=1,jcmb-1
         n=jcmb-i
         if(event(n)%en.gt.e0)then
            event(n+1)=act                     ! add particle to stack
            event(n+1)%w=event(n+1)%w/aweight  ! total weight of the particle
            return
         else
            event(n+1)=event(n)                ! re-arrange stack
         endif
      enddo
   endif
   event(1)=act                        ! 1st particle in stack
   event(1)%w=event(1)%w/aweight       ! total weight of the particle

end subroutine store_particle
!=============================================================================!
!=============================================================================!
!        get a particle from stack                                            !
!=============================================================================!
subroutine get_particle(y,e0,zz,weight,s,dt,icq)
!-----------------------------------------------------------------------------
! output:
!        e0       - particle energy;
!        y        - particle position and velocity;
!        zz       - particle redshift position;
!        weight   - particle weight;
!        s        - total traveled distance for the particle;
!        dt       - particle time delay;
!        icq      - particle type (0 - photon, +/-1 - electron/positron);
!-----------------------------------------------------------------------------
   use stack
   implicit none
   integer icq
   double precision y(6),e0,zz,weight,s,dt

   y=event(jcmb)%y
   e0=event(jcmb)%en
   zz=event(jcmb)%z
   weight=event(jcmb)%w
   s=event(jcmb)%s
   dt=event(jcmb)%dt
   icq=event(jcmb)%icq

   jcmb=jcmb-1

end subroutine get_particle
!=============================================================================!
!=============================================================================!
!            subroutines handling single interaction                          !
!=============================================================================!
!=============================================================================!
! sample c.m. energy for gamma-gamma_EBL interaction                          !
!=============================================================================!
subroutine sample_photon(e0,zz,sgam,ierr)
!-----------------------------------------------------------------------------
! calls: w_EBL_density, sigpair, psran
! input:
!        e0 - photon energy;
!        zz - photon redshift position
! output:
!        sgam - c.m. energy for gamma-gamma interaction;
!        ierr - error code (0 - o.k., 1 - error)
!
! NB: rejection method used - optimized for the EBL models of Kneiske & Doll;
! the procedure may have to be re-adjusted for different (higher EBL) models
!-----------------------------------------------------------------------------
   use EBL_fit
   use constants
   use user_variables, only : ethr,model
   use internal, only : debug
   implicit none
   integer ierr,nrej,nrejmax,nrenorm
   double precision e0,zz,sgam,de,emin,emin0 &
   ,etrans1,etrans2,aw1,aw2,aw3,gb,gb0,gbmax,gbnorm,rrr,gnorm1,gnorm2,gnorm3
   double precision psran,sigpair,w_EBL_density

   de=4.d0*e0  
   emin0=ame**2/e0                  ! minimal required energy for EBL photon
   if (emin0.ge.eirmax) then        ! wrong kinematics
      !!!     write(*,*)'photon:',emin0,eirmax,e0
      ierr=1
      return
   end if

   nrej=0
   nrejmax=3000                     ! user-defined limit on the N of rejections
   nrenorm=0
   gbmax=0.d0
   gbnorm=2.5d0                      ! normalization factor for rejection

   etrans1=1.d-6                    ! parameters for 'proposal function'
   etrans2=eirmin*(1.d0+zz)**1.25

   ! partial weights for different energy intervals for 'proposal function'  
   if(emin0.lt.etrans1)then
      gnorm1=w_EBL_density(etrans1,zz)*etrans1*sigpair(etrans1*de)
      aw1=gnorm1*(etrans1/emin0-1.d0)*etrans1
   else
      aw1=0.d0
   endif

   if(emin0.lt.etrans2)then
      gnorm2=w_EBL_density(max(etrans1,2.d0*emin0),zz)*max(etrans1,2.d0*emin0) &
      *sigpair(max(etrans1,2.d0*emin0)*de)
      aw2=gnorm2*(exp((max(etrans1,2.d0*emin0)-max(etrans1,emin0)) &
      /akt/(1.d0+zz)/1.5)-exp((max(etrans1,2.d0*emin0)-etrans2) &
      /akt/(1.d0+zz)/1.5))*akt*(1.d0+zz)*1.5d0
   else
      aw2=0.d0
   endif

   gnorm3=w_EBL_density(eirmax/2.d0,zz)*eirmax*sigpair(eirmax*de)
   aw3=gnorm3*((eirmax/max(emin0,etrans2))**2.5d0-1.d0)/2.5d0*eirmax

   1 rrr=psran()*(aw1+aw2+aw3)

   ! sample emin (= sgam/de) according to the 'proposal function';
   ! define 'rejection function' ( gb = f(emin) / f_proposal(emin) )
   if(rrr.lt.aw1)then
      emin=etrans1/(1.d0+psran()*(etrans1/emin0-1.d0))
      gb0=gnorm1*(etrans1/emin)**2
   elseif(rrr.lt.aw1+aw2)then
      emin=etrans2-akt*(1.d0+zz)*1.5d0*dlog(1.d0-psran() &
      *(1.d0-exp((etrans2-max(emin0,etrans1))/akt/(1.d0+zz)/1.5d0)))
      gb0=gnorm2*exp((max(etrans1,2.d0*emin0)-emin)/akt/(1.d0+zz)/1.5d0)
   else
      emin=eirmax/(1.d0-psran()*(1.d0-(eirmax/max(emin0,etrans2))**2.5d0))**.4d0
      gb0=gnorm3*(eirmax/emin)**3.5d0
   endif
   gb0=gb0*gbnorm

   if(model.eq.1)then
      gb0=gb0*2.d0
      if(emin0/etrans2.gt.1.d0)gb0=gb0/10.d0
   elseif(model.eq.2)then
      gb0=gb0*20.d0
      if(emin0/etrans2.lt.1.d0)gb0=gb0*30.d0
   elseif(model.eq.3)then
      gb0=gb0*2.d0
      if(emin0/etrans2.gt.1.d0)gb0=gb0*2.2d1
      if(emin0/etrans2.gt.1.d2)gb0=gb0*6.d0
   elseif(model.eq.4)then
      gb0=gb0*4.5d0
      if(emin0/etrans2.gt.1.d0)gb0=gb0/20.d0
      if(emin0/etrans2.gt.1.d2)gb0=gb0*1.5d0
   elseif(model.eq.5)then
      gb0=gb0*2.d0
      if(emin0/etrans2.gt.1.d0)gb0=gb0/20.d0
      if(emin0/etrans2.gt.1.d2)gb0=gb0*1.5d0
   elseif(model.eq.6)then
      gb0=gb0*2.d0
      if(emin0/etrans2.gt.1.d0)gb0=gb0/10.d0
      if(emin0/etrans2.gt.1.d2)gb0=gb0*1.2d0
   endif
   !

   sgam=emin*de                 ! c.m. energy for gamma-gamma interaction
   gb=w_EBL_density(emin,zz)*sigpair(sgam)*emin/gb0

   !  if (gb.gt.1.d0.and.nrenorm.eq.0) write(*,*)'sample_cmb(photon): gb=' &
   !  ,gb,nrenorm,emin,emin0,emin0/etrans2  !/1.d3

   if (psran().gt.gb) then       ! rejection
      nrej=nrej+1                ! total number of rejections for current sampling
      gbmax=max(gbmax,gb)        ! maximal value for rejection function
      if(nrej.gt.nrejmax)then    ! too many rejections
         if(gbmax.le.0.d0)then     ! wrong kinematics
            write(*,*)'photon: gbmax=0!!!'
            ierr=1
            return
         else
            !       write(*,*)'nrej(gamma)>nrejmax',nrej,emin0/etrans2,nrenorm,e0/1.d12,gbmax
            gbnorm=gbnorm*gbmax*2.d0 ! change normalization for the rejection function
            gbmax=0.d0
            nrenorm=nrenorm+1
            nrej=0
         endif
      endif
      goto 1                     ! new try
   end if

end subroutine sample_photon

!=============================================================================!
!=============================================================================!
! sample c.m. energy for e(+-)-gamma_EBL interaction                          !
!=============================================================================!
subroutine sample_electron(e0,zz,sgam,ierr)
!-----------------------------------------------------------------------------
! calls: w_EBL_density, sigics, psran
! input:
!        e0 - electron/positron energy;
!        zz - electron/positron redshift position
! output:
!        sgam - c.m. energy for e(+-)-gamma interaction;
!        ierr - error code (0 - o.k., 1 - error)
!
! NB: rejection method used - optimized for the EBL models of Kneiske & Doll;
! the procedure may have to be re-adjusted for different (higher EBL) models
!-----------------------------------------------------------------------------
   use EBL_fit
   use constants
   use user_variables, only : ethr,model
   use internal, only : debug
   implicit none
   integer nrej,ierr,nrejmax,nrenorm
   double precision e0,zz,sgam,de,emin,emin0,gb,gb0,gbmax,gbnorm   &
   & ,etrans1,etrans2,etrans3,aw1,aw2,aw3,aw4,rrr,gnorm1,gnorm2,gnorm3,gnorm4
   double precision psran,sigics,w_EBL_density

   de=2.d0*e0*(1.d0+dsqrt(max(0.d0,1.d0-ame**2/e0**2)))
   emin0=ame**2*ethr/(e0-ethr)/de    ! minimal required energy for EBL photon

   if (emin0.ge.eirmax) then         ! wrong kinematics
      !!!     write(*,*)'electron:',emin0,eirmax,e0
      ierr=1
      return
   end if

   nrej=0
   nrejmax=3000                      ! user-defined limit on the N of rejections
   nrenorm=0
   gbmax=0.d0
   gbnorm=3.d0                       ! rejection normalization factor

   etrans1=2.d-10                    ! parameters for 'proposal function'
   etrans2=4.d-7
   etrans3=eirmin*(1.d0+zz)**1.2

   ! partial weights for different energy intervals for 'proposal function'  
   if(emin0.lt.etrans1)then
      gnorm1=emin0*w_EBL_density(emin0,zz)*sigics(e0,ame**2+2.d0*emin0*de)
      aw1=gnorm1*((etrans1/emin0)**(brad1+1.d0)-1.d0)/(brad1+1.d0)*emin0
   else
      aw1=0.d0
   endif

   if(emin0.lt.etrans2)then
      gnorm2=etrans2*w_EBL_density(etrans2,zz)*sigics(e0,ame**2+etrans2*de)
      if(etrans2*de.gt.10.d0*ame**2)then
         gnorm2=max(gnorm2,w_EBL_density(etrans1,zz)*sigics(e0,ame**2 &
         +etrans1*de)*etrans1*(etrans1/etrans2)**(1.d0-brad2))
         aw2=gnorm2*(1.d0-(max(emin0,etrans1)/etrans2)**brad2)/brad2*etrans2
      else
         gnorm2=max(gnorm2,w_EBL_density(etrans1,zz) &
         *sigics(e0,ame**2+etrans1*de)*etrans1*(etrans2/etrans1)**brad2)
         aw2=gnorm2*(1.d0-(max(emin0,etrans1)/etrans2)**(brad2+1.d0)) &
         /(brad2+1.d0)*etrans2
      endif
   else
      aw2=0.d0
   endif

   if(emin0.lt.etrans3)then
      gnorm3=etrans3*w_EBL_density(etrans3,zz)*sigics(e0,ame**2+etrans3*de)
      if(etrans3*de.gt.10.d0*ame**2)then
         gnorm3=max(gnorm3,w_EBL_density(etrans2,zz)*sigics(e0,ame**2 &
         +etrans2*de)*etrans2*exp(-(etrans3-etrans2)/akt/(1.d0+zz)))
         aw3=gnorm3*(exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))-1.d0) &
         *akt*(1.d0+zz)
      else
         gnorm3=max(gnorm3,w_EBL_density(max(2.d0*emin0,etrans2),zz) &
         *sigics(e0,ame**2+max(2.d0*emin0,etrans2)*de)*etrans3 &
         *exp(-(etrans3-max(2.d0*emin0,etrans2))/akt/(1.d0+zz)))
         aw3=gnorm3*akt*(1.d0+zz)*(max(emin0,etrans2)/etrans3 &
         *exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))-1.d0 &
         +akt*(1.d0+zz)/etrans3 &
         *(exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))-1.d0))
      endif
   else
      aw3=0.d0
   endif

   gnorm4=eirmax*w_EBL_density(eirmax/2.d0,zz)*sigics(e0,ame**2+eirmax*de)
   if(etrans3*de.gt.10.d0*ame**2)then
      gnorm4=max(gnorm4,w_EBL_density(etrans3,zz) &
      *sigics(e0,ame**2+etrans3*de)*etrans3*(etrans3/eirmax)**3)
      aw4=gnorm4*((eirmax/max(emin0,etrans3))**2-1.d0)/2.d0*eirmax
   else
      gnorm4=max(gnorm4,w_EBL_density(etrans3,zz) &
      *sigics(e0,ame**2+etrans3*de)*etrans3*(etrans3/eirmax)**2)
      aw4=gnorm4*(eirmax/max(emin0,etrans3)-1.d0)*eirmax
   endif

   1 rrr=psran()*(aw1+aw2+aw3+aw4)

   ! sample emin (= (sgam-m_e^2)/de) according to the 'proposal function';
   ! define 'rejection function' ( gb = f(emin) / f_proposal(emin) )
   if(rrr.lt.aw1)then
      emin=etrans1*(1.d0-psran()*(1.d0-(emin0/etrans1)**(brad1+1.d0))) &
      **(1.d0/(brad1+1.d0))
      gb0=gnorm1*(emin/emin0)**brad1
   elseif(rrr.lt.aw1+aw2)then
      if(etrans2*de.gt.10.d0*ame**2)then
         emin=etrans2*(1.d0-psran() &
         *(1.d0-(max(emin0,etrans1)/etrans2)**brad2))**(1.d0/brad2)
         gb0=gnorm2*(emin/etrans2)**(brad2-1.d0)
      else
         emin=etrans2*(1.d0-psran()*(1.d0-(max(emin0,etrans1)/etrans2) &
         **(1.d0+brad2)))**(1.d0/(1.d0+brad2))
         gb0=gnorm2*(emin/etrans2)**brad2
      endif
   elseif(rrr.lt.aw1+aw2+aw3)then
      2  emin=etrans3-akt*(1.d0+zz)*dlog(1.d0-psran() &
      *(1.d0-exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))))
      gb0=gnorm3*exp((etrans3-emin)/akt/(1.d0+zz))

      if(etrans3*de.lt.10.d0*ame**2)then
         if(psran().gt.emin/etrans3)goto 2
         gb0=gb0*emin/etrans3
      endif
   else
      if(etrans3*de.gt.10.d0*ame**2)then
         emin=eirmax/dsqrt(1.d0-psran()*(1.d0-(eirmax/max(emin0,etrans3))**2))
         gb0=gnorm4*(eirmax/emin)**3
      else
         emin=eirmax/(1.d0-psran()*(1.d0-eirmax/max(emin0,etrans3)))
         gb0=gnorm4*(eirmax/emin)**2
      endif
   endif

   gb0=gb0*gbnorm

   if(model.eq.1)then
      if(emin0/etrans3.lt..3d0)then
         gb0=gb0*1.2d0
      elseif(emin0/etrans3.gt..3d0.and.emin0/etrans3.lt..55d0)then
         gb0=gb0*20.d0
      endif
   elseif(model.eq.4)then
      if(emin0/etrans3.lt.1.d-2)gb0=gb0*2.d0
   elseif(model.eq.3)then
      gb0=gb0*1.1d0
      if(emin0/etrans3.gt.15.d0)then
         gb0=gb0*220.
      elseif(emin0/etrans3.gt..8d0)then
         gb0=gb0*60.
      endif
   endif

   if(.not.(gb0.gt.0.d0.and.gb0.lt.1.d60))then
      write(*,*)'electron: gb0=0',gb0,gbnorm
      stop
      ierr=1
      return
   endif

   sgam=ame**2+emin*de                 ! c.m. energy for e(+-)-gamma interaction
   gb=w_EBL_density(emin,zz)*sigics(e0,sgam)*emin/gb0

   !if (gb.gt.1.d0.and.nrenorm.eq.0) write(*,*)'sample_cmb(electron): gb='  &
   !   ,gb,nrej,nrenorm,emin0/etrans3

   if (psran().gt.gb) then       ! rejection
      nrej=nrej+1                ! total number of rejections for current sampling
      gbmax=max(gbmax,gb)        ! maximal value for rejection function
      if(nrej.gt.nrejmax)then    ! too many rejections
         if(gbmax.le.0.d0)then     ! wrong kinematics
            write(*,*)'electron: gbmax=0!!!'
            ierr=1
            return
         else
   !       write(*,*)'nrej(e)>nrejmax'  &
   !       ,nrej,emin0/etrans3,nrenorm,gbmax,e0/1.d9,emin0/etrans2  !,gb,gbnorm,gb0
            gbnorm=gbnorm*gbmax*2.d0 ! change norm for the rejection function
            gbmax=0.d0
            nrenorm=nrenorm+1
            if(nrenorm.gt.100)stop
            nrej=0
         endif
      endif
      goto 1                     ! new try
   end if

end subroutine sample_electron
!=============================================================================!
!=============================================================================!
!     weighted background photon density Eq. (9)    (interpolation)           ! 
!=============================================================================!
double precision function w_EBL_density(emin,zz)
!-----------------------------------------------------------------------------
! input:
!        emin - minimal required energy for a background photon
!        zz - current redshift
!-----------------------------------------------------------------------------
   use EBL_fit
   use xsec
   implicit none
   integer k,k1,jz,l1
   double precision emin,zz,dz,wk(3),wz(3),yl

   w_EBL_density=0.d0
   if (emin.ge.eirmax) return

   if (emin.le.erad1) then
      yl=log(emin*1.d15)/dlog(erad1*1.d15)*10.d0+1.d0
   else
      yl=log(emin/erad1)/dlog(eirmax/erad1)*40.d0+11.d0
   endif
   k=min(int(yl),48)
   k=max(k,1)
   if (k.eq.10) k=9
   wk(2)=yl-k
   wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
   wk(1)=1.d0-wk(2)+wk(3)
   wk(2)=wk(2)-2.d0*wk(3)

   dz=10.d0*zz+1.d0
   jz=min(49,int(dz))
   wz(2)=dz-jz
   wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
   wz(1)=1.d0-wz(2)+wz(3)
   wz(2)=wz(2)-2.d0*wz(3)

   do k1=1,3
      do l1=1,3
         w_EBL_density=w_EBL_density+denc(k+k1-1,jz+l1-1)*wk(k1)*wz(l1)
      end do
   end do
   w_EBL_density=exp(w_EBL_density-(1.d0-brad1)*dlog(emin))

end function w_EBL_density
!=============================================================================!
!=============================================================================!
!  interaction length for interactions on EBL photons    (interpolation)      !
!=============================================================================!
double precision function int_length(e0,zz,icq)
!-----------------------------------------------------------------------------
! calls: rate_bb, psran
! input:
!        e0 - particle energy;
!        zz - particle redshift position;
!        icq - particle type (0 - photon, +/-1 - electron/positron)
!-----------------------------------------------------------------------------
   use internal, only : rcmb
   implicit none
   double precision e0,zz,sig
   double precision rate_bb,psran
   integer icq

   sig = rate_bb(e0,zz,icq)           ! inverse mean free pass
   if (sig.le.0.d0) then              ! particle energy below threshold
      int_length=1.1d0*rcmb           ! -> no interaction
   else
      int_length=-log(psran())/sig    ! interaction lengt 
   endif

end function int_length
!=============================================================================!
!=============================================================================!
!        pair production x-section                                            !
!=============================================================================!
double precision function sigpair(sgam)
!-----------------------------------------------------------------------------
! input:
!        sgam - c.m. energy for gamma-gamma interaction
!-----------------------------------------------------------------------------
   use constants, only : ame,sigtmp
   implicit none
   double precision sgam,bet

   bet=1.d0-4.d0*ame**2/sgam
   if (bet.le.0.d0) then
      sigpair=0.d0
   else
      bet=sqrt(max(0.d0,bet))
      sigpair=sigtmp*3.d0/4.d0*ame**2/sgam*((3.d0-bet**4) &
      & *log((1.d0+bet)/(1.d0-bet))-2.d0*bet*(2.d0-bet**2))
   end if

end function sigpair
!=============================================================================!
!=============================================================================!
!        energy partition for pair production                                 !
!=============================================================================!
double precision function zpair(sgam)
!-----------------------------------------------------------------------------
! calls: psran
! input:
!        sgam - c.m. energy for gamma-gamma interaction
!-----------------------------------------------------------------------------
   use constants, only : ame
   use internal, only : debug
   implicit none
   double precision sgam,bet,zmin,z,gb
   double precision psran

   bet=sqrt(max(0.d0,1.d0-4.d0*ame**2/sgam))
   zmin=(1.d0-bet)/2.d0
   do 
      z=.5d0*(2.d0*zmin)**psran()
      gb=(z**2/(1.d0-z)+1.d0-z+(1.d0-bet**2)/(1.d0-z)-(1.d0-bet**2)**2 &
      &/4.d0/z/(1.d0-z)**2)/(1.d0+2.d0*bet**2*(1.d0-bet**2))
      if (debug.gt.0.and.gb.gt.1.d0) write(*,*)'zpair: gb=',gb
      if (psran().lt.gb) exit
   end do
   if (psran().gt.0.5d0) z=1.d0-z
   zpair=z

end function zpair
!=============================================================================!
!=============================================================================!
!        (inverse) Compton x-section                                          !
!=============================================================================!
double precision function sigics(e0,sgam)
!-----------------------------------------------------------------------------
! input:
!        e0 - electron/positron energy;
!        sgam - c.m. energy for e(+-)-gamma interaction
!-----------------------------------------------------------------------------
   use constants, only : ame,sigtmp
   use user_variables, only : ethr
   implicit none
   double precision e0,sgam,zmax,zm,t

   zmax=1.d0-ethr/e0
   zm=ame**2/sgam
   zmax=max(zmax,zm)
   if (1.d0-zm.le.3.d-3) then              ! use series expansion near threshold
      t=min(1.d0,(zmax-zm)/(1.d0-zm))
      sigics=.75d0*sigtmp*t*(4.d0*t*(t*(.5d0/zmax-(1.d0+zm)/3.d0 &
      &  /zm**2)-min(1.d0,(1.d0-zmax)/(1.d0-zm))*zm/zmax/2.d0) &
      &  +(zmax+zm)*zm/2.d0+1.d0-(zmax-zm)/zm/2.d0 &
      &  +(zmax-zm)**2/zm**2/3.d0)
   else
      sigics=.75d0*sigtmp*zm*min(1.d0,(zmax-zm)/(1.d0-zm)) &
      &  *(dlog(zmax/zm)/(zmax-zm)*(1.d0-4.d0*zm*(1.d0+zm)/(1.d0-zm)**2) &
      &  +4.d0*(zm/zmax+zm)/(1.d0-zm)**2+(zmax+zm)/2.d0)
   end if
   if(sigics.le.0.d0)then
      ! write(*,*)'sigics:e0,sgam,zm,zmax',e0,sgam,zm,zmax,sigics
      sigics=0
   endif
end function sigics
!=============================================================================!
!=============================================================================!
!        energy partition for inverse Compton    (E-fraction taken by e+-)    !
!=============================================================================!
double precision function zics(e0,sgam)
!-----------------------------------------------------------------------------
! calls: psran
! input:
!        sgam - c.m. energy for e(+-)-gamma interaction
!-----------------------------------------------------------------------------
   use constants, only : ame
   use user_variables, only : ethr
   use internal, only : debug
   implicit none
   double precision e0,sgam,zmin,zmax,z,gb
   double precision psran

   zmax=1.d0-ethr/e0
   zmin=ame**2/sgam

   if (zmin.ge.zmax) then
      if (debug.gt.0) call error('zmin>zmax in zics',1)
      zics=zmin
   else
      do
         z=zmin*(zmax/zmin)**psran()
         gb=(1.d0+z*z)/2.d0-2.d0*zmin/z*(z-zmin)*(1.d0-z)/(1.d0-zmin)**2
         !if (debug>0.and.gb.gt.1.d0) write(*,*)'zics: gb=',gb,z,zmin,zmax
         if (psran().lt.gb) exit
      end do
      zics=z
   endif

end function zics
!=============================================================================!
!=============================================================================!
!        cross section * energy loss (under threshold) for ICS
!=============================================================================!
double precision function zsigics(e0,sgam) 
!-----------------------------------------------------------------------------
! calls: error
! input:
!        e0 - electron/positron energy;
!        sgam - c.m. energy for e(+-)-gamma interaction
!-----------------------------------------------------------------------------
   use constants, only : ame,sigtmp
   use user_variables, only : ethr
   implicit none
   double precision e0,sgam,zmax,zm

   zm=ame**2/sgam
   zmax=max(1.d0-ethr/e0,zm)

   if(zm.lt..5d0)then
      zsigics=.75d0*sigtmp*zm*(1.d0-zmax)/(1.d0-zm)                        &
      &  *((-dlog(zmax)/(1.d0-zmax)-1.d0)*(1.d0-4.d0*zm*(1.d0+2.d0*zm)     &
      &  /(1.d0-zm)**2)+(1.d0-zmax)*(1.d0+2.d0*zmax)/6.d0                  &
      &  +2.d0*zm*(1.d0-zmax)/(1.d0-zm)**2*(1.d0+2.d0*zm/zmax))
   else if (zm.eq.1.d0) then! added in v203 for limiting case  <--------
      zsigics = 3.d0*sigtmp ! added in v203 for limiting ccase <--------
   else
      zsigics=.75d0*sigtmp*zm*(1.d0-zmax)**2/(1.d0-zm)                     &
      &  *(1.d0+.25d0*(1.d0-zmax)**2+4.d0*zm*(1.d0-zmax)/(1.d0-zm)**2      &
      &  *(zm/zmax-(1.d0+2.d0*zm)/3.d0*(1.d0+.75d0*(1.d0-zmax))))
   endif


   if(zsigics.lt.0.d0)then
   call error('zsigics<0',1)
   zsigics=0.d0
   endif

end function zsigics
!=============================================================================!
!=============================================================================!
!        relative energy loss (per cm) (due to under-threshold ICS photons)   !
!=============================================================================!
double precision function zloss(e0,zz)
!-----------------------------------------------------------------------------
! input:
!        e0 - electron/positron energy;
!        zz - current redshift
!-----------------------------------------------------------------------------
   use constants
   use EBL_fit
   use user_variables, only : ethr,egmax
   use xsec
   implicit none
   integer jz,k,k1,l1
   double precision e0,zz,yl,dz,wk(3),wz(3)

   zloss=0.d0
   yl=dlog(e0/ame)/dlog(egmax/ame)*50.d0+1.d0
   k=min(int(yl),49)
   k=max(k,2)
   wk(2)=yl-k
   wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
   wk(1)=1.d0-wk(2)+wk(3)
   wk(2)=wk(2)-2.d0*wk(3)

   dz=10.d0*zz+1.d0
   jz=min(49,int(dz))
   wz(2)=dz-jz
   wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
   wz(1)=1.d0-wz(2)+wz(3)
   wz(2)=wz(2)-2.d0*wz(3)

   do k1=1,3
      do l1=1,3
         zloss=zloss+zsigc(k+k1-1,jz+l1-1)*wk(k1)*wz(l1)
      enddo
   enddo
   zloss=exp(zloss)

end function zloss
!=============================================================================!
!=============================================================================!
!       interpolates interaction rate/cm on photon background                 !
!=============================================================================!
double precision function rate_bb(e0,zz,icq)
!-----------------------------------------------------------------------------
! input:
!        e0 - particle energy;
!        zz - particle redshift position;
!        icq - particle type (0 - photon, +/-1 - electron/positron)
!-----------------------------------------------------------------------------
   use constants
   use EBL_fit
   use user_variables, only : ethr,egmax
   use xsec
   implicit none
   integer icq,ica,jz,k,k1,l1
   double precision e0,zz,emin,yl,dz,wk(3),wz(3)

   rate_bb=0.d0
   if (icq.eq.0) then
      emin=ame**2/eirmax
   else
      emin=ethr+.5d0*ame**2/eirmax/(1.d0+dsqrt(1.d0+ame**2/ethr**2/eirmax*(ethr-eirmax)))
   end if
   if (e0.le.emin) return

   ica=iabs(icq)+1
   yl=log(e0/emin)/log(egmax/emin)*50.d0+1.d0
   k=min(int(yl),49)
   k=max(k,2)
   wk(2)=yl-k
   wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
   wk(1)=1.d0-wk(2)+wk(3)
   wk(2)=wk(2)-2.d0*wk(3)

   dz=10.d0*zz+1.d0
   jz=min(49,int(dz))
   wz(2)=dz-jz
   wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
   wz(1)=1.d0-wz(2)+wz(3)
   wz(2)=wz(2)-2.d0*wz(3)

   do k1=1,3
      do l1=1,3
         rate_bb=rate_bb+sigc(k+k1-1,jz+l1-1,ica)*wk(k1)*wz(l1)
      enddo
   enddo
   rate_bb=exp(rate_bb)

end function rate_bb
!=============================================================================!
!=============================================================================!
!        electron synchrotron energy loss (eV/cm)                     !
!=============================================================================!
double precision function eloss(e0,begmf)
!-----------------------------------------------------------------------------
! input:
!        e0 - electron/positron energy;
!        begmf - strength of (transverse) extragalactic B-field
!-----------------------------------------------------------------------------
   use constants
   implicit none
   double precision e0,begmf,chi

   chi=sqrt(max(0.d0,(e0/ame)**2-1.d0))*begmf/bcr
   eloss = chi**2/(1.d0+4.8d0*(1.d0+chi)*log(1.d0+1.7d0*chi) +  &
   3.44d0*chi**2)**(2.d0/3.d0)*ame**2/137.d0/1.5d0*1.d7/197.d0 
end function eloss
!=============================================================================!
!=============================================================================!
!        electron/positron deflection angle (rad/cm) by EGMF                  !
!=============================================================================!
double precision function themf(e0,begmf)
!-----------------------------------------------------------------------------
! input:
!        e0    - electron/positron energy;
!        begmf - strength of (transverse) extragalactic B-field;
!-----------------------------------------------------------------------------
  implicit none
  double precision e0,begmf
  themf=294.d0*begmf/e0
end function themf
!=============================================================================!
!=============================================================================!