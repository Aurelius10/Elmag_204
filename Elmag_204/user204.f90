!=============================================================================!
!=============================================================================!
module bfield_parameters
   implicit none
   save
   integer, parameter :: model_b = 2           ! 1 uniform, 2 turbulent
   integer, parameter :: n_k=100               ! Number of modes
   double precision, parameter ::           &
      B_rms = 1.d-14,                       &  ! in Gauss
      gamma = 5.0d0/3.d0,                   &  ! power-law fluctuations
      k_min = 1/(3.086d18*.75d6),            &
      k_max = 1/(3.086d18*.75d4),            &
      h = 0.5d0               ! helicity of bfield 1: right handed, 0: left handed
   double precision :: L_c = 3.086d18*1.d6     ! 1Mpc
   double precision O_min
end module bfield_parameters
!=============================================================================!
!=============================================================================!
module user_variables
   use constants, only : degree_rad
   implicit none
   save
   integer,parameter :: old_calc = 0   ! 0: calculate high energy (>e10 eV)
                                       ! 1: calculate low energy  (<e10 eV)
   integer,parameter :: model = 1      ! EBL model: 1 best fit, 2 lower limit,
                                       ! 3: Franceschini arXiv:0805.1841
                                       ! 4: Model C Finke et al. arXiv:0905.1115
                                       ! 5: Gilmore et al. arXiv:1104.0671
                                       ! 6: Dominguez et al. arXiv:1007.1459
   integer,parameter :: nmax=100000   ! # injected photons
   
   double precision,parameter ::    &
      ethr=1.d6,                    & ! energy threshold, eV
      egmax=2.d13,                  & ! maximal gamma energy, eV
      ir_rescale=1.d0,              & ! multiplicative rescaling IR background
      th_jet=6.d0,                  & ! jet opening angle/degrees
      th_jetx=6.d0,                 & ! jet misalignment angle towards earth
      a_smp=0.0d0,                  & ! 1 max. weighted sampling, 0-no sampling
      vx0=sin(th_jetx*degree_rad),  &
      vz0=cos(th_jetx*degree_rad)
end module user_variables
!=============================================================================!
!=============================================================================!
module user_result
   implicit none
   save
   integer, parameter ::            &
      n_bin = 81,                   &  ! # of bins in agamma
      n_bin2 = 2,                   &  ! # angular bin
      n_bin3 = 7,                   &  ! # time bins old
      n_bint = 5,                   &  ! # time bins new
      n_binx = 541,                 &  ! # x-dir resolution
      n_biny = 421,                 &  ! # y-dir resolution
      n_bind = 30,                  &  ! # bins per degree
      shiftx = 3                       ! # telling plotter to shift the x-axs
   double precision spec(n_bin,0:1),spec_tot(n_bin,0:1)     ! diffuse spectrum
   double precision gam_th(n_bin,n_bin2),gam_th_tot(n_bin,n_bin2) ! spectrum
   double precision spec_t(n_bin,n_bin3),spec_t_tot(n_bin,n_bin3) ! spectrum
   double precision cum_t(n_bin,n_bin3)                     ! cum. spectrum
   double precision anobs(n_bint,n_binx,n_biny),anobs_tot(n_bint,n_binx,n_biny)
   integer     :: n_reg(n_bint)=0      ! Counting particles in each time bin
   integer     :: n_kc=0               ! Counting bad particles
   character*5 :: filename = '_mpi'         ! name in output
end module user_result
!=============================================================================!
!=============================================================================!
subroutine user_main(myid,nmax)
!-----------------------------------------------------------------------------!
!           calls: plots,initial_particle,cascade                             !
!-----------------------------------------------------------------------------!
   use bfield_parameters, only : h
   use user_result
   implicit none
   integer myid,nl,icq,nmax
   double precision z,e0,weight

   if (myid==0) then
      call plots
      write(*,*)'start with nmax = ',nmax
   end if
   spec=0.d0
   gam_th=0.d0
   spec_t=0.d0
   cum_t=0.d0
   anobs=0.d0

   z = 0.14d0                                      ! initial redshift
   do nl=1,nmax
      call initial_particle(e0,weight)             ! generate initial energy
      icq = 0                                      ! (0 - gamma, +-1 - e+-)
      call cascade(icq,e0,weight,z)                ! starts EM-cascade
      if (myid.eq.0.and.mod(nl*100,nmax).eq.0) write(*,*)'nl=',nl
   enddo
end subroutine user_main
!=============================================================================!
!        sample initial photon/electron from the input spectra                !
!=============================================================================!
subroutine initial_particle(e0,weight)
!-----------------------------------------------------------------------------
! output:
!        e0       - particle energy;
!        weight   - initial particle weight;
!-----------------------------------------------------------------------------
   use user_variables
   implicit none
   double precision e0,weight,emin,gam,ebreak
   double precision psran

   !gam=-2.d0
   gam=-0.67d0
   emin=1.d12
   ebreak=2.d12
   e0=emin*(egmax/emin)**psran()                      !energy: uniform in ln E 
   weight=(e0/ebreak)**(gam+1.d0)*log(egmax/emin)

end subroutine initial_particle
!=============================================================================!
!        Additional weighting according to a jet distribution                 !
!=============================================================================!
subroutine jet_distribution(the_s,w_jet)
!-----------------------------------------------------------------------------!
!        input    :  the_s    -  angle from jet center
!        output   :  w_jet    -  weigth given jet distribution
!-----------------------------------------------------------------------------!
   use user_variables, only : th_jet
   implicit none
   double precision the_s,w_jet
   ! radial gaussian distribution with weight 1 at the center and
   w_jet=exp(-the_s**2/th_jet**2)  ! gaussian distribution
   !w_jet=1                       ! uniform
end subroutine jet_distribution
!=============================================================================!
!=============================================================================!
subroutine register(e0,thex,they,weight,dt,icq)
   use user_variables
   use user_result
   use constants, only : pi,degree_rad
   implicit none
   integer icq,i,j,k
   double precision e0,weight,thex,they,dt
   double precision thereg_en,theta,theta95
  
   ! diffuse energy spectrum:
   i=min(n_bin,int(log(e0/ethr)/log(egmax/ethr)*(n_bin-1))+1)
   i=max(i,1)
   spec(i,abs(icq))=spec(i,abs(icq))+weight*e0/log(egmax/ethr)*(n_bin-1)

   if (icq.ne.0) return                                ! forget electrons
   if (old_calc.eq.0) call psf_spread(e0,thex,they,weight,dt)

   theta=sqrt(thex**2+they**2)
   theta95 = thereg_en(e0)
   if (theta<theta95) then
      j=1
   else
      j=2
   end if
   gam_th(i,j)=gam_th(i,j)+weight*e0/log(egmax/ethr)*(n_bin-1)
   if (j.eq.1) then                                ! only gamma's inside PSF
      if (theta.eq.0.d0.or.dt.eq.0) then
         k=1
      else
         k = int(log10(dt))+1
      end if
      if (k<1) k=1
      if (k>n_bin3) k=n_bin3
      spec_t(i,k)=spec_t(i,k)+weight*e0/log(egmax/ethr)*(n_bin-1)
   end if

end subroutine register
!=============================================================================!
!        distributing detected particle according to a PSF function           !
!=============================================================================!
subroutine psf_spread(e0,thex,they,weight,dt)
   use user_result
   use constants, only : two_pi,rad_degree
   implicit none
   integer i,j,k
   double precision e0,thex,they,weight,dt
   double precision rx,ry,r,sig,thereg_en

   sig=thereg_en(e0)
   if (dt.lt.1.d-4) then
      k=5                  ! No-interacting particles stored in own array
   elseif (dt.lt.1.d5) then
      k=1
   elseif (dt.lt.1.d6) then
      k=2
   elseif (dt.lt.3.d6) then
      k=3
   elseif (dt.lt.1.d7) then
      k=4
   else
      return
   end if
   n_reg(k) = n_reg(k) + 1

   do i=1,n_binx
      rx=(i-n_binx/2.d0)/real(n_bind) + shiftx
      do j=1,n_biny
         ry=(j-n_biny/2.d0)/real(n_bind)
         r=sqrt((rx-thex)**2 + (ry-they)**2)
         ! homogeneous distribution of PSF
         if (r.lt.sig) anobs(k,i,j) = anobs(k,i,j) + weight*e0/sig**2
      end do
   end do
end subroutine psf_spread
!=============================================================================!
!        theta_PSF/degree (95%) from Fermi and HESS                           !
!=============================================================================!
double precision function thereg_en(en)
   implicit none
   double precision en

   if (en.lt.1.d6) then                            ! avoid theta>360 degrees 
      thereg_en = 359.d0
      return
   end if

   if (en.le.3e11) then                                             ! Fermi
      thereg_en = 10.d0**7.12d0*en**(-0.766d0)             
      if (en>1d9) thereg_en = thereg_en + 0.2d0*exp(-1.d10/en)
   else                                                             ! HESS
      thereg_en = 0.11d0                                      
   endif
end function thereg_en
!=============================================================================!
!=============================================================================!
subroutine user_output(n_max,n_proc) 
   use user_variables
   use user_result
   implicit none
   integer n_max,n_proc,i,j
   character*5 str

   spec_tot = spec_tot/dble(n_proc*n_max)
   gam_th_tot = gam_th_tot/dble(n_proc*n_max)
   spec_t_tot = spec_t_tot/dble(n_proc*n_max)
   anobs_tot = anobs_tot/dble(n_proc*n_max)
   
   cum_t(:,1) = spec_t_tot(:,1)
   do i=2,n_bin3
      cum_t(:,i) = cum_t(:,i-1)+spec_t_tot(:,i)
   end do

   open(unit=31,file='Output/spec_diff'//filename) ! energy spectrum g,e
   open(unit=32,file='Output/spec_95'//filename)   ! energy spectrum g in PSF
   open(unit=33,file='Output/spec_t'//filename)    ! energy-time spectra g in PSF
   open(unit=34,file='Output/spec_c'//filename)    ! energy-time spectra g in PSF
   do i=1,n_bin-1
      write(31,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
                  real(spec_tot(i,0)),real(spec_tot(i,1))
      write(32,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
                  real(gam_th_tot(i,1)),real(gam_th_tot(i,2))
      write(33,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
                  real(spec_t_tot(i,:))
      write(34,*) real(ethr*(egmax/ethr)**((i-.5d0)/(n_bin-1))), &
                  real(cum_t(i,:))
   end do
   close(31)
   close(32)
   close(33)
   close(34)

   do j=1,n_bint
      write (str,'(I1.1)') j
      open(unit=35,file='Output/AngRes/angle_matrix'//str)
      write(35,*) th_jetx, n_binx, n_biny, shiftx, n_bind
      do i=1,n_binx
         write(35,*) real(anobs_tot(j,i,:))
      end do
      close(35)
   end do
   write(*,*) "jet off by", th_jetx, "degrees"

end subroutine user_output
!=============================================================================!
!=============================================================================!
subroutine plots
   implicit none
   integer icq
   double precision E,x(3),z0,rate(2),psf95
   double precision rate_EBL_tab,thereg_en

   x(1) = 125.*3.086d24
   x(2) = 808.*3.086d24
   x(3) = 1500.*3.086d24

   z0 = 0.d0
   icq = 0 

   E = 1.d8
   open(41,file='Output/lint')
   open(42,file='Output/tau')
   open(43,file='Output/psf95')
   do
      E=E*1.1d0
      rate(1) = rate_EBL_tab(E,z0,icq)
      rate(2) = rate_EBL_tab(E,z0,1)
      write(41,*) real(E),rate(1)*3.086d24,rate(2)*3.086d24
      psf95 = thereg_en(E)
      write(43,*) real(E),psf95
      if (E>1d21) exit
   end do
   close(41)
   close(42)
   close(43)
end subroutine plots
