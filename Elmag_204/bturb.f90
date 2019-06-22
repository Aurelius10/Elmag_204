!=============================================================================!
!=============================================================================!
module random_phases
  use bfield_parameters
  implicit none
  save
  double precision ca(n_k),sa(n_k),cp(n_k),b(n_k),sp(n_k),ct(n_k), &
&   st(n_k),s(n_k)
end module random_phases
!=============================================================================!
!=============================================================================!
subroutine init_turb(myid)
  use constants, only : two_pi
  use bfield_parameters
  implicit none
  integer myid
  double precision L_min,L_max,z

  if (model_b.eq.1) return
  call init_phases
  call init_norm

  L_min = two_pi/k_max
  L_max = two_pi/k_min

  if (L_min==L_max) then
    L_c = L_max/2.d0
  elseif (gamma.lt.1.01d0) then
    z = L_min/L_max
    L_c = L_max/2.d0*(1-z)/log(1/z)
  else
    z = L_min/L_max
    L_c = L_max/2.d0*(gamma-1.d0)/gamma*(1.d0-z**gamma)/(1.d0-z**(gamma-1.d0))
  end if
  if (myid.eq.0) then
    write(*,*) "!  L_coh = ", L_c
    write(*,*) "!  B_rms = ", B_rms
  end if
end subroutine init_turb

!=============================================================================!
!=============================================================================!

subroutine generate_B(x, Omega)
  use constants
  use bfield_parameters
  use random_phases
  implicit none
  integer i
  double precision k,l,O_k,d_O(3),Omega(3),x(3)
  complex*16 v(3),ar
  complex*16, parameter :: im = (0.d0,1.d0)

  if (model_b.eq.1) then
    Omega(1) = B_rms
    Omega(2) = 0.d0
    Omega(3) = 0.d0
    return
  end if

  l = (k_max/k_min)**(1.d0/n_k) 
  d_O = 0.d0
  k = k_min
  
  do i=1,n_k

     O_k = O_min * (k/k_min)**(-gamma/2.d0)
     v(1) = -ca(i)*sp(i) - im*s(i)*sa(i)*st(i)*cp(i)
     v(2) =  ca(i)*cp(i) - im*s(i)*sa(i)*st(i)*sp(i)
     v(3) =                im*s(i)*sa(i)*ct(i)
     ar = k*( ct(i)*(cp(i)*x(1)+sp(i)*x(2))+st(i)*x(3) ) + b(i) 
     d_O = d_O + O_k * real( v*exp(im*ar) )

     k=l*k

  end do
  Omega = d_O    
      
end subroutine generate_B
!=============================================================================!
!=============================================================================!
subroutine init_phases
  use random_phases
  use constants
  use bfield_parameters, only : h
  implicit none
  integer i
  double precision r,cc,ss
  double precision ran0

  do i=1,n_k
     
     call cs(cc,ss) 
     ca(i) = cc            
     sa(i) = ss            
          
     call cs(cc,ss) 
     cp(i) = cc                           ! 2-dim: p(i) = 0.d0
     sp(i) = ss
     b(i) = two_pi * ran0()! 0.2!ran0()
     r = pi * ran0()! 0.7!ran0()
     ct(i) = cos(r)
     st(i) = sin(r)
     s(i)=-1
     if (ran0().lt.h) s(i)=1

  enddo

end subroutine init_phases
!=============================================================================!
!=============================================================================!
subroutine cs(cc,ss)
  implicit none
  double precision cc,ss,s1,s2,s3
  double precision ran0
  
  do
     s1=2.d0*ran0()-1.d0 !2.d0*0.82-1.d0
     s2=2.d0*ran0()-1.d0 !2.d0*0.41-1.d0
     s3=s1*s1+s2*s2
     if(s3<=1.d0) exit
  enddo
  s3=dsqrt(s3)
  cc=s1/s3
  ss=s2/s3      
  
end subroutine cs
!=============================================================================!
!=============================================================================!
subroutine init_norm
  use constants
  use bfield_parameters
  implicit none
  integer i
  double precision k,N,v

  N = 0.d0
  v = (k_max/k_min)**(1.d0/n_k)
  k = k_min
  do i=1,n_k
     N = N+(k/k_min)**(-gamma)   
     k = v*k
  enddo
  O_min = B_rms/sqrt(N)

end subroutine init_norm

!=============================================================================!
!=============================================================================!
