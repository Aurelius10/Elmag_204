program Elmagcasc
!-----------------------------------------------------------------------------
! calls: init,user_main,user_output,banner
!-----------------------------------------------------------------------------
  use user_variables, only : nmax,old_calc
  use user_result
  implicit none
  integer i
  integer myid,n_proc

! non-MPI values
  myid = 0
  n_proc = 1 

  call init(myid,n_proc)
  if (old_calc.eq.0) call init_turb(myid)
  call user_main(myid,nmax)

! non-MPI values
  spec_tot = spec
  gam_th_tot = gam_th
  spec_t_tot = spec_t
  anobs_tot = anobs
  
  do i=1,n_bint
    write(*,*) n_reg(i), 'detectable photons in time bin', i
  end do
  write(*,*) n_kc, 'bad particles'

  if (myid==0) then
     call user_output(nmax,n_proc)  ! make plots
     call banner(n_proc,1)
  end if
  close(99)               

end program Elmagcasc