program Elmagcasc
   use mpi
   use user_variables, only : nmax,old_calc
   use user_result
   use random_phases
   use bfield_parameters, only : n_k
   implicit none
   integer i
   integer ierr,n_procs,myid,n_array,n_bad,n_ok(n_bint)
   n_ok=0
   n_bad=0
   
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,n_procs,ierr)
   
   call init(myid,n_procs)
   if (old_calc.eq.0) then
      call init_turb(myid)
   end if
   
   ! Same B-field on all processors
   CALL MPI_BCAST(ca, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(sa, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(cp, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(b, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(sp, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(ct, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(st, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST(s, n_k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

   call user_main(myid,nmax)
   
   n_array = 2*n_bin
   call MPI_REDUCE(spec,spec_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr) ! sum individal arrays spec
   call MPI_REDUCE(gam_th,gam_th_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr) ! sum individal arrays gam_th
   n_array = n_bin*n_bin3
   call MPI_REDUCE(spec_t,spec_t_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr) ! sum individal arrays spec
   
   n_array = n_bint*n_binx*n_biny
   call MPI_REDUCE(anobs,anobs_tot,n_array,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(n_kc,n_bad,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call MPI_REDUCE(n_reg,n_ok,n_bint,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   if (myid==0) then
      call user_output(nmax,n_procs) ! make plots
      do i=1,n_bint
         write(*,*) n_ok(i), 'detectable photons in time bin', i
      end do
      call banner(n_procs,1)
   end if
   
   close(99) 
   
   call MPI_Finalize(ierr) 
   
 end program Elmagcasc