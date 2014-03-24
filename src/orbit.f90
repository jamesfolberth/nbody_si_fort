! orbit.f90

program orbit
   
   use utils ! helper routines, data types
   use orbital_data
   use coords 
   use kepler
   use symplectic

   implicit none


   ! run config parameters
   character (len=256), parameter :: savefile='../data/orbit.h5'
   integer (kind=intk), parameter :: N_record_int = 1000 ! record state every N_record_int time steps
   integer (kind=intk), parameter :: N_saves = 20
   
   real (kind=dblk), parameter :: t0 = 1941.+6./365.25 ! JD=2430000.5
   real (kind=dblk), parameter :: t1 = t0+10**5 ! JD=2430000.5
   real (kind=dblk), parameter :: dt = 1.0_dblk ! time step
   integer (kind=intk), parameter :: N_records = ceiling((t1-t0)/N_record_int)
   integer (kind=intk), parameter :: N_save_int = ceiling(dble(N_records)/dble(N_saves))
   
   ! variables
   real (kind=dblk), allocatable :: t(:)
   real (kind=dblk) :: P_wrk(3*n_masses),Q_wrk(3*n_masses),&
      Pjac_wrk(3*n_masses),Qjac_wrk(3*n_masses)
   real (kind=dblk), allocatable :: P(:,:),Q(:,:),Pjac(:,:),Qjac(:,:)
   real (kind=dblk) :: jacQ(3*n_masses,3*n_masses),&
      LUjacQ(3*n_masses,3*n_masses),jacP(3*n_masses,3*n_masses),&
      LUjacP(3*n_masses,3*n_masses),jacT(3*n_masses,3*n_masses),&
      eta(n_masses),m_vec_jac(n_masses), g_param_jac(n_masses)
   integer (kind=intk):: PjacQ(3*n_masses), PjacP(3*n_masses)

   ! variables 'local' to kepler_step
   real (kind=dblk) :: kep_r(n_masses), kep_v(n_masses), kep_u(n_masses),&
      kep_a(n_masses), kep_n(n_masses), kep_EC(n_masses), kep_ES(n_masses),&
      kep_e(n_masses), kep_dE(n_masses), kep_dtv(n_masses),kep_C(n_masses),&
      kep_S(n_masses),kep_f(n_masses), kep_g(n_masses), kep_aor(n_masses),&
      kep_fp(n_masses),kep_gp(n_masses)

   ! variables 'local' to symp_step
   real (kind=dblk) :: symp_interdv(3*n_masses),symp_interdvjac(3*n_masses),&
      symp_Q_wrk(3*n_masses),symp_qimqj(3*n_masses,3*n_masses),&
      symp_qimqjnrm(n_masses,n_masses),symp_ind_wrk1(3*n_masses),&
      symp_ind_wrk2(3*n_masses)

   integer (kind=intk) :: i,j,save_counter,clock
   real (kind=dblk) :: ti


   ! initialize (allocatable) variables
   allocate(t(N_records))
   allocate(P(3*n_masses, N_records), Q(3*n_masses, N_records))
   allocate(Pjac(3*n_masses, N_records), Qjac(3*n_masses, N_records))

   P=0; Q=0; Pjac=0; Qjac=0;


   ! Set the initial conditions
   call WH_initial_data(P,Q)
   !call print_vector(Q(:,1))

   ! set up jacobi coordinate functions
   call jacobi_setup(jacQ,jacP,jact,PjacQ,LUjacQ,PjacP,LUjacP,eta,m_vec_jac,g_param_jac)

   Q_wrk = Q(:,1); P_wrk = P(:,1);
   call apply_jacobi(Q(:,1),P(:,1),Qjac_wrk,Pjac_wrk,jacQ,jacP)
   Qjac(:,1) = Qjac_wrk; Pjac(:,1) = Pjac_wrk;

   ! Save data initially (useful for debugging stuff)
   if (debug) then
      call save_orbit(savefile,t,Q,P,Qjac,Pjac,jacQ,jacP,jact,PjacQ,LUjacQ,PjacP,LUjacP,m_vec,m_vec_jac,g_const,g_param,g_param_jac)
   end if

  
   ! main loop
   write (*,"(A,10ES13.4)"), "Num iterations (approx) = ",(t1-t0)/dt

   i = 0
   save_counter = 0
   ti = t0
   
   ! Do half a step of kepler
   ! TODO should ti also be moved forward by 0.5dt?
   call kepler_step(Qjac_wrk, Pjac_wrk, 0.5_dblk*dt, kep_r,kep_v,kep_u,kep_a,&
      kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,kep_C,kep_S,&
      kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)


   call tic(clock)
   main: do while (ti < t1)

      integrate: do j=1,N_record_int
        
         ! second-order SI method (with 0.5*dt from above and below)
         ! Do a full step of SI
         call symp_step(Qjac_wrk,Pjac_wrk,dt,symp_interdv,symp_interdvjac,&
            symp_Q_wrk,symp_qimqj,symp_qimqjnrm,jacQ,PjacQ,LUjacQ,jacT,&
            symp_ind_wrk1,symp_ind_wrk2,m_vec,m_vec_jac,g_const)

         ! Do full step of kepler
         call kepler_step(Qjac_wrk, Pjac_wrk, dt, kep_r,kep_v,kep_u,kep_a,&
            kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,kep_C,kep_S,&
            kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)
         
         ti = ti + dt;
      end do integrate

      ! Record state of system
      i=i+1; 
      t(i)=ti;
      Qjac(:,i) = Qjac_wrk; Pjac(:,i) = Pjac_wrk;
      call apply_jacobi_inv(Qjac_wrk,Pjac_wrk,Q(:,i),P(:,i),&
         PjacQ,LUjacQ,PjacP,LUjacP)


      if (((t1-t0)/dt > 10**7) .and. &
         (mod(i,floor(dble(N_save_int)/100_dblk)) == 0)) then
         print *, "Percent complete: ", ti/t1*100_dblk
         
         ! save data every N_save_int iterations of i
         if (mod(i,N_save_int) == 0) then
            call save_orbit(savefile,t,Q,P,Qjac,Pjac,jacQ,jacP,jact,&
               PjacQ,LUjacQ,PjacP,LUjacP,m_vec,m_vec_jac,g_const,&
               g_param,g_param_jac)

         end if
      end if

   end do main
   call toc(clock)

   ! Finish the integration with a half-step (0.5dt) of kepler
   call kepler_step(Qjac_wrk, Pjac_wrk, 0.5_dblk*dt, kep_r,kep_v,kep_u,kep_a,&
      kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,kep_C,kep_S,&
      kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)

   Qjac(:,i) = Qjac_wrk; Pjac(:,i) = Pjac_wrk;
   call apply_jacobi_inv(Qjac_wrk,Pjac_wrk,Q(:,i),P(:,i),&
      PjacQ,LUjacQ,PjacP,LUjacP)
 
   call save_orbit(savefile,t,Q,P,Qjac,Pjac,jacQ,jacP,jact,&
      PjacQ,LUjacQ,PjacP,LUjacP,m_vec,m_vec_jac,g_const,&
      g_param,g_param_jac)

   print *, "Computation complete!"
end program orbit
