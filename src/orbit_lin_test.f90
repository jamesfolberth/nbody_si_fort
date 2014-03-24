! orbit_lin_test.f90
! Test the linearization by integrating the test trajectory and reference trajectory forward one time step

program orbit_lin_test
   
   use utils ! helper routines, data types
   use orbital_data
   use coords 
   use kepler
   use symplectic

   implicit none

   ! Interface for F77 DLSODE
   interface
      subroutine dlsode(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,&
            lrw,iwork,liw,jac,mf)
         implicit none
         external :: f,jac
         integer, intent(in) :: itol,itask,iopt,lrw,liw,mf
         integer, intent(inout) :: istate
         integer, intent(in),dimension(*) :: neq
         integer, intent(inout) :: iwork(liw)
         double precision, intent(in) :: tout
         double precision, intent(inout) :: t
         double precision, intent(inout),dimension(*) :: y,rwork,rtol,atol
      end subroutine dlsode
   end interface


   ! run config parameters
   character (len=256), parameter :: savefile='../data/orbit_lin_test.h5'
   integer (kind=intk), parameter :: N_record_int = 100 ! record state every N_record_int time steps
   integer (kind=intk), parameter :: N_saves = 100
   
   real (kind=dblk), parameter :: t0 = 1941.+6./365.25 ! JD=2430000.5
   real (kind=dblk), parameter :: t1 = t0+10**5 ! JD=2430000.5
   !real (kind=dblk), parameter :: dt = 1.0_dblk ! time step
   real (kind=dblk), parameter :: dt = 0.1_dblk ! time step
   real (kind=dblk), parameter :: traj_pert = 1.0E-4_dblk
   integer (kind=intk), parameter :: N_records = ceiling((t1-t0)/N_record_int)
   integer (kind=intk), parameter :: N_save_int = ceiling(dble(N_records)/dble(N_saves))
   integer (kind=intk) :: N_iterations


   ! variables
   real (kind=dblk), allocatable :: t(:)

      ! reference trajectory
   real (kind=dblk) :: P_wrk(3*n_masses),Q_wrk(3*n_masses),&
      Pjac_wrk(3*n_masses),Qjac_wrk(3*n_masses)
   real (kind=dblk), allocatable :: P(:,:),Q(:,:),Pjac(:,:),Qjac(:,:)
      ! test trajectory
   real (kind=dblk) :: P_tst_wrk(3*n_masses),Q_tst_wrk(3*n_masses),&
      Pjac_tst_wrk(3*n_masses),Qjac_tst_wrk(3*n_masses)
   real (kind=dblk), allocatable :: P_tst(:,:),Q_tst(:,:),Pjac_tst(:,:),Qjac_tst(:,:)


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

   integer (kind=intk) :: i,j,clock,iter
   real (kind=dblk) :: ti

   ! Phase space distances
   real (kind=dblk), allocatable :: ps_dist_rat(:)
   real (kind=dblk) :: ps_dist_rat_avg,dxnrm
   real (kind=dblk) :: dPjac(3*n_masses),dQjac(3*n_masses)
   real (kind=dblk) :: dP(3*n_masses),dQ(3*n_masses)

   ! initialize (allocatable) variables
   allocate(t(N_records),ps_dist_rat(N_records))
   allocate(P(3*n_masses, N_records), Q(3*n_masses, N_records))
   allocate(Pjac(3*n_masses, N_records), Qjac(3*n_masses, N_records))
   allocate(P_tst(3*n_masses, N_records), Q_tst(3*n_masses, N_records))
   allocate(Pjac_tst(3*n_masses, N_records), Qjac_tst(3*n_masses, N_records))
   !allocate(P(3*n_masses,1), Q(3*n_masses,1))
   !allocate(Pjac(3*n_masses,1), Qjac(3*n_masses,1))
   !allocate(P_tst(3*n_masses,1), Q_tst(3*n_masses,1))
   !allocate(Pjac_tst(3*n_masses,1), Qjac_tst(3*n_masses,1))


   P=0; Q=0; Pjac=0; Qjac=0;
   P_tst=0; Q_tst=0; Pjac_tst=0; Qjac_tst=0;


   ! Set the initial conditions
   call WH_initial_data(P,Q)
   call WH_initial_data(P_tst,Q_tst)

   
   ! Perturb the initial data for the test trajectory (but only for pluto)
   Q_tst(4,1) = Q_tst(4,1) + traj_pert
   !call print_vector(Q(:,1))
   !call print_vector(P(:,1))

   
   ! set up jacobi coordinate functions
   call jacobi_setup(jacQ,jacP,jact,PjacQ,LUjacQ,PjacP,LUjacP,eta,m_vec_jac,g_param_jac)

   Q_wrk = Q(:,1); P_wrk = P(:,1);
   Q_tst_wrk = Q_tst(:,1); P_tst_wrk = P_tst(:,1);

   call apply_jacobi(Q(:,1),P(:,1),Qjac_wrk,Pjac_wrk,jacQ,jacP)
   Qjac(:,1) = Qjac_wrk; Pjac(:,1) = Pjac_wrk;
   call apply_jacobi(Q_tst(:,1),P_tst(:,1),Qjac_tst_wrk,Pjac_tst_wrk,jacQ,jacP)
   Qjac_tst(:,1) = Qjac_tst_wrk; Pjac_tst(:,1) = Pjac_tst_wrk;

   
   ! Do half a step of kepler
   call kepler_step(Qjac_wrk, Pjac_wrk, 0.5_dblk*dt,&
      kep_r,kep_v,kep_u,kep_a,kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,&
      kep_C,kep_S,kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)
   call kepler_step(Qjac_tst_wrk, Pjac_tst_wrk, 0.5_dblk*dt,&
      kep_r,kep_v,kep_u,kep_a,kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,&
      kep_C,kep_S,kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)
  
   ! One step of SI
   call symp_step(Qjac_wrk,Pjac_wrk,dt,symp_interdv,&
      symp_interdvjac,symp_Q_wrk,symp_qimqj,symp_qimqjnrm,jacP,PjacQ,&
      LUjacQ,jacT,symp_ind_wrk1,symp_ind_wrk2,m_vec,m_vec_jac,g_const)
   call symp_step(Qjac_tst_wrk,Pjac_tst_wrk,dt,symp_interdv,&
      symp_interdvjac,symp_Q_wrk,symp_qimqj,symp_qimqjnrm,jacP,PjacQ,&
      LUjacQ,jacT,symp_ind_wrk1,symp_ind_wrk2,m_vec,m_vec_jac,g_const)


   ! Finish the integration with a half-step (0.5dt) of kepler
   call kepler_step(Qjac_wrk, Pjac_wrk, 0.5_dblk*dt,&
      kep_r,kep_v,kep_u,kep_a,kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,&
      kep_C,kep_S,kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)
   call kepler_step(Qjac_tst_wrk, Pjac_tst_wrk, 0.5_dblk*dt,&
      kep_r,kep_v,kep_u,kep_a,kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,&
      kep_C,kep_S,kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)
   
   Qjac(:,2) = Qjac_wrk; Pjac(:,2) = Pjac_wrk;
   call apply_jacobi_inv(Qjac_wrk,Pjac_wrk,Q(:,2),P(:,2),&
      PjacQ,LUjacQ,PjacP,LUjacP)

   Qjac_tst(:,2) = Qjac_tst_wrk; Pjac_tst(:,2) = Pjac_tst_wrk;
   call apply_jacobi_inv(Qjac_tst_wrk,Pjac_tst_wrk,Q_tst(:,2),P_tst(:,2),&
      PjacQ,LUjacQ,PjacP,LUjacP)


   call save_orbit_lyapunov(savefile,t,Q,P,Qjac,Pjac,Q_tst,P_tst,&
      Qjac_tst,Pjac_tst,jacQ,jacP,jact,PjacQ,LUjacQ,PjacP,LUjacP,&
      m_vec,m_vec_jac,g_const,g_param,g_param_jac,ps_dist_rat)
 
   print *, "Computation complete!"

end program orbit_lin_test
