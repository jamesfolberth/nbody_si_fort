! lyapunov.f90

! using the "variational method" of Sussman/Wisdom 1988, estimate
! the largest positive lyapunov exponent of Pluto.
! Due to the large number of calls to DLSODE, this code is considerably
! slower than orbit.f90

program lyapunov
   
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
   character (len=256), parameter :: savefile='../data/lyapunov_dt_0.1.h5'
   integer (kind=intk), parameter :: N_record_int = 1000 ! record state every N_record_int time steps
   integer (kind=intk), parameter :: N_saves = 1000
   
   real (kind=dblk), parameter :: t0 = 1941.+6./365.25 ! JD=2430000.5
   real (kind=dblk), parameter :: t1 = t0+10**9 ! JD=2430000.5
   real (kind=dblk), parameter :: dt = 0.50_dblk ! time step
   real (kind=dblk), parameter :: traj_pert = 1.0E0_dblk
   integer (kind=intk), parameter :: N_records = ceiling((t1-t0)/(dt*N_record_int))
   integer (kind=intk), parameter :: N_save_int = ceiling(dble(N_records)/dble(N_saves))


   ! DLSODE setup
   external :: lin_func, lin_jac
   integer (kind=intk) :: dlsode_neq(1) = (/6*n_masses/)
   integer (kind=intk) :: dlsode_itol = 1
   real (kind=dblk) :: dlsode_rtol(1)=(/1E-10_dblk/)
   real (kind=dblk) :: dlsode_atol(1)=(/1E-10_dblk/)
   integer (kind=intk) :: dlsode_itask = 1
   integer (kind=intk) :: dlsode_istate = 1
   integer (kind=intk) :: dlsode_iopt = 0
   real (kind=dblk) :: dlsode_rwork(20+16*6*n_masses)
   integer (kind=intk) :: dlsode_lrw = 20+16*6*n_masses
   integer (kind=intk) :: dlsode_iwork(20)
   integer (kind=intk) :: dlsode_liw = 20
   integer (kind=intk) :: dlsode_mf = 10

   real (kind=dblk) :: dlsode_t,temp
   real (kind=dblk) :: dqp_lin(2*6*n_masses+n_masses+2)
   real (kind=dblk), allocatable :: dqp_norm(:),dqp_norm_wrk,dqp_norm_prev,&
      lambda(:)
   ! dqp_lin storage
   ! [ dq dp | q_0 p_0 m_vec g_const build_d2Hdq2_flag ]^T
   !    neq  | other variables


   ! variables
   real (kind=dblk), allocatable :: t(:)

   ! reference trajectory (similar to orbit.f90, except that we don't
   ! store the position and momentum at each time step)
   real (kind=dblk) :: P_wrk(3*n_masses),Q_wrk(3*n_masses),&
      Pjac_wrk(3*n_masses),Qjac_wrk(3*n_masses)
   real (kind=dblk), allocatable :: P(:,:),Q(:,:),Pjac(:,:),Qjac(:,:)
  
   ! coordinate transformation variables
   real (kind=dblk) :: jacQ(3*n_masses,3*n_masses),&
      LUjacQ(3*n_masses,3*n_masses),jacP(3*n_masses,3*n_masses),&
      LUjacP(3*n_masses,3*n_masses),jacT(3*n_masses,3*n_masses),&
      eta(n_masses),m_vec_jac(n_masses), g_param_jac(n_masses)
   integer (kind=intk):: PjacQ(3*n_masses), PjacP(3*n_masses)

   ! variables 'local' to kepler_step
   ! see the note in orbit.f90
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

   ! initialize (allocatable) variables
   allocate(t(N_records),dqp_norm(N_records),lambda(N_records))
   allocate(P(3*n_masses,1), Q(3*n_masses,1))
   allocate(Pjac(3*n_masses,1), Qjac(3*n_masses,1))
   P=0.0_dblk; Q=0.0_dblk; Pjac=0.0_dblk; Qjac=0.0_dblk;


   ! Set the initial conditions
   call WH_initial_data(P,Q)

   
   ! set up jacobi coordinate functions
   call jacobi_setup(jacQ,jacP,jact,PjacQ,LUjacQ,PjacP,LUjacP,eta,m_vec_jac,g_param_jac)

   Q_wrk = Q(:,1); P_wrk = P(:,1);

   call apply_jacobi(Q(:,1),P(:,1),Qjac_wrk,Pjac_wrk,jacQ,jacP)
   Qjac(:,1) = Qjac_wrk; Pjac(:,1) = Pjac_wrk;


   dqp_lin = 0.0_dblk
   dqp_lin(4) = traj_pert
   dqp_lin(5) = -0.5*traj_pert
   dqp_lin(6) = -0.239*traj_pert
   dqp_lin(22) = traj_pert
   dqp_lin(23) = 0.52333423*traj_pert
   dqp_norm_prev = dsqrt(norm2(dqp_lin(4:6))**2+norm2(dqp_lin(22:24))**2)
   dqp_lin(1:6*n_masses) = dqp_lin(1:6*n_masses)/dqp_norm_prev
   dqp_lin(37:54) = Q_wrk
   dqp_lin(55:72) = P_wrk
   dqp_lin(73:78) = m_vec
   dqp_lin(79) = g_const
   dqp_lin(80) = -1.0_dblk
   
   
   ! Intgrate dqp
   dqp_lin(80) = -1.0_dblk
   dlsode_t = 0.0_dblk
   dlsode_istate=1

   call dlsode(lin_func,dlsode_neq,dqp_lin,dlsode_t,0.5_dblk*dt,&
      dlsode_itol,dlsode_rtol,dlsode_atol,dlsode_itask,dlsode_istate,&
      dlsode_iopt,dlsode_rwork,dlsode_lrw,dlsode_iwork,dlsode_liw,&
      lin_jac,dlsode_mf)
 
   dqp_norm = 1.0_dblk
   dqp_norm_wrk = hypot(norm2(dqp_lin(4:6)),norm2(dqp_lin(22:24)))&
      /dqp_norm_prev
   
   dqp_lin(1:6*n_masses) = dqp_lin(1:6*n_masses)/dqp_norm_wrk

   !print *,dlsode_istate ! istate=2 is success
   !call print_vector(dqp_lin(1:6*n_masses))

   i = 1
   iter = 1
   ti = t0
 
   ! Do half a step of kepler
   ! TODO should ti also be moved forward by 0.5dt?
   call kepler_step(Qjac_wrk, Pjac_wrk, 0.5_dblk*dt,&
      kep_r,kep_v,kep_u,kep_a,kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,&
      kep_C,kep_S,kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)
  

   call tic(clock)
   main: do while (ti < t1)

      integrate: do j=1,N_record_int

         ! Integrate the variational problem (with the old data)
         call apply_jacobi_inv(Qjac_wrk,Pjac_wrk,Q_wrk,P_wrk,&
            PjacQ,LUjacQ,PjacP,LUjacP)

         ! renormalize length of dqp_lin vector
         !temp = dsqrt(norm2(dqp_lin(4:6))**2+norm2(dqp_lin(22:24))**2)
         temp = hypot(norm2(dqp_lin(4:6)),norm2(dqp_lin(22:24)))
         dqp_lin(1:6*n_masses) = dqp_lin(1:6*n_masses) / temp
         dqp_norm_wrk = dqp_norm_wrk * temp

         ! store reference trajectory data (computed with SI method)
         dqp_lin(37:54) = Q_wrk
         dqp_lin(55:72) = P_wrk
         dqp_lin(80) = -1.0_dblk ! build flag
         dlsode_istate=1

         call dlsode(lin_func,dlsode_neq,dqp_lin,dlsode_t,dlsode_t+dt,&
            dlsode_itol,dlsode_rtol,dlsode_atol,dlsode_itask,dlsode_istate,&
            dlsode_iopt,dlsode_rwork,dlsode_lrw,dlsode_iwork,dlsode_liw,&
            lin_jac,dlsode_mf)
      
         !call print_vector(dqp_lin(1:36))
     

         ! second-order SI method (with 0.5*dt from above and below)
         ! Do a full step of SI
         call symp_step(Qjac_wrk,Pjac_wrk,dt,symp_interdv,&
            symp_interdvjac,symp_Q_wrk,symp_qimqj,symp_qimqjnrm,jacP,PjacQ,&
            LUjacQ,jacT,symp_ind_wrk1,symp_ind_wrk2,m_vec,m_vec_jac,g_const)

         ! Do full step of kepler
         call kepler_step(Qjac_wrk, Pjac_wrk,dt,&
            kep_r,kep_v,kep_u,kep_a,kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,&
            kep_C,kep_S,kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)

         ti = ti + dt

      end do integrate
   
      !call toc(clock)
      !stop('terminating integration after one loop of integration (nothing is wrong)')

      !print *, ti,dlsode_t,norm2(dqp_lin(1:6*n_masses))
      if (dlsode_istate /= 2) stop

      dqp_norm(i) = dqp_norm_wrk
      lambda(i) = 1.0_dblk/(ti-t0)*dlog(dqp_norm(i))

      t(i)=ti;
      i=i+1; 

      ! Record state of system every once in a while
      ! note that we don't store the reference trajector; only enough
      ! information to analyze the lyapunov exponent data
      if (mod(i,floor(dble(N_save_int)/100_dblk)) == 0) then
         print *, "Percent complete: ", ti/t1*100_dblk
         print *, temp,dqp_norm(i-1),lambda(i-1)

         ! save data every N_save_int iterations of i
         if (mod(i,N_save_int) == 0) then
            print *, "saving"
            call save_orbit_lyapunov(savefile,t,N_record_int,dqp_norm,lambda)
 
         end if
      end if

   end do main
   call toc(clock)

   ! Finish the integration with a half-step (0.5dt) of kepler
   ! (this isn't needed for lyapunov calculations)
   call kepler_step(Qjac_wrk, Pjac_wrk, 0.5_dblk*dt,&
      kep_r,kep_v,kep_u,kep_a,kep_n,kep_EC,kep_ES,kep_e,kep_dE,kep_dtv,&
      kep_C,kep_S,kep_f,kep_g,kep_aor,kep_fp,kep_gp,m_vec_jac,g_param_jac)
   
   call save_orbit_lyapunov(savefile,t,N_record_int,dqp_norm,lambda)
 
   print *, "Computation complete!"

   deallocate(t,dqp_norm,lambda)

end program lyapunov


! This is the linearization of the Hamiltonian ODEs
subroutine lin_func(neq,t,y,ydot)
   integer (kind=4) :: neq
   real (kind=8) :: t,y(neq),ydot(neq)

   integer (kind=4) :: i,n_masses,m,k
   real (kind=8) :: q_0(18),p_0(18),m_vec(6),g_const
   real (kind=8) :: qmmqk(3),rmk,temp

   real (kind=8), save :: d2Hdq2(18,18)

   ! dqp_lin = y storage
   ! [ dq dp | q_0 p_0 m_vec g_const build_d2Hdq2_flag ]^T
   !    neq  | other variables  <-- sizes

   ! if build_d2Hdq2 <= 0.0d0, then we build d2Hdq2
   ! We can probably get away with only building it once (due to short 
   ! time steps)
   ! TODO This is probably true, but test it to verify.

   n_masses = 6
   g_const = y(79)
   m_vec = y(73:78)
   p_0 = y(55:72)
   q_0 = y(37:54)

   ! dq' = d^2H/dp^2 * dp
   do i=0,n_masses-1
      ydot(3*i+1) = 1/m_vec(i+1)*y(3*n_masses+3*i+1)
      ydot(3*i+2) = 1/m_vec(i+1)*y(3*n_masses+3*i+2)
      ydot(3*i+3) = 1/m_vec(i+1)*y(3*n_masses+3*i+3)
   end do
   !ydot(1) = 1/m_vec(2)*y(4)
   !ydot(2) = 1/m_vec(2)*y(5)
   !ydot(3) = 1/m_vec(2)*y(6)


   ! dp' = -d^2H/dq^2 * dq
   ! Build d2Hdq2
   if (y(80) <= 0.0d0) then  ! this is the build_d2Hdq2 flag
      d2Hdq2 = 0.0d0 ! local variables are saved after each call, so 
                     ! be sure to init to zero!
      do m=0,n_masses-1
         do k=0,m-1
            qmmqk = q_0(3*m+1:3*m+3) - q_0(3*k+1:3*k+3)
            rmk = norm2(qmmqk)
            temp = g_const*m_vec(m+1)*m_vec(k+1)/rmk**3
            d2Hdq2(3*m+1:3*m+3,3*k+1) = temp*((/1.0d0,0.0d0,0.0d0/)&
               -3.0d0/(rmk**2)*qmmqk(1)*qmmqk)
            d2Hdq2(3*m+1:3*m+3,3*k+2) = temp*((/0.0d0,1.0d0,0.0d0/)&
               -3.0d0/(rmk**2)*qmmqk(2)*qmmqk)
            d2Hdq2(3*m+1:3*m+3,3*k+3) = temp*((/0.0d0,0.0d0,1.0d0/)&
               -3.0d0/(rmk**2)*qmmqk(3)*qmmqk)

            d2Hdq2(3*k+1:3*k+3,3*m+1) = temp*((/1.0d0,0.0d0,0.0d0/)&
               -3.0d0/(rmk**2)*qmmqk(1)*qmmqk)
            d2Hdq2(3*k+1:3*k+3,3*m+2) = temp*((/0.0d0,1.0d0,0.0d0/)&
               -3.0d0/(rmk**2)*qmmqk(2)*qmmqk)
            d2Hdq2(3*k+1:3*k+3,3*m+3) = temp*((/0.0d0,0.0d0,1.0d0/)&
               -3.0d0/(rmk**2)*qmmqk(3)*qmmqk)
         end do
      end do

      ! Need to init diagonal blocks to zero first!!
      do m=0,n_masses-1
         do k=0,n_masses-1
            if (m /= k) then
               d2Hdq2(3*m+1:3*m+3,3*m+1:3*m+3) = &
                  d2Hdq2(3*m+1:3*m+3,3*m+1:3*m+3)&
                  - d2Hdq2(3*k+1:3*k+3,3*m+1:3*m+3)
            end if
         end do
      end do

      y(80) = -1.0d0 ! Change the build flag so we don't build again for this
                    ! call
   end if

   ! d2Hdq2 is dense
   ! N = 18 is near the cutoff for using matmul vs. dgemv
   call dgemv('N',18,18,1.0d0,d2Hdq2,18,y(1),1,0.0d0,ydot(19),1)

end subroutine lin_func

! Jacobian for lin_func; will not be called for dlsode_mf=10 (Adam's 
! multistep method)
subroutine lin_jac(neq,t,y,ml,mu,pd,nrowpd)
   integer (kind=4) :: neq,ml,mu,nrowpd
   real (kind=8) :: t,y(neq),pd(nrowpd,neq)

   neq = 0;ml=0;mu=0;nrowpd=0;
   t=0d0;y=0d0;pd=0d0;
   print *, "calling lin_jac"
end subroutine
