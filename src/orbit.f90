! orbit.f90

program main
   
   use HDF5 ! HDF5 for saving data
   use utils ! helper routines, data types
   use orbital_data
   use coords

   implicit none


   ! run config parameters
   character (len=256), parameter :: savefile='../data/orbit.h5'
   integer (kind=intk), parameter :: N_record_int = 1000 ! record state every N_record_int time steps
   integer (kind=intk), parameter :: N_saves = 10
   
   real (kind=dblk), parameter :: t0 = 1941.+6./365.25 ! JD=2430000.5
   real (kind=dblk), parameter :: t1 = t0+10**5 ! JD=2430000.5
   real (kind=dblk), parameter :: dt = 1 ! time step
   integer (kind=intk), parameter :: N_records = ceiling((t1-t0)/N_record_int)
   integer (kind=intk), parameter :: N_save_int = ceiling(dble(N_records)/dble(N_saves))
   
   ! variables
   real (kind=dblk), allocatable :: t(:)
   real (kind=dblk), allocatable :: P(:,:),Q(:,:),Pjac(:,:),Qjac(:,:)
   real (kind=dblk), allocatable :: jacQ(:,:), LUjacQ(:,:), &
                                    jacP(:,:), LUjacP(:,:), &
                                    jacT(:,:), & ! used in interaction_dV 
                                    eta(:),m_vec_jac(:), g_param_jac(:)
   integer (kind=intk), allocatable :: PjacQ(:), PjacP(:)

   integer (kind=intk) :: i,j, save_counter
   real (kind=dblk) :: ti


   ! initialize variables
   allocate(t(N_records))
   allocate(P(3*n_masses, N_records), Q(3*n_masses, N_records))
   allocate(Pjac(3*n_masses, N_records), Qjac(3*n_masses, N_records))
   allocate(jacQ(3*n_masses, 3*n_masses), LUjacQ(3*n_masses,3*n_masses))
   allocate(jacP(3*n_masses, 3*n_masses), LUjacP(3*n_masses,3*n_masses))
   allocate(PjacQ(3*n_masses), PjacP(3*n_masses))
   allocate(jacT(3*n_masses,3*n_masses))
   allocate(eta(n_masses), m_vec_jac(n_masses), g_param_jac(n_masses))

   P=0; Q=0; Pjac=0; Qjac=0;

   call WH_initial_data(P,Q)



   ! set up jacobi coordinate functions
   call jacobi_setup(jacQ,jacP,jact,PjacQ,LUjacQ,PjacP,LUjacP,eta,m_vec_jac,g_param_jac)

   call apply_jacobi(Q(:,1),P(:,1),Qjac(:,1),Pjac(:,1),jacQ,jacP)

   call save_data(savefile,t,Q,P,Qjac,Pjac,jacQ,jacP,jact,PjacQ,LUjacQ,PjacP,LUjacP,m_vec,m_vec_jac,g_const,g_param)

  
   ! main loop
   if (debug) then
      write (*,"(A,10ES13.4)"), "Num iterations (approx) = ",(t1-t0)/dt
   end if

   i = 0
   save_counter = 0
   ti = t0
   do while (ti < t1)

      do j=1,N_record_int
         
         ti = ti + dt;
      end do

      ! Record state of system
      i=i+1; t(i)=ti;
      
      ! save data every N_save_int iterations of i
      if (mod(i,N_save_int) == 0) then
         !if (debug) print *, "saving data (but not really)"
      end if
   end do


end program main
