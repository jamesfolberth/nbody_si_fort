! orbit.f90

program main
   
   use HDF5 ! HDF5 for saving data
   use utils ! helper routines, data types
   use orbital_data


   ! run config parameters
   character (len=255), parameter :: savefile='data/orbit.h5'
   integer (kind=intk), parameter :: N_save_int = 1000 ! save data every N_save_int times iterations
   
   real (kind=dblk), parameter :: t0 = 1941.+6./365.25 ! JD=2430000.5
   real (kind=dblk), parameter :: t1 = t0+10**4 ! JD=2430000.5
   real (kind=dblk), parameter :: dt = 1 ! time step
   integer (kind=intk), parameter :: N_saves = ceiling((t1-t0)/N_save_int)
   
   ! variables
   real (kind=dblk), allocatable :: t(:)
   real (kind=dblk), allocatable :: P(:,:),Q(:,:),Pjac(:,:),Qjac(:,:)

   integer (kind=intk) :: i,j
   real (kind=dblk) :: ti


   ! initialize 
   ! TODO time handling should be done differently
   !      write time at each save, not before iterations   
   !t = linspace(t0,t1,floor((t1-t0)/(dt_save_init)))
   !dt_save = t(2)-t(1)
   !dt = dt_save/N_step
   allocate(t(N_saves))

   if (debug) then
      print *, "Num iterations (approx) = ", floor((t1-t0)/dt)
   end if

   allocate(P(3*n_masses, size(t,1)), Q(3*n_masses, size(t,1)))
   allocate(Pjac(3*n_masses, size(t,1)), Qjac(3*n_masses, size(t,1)))
   P=0; Q=0; Pjac=0; Qjac=0;

  
   ! main loop
   i = 0
   ti = t0
   do while (ti < t1)

      do j=1,N_save_int
         
         ti = ti + dt;
      end do

      i=i+1; t(i)=ti;
   end do


end program main
