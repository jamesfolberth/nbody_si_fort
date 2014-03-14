! utils.f90

module utils
   
   use HDF5
   implicit none

   integer, parameter :: intk = kind(1) ! get kind numbers of various types
   integer, parameter :: dblk = kind(1d0)
   real (kind=dblk), parameter :: pi = 4.0_dblk*atan(1.0_dblk)
   logical, parameter :: debug = .true.


   public write_dset
   private dwrite_dset_rank1, iwrite_dset_rank1,&
           dwrite_dset_rank2, dwrite_dset_rank0

   interface write_dset
      module procedure dwrite_dset_rank1, iwrite_dset_rank1,&
                       dwrite_dset_rank2, dwrite_dset_rank0
   end interface write_dset


   contains

   !!!!!!!!!!!!!!!!!!!!
   ! Mat/Vec routines !
   !!!!!!!!!!!!!!!!!!!!
   ! {{{
   ! linspace
   !  behaves like Matlab's linspace
   function linspace(x_start, x_end, N)
      real (kind=dblk), intent(in) :: x_start, x_end
      integer (kind=intk), intent(in) :: N

      real (kind=dblk), allocatable :: linspace(:)
      real (kind=dblk) :: h
      integer (kind=intk) :: i

      allocate(linspace(N))
      h = (x_end-x_start) / dble(N-1)

      linspace(1) = x_start
      do i=1+1,N-1
         linspace(i) = linspace(i-1)+h
      end do
      linspace(N)=x_end
   end function linspace
   ! }}}


   !!!!!!!!!!
   ! Timing !
   !!!!!!!!!!
   ! {{{
   ! behaves like Octave's tic()
   subroutine tic(t)
      integer, intent(out) :: t
      call system_clock(t)
   end subroutine tic
   
   ! returns time in seconds from now to time described by t
   real function toc_return(t)
      integer, intent(in) :: t
      integer :: now, clock_rate
   
      call system_clock(now, clock_rate)
      toc_return = real(now - t)/real(clock_rate)
   end function toc_return
   
   ! prints time in seconds from now to time described by t
   subroutine toc(t)
      integer, intent(in) :: t
      real (kind=8) :: time
      integer :: now, clock_rate
   
      call system_clock(now, clock_rate)
      time = real(now - t)/real(clock_rate)
      print *,"Elapsed time is ",time," seconds."
   end subroutine toc
   ! }}}


   !!!!!!!!!!
   ! Stdout !
   !!!!!!!!!!
   ! {{{
   ! print_vector
   !  print a rank 1 array
   subroutine print_vector(v)
      real (kind=dblk), intent(in) :: v(:)
      character (len=30) :: rowfmt

      integer (kind=intk) :: i
   
      write(rowfmt, "(A,I4,A)") "(",1,"(1X,SS,10Es13.4))"
      print *,
      row_print: do i=1,size(v,1)
         write(*, fmt=rowfmt) v(i)
      end do row_print
   end subroutine print_vector

   ! print_array
   !  print a rank 2 array
   subroutine print_array(A)
      real (kind=dblk), intent(in) :: A(:,:)
      character (len=30) :: rowfmt

      integer (kind=intk) :: i,j
   
      write(rowfmt, "(A,I4,A)") "(",size(A,2),"(1X,SS,10Es13.4))"
      print *,
      row_print: do i=1,size(A,1)
         write(*, fmt=rowfmt) (A(i,j), j=1,size(A,2))
      end do row_print
   end subroutine print_array
   ! }}}


   !!!!!!!!
   ! HDF5 !
   !!!!!!!!
   ! {{{
   subroutine save_data(savefile,t,Q,P,Qjac,Pjac,jacQ,jacP,jacT,&
                        PjacQ,LUjacQ,PjacP,LUjacP,m_vec,m_vec_jac,&
                        g_const,g_param,g_param_jac)
      character (len=256) :: savefile
      real (kind=dblk) :: t(:),Q(:,:),P(:,:),Qjac(:,:),Pjac(:,:),&
                          jacQ(:,:),jacP(:,:),jacT(:,:),&
                          LUjacQ(:,:),LUjacP(:,:),&
                          m_vec(:),m_vec_jac(:),&
                          g_const,g_param,g_param_jac(:)
      integer (kind=intk) :: PjacQ(:),PjacP(:)

      integer (kind=intk) :: h5error, file_id

      call h5open_f(h5error)
      call h5fcreate_f(savefile, H5F_ACC_TRUNC_F, file_id, h5error)

      ! start writing data
      call write_dset(file_id, t, "t")
      call write_dset(file_id, Q, "Q")
      call write_dset(file_id, P, "P")
      call write_dset(file_id, Qjac, "Qjac")
      call write_dset(file_id, Pjac, "Pjac")
      call write_dset(file_id, jacQ, "jacQ")
      call write_dset(file_id, jacP, "jacP")
      call write_dset(file_id, jacT, "jacT")
      call write_dset(file_id, PjacQ, "PjacQ")
      call write_dset(file_id, LUjacQ, "LUjacQ")
      call write_dset(file_id, PjacP, "PjacP")
      call write_dset(file_id, LUjacP, "LUjacP")
      call write_dset(file_id, m_vec, "m_vec")
      call write_dset(file_id, m_vec_jac, "m_vec_jac")
      call write_dset(file_id, g_const, "g_const")
      call write_dset(file_id, g_param, "g_param")
      call write_dset(file_id, g_param_jac, "g_param_jac")
    
      call h5fclose_f(file_id,h5error)
      call h5close_f(h5error)

   end subroutine save_data

  
   ! These routines make the dspace, dset, and write data
   ! they all interface to the 'write_dset' generic routine
   subroutine dwrite_dset_rank0(file_id, scalar, dsetname)
      integer (kind=intk) :: file_id
      real (kind=dblk) :: scalar
      character (len=*) :: dsetname

      integer (kind=intk) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(1)

      dims(1) = 0
      call h5screate_simple_f(0,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,scalar,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)

   end subroutine dwrite_dset_rank0

   subroutine dwrite_dset_rank1(file_id, array, dsetname)
      integer (kind=intk) :: file_id
      real (kind=dblk) :: array(:)
      character (len=*) :: dsetname

      integer (kind=intk) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(1)

      dims(1) = size(array)
      call h5screate_simple_f(1,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)

   end subroutine dwrite_dset_rank1

   subroutine iwrite_dset_rank1(file_id, array, dsetname)
      integer (kind=intk) :: file_id
      integer (kind=intk) :: array(:)
      character (len=*) :: dsetname

      integer (kind=intk) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(1)

      dims(1) = size(array)
      call h5screate_simple_f(1,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_INTEGER,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,array,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)

   end subroutine iwrite_dset_rank1


   subroutine dwrite_dset_rank2(file_id, array, dsetname)
      integer (kind=intk) :: file_id
      real (kind=dblk) :: array(:,:)
      character (len=*) :: dsetname

      integer (kind=intk) :: dspace_id, dset_id, h5error
      integer (kind=HSIZE_T) :: dims(2)

      dims(1) = size(array,1); dims(2) = size(array,2)
      call h5screate_simple_f(2,dims,dspace_id,h5error)
      call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,h5error)
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,array,dims,h5error)
      call h5dclose_f(dset_id,h5error)
      call h5sclose_f(dspace_id,h5error)

   end subroutine dwrite_dset_rank2
   ! }}}


end module utils

