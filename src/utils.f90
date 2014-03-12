! utils.f90

module utils

   use HDF5

   integer, parameter :: intk = kind(1) ! get kind numbers of various types
   integer, parameter :: dblk = kind(1d0)
   logical, parameter :: debug = .true.

 
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
      real (kind=8), intent(in) :: v(:)
      character (len=30) :: rowfmt
   
      write(rowfmt, "(A,I4,A)") "(",1,"(1X,SS,10Es13.4))"
      print *,
      row_print: do i=1,size(v,1)
         write(*, fmt=rowfmt) v(i)
      end do row_print
   end subroutine print_vector

   ! print_array
   !  print a rank 2 array
   subroutine print_array(A)
      real (kind=8), intent(in) :: A(:,:)
      character (len=30) :: rowfmt
   
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

end module utils

