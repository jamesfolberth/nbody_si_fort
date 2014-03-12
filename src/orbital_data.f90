! orbital_data.f90

module orbital_data
   
   use utils
  
   ! Various solar system constants
   real (kind=dblk), parameter :: AU = 1.4959787066E11
   real (kind=dblk), parameter :: m_earth = 5.9736E24
   real (kind=dblk), parameter :: year_day = 365.25636004
   real (kind=dblk), parameter :: year_sec = year_day*24*3600
   real (kind=dblk), parameter :: m_sun = 1.9891E30/m_earth
   real (kind=dblk), parameter :: g_const = 6.67428E-11 / (AU**3) * m_earth * year_sec**2
   real (kind=dblk), parameter :: g_param = g_const*m_sun
  
   ! masses from Carroll and Ostlie
   ! Sun Pluto Jupiter ... Neptune (in Earth masses)
   real (kind=dblk), parameter, dimension(6) :: m_vec = (/ m_sun, 0.00218d0, 317.83d0, 95.159d0, 14.536d0, 17.147d0 /) ! earth masses
   integer (kind=intk), parameter :: n_masses = size(m_vec,1)

   contains

      ! Wisdom/Holman initial data 
      ! TODO which paper?
      subroutine WH_initial_data(P,Q)
         real (kind=dblk), intent(inout) :: P(:,:), Q(:,:)

         integer (kind=intk) :: i, ind0, ind1
         real (kind=dblk) :: M
         real (kind=dblk) :: qc(3), vc(3)

         ! Initial data
         ! {{{
         ! Sun
         i=0; ind0=3*i+1; ind1=3*i+3;
         ! need to convert from AU/year to AU/day
         P(ind0:ind1,1) = m_vec(i+1)*year_day*& 
            (/ 6.69048890636161E-6,&
               -6.33922479583593E-6,&
               -3.13202145590767E-9 /)
         Q(ind0:ind1,1) = &
            (/ -4.06428567034226E-3,&
               -6.08813756435987E-3,&
               -1.66162304225834E-6 /)

         ! Pluto
         i=1; ind0=3*i+1; ind1=3*i+3;
         P(ind0:ind1,1) = m_vec(i+1)*year_day*& 
            (/ -1.76936577252484E-3,&
               -2.06720938381724E-3,&
               6.58091931493844E-4 /)
         Q(ind0:ind1,1) = &
            (/ -2.13858977531573E1,&
               3.20719104739886E1,&
               2.49245689556096E0 /)

         ! Jupiter
         i=2; ind0=3*i+1; ind1=3*i+3;
         P(ind0:ind1,1) = m_vec(i+1)*year_day*& 
            (/ -5.5979796931066E-3,&
               5.51815399480116E-3,&
               -2.66711392865591E-6 /)
         Q(ind0:ind1,1) = &
            (/ 3.40546614227466E0,&
               3.62978190075864E0,&
               3.42386261766577E-2 /)

         ! Saturn
         i=3; ind0=3*i+1; ind1=3*i+3;
         P(ind0:ind1,1) = m_vec(i+1)*year_day*& 
            (/ -4.17354020307064E-3,&
               3.99723751748116E-3,&
               1.67206320571441E-5 /)
         Q(ind0:ind1,1) = &
            (/ 6.60801554403466E0,&
               6.38084674585064E0,&
               -1.3614596372452E-1 /)

         ! Uranus
         i=4; ind0=3*i+1; ind1=3*i+3;
         P(ind0:ind1,1) = m_vec(i+1)*year_day*& 
            (/ -3.25884806151064E-3,&
               2.06438412905916E-3,&
               -2.17699042180559E-5 /)
         Q(ind0:ind1,1) = &
            (/ 1.11636331405597E1,&
               1.60373479057256E1,&
               3.61783279369958E-1 /)

         ! Neptune
         i=5; ind0=3*i+1; ind1=3*i+3;
         P(ind0:ind1,1) = m_vec(i+1)*year_day*& 
            (/ -2.17471785045538E-4,&
               -3.11361111025884E-3,&
               3.58344705491441E-5 /)
         Q(ind0:ind1,1) = &
            (/ -3.01777243405203E1,&
               1.9115531499806E0,&
               -1.53887595621042E-1 /)

         ! }}}
     
         ! Subtract off motion of center of mass of solar system (this is an integral of motion)
         M = sum(m_vec)
         qc = 0; vc = 0;
         do i=1,n_masses
            qc = qc + m_vec(i)*Q(3*(i-1)+1:3*(i-1)+3,1)
            vc = vc + P(3*(i-1)+1:3*(i-1)+3,1)
         end do

         qc = qc/M; vc = vc/M;

         do i=0,n_masses-1
            Q(3*i+1:3*i+3,1) = Q(3*i+1:3*i+3,1) - qc
            P(3*i+1:3*i+3,1) = m_vec(i+1)*(P(3*i+1:3*i+3,1)/m_vec(i+1) - vc)
         end do

         if (debug) then
            if (size(Q,1) /= 3*n_masses) then
               print *, "error: orbital_data.f90: number of variables in Q does not match number of masses in m_vec"
            end if
         end if
            

 



 



      end subroutine WH_initial_data

end module orbital_data
