! symplectic.f90

module symplectic

   use orbital_data
   use utils
   use coords

   implicit none

   contains

   ! Driver routine for symplectic integration of interaction terms
   subroutine symp_step(Qjac_wrk,Pjac_wrk,tau,interdv,interdvjac,Q_wrk,qimqj,&
         qimqjnrm,jacQ,PjacQ,LUjacQ,jacT,ind_wrk1,ind_wrk2,m_vec,m_vec_jac,g_const)
      real (kind=dblk), intent(inout) :: Qjac_wrk(:),Pjac_wrk(:),interdv(:),&
         interdvjac(:),Q_wrk(:),qimqj(:,:),qimqjnrm(:,:),jacQ(:,:),&
         LUjacQ(:,:),jacT(:,:), ind_wrk1(:),ind_wrk2(:)
      real (kind=dblk), intent(in) :: tau,m_vec(:),m_vec_jac(:),g_const
      integer (kind=intk) :: PjacQ(:)


      integer (kind=intk) :: i,j
      real (kind=dblk) :: temp(3), qhnrm(n_masses),qjacnrm(n_masses)

      ! first step of Ruth's 2nd order
      ! only do center of mass term (see eqn. 4.22 in JMF senior thesis)
      Qjac_wrk(1:3) = Qjac_wrk(1:3) + tau/2.0_dblk*Pjac_wrk(1:3)/m_vec_jac(1)
      
      ! Calculate the gradient of the potential for the "kick" to Pjac
      ! Convert Qjac_wrk from Jacobi
      call apply_jacobi_invqp(Qjac_wrk,Q_wrk,PjacQ,LUjacQ)


      ! Direct terms
      !qimqj = 0.0_dblk
      !qimqjnrm = 0.0_dblk
      interdv = 0.0_dblk
      interdvjac = 0.0_dblk
     
      ! r_i - r_j in heliocentric cartesian coordinates
      ! store in lower triangle
      ! qimqj(i,j) = - qimqj(i,j)
      do j=1,n_masses-1 ! don't need sun terms
         temp = Q_wrk(3*j+1:3*j+3)
         do i=j,n_masses-1
            qimqj(3*i+1:3*i+3,j+1) = Q_wrk(3*i+1:3*i+3) - temp
            qimqjnrm(i+1,j+1) = norm2(qimqj(3*i+1:3*i+3,j+1))
         end do
      end do
      !call print_array(qimqj)
      !call print_array(qimqjnrm)


      ! find dv_i/dt from Saha 1994
      ! first entry
      do i=2,n_masses-1
         interdv(4:6) = interdv(4:6) - g_const*m_vec(i+1)*&
            qimqj(3*i+1:3*i+3,2) / (qimqjnrm(i+1,2)**3)
      end do

      ! 'middle' entries
      ! This was kind of a 'bear' to work out...
      do j=2,n_masses-2
        
         ! sweep through interdv
         ! add left parts
         do i=j,n_masses-2
            interdv(3*i+1:3*i+3) = interdv(3*i+1:3*i+3) + g_const*m_vec(j)*&
               qimqj(3*i+1:3*i+3,j) / (qimqjnrm(i+1,j)**3)
         end do
        
         ! add right parts (but sweep down to do it, using anti-symmetry of
         ! qimqj)
         do i=j+1,n_masses-1
            interdv(3*j+1:3*j+3) = interdv(3*j+1:3*j+3) - g_const*m_vec(i+1)*&
               qimqj(3*i+1:3*i+3,j+1) / (qimqjnrm(i+1,j+1)**3)
         end do
      end do

      ! final entry
      do i=1,n_masses-2
         interdv(3*n_masses-2:3*n_masses) = interdv(3*n_masses-2:3*n_masses) &
            + g_const*m_vec(i+1)*qimqj(3*n_masses-2:3*n_masses,i+1)&
            / (qimqjnrm(n_masses,i+1)**3)
      end do

      !call print_vector(interdv)
     
      ! convert to Jacobi coords
      call apply_jacobiqp(interdv,interdvjac,jacQ)
      ! make it a momentum
      do i=1,n_masses-1
         interdvjac(3*i+1:3*i+3) = interdvjac(3*i+1:3*i+3) * m_vec_jac(i+1) 
      end do

      ! Indirect terms
      ! Saha 1994; Encke's method for a better way of doing r'/r'^3 - r/r^3
      ! Encke's method is discussed in Danby - Fundamentals of Celestial Mech.
      temp = Q_wrk(1:3) ! Sun's position (we're moving to helio-centric)
      do i=1,n_masses-1
         ind_wrk1(3*i+1:3*i+3) = Q_wrk(3*i+1:3*i+3) - temp ! heliocentric
         qhnrm(i+1) = norm2(ind_wrk1(3*i+1:3*i+3))
         qjacnrm(i+1) = norm2(Qjac_wrk(3*i+1:3*i+3))
         ind_wrk1(3*i+1:3*i+3) = ind_wrk1(3*i+1:3*i+3) / (qhnrm(i+1)**3)
      end do
      ! we don't use 1st entry of qhnrm or qjacnrm
      ind_wrk1(1:3) = 0.0_dblk

      call apply_jacT(ind_wrk1,ind_wrk2,jacT)
     
     
      do i=1,n_masses-1
         ind_wrk1(3*i+1:3*i+3) = -g_const*m_vec(1)*m_vec(i+1) &
            / (qjacnrm(i+1)**3) * (Qjac_wrk(3*i+1:3*i+3) - qjacnrm(i+1)**3*&
            ind_wrk2(3*i+1:3*i+3))
      end do
      ind_wrk1(1:3) = 0.0_dblk
      
      interdvjac = interdvjac + ind_wrk1

      ! Apply the rest of Ruth's 2nd order
      ! "Kick" Pjac
      Pjac_wrk = Pjac_wrk - tau*interdvjac

      ! center of mass
      Qjac_wrk(1:3) = Qjac_wrk(1:3) + tau/2.0_dblk*Pjac_wrk(1:3)/m_vec_jac(1)

   end subroutine symp_step

end module symplectic
