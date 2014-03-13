! coords.f90

module coords

   use utils
   use orbital_data

   contains

      ! Set up Jacobi coordinate transformation stuff
      subroutine jacobi_setup(jacQ,jacP,jacT,PjacQ,LUjacQ,PjacP,LUjacP,eta,m_vec_jac, g_param_jac)
         real (kind=dblk) :: jacQ(:,:), jacP(:,:), jacT(:,:), &
                             LUjacQ(:,:), LUjacP(:,:), &
                             eta(:), m_vec_jac(:), g_param_jac(:)
         integer (kind=intk) :: PjacQ(:), PjacP(:)

         integer (kind=intk) :: i, info

         ! cumulative masses
         eta(1) = m_vec(1)
         do i=2,n_masses
            eta(i) = eta(i-1) + m_vec(i)
         end do

         ! Jacobi masses
         m_vec_jac(1) = eta(n_masses)
         g_param_jac(1) = eta(1)*g_param
         do i=2,n_masses
            m_vec_jac(i) = eta(i-1)/eta(i)*m_vec(i)
            g_param_jac = eta(i)/eta(i-1)*g_param
         end do

         ! Qjac = jacQ*Q 
         ! Transform position to jacobi coordinates
         jacQ = 0
         ! First row
         jacQ(1,1:3*n_masses-2:3) = m_vec/eta(n_masses)
         jacQ(2,2:3*n_masses-1:3) = m_vec/eta(n_masses)
         jacQ(3,3:3*n_masses-0:3) = m_vec/eta(n_masses)
         do i=1,n_masses-1
            jacQ(3*i+1, 1:3*i:3) = -m_vec(1:i)/eta(i)
            jacQ(3*i+2, 2:3*i+1:3) = -m_vec(1:i)/eta(i)
            jacQ(3*i+3, 3:3*i+2:3) = -m_vec(1:i)/eta(i)

            jacQ(3*i+1,3*i+1) = 1d0
            jacQ(3*i+2,3*i+2) = 1d0
            jacQ(3*i+2,3*i+3) = 1d0
         end do

         ! Compute PLU decomp of jacQ (so we can solve later)
         LUjacQ = jacQ
         call dgetrf(3*n_masses, 3*n_masses, LUjacQ, 3*n_masses, PjacQ, info)
         if (info /= 0) then
            print *, "error: coords.f90: jacobi_setup: PLU of jacQ failed"
            stop
         end if

         ! Transform momentum to jacobi coordinates
         jacP = 0
         jacP(1,1:3*n_masses-2:3) = (/ (1d0,i=1,n_masses) /)
         jacP(2,2:3*n_masses-1:3) = (/ (1d0,i=1,n_masses) /)
         jacP(3,3:3*n_masses-0:3) = (/ (1d0,i=1,n_masses) /)
         do i=1,n_masses-1
            jacP(3*i+1,1:3*i+0:3) = -m_vec(i+1)/eta(i+1)
            jacP(3*i+2,2:3*i+1:3) = -m_vec(i+1)/eta(i+1)
            jacP(3*i+3,3:3*i+2:3) = -m_vec(i+1)/eta(i+1)

            jacP(3*i+1,3*i+1) = eta(i)/eta(i+1)
            jacP(3*i+2,3*i+2) = eta(i)/eta(i+1)
            jacP(3*i+3,3*i+3) = eta(i)/eta(i+1)
         end do

         ! Compute PLU decomp of jacP (so we can solve later)
         LUjacP = jacP
         call dgetrf(3*n_masses, 3*n_masses, LUjacP, 3*n_masses, PjacP, info)
         if (info /= 0) then
            print *, "error: coords.f90: jacobi_setup: PLU of jacP failed"
            stop
         end if

         ! Some other transform (Saha 199?)
         jacT = 0
         do i=1,3*n_masses
            jacT(i,i) = 1d0
         end do
         
         do i=1,n_masses-1
            jacT(3*i+1, 3*(i+1)+1:3*(n_masses-1)+1:3) = m_vec(i+2:n_masses)/eta(i+1)
            jacT(3*i+2, 3*(i+1)+2:3*(n_masses-1)+2:3) = m_vec(i+2:n_masses)/eta(i+1)
            jacT(3*i+3, 3*(i+1)+3:3*(n_masses-1)+3:3) = m_vec(i+2:n_masses)/eta(i+1)
         end do

      end subroutine jacobi_setup


      subroutine apply_jacobi(q,p,qjac,pjac,jacQ,jacP)
         real (kind=dblk) :: q(:), p(:), qjac(:), pjac(:),&
                             jacQ(:,:), jacP(:,:)

         call dgemv("N",3*n_masses,3*n_masses,1d0,jacQ,3*n_masses,q,1,0d0,qjac,1)
         call dgemv("N",3*n_masses,3*n_masses,1d0,jacP,3*n_masses,p,1,0d0,pjac,1)

      end subroutine apply_jacobi

      subroutine apply_jacobi_inv(qjac,pjac,q,p,LUjacQ,PjacQ,LUjacP,PjacP)
         real (kind=dblk) :: qjac(:),pjac(:),q(:),p(:),&
                             LUjacQ(:,:), LUjacP(:,:)
         integer (kind=intk) :: PjacQ(:), PjacP(:)
      
      end subroutine apply_jacobi_inv


end module coords
