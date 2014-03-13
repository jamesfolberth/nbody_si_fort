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
            jacQ(3*i+3,3*i+3) = 1d0
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


      ! convert to Jacobi coordinates (uses sparse mat-vec)
      subroutine apply_jacobi(q,p,qjac,pjac,jacQ,jacP)
         real (kind=dblk) :: q(:), p(:), qjac(:), pjac(:),&
                             jacQ(:,:), jacP(:,:)

         integer (kind=intk) :: i,j
         !real (kind=dblk) :: qjt(18),pjt(18)

         !do j=1,3*n_masses
         !   do i=1,3*n_masses
         !      qjac(i) = qjac(i) + jacQ(i,j)*q(j);
         !      pjac(i) = pjac(i) + jacP(i,j)*p(j);
         !   end do
         !end do

         !call print_vector(qjac)
         !call print_vector(pjac)

         !qjt = qjac
         !pjt = pjac
         !qjac = 0
         !pjac = 0

         ! qjac <- jacQ*q
         ! pjac <- jacP*p
         ! Note that jacQ and jacP have the same sparsity pattern
         ! This is about 3 times faster than the above (dense) mat-vec
         do j = 1,3*n_masses,3
            ! multiply against first three rows of jacQ, jacP
            qjac(1) = qjac(1) + jacQ(1,j)*q(j)
            qjac(2) = qjac(2) + jacQ(2,j+1)*q(j+1)
            qjac(3) = qjac(3) + jacQ(3,j+2)*q(j+2)

            pjac(1) = pjac(1) + jacP(1,j)*p(j)
            pjac(2) = pjac(2) + jacP(2,j+1)*p(j+1)
            pjac(3) = pjac(3) + jacP(3,j+2)*p(j+2)
         end do

         do i=4,3*n_masses,3
            qjac(i) = qjac(i) + jacQ(i,1)*q(1)
            qjac(i+1) = qjac(i+1) + jacQ(i+1,2)*q(2)
            qjac(i+2) = qjac(i+2) + jacQ(i+2,3)*q(3)

            pjac(i) = pjac(i) + jacP(i,1)*p(1)
            pjac(i+1) = pjac(i+1) + jacP(i+1,2)*p(2)
            pjac(i+2) = pjac(i+2) + jacP(i+2,3)*p(3)
         end do
 
         do j=4,3*n_masses,3
            do i = j,3*n_masses,3
               qjac(i) = qjac(i) + jacQ(i,j)*q(j)
               qjac(i+1) = qjac(i+1) + jacQ(i+1,j+1)*q(j+1)
               qjac(i+2) = qjac(i+2) + jacQ(i+2,j+2)*q(j+2)

               pjac(i) = pjac(i) + jacP(i,j)*p(j)
               pjac(i+1) = pjac(i+1) + jacP(i+1,j+1)*p(j+1)
               pjac(i+2) = pjac(i+2) + jacP(i+2,j+2)*p(j+2)
            end do
         end do

         !call print_vector(qjac)
         !call print_vector(pjac)
         !print *, norm2(qjac-qjt)
         !print *, norm2(pjac-pjt)

      end subroutine apply_jacobi


      ! Convert from Jacobi coordinates
      subroutine apply_jacobi_inv(qjac,pjac,q,p,PjacQ,LUjacQ,PjacP,LUjacP)
         real (kind=dblk) :: qjac(:),pjac(:),q(:),p(:),&
                             LUjacQ(:,:), LUjacP(:,:)
         integer (kind=intk) :: PjacQ(:), PjacP(:)

         ! TODO write sparse FW/BW subs routine

      
      end subroutine apply_jacobi_inv



end module coords
