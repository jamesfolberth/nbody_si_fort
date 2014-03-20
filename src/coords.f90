! coords.f90

module coords
   
   use utils
   use orbital_data

   implicit none

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
            g_param_jac(i) = eta(i)/eta(i-1)*g_param
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
            jacT(3*i+1, 3*(i+1)+1:3*(n_masses-1)+1:3) = m_vec(i+2:n_masses)&
               /eta(i+1)
            jacT(3*i+2, 3*(i+1)+2:3*(n_masses-1)+2:3) = m_vec(i+2:n_masses)&
               /eta(i+1)
            jacT(3*i+3, 3*(i+1)+3:3*(n_masses-1)+3:3) = m_vec(i+2:n_masses)&
               /eta(i+1)
         end do

      end subroutine jacobi_setup


      ! convert to Jacobi coordinates (uses sparse mat-vec)
      pure subroutine apply_jacobi(q,p,qjac,pjac,jacQ,jacP)
         real (kind=dblk),intent(in) :: q(:), p(:),jacQ(:,:), jacP(:,:)
         real (kind=dblk),intent(out) :: qjac(:),pjac(:)

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
         qjac = 0.0_dblk
         pjac = 0.0_dblk
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
      ! Sparse routine is ~3.3 times faster than dense dgetrs
      subroutine apply_jacobi_inv(qjac,pjac,q,p,PjacQ,LUjacQ,PjacP,LUjacP)
         real (kind=dblk),intent(in):: qjac(:),pjac(:),LUjacQ(:,:), LUjacP(:,:)
         integer (kind=intk),intent(in) :: PjacQ(:), PjacP(:)
         real (kind=dblk),intent(out) :: q(:),p(:)

         !integer (kind=intk) :: info

         q = qjac
         !call dgetrs('N',3*n_masses,1,LUjacQ,3*n_masses,PjacQ,q,3*n_masses,info)
         call apply_jacobi_invqp(qjac,q,PjacQ,LUjacQ)
         
         p = pjac
         !call dgetrs('N',3*n_masses,1,LUjacP,3*n_masses,PjacP,p,3*n_masses,info)
         call apply_jacobi_invqp(pjac,p,PjacP,LUjacP)
      
      end subroutine apply_jacobi_inv

      
      ! Convert Qjac from Jacobi coordinates
      ! Sparse routine is ~3.3 times faster than dense dgetrs
      pure subroutine apply_jacobi_invqp(qjac,q,PjacQ,LUjacQ)
         real (kind=dblk), intent(in) :: qjac(:),LUjacQ(:,:)
         integer (kind=intk), intent(in) :: PjacQ(:)
         real (kind=dblk), intent(out) :: q(:)

         integer (kind=intk) :: i,j,ip
         real (kind=dblk) :: temp


         q = qjac
         ! Apply permutation vector to RHS
         !call dlaswp(3*n_masses,qjac,3*n_masses,1,3*n_masses,PjacQ,1)
         do i=1,3*n_masses
            ip = PjacQ(i)
            temp = q(i)
            q(i) = q(ip)
            q(ip) = temp
         end do

         ! Apply FW subs to system with LUjacQ/LUjacP
         ! Dense
         !do j=1,3*n_masses
         !   ! it is assumed that L(j,j) = 1

         !   do i=j+1,3*n_masses
         !      q(i) = q(i) - LUjacQ(i,j)*q(j)
         !   end do
         !end do
         ! Sparse
         do j=0,n_masses-1
            ! it is assumed that L(j,j) = 1

            do i=j+1,n_masses-1
               q(3*i+1) = q(3*i+1) - LUjacQ(3*i+1,3*j+1)*q(3*j+1)
               q(3*i+2) = q(3*i+2) - LUjacQ(3*i+2,3*j+2)*q(3*j+2)
               q(3*i+3) = q(3*i+3) - LUjacQ(3*i+3,3*j+3)*q(3*j+3)
            end do
         end do

         ! Apply BW subs
         ! Dense
         !do j=3*n_masses,1,-1
         !   ! we should never have a zero on the diagonal of LUjacQ/LUjacP

         !   q(j) = q(j)/LUjacQ(j,j)
         !   do i=j-1,1,-1
         !      q(i) = q(i) - LUjacQ(i,j)*q(j)
         !   end do
         !end do
         ! Sparse
         do j=n_masses-1,0,-1
            ! we should never have zeros on the diagonal

            q(3*j+1) = q(3*j+1) / LUjacQ(3*j+1,3*j+1)
            q(3*j+2) = q(3*j+2) / LUjacQ(3*j+2,3*j+2)
            q(3*j+3) = q(3*j+3) / LUjacQ(3*j+3,3*j+3)

            do i=j-1,0,-1
               q(3*i+1) = q(3*i+1) - LUjacQ(3*i+1,3*j+1)*q(3*j+1)
               q(3*i+2) = q(3*i+2) - LUjacQ(3*i+2,3*j+2)*q(3*j+2)
               q(3*i+3) = q(3*i+3) - LUjacQ(3*i+3,3*j+3)*q(3*j+3)
            end do
         end do
      
      end subroutine apply_jacobi_invqp


      
      ! convert Q or P to Jacobi coordinates (uses sparse mat-vec)
      ! used in symplectic.f90
      pure subroutine apply_jacobiqp(p,pjac,jacP)
         real (kind=dblk), intent(in) :: p(:),jacP(:,:)
         real (kind=dblk), intent(out) :: pjac(:)

         integer (kind=intk) :: i,j
         ! pjac <- jacP*p
         ! Note that jacQ and jacP have the same sparsity pattern
         ! This is about 3 times faster than a dense mat-vec
         do j = 1,3*n_masses,3
            ! multiply against first three rows of jacP
            pjac(1) = pjac(1) + jacP(1,j)*p(j)
            pjac(2) = pjac(2) + jacP(2,j+1)*p(j+1)
            pjac(3) = pjac(3) + jacP(3,j+2)*p(j+2)
         end do

         do i=4,3*n_masses,3
            pjac(i) = pjac(i) + jacP(i,1)*p(1)
            pjac(i+1) = pjac(i+1) + jacP(i+1,2)*p(2)
            pjac(i+2) = pjac(i+2) + jacP(i+2,3)*p(3)
         end do
 
         do j=4,3*n_masses,3
            do i = j,3*n_masses,3
               pjac(i) = pjac(i) + jacP(i,j)*p(j)
               pjac(i+1) = pjac(i+1) + jacP(i+1,j+1)*p(j+1)
               pjac(i+2) = pjac(i+2) + jacP(i+2,j+2)*p(j+2)
            end do
         end do

      end subroutine apply_jacobiqp


      ! Apply qt = jacT*q
      ! jacT is a sum from Saha 1994, represented as mat-vec
      ! This is essentially sparse mat-vec
      pure subroutine apply_jacT(q,qt,jacT)
         real (kind=dblk), intent(in) :: q(:),jacT(:,:)
         real (kind=dblk), intent(out) :: qt(:)

         integer (kind=intk) :: i,j

         qt = 0.0_dblk
         qt(1:3) = q(1:3)
         do j=4,3*n_masses,3
            do i=4,j,3
               qt(i) = qt(i) + jacT(i,j)*q(j)
               qt(i+1) = qt(i+1) + jacT(i+1,j+1)*q(j+1)
               qt(i+2) = qt(i+2) + jacT(i+2,j+2)*q(j+2)
            end do
         end do

      end subroutine apply_jacT


end module coords
