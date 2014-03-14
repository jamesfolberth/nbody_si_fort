! kepler.f90

module kepler
   
   use utils
   use orbital_data

   contains

   ! Initial value problem I from Danby
   subroutine kepler_step(Qjac_wrk, Pjac_wrk,tau,r,v,u,a,n,EC,ES,&
                          e,dE,dtv,C,S,f,g,aor,fp,gp,m_vec_jac,g_param_jac)
      real (kind=dblk) :: Qjac_wrk(:), Pjac_wrk(:),tau, &
                          r(:), v(:), u(:), a(:), n(:), EC(:), ES(:),&
                          e(:), dE(:), dtv(:), C(:), S(:),&
                          f(:),g(:),aor(:),fp(:),gp(:),m_vec_jac(:),&
                          g_param_jac(:)

      integer (kind=intk) :: i
      real (kind=dblk) :: bigR(3), bigV(3), temp(3*n_masses)

      r = 0.0_dblk
      v = 0.0_dblk
      u = 0.0_dblk
      do i=1,n_masses-1
         bigV = Pjac_wrk(3*i+1:3*i+3)/m_vec_jac(i+1)
         bigR = Qjac_wrk(3*i+1:3*i+3)
         v(i+1) = norm2(bigV)
         r(i+1) = norm2(bigR)
         u(i+1) = bigR(1)*bigV(1)+bigR(2)*bigV(2)+bigR(3)*bigV(3)
      end do

      a = 1.0_dblk / (2.0_dblk / r - v*v / g_param_jac)

      if (debug) then
         if (any(a(2:n_masses) <= 0)) then
            print *, "error: kepler.f90: kepler_step: at least one semi-major axis went negative"
            stop
         end if
      end if

      n = sqrt(g_param_jac / (a * a * a))

      EC = 1.0_dblk - r/a
      ES = u / (n*a*a)
      e = sqrt(EC*EC+ES*ES)

      dE = 0.0_dblk
      dtv = 0.0_dblk
      call kepeq2(dE, dtv, tau, e, n, EC, ES)
      
      C = dcos(dE); S = dsin(dE);

      ! Lagrange f,g functions
      f = a / r * (C-1.0_dblk) + 1.0_dblk
      g = dtv + 1.0_dblk / n * (S-dE)
      aor = 1.0_dblk / (1.0_dblk - EC*C+ES*S)
      fp = -aor*a/r*n*S
      gp = (C-1.0_dblk)*aor+1.0_dblk

      temp = Qjac_wrk
      do i=1,n_masses-1
         Qjac_wrk(3*i+1:3*i+3) = f(i+1)*Qjac_wrk(3*i+1:3*i+3)+g(i+1)*&
            Pjac_wrk(3*i+1:3*i+3)/m_vec_jac(i+1)
         Pjac_wrk(3*i+1:3*i+3) = m_vec_jac(i+1)*(fp(i+1)*temp(3*i+1:3*i+3))&
            +gp(i+1)*Pjac_wrk(3*i+1:3*i+3)
      end do


   end subroutine kepler_step


   ! Solve Kepler's equation for IVP I in Danby
   subroutine kepeq2(dE, dtv, tau, e, n, EC, ES)
      real (kind=dblk) :: dE(:), dtv(:), tau, e(:), n(:), EC(:), ES(:)

      integer (kind=dblk) :: i,counter
      real (kind=dblk) :: ndt,y,sigma,dx,c,s,f,fp,fpp,fppp

      do i=2,n_masses ! Don't do the sun
         dtv(i) = tau - floor(n(i)*tau/(2*pi))*2*pi/n(i)
         ndt = n(i)*dtv(i)

         if (e(i) < 0.1_dblk) then
            dE(i) = ndt
         else
            y = ndt - ES(i)
            sigma = sign(1.0_dblk,ES(i)*dcos(y)+EC(i)*dsin(y))
            dE(i) = y + 0.85_dblk*sigma*e(i) - ES(i)
         end if

         dx = 1.0_dblk
         if (debug) counter = 0
         do while (dabs(dx) > 10.0E-14_dblk .and. counter < 10)
            c = dcos(dE(i)); s = dsin(dE(i));
            ! Halley-type quartic method; Danby 6.6.7
            f = dE(i) - EC(i)*s+ES(i)*(1.0_dblk-c)-ndt
            fp = 1.0_dblk - EC(i)*c+ES(i)*s
            fpp = EC(i)*s+ES(i)*c
            fppp = EC(i)*c-ES(i)*s
            dx = -f/fp
            dx = -f/(fp+dx*fpp/2.0_dblk)
            dx = -f/(fp+dx*fpp/2.0_dblk+dx*dx*fppp/6.0_dblk)
            dE(i) = dE(i) + dx
            if (debug) counter = counter + 1
         end do

         if (debug) then
            if (counter >= 10) then 
               print *, "error: kepler.f90: kepeq2: Took too long in kepeq2"
               stop
            end if
         end if

         dE(i) = mod(dE(i), 2*pi)

      end do

   end subroutine kepeq2


end module kepler
