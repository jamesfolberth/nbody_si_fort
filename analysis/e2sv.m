function [rv_cent, vv_cent] = e2sv(orb, params)
% Orbital elements to state vectors
% Algorithm 4.5 of Curtis

% perifocal coordinates
%fprintf(2, 'Might need to do repmat here\n');
rv_peri = (orb.h).^2/params.mu./(1+orb.e.*cos(orb.theta)).*[cos(orb.theta); sin(orb.theta); 0.0];
vv_peri = params.mu./orb.h.*[-sin(orb.theta); orb.e + cos(orb.theta); 0.0];

Q_peri_cent = zeros([3 3]);
Q_peri_cent(1:3) = [-sin(orb.Omega)*cos(orb.i)*sin(orb.omega)+cos(orb.Omega)*cos(orb.omega);
                    cos(orb.Omega)*cos(orb.i)*sin(orb.omega)+sin(orb.Omega)*cos(orb.omega);
                    sin(orb.i)*sin(orb.omega)];

Q_peri_cent(4:6) = [-sin(orb.Omega)*cos(orb.i)*cos(orb.omega)-cos(orb.Omega)*sin(orb.omega);
                    cos(orb.Omega)*cos(orb.i)*cos(orb.omega)-sin(orb.Omega)*sin(orb.omega);
                    sin(orb.i)*cos(orb.omega)];

Q_peri_cent(7:9) = [sin(orb.Omega)*sin(orb.i);
                    -cos(orb.Omega)*sin(orb.i);
                    cos(orb.i)];

rv_cent = Q_peri_cent*rv_peri;
vv_cent = Q_peri_cent*vv_peri;

end
