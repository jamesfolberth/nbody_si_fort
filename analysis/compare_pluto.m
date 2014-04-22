
function [] = compare_pluto()

close all;

%plot_orbit('../data/orbit_test_control.h5', 'orb.h');
load('../data/orbit_test_control.h5')

plot_planet = 1;
vars = (3*(plot_planet)+1):(3*(plot_planet)+3);
orb = sv2e(Q(vars,:), P(vars,:)/m_vec(plot_planet+1), g_param);
rp = (orb.h).^2./g_param./(1+orb.e);
ra = (orb.h).^2./g_param./(1-orb.e);
a = 0.5*(rp+ra);
figure();
%plot(t, Q(vars(1), :));
%plot(t, a)
plot(t, orb.i*180/pi,'b.')
%plot(t, orb.h)
%plot(t, orb.e)
%plot(t,orb.e.*sin(orb.omega+orb.Omega),'b.')
axis([1960 10^9+1960]);


%plot_orbit('../data/orbit_test_small_pluto.h5', 'orb.h');
%load('../data/orbit_test_small_pluto.h5')
load('../data/orbit_test_small_pluto_10e6.h5')

plot_planet = 1;
vars = (3*(plot_planet)+1):(3*(plot_planet)+3);
orb = sv2e(Q(vars,:), P(vars,:)/m_vec(plot_planet+1), g_param);
rp = (orb.h).^2./g_param./(1+orb.e);
ra = (orb.h).^2./g_param./(1-orb.e);
a = 0.5*(rp+ra);
hold on;
plot(t, orb.i*180/pi,'r.')
%plot(t,orb.e.*sin(orb.omega+orb.Omega),'r.')
hold off;

end 
