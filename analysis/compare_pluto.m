% GNU Octave

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
%plot(t(1:10:end), (orb.i*180/pi)(1:10:end),'b')
%plot(t, orb.h)
%plot(t, orb.e)
%plot(t(1:50:end), (orb.i*180/pi)(1:50:end),'b')
plot(t(1:50:end),(orb.e.*sin(orb.omega+orb.Omega))(1:50:end),'b')
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
%plot(t(1:50:end), (orb.i*180/pi)(1:50:end),'r')
plot(t(1:50:end),(orb.e.*sin(orb.omega+orb.Omega))(1:50:end),'r')
hold off;

xlabel('$t$ ($\mathrm{yr}$)','Interpreter','tex');
%ylabel('$h$ ()','Interpreter','tex');
%ylabel('$i$ ()','Interpreter','tex');
%print('../temp_plots/figures/comp_pluto_h.tikz','-dtikz','-S640,480');
%print('../temp_plots/figures/comp_pluto_i.tikz','-dtikz','-S640,480');




end 
