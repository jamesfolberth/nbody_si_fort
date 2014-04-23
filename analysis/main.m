
clear;

%plot_orbit('../data/orbit.h5','herr')
%plot_orbit('../data/orbit.h5','orb.h')
%plot_orbit('../data/orbit_test_control.h5','orb.h')
%plot_orbit('../data/orbit_test_small_pluto.h5','orb.h')

%plot_orbit('../data/orbit_2014_03_14_2130.h5','herr')
%plot_orbit('../data/orbit_2014_03_14_2130.h5','orb.h')

% orbit_control.tikz
plot_orbit('../data/orbit_test_control.h5','orb.h');
print('../temp_plots/figures/orbit_control.tikz','-dtikz','-S640,480');
