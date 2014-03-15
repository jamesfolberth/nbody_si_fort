function [t,Q,P,Qjac,Pjac,jacQ,jacP,jacT,PjacQ,LUjacQ,PjacP,LUjacP,m_vec,m_vec_jac,g_const,g_param,g_param_jac] = load_orbit_data(savefile)
% load the HDF5 file 'savefile', read the data, and output

if exist(savefile)
   t = h5read(savefile, '/t');
   Q = h5read(savefile, '/Q');
   P = h5read(savefile, '/P');
   Qjac = h5read(savefile, '/Qjac');
   Pjac = h5read(savefile, '/Pjac');
   jacQ = h5read(savefile, '/jacQ');
   jacP = h5read(savefile, '/jacP');
   jacT = h5read(savefile, '/jacT');
   PjacQ = h5read(savefile, '/PjacQ');
   LUjacQ = h5read(savefile, '/LUjacQ');
   PjacP = h5read(savefile, '/PjacP');
   LUjacP = h5read(savefile, '/LUjacQ');
   m_vec = h5read(savefile, '/m_vec');
   m_vec_jac = h5read(savefile, '/m_vec_jac');
   g_const = h5read(savefile, '/g_const');
   g_param = h5read(savefile, '/g_param');
   g_param_jac = h5read(savefile, '/g_param_jac');

else
   error(sprintf('load_orbit_data.m: cannot load HDF5 file: %s', savefile));
end

end 
