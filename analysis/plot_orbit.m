
function [] = plot_orbit(savefile,plot_type)

if exist(savefile)
   [t,Q,P,Qjac,Pjac,jacQ,jacP,jacT,PjacQ,LUjacQ,PjacP,LUjacP,m_vec,m_vec_jac,g_const,g_param,g_param_jac] = load_orbit_data(savefile);
else
   error(sprintf('plot_stuff_orbit.m: data file %s not found',savefile));
end


if strcmp(plot_type,'herr')
   %% Hamiltonian %%
   H = zeros(size(Pjac(1,:)));
   for plot_planet=1:numel(m_vec)-1
      vars = (3*(plot_planet)+1):(3*(plot_planet)+3);
      
      % Kepler part
      H = H + sum(Pjac(vars,:).*Pjac(vars,:))/(2*m_vec_jac(plot_planet+1));
      H = H - g_param*m_vec(plot_planet+1)./sqrt(sum((Qjac(vars,:)).^2));
   
      % Interaction terms
      H = H + g_const*m_vec(1)*m_vec(plot_planet+1).*(1./sqrt(sum(Qjac(vars,:).^2))-1./sqrt(sum((Q(vars,:)-Q(1:3,:)).^2)));
      
      for j=plot_planet+1:numel(m_vec)-1
         H = H - g_const*m_vec(plot_planet+1)*m_vec(j+1)./sqrt(sum((Q(vars,:)-Q(3*j+1:3*j+3,:)).^2));
      end
   end 
   
   plot_planet = 0;
   vars = (3*(plot_planet)+1):(3*(plot_planet)+3);
   H = H + sum(Pjac(vars,:).*Pjac(vars,:))/(2*m_vec_jac(plot_planet+1));
   
   figure()
   Herrr = (H-H(1))/H(1);
   plot(t(1:end), Herrr(1:end))
end

end 
