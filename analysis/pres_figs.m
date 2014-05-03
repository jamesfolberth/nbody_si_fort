% GNU Octave
function [] = pres_figs(fig)
   
   SIZE = '-S640,480';

   switch fig
      % relative error in Hamiltonian
      case {1,"1"}
         % {{{
         savefile = '../data/orbit_test_control.h5';
         
         load(savefile);
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
         
         figure('visible','off')
         Herrr = (H-H(1))/H(1);
         plot(t(1:50:end), Herrr(1:50:end),'-')
         
         xlabel('$t$ (yr)','Interpreter','tex');
         ylabel('Relative Error in Hamiltonian')
         print('figures/rel_err_pres.tikz','-dtikz',SIZE);
         %print('figures/rel_err.png','-dpng','-S640,480');
         % }}}

      case {2,"2"}
         % {{{
         savefile = '../data/orbit_test_control.h5';
         
         load(savefile);

         plot_planet = 1; % but Pluto isn't a planet! :(
         vars = (3*(plot_planet)+1):(3*(plot_planet)+3);
         orb = sv2e(Q(vars,:), P(vars,:)/m_vec(plot_planet+1), g_param);
         rp = (orb.h).^2./g_param./(1+orb.e);
         ra = (orb.h).^2./g_param./(1-orb.e);
         a = 0.5*(rp+ra);
         figure('visible','off');
         tdat = t(1:50:end);
         dat = orb.e.*sin(orb.omega+orb.Omega);
         %dat = orb.i*180/pi;
         dat = dat(1:50:end);
         plot(tdat,dat)
         %axis([1960 10^9+1960]);
         %axis([t(1),t(end)]);
         %title('$h=e\sin(\omega+\Omega)$','Interpreter','tex');
         xlabel('$t$ ($\mathrm{yr}$)','Interpreter','tex');
         ylabel('$h=e\sin(\omega+\Omega)$','Interpreter','tex');

         print('figures/orb_h_pres.tikz','-dtikz',SIZE);

         % }}}

      case {3,"3"}
         % {{{
         savefile = '../data/lyapunov_dt_0.1_2014_03_26_1300.h5';
         load(savefile);

         ZERO_TOL = 10^(-8);

         %ind0 = 3*10^6;
         ind0 = 1;
         ind1 = 10^7;
         
         % Isolate probably convergent data in tail of lyapunov
         tdat = t(ind0:ind1);
         ldat = lambda(ind0:ind1);
      
         % strip too small values
         bad_inds = ldat < ZERO_TOL;
         good_inds = not(bad_inds);
         %tdat = tdat(good_inds)(1:500:end);
         %ldat = ldat(good_inds)(1:500:end);
         tdat = tdat(good_inds);
         ldat = ldat(good_inds);
         tdat = cat(2,tdat(1:5:10^3),tdat(10^3+50:50:10^4),...
                  tdat(10^4+250:250:10^5),tdat(10^5+500:500:10^6),...
                  tdat(10^6+2500:2500:end));
         ldat = cat(2,ldat(1:5:10^3),ldat(10^3+50:50:10^4),...
                  ldat(10^4+250:250:10^5),ldat(10^5+500:500:10^6),...
                  ldat(10^6+2500:2500:end));

         % plot interdasting region
         figure('visible','off')
         loglog(tdat,ldat,'b.');
         axis([10^3 10^9]);
        
         % lyap_full: change ind0 to 1, unset axis
         xlabel('t ($\mathrm{yr}$)','Interpreter','tex')
         ylabel('$\mu_{max}$ ($\mathrm{yr}^{-1}$)','Interpreter','tex');
         print('figures/lyap_full_pres.tikz','-dtikz',SIZE);

         % }}}

      case {4,"4"}
         % {{{
         savefile = '../data/lyapunov_dt_0.1_2014_03_26_1300.h5';
         load(savefile);

         pkg load optim;

         ZERO_TOL = 10^(-8);

         ind0 = 3*10^6;
         ind1 = 10^7;
         
         % Isolate probably convergent data in tail of lyapunov
         tdat = t(ind0:ind1);
         ldat = lambda(ind0:ind1);
      
         % strip too small values
         bad_inds = ldat < ZERO_TOL;
         good_inds = not(bad_inds);
         tdat = tdat(good_inds)(1:500:end);
         ldat = ldat(good_inds)(1:500:end);
        
         % plot interdasting region
         figure('visible','off');
         loglog(t(10^6:250:end),lambda(10^6:250:end),'b.'); % plot more data for zoom_fit plot
         axis([10^8 5*10^9, 10^(-7.3), 10^-6.8]);
        
         % Levenberg-Marquardt to fit curve
         % TODO bad values in covp
         pin = [10^-7; 10^-9; 10^-9; tdat(1)];
         tic
         [f,p,cvg,iter,corp,covp] = leasqr(tdat,ldat,pin,@(x_,p_) fit_fcn(x_,p_));
         %p = [ 9.9916e-08;
         %      1.5883e-09;
         %     -3.6356e-09;
         %      7.9484e+08];
      
         %% Try using curvefit_stat?
         %settings = optimset();
         %settings.ret_covp = true;
         %%residmin_stat(@(x_,p_) fit_fcn(x_,p_),p,settings)
         %curvefit_stat(@(p_,x_) fit_fcn(x_,p_),p,tdat,ldat,settings)
       
         % Plot fit curve
         hold on;
         lfit = fit_fcn(tdat,p);
         loglog(tdat(1:10:end),lfit(1:10:end),'r','linewidth',3)
         axis([10^8 10^9, 10^(-7.3), 10^-6.8]);
         hold off;
      
         % lyapunov exponent is limit as t->infinity
         lyap_exp = p(1) % 1/years
         lyap_time = 1/p(1) % years
      
         xlabel('t ($\mathrm{yr}$)','Interpreter','tex')
         ylabel('$\mu_{max}$ ($\mathrm{yr}^{-1}$)','Interpreter','tex');
         legend('data','fit');
         print('figures/lyap_zoom_fit_pres.tikz','-dtikz',SIZE);

         % }}}

      otherwise
         error('you dun goofed');
   end
end

% shifted exponential decay y = p(1)+p(2)*exp(p(3)*(x-p(4)))
function [y] = fit_fcn(x,p)
   y = p(1)+p(2)*exp(p(3)*(x-p(4)));
end
