
function [] = lyap()
   
   pkg load optim % GNU Octave optimization package ('leasqr')

   ZERO_TOL = 10^(-8);

   load('../data/lyapunov_dt_0.1_2014_03_26_1300.h5');
   
   % look towards the end, near index 3*10^6
   % loglog(1:10^7,lambda,'.');
   
   ind0 = 3*10^6;
   %ind0 = 1;
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
   %loglog(tdat,ldat,'b');
   loglog(t(10^6:500:end),lambda(10^6:500:end),'b'); % plot more data for zoom_fit plot
   axis([10^8 5*10^9, 10^(-7.3), 10^-6.8]);
  
   % lyap_full: change ind0 to 1, unset axis
   %xlabel('t ($\mathrm{yr}$)')
   %ylabel('$\lambda_1$ ($\mathrm{yr}^{-1}$)','Interpreter','tex');
   %print('../temp_plots/figures/lyap_full.tikz','-dtikz','-S640,480');
   %error('stuff')

   % Levenberg-Marquardt to fit curve
   % TODO bad values in covp
   pin = [10^-7; 10^-9; 10^-9; tdat(1)];
   tic
   [f,p,cvg,iter,corp,covp] = leasqr(tdat,ldat,pin,@(x_,p_) fit_fcn(x_,p_));
   toc
   p
   covp
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
   loglog(tdat(1:10:end),lfit(1:10:end),'r')
   axis([10^8 5*10^9, 10^(-7.3), 10^-6.8]);
   hold off;

   % lyapunov exponent is limit as t->infinity
   lyap_exp = p(1)
   lyap_time = 1/p(1) % years

   xlabel('t ($\mathrm{yr}$)')
   ylabel('$\lambda_1$ ($\mathrm{yr}^{-1}$)','Interpreter','tex');
   print('../temp_plots/figures/lyap_zoom_fit.tikz','-dtikz','-S640,480');
 

end


% shifted exponential decay y = p(1)+p(2)*exp(p(3)*(x-p(4)))
function [y] = fit_fcn(x,p)
   y = p(1)+p(2)*exp(p(3)*(x-p(4)));
end
