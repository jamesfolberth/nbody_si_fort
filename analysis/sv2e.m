function [orb] = sv2e(rv, vv, g_param)
% State vectors to elements
% Algorithm 4.2 of Curtis

orb = {};

%r = norm(rv,2,'columns');
r = norm_cols(rv,2);
%v = norm(vv,2,'columns');
v = norm_cols(vv,2);
vr = dot(rv,vv)./r;

hv = cross(rv,vv);
%orb.h = norm(hv,2,'columns');
orb.h = norm_cols(hv,2);

orb.i = acos(hv(3,:)./orb.h);

Nv = cross(repmat([0;0;1], 1, size(hv, 2)), hv);
%N = norm(Nv,2,'columns');
N = norm_cols(Nv,2);

orb.Omega = (acos(Nv(1,:)./N)).*(Nv(2,:)>=0.0) ...
           +(2*pi-acos(Nv(1,:)./N)).*(Nv(2,:)<0.0);
%if Nv(2,:) >= 0.0
%    orb.Omega = acos(Nv(1,:)./N);
%else
%    orb.Omega = 2*pi-acos(Nv(1,:)/N);
%end

ev = 1/g_param*(repmat(v.^2-(g_param)./r,3,1).*rv-repmat(r.*vr,3,1).*vv);
%orb.e = norm(ev,2,'columns');
orb.e = norm_cols(ev,2);

orb.omega = (acos(dot(Nv, ev)./(N.*orb.e))).*(ev(3,:)>=0.0) ...
           +(2*pi-acos(dot(Nv, ev)./(N.*orb.e))).*(ev(3,:)<0.0);
%if ev(3) >= 0.0
%    orb.omega = acos(dot(Nv, ev)./(N.*orb.e));
%else
%    orb.omega = 2*pi - acos(dot(Nv, ev)./(N.*orb.e));
%end

orb.theta = (acos(dot(ev, rv)./((orb.e).*r))).*(vr>=0.0) ...
           +(2*pi-acos(dot(ev, rv)./((orb.e).*r))).*(vr<0.0);
%if vr >= 0
%    orb.theta = acos(dot(ev, rv)/(orb.e*r));
%else
%    orb.theta = 2*pi - acos(dot(ev, rv)/(orb.e*r));
%end

%rp = (orb.h)^2/params.mu/(1+orb.e);
%ra = (orb.h)^2/params.mu/(1-orb.e);
%a = 0.5*(rp+ra);

end

% Hack to get norm(x,p,'columns') to work in matlab
function [x_nrm] = norm_cols(x,p)
   x_nrm = zeros([1,size(x,2)]);
   x_nrm = sqrt(sum(bsxfun(@times,x,x),1));

end
