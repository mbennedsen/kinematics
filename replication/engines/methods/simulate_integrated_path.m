function S = simulate_integrated_path(C,d,rep)
% Integrate the zero, ramp or smoothed NASA acceleration path.
validateattributes(d,{'numeric'},{'scalar','integer'});
assert(ismember(d,[1 2 4]),'Unknown integrated-path design.');
validateattributes(rep,{'numeric'},{'scalar','integer','>=',1,'<=',C.cfg.nRep});
s = C.sample; p = C.est.paramsHat; cfg = C.cfg;
n = s.T; dt = 1/12; u = s.t-s.t(1); L = u(end);
alpha = zeros(3,n); alpha(:,1) = [C.x0;C.v0;C.aNASA(1)];
switch d
    case 1
        a = zeros(n,1);
    case 2
        a = cfg.rampEnd*u/L;
    case 4
        a = C.aNASA;
end
alpha(3,:) = a';
% Piecewise-linear interpolation of a gives exact cubic x and quadratic v.
for k = 2:n
    da = a(k)-a(k-1);
    alpha(2,k) = alpha(2,k-1)+dt*(a(k-1)+da/2);
    alpha(1,k) = alpha(1,k-1)+dt*alpha(2,k-1)+dt^2*(a(k-1)/2+da/6);
end
% Common observation-noise draws across DGPs, independent of latent shocks.
% Burn in ARMA with the first era and the recurring monthly variance pattern.
rng(cfg.seed+1000000+rep,'twister');
e = 0; lastU = 0; noise = zeros(n,1);
for k = (1-cfg.burninMonths):n
    if k < 1
        month = mod(s.month(1)+k-2,12)+1; era = s.era(1);
    else
        month = s.month(k); era = s.era(k);
    end
    scale = 1;
    if C.est.setup.useMonthlyVariance, scale = scale*exp(p.gamma(month)); end
    if numel(p.delta) >= era, scale = scale*exp(p.delta(era)); end
    innovation = sqrt(max(p.sig2u*scale,C.est.setup.minVariance))*randn;
    e = p.ar(1)*e + innovation + p.maCoeffs(1)*lastU;
    lastU = innovation;
    if k >= 1, noise(k) = e; end
end
S = struct('alpha',alpha,'noise',noise);
S.y = alpha(1,:)' + s.W*p.beta + noise;
S.truth = [alpha(1,end),10*alpha(2,end),100*alpha(3,end), ...
    100*(alpha(2,end)-alpha(2,C.index0))/(s.t(end)-C.actualT0)];
end
