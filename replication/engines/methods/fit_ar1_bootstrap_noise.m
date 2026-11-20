function noise=fit_ar1_bootstrap_noise(u)
% Exact stationary zero-mean Gaussian AR(1) MLE; same likelihood engine as ARMA.
u=u(:);assert(all(isfinite(u)) && numel(u)>10);
opt=optimset('Display','off','TolX',1e-10,'TolFun',1e-10,'MaxIter',1000,'MaxFunEvals',2000);
[z,nll,flag]=fminsearch(@(z)arma11_profile_nll(u,tanh(z),0),0,opt);phi=tanh(z);
assert(flag>0 && abs(phi)<.9999,'AR(1) noise fit failed or reached stationarity boundary.');
[~,s2,innov,F]=arma11_profile_nll(u,phi,0);
noise=struct('phi',phi,'theta',0,'innovationVariance',s2,'negativeLogLikelihood',nll, ...
 'optimizerExitFlag',flag,'innovations',innov,'unitInnovationVariances',F, ...
 'initialization','exact stationary','varianceModel','constant innovation variance','mean',0,'nearMABoundary',false);
end
