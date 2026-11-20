function noise = fit_arma11_bootstrap_noise(u)
% Zero-mean stationary/invertible Gaussian ARMA(1,1) MLE on fitted residuals.
u=u(:);assert(all(isfinite(u)) && numel(u)>10);
rho=(u(1:end-1)'*u(2:end))/(u(1:end-1)'*u(1:end-1));rho=max(-.95,min(.95,rho));
starts=[rho 0;0 0;.8 -.4;.3 .3;-.3 .3];fits=nan(size(starts,1),5);
opt=optimset('Display','off','TolX',1e-9,'TolFun',1e-9,'MaxIter',2000,'MaxFunEvals',5000);
for i=1:size(starts,1)
 [z,f,flag]=fminsearch(@(z)arma11_profile_nll(u,tanh(z(1)),tanh(z(2))),atanh(starts(i,:)),opt);
 fits(i,:)=[tanh(z(:))' f flag i];
end
ok=fits(:,4)>0 & all(isfinite(fits),2) & max(abs(fits(:,1:2)),[],2)<1;
assert(any(ok),'ARMA likelihood optimization failed.');values=fits(:,3);values(~ok)=Inf;[~,j]=min(values);
phi=fits(j,1);theta=fits(j,2);[nll,s2,innov,F]=arma11_profile_nll(u,phi,theta);
assert(abs(phi)<.9999,'AR estimate too close to the stationarity boundary; inspect.');
% A near-unit MA root is stationary and simulable; flag its weak identification.
nearMABoundary=abs(theta)>.999;
noise=struct('phi',phi,'theta',theta,'innovationVariance',s2,'negativeLogLikelihood',nll, ...
 'optimizerExitFlag',fits(j,4),'multiStartFits',fits,'innovations',innov,'unitInnovationVariances',F, ...
 'initialization','exact stationary','varianceModel','constant innovation variance','mean',0,'nearMABoundary',nearMABoundary);
end
