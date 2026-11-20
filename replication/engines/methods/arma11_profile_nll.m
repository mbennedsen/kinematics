function [nll,sigma2,innov,F] = arma11_profile_nll(u,phi,theta)
% Exact stationary Gaussian likelihood, concentrating out innovation variance.
% u_t = phi*u_(t-1) + e_t + theta*e_(t-1); Var(e_t)=sigma2.
u=u(:);n=numel(u);nll=Inf;sigma2=NaN;innov=nan(n,1);F=innov;
if abs(phi)>=1 || abs(theta)>=1 || any(~isfinite(u)),return;end
P=[(1+theta^2+2*phi*theta)/(1-phi^2),1;1,1];a=zeros(2,1);A=[phi theta;0 0];Q=ones(2);
for k=1:n
 f=P(1,1);if ~(f>0 && isfinite(f)),return;end
 v=u(k)-a(1);innov(k)=v;F(k)=f;K=P(:,1)/f;
 a=A*(a+K*v);P=A*(P-K*P(1,:))*A'+Q;P=(P+P')/2;
end
sigma2=mean(innov.^2./F);
if sigma2>0,nll=.5*(n*(log(2*pi)+1+log(sigma2))+sum(log(F)));end
end
