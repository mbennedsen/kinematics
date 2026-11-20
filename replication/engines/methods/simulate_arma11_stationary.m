function U=simulate_arma11_stationary(noise,n,seeds)
% Exact stationary starts; independent reproducible stream for each replication.
phi=noise.phi;theta=noise.theta;s2=noise.innovationVariance;
P=s2*[(1+theta^2+2*phi*theta)/(1-phi^2),1;1,1];L=chol(P,'lower');
U=zeros(n,numel(seeds));
for b=1:numel(seeds)
 stream=RandStream('mt19937ar','Seed',seeds(b));z=randn(stream,n+2,1);initial=L*z(1:2);e=sqrt(s2)*z(3:end);
 last=initial(1);lastE=initial(2);
 for k=1:n,last=phi*last+theta*lastE+e(k);U(k,b)=last;lastE=e(k);end
end
end
