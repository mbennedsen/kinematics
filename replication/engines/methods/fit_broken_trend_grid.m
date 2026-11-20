function out=fit_broken_trend_grid(t,Z,Y,candidateDates)
% Full least-squares breakpoint search, using FWL to reuse a fixed design.
% Y columns are separate samples. All trend and covariate coefficients refit.
t=t(:);tRef=mean(t);s=t-tRef;n=numel(t);assert(size(Y,1)==n && size(Z,1)==n);
X0=[ones(n,1),s,Z];[Q,R]=qr(X0,0);assert(rcond(R)>1e-12,'Rank-deficient base design.');
H=max(t-candidateDates(:)',0);Hr=H-Q*(Q'*H);d=sum(Hr.^2,1)';assert(all(d>1e-12));
Yr=Y-Q*(Q'*Y);cross=Hr'*Yr;score=cross.^2./d;[~,index]=max(score,[],1);
lin=sub2ind(size(cross),index,1:size(Y,2));jump=cross(lin)./d(index)';
base=R\(Q'*(Y-H(:,index).*jump));coeff=[base(1:2,:);jump;base(3:end,:)];
out=struct('coeff',coeff,'breakDate',reshape(candidateDates(index),1,[]),'breakIndex',index, ...
 'tRef',tRef,'sse',sum(Yr.^2,1)-score(lin));
end
