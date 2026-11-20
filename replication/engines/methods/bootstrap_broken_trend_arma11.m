function result=bootstrap_broken_trend_arma11(E,opts)
% Parametric percentile bootstrap of the entire existing OLS break search.
% Keeps the empirical design/covariates fixed; refits ALL regression coefficients
% and reselects the breakpoint in every draw. ARMA is a bootstrap noise model,
% not a change to the OLS point estimator. No null-break test is performed.
t=E.sample.t(:);Z=E.sample.Z;y=E.sample.y(:);n=numel(t);B=opts.B;
assert(B>=19 && all(isfinite([t y]),'all') && all(isfinite(Z),'all'));
point=fit_broken_trend_grid(t,Z,y,E.candidateBreakDates);
assert(abs(point.breakDate-E.breakDate)<1e-8 && max(abs(point.coeff-E.coeff))<1e-8,'Point estimator mismatch.');
s=t-E.tRef;X=[ones(n,1),s,max(t-E.breakDate,0),Z];mu=X*E.coeff;u=y-mu;
assert(max(abs(u-E.resid))<1e-8);
if strcmp(opts.noiseModel,'AR1'),noise=fit_ar1_bootstrap_noise(u);
elseif strcmp(opts.noiseModel,'ARMA11'),noise=fit_arma11_bootstrap_noise(u);
else,error('Unsupported bootstrap noise model');end
seeds=opts.seed+(1:B);U=simulate_arma11_stationary(noise,n,seeds);Y=mu+U;
fitted=fit_broken_trend_grid(t,Z,Y,E.candidateBreakDates);b=fitted.coeff;tau=fitted.breakDate;
% Independent reference enumerations guard the optimized projection/search.
for k=unique([1,min(7,B),B])
 best=Inf;bc=[];tc=NaN;
 for j=1:numel(E.candidateBreakDates)
  tj=E.candidateBreakDates(j);XX=[ones(n,1),s,max(t-tj,0),Z];bb=XX\Y(:,k);ss=sum((Y(:,k)-XX*bb).^2);
  if ss<best,best=ss;bc=bb;tc=tj;end
 end
 assert(abs(tc-tau(k))<1e-8 && max(abs(bc-b(:,k)))<1e-8,'Search disagrees with independent enumeration.');
end
est=[E.position.endpoint,E.velocity.endpoint,0];
draws=[b(1,:)+b(2,:)*s(end)+b(3,:).*(t(end)-tau);10*(b(2,:)+b(3,:));zeros(1,B)]';
% Use right-hand slopes at a break: a break exactly at T0 is excluded.
windowT0=opts.windowT0(:);windowNames=string(opts.windowNames(:));nw=numel(windowT0);
assert(nw==numel(windowNames) && all(windowT0<t(end)) && all(windowT0>=t(1)));
wd=zeros(B,nw);we=zeros(nw,1);
for j=1:nw
 t0=windowT0(j);wd(:,j)=(100*b(3,:).*(double(t(end)>=tau)-double(t0>=tau))/(t(end)-t0))';
 we(j)=100*E.coeff(3)*(double(t(end)>=E.breakDate)-double(t0>=E.breakDate))/(t(end)-t0);
end
pcts=100*[opts.alpha/2,1-opts.alpha/2];ci=prctile([draws,wd(:,1)],pcts,1)';ci(3,:)=NaN;
W=table(windowNames,windowT0,repmat(t(end),nw,1),we,prctile(wd,pcts(1),1)',prctile(wd,pcts(2),1)', ...
 mean(tau>windowT0,2),repmat(B,nw,1),'VariableNames', ...
 {'Window','ActualT0','EndTime','Estimate','CILow','CIHigh','BreakAfterStartFraction','ValidDraws'});
% Paired fixed-break bootstrap, for diagnosing the effect of selection itself.
bfixed=X\Y;fixed=[bfixed(1,:)+bfixed(2,:)*s(end)+bfixed(3,:)*(t(end)-E.breakDate);10*(bfixed(2,:)+bfixed(3,:))]';
result=struct('options',opts,'point',E,'noise',noise,'estimates',[est we(1)], ...
 'intervals',ci,'draws',[draws,wd(:,1)],'allWindowDraws',wd,'windows',W, ...
 'coefficients',b,'breakDates',tau,'seeds',seeds,'fixedBreakDraws',fixed, ...
 'fixedBreakIntervals',prctile(fixed,pcts,1)','slopeChangeDraws',10*b(3,:)', ...
 'slopeChangeInterval',prctile(10*b(3,:),pcts),'breakDateInterval',prctile(tau,pcts), ...
 'validDraws',B,'referenceSearchChecks',numel(unique([1,min(7,B),B])));
assert(all(isfinite(result.draws),'all') && all(isfinite(b),'all'));
end
