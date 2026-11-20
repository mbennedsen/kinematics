function collect_empirical(engine,dest,kind,name,methods,B)
% Export fresh native fits to exactly the numerical interface used by reporting.
if nargin<6,B=1999;end
out=fullfile(dest,'empirical');paths=fullfile(out,'reporting');ssdir=fullfile(out,'validated_monthly');
for d={out,paths,ssdir},if ~isfolder(d{1}),mkdir(d{1});end,end
for mm=1:numel(methods)
 m=methods{mm};folder=fullfile(engine,'results',kind,name,sprintf('%s_B%d',m,B));
 R=readtable(fullfile(folder,'Summary.csv'),'TextType','string');
 merge_rows(fullfile(out,'Empirical_Results.csv'),R,{'Dataset','Method','Specification'});
 if strcmp(kind,'annual'),merge_rows(fullfile(out,'Annual_Results.csv'),R,{'Dataset','Method'});end
 sub='';if strcmp(m,'changepoint')||(strcmp(m,'ss')&&strcmp(kind,'annual')),sub='point';end
 L=load(fullfile(folder,sub,'CompletedFit.mat'),'fit');E=L.fit.detail.(name);
 if strcmp(m,'changepoint')
  Lb=load(fullfile(folder,'Bootstrap.mat'),'result');z=Lb.result;
  q=prctile(z.breakDates,[5 95]);T0=R.ActualT0(1);T=R.EndTime(1);
  diag=table(string(kind),string(name),z.validDraws,z.noise.phi,z.noise.theta,sqrt(z.noise.innovationVariance), ...
   z.noise.optimizerExitFlag,q(1),q(2),mean(z.breakDates>T0&z.breakDates<T), ...
   'VariableNames',{'Specification','Dataset','ValidDraws','AR1','MA1','InnovationSD','NoiseExitFlag','BreakDateLow','BreakDateHigh','BreakInsideDecadeFraction'});
  diag.RateIntervalWidthRatio=diff(z.intervals(2,:))/diff(z.fixedBreakIntervals(2,:));
  merge_rows(fullfile(out,'Changepoint_Diagnostics.csv'),diag,{'Specification','Dataset'});
 end
 if ~strcmp(kind,'monthly'),continue;end
 if strcmp(m,'ss')
  t=E.sampleData.t(:);y=E.sampleData.y(:);v=kinematic_path_values(E.pointEstimate.smoother.alphaSmooth);
  lo=E.bootstrap.pathBands.low;hi=E.bootstrap.pathBands.high;
  Z=array2table([t y v lo(:,1) hi(:,1) lo(:,2) hi(:,2) lo(:,3) hi(:,3)], ...
   'VariableNames',{'t','y','x','v','a','average_a','xLow','xHigh','vLow','vHigh','aLow','aHigh'});
  writetable(Z,fullfile(ssdir,[name '_ss_paths.csv']));
 else
  t=E.sample.t(:);
  if strcmp(m,'kernel')
   h=L.fit.options.kernel.bandwidthYears;y=E.yAdj(:);v=zeros(numel(t),3);
   for k=1:numel(t)
    u=(t-t(k))/h;ix=abs(u)<1;u=u(ix);w=sqrt((1-abs(u).^3).^3);X=[ones(size(u)),u,u.^2];b=(X.*w)\(y(ix).*w);
    v(k,:)=[b(1),10*b(2)/h,200*b(3)/h^2];
   end
   assert(max(abs(v(:,1)-E.xFull(:)))<1e-8,'Kernel path mismatch.');
  else
   b=E.coeff;
   if strcmp(m,'quadratic'),tt=t-mean(t);v=[b(1)+b(2)*tt+b(3)*tt.^2,10*(b(2)+2*b(3)*tt),repmat(200*b(3),numel(t),1)];
   else,tt=t-E.tRef;v=[b(1)+b(2)*tt+b(3)*max(t-E.breakDate,0),10*(b(2)+b(3)*(t>=E.breakDate)),zeros(numel(t),1)];end
  end
  for k=1:3,objects=["x","v","a"];assert(abs(v(end,k)-R.Estimate(R.Object==objects(k)))<1e-8);end
  writetable(array2table([t v],'VariableNames',{'t','x','v','a'}),fullfile(paths,[name '_' m '_paths.csv']));
 end
end
end
