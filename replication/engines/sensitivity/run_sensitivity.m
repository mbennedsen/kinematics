function run_sensitivity(d,base)
if nargin<2,base=fileparts(mfilename('fullpath'));end
root=fileparts(base);
snapshot=verify_engine_inputs(base);
names={'NASA','NOAA','HadCRU','Berk','ERA5'};name=names{d};
input=fullfile(root,'data','Monthly_HadISST.csv');
L=load(fullfile(root,'empirical','results','monthly',name,'ss_B1999','CompletedFit.mat'),'fit');
E=L.fit.detail.(name);ref=L.fit.options;t=E.sampleData.t;T0=t(end-120);assert(abs(t(end)-T0-10)<1e-8);
R=table();out=fullfile(base,'results_robustness',name);
for mult=[.5 .75 1 1.5 2]
 o=ref;o.bootstrap.B=0;o.bootstrap.useParallel=false;o.bootstrap.nWorkers=1;o.bootstrap.pathBands=false;
 o.stateSpace.fixedQ=E.pointEstimate.paramsHat.q*mult^2;
 o.outputDir=fullfile(out,sprintf('ss_sigma_%g',mult));
 fit=cached_empirical_fit(@estimate_monthly_state_space,input,o,snapshot);
 G=fit.detail.(name);assert(G.pointEstimate.exitflag>0);
 row=empirical_summary_rows(name,'ss',G,0);row.Variant=repmat("sigma_"+string(mult),4,1);R=[R;row]; %#ok<AGROW>
end
for h=15:5:35
 o=struct('datasets',{{name}},'seed',550000+d,'rebaseData',true, ...
 'outputDir',fullfile(out,sprintf('kernel_h%d',h)), ...
 'kernel',struct('bandwidthYears',h,'averageStartYear',T0,'trailingMonths',60), ...
 'bootstrap',struct('B',0,'ciAlpha',.1,'useParallel',false,'verboseEvery',0));
 fit=cached_empirical_fit(@estimate_local_quadratic,input,o,snapshot);
 row=empirical_summary_rows(name,'kernel',fit.detail.(name),0);row.Variant=repmat("bandwidth_"+string(h),4,1);R=[R;row]; %#ok<AGROW>
end
for year=[1950 1960 1970 1980 1990]
 for kind={'quadratic','changepoint'}
  method=kind{1};o=struct('datasets',{{name}},'rebaseData',true,'inference',struct('ciLevel',.9), ...
   'outputDir',fullfile(out,sprintf('%s_start%d',method,year)));
  if strcmp(method,'quadratic'),o.quadratic=struct('windowStartYear',year,'hacLagMonths',24);fun=@estimate_quadratic_benchmark;
  else,o.broken=struct('windowStartYear',year,'hacLagMonths',24,'minSegmentMonths',60,'candidateStepMonths',1);fun=@estimate_broken_trend_benchmark;end
  fit=cached_empirical_fit(fun,input,o,snapshot);row=empirical_summary_rows(name,method,fit.detail.(name),0);
  row.Variant=repmat("start_"+string(year),4,1);R=[R;row]; %#ok<AGROW>
 end
end
writetable(R,fullfile(out,'Summary.csv'));
o=struct('datasets',{{name}},'rebaseData',true,'outputDir',fullfile(out,'velocity_rw'), ...
 'covariates',ref.covariates,'figure',struct('visible','off'));
cached_empirical_fit(@estimate_velocity_rw_robustness,input,o,snapshot);
% Alternative windows: point contrasts from the same fitted state-space path.
z=E.pointEstimate.smoother.alphaSmooth;A=table();
for year=[2010 2012 2014 2015 2016 2018 2020]
 i=find(t>=year,1);a=100*(z(2,end)-z(2,i))/(t(end)-t(i));
 A=[A;table(string(name),year,t(i),t(end),a,'VariableNames',{'Dataset','RequestedT0','ActualT0','EndTime','Estimate'})]; %#ok<AGROW>
end
writetable(A,fullfile(out,'AlternativeWindows_SS_Points.csv'));
end
