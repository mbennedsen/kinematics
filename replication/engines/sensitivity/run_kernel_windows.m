function run_kernel_windows(d,base)
if nargin<2,base=fileparts(mfilename('fullpath'));end
root=fileparts(base);
snapshot=verify_engine_inputs(base);
input=fullfile(root,'data','Monthly_HadISST.csv');names={'NASA','NOAA','HadCRU','Berk','ERA5'};rows=table();
for d=d
 name=names{d};L=load(fullfile(root,'empirical','results','monthly',name,'kernel_B1999','CompletedFit.mat'),'fit');
 for year=[2010 2012 2014 2015 2016 2018 2020]
  o=L.fit.options;o.kernel.averageStartYear=year;o.outputDir=fullfile(base,'results_kernel_windows',name,sprintf('since_%d',year));
  fit=cached_empirical_fit(@estimate_local_quadratic,input,o,snapshot);E=fit.detail.(name);R=L.fit.detail.(name);
  assert(max(abs(E.bootstrapStats-R.bootstrapStats),[],'all')<1e-10,'Endpoint draws changed when only the contrast changed.');
  A=E.averageAcceleration;
  rows=[rows;table(string(name),year,E.actualT0,E.sample.t(end),A.estimate,A.ci_low,A.ci_high, ...
   'VariableNames',{'Dataset','RequestedT0','ActualT0','EndTime','Estimate','CILow','CIHigh'})]; %#ok<AGROW>
 end
end
writetable(rows,fullfile(base,['AlternativeWindows_Kernel_' name '.csv']));
end
