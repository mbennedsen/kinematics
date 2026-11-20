function collect_robustness(rob,emp,dest)
out=fullfile(dest,'empirical');S=table();C=table();V=table();names={'NASA','NOAA','HadCRU','Berk','ERA5'};
for d=1:5
 name=names{d};folder=fullfile(rob,'results_robustness',name);R=readtable(fullfile(folder,'Summary.csv'),'TextType','string');S=[S;R]; %#ok<AGROW>
 for year=[1950 1960 1970 1980 1990]
  L=load(fullfile(folder,sprintf('changepoint_start%d',year),'CompletedFit.mat'),'fit');E=L.fit.detail.(name);
  Z=table(repmat(string(name),2,1),repmat("changepoint",2,1),["slope";"break"],[10*E.coeff(3);E.breakDate],repmat("start_"+year,2,1), ...
   'VariableNames',{'Dataset','Method','Object','Estimate','Variant'});C=[C;Z]; %#ok<AGROW>
 end
 L=load(fullfile(emp,'results','monthly',name,'ss_B1999','CompletedFit.mat'),'fit');base=L.fit.detail.(name).pointEstimate.smoother.alphaSmooth;
 L=load(fullfile(folder,'velocity_rw','CompletedFit.mat'),'fit');E=L.fit.detail.(name);a=E.estimate.smoother.alphaSmooth;
 assert(E.estimate.exitflag>0,'Velocity alternative optimizer did not converge.');
 Z=table(string(name),base(1,end),a(1,end),10*base(2,end),10*a(2,end),E.estimate.exitflag, ...
  'VariableNames',{'Dataset','BaselineLevel','AlternativeLevel','BaselineRate','AlternativeRate','OptimizerExitFlag'});V=[V;Z]; %#ok<AGROW>
end
writetable(S,fullfile(out,'Sensitivity_Results.csv'));writetable(C,fullfile(out,'ChangepointSensitivity.csv'));
writetable(V,fullfile(out,'reporting','VelocityModel_Comparison.csv'));
end
