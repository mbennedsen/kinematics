function row = run_no_covariates(d,B,nWorkers)
% Table S16 uses the adjusted fits' observation mask and no climate controls.
root=fileparts(mfilename('fullpath'));
snapshot=verify_engine_inputs(root);
names={'NASA','NOAA','HadCRU','Berk','ERA5','Mean'};
name=names{d};
input=fullfile(fileparts(root),'data','Monthly_NoCovariates.csv');
folder=fullfile(root,'results','no_covariates',name,sprintf('ss_B%d',B));
o=struct('datasets',{{name}},'seed',550000+6+d,'rebaseData',true,'outputDir',folder);
o.covariates=struct('names',{{}},'maxLags',[],'almonDegrees',[], ...
    'standardize',false,'includeMonthlyDummies',false);
o.averageAcceleration=struct('enabled',true,'T0',2015);
o.bootstrap=struct('B',B,'ciAlpha',.1,'useParallel',nWorkers>1,'nWorkers',nWorkers, ...
    'pathBands',true,'method','posteriorResidual','inferenceMode','bootstrapError', ...
    'includeParameterUncertainty',true,'startAtOriginalTheta',true,'maxRefitTries',2);
fit=cached_empirical_fit(@estimate_monthly_state_space,input,o,snapshot);
E=fit.detail.(name);
assert(E.pointEstimate.exitflag>0 && isempty(E.pointEstimate.paramsHat.beta));
velocity_window_intervals(E,name,folder);
row=empirical_summary_rows(name,'ss',E,B);
row.Specification=repmat("no_covariates",height(row),1);
writetable(row,fullfile(folder,'Summary.csv'));
end
