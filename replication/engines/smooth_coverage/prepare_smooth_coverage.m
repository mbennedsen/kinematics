function [root,C,S]=prepare_smooth_coverage(mode)
% Read-only verification; safe to call independently from every cluster array.
if nargin<1,mode='production';end
assert(ismember(string(mode),["production","smoke"]));
root=fileparts(mfilename('fullpath'));
snapshot=verify_engine_inputs(root);
L=load(fullfile(fileparts(root),'data','Simulation','Calibration.mat'),'C');C=L.C;
assert(C.schemaVersion==57 && C.sample.T==1748 && C.endTime==2026.625 && C.actualT0==2016.625);
assert(numel(C.selected)==5 && C.selected(3)==4 && strcmp(C.scenarios{C.selected(3)},'NASA_smoothed'));
assert(isequal(C.fitOptions.covariates.maxLags,[12 12 12]) && isequal(C.fitOptions.covariates.almonDegrees,[2 1 1]));
s=C.fitOptions.stateSpace;assert(s.arOrder==1 && s.maOrder==1 && s.useMonthlyVariance && isempty(s.varianceBreakYears) && isempty(s.fixedQ));
o=C.empiricalOptions.ss;
assert(o.bootstrap.ciAlpha==.1 && strcmpi(o.bootstrap.inferenceMode,'bootstrapError') && strcmpi(o.bootstrap.method,'posteriorResidual'));
assert(C.empiricalOptions.kernel.bootstrap.ciAlpha==.1 && C.empiricalOptions.quadratic.inference.ciLevel==.9);
assert(C.empiricalOptions.quadratic.quadratic.windowStartYear==1970 && C.empiricalOptions.changepoint.broken.windowStartYear==1970);
assert(strcmp(which('estimate_monthly_state_space'),fullfile(fileparts(root),'methods','estimate_monthly_state_space.m')),'Estimator shadowed by another folder. Start MATLAB in this package.');
S=struct('mode',char(mode),'coverageR',400,'B',399,'snapshot',snapshot, ...
 'manuscriptDGP',3,'calibrationIndex',C.selected(3),'rebaseData',false);
prefix=root;if strcmp(mode,'smoke'),S.coverageR=1;S.B=19;prefix=fullfile(root,'validation');end
S.coverageDir=fullfile(prefix,sprintf('dgp3_coverage_R%d_B%d',S.coverageR,S.B));
S.summaryDir=fullfile(prefix,'summary');
fprintf('PASS: manuscript DGP 3, %s R=%d B=%d, no individual rebasing.\n',mode,S.coverageR,S.B);
end
