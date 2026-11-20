function [root,snapshot,C,S]=prepare_accuracy(mode)
% Read-only verification. Preparation is already frozen; safe in every array.
if nargin<1,mode='production';end
root=fileparts(mfilename('fullpath'));
assert(ismember(string(mode),["production","smoke"]),'Use production or smoke.');
snapshot=verify_engine_inputs(root);
L=load(fullfile(fileparts(root),'data','Simulation','Calibration.mat'),'C');C=L.C;
assert(C.schemaVersion==57 && C.sample.T==1748 && abs(C.endTime-C.actualT0-10)<1e-8);
assert(isequal(C.fitOptions.covariates.maxLags,[12 12 12]) && isequal(C.fitOptions.covariates.almonDegrees,[2 1 1]));
s=C.fitOptions.stateSpace;assert(s.arOrder==1 && s.maOrder==1 && s.useMonthlyVariance && isempty(s.varianceBreakYears) && isempty(s.fixedQ));
S=struct('mode',char(mode),'accuracyR',1000,'snapshot',snapshot);
prefix=root;
if strcmp(mode,'smoke'),S.accuracyR=1;prefix=fullfile(root,'validation');end
S.accuracyDir=fullfile(prefix,sprintf('accuracy_R%d',S.accuracyR));
end
