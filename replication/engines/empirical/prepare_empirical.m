function [root,snapshot,C]=prepare_empirical
% Read-only preflight. Never regenerates inputs or hashes user launch scripts.
root=fileparts(mfilename('fullpath'));
snapshot=verify_engine_inputs(root);
data=fullfile(fileparts(root),'data');
L=load(fullfile(data,'Configuration.mat'),'C');C=L.C;
M=readtable(fullfile(data,'Monthly_HadISST.csv'));A=readtable(fullfile(data,'Annual_HadISST.csv'));
assert(max(abs(diff(M.t)-1/12))<1e-8 && all(diff(A.year)==1) && A.year(end)==2025);
for j=1:5
 name=C.names{j};assert(abs(mean(A.(name)(A.year>=1991 & A.year<=2020))-.88)<1e-12);
 assert(~C.annual.(name).ss.rebaseData);
 assert(C.annual.(name).ss.stateSpace.fixedQ==C.monthly.(name).fixedQReference);
end
assert(contains(C.indexConvention,'Ordinary HadISST'));
fprintf('PASS: verified empirical inputs; 6 monthly and 5 annual records.\n');
end
