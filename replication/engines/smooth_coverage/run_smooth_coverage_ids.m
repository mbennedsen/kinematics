function run_smooth_coverage_ids(ids,nWorkers,mode)
% Explicit-ID recovery; IDs and RNG seeds do not depend on worker count.
if nargin<2,nWorkers=20;end
if nargin<3,mode='production';end
[~,C,S]=prepare_smooth_coverage(mode);
validateattributes(nWorkers,{'numeric'},{'scalar','integer','positive'});
if isempty(ids),fprintf('No requested replications.\n');return;end
validateattributes(ids,{'numeric'},{'vector','integer','>=',1,'<=',S.coverageR});
assert(numel(unique(ids))==numel(ids));
if ~isfolder(S.coverageDir),mkdir(S.coverageDir);end
if nWorkers>1
 pool=gcp('nocreate');if ~isempty(pool)&&pool.NumWorkers~=nWorkers,delete(pool);pool=[];end
 if isempty(pool)
  c=parcluster('local');storage=tempname;mkdir(storage);c.JobStorageLocation=storage;
  parpool(c,nWorkers);
 end
end
for k=1:numel(ids),smooth_coverage_replication(C,S,ids(k),nWorkers);end
end
