function run_accuracy(arrayIndex,nArrays,nWorkers,mode)
if nargin<2,nArrays=12;end
if nargin<3,nWorkers=20;end
if nargin<4,mode='production';end
validateattributes(nArrays,{'numeric'},{'scalar','integer','positive'});
validateattributes(arrayIndex,{'numeric'},{'scalar','integer','>=',1,'<=',nArrays});
[~,~,C,S]=prepare_accuracy(mode);
ids=arrayIndex:nArrays:(numel(C.selected)*S.accuracyR);
validateattributes(nWorkers,{'numeric'},{'scalar','integer','positive'});
if ~isfolder(S.accuracyDir),mkdir(S.accuracyDir);end
if nWorkers>1
 pool=gcp('nocreate');
 if ~isempty(pool)&&pool.NumWorkers~=nWorkers,delete(pool);pool=[];end
 if isempty(pool)
  cluster=parcluster('local');storage=tempname;mkdir(storage);
  cluster.JobStorageLocation=storage;parpool(cluster,nWorkers);
 end
 parfor k=1:numel(ids),accuracy_replication(C,S,ids(k));end
else
 for k=1:numel(ids),accuracy_replication(C,S,ids(k));end
end
end
