function run_smooth_coverage(arrayIndex,nArrays,nWorkers,mode)
if nargin<2,nArrays=12;end
if nargin<3,nWorkers=20;end
if nargin<4,mode='production';end
validateattributes(nArrays,{'numeric'},{'scalar','integer','positive'});
validateattributes(arrayIndex,{'numeric'},{'scalar','integer','>=',1,'<=',nArrays});
[~,~,S]=prepare_smooth_coverage(mode);
run_smooth_coverage_ids(arrayIndex:nArrays:S.coverageR,nWorkers,mode);
end
