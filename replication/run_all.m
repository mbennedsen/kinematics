function run_all(mode,options)
% RUN_ALL Reproduce all 13 figures and 23 supplementary tables.
% run_all('saved') is fast. run_all('reestimate') is a large research workload.
% A single fresh run is shared across items to avoid repeating fitted models.
if nargin<1,mode='saved';end
if nargin<2,options=struct();end
if strcmpi(mode,'reestimate') && ~isfield(options,'RunID')
 options.RunID=['all_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))];
end
[items,~]=item_catalog;
for k=1:numel(items),feval(items{k},mode,options);end
end
