function result=replicate_item(item,mode,options)
% One interface for both sources; neither plotting nor table code fits models.
[items,titles]=item_catalog;ix=find(strcmp(items,item));assert(isscalar(ix));
mode=validatestring(mode,{'saved','reestimate'});
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~isfield(options,'Workers'),options.Workers=1;end
validateattributes(options.Workers,{'numeric'},{'scalar','integer','positive'});
if ~isfield(options,'Visible'),options.Visible='off';end
oldpath=path;restore=onCleanup(@()path(oldpath));
if strcmp(mode,'saved')
 source=fullfile(root,'outputs','saved');out=fullfile(root,'outputs','reproduced','saved');
else
 if ~isfield(options,'RunID'),options.RunID=['run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))];end
 assert(~isempty(regexp(options.RunID,'^[A-Za-z0-9_-]+$','once')),'RunID may contain letters, numbers, underscores and hyphens only.');
 run=fullfile(root,'outputs','reestimated',options.RunID);source=fullfile(run,'results');out=fullfile(run,'reproduced');
 fprintf('Re-estimation run: %s\n',run);
 estimate_item(root,run,item,options);
end
if ~isfolder(out),mkdir(out);end
if startsWith(item,'Table')
 [values,panels]=paper_table(item,source);
 emit_table(item,titles{ix},values,panels,out);
 result=values;
else
 result=paper_figure(item,source,fullfile(root,'data'),options.Visible);
 exportgraphics(result,fullfile(out,[item '.pdf']),'ContentType','vector');
 exportgraphics(result,fullfile(out,[item '.png']),'Resolution',200);
 savefig(result,fullfile(out,[item '.fig']));
 if strcmp(options.Visible,'off'),close(result);result=fullfile(out,[item '.pdf']);end
end
provenance=struct('Item',item,'Title',titles{ix},'Mode',mode,'Source',source, ...
 'MATLAB',version,'Created',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss z')));
save(fullfile(out,[item '_provenance.mat']),'provenance');
fprintf('%s reproduced: %s\n',item,out);
end
