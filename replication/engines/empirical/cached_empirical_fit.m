function fit=cached_empirical_fit(fun,input,o,snapshot)
if ~isfolder(o.outputDir),mkdir(o.outputDir);end
file=fullfile(o.outputDir,'CompletedFit.mat');
key=struct('inputSHA',file_sha256(input),'options',o, ...
 'sourceSHA',file_sha256(which(func2str(fun))),'packageSHA',snapshot);
if isfile(file)
 L=load(file,'key','fit');[ok,d]=compare_fit_settings(key,L.key);
 assert(ok,'Cached inputs/options differ. Preserve output; use a new package/folder.\n%s',strjoin(d,newline));
 fit=L.fit;fprintf('Reusing completed fit: %s\n',file);return;
end
diary(fullfile(o.outputDir,'Estimation.log'));cleanup=onCleanup(@()diary('off'));started=tic;
fit=fun(input,o);seconds=toc(started);
save(fullfile(o.outputDir,'CompletedFit.tmp.mat'),'key','fit','seconds','-v7.3');
movefile(fullfile(o.outputDir,'CompletedFit.tmp.mat'),file,'f');
end
