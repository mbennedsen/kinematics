function write_parameter_tables
base=fileparts(mfilename('fullpath'));
root=fileparts(base);
names={'NASA','NOAA','HadCRU','Berk','ERA5'};detail=struct();
for i=1:5
 name=names{i};L=load(fullfile(root,'empirical','results','monthly',name,'ss_B1999','CompletedFit.mat'),'fit');
 E=L.fit.detail.(name);E=rmfield(E,'bootstrap');detail.(name)=E;opts=L.fit.options;
end
opts.datasets=names;folder=fullfile(base,'validated_monthly','parameters');if ~isfolder(folder),mkdir(folder);end
input=fullfile(folder,'Combined_Monthly_Point_Fits.mat');save(input,'detail','opts','-v7.3');
write_ss_parameter_appendix_tables(input,struct('outputDir',folder));
write_diagnostics(input,fullfile(folder,'ResidualDiagnostics.csv'));
fprintf('Updated all baseline SS parameter and residual-diagnostic tables.\n');
end
