function stage_data(root,work)
% Each estimation workspace holds one complete copy of the frozen inputs.
source=fullfile(root,'data');
folder=fullfile(work,'data');mkdir(folder);
files={'Monthly_HadISST.csv','Annual_HadISST.csv','Monthly_NoCovariates.csv', ...
    'Annual_Baseline_Offsets.csv','Monthly_Smoothing_References.csv','Configuration.mat'};
for j=1:numel(files)
    copyfile(fullfile(source,files{j}),fullfile(folder,files{j}));
end
copyfile(fullfile(source,'Simulation'),fullfile(folder,'Simulation'));
end
