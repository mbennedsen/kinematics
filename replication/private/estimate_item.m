function estimate_item(root,run,item,opt)
% Stage the exact inputs and engines once per RunID.
work=fullfile(run,'work');dest=fullfile(run,'results');
source=fullfile(root,'replication','engines');
if ~isfolder(run),mkdir(run);end
if ~isfolder(work)
    pending=tempname(run);mkdir(pending);
    clean=onCleanup(@() remove_staging(pending));
    copyfile(fullfile(source,'*'),pending);
    stage_data(root,pending);
    movefile(pending,work);
end
manifest='SHA256SUMS.tsv';
assert(isfile(fullfile(work,manifest)) && ...
    strcmp(fileread(fullfile(source,manifest)),fileread(fullfile(work,manifest))), ...
    'The RunID belongs to a different release. Choose a new RunID.');
if ~isfolder(dest),mkdir(dest);end
emp=fullfile(work,'empirical');rob=fullfile(work,'sensitivity');sim=fullfile(work,'accuracy');
addpath(fullfile(work,'methods'),emp,rob,sim,fullfile(work,'smooth_coverage'), ...
    fullfile(work,'reference_coverage'),'-begin');
names={'NASA','NOAA','HadCRU','Berk','ERA5','Mean'};methods={'ss','kernel','quadratic','changepoint'};
monthly=[];annual=false;robust=false;windows=false;params=false;noCov=false;accuracy=false;coverage=false;reference=false;
if ismember(item,{'figure1','figure3','FigureS9'}),monthly=6;
elseif strcmp(item,'figure2'),monthly=1:6;
elseif startsWith(item,'FigureS')
 n=str2double(extractAfter(item,'FigureS'));if n>=3&&n<=8,order=[6 1 2 3 4 5];monthly=order(n-2);end
else
 n=str2double(extractAfter(item,'TableS'));
 switch n
  case 1,monthly=1:6;
  case {2,3},monthly=6;
  case {5,6},accuracy=true;
  case 9,annual=true;
  case {10,11,12,13,14},monthly=1:5;methods={'ss'};robust=true;
  case 15,monthly=1:5;windows=true;
  case 16,noCov=true;
  case 7,coverage=true;
  case 8,reference=true;
  case {17,18,19,20,21,22},monthly=1:5;methods={'ss'};params=true;
  case 23,monthly=1:6;methods={'changepoint'};annual=true;
 end
end
% Data-only illustrations and design truth have no fitted estimator to rerun.
if strcmp(item,'FigureS10')||strcmp(item,'TableS4')
 make_design(dest);
end
for d=monthly
 run_empirical(d,1999,opt.Workers,methods);
 collect_empirical(emp,dest,'monthly',names{d},methods);
end
if annual
 am=methods;if strcmp(item,'TableS9'),am={'ss','kernel','quadratic','changepoint'};end
 for d=1:5
  run_empirical(d+6,1999,opt.Workers,am);
  collect_empirical(emp,dest,'annual',names{d},am);
 end
end
if robust
 for d=1:5,run_sensitivity(d,rob);end
 collect_robustness(rob,emp,dest);
end
if windows
 for d=1:5,run_kernel_windows(d,rob);end
 collect_windows(rob,emp,dest);
end
if params
 write_parameter_tables;
 target=fullfile(dest,'empirical','validated_monthly','parameters');if ~isfolder(target),mkdir(target);end
 copyfile(fullfile(rob,'validated_monthly','parameters','*'),target);
end
if noCov
 T=table();
 for d=1:6
  row=run_no_covariates(d,1999,opt.Workers);T=[T;row]; %#ok<AGROW>
 end
 target=fullfile(dest,'no_covariates');if ~isfolder(target),mkdir(target);end
 writetable(T,fullfile(target,'Endpoint_Summary.csv'));
end
if accuracy
 run_accuracy(1,1,opt.Workers);
 status=collect_accuracy();
 assert(status.accuracyComplete,'Incomplete accuracy experiment.');
 target=fullfile(dest,'simulations','accuracy'); if ~isfolder(target),mkdir(target);end
 copyfile(fullfile(sim,'summary','AccuracySummary.csv'),target);
end
if coverage
 study = fullfile(work,'smooth_coverage');
 run_smooth_coverage(1,1,opt.Workers);
 status = collect_smooth_coverage(1);
 assert(status.coverageComplete,'Incomplete smooth-path experiment.');
 target = fullfile(dest,'simulations','smooth_coverage');
 if ~isfolder(target), mkdir(target); end
 copyfile(fullfile(study,'summary','CoverageSummary.csv'),target);
end
if reference
 study = fullfile(work,'reference_coverage');
 for design = ["quadratic_truth","changepoint_truth"]
  run_parametric_reference_study(design,1:400,opt.Workers);
 end
 target = fullfile(dest,'simulations','reference_coverage');
 if ~isfolder(target), mkdir(target); end
 results = collect_reference_coverage(study,1:400);
 writetable(results,fullfile(target,'Summary.csv'));
end
end

function remove_staging(folder)
if isfolder(folder),rmdir(folder,'s');end
end
