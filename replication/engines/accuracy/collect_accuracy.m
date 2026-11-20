function out=collect_accuracy(mode)
% Partial results are reported explicitly; missing/failed IDs can be rerun.
if nargin<1,mode='production';end
[root,~,C,S]=prepare_accuracy(mode);
methods={'ss','kernel','quadratic','changepoint'};objects={'x_T','v_T','a_T','average_a_past_decade'};
Arows=struct([]);Asummary=struct([]);missingA=[];retryA=[];Aok=false(5*S.accuracyR,1);
for di=1:5
 errors=nan(S.accuracyR,4,4);attempted=0;
 for r=1:S.accuracyR
  job=(di-1)*S.accuracyR+r;f=fullfile(S.accuracyDir,sprintf('job_%05d.mat',job));
  if ~isfile(f),missingA(end+1)=job;continue;end %#ok<AGROW>
  L=load(f,'R');R=L.R;assert(strcmp(R.snapshot,S.snapshot)&&strcmp(R.mode,S.mode)&&R.job==job&&R.rep==r&&R.dgp==C.selected(di));
  attempted=attempted+1;Aok(job)=all(R.success);if ~Aok(job),retryA(end+1)=job;end %#ok<AGROW>
  for m=1:4
   if R.success(m),errors(r,m,:)=R.estimates(m,:)-R.truth;end
   for k=1:4
    z=struct('JobID',job,'DGP',string(C.scenarios{C.selected(di)}),'Replication',r,'Method',string(methods{m}), ...
     'Estimand',string(objects{k}),'Success',R.success(m),'Truth',R.truth(k),'Estimate',R.estimates(m,k),'ErrorMessage',R.message(m));
    if isempty(Arows),Arows=z;else,Arows(end+1)=z;end %#ok<AGROW>
   end
  end
 end
 for m=1:4
  for k=1:4
   e=errors(:,m,k);e=e(isfinite(e));n=numel(e);rmse=sqrt(mean(e.^2));
   z=struct('DGP',string(C.scenarios{C.selected(di)}),'Method',string(methods{m}),'Estimand',string(objects{k}), ...
    'Planned',S.accuracyR,'Returned',attempted,'Valid',n,'Failed',attempted-n,'Bias',mean(e), ...
    'BiasMCSE',std(e)/sqrt(n),'RMSE',rmse,'RMSE_MCSE',std(e.^2)/sqrt(n)/(2*max(rmse,eps)));
   if isempty(Asummary),Asummary=z;else,Asummary(end+1)=z;end %#ok<AGROW>
  end
 end
end
summary=fullfile(root,'summary');if strcmp(mode,'smoke'),summary=fullfile(root,'validation','summary');end
if ~isfolder(summary),mkdir(summary);end
writetable(struct2table(Asummary),fullfile(summary,'AccuracySummary.csv'));
if ~isempty(Arows),writetable(struct2table(Arows),fullfile(summary,'AccuracyReplications.csv'));end
out=struct('missingAccuracy',missingA,'retryAccuracy',retryA, ...
 'accuracyComplete',all(Aok),'snapshot',S.snapshot,'mode',mode);
for label={'missingAccuracy','retryAccuracy'}
 ids=out.(label{1});writetable(table(ids(:),'VariableNames',{'JobID'}),fullfile(summary,[label{1} '.csv']));
end
fprintf('ACCURACY: %d/%d complete; %d missing; %d need retry.\n',sum(Aok),numel(Aok),numel(missingA),numel(retryA));
end
