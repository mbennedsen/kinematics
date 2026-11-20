function out=collect_smooth_coverage(nArrays,mode)
% Common-DGP summaries, completeness accounting, and saved-interval audit.
if nargin<1,nArrays=12;end
if nargin<2,mode='production';end
[root,C,S]=prepare_smooth_coverage(mode);
validateattributes(nArrays,{'numeric'},{'scalar','integer','positive'});
methods={'ss','kernel','quadratic','changepoint'};objects={'x_T','v_T','a_T','average_a_past_decade'};
[~,fixedTruth,fixedInfo]=simulate_smooth_record(C,1);
Crows=struct([]);missingC=[];retryC=[];Cok=false(S.coverageR,1);Cpresent=Cok;
for r=1:S.coverageR
 f=fullfile(S.coverageDir,sprintf('rep_%04d',r),'Results.mat');
 if ~isfile(f),missingC(end+1)=r;continue;end %#ok<AGROW>
 L=load(f,'R');R=L.R;assert(strcmp(R.snapshot,S.snapshot)&&strcmp(R.mode,S.mode)&&R.rep==r&&R.B==S.B);
 assert(R.manuscriptDGP==3 && R.calibrationIndex==C.selected(3) && ~R.rebaseData);
 assert(max(abs(R.truth-fixedTruth))<1e-10);
 input=fullfile(fileparts(f),'OuterSample.csv');assert(strcmp(R.inputSHA,file_sha256(input)));
 L=load(fullfile(fileparts(f),'OuterTruth.mat'),'R');assert(isequaln(R.truth,L.R.truth));
 assert(L.R.info.rebaseShift==0 && isequaln(L.R.info.alpha,fixedInfo.alpha));
 check_smooth_intervals(R,S);
 Cpresent(r)=true;Cok(r)=all(R.success);if ~Cok(r),retryC(end+1)=r;end %#ok<AGROW>
 for m=1:4
  for k=1:4
   structural=m==4&&k==3;ok=R.success(m)&&~structural;
   cov=NaN;lm=NaN;um=NaN;len=NaN;
   if ok,cov=double(R.low(m,k)<=R.truth(k)&&R.truth(k)<=R.high(m,k));lm=double(R.low(m,k)>R.truth(k));um=double(R.high(m,k)<R.truth(k));len=R.high(m,k)-R.low(m,k);end
   z=struct('Replication',r,'Method',string(methods{m}),'Estimand',string(objects{k}),'Success',R.success(m),'StructuralZero',structural, ...
    'Truth',R.truth(k),'Point',R.point(m,k),'Low',R.low(m,k),'High',R.high(m,k),'Coverage',cov,'LowerMiss',lm,'UpperMiss',um, ...
    'IntervalLength',len,'ValidInnerDraws',R.validDraws(m,k),'ErrorMessage',R.message(m));
   if isempty(Crows),Crows=z;else,Crows(end+1)=z;end %#ok<AGROW>
  end
 end
end
cs=struct([]);if isempty(Crows),T=table();else,T=struct2table(Crows);end
for m=1:4
 for k=1:4
  if isempty(T),z=T;else,z=T(T.Method==methods{m}&T.Estimand==objects{k},:);end
  structural=m==4&&k==3;n=0;np=0;bias=NaN;len=NaN;minB=NaN;maxB=NaN;est=nan(1,3);se=est;wl=est;wh=est;
  if ~isempty(z)
   good=z.Success;np=sum(good);bias=mean(z.Point(good)-z.Truth(good));
   q=z(good&~z.StructuralZero,:);n=height(q);
   if n>0
    val=[q.Coverage,q.LowerMiss,q.UpperMiss];est=mean(val,1);se=sqrt(est.*(1-est)/n);[wl,wh]=wilson_local(sum(val,1),n);
    len=mean(q.IntervalLength);minB=min(q.ValidInnerDraws);maxB=max(q.ValidInnerDraws);
   end
  end
  r=struct('Method',string(methods{m}),'Estimand',string(objects{k}),'Planned',S.coverageR,'Returned',height(z), ...
   'ValidPoints',np,'Failed',height(z)-np,'ValidIntervals',n,'StructuralZero',structural,'Coverage',est(1),'CoverageMCSE',se(1), ...
   'CoverageWilsonLow',wl(1),'CoverageWilsonHigh',wh(1),'LowerMiss',est(2),'LowerMissMCSE',se(2), ...
   'LowerMissWilsonLow',wl(2),'LowerMissWilsonHigh',wh(2),'UpperMiss',est(3),'UpperMissMCSE',se(3), ...
   'UpperMissWilsonLow',wl(3),'UpperMissWilsonHigh',wh(3),'Bias',bias,'MeanLength',len,'B',S.B,'MinValidInner',minB,'MaxValidInner',maxB);
  if isempty(cs),cs=r;else,cs(end+1)=r;end %#ok<AGROW>
 end
end
summary=S.summaryDir;if ~isfolder(summary),mkdir(summary);end
writetable(struct2table(cs),fullfile(summary,'CoverageSummary.csv'));
if ~isempty(Crows),writetable(T,fullfile(summary,'CoverageReplications.csv'));end
out=struct('missingCoverage',missingC,'retryCoverage',retryC, ...
 'coverageComplete',all(Cok),'snapshot',S.snapshot,'mode',mode,'manuscriptDGP',3, ...
 'planned',S.coverageR,'returned',sum(Cpresent),'complete',sum(Cok),'B',S.B);
for label={'missingCoverage','retryCoverage'}
 ids=out.(label{1});writetable(table(ids(:),'VariableNames',{'Replication'}),fullfile(summary,[label{1} '.csv']));
end
P=table((1:nArrays)','VariableNames',{'ArrayIndex'});
for name={'Assigned','Missing','Retry'},P.(name{1})=zeros(nArrays,1);end
for a=1:nArrays
 ic=a:nArrays:S.coverageR;P.Assigned(a)=numel(ic);P.Missing(a)=sum(~Cpresent(ic));P.Retry(a)=sum(Cpresent(ic)&~Cok(ic));
end
writetable(P,fullfile(summary,'ArrayProgress.csv'));save(fullfile(summary,'Collection.mat'),'out','S');
fprintf('DGP3 COVERAGE: %d/%d complete; %d missing; %d need retry; B=%d.\n',sum(Cok),S.coverageR,numel(missingC),numel(retryC),S.B);
if ~all(Cok),fprintf('INCOMPLETE: do not treat partial summaries as final.\n');else,fprintf('COMPLETE: all four methods and interval checks passed.\n');end
disp(P);
end
function [lo,hi]=wilson_local(k,n)
z=1.95996398454005;p=k/n;den=1+z^2/n;c=(p+z^2/(2*n))/den;h=z*sqrt(p.*(1-p)/n+z^2/(4*n^2))/den;lo=c-h;hi=c+h;
end
