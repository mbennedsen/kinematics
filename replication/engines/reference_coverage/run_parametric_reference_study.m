function run_parametric_reference_study(design,ids,nWorkers)
% Correct trend specification: matching benchmark, NASA irregular calibration.
if nargin<2,ids=1:400;end
if nargin<3,nWorkers=1;end
assert(ismember(string(design),["quadratic_truth","changepoint_truth"]));
base=fileparts(mfilename('fullpath'));
root=fullfile(base,char(design));if ~isfolder(root),mkdir(root);end
validateattributes(ids,{'numeric'},{'vector','integer','>=',1,'<=',400});
assert(numel(unique(ids))==numel(ids));
snapshot=verify_engine_inputs(base);
L=load(fullfile(fileparts(base),'data','Simulation','Calibration.mat'),'C');C=L.C;
C.referenceMethod=1+strcmp(design,'changepoint_truth');
assert(C.schemaVersion==57 && C.sample.T==1748 && abs(C.endTime-C.actualT0-10)<1e-8);
assert(C.empiricalOptions.quadratic.quadratic.windowStartYear==1970);
assert(C.empiricalOptions.changepoint.broken.windowStartYear==1970);
alpha=reference_path(C);
truth=[alpha(1,end),10*alpha(2,end),100*alpha(3,end), ...
 100*(alpha(2,end)-alpha(2,C.index0))/(C.endTime-C.actualT0)];
if C.referenceMethod==1,assert(truth(3)>0 && abs(truth(3)-truth(4))<1e-10);
else,assert(all(truth(3:4)==0));end
P=table(C.sample.t,alpha(1,:)',10*alpha(2,:)',100*alpha(3,:)', ...
 'VariableNames',{'t','x','v','a'});writetable(P,fullfile(root,'FixedPath.csv'));
D=struct('snapshot',snapshot,'truth',truth,'endTime',C.endTime,'averageStart',C.actualT0, ...
 'outerR',400,'innerB',399,'definition','Fitted parametric trend and covariate component held fixed; matching benchmark fitted','design',char(design), ...
 'quadraticOptions',C.empiricalOptions.quadratic,'changepointOptions',C.empiricalOptions.changepoint);
save_checkpoint(fullfile(root,'Design.mat'),D);
fprintf('Fixed reference path: x=%.9f, v=%.9f, a=%.9f, average a=%.9f; R=400, B=399.\n',truth);
if nWorkers>1
 pool=gcp('nocreate');if isempty(pool),parpool('Processes',nWorkers);else,assert(pool.NumWorkers==nWorkers);end
 parfor j=1:numel(ids),one_rep(C,root,snapshot,ids(j));end
else
 for j=1:numel(ids),one_rep(C,root,snapshot,ids(j));end
end
fprintf('Requested parametric reference replications finished.\n');
end

function one_rep(C,root,snapshot,rep)
maxNumCompThreads(1);
folder=fullfile(root,'results',sprintf('rep_%04d',rep));if ~isfolder(folder),mkdir(folder);end
file=fullfile(folder,'Results.mat');input=fullfile(folder,'OuterSample.csv');
if isfile(file)
 L=load(file,'R');R=L.R;assert(strcmp(R.snapshot,snapshot)&&R.rep==rep&&R.B==399);
 assert(strcmp(R.inputSHA,file_sha256(input)));
 if R.success(C.referenceMethod),return;end
else
 [raw,truth,info]=outer(C,rep);writetable(raw,input);
 R=struct('snapshot',snapshot,'rep',rep,'B',399,'inputSHA',file_sha256(input), ...
  'truth',truth,'outerInfo',info,'testedMethod',C.referenceMethod,'point',nan(2,4),'low',nan(2,4),'high',nan(2,4), ...
  'success',false(2,1),'validDraws',nan(2,4),'seconds',nan(2,1),'message',strings(2,1));
end
methods={'quadratic','changepoint'};
for m=C.referenceMethod
 if R.success(m),continue;end
 work=tempname(folder);mkdir(work);clean=onCleanup(@()rmdir(work,'s'));tick=tic;
 try
  o=C.empiricalOptions.(methods{m});o.datasets={'NASA'};o.rebaseData=true;o.outputDir=work;
  o.seed=591000000+10000*rep+1000*m;
  if m==1
   evalc('fit=estimate_quadratic_benchmark(input,o);');E=fit.detail.NASA;
   fields={'position','velocity','acceleration'};point=zeros(1,4);lo=point;hi=point;
   for k=1:3,z=E.(fields{k});point(k)=z.endpoint;lo(k)=z.endpoint_ci_low;hi(k)=z.endpoint_ci_high;end
   point(4)=point(3);lo(4)=lo(3);hi(4)=hi(3);
   R.quadraticCoefficients=E.coeff;R.quadraticOptions=o;
  else
   evalc('fit=estimate_broken_trend_benchmark(input,o);');E=fit.detail.NASA;
   opt=struct('noiseModel','ARMA11','B',399,'alpha',.1,'seed',o.seed, ...
    'windowT0',C.actualT0,'windowNames',{{'past_decade'}});
   b=bootstrap_broken_trend_arma11(E,opt);
   point=b.estimates;lo=b.intervals(:,1)';hi=b.intervals(:,2)';
   R.draws=b.draws;R.validDraws(m,:)=sum(isfinite(b.draws),1);R.validDraws(m,3)=NaN;
   R.breakDate=E.breakDate;R.breakDates=b.breakDates;R.cpNoise=b.noise;R.cpCoefficients=b.coefficients;
   R.candidateBreakDates=E.candidateBreakDates;R.referenceSearchChecks=b.referenceSearchChecks;
   R.changepointOptions=o;R.bootstrapOptions=opt;
   assert(point(3)==0 && isnan(lo(3)) && isnan(hi(3)));
   assert(all(R.validDraws(m,[1 2 4])==399));
  end
  assert(abs(E.sample.t(end)-C.endTime)<1e-8);
  assert(E.sample.t(1)>=1970 && E.sample.t(1)<1970.1);
  keep=1:4;if m==2,keep=[1 2 4];end
  assert(all(isfinite([point lo(keep) hi(keep)]))&&all(lo(keep)<=hi(keep)));
  R.point(m,:)=point;R.low(m,:)=lo;R.high(m,:)=hi;R.success(m)=true;R.message(m)="";
 catch ME
  R.message(m)=string(getReport(ME,'extended','hyperlinks','off'));
 end
 R.seconds(m)=toc(tick);R.matlabVersion=version;save_checkpoint(file,R);clear clean
 fprintf('rep %d/400 %s success=%d %.2fs\n',rep,methods{m},R.success(m),R.seconds(m));
 if ~R.success(m),fprintf('%s\n',R.message(m));end
end
end

function [raw,truth,info]=outer(C,rep)
s=C.sample;p=C.est.paramsHat;n=s.T;alpha=reference_path(C);
% Same independent outer noise streams as the two completed fixed-null designs.
seed=570110300+rep;rng(seed,'twister');noise=zeros(n,1);e=0;lastU=0;
for k=(1-C.cfg.burninMonths):n
 if k<1,month=mod(s.month(1)+k-2,12)+1;era=s.era(1);else,month=s.month(k);era=s.era(k);end
 scale=1;
 if C.est.setup.useMonthlyVariance,scale=scale*exp(p.gamma(month));end
 if numel(p.delta)>=era,scale=scale*exp(p.delta(era));end
 u=sqrt(max(p.sig2u*scale,C.est.setup.minVariance))*randn;
 e=p.ar(1)*e+u+p.maCoeffs(1)*lastU;lastU=u;if k>=1,noise(k)=e;end
end
methods={'quadratic','changepoint'};P=C.parametric.(methods{C.referenceMethod});
y=alpha(1,:)'+P.covariateContribution+noise;
raw=C.raw;raw.NASA(:)=NaN;raw.NASA(C.sampleIndex)=y;
ref=s.t>=1991&s.t<2021;shift=.88-mean(y(ref));
truth=[alpha(1,end)+shift,10*alpha(2,end),100*alpha(3,end), ...
 100*(alpha(2,end)-alpha(2,C.index0))/(C.endTime-C.actualT0)];
info=struct('seed',seed,'rebaseShift',shift,'noise',noise,'averageStart',C.actualT0,'endTime',C.endTime);
end

function alpha=reference_path(C)
t=C.sample.t(:);alpha=zeros(3,numel(t));methods={'quadratic','changepoint'};
P=C.parametric.(methods{C.referenceMethod});b=P.coeff;u=t-P.tRef;
if C.referenceMethod==1
 alpha(1,:)=(b(1)+b(2)*u+b(3)*u.^2)';
 alpha(2,:)=(b(2)+2*b(3)*u)';alpha(3,:)=2*b(3);
else
 alpha(1,:)=(b(1)+b(2)*u+b(3)*max(0,t-P.breakDate))';
 alpha(2,:)=(b(2)+b(3)*double(t>=P.breakDate))';
 alpha(3,abs(t-P.breakDate)<1e-9)=NaN;
 assert(P.breakDate<C.actualT0);
end
end
