function summary=run_empirical(jobIndex,B,nWorkers,methods)
% Jobs 1:6 monthly NASA/NOAA/HadCRU/Berk/ERA5/Mean; 7:11 annual first five.
% Whole-record jobs, NOT chunks. All four methods by default.
if nargin<2,B=1999;end
if nargin<3,nWorkers=20;end
if nargin<4,methods={'quadratic','changepoint','kernel','ss'};end
validateattributes(jobIndex,{'numeric'},{'scalar','integer','>=',1,'<=',11});
validateattributes(B,{'numeric'},{'scalar','integer','>=',1});
validateattributes(nWorkers,{'numeric'},{'scalar','integer','>=',1});
assert(iscell(methods)&&all(ismember(methods,{'ss','kernel','quadratic','changepoint'})));
assert(~ismember('changepoint',methods)||B>=19,'Changepoint smoke tests require B >= 19.');
[root,snapshot,C]=prepare_empirical;old=pwd;cleanup=onCleanup(@()cd(old));cd(root);
annual=jobIndex>6;d=jobIndex-6*annual;name=C.names{d};kind='monthly';dataName='Monthly_HadISST.csv';lag=120;
if annual,kind='annual';dataName='Annual_HadISST.csv';lag=10;end
input=fullfile(fileparts(root),'data',dataName);inputSHA=file_sha256(input);
raw=readtable(input);last=find(isfinite(raw.(name)),1,'last');T0=raw.t(last-lag);
assert(abs(raw.t(last)-T0-10)<1e-8);
jobFolder=fullfile(root,'results',kind,name);if ~isfolder(jobFolder),mkdir(jobFolder);end
summary=table();
for mm=1:numel(methods)
 method=methods{mm};folder=fullfile(jobFolder,sprintf('%s_B%d',method,B));
 if ~isfolder(folder),mkdir(folder);end
 o=C.(kind).(name).(method);o.datasets={name};o.outputDir=folder;
 fprintf('\nJOB %d: %s / %s / %s / B=%d / workers=%d\n',jobIndex,kind,name,method,B,nWorkers);
 if strcmp(method,'ss') && annual
  o.outputDir=fullfile(folder,'point');o.bootstrap.B=0;o.bootstrap.useParallel=false;
  fit=cached_empirical_fit(@estimate_annual_state_space,input,o,snapshot);
  E=fit.detail.(name);assert(E.pointEstimate.exitflag>0,'Annual point optimization failed.');
  record=bootstrap_annual_state_space('record',input,o,E.pointEstimate,name,20261008+100000*d,inputSHA);
  if nWorkers>1,ensure_pool(nWorkers);end
  nw=nWorkers;if nw==1,nw=0;end
  z=bootstrap_annual_state_space('run',record,B,nw,folder,100);
  idx=[1 2 3 7];p=z.summary.Estimate(idx)';lo=z.summary.CI_Low(idx)';hi=z.summary.CI_High(idx)';
  valid=z.summary.SuccessfulDraws(1);t=record.sample.t;flag=E.pointEstimate.exitflag;br=NaN;
  inference="Bootstrap error, conditional on corrected monthly sigma";
 elseif strcmp(method,'ss')
  o.bootstrap.B=B;o.bootstrap.useParallel=nWorkers>1;o.bootstrap.nWorkers=nWorkers;
  o.averageAcceleration.T0=T0;
  fit=cached_empirical_fit(@estimate_monthly_state_space,input,o,snapshot);E=fit.detail.(name);
  assert(E.pointEstimate.exitflag>0,'Monthly point optimization failed.');
  t=E.sampleData.t;P=kinematic_path_values(E.pointEstimate.smoother.alphaSmooth);p=P(end,:);
  bands=E.bootstrap.pathBands;lo=bands.low(end,:);hi=bands.high(end,:);valid=bands.validDraws;
  fields={'position','velocity','acceleration'};
  for kk=1:3
   lo(kk)=E.bootstrap.(fields{kk}).endpoint_ci_low;
   hi(kk)=E.bootstrap.(fields{kk}).endpoint_ci_high;
  end
  assert(valid>=ceil(.95*B),'Fewer than 95%% valid paired bootstrap paths. Inspect saved fit.');
  flag=E.pointEstimate.exitflag;br=NaN;inference="Bootstrap error";
  velocity_window_intervals(E,name,folder);
  assert(abs(p(4)-E.pointAverageAcceleration.estimate)<1e-8);
 else
  switch method
   case 'kernel'
    o.bootstrap.B=B;o.bootstrap.useParallel=false;o.kernel.averageStartYear=T0;
    fun=@estimate_local_quadratic;
   case 'quadratic',fun=@estimate_quadratic_benchmark;
   case 'changepoint',fun=@estimate_broken_trend_benchmark;
  end
  if strcmp(method,'changepoint'),o.outputDir=fullfile(folder,'point');end
  fit=cached_empirical_fit(fun,input,o,snapshot);E=fit.detail.(name);t=E.sample.t;flag=NaN;br=NaN;
  if strcmp(method,'changepoint')
   noise='ARMA11';if annual,noise='AR1';end
   starts=T0;labels={'past_decade'};
   for yr=[2010 2012 2014 2015 2016 2018 2020]
    starts(end+1)=t(find(t>=yr,1));labels{end+1}=sprintf('since_%d',yr); %#ok<AGROW>
   end
   opt=struct('noiseModel',noise,'B',B,'alpha',.1,'seed',558000+10000*(1+2*annual)+3000*d, ...
    'windowT0',starts,'windowNames',{labels});
   key=struct('inputSHA',inputSHA,'options',opt,'packageSHA',snapshot);file=fullfile(folder,'Bootstrap.mat');
   if isfile(file)
    L=load(file,'key','result');assert(isequaln(key,L.key),'Changepoint inputs/options changed.');z=L.result;
   else
    result=bootstrap_broken_trend_arma11(E,opt);z=result;save(file,'key','result','-v7.3');
   end
   p=z.estimates;lo=z.intervals(:,1)';hi=z.intervals(:,2)';valid=z.validDraws;br=E.breakDate;
   inference=string(noise)+" full-search bootstrap percentile";
   writetable(z.windows,fullfile(folder,'AverageAcceleration_Windows.csv'));
   draws=table((1:B)',z.seeds',z.breakDates',z.draws(:,1),z.draws(:,2),z.draws(:,3),z.draws(:,4),z.slopeChangeDraws, ...
    'VariableNames',{'Replication','Seed','BreakDate','x','v','a','average_a','SlopeChange'});
   writetable(draws,fullfile(folder,'Draws.csv'));
  else
   p=[E.position.endpoint,E.velocity.endpoint,E.acceleration.endpoint,NaN];
   lo=[E.position.endpoint_ci_low,E.velocity.endpoint_ci_low,E.acceleration.endpoint_ci_low,NaN];
   hi=[E.position.endpoint_ci_high,E.velocity.endpoint_ci_high,E.acceleration.endpoint_ci_high,NaN];
   if strcmp(method,'quadratic')
    p(4)=p(3);lo(4)=lo(3);hi(4)=hi(3);valid=NaN;inference="Newey-West HAC";
   else
    p(4)=E.averageAcceleration.estimate;lo(4)=E.averageAcceleration.ci_low;hi(4)=E.averageAcceleration.ci_high;
    valid=sum(isfinite(E.averageAcceleration.bootstrap));inference="Moving-block bootstrap percentile";
    assert(valid>=ceil(.95*B),'Fewer than 95%% valid kernel draws.');
   end
  end
 end
 assert(all(isfinite(p)) && abs(t(end)-T0-10)<1e-8);
 keep=1:4;if strcmp(method,'changepoint'),keep=[1 2 4];assert(p(3)==0 && isnan(lo(3)) && isnan(hi(3)));end
 assert(all(isfinite([lo(keep),hi(keep)])) && all(lo(keep)<=hi(keep)));
 draws=B;if strcmp(method,'quadratic'),draws=0;end
 R=table(repmat(string(name),4,1),repmat(string(method),4,1),["x";"v";"a";"average_a"],p',lo',hi', ...
  repmat(t(1),4,1),repmat(t(end),4,1),repmat(T0,4,1),repmat(draws,4,1),repmat(valid,4,1), ...
  repmat(br,4,1),repmat(flag,4,1),repmat(inference,4,1),repmat(string(kind),4,1), ...
  'VariableNames',{'Dataset','Method','Object','Estimate','CILow','CIHigh','StartTime','EndTime', ...
  'ActualT0','BootstrapB','ValidDraws','BreakDate','OptimizerExitFlag','Inference','Specification'});
 writetable(R,fullfile(folder,'Summary.csv'));summary=[summary;R]; %#ok<AGROW>
 fprintf('PASS %s/%s/%s: updated intervals saved.\n',kind,name,method);
end
allMethods={'ss','kernel','quadratic','changepoint'};hashes=cell(4,1);complete=true;
for j=1:4
 file=fullfile(jobFolder,sprintf('%s_B%d',allMethods{j},B),'Summary.csv');
 if ~isfile(file),complete=false;break;end
 hashes{j}=file_sha256(file);
end
if complete
 save(fullfile(jobFolder,sprintf('JobComplete_B%d.mat',B)),'jobIndex','B','snapshot','allMethods','hashes');
 fprintf('COMPLETE: job %d, all four methods, B=%d.\n',jobIndex,B);
end
end

function ensure_pool(n)
p=gcp('nocreate');if ~isempty(p)&&p.NumWorkers~=n,delete(p);p=[];end
if isempty(p),parpool('local',n);end
end
