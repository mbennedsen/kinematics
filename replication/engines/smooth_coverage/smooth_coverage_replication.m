function smooth_coverage_replication(C,S,rep,nWorkers)
folder=fullfile(S.coverageDir,sprintf('rep_%04d',rep));
if ~isfolder(folder),mkdir(folder);end
file=fullfile(folder,'Results.mat');sampleFile=fullfile(folder,'OuterSample.csv');outerFile=fullfile(folder,'OuterTruth.mat');
key=[S.snapshot ':' S.mode ':' num2str(rep)];
if isfile(outerFile)
 L=load(outerFile,'R');O=L.R;assert(strcmp(O.key,key));
 assert(strcmp(O.inputSHA,file_sha256(sampleFile)),'Outer sample changed.');
else
 [raw,truth,info]=simulate_smooth_record(C,rep);writetable(raw,sampleFile);
 O=struct('key',key,'truth',truth,'info',info,'seed',info.seed,'inputSHA',file_sha256(sampleFile));
 save_checkpoint(outerFile,O);
end
R=struct('snapshot',S.snapshot,'mode',S.mode,'rep',rep,'B',S.B,'inputSHA',O.inputSHA, ...
 'truth',O.truth,'manuscriptDGP',3,'calibrationIndex',C.selected(3),'rebaseData',false,'point',nan(4,4),'low',nan(4,4),'high',nan(4,4), ...
 'success',false(4,1),'validDraws',nan(4,4),'draws',{cell(4,1)}, ...
 'seconds',nan(4,1),'message',strings(4,1));
if isfile(file)
 L=load(file,'R');R=L.R;
 assert(strcmp(R.snapshot,S.snapshot)&&strcmp(R.mode,S.mode)&&R.B==S.B&&strcmp(R.inputSHA,O.inputSHA));
 if all(R.success),return;end
end
methods={'ss','kernel','quadratic','changepoint'};
for m=1:4
 if R.success(m),continue;end
 tick=tic;work=tempname(folder);mkdir(work);clean=onCleanup(@()rmdir(work,'s'));
 try
  o=C.empiricalOptions.(methods{m});o.datasets={'NASA'};o.rebaseData=false;o.outputDir=work;
  o.seed=59260820+100000*rep+1000*m;
  switch m
   case 1
    inner=fullfile(folder,'ss_inner');if ~isfolder(inner),mkdir(inner);end
    o.averageAcceleration=struct('enabled',true,'T0',C.actualT0);
    o.bootstrap.B=S.B;o.bootstrap.useParallel=nWorkers>1;o.bootstrap.nWorkers=nWorkers;o.bootstrap.pathBands=false;
    o.bootstrap.checkpointDir=inner;o.bootstrap.checkpointKey=[key ':' O.inputSHA];
    retry_failed_draws(inner,o.bootstrap.checkpointKey,S.B);
    evalc('fit=estimate_monthly_state_space(sampleFile,o);');E=fit.detail.NASA;
    assert(E.pointEstimate.exitflag>0);assert(max(abs(E.sampleData.t-C.sample.t))<1e-8);
    z=E.pointEstimate.smoother.alphaSmooth;
    point=[z(1,end),10*z(2,end),100*z(3,end),E.pointAverageAcceleration.estimate];
    lo=zeros(1,4);hi=lo;valid=lo;fields={'position','velocity','acceleration'};
    for k=1:3
     b=E.bootstrap.(fields{k});lo(k)=b.endpoint_ci_low;hi(k)=b.endpoint_ci_high;
    end
    b=E.bootstrap.averageAcceleration;lo(4)=b.endpoint_ci_low;hi(4)=b.endpoint_ci_high;
    % Read paired, unfiltered errors from restart files, preserving failed draws.
    draws=nan(S.B,4);
    for b=1:S.B
     d=load(fullfile(inner,sprintf('draw_%05d.mat',b)));
     assert(strcmp(d.key,o.bootstrap.checkpointKey));draws(b,:)=[d.xRow(1),d.vRow(1),d.aRow(1),d.averageAccelerationRow];
    end
    valid=sum(isfinite(draws),1);
    R.ssTheta=E.pointEstimate.thetaHat;R.ssLogLik=E.pointEstimate.logLik;R.ssExitflag=E.pointEstimate.exitflag;
    R.ssPath=z(1:3,:);R.pathTimes=E.sampleData.t;
   case 2
    o.kernel.averageStartYear=C.actualT0;o.bootstrap.B=S.B;o.bootstrap.useParallel=false;
    evalc('fit=estimate_local_quadratic(sampleFile,o);');E=fit.detail.NASA;
    assert(abs(E.actualT0-C.actualT0)<1e-8 && max(abs(E.sample.t-C.sample.t))<1e-8);
    point=zeros(1,4);lo=point;hi=point;fields={'position','velocity','acceleration'};
    for k=1:3,b=E.(fields{k});point(k)=b.endpoint;lo(k)=b.endpoint_ci_low;hi(k)=b.endpoint_ci_high;end
    b=E.averageAcceleration;point(4)=b.estimate;lo(4)=b.ci_low;hi(4)=b.ci_high;
    draws=[E.bootstrapStats(:,[1 3 5]),b.bootstrap];valid=sum(isfinite(draws),1);
   case 3
    evalc('fit=estimate_quadratic_benchmark(sampleFile,o);');E=fit.detail.NASA;
    point=zeros(1,4);lo=point;hi=point;fields={'position','velocity','acceleration'};
    for k=1:3,b=E.(fields{k});point(k)=b.endpoint;lo(k)=b.endpoint_ci_low;hi(k)=b.endpoint_ci_high;end
    point(4)=point(3);lo(4)=lo(3);hi(4)=hi(3);draws=[];valid=nan(1,4);
   case 4
    evalc('fit=estimate_broken_trend_benchmark(sampleFile,o);');E=fit.detail.NASA;
    opt=struct('noiseModel','ARMA11','B',S.B,'alpha',.1,'seed',o.seed, ...
     'windowT0',C.actualT0,'windowNames',{{'past_decade'}});
    b=bootstrap_broken_trend_arma11(E,opt);point=b.estimates;lo=b.intervals(:,1)';hi=b.intervals(:,2)';
    draws=b.draws;valid=sum(isfinite(draws),1);valid(3)=NaN;
    R.breakDate=E.breakDate;R.breakDates=b.breakDates;R.cpNoise=b.noise;R.cpCoefficients=b.coefficients;
    R.referenceSearchChecks=b.referenceSearchChecks;
  end
  if m>1,assert(abs(E.sample.t(end)-C.endTime)<1e-8);end
  keep=1:4;if m==4,keep=[1 2 4];assert(point(3)==0 && isnan(lo(3)) && isnan(hi(3)));end
  assert(all(isfinite([point,lo(keep),hi(keep)])) && all(lo(keep)<=hi(keep)));
  if m~=3,assert(all(valid(keep)>=ceil(.95*S.B)),'Fewer than 95%% valid inner draws.');end
  R.point(m,:)=point;R.low(m,:)=lo;R.high(m,:)=hi;R.draws{m}=draws;R.validDraws(m,:)=valid;
  R.success(m)=true;R.message(m)="";
 catch ME
  R.message(m)=string(getReport(ME,'extended','hyperlinks','off'));
 end
 R.seconds(m)=toc(tick);R.matlabVersion=version;save_checkpoint(file,R);clear clean
 fprintf('DGP3 COVERAGE rep %d/%d, %s: success=%d, %.1f min.\n',rep,S.coverageR,methods{m},R.success(m),R.seconds(m)/60);
end
end
