function accuracy_replication(C,S,job)

d=C.selected(floor((job-1)/S.accuracyR)+1);rep=mod(job-1,S.accuracyR)+1;
file=fullfile(S.accuracyDir,sprintf('job_%05d.mat',job));
R=struct('snapshot',S.snapshot,'mode',S.mode,'job',job,'dgp',d,'rep',rep, ...
 'truth',nan(1,4),'estimates',nan(4,4),'success',false(4,1),'seconds',nan(4,1), ...
 'message',strings(4,1),'ssExitflag',NaN,'breakDate',NaN);
if isfile(file)
 L=load(file,'R');R=L.R;assert(strcmp(R.snapshot,S.snapshot)&&strcmp(R.mode,S.mode)&&R.job==job);
 if all(R.success),return;end
end
sim=simulate_temperature_record(C,d,rep);
if any(R.success),assert(max(abs(R.truth-sim.truth))<1e-10);end
R.truth=sim.truth;raw=C.raw;raw.NASA(:)=NaN;raw.NASA(C.sampleIndex)=sim.y;
work=tempname(S.accuracyDir);mkdir(work);clean=onCleanup(@()rmdir(work,'s'));
input=fullfile(work,'Sample.csv');writetable(raw,input);R.inputSHA=file_sha256(input);
methods={'ss','kernel','quadratic','changepoint'};
for m=1:4
 if R.success(m),continue;end
 tick=tic;
 try
  o=C.empiricalOptions.(methods{m});o.datasets={'NASA'};o.rebaseData=false;
  o.outputDir=fullfile(work,methods{m});o.seed=C.cfg.seed+3000000+10000*m+rep;
  switch m
   case 1
    o.bootstrap.B=0;o.bootstrap.useParallel=false;o.bootstrap.pathBands=false;
    o.averageAcceleration=struct('enabled',true,'T0',C.actualT0);
    evalc('fit=estimate_monthly_state_space(input,o);');E=fit.detail.NASA;P=E.pointEstimate;
    assert(P.exitflag>0,'SS point optimizer did not converge.');
    assert(numel(E.sampleData.t)==C.sample.T && max(abs(E.sampleData.t-C.sample.t))<1e-8);
    z=P.smoother.alphaSmooth;values=[z(1,end),10*z(2,end),100*z(3,end),100*(z(2,end)-z(2,C.index0))/10];
    R.ssExitflag=P.exitflag;R.ssTheta=P.thetaHat;R.ssLogLik=P.logLik;R.ssPath=z(1:3,:);R.pathTimes=E.sampleData.t;
   case 2
    o.kernel.averageStartYear=C.actualT0;o.bootstrap.B=0;o.bootstrap.useParallel=false;
    evalc('fit=estimate_local_quadratic(input,o);');E=fit.detail.NASA;
    values=[E.position.endpoint,E.velocity.endpoint,E.acceleration.endpoint,E.averageAcceleration.estimate];
    assert(abs(E.actualT0-C.actualT0)<1e-8);
   case 3
    evalc('fit=estimate_quadratic_benchmark(input,o);');E=fit.detail.NASA;
    values=[E.position.endpoint,E.velocity.endpoint,E.acceleration.endpoint,E.acceleration.endpoint];
   case 4
    evalc('fit=estimate_broken_trend_benchmark(input,o);');E=fit.detail.NASA;R.breakDate=E.breakDate;
    % Match the empirical bootstrap's right-hand slope convention at T0.
    av=100*E.coeff(3)*(double(C.endTime>=E.breakDate)-double(C.actualT0>=E.breakDate))/10;
    values=[E.position.endpoint,E.velocity.endpoint,0,av];
  end
  if m>1,assert(abs(E.sample.t(end)-C.endTime)<1e-8);end
  assert(all(isfinite(values)));R.estimates(m,:)=values;R.success(m)=true;R.message(m)="";
 catch ME
  R.message(m)=string(getReport(ME,'extended','hyperlinks','off'));
 end
 R.seconds(m)=toc(tick);R.matlabVersion=version;save_checkpoint(file,R);
end
fprintf('ACCURACY job %d: %s rep %d, %d/4 methods successful.\n',job,C.scenarios{d},rep,sum(R.success));
end
