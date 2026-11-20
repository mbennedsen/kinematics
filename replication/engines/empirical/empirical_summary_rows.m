function R=empirical_summary_rows(name,method,E,B)
lo=nan(1,4);hi=lo;breakDate=NaN;exitflag=NaN;validB=NaN;
if strcmp(method,'ss')
 t=E.sampleData.t;z=kinematic_path_values(E.pointEstimate.smoother.alphaSmooth);p=z(end,:);
 exitflag=E.pointEstimate.exitflag;
 if B>0
  F={'position','velocity','acceleration'};
  for j=1:3,lo(j)=E.bootstrap.(F{j}).endpoint_ci_low;hi(j)=E.bootstrap.(F{j}).endpoint_ci_high;end
  lo(4)=E.bootstrap.pathBands.low(end,4);hi(4)=E.bootstrap.pathBands.high(end,4);validB=E.bootstrap.pathBands.validDraws;
 end
 inference="Bootstrap error";
else
 t=E.sample.t;p=[E.position.endpoint E.velocity.endpoint NaN NaN];
 lo(1:2)=[E.position.endpoint_ci_low E.velocity.endpoint_ci_low];hi(1:2)=[E.position.endpoint_ci_high E.velocity.endpoint_ci_high];
 if strcmp(method,'changepoint')
  breakDate=E.breakDate;g=zeros(numel(E.coeff),1);g(3)=10*(double(t(end)>breakDate)-double(t(end-120)>breakDate));
  p(4)=g'*E.coeff;se=sqrt(max(0,g'*E.Vhac*g));
  if any(g),lo(4)=p(4)-1.64485362695147*se;hi(4)=p(4)+1.64485362695147*se;end
  if abs(t(end)-breakDate)>1e-9,p(3)=0;end
  inference="HAC conditional on selected break";B=0;
 else
  p(3)=E.acceleration.endpoint;lo(3)=E.acceleration.endpoint_ci_low;hi(3)=E.acceleration.endpoint_ci_high;
  if strcmp(method,'quadratic')
   p(4)=p(3);lo(4)=lo(3);hi(4)=hi(3);B=0;inference="Newey-West HAC";
  else
   assert(abs(E.actualT0-t(end-120))<1e-8);p(4)=E.averageAcceleration.estimate;
   lo(4)=E.averageAcceleration.ci_low;hi(4)=E.averageAcceleration.ci_high;
   validB=sum(isfinite(E.averageAcceleration.bootstrap));inference="Moving-block bootstrap percentile";
  end
 end
end
assert(numel(t)>120 && abs(t(end)-t(end-120)-10)<1e-8 && all(isfinite(p)));
R=table(repmat(string(name),4,1),repmat(string(method),4,1),["x";"v";"a";"average_a"], ...
 p',lo',hi',repmat(t(1),4,1),repmat(t(end),4,1),repmat(t(end-120),4,1), ...
 repmat(B,4,1),repmat(validB,4,1),repmat(breakDate,4,1),repmat(exitflag,4,1),repmat(inference,4,1), ...
 'VariableNames',{'Dataset','Method','Object','Estimate','CILow','CIHigh','StartTime','EndTime','ActualT0','BootstrapB','ValidDraws','BreakDate','OptimizerExitFlag','Inference'});
end
