function [r,P]=paper_table(item,source)
% Numerical tables and display panels shared by saved and re-estimated modes.
n=str2double(extractAfter(item,'TableS'));e=fullfile(source,'empirical');
rd=@(file)readtable(file,'TextType','string');
ds=["NASA","NOAA","HadCRU","Berk","ERA5"];ms=["ss","kernel","quadratic","changepoint"];
ob=["x","v","a","average_a"];ot=["Warming level (degrees C)","Warming rate (degrees C / decade)","Acceleration (degrees C / decade^2)","Past-decade average acceleration (degrees C / decade^2)"];
P=cell(0,3);
switch n
 case {1,9,16}
  if n==1,r=rd(fullfile(e,'Empirical_Results.csv'));r=r(r.Specification=="monthly",:);names=[ds,"Mean"];
  elseif n==9,r=rd(fullfile(e,'Annual_Results.csv'));names=ds;
  else,r=rd(fullfile(source,'no_covariates','Endpoint_Summary.csv'));names=[ds,"Mean"];ms="ss";end
  for k=1:4
   z=r(r.Object==ob(k),:);P(end+1,:)={ot(k),method_header(ms),method_rows(z,names,ms,n~=16)}; %#ok<AGROW>
  end
 case {2,3}
  threshold = 1.5;
  if n == 3, threshold = 2; end
  r = continuation_results(e,threshold);
  z = r(:,{'Continuation','Level','Rate','Acceleration','CrossingYear'});
  rows = tabcells(z,3);
  for j = 1:height(z), rows{j,5} = fmt(z.CrossingYear(j),1); end
  P = {"",{'Continuation','Level','Rate','Acceleration','Crossing year'},rows};
 case 4
  r=rd(fullfile(source,'simulations','accuracy','DGPTruth.csv'));
  P={"",{'DGP','Level','Rate','Acceleration','Average acceleration'},tabcells(r,3)};
 case {5,6}
  r=rd(fullfile(source,'simulations','accuracy','AccuracySummary.csv'));assert(all(r.Valid==1000),'Accuracy table requires all 1000 replications per DGP.');
  keys=["x_T","v_T","a_T","average_a_past_decade"];dgps=["zero","ramp","NASA_smoothed","NASA_quadratic","NASA_changepoint"];
  for k=1:4
   rows=cell(5,5);
   for j=1:5
    rows{j,1}=char(dgps(j));v=r(r.DGP==dgps(j)&r.Estimand==keys(k),:);
    for m=1:4
     z=one(v,'Method',ms(m));
     if n==5,rows{j,m+1}=fmt(z.RMSE);if z.RMSE==min(v.RMSE),rows{j,m+1}=['*' rows{j,m+1}];end
     else,rows{j,m+1}=[fmt(z.Bias) ' (' fmt(z.BiasMCSE) ')'];end
    end
   end
   P(end+1,:)={ot(k),{'DGP','Semiparametric','Nonparametric','Quadratic','Changepoint'},rows}; %#ok<AGROW>
  end
 case {7,8}
  if n == 7
   r = rd(fullfile(source,'simulations','smooth_coverage','CoverageSummary.csv'));
   assert(all(r.Returned==400 & r.B==399) && height(r)==16);
   count = r.ValidIntervals; methods = ms;
   lo = 'LowerMiss'; hi = 'UpperMiss';
  else
   r = rd(fullfile(source,'simulations','reference_coverage','Summary.csv'));
   assert(all(r.Returned==400 & r.Failed==0) && height(r)==8);
   count = r.Intervals; methods = ["quadratic","changepoint"];
   lo = 'LowerAboveTruth'; hi = 'UpperBelowTruth';
   r = r(:,{'Design','Method','Estimand','Returned','Failed','Intervals','Coverage','CoverageMCSE',lo,hi});
  end
  structural = r.Method=="changepoint" & r.Estimand=="a_T";
  assert(sum(structural)==1 && all(count(~structural)==400) && count(structural)==0);
  for m=methods
   z=r(r.Method==m,{'Estimand','Coverage','CoverageMCSE',lo,hi});
   assert(height(z)==4 && numel(unique(z.Estimand))==4);
   P(end+1,:)={method_name(m),{'Estimand','Coverage','MC s.e.','Below','Above'},tabcells(z,3)};
  end
 case {10,11,12,13}
  r=rd(fullfile(e,'Sensitivity_Results.csv'));
  switch n
   case 10,m="quadratic";variants="start_"+[1950 1960 1970 1980 1990];heads=string([1950 1960 1970 1980 1990]);objs=["x","v","a"];
   case 11,m="changepoint";variants="start_"+[1950 1960 1970 1980 1990];heads=string([1950 1960 1970 1980 1990]);objs=["x","v","slope","break"];
   case 12,m="kernel";variants="bandwidth_"+[15 20 25 30 35];heads=string([15 20 25 30 35]);objs=ob;
   case 13,m="ss";variants="sigma_"+[.5 .75 1 1.5 2];heads=string([.5 .75 1 1.5 2]);objs=ob;
  end
  r=r(r.Method==m,{'Dataset','Method','Object','Estimate','Variant'});
  if n==11,z=rd(fullfile(e,'ChangepointSensitivity.csv'));r=[r;z(:,r.Properties.VariableNames)];end
  r=r(ismember(r.Object,objs),:);
  for k=1:numel(objs)
   cross=objs(k)~="break";rows=cell(5+cross,6);
   for j=1:5
    rows{j,1}=dataset_name(ds(j));
    for h=1:5
     z=r(r.Dataset==ds(j)&r.Object==objs(k)&r.Variant==variants(h),:);assert(height(z)==1);
     if cross,rows{j,h+1}=fmt(z.Estimate);else,y=floor(z.Estimate);mo=round((z.Estimate-y)*12+.5);rows{j,h+1}=sprintf('%d-%02d',y,mo);end
    end
   end
   if cross,rows{6,1}='Cross-dataset mean';for h=1:5,z=r(r.Object==objs(k)&r.Variant==variants(h),:);rows{6,h+1}=fmt(mean(z.Estimate));end,end
   title=objs(k);if ismember(title,ob),title=ot(ob==title);end
   P(end+1,:)={title,[{'Dataset'},cellstr(heads)],rows}; %#ok<AGROW>
  end
 case 14
  r=rd(fullfile(e,'reporting','VelocityModel_Comparison.csv'));assert(all(r.OptimizerExitFlag>0));
  rows=cell(6,5);for j=1:5,rows(j,:)=[{dataset_name(r.Dataset(j))},arrayfun(@fmt,r{j,2:5},'UniformOutput',false)];end
  rows(6,:)=[{'Cross-dataset mean'},arrayfun(@fmt,mean(r{:,2:5},1),'UniformOutput',false)];
  P={"",{'Dataset','Baseline level','Alternative level','Baseline rate','Alternative rate'},rows};
 case 15
  r=rd(fullfile(e,'reporting','AlternativeWindows_AllMethods.csv'));
  for year=[2010 2012 2014 2015 2016 2018 2020],P(end+1,:)={"Starting January "+year,method_header(ms),method_rows(r(r.Year==year,:),ds,ms,true)};end %#ok<AGROW>
 case {17,18,19,20,21,22}
  p=fullfile(e,'validated_monthly','parameters');
  switch n
   case 17
    r=rd(fullfile(p,'SS_ParameterAppendix_Dynamics.csv'));
    z=r(:,{'Dataset','q','MonthlyAccelerationShockSD_decade2','ARMA_AR1','ARMA_MA1','BaseResidualVariance','BaseResidualSD'});
    heads={'Dataset','sigma squared','Monthly shock SD','AR(1)','MA(1)','Residual variance','Residual SD'};rows=tabcells(z,4);
    for j=1:height(z),rows{j,2}=sprintf('%.3e',z.q(j));rows{j,6}=sprintf('%.3e',z.BaseResidualVariance(j));end
   case {18,21}
    suffix='BetaBasis';if n==21,suffix='ImpliedLagWeights';end
    r=rd(fullfile(p,['SS_ParameterAppendix_' suffix '.csv']));z=r(:,[{'Label'},cellstr(ds)]);heads=[{'Parameter'},cellstr(ds)];rows=tabcells(z,4);
    for j=1:height(z),for k=2:width(z),x=z{j,k};if abs(x)>0&&abs(x)<.00005,rows{j,k}=sprintf('%.3e',x);end,end,end
   case {19,20}
    r=rd(fullfile(p,'SS_ParameterAppendix_MonthlyVarianceScales.csv'));suffix="_VarianceScale";digits=3;if n==20,suffix="_ResidualVariance";digits=4;end
    z=r(:,[{'MonthName'},cellstr(ds+suffix)]);heads=[{'Month'},cellstr(ds)];rows=tabcells(z,digits);
   case 22
    r=rd(fullfile(p,'ResidualDiagnostics.csv'));
    cols={'Dataset','InnovMean','InnovStd','InnovSkew','InnovKurt','JB_p','AR1_Coeff','AR1_p','DW','LB1_p','LB2_p','LB6_p','LB12_p','LB24_p','LBsq1_p','LBsq2_p','LBsq6_p','LBsq12_p','LBsq24_p','ARCH_p'};
    z=r(:,cols);heads=cols;rows=tabcells(z,3);
    for j=1:height(z),for k=2:width(z),if endsWith(cols{k},'_p')&&z{j,k}<1e-4,rows{j,k}='<0.0001';end,end,end
  end
  P={"",heads,rows};
 case 23
  r=rd(fullfile(e,'Changepoint_Diagnostics.csv'));
  for kind=["monthly","annual"]
   z=r(r.Specification==kind,:);rows=cell(height(z),6);
   for j=1:height(z)
    rows(j,:)={dataset_name(z.Dataset(j)),fmt(z.AR1(j)),fmt(z.MA1(j)),fmt(z.InnovationSD(j)), ...
     ['[' fmt(z.BreakDateLow(j),2) ', ' fmt(z.BreakDateHigh(j),2) ']'],fmt(100*z.BreakInsideDecadeFraction(j),1)};
   end
   P(end+1,:)={kind,{'Dataset','AR(1)','MA(1)','Innovation SD','90% breakpoint range','Inside decade (%)'},rows}; %#ok<AGROW>
  end
 otherwise,error('Unknown table.');
end
end
function heads=method_header(ms)
heads={'Dataset'};for m=ms,heads=[heads,{char(method_name(m))+" estimate",char(method_name(m))+" 90% interval"}];end
end
function rows=method_rows(T,names,ms,cross)
rows=cell(numel(names)+cross,1+2*numel(ms));
for j=1:numel(names)
 rows{j,1}=dataset_name(names(j));
 for k=1:numel(ms)
  z=T(T.Dataset==names(j)&T.Method==ms(k),:);assert(height(z)==1);
  rows{j,2*k}=fmt(z.Estimate);rows{j,2*k+1}=interval(z);
 end
end
if cross
 rows{end,1}='Cross-dataset mean';
 for k=1:numel(ms)
  z=T(T.Method==ms(k)&ismember(T.Dataset,["NASA","NOAA","HadCRU","Berk","ERA5"]),:);assert(height(z)==5);
  rows{end,2*k}=fmt(mean(z.Estimate));rows{end,2*k+1}='--';
 end
end
end
function x=one(T,key,value),x=T(T.(key)==value,:);assert(height(x)==1);end
function rows=tabcells(T,d)
rows=cell(height(T),width(T));
for j=1:height(T),for k=1:width(T),v=T{j,k};if isnumeric(v)||islogical(v),rows{j,k}=fmt(v,d);else,rows{j,k}=char(string(v));end,end,end
end
function s=fmt(v,d)
if nargin<2,d=3;end
if ~isfinite(v),s='--';return;end
v=round(v,d);if v==0,v=0;end;s=sprintf('%.*f',d,v);
end
function s=interval(z)
if ~isfinite(z.CILow)||~isfinite(z.CIHigh),s='--';return;end
lo=fmt(z.CILow);hi=fmt(z.CIHigh);
if abs(z.CILow)>0&&abs(z.CILow)<.0005,lo=fmt(z.CILow,4);end
if abs(z.CIHigh)>0&&abs(z.CIHigh)<.0005,hi=fmt(z.CIHigh,4);end
s=['[' lo ', ' hi ']'];
end
function s=method_name(m)
keys=["ss","kernel","quadratic","changepoint"];labels=["Semiparametric","Nonparametric","Quadratic","Changepoint"];s=labels(keys==m);
end
function s=dataset_name(d)
keys=["NASA","NOAA","HadCRU","Berk","ERA5","Mean"];labels=["NASA","NOAA","HadCRUT5","Berkeley Earth","ERA5","Arithmetic mean"];s=char(labels(keys==d));
end
