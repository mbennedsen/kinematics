function validate_package(mode)
% VALIDATE_PACKAGE('saved'): render and check all 36 published items.
% VALIDATE_PACKAGE('smoke'): also test clean-workspace setup, DGP reconstruction,
% and fresh kernel/quadratic fits (3 bootstrap draws; not production intervals).
if nargin<1,mode='saved';end
mode=validatestring(mode,{'saved','smoke'});
root=fileparts(fileparts(mfilename('fullpath')));old=path;c=onCleanup(@()path(old));
verify_snapshot;
if strcmp(mode,'saved')
 run_all('saved');
 out=fullfile(root,'outputs','reproduced','saved');
 [items,~]=item_catalog;
 for k=1:numel(items)
  item=items{k};
  if startsWith(item,'Table'), formats={'csv','mat','tex','html'};
  else, formats={'fig','pdf','png'}; end
  for j=1:numel(formats)
   written=dir(fullfile(out,[item '.' formats{j}]));
   assert(isscalar(written) && written.bytes>0,'Missing output for %s.',item);
  end
 end
 % Independent manuscript values at its printed precision.
 for k=2:3
  tableData=readtable(fullfile(out,sprintf('TableS%d.csv',k)));
  dates=[2030.9 2028.7 2028.5 2028.4 2027.5 2027.5;2046.7 2042.9 2039.7 2037.5 2037.5 2035.7];
  assert(max(abs(tableData.CrossingYear(:)'-dates(k-1,:)))<.05);
 end
 for k=7:8
  tableData=readtable(fullfile(out,sprintf('TableS%d.csv',k)),'TextType','string');
  rate=tableData(tableData.Estimand=="v_T",:);
  methods=["quadratic","changepoint"]; expected=[0 .495;.79 .8475];
  for j=1:2
   assert(abs(rate.Coverage(rate.Method==methods(j))-expected(k-6,j))<1e-10);
  end
 end
end
if strcmp(mode,'smoke')
 opt=struct('Workers',1,'RunID',['validation_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
 FigureS1('reestimate',opt);FigureS10('reestimate',opt);TableS4('reestimate',opt);
 run=fullfile(root,'outputs','reestimated',opt.RunID);work=fullfile(run,'work');
 eng=fullfile(work,'empirical');
 addpath(eng,fullfile(work,'methods'),'-begin');prepare_empirical;
 R=run_empirical(5,3,1,{'quadratic','kernel'});
 collect_empirical(eng,fullfile(run,'smoke_exports'),'monthly','ERA5',{'quadratic','kernel'},3);
 saved=readtable(fullfile(root,'outputs','saved','empirical','Empirical_Results.csv'),'TextType','string');
 for j=1:height(R)
  z=saved(saved.Dataset==R.Dataset(j)&saved.Method==R.Method(j)&saved.Object==R.Object(j)&saved.Specification=="monthly",:);
  assert(height(z)==1&&abs(z.Estimate-R.Estimate(j))<1e-7,'Fresh point estimate does not match paper.');
 end
 a=readtable(fullfile(run,'results','simulations','accuracy','DGPTruth.csv'),'TextType','string');
 b=readtable(fullfile(root,'outputs','saved','simulations','accuracy','DGPTruth.csv'),'TextType','string');
 assert(isequal(a.DGP,b.DGP)&&max(abs(a{:,2:5}-b{:,2:5}),[],'all')<1e-10,'DGP truth mismatch.');
 smooth=fullfile(work,'smooth_coverage');
 addpath(smooth,'-begin');
 if ~isfolder(fullfile(smooth,'validation')), mkdir(fullfile(smooth,'validation')); end
 validate_smooth_design;
 study=fullfile(work,'reference_coverage');
 addpath(study,'-begin');
 for design=["quadratic_truth","changepoint_truth"]
  run_parametric_reference_study(design,1,1);
 end
 reference=collect_reference_coverage(study,1);
 assert(height(reference)==8 && sum(reference.Intervals)==7);
 fprintf('SMOKE PASS: both fixed parametric paths and DGP 3 construction.\n');
 fprintf('SMOKE PASS: clean setup, deterministic designs, and fresh kernel/quadratic point estimates.\n');
 fprintf('Three-draw smoke intervals are not production intervals.\n');
end
fprintf('VALIDATION PASS: %s\n',mode);
end
