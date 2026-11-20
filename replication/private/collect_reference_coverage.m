function summary = collect_reference_coverage(study,ids)
% Recompute the published Table S8 columns from the two independent sets of intervals.
methods=["quadratic","changepoint"]; estimands=["x_T","v_T","a_T","average_a_past_decade"];
assert(~isempty(ids) && numel(unique(ids))==numel(ids));
summary=table();
for m=1:2
    design=methods(m)+"_truth"; below=false(numel(ids),4); above=below;
    hasInterval=true(numel(ids),4);
    if m==2, hasInterval(:,3)=false; end
    for k=1:numel(ids)
        folder=fullfile(study,design,'results',sprintf('rep_%04d',ids(k)));
        loaded=load(fullfile(folder,'Results.mat'),'R'); r=loaded.R;
        assert(r.rep==ids(k) && r.B==399 && r.success(m) && r.testedMethod==m);
        lo=r.low(m,:); hi=r.high(m,:); truth=r.truth(:)';
        assert(isequal(isfinite(lo)&isfinite(hi),hasInterval(k,:)) && all(isfinite(truth)));
        assert(all(lo(hasInterval(k,:))<=hi(hasInterval(k,:))));
        below(k,:)=truth<lo; above(k,:)=truth>hi;
    end
    for j=1:4
        n=sum(hasInterval(:,j)); p=NaN; se=NaN; left=NaN; right=NaN;
        if n>0
            left=mean(below(:,j)); right=mean(above(:,j));
            p=mean(~below(:,j)&~above(:,j)); se=sqrt(p*(1-p)/n);
        end
        row=table(design,methods(m),estimands(j),numel(ids),0,n,p,se,left,right, ...
            'VariableNames',{'Design','Method','Estimand','Returned','Failed','Intervals', ...
            'Coverage','CoverageMCSE','LowerAboveTruth','UpperBelowTruth'});
        summary=[summary;row]; %#ok<AGROW>
    end
end
end
