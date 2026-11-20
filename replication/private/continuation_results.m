function r = continuation_results(empirical,threshold)
% Conditional crossing of the mean record starting at each model's own endpoint.
fit = readtable(fullfile(empirical,'Empirical_Results.csv'),'TextType','string');
fit = fit(fit.Specification=="monthly" & fit.Dataset=="Mean",:);
methods = ["quadratic";"changepoint";"ss";"ss";"kernel";"kernel"];
labels = ["Quadratic continuation";"Changepoint continuation"; ...
    "Semiparametric, constant rate";"Semiparametric, constant acceleration"; ...
    "Nonparametric, constant rate";"Nonparametric, constant acceleration"];
r = table(labels,methods,'VariableNames',{'Continuation','Method'});
for j = 1:height(r)
    rows = fit(fit.Method==methods(j),:);
    assert(height(rows)==4 && numel(unique(rows.Object))==4 && numel(unique(rows.EndTime))==1);
    x = rows.Estimate(rows.Object=="x");
    v = rows.Estimate(rows.Object=="v");
    a = rows.Estimate(rows.Object=="a");
    if endsWith(labels(j),"constant rate"), a = 0; end
    endpoint = rows.EndTime(1);
    gap = threshold-x;
    discriminant = v^2 + 2*a*gap;
    year = NaN;
    if gap <= 0
        year = endpoint;
    elseif discriminant >= 0 && v+sqrt(discriminant)>0
        year = endpoint + 20*gap/(v+sqrt(discriminant));
    end
    r.Endpoint(j) = endpoint;
    r.Level(j) = x;
    r.Rate(j) = v;
    r.Acceleration(j) = a;
    r.CrossingYear(j) = year;
end
end
