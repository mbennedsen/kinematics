function write_diagnostics(inputFile,outFile)
s=load(inputFile,'detail');names={'NASA','NOAA','HadCRU','Berk','ERA5'};rows=cell(5,1);
for i=1:5
 d=innovation_diagnostics_local(s.detail.(names{i}).pointEstimate.filter.vStd,struct('lbLags',[1 2 6 12 24],'archLags',12));
 r=struct('Dataset',string(names{i}),'InnovMean',d.mean,'InnovStd',d.std,'InnovSkew',d.skew,'InnovKurt',d.kurt,'JB_p',d.JB_p,'AR1_Coeff',d.AR1_coeff,'AR1_p',d.AR1_p,'DW',d.DW);
 for h=[1 2 6 12 24],for prefix={'LB','LBsq'},key=sprintf('%s%d_p',prefix{1},h);r.(key)=d.(key);end,end
 r.ARCH_p=d.ARCH_p;rows{i}=r;
end
writetable(struct2table(vertcat(rows{:})),outFile);
end
function out = innovation_diagnostics_local(x, diagOpts)

    x = x(:);
    xObs = x(isfinite(x));

    out = struct();
    out.nObserved = numel(xObs);
    out.mean = nanmean_local(xObs);
    out.std = nanstd_local(xObs);
    out.skew = skewness_local(xObs);
    out.kurt = kurtosis_local(xObs);
    [out.JB_stat, out.JB_p] = jarque_bera_local(xObs);
    [out.AR1_coeff, out.AR1_p] = ar1_test_local(x);
    out.DW = durbin_watson_local(x);

    for h = diagOpts.lbLags(:)'
        [Q, p] = ljung_box_local(x, h);
        out.(sprintf('LB%d_stat', h)) = Q;
        out.(sprintf('LB%d_p', h)) = p;
        [Q2, p2] = ljung_box_local(x.^2, h);
        out.(sprintf('LBsq%d_stat', h)) = Q2;
        out.(sprintf('LBsq%d_p', h)) = p2;
    end

    [out.ARCH_stat, out.ARCH_p] = arch_lm_local(x, diagOpts.archLags);
end

%% ========================================================================
function [stat, p] = jarque_bera_local(x)
    x = x(isfinite(x));
    n = numel(x);
    if n < 8
        stat = NaN; p = NaN; return;
    end
    s = skewness_local(x);
    k = kurtosis_local(x);
    stat = n/6 * (s^2 + 0.25*(k - 3)^2);
    p = 1 - chi2cdf_local(stat, 2);
end

%% ========================================================================
function [phi, p] = ar1_test_local(x)

    x0 = x(1:end-1);
    x1 = x(2:end);
    idx = isfinite(x0) & isfinite(x1);
    y = x1(idx);
    z = x0(idx);

    n = numel(y);
    if n < 10
        phi = NaN; p = NaN; return;
    end

    X = [ones(n,1), z];
    b = X \ y;
    e = y - X * b;
    s2 = sum(e.^2) / max(n - size(X,2), 1);
    XXinv = pinv(X' * X);
    se = sqrt(diag(s2 * XXinv));

    phi = b(2);
    zstat = phi / max(se(2), eps);
    p = 2 * (1 - normcdf_local(abs(zstat)));
end

%% ========================================================================
function dw = durbin_watson_local(x)

    x = x(:);
    idx = isfinite(x(2:end)) & isfinite(x(1:end-1));
    d = x(2:end) - x(1:end-1);
    num = sum(d(idx).^2);
    den = sum(x(isfinite(x)).^2);
    if den <= 0
        dw = NaN;
    else
        dw = num / den;
    end
end

%% ========================================================================
function [Q, p] = ljung_box_local(x, m)

    x = x(:);
    xObs = x(isfinite(x));
    n = numel(xObs);
    if n <= m + 1
        Q = NaN; p = NaN; return;
    end

    rho = nan(m, 1);
    mu = mean(xObs);

    for h = 1:m
        xLead = x(1+h:end);
        xLag  = x(1:end-h);
        idx = isfinite(xLead) & isfinite(xLag);
        if sum(idx) < 3
            continue;
        end
        a = xLead(idx) - mu;
        b = xLag(idx) - mu;
        den = sqrt(sum(a.^2) * sum(b.^2));
        if den > 0
            rho(h) = sum(a .* b) / den;
        end
    end

    good = isfinite(rho);
    if ~any(good)
        Q = NaN; p = NaN; return;
    end

    hVec = (1:m)';
    Q = n * (n + 2) * sum((rho(good).^2) ./ max(n - hVec(good), 1));
    p = 1 - chi2cdf_local(Q, sum(good));
end

%% ========================================================================
function [lmStat, p] = arch_lm_local(x, q)

    x = x(:);
    y = x.^2;
    [Y, X] = build_lag_matrix_local(y, q);

    nEff = size(Y, 1);
    if nEff <= q + 5
        lmStat = NaN; p = NaN; return;
    end

    b = X \ Y;
    e = Y - X * b;
    sse = sum(e.^2);
    sst = sum((Y - mean(Y)).^2);

    if sst <= 0
        lmStat = NaN; p = NaN; return;
    end

    R2 = 1 - sse / sst;
    lmStat = nEff * R2;
    p = 1 - chi2cdf_local(lmStat, q);
end

%% ========================================================================
function [Y, X] = build_lag_matrix_local(x, q)

    n = numel(x);
    rows = false(n, 1);
    Xfull = nan(n, q + 1);

    for t = (q+1):n
        vals = x((t-q):t);
        if all(isfinite(vals))
            rows(t) = true;
            Xfull(t, :) = [1, x(t-1:-1:t-q)'];
        end
    end

    Y = x(rows);
    X = Xfull(rows, :);
end

%% ========================================================================
function S = regularize_cov_local(S)

    S = symmetrize_local(S);
    jitter = 1e-10;
    for k = 1:8
        S = S + jitter * eye(size(S));
        [~, p] = chol(S);
        if p == 0
            return;
        end
        jitter = jitter * 10;
    end
end

%% ========================================================================
function A = symmetrize_local(A)
    A = 0.5 * (A + A');
end

%% ========================================================================
function x = coerce_numeric_local(x)

    if isnumeric(x)
        x = x(:); return;
    end
    if islogical(x)
        x = double(x(:)); return;
    end
    if iscell(x)
        x = string(x);
    end
    if isstring(x) || ischar(x) || iscategorical(x)
        x = str2double(string(x));
        x = x(:); return;
    end
    error('Could not coerce column to numeric.');
end

%% ========================================================================
function m = nanmean_local(x)
    x = x(isfinite(x));
    if isempty(x), m = NaN; else, m = mean(x); end
end

%% ========================================================================
function s = nanstd_local(x)
    x = x(isfinite(x));
    if numel(x) < 2, s = NaN; else, s = std(x, 0); end
end

%% ========================================================================
function v = nanvar_local(x)
    x = x(isfinite(x));
    if numel(x) < 2, v = NaN; else, v = var(x, 0); end
end

%% ========================================================================
function s = skewness_local(x)
    x = x(isfinite(x));
    if numel(x) < 3, s = NaN; return; end
    mu = mean(x); sd = std(x, 0);
    if sd <= 0, s = NaN; return; end
    s = mean(((x - mu) / sd).^3);
end

%% ========================================================================
function k = kurtosis_local(x)
    x = x(isfinite(x));
    if numel(x) < 4, k = NaN; return; end
    mu = mean(x); sd = std(x, 0);
    if sd <= 0, k = NaN; return; end
    k = mean(((x - mu) / sd).^4);
end

%% ========================================================================
function p = normcdf_local(z)
    p = 0.5 * erfc(-z ./ sqrt(2));
end

%% ========================================================================
function y = chi2cdf_local(x, k)
    y = gammainc(max(x, 0) / 2, k / 2, 'lower');
end

%% ========================================================================
function m = nanmean_rows_local(X)

    m = nan(size(X,1),1);
    for i = 1:size(X,1)
        xi = X(i,:);
        xi = xi(isfinite(xi));
        if isempty(xi)
            m(i) = NaN;
        else
            m(i) = mean(xi);
        end
    end
end

%% ========================================================================
function s = nansum_local(x)
    x = x(isfinite(x));
    if isempty(x), s = NaN; else, s = sum(x); end
end

%% ========================================================================
function out = merge_structs_local(a, b)
% Recursive merge: fields in b overwrite fields in a.

    out = a;
    if isempty(b)
        return;
    end

    bFields = fieldnames(b);
    for i = 1:numel(bFields)
        f = bFields{i};
        if isstruct(b.(f)) && isfield(out, f) && isstruct(out.(f))
            out.(f) = merge_structs_local(out.(f), b.(f));
        else
            out.(f) = b.(f);
        end
    end
end
