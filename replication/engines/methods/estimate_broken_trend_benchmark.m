function out = estimate_broken_trend_benchmark(dataFile, userOpts)
% ESTIMATE_BROKEN_TREND_BENCHMARK
% Parametric broken-trend benchmark for the temperature-kinematics paper.
%
% Model:
%   y_t = alpha + beta * ttil + dBeta * max(ttil - tc_ttil, 0) + gamma' z_t + u_t
%
% where tc is selected by minimizing the residual sum of squares (SSE)
% over an admissible grid of candidate break dates.
%
% Notes:
%   - Uses the same Almon distributed-lag covariate structure as the SS model
%   - Uses Newey-West HAC inference conditional on the selected break date
%   - Reports:
%         (i)   endpoint and trailing-average position x(t)
%         (ii)  endpoint and trailing-average warming rate v(t)
%         (iii) discrete slope change Delta v at the break
%   - For a broken-trend model, smooth local acceleration a(t) is zero away
%     from the break and undefined at the break. Hence this script reports
%     the slope change instead of an "acceleration" table.
%
% Output CSV files:
%   BrokenTrend_Table1_Position.csv
%   BrokenTrend_Table2_Velocity.csv
%   BrokenTrend_Table3_SlopeChange.csv
%   BrokenTrend_Table4_BreakDate.csv

    assert(nargin >= 1 && ~isempty(dataFile),'Provide the path to the analysis panel.');
    if nargin < 2
        userOpts = struct();
    end

    opts = default_options_local();
    opts = merge_structs_local(opts, userOpts);

    if ~exist(dataFile, 'file')
        error('Data file not found: %s', dataFile);
    end
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end

    raw = readtable(dataFile);
    data = prepare_base_data_local(raw, opts);

    datasets = opts.datasets(:)';
    nD = numel(datasets);

    posRows   = cell(nD, 1);
    velRows   = cell(nD, 1);
    jumpRows  = cell(nD, 1);
    breakRows = cell(nD, 1);
    detail = struct();

    fprintf('\n============================================================\n');
    fprintf('Estimating broken-trend benchmark model\n');
    fprintf('============================================================\n');
    fprintf('Window start year: %.0f\n', opts.broken.windowStartYear);
    fprintf('Minimum segment length: %d months\n', opts.broken.minSegmentMonths);
    fprintf('Candidate break step: %d months\n', opts.broken.candidateStepMonths);
    fprintf('HAC lag truncation L: %d months\n', opts.broken.hacLagMonths);

    for i = 1:nD
        datasetName = datasets{i};
        fprintf('\n------------------------------------------------------------\n');
        fprintf('Dataset: %s\n', datasetName);
        fprintf('------------------------------------------------------------\n');

        sample = build_sample_local(data, datasetName, opts);
        est = fit_broken_trend_model_local(sample, opts);

        fprintf('  break date = %s\n', est.breakDateString);
        fprintf('  x(T)       = %10.4f\n', est.position.endpoint);
        fprintf('  xbar(60)   = %10.4f\n', est.position.trailing60);
        fprintf('  v(T)       = %10.4f degC/decade\n', est.velocity.endpoint);
        fprintf('  vbar(60)   = %10.4f degC/decade\n', est.velocity.trailing60);
        fprintf('  Delta v    = %10.4f degC/decade\n', est.slopeChange.estimate);

        posRows{i}   = make_row_local(datasetName, est.position);
        velRows{i}   = make_row_local(datasetName, est.velocity);
        jumpRows{i}  = make_jump_row_local(datasetName, est.slopeChange);
        breakRows{i} = make_break_row_local(datasetName, est);

        detail.(datasetName) = est;
    end

    posTable   = struct2table(vertcat(posRows{:}));
    velTable   = struct2table(vertcat(velRows{:}));
    jumpTable  = struct2table(vertcat(jumpRows{:}));
    breakTable = struct2table(vertcat(breakRows{:}));

    writetable(posTable,   fullfile(opts.outputDir, 'BrokenTrend_Table1_Position.csv'));
    writetable(velTable,   fullfile(opts.outputDir, 'BrokenTrend_Table2_Velocity.csv'));
    writetable(jumpTable,  fullfile(opts.outputDir, 'BrokenTrend_Table3_SlopeChange.csv'));
    writetable(breakTable, fullfile(opts.outputDir, 'BrokenTrend_Table4_BreakDate.csv'));
    save(fullfile(opts.outputDir, 'BrokenTrend_Benchmark_Output.mat'), 'detail', 'opts', '-v7.3');

    out = struct();
    out.options = opts;
    out.positionTable = posTable;
    out.velocityTable = velTable;
    out.slopeChangeTable = jumpTable;
    out.breakTable = breakTable;
    out.detail = detail;

    fprintf('\nSaved broken-trend benchmark output to:\n  %s\n', opts.outputDir);
end

%% ========================================================================
function opts = default_options_local()

    opts = struct();
    opts.outputDir = fullfile(pwd, 'results_broken_trend_benchmark');
    opts.datasets = {'NASA','NOAA','HadCRU','Berk','ERA5'};
    opts.rebaseData = true;

    opts.covariates = struct();
    opts.covariates.names = {'NINO34','volc2','ssn'};
    opts.covariates.standardize = true;
    opts.covariates.maxLags = [12 12 12];
    opts.covariates.almonDegrees = [2 1 1];

    opts.broken = struct();
    opts.broken.windowStartYear = 1970;   % More natural for a broken-trend benchmark
    opts.broken.trailingMonths = 60;
    opts.broken.hacLagMonths = 24;
    opts.broken.minSegmentMonths = 60;    % Minimum observations on each side of break
    opts.broken.candidateStepMonths = 1;  % Search every month

    opts.inference = struct();
    opts.inference.ciLevel = 0.95;
end

%% ========================================================================
function data = prepare_base_data_local(raw, opts)

    data = struct();
    data.t = coerce_numeric_column_local(raw.t);

    for i = 1:numel(opts.datasets)
        nm = opts.datasets{i};
        data.temp.(nm) = coerce_numeric_column_local(raw.(nm));
    end

    if opts.rebaseData
        data.temp = rebase_temperature_series_local(data.temp, data.t);
    end

    covNames = opts.covariates.names;
    for i = 1:numel(covNames)
        x = coerce_numeric_column_local(raw.(covNames{i}));
        if opts.covariates.standardize
            mu = nanmean_local(x);
            sd = nanstd_local(x);
            if ~(isfinite(sd) && sd > 0)
                sd = 1;
            end
            data.covariates.(covNames{i}) = (x - mu) ./ sd;
        else
            data.covariates.(covNames{i}) = x;
        end
    end
end

%% ========================================================================
function sample = build_sample_local(data, datasetName, opts)

    y = data.temp.(datasetName);
    [Z, zNames] = build_almon_design_local(data.covariates, opts.covariates.names, ...
        opts.covariates.maxLags, opts.covariates.almonDegrees);

    valid = isfinite(y) & all(isfinite(Z), 2) & (data.t >= opts.broken.windowStartYear);
    idx = find(valid);

    if isempty(idx)
        error('No valid observations for dataset %s in broken-trend window.', datasetName);
    end

    use = idx(1):idx(end);

    sample = struct();
    sample.name = datasetName;
    sample.t = data.t(use);
    sample.y = y(use);
    sample.Z = Z(use, :);
    sample.ZNames = zNames;
    sample.n = numel(sample.y);

    if sample.n < max(opts.broken.trailingMonths, 2 * opts.broken.minSegmentMonths + 1)
        error('Broken-trend sample too short for dataset %s.', datasetName);
    end
end

%% ========================================================================
function [Z, zNames] = build_almon_design_local(covStruct, names, maxLags, degrees)

    nCov = numel(names);
    if numel(maxLags) ~= nCov || numel(degrees) ~= nCov
        error('maxLags and almonDegrees must have one entry per covariate.');
    end

    n = numel(covStruct.(names{1}));
    totalCols = sum(degrees + 1);
    Z = nan(n, totalCols);
    zNames = cell(1, totalCols);

    k = 0;
    for i = 1:nCov
        x = covStruct.(names{i});
        L = maxLags(i);
        d = degrees(i);

        lagged = nan(n, L + 1);
        for ell = 0:L
            lagged(:, ell + 1) = lag_series_local(x, ell);
        end

        lagIdx = (0:L)';
        for p = 0:d
            k = k + 1;
            Z(:, k) = lagged * (lagIdx .^ p);
            zNames{k} = sprintf('%s_AlmonDeg%d', names{i}, p);
        end
    end
end

%% ========================================================================
function xLag = lag_series_local(x, L)
    xLag = nan(size(x));
    if L == 0
        xLag = x;
    else
        xLag((L+1):end) = x(1:(end-L));
    end
end

%% ========================================================================
function est = fit_broken_trend_model_local(sample, opts)

    n = sample.n;
    tRef = mean(sample.t);
    ttil = sample.t - tRef;

    minSeg = opts.broken.minSegmentMonths;
    step = opts.broken.candidateStepMonths;

    candIdx = (minSeg + 1):step:(n - minSeg);
    nCand = numel(candIdx);

    if nCand < 1
        error('No admissible candidate break dates for dataset %s.', sample.name);
    end

    sse = nan(nCand, 1);
    bAll = cell(nCand, 1);

    for j = 1:nCand
        jc = candIdx(j);
        tcTil = ttil(jc);
        Xtrend = [ones(n, 1), ttil, max(0, ttil - tcTil)];
        X = [Xtrend, sample.Z];

        b = X \ sample.y;
        u = sample.y - X * b;

        sse(j) = sum(u.^2);
        bAll{j} = b;
    end

    [~, jStar] = min(sse);
    idxStar = candIdx(jStar);
    tc = sample.t(idxStar);
    tcTil = ttil(idxStar);

    Xtrend = [ones(n, 1), ttil, max(0, ttil - tcTil)];
    X = [Xtrend, sample.Z];
    b = X \ sample.y;
    u = sample.y - X * b;

    Vhac = newey_west_local(X, u, opts.broken.hacLagMonths);

    % Endpoint and trailing objects
    tT = ttil(end);
    tailIdx = (n - opts.broken.trailingMonths + 1):n;
    tTail = ttil(tailIdx);

    hingeT = max(0, tT - tcTil);
    hingeTail = max(0, tTail - tcTil);

    % x(T) and trailing-average x
    g_xT = [1; tT; hingeT; zeros(size(sample.Z, 2), 1)];
    g_xbar = [1; mean(tTail); mean(hingeTail); zeros(size(sample.Z, 2), 1)];

    % v(T) and trailing-average v
    % v(t) = beta + dBeta * 1{t > tc}
    iT = double(tT > tcTil);
    iTail = double(tTail > tcTil);

    g_vT = [0; 1; iT; zeros(size(sample.Z, 2), 1)];
    g_vbar = [0; 1; mean(iTail); zeros(size(sample.Z, 2), 1)];

    % Discrete slope change Delta v = dBeta
    g_jump = [0; 0; 1; zeros(size(sample.Z, 2), 1)];

    est = struct();
    est.sample = sample;
    est.tRef = tRef;
    est.breakIndex = idxStar;
    est.breakDate = tc;
    est.breakDateString = decimal_year_to_string_local(tc);
    est.candidateBreakDates = sample.t(candIdx);
    est.sseProfile = sse;
    est.coeff = b;
    est.resid = u;
    est.Vhac = Vhac;

    est.position = summarize_delta_local(b, Vhac, g_xT, g_xbar, 1, 1, opts.inference.ciLevel);
    est.velocity = summarize_delta_local(b, Vhac, g_vT, g_vbar, 10, 10, opts.inference.ciLevel);
    est.slopeChange = summarize_jump_local(b, Vhac, g_jump, 10, opts.inference.ciLevel);
end

%% ========================================================================
function s = summarize_delta_local(b, V, gEndpoint, gTrailing, endpointScale, trailingScale, ciLevel)

    zCrit = zcrit_from_level_local(ciLevel);

    estE = endpointScale * (gEndpoint' * b);
    varE = endpointScale^2 * (gEndpoint' * V * gEndpoint);
    seE = sqrt(max(varE, 0));

    estT = trailingScale * (gTrailing' * b);
    varT = trailingScale^2 * (gTrailing' * V * gTrailing);
    seT = sqrt(max(varT, 0));

    s = struct();
    s.endpoint = estE;
    s.endpoint_ci_low = estE - zCrit * seE;
    s.endpoint_ci_high = estE + zCrit * seE;
    s.endpoint_se = seE;

    s.trailing60 = estT;
    s.trailing60_ci_low = estT - zCrit * seT;
    s.trailing60_ci_high = estT + zCrit * seT;
    s.trailing60_se = seT;
end

%% ========================================================================
function s = summarize_jump_local(b, V, g, scale, ciLevel)

    zCrit = zcrit_from_level_local(ciLevel);

    est = scale * (g' * b);
    varEst = scale^2 * (g' * V * g);
    se = sqrt(max(varEst, 0));

    s = struct();
    s.estimate = est;
    s.ci_low = est - zCrit * se;
    s.ci_high = est + zCrit * se;
    s.se = se;
end

%% ========================================================================
function z = zcrit_from_level_local(ciLevel)
    alpha = 1 - ciLevel;
    z = sqrt(2) * erfinv(2 * (1 - alpha / 2) - 1);
end

%% ========================================================================
function V = newey_west_local(X, u, L)

    [n, k] = size(X);
    XXinv = inv(X' * X);

    S = zeros(k, k);
    for t = 1:n
        xt = X(t, :)';
        S = S + (u(t)^2) * (xt * xt');
    end

    for ell = 1:L
        w = 1 - ell / (L + 1);
        G = zeros(k, k);
        for t = (ell + 1):n
            xt = X(t, :)';
            xl = X(t - ell, :)';
            G = G + u(t) * u(t - ell) * (xt * xl');
        end
        S = S + w * (G + G');
    end

    V = XXinv * S * XXinv;
    V = 0.5 * (V + V');
end

%% ========================================================================
function row = make_row_local(datasetName, s)
    row = struct();
    row.Dataset = string(datasetName);
    row.Endpoint_Est = s.endpoint;
    row.Endpoint_CI_Low = s.endpoint_ci_low;
    row.Endpoint_CI_High = s.endpoint_ci_high;
    row.Trailing60_Est = s.trailing60;
    row.Trailing60_CI_Low = s.trailing60_ci_low;
    row.Trailing60_CI_High = s.trailing60_ci_high;
end

%% ========================================================================
function row = make_jump_row_local(datasetName, s)
    row = struct();
    row.Dataset = string(datasetName);
    row.SlopeChange_Est = s.estimate;
    row.SlopeChange_CI_Low = s.ci_low;
    row.SlopeChange_CI_High = s.ci_high;
end

%% ========================================================================
function row = make_break_row_local(datasetName, est)
    row = struct();
    row.Dataset = string(datasetName);
    row.BreakDate_DecimalYear = est.breakDate;
    row.BreakDate_YearMonth = string(est.breakDateString);
    row.SSE_At_Break = est.sseProfile(est.breakIndex - opts_index_shift_local(est));
end

%% ========================================================================
function idx0 = opts_index_shift_local(est)
% Helper to map breakIndex back to the candidate-break vector index.
    idx0 = find(abs(est.candidateBreakDates - est.breakDate) < 1e-10, 1, 'first');
    if isempty(idx0)
        idx0 = NaN;
    end
end

%% ========================================================================
function s = decimal_year_to_string_local(t)

    yr = floor(t);
    frac = t - yr;

    % Approximate month from decimal year
    mo = floor(frac * 12) + 1;
    mo = max(1, min(12, mo));

    s = sprintf('%04d-%02d', yr, mo);
end

%% ========================================================================
function x = coerce_numeric_column_local(x)

    if isnumeric(x)
        x = x(:);
        return;
    end

    if islogical(x)
        x = double(x(:));
        return;
    end

    if iscell(x)
        x = string(x);
    end

    if isstring(x) || ischar(x) || iscategorical(x)
        x = str2double(string(x));
        x = x(:);
        return;
    end

    error('Could not coerce column to numeric vector.');
end

%% ========================================================================
function temp = rebase_temperature_series_local(temp, t)

% Foster-Rahmstorf preindustrial translation: each dataset is shifted by a
% constant so its 1991-2020 average equals 0.88 degC.

    targetMean = 0.88;
    maskRef = (t >= 1991) & (t < 2021);
    names = fieldnames(temp);

    for i = 1:numel(names)
        nm = names{i};
        mu = nanmean_local(temp.(nm)(maskRef));
        if ~isfinite(mu)
            error('Cannot rebase %s: no finite values over 1991-2020.', nm);
        end
        temp.(nm) = temp.(nm) + (targetMean - mu);
    end
end

%% ========================================================================
function m = nanmean_local(x)
    x = x(isfinite(x));
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end

%% ========================================================================
function s = nanstd_local(x)
    x = x(isfinite(x));
    if numel(x) < 2
        s = NaN;
    else
        s = std(x, 0);
    end
end

%% ========================================================================
function m = nanmean_rows_local(X)
    m = nan(size(X, 1), 1);
    for i = 1:size(X, 1)
        xi = X(i, :);
        xi = xi(isfinite(xi));
        if isempty(xi)
            m(i) = NaN;
        else
            m(i) = mean(xi);
        end
    end
end

%% ========================================================================
function out = merge_structs_local(a, b)
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
