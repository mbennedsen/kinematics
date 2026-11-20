
function out = estimate_quadratic_benchmark(dataFile, userOpts)
% ESTIMATE_QUADRATIC_BENCHMARK
% Parametric quadratic benchmark for the temperature-kinematics paper.
%
% Fits a quadratic time trend with climate covariates.
%   - fit y_t = alpha + gamma1 * ttilde + gamma2 * ttilde^2 + beta' z_t + u_t
%   - use the same Almon distributed-lag covariate structure as the SS model
%   - estimate endpoint and 5-year trailing-average x(t), v(t), a(t)
%   - use Newey-West HAC inference with Bartlett kernel and lag truncation L
%
% The script writes CSV files with the information needed for Tables 1-3:
%   Quadratic_Table1_Position.csv
%   Quadratic_Table2_Velocity.csv
%   Quadratic_Table3_Acceleration.csv

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

    posRows = cell(nD, 1);
    velRows = cell(nD, 1);
    accRows = cell(nD, 1);
    detail = struct();

    fprintf('\n============================================================\n');
    fprintf('Estimating quadratic benchmark model\n');
    fprintf('============================================================\n');
    fprintf('Window start year t0: %.0f\n', opts.quadratic.windowStartYear);
    fprintf('HAC lag truncation L: %d months\n', opts.quadratic.hacLagMonths);

    for i = 1:nD
        datasetName = datasets{i};
        fprintf('\n------------------------------------------------------------\n');
        fprintf('Dataset: %s\n', datasetName);
        fprintf('------------------------------------------------------------\n');

        sample = build_sample_local(data, datasetName, opts);
        est = fit_quadratic_model_local(sample, opts);

        fprintf('  x(T)      = %10.4f\n', est.position.endpoint);
        fprintf('  xbar(60)  = %10.4f\n', est.position.trailing60);
        fprintf('  v(T)      = %10.4f degC/decade\n', est.velocity.endpoint);
        fprintf('  vbar(60)  = %10.4f degC/decade\n', est.velocity.trailing60);
        fprintf('  a(T)      = %10.4f degC/decade^2\n', est.acceleration.endpoint);
        fprintf('  abar(60)  = %10.4f degC/decade^2\n', est.acceleration.trailing60);

        posRows{i} = make_row_local(datasetName, est.position);
        velRows{i} = make_row_local(datasetName, est.velocity);
        accRows{i} = make_row_local(datasetName, est.acceleration);

        detail.(datasetName) = est;
    end

    posTable = struct2table(vertcat(posRows{:}));
    velTable = struct2table(vertcat(velRows{:}));
    accTable = struct2table(vertcat(accRows{:}));

    writetable(posTable, fullfile(opts.outputDir, 'Quadratic_Table1_Position.csv'));
    writetable(velTable, fullfile(opts.outputDir, 'Quadratic_Table2_Velocity.csv'));
    writetable(accTable, fullfile(opts.outputDir, 'Quadratic_Table3_Acceleration.csv'));
    save(fullfile(opts.outputDir, 'Quadratic_Benchmark_Output.mat'), 'detail', 'opts', '-v7.3');

    out = struct();
    out.options = opts;
    out.positionTable = posTable;
    out.velocityTable = velTable;
    out.accelerationTable = accTable;
    out.detail = detail;

    fprintf('\nSaved quadratic benchmark output to:\n  %s\n', opts.outputDir);
end

%% ========================================================================
function opts = default_options_local()

    opts = struct();
    opts.outputDir = fullfile(pwd, 'results_quadratic_benchmark');
    opts.datasets = {'NASA','NOAA','HadCRU','Berk','ERA5'};
    opts.rebaseData = true;

    opts.covariates = struct();
    opts.covariates.names = {'NINO34','volc2','ssn'};
    opts.covariates.standardize = true;
    opts.covariates.maxLags = [12 12 12];
    opts.covariates.almonDegrees = [2 1 1];

    opts.quadratic = struct();
    opts.quadratic.windowStartYear = 2005;
    opts.quadratic.trailingMonths = 60;
    opts.quadratic.hacLagMonths = 24;

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

    valid = isfinite(y) & all(isfinite(Z), 2) & (data.t >= opts.quadratic.windowStartYear);
    idx = find(valid);

    if isempty(idx)
        error('No valid observations for dataset %s in quadratic window.', datasetName);
    end

    use = idx(1):idx(end);

    sample = struct();
    sample.name = datasetName;
    sample.t = data.t(use);
    sample.y = y(use);
    sample.Z = Z(use, :);
    sample.ZNames = zNames;
    sample.n = numel(sample.y);

    if sample.n < opts.quadratic.trailingMonths
        error('Quadratic sample too short for trailing %d months in dataset %s.', ...
            opts.quadratic.trailingMonths, datasetName);
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
function est = fit_quadratic_model_local(sample, opts)

    tbar = mean(sample.t);
    ttil = sample.t - tbar;
    Xtrend = [ones(sample.n, 1), ttil, ttil.^2];
    X = [Xtrend, sample.Z];

    b = X \ sample.y;
    u = sample.y - X * b;

    Vhac = newey_west_local(X, u, opts.quadratic.hacLagMonths);

    tT = ttil(end);
    tailIdx = (sample.n - opts.quadratic.trailingMonths + 1):sample.n;
    tTail = ttil(tailIdx);

    g_xT = [1; tT; tT.^2; zeros(size(sample.Z, 2), 1)];
    g_xbar = [1; mean(tTail); mean(tTail.^2); zeros(size(sample.Z, 2), 1)];

    g_vT = [0; 1; 2*tT; zeros(size(sample.Z, 2), 1)];
    g_vbar = [0; 1; 2*mean(tTail); zeros(size(sample.Z, 2), 1)];

    g_a = [0; 0; 2; zeros(size(sample.Z, 2), 1)];

    est = struct();
    est.sample = sample;
    est.coeff = b;
    est.resid = u;
    est.Vhac = Vhac;

    est.position = summarize_delta_local(b, Vhac, g_xT, g_xbar, 1, 1, opts.inference.ciLevel);
    est.velocity = summarize_delta_local(b, Vhac, g_vT, g_vbar, 10, 10, opts.inference.ciLevel);
    est.acceleration = summarize_delta_local(b, Vhac, g_a, g_a, 100, 100, opts.inference.ciLevel);
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
    s.endpoint_p = one_sided_p_local(estE, seE);

    s.trailing60 = estT;
    s.trailing60_ci_low = estT - zCrit * seT;
    s.trailing60_ci_high = estT + zCrit * seT;
    s.trailing60_p = one_sided_p_local(estT, seT);
end

%% ========================================================================
function z = zcrit_from_level_local(ciLevel)
    alpha = 1 - ciLevel;
    z = sqrt(2) * erfinv(2 * (1 - alpha / 2) - 1);
end

%% ========================================================================
function p = one_sided_p_local(est, se)

    if ~(isfinite(se) && se > 0)
        p = NaN;
        return;
    end
    z = est / se;
    p = 1 - normcdf_local(z);
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
    row.Endpoint_P = s.endpoint_p;
    row.Trailing60_Est = s.trailing60;
    row.Trailing60_CI_Low = s.trailing60_ci_low;
    row.Trailing60_CI_High = s.trailing60_ci_high;
    row.Trailing60_P = s.trailing60_p;
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
function p = normcdf_local(z)
    p = 0.5 * erfc(-z / sqrt(2));
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
