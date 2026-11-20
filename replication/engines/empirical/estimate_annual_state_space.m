
function out = estimate_annual_state_space(dataFile, userOpts)
% ESTIMATE_ANNUAL_STATE_SPACE
% Estimate the kinematic state-space model for annual observations.
%
% The annual model uses the following specification.
%   * annual dt = 1
%   * no monthly variance (useMonthlyVariance = false)
%   * free lag basis instead of Almon polynomial
%   * annual covariate lag structure
%   * trailing years instead of trailing months
%
% Output CSV files:
%   SS_Annual_Table1_Position.csv
%   SS_Annual_Table2_Velocity.csv
%   SS_Annual_Table3_Acceleration.csv

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

    rng(opts.seed, 'twister');
    if opts.bootstrap.useParallel
        p = gcp('nocreate');
        if isempty(p)
            if isempty(opts.bootstrap.nWorkers)
                parpool('local');
            else
                parpool('local', opts.bootstrap.nWorkers);
            end
        elseif ~isempty(opts.bootstrap.nWorkers) && p.NumWorkers ~= opts.bootstrap.nWorkers
            delete(p);
            parpool('local', opts.bootstrap.nWorkers);
        end
    end


    raw = readtable(dataFile);
    data = prepare_data_local(raw, opts);

    nD = numel(opts.datasets);
    posRows = cell(nD, 1);
    velRows = cell(nD, 1);
    accRows = cell(nD, 1);
    detail = struct();

    fprintf('\n============================================================\n');
    fprintf('Estimating annual state-space estimates\n');
    fprintf('============================================================\n');
    fprintf('Datasets: %s\n', strjoin(opts.datasets, ', '));
    fprintf('Bootstrap replications: %d\n', opts.bootstrap.B);
    fprintf('Include parameter uncertainty: %s\n', tf_to_yesno_local(opts.bootstrap.includeParameterUncertainty));
    fprintf('Bootstrap CI level: %.1f%%\n', 100 * (1 - opts.bootstrap.ciAlpha));
    fprintf('Use parfor: %s\n', tf_to_yesno_local(opts.bootstrap.useParallel));
    if opts.bootstrap.useParallel && ~isempty(opts.bootstrap.nWorkers)
        fprintf('Requested workers: %d\n', opts.bootstrap.nWorkers);
    end

    for j = 1:nD
        datasetName = opts.datasets{j};

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Dataset: %s\n', datasetName);
        fprintf('------------------------------------------------------------\n');

        sample = build_sample_local(data, datasetName, opts);

        fprintf('start year = %d\n', sample.t(1));
        fprintf('end year   = %d\n', sample.t(end));
        fprintf('T = %d, observed y = %d, internal missing y = %d\n', ...
            sample.T, sum(sample.yObserved), sum(~sample.yObserved));
        fprintf('Using maxLags      = [%s]\n', num2str(opts.covariates.maxLags));
        fprintf('Covariate basis structure:\n');
        for i = 1:numel(opts.covariates.names)
            fprintf('  %s: lags 0 to %d\n', opts.covariates.names{i}, opts.covariates.maxLags(i));
        end

        pointEst = estimate_individual_model_local(sample, opts.stateSpace, []);
        pointStates = pointEst.smoother.alphaSmooth(1:3, :);
        pointStatsRaw = make_estimands_local(pointStates, opts.trailingYears);
        pointStats = scale_estimands_local(pointStatsRaw);

        boot = run_bootstrap_local(sample, pointEst, opts);

        posRows{j} = make_table_row_local(datasetName, pointStats.x, boot.position);
        velRows{j} = make_table_row_local(datasetName, pointStats.v, boot.velocity);
        accRows{j} = make_table_row_local(datasetName, pointStats.a, boot.acceleration);

        detail.(datasetName) = struct();
        detail.(datasetName).sample = sample.meta;
        detail.(datasetName).pointEstimate = pointEst;
        detail.(datasetName).pointEstimandsRaw = pointStatsRaw;
        detail.(datasetName).pointEstimandsScaled = pointStats;
        detail.(datasetName).bootstrap = boot;

        fprintf('  x(T)      = %10.4f\n', pointStats.x.endpoint);
        fprintf('  xbar(%d)  = %10.4f\n', opts.trailingYears, pointStats.x.trailingN);
        fprintf('  v(T)      = %10.4f degC/decade\n', pointStats.v.endpoint);
        fprintf('  vbar(%d)  = %10.4f degC/decade\n', opts.trailingYears, pointStats.v.trailingN);
        fprintf('  a(T)      = %10.4f degC/decade^2\n', pointStats.a.endpoint);
        fprintf('  abar(%d)  = %10.4f degC/decade^2\n', opts.trailingYears, pointStats.a.trailingN);

        if opts.bootstrap.B > 0
            fprintf('  %.0f%% CI x(T):     [%10.4f, %10.4f]\n', 100 * (1 - opts.bootstrap.ciAlpha), boot.position.endpoint_ci_low, boot.position.endpoint_ci_high);
            fprintf('  %.0f%% CI v(T):     [%10.4f, %10.4f]\n', 100 * (1 - opts.bootstrap.ciAlpha), boot.velocity.endpoint_ci_low, boot.velocity.endpoint_ci_high);
            fprintf('  %.0f%% CI a(T):     [%10.4f, %10.4f]\n', 100 * (1 - opts.bootstrap.ciAlpha), boot.acceleration.endpoint_ci_low, boot.acceleration.endpoint_ci_high);
        end
    end

    posTable = struct2table(vertcat(posRows{:}));
    velTable = struct2table(vertcat(velRows{:}));
    accTable = struct2table(vertcat(accRows{:}));

    writetable(posTable, fullfile(opts.outputDir, 'SS_Annual_Table1_Position.csv'));
    writetable(velTable, fullfile(opts.outputDir, 'SS_Annual_Table2_Velocity.csv'));
    writetable(accTable, fullfile(opts.outputDir, 'SS_Annual_Table3_Acceleration.csv'));
    save(fullfile(opts.outputDir, 'SS_Annual_Table_Bootstrap_Output.mat'), 'detail', 'opts', '-v7.3');

    out = struct();
    out.options = opts;
    out.positionTable = posTable;
    out.velocityTable = velTable;
    out.accelerationTable = accTable;
    out.detail = detail;

    fprintf('\nSaved SS annual table outputs to:\n  %s\n', opts.outputDir);
end

%% ========================================================================
function opts = default_options_local()

    opts = struct();
    opts.seed = 12345;
    opts.outputDir = fullfile(pwd, 'results_annual_state_space');
    opts.datasets = {'NASA','NOAA','HadCRU','Berk','ERA5'};
    opts.rebaseData = true;
    opts.trailingYears = 5;

    opts.covariates = struct();
    opts.covariates.names = {'NINO34','volc2','ssn'};
    opts.covariates.standardize = true;
    opts.covariates.maxLags = [0 0 0];

    opts.stateSpace = struct();
    opts.stateSpace.dt = 1;
    opts.stateSpace.rhoBound = 0.995;
    opts.stateSpace.maBound = 0.995;
    opts.stateSpace.errorModel = 'ARMA11'; % 'ARMA11', 'AR1', or 'IID'
    opts.stateSpace.useMonthlyVariance = false;
    opts.stateSpace.varianceBreakYears = [];
    opts.stateSpace.jitter = 1e-8;
    opts.stateSpace.minVariance = 1e-10;
    opts.stateSpace.qLowerBound = 0;
    opts.stateSpace.fixedQ = []; % Optional fixed diffusion variance sigma^2.
    opts.stateSpace.optimizerMaxIter = 600;
    opts.stateSpace.optimizerTol = 1e-6;
    opts.stateSpace.verbose = false;

    opts.bootstrap = struct();
    opts.bootstrap.B = 499;
    opts.bootstrap.includeParameterUncertainty = true;
    opts.bootstrap.drawInitialState = true;
    opts.bootstrap.startAtOriginalTheta = true;
    opts.bootstrap.ciAlpha = 0.05;
    opts.bootstrap.maxRefitTries = 2;
    opts.bootstrap.useParallel = false;
    opts.bootstrap.nWorkers = [];

end

%% ========================================================================
function boot = run_bootstrap_local(sample, pointEst, opts)

    B = opts.bootstrap.B;
    if B <= 0
        boot = summarize_bootstrap_local([], [], [], opts.bootstrap.ciAlpha);
        return;
    end

    fprintf('  Bootstrapping (%d reps', B);
    if opts.bootstrap.includeParameterUncertainty
        fprintf(', re-estimating parameters');
    else
        fprintf(', fixed parameters via simulation smoother');
    end
    if opts.bootstrap.useParallel
        fprintf(', parfor');
    end
    fprintf(') ...\n');

    if ~opts.bootstrap.includeParameterUncertainty
        [repsX, repsV, repsA] = run_conditional_uncertainty_local(thetaHat_from_pointest_local(pointEst), sample, pointEst.setup, opts);
        boot = summarize_bootstrap_local(repsX, repsV, repsA, opts.bootstrap.ciAlpha);
        return;
    end

    repsX = nan(B, 2);
    repsV = nan(B, 2);
    repsA = nan(B, 2);

    if opts.bootstrap.useParallel
        parfor b = 1:B
            [xRow, vRow, aRow] = run_parameter_bootstrap_rep_local(sample, pointEst, opts);
            repsX(b, :) = xRow;
            repsV(b, :) = vRow;
            repsA(b, :) = aRow;

            fprintf('    %d / %d\n', b, B);
        end
    else
        for b = 1:B
            fprintf('    %d / %d\n', b, B);

            [xRow, vRow, aRow] = run_parameter_bootstrap_rep_local(sample, pointEst, opts);
            repsX(b, :) = xRow;
            repsV(b, :) = vRow;
            repsA(b, :) = aRow;
        end
    end

    repsX = repsX(all(isfinite(repsX),2), :);
    repsV = repsV(all(isfinite(repsV),2), :);
    repsA = repsA(all(isfinite(repsA),2), :);

    boot = summarize_bootstrap_local(repsX, repsV, repsA, opts.bootstrap.ciAlpha);
end


%% ========================================================================
function [repsX, repsV, repsA] = run_conditional_uncertainty_local(thetaHat, sample, setup, opts)

    B = opts.bootstrap.B;
    trailingYears = opts.trailingYears;

    [~, filt, pred, params] = negloglik_local(thetaHat, sample, setup, opts.stateSpace, false);
    F = state_transition_local(params, setup);
    sm = rts_smoother_local(filt, pred, F);

    repsX = nan(B, 2);
    repsV = nan(B, 2);
    repsA = nan(B, 2);

    if opts.bootstrap.useParallel
        parfor b = 1:B
            alphaDraw = simulation_smoother_draw_local(filt, sm, F);
            bsStats = scale_estimands_local(make_estimands_local(alphaDraw(1:3, :), trailingYears));
            repsX(b, :) = [bsStats.x.endpoint, bsStats.x.trailingN];
            repsV(b, :) = [bsStats.v.endpoint, bsStats.v.trailingN];
            repsA(b, :) = [bsStats.a.endpoint, bsStats.a.trailingN];
        end
    else
        for b = 1:B
            if mod(b, max(1, round(B/10))) == 0 || b == 1 || b == B
                fprintf('    %d / %d\n', b, B);
            end

            alphaDraw = simulation_smoother_draw_local(filt, sm, F);
            bsStats = scale_estimands_local(make_estimands_local(alphaDraw(1:3, :), trailingYears));

            repsX(b, :) = [bsStats.x.endpoint, bsStats.x.trailingN];
            repsV(b, :) = [bsStats.v.endpoint, bsStats.v.trailingN];
            repsA(b, :) = [bsStats.a.endpoint, bsStats.a.trailingN];
        end
    end
end


%% ========================================================================
function bsSample = simulate_bootstrap_sample_local(sample, pointEst, opts)

    params = pointEst.paramsHat;
    setup = pointEst.setup;
    F = state_transition_local(params, setup);
    H = setup.H(:)';
    T = sample.T;
    m = setup.stateDim;

    alpha = nan(m, T);
    y = nan(T, 1);

    a1 = pointEst.smoother.alphaSmooth(:, 1);
    P1 = regularize_cov_local(pointEst.smoother.PSmooth(:, :, 1));

    if opts.bootstrap.drawInitialState
        alpha(:, 1) = mvnrnd_local(a1, P1);
    else
        alpha(:, 1) = a1;
    end

    if sample.yObserved(1)
        y(1) = H * alpha(:, 1) + sample.W(1, :) * params.beta;
    end

    for t = 2:T
        shock = draw_state_noise_local(params, setup, t);
        alpha(:, t) = F * alpha(:, t-1) + shock;

        if sample.yObserved(t)
            y(t) = H * alpha(:, t) + sample.W(t, :) * params.beta;
        end
    end

    bsSample = sample;
    bsSample.y = y;
    bsSample.yObserved = ~isnan(y);
end

%% ========================================================================
function shock = draw_state_noise_local(params, setup, t)

    qxva = params.q * setup.SigmaXVA;
    xvaShock = mvnrnd_local(zeros(3,1), qxva);

    scale = 1;
    if numel(params.delta) >= setup.era(t)
        scale = scale * exp(params.delta(setup.era(t)));
    end
    Rt = max(params.sig2u * scale, setup.minVariance);
    uShock = sqrt(Rt) * randn;

    shock = [xvaShock; uShock; uShock];
end

%% ========================================================================
function out = evaluate_given_theta_local(sample, stateOpts, thetaHat)

    setup = build_setup_local(sample, stateOpts);
    [~, filt, pred, params] = negloglik_local(thetaHat, sample, setup, stateOpts, false);
    F = state_transition_local(params, setup);
    sm = rts_smoother_local(filt, pred, F);

    out = struct();
    out.thetaHat = thetaHat;
    out.logLik = NaN;
    out.exitflag = NaN;
    out.optimizerOutput = struct();
    out.filter = filt;
    out.pred = pred;
    out.smoother = sm;
    out.paramsHat = params;
    out.setup = setup;
end

%% ========================================================================
function alphaDraw = simulation_smoother_draw_local(filt, sm, F)

    T = size(filt.aFilt, 2);
    m = size(filt.aFilt, 1);

    alphaDraw = nan(m, T);
    alphaDraw(:, T) = mvnrnd_local(sm.alphaSmooth(:, T), sm.PSmooth(:, :, T));

    for t = T-1:-1:1
        Ptt = regularize_cov_local(filt.PFilt(:, :, t));
        Ptpred = regularize_cov_local(filt.PPred(:, :, t+1));
        Jt = (Ptt * F') / Ptpred;
        mt = filt.aFilt(:, t) + Jt * (alphaDraw(:, t+1) - filt.aPred(:, t+1));
        Ct = Ptt - Jt * Ptpred * Jt';
        Ct = regularize_cov_local(Ct);
        alphaDraw(:, t) = mvnrnd_local(mt, Ct);
    end
end

%% ========================================================================
function th = thetaHat_from_pointest_local(pointEst)
    th = pointEst.thetaHat;
end

%% ========================================================================
function [xRow, vRow, aRow] = run_parameter_bootstrap_rep_local(sample, pointEst, opts)

    xRow = [NaN, NaN];
    vRow = [NaN, NaN];
    aRow = [NaN, NaN];

    success = false;
    thetaStart = pointEst.thetaHat;

    for tryNum = 1:opts.bootstrap.maxRefitTries
        bsSample = simulate_bootstrap_sample_local(sample, pointEst, opts);

        try
            if opts.bootstrap.startAtOriginalTheta
                thetaStartUse = thetaStart;
            else
                thetaStartUse = [];
            end
            bsEst = estimate_individual_model_local(bsSample, opts.stateSpace, thetaStartUse);
            thetaStar = bsEst.thetaHat;
            success = all(isfinite(thetaStar));
            if success
                break;
            end
        catch
            success = false;
        end
    end

    if ~success
        return;
    end

    try
        condEst = evaluate_given_theta_local(sample, opts.stateSpace, thetaStar);
        F = state_transition_local(condEst.paramsHat, condEst.setup);
        alphaDraw = simulation_smoother_draw_local(condEst.filter, condEst.smoother, F);
        bsStats = scale_estimands_local(make_estimands_local(alphaDraw(1:3, :), opts.trailingYears));

        xRow = [bsStats.x.endpoint, bsStats.x.trailingN];
        vRow = [bsStats.v.endpoint, bsStats.v.trailingN];
        aRow = [bsStats.a.endpoint, bsStats.a.trailingN];
    catch
        % leave NaN rows
    end
end

%% ========================================================================
function summary = summarize_bootstrap_local(repsX, repsV, repsA, alpha)

    summary = struct();
    summary.position = summarize_two_stats_local(repsX, alpha);
    summary.velocity = summarize_two_stats_local(repsV, alpha);
    summary.acceleration = summarize_two_stats_local(repsA, alpha);

    summary.replicates = struct();
    summary.replicates.position = repsX;
    summary.replicates.velocity = repsV;
    summary.replicates.acceleration = repsA;
end

%% ========================================================================
function s = summarize_two_stats_local(X, alpha)

    if isempty(X)
        s = struct();
        s.endpoint_ci_low = NaN;
        s.endpoint_ci_high = NaN;
        s.endpoint_p = NaN;
        s.trailingN_ci_low = NaN;
        s.trailingN_ci_high = NaN;
        s.trailingN_p = NaN;
        return;
    end

    X = X(all(isfinite(X), 2), :);
    if isempty(X)
        s = struct();
        s.endpoint_ci_low = NaN;
        s.endpoint_ci_high = NaN;
        s.endpoint_p = NaN;
        s.trailingN_ci_low = NaN;
        s.trailingN_ci_high = NaN;
        s.trailingN_p = NaN;
        return;
    end

    lo = 100 * (alpha / 2);
    hi = 100 * (1 - alpha / 2);

    endpoint = X(:, 1);
    trailing = X(:, 2);

    s = struct();
    s.endpoint_ci_low = prctile(endpoint, lo);
    s.endpoint_ci_high = prctile(endpoint, hi);
    s.endpoint_p = empirical_one_sided_p_local(endpoint, 0);
    s.trailingN_ci_low = prctile(trailing, lo);
    s.trailingN_ci_high = prctile(trailing, hi);
    s.trailingN_p = empirical_one_sided_p_local(trailing, 0);
end

%% ========================================================================
function p = empirical_one_sided_p_local(x, nullValue)

    x = x(isfinite(x));
    if isempty(x)
        p = NaN;
        return;
    end
    p = (1 + sum(x <= nullValue)) / (numel(x) + 1);
end

%% ========================================================================
function row = make_table_row_local(datasetName, pointObj, bootObj)

    row = struct();
    row.Dataset = string(datasetName);
    row.Endpoint_Est = pointObj.endpoint;
    row.Endpoint_CI_Low = bootObj.endpoint_ci_low;
    row.Endpoint_CI_High = bootObj.endpoint_ci_high;
    row.Endpoint_P = bootObj.endpoint_p;
    row.TrailingN_Est = pointObj.trailingN;
    row.TrailingN_CI_Low = bootObj.trailingN_ci_low;
    row.TrailingN_CI_High = bootObj.trailingN_ci_high;
    row.TrailingN_P = bootObj.trailingN_p;
end

%% ========================================================================
function scaled = scale_estimands_local(raw)

    scaled = struct();

    scaled.x = struct();
    scaled.x.endpoint = raw.x_T;
    scaled.x.trailingN = raw.x_bar;

    scaled.v = struct();
    scaled.v.endpoint = 10 * raw.v_T;
    scaled.v.trailingN = 10 * raw.v_bar;

    scaled.a = struct();
    scaled.a.endpoint = 100 * raw.a_T;
    scaled.a.trailingN = 100 * raw.a_bar;
end

%% ========================================================================
function estimands = make_estimands_local(alphaPath, trailingYears)

    x = alphaPath(1, :);
    v = alphaPath(2, :);
    a = alphaPath(3, :);

    L = min(trailingYears, size(alphaPath, 2));
    idx = (size(alphaPath, 2) - L + 1):size(alphaPath, 2);

    estimands = struct();
    estimands.x_T = x(end);
    estimands.x_bar = mean(x(idx));
    estimands.v_T = v(end);
    estimands.v_bar = mean(v(idx));
    estimands.a_T = a(end);
    estimands.a_bar = mean(a(idx));
end

%% ========================================================================
function data = prepare_data_local(raw, opts)

    requiredCols = [{'year'}, opts.datasets, opts.covariates.names];
    for i = 1:numel(requiredCols)
        if ~ismember(requiredCols{i}, raw.Properties.VariableNames)
            error('Missing required column in CSV: %s', requiredCols{i});
        end
    end

    data = struct();
    data.t = coerce_numeric_local(raw.year);
    data.n = numel(data.t);
    data.calendarYear = data.t;

    for i = 1:numel(opts.datasets)
        nm = opts.datasets{i};
        data.temp.(nm) = coerce_numeric_local(raw.(nm));
    end

    if opts.rebaseData
        data.temp = rebase_temperature_local(data.temp, data.t);
    end

    covNames = opts.covariates.names;
    for i = 1:numel(covNames)
        nm = covNames{i};
        x = coerce_numeric_local(raw.(nm));
        if opts.covariates.standardize
            mu = nanmean_local(x);
            sd = nanstd_local(x);
            if ~(isfinite(sd) && sd > 0)
                error('Covariate %s has zero or undefined standard deviation.', nm);
            end
            data.covariates.(nm) = (x - mu) ./ sd;
        else
            data.covariates.(nm) = x;
        end
    end
end

%% ========================================================================
function temp = rebase_temperature_local(temp, t)
% NASA / NOAA / HadCRU / Berk: mean over 1850-1900 (available portion)
% ERA5: shifted to match mean-of-four over 1981-2010

    mainFour = {'NASA','NOAA','HadCRU','Berk'};
    mask1850 = (t >= 1850) & (t <= 1900);

    for i = 1:numel(mainFour)
        nm = mainFour{i};
        if isfield(temp, nm)
            mu = nanmean_local(temp.(nm)(mask1850));
            temp.(nm) = temp.(nm) - mu;
        end
    end

    if isfield(temp, 'ERA5')
        mask8110 = (t >= 1981) & (t <= 2010);
        X = nan(sum(mask8110), numel(mainFour));
        for i = 1:numel(mainFour)
            nm = mainFour{i};
            if isfield(temp, nm)
                X(:, i) = temp.(nm)(mask8110);
            end
        end
        mu4 = nanmean_local(nanmean_rows_local(X));
        muE = nanmean_local(temp.ERA5(mask8110));
        temp.ERA5 = temp.ERA5 + (mu4 - muE);
    end
end

%% ========================================================================
function sample = build_sample_local(data, datasetName, opts)

    lagBasis = build_lag_basis_local(data.covariates, opts.covariates.names, ...
        opts.covariates.maxLags);

    yFull = data.temp.(datasetName);
    validW = all(isfinite(lagBasis.W), 2);
    obsY = ~isnan(yFull);

    validObs = validW & obsY;
    idx = find(validObs);
    if isempty(idx)
        error('No valid observations for %s after building lag basis.', datasetName);
    end

    startIdx = idx(1);
    endIdx = idx(end);

    use = false(size(yFull));
    use(startIdx:endIdx) = validW(startIdx:endIdx);

    sample = struct();
    sample.name = datasetName;
    sample.y = yFull(use);
    sample.W = lagBasis.W(use, :);
    sample.WNames = lagBasis.WNames;
    sample.t = data.t(use);
    sample.era = make_era_index_local(data.calendarYear(use), opts.stateSpace.varianceBreakYears);
    sample.yObserved = ~isnan(sample.y);
    sample.T = numel(sample.y);
    sample.lagInfo = lagBasis;

    sample.meta = struct();
    sample.meta.dataset = datasetName;
    sample.meta.startYear = sample.t(1);
    sample.meta.endYear = sample.t(end);
    sample.meta.T = sample.T;
    sample.meta.observedY = sum(sample.yObserved);
    sample.meta.internalMissingY = sum(~sample.yObserved);
    sample.meta.WNames = sample.WNames;
end

%% ========================================================================
function lagBasis = build_lag_basis_local(covStruct, covNames, maxLags)
% Build free lag basis: for each covariate with maxLag L, create L+1 columns
% (lag 0, lag 1, ..., lag L). No polynomial basis.

    nCov = numel(covNames);
    n = numel(covStruct.(covNames{1}));
    K = sum(maxLags + 1);

    W = nan(n, K);
    WNames = cell(1, K);
    coeffInfo = cell(1, nCov);

    k = 0;
    for i = 1:nCov
        x = covStruct.(covNames{i});
        L = maxLags(i);

        coeffInfo{i} = struct();
        coeffInfo{i}.covariate = covNames{i};
        coeffInfo{i}.maxLag = L;
        coeffInfo{i}.lagIndices = 0:L;

        for lag = 0:L
            col = nan(n, 1);
            for t = (L+1):n
                val = x(t - lag);
                if isfinite(val)
                    col(t) = val;
                end
            end
            k = k + 1;
            W(:, k) = col;
            WNames{k} = sprintf('%s_lag%d', covNames{i}, lag);
        end
    end

    lagBasis = struct();
    lagBasis.W = W;
    lagBasis.WNames = WNames;
    lagBasis.coeffInfo = coeffInfo;
end

%% ========================================================================
function era = make_era_index_local(calendarYear, breakYears)

    if isempty(breakYears)
        era = ones(size(calendarYear));
        return;
    end

    breaks = unique(round(breakYears(:)'));
    era = ones(size(calendarYear));
    for i = 1:numel(breaks)
        era(calendarYear >= breaks(i)) = i + 1;
    end
end

%% ========================================================================
function est = estimate_individual_model_local(sample, opts, thetaStart)

    setup = build_setup_local(sample, opts);
    if nargin < 3 || isempty(thetaStart)
        thetaStart = initial_theta_local(sample, setup);
    end

    [thetaHat, logLik, exitflag, out] = optimize_local(thetaStart, sample, setup, opts);
    [~, filt, pred, params] = negloglik_local(thetaHat, sample, setup, opts, false);
    F = state_transition_local(params, setup);
    sm = rts_smoother_local(filt, pred, F);

    est = struct();
    est.thetaHat = thetaHat;
    est.logLik = logLik;
    est.exitflag = exitflag;
    est.optimizerOutput = out;
    est.filter = filt;
    est.smoother = sm;
    est.paramsHat = params;
    est.setup = setup;
end

%% ========================================================================
function setup = build_setup_local(sample, opts)

    setup = struct();
    setup.dt = opts.dt;
    setup.T = sample.T;
    setup.K = size(sample.W, 2);
    setup.era = sample.era(:);
    setup.nEra = max(sample.era(:));
    setup.rhoBound = opts.rhoBound;
    setup.maBound = opts.maBound;
    setup.errorModel = upper(string(opts.errorModel));
    if ~(setup.errorModel == "ARMA11" || setup.errorModel == "AR1" || setup.errorModel == "IID")
        error('Unsupported stateSpace.errorModel: %s. Use ARMA11, AR1, or IID.', opts.errorModel);
    end
    setup.estimateRho = setup.errorModel ~= "IID";
    setup.estimateMA = setup.errorModel == "ARMA11";
    setup.jitter = opts.jitter;
    setup.minVariance = opts.minVariance;
    setup.qLowerBound = opts.qLowerBound;
    setup.fixedQ = opts.fixedQ;
    if ~isempty(setup.fixedQ)
        validateattributes(setup.fixedQ, {'numeric'}, {'scalar','real','finite','nonnegative'});
    end

    dt = setup.dt;
    setup.SigmaXVA = [dt^5/20, dt^4/8, dt^3/6; ...
                      dt^4/8,  dt^3/3, dt^2/2; ...
                      dt^3/6,  dt^2/2, dt      ];

    setup.numGammaFree = 0;
    setup.numDeltaFree = max(setup.nEra - 1, 0);

    setup.stateDim = 5;
    setup.H = [1, 0, 0, 1, 0];
end

%% ========================================================================
function theta0 = initial_theta_local(sample, setup)

    obs = sample.yObserved;
    y = sample.y(obs);
    W = sample.W(obs, :);

    beta0 = W \ y;
    u = y - W * beta0;

    rho0 = 0.5;
    if numel(u) > 2
        c = corrcoef(u(2:end), u(1:end-1));
        if numel(c) >= 4 && all(isfinite(c(:)))
            rho0 = max(min(c(1,2), 0.9), -0.9);
        end
    end
    q0 = max(1e-8, 1e-4 * nanvar_local(u));
    sig2u0 = max(1e-6, 0.25 * nanvar_local(u));

    delta0 = zeros(setup.numDeltaFree, 1);

    theta0 = beta0(:);
    if setup.estimateRho
        theta0 = [theta0; atanh(rho0 / setup.rhoBound)];
    end
    if setup.estimateMA
        ma0 = 0.1;
        theta0 = [theta0; atanh(ma0 / setup.maBound)];
    end
    if isempty(setup.fixedQ)
        theta0 = [theta0; log(max(q0, 1e-10))];
    end
    theta0 = [theta0; log(max(sig2u0, 1e-10)); delta0];
end

%% ========================================================================
function [thetaHat, logLik, exitflag, out] = optimize_local(theta0, sample, setup, opts)

    obj = @(th) negloglik_local(th, sample, setup, opts, true);

    hasFminunc = exist('fminunc', 'file') == 2;
    if hasFminunc
        optimOpts = optimoptions('fminunc', ...
            'Algorithm', 'quasi-newton', ...
            'Display', 'off', ...
            'MaxIterations', opts.optimizerMaxIter, ...
            'MaxFunctionEvaluations', 50000, ...
            'OptimalityTolerance', opts.optimizerTol, ...
            'StepTolerance', 1e-10);
        [thetaHat, fval, exitflag, out] = fminunc(obj, theta0, optimOpts);
    else
        optimOpts = optimset('Display', 'off', ...
            'MaxIter', opts.optimizerMaxIter, ...
            'MaxFunEvals', 50000, ...
            'TolX', opts.optimizerTol, ...
            'TolFun', opts.optimizerTol);
        [thetaHat, fval, exitflag, out] = fminsearch(obj, theta0, optimOpts);
    end

    logLik = -fval;
end

%% ========================================================================
function [nll, filt, pred, params] = negloglik_local(thetaRaw, sample, setup, opts, onlyNll)

    if nargin < 5
        onlyNll = false;
    end

    params = unpack_params_local(thetaRaw, setup);
    F = state_transition_local(params, setup);
    H = setup.H(:)';

    m = setup.stateDim;
    T = sample.T;
    y = sample.y;
    W = sample.W;

    firstObs = find(sample.yObserved, 1, 'first');
    if isempty(firstObs)
        error('Sample contains no observed y values.');
    end

    aPrev = zeros(m, 1);
    aPrev(1) = y(firstObs) - W(firstObs, :) * params.beta;

    PPrev = 100 * eye(m);
    PPrev(4,4) = 10;
    PPrev(5,5) = 10;

    nll = 0;

    if ~onlyNll
        filt = struct();
        filt.aPred = nan(m, T);
        filt.PPred = nan(m, m, T);
        filt.aFilt = nan(m, T);
        filt.PFilt = nan(m, m, T);
        filt.v = nan(T, 1);
        filt.F = nan(T, 1);
        filt.K = nan(m, T);
        filt.vStd = nan(T, 1);
        pred = struct();
        pred.a0 = aPrev;
        pred.P0 = PPrev;
    else
        filt = [];
        pred = [];
    end

    for t = 1:T
        if t == 1
            aPred = aPrev;
            PPred = PPrev;
        else
            Qt = state_noise_local(params, setup, t);
            aPred = F * aPrev;
            PPred = F * PPrev * F' + Qt;
        end
        PPred = regularize_cov_local(PPred);

        muObs = H * aPred + sample.W(t, :) * params.beta;

        if sample.yObserved(t)
            vt = y(t) - muObs;
            Ft = H * PPred * H' + setup.jitter;
            Ft = max(Ft, setup.jitter);

            Kt = (PPred * H') / Ft;
            aFilt = aPred + Kt * vt;
            PFilt = PPred - Kt * Ft * Kt';
            PFilt = regularize_cov_local(PFilt);

            nll = nll + 0.5 * (log(2*pi) + log(Ft) + vt^2 / Ft);
        else
            vt = NaN;
            Ft = NaN;
            Kt = zeros(m, 1);
            aFilt = aPred;
            PFilt = PPred;
        end

        aPrev = aFilt;
        PPrev = PFilt;

        if ~onlyNll
            filt.aPred(:, t) = aPred;
            filt.PPred(:, :, t) = PPred;
            filt.aFilt(:, t) = aFilt;
            filt.PFilt(:, :, t) = PFilt;
            filt.v(t) = vt;
            filt.F(t) = Ft;
            filt.K(:, t) = Kt;
            if isfinite(vt) && isfinite(Ft) && Ft > 0
                filt.vStd(t) = vt / sqrt(Ft);
            end
        end
    end

    if ~isfinite(nll)
        nll = 1e12;
    end
end

%% ========================================================================
function params = unpack_params_local(thetaRaw, setup)

    k = 0;
    K = setup.K;
    params = struct();

    params.beta = thetaRaw(1:K);
    k = K;

    if setup.estimateRho
        params.rho = setup.rhoBound * tanh(thetaRaw(k+1));
        k = k + 1;
    else
        params.rho = 0;
    end

    if setup.estimateMA
        params.ma = setup.maBound * tanh(thetaRaw(k+1));
        k = k + 1;
    else
        params.ma = 0;
    end

    if ~isempty(setup.fixedQ)
        params.q = setup.fixedQ;
    elseif setup.qLowerBound > 0
        params.q = setup.qLowerBound + exp(thetaRaw(k+1));
        k = k + 1;
    else
        params.q = exp(thetaRaw(k+1));
        k = k + 1;
    end

    params.sig2u = exp(thetaRaw(k+1));
    k = k + 1;

    nDelta = setup.numDeltaFree;
    if nDelta > 0
        deltaFree = thetaRaw(k+1:k+nDelta);
        params.delta = [0; deltaFree(:)];
    else
        params.delta = 0;
    end
end

%% ========================================================================
function F = state_transition_local(params, setup)

    dt = setup.dt;
    F = [1, dt, 0.5*dt^2, 0,         0; ...
         0, 1,  dt,       0,         0; ...
         0, 0,  1,        0,         0; ...
         0, 0,  0, params.rho, params.ma; ...
         0, 0,  0, 0,         0];
end

%% ========================================================================
function Qt = state_noise_local(params, setup, t)

    Qxva = params.q * setup.SigmaXVA;

    scale = 1;
    if numel(params.delta) >= setup.era(t)
        scale = scale * exp(params.delta(setup.era(t)));
    end

    Rt = max(params.sig2u * scale, setup.minVariance);
    Qu = Rt * [1 1; 1 1];

    Qt = blkdiag(Qxva, Qu);
    Qt = symmetrize_local(Qt);
end

%% ========================================================================
function sm = rts_smoother_local(filt, pred, F)

    %#ok<INUSD> pred retained for compatibility

    T = size(filt.aFilt, 2);
    m = size(filt.aFilt, 1);

    alphaSmooth = nan(m, T);
    PSmooth = nan(m, m, T);

    alphaSmooth(:, T) = filt.aFilt(:, T);
    PSmooth(:, :, T) = filt.PFilt(:, :, T);

    for t = T-1:-1:1
        Ptt = regularize_cov_local(filt.PFilt(:, :, t));
        Ptpred = regularize_cov_local(filt.PPred(:, :, t+1));
        Jt = (Ptt * F') / Ptpred;
        alphaSmooth(:, t) = filt.aFilt(:, t) + Jt * (alphaSmooth(:, t+1) - filt.aPred(:, t+1));
        PSmooth(:, :, t) = Ptt + Jt * (PSmooth(:, :, t+1) - Ptpred) * Jt';
        PSmooth(:, :, t) = regularize_cov_local(PSmooth(:, :, t));
    end

    sm = struct();
    sm.alphaSmooth = alphaSmooth;
    sm.PSmooth = PSmooth;
end

%% ========================================================================
function S = regularize_cov_local(S)
    S = symmetrize_local(S);
    [V, D] = eig(S);
    d = diag(D);
    d(~isfinite(d)) = 0;
    d = max(d, 1e-10);
    S = V * diag(d) * V';
    S = symmetrize_local(S);
end

%% ========================================================================
function A = symmetrize_local(A)
    A = 0.5 * (A + A');
end

%% ========================================================================
function x = coerce_numeric_local(x)

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
function v = nanvar_local(x)
    x = x(isfinite(x));
    if numel(x) < 2
        v = NaN;
    else
        v = var(x, 0);
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
function x = mvnrnd_local(mu, Sigma)

    mu = mu(:);
    Sigma = regularize_cov_local(Sigma);
    C = chol(Sigma, 'lower');
    x = mu + C * randn(numel(mu), 1);
end

%% ========================================================================
function s = tf_to_yesno_local(tf)
    if tf
        s = 'yes';
    else
        s = 'no';
    end
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
