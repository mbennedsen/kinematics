function out = estimate_velocity_rw_robustness(dataFile, userOpts)
% ESTIMATE_VELOCITY_RW_ROBUSTNESS
% Robustness check for the kinematic state-space hierarchy without
% acceleration:
%
%   dx(t) = v(t) dt
%   dv(t) = sigma dW(t)
%
% The observation equation and nuisance structure mirror the individual
% x-v-a state-space model: Almon-lag covariates plus an ARMA(1,1)
% irregular component. The main outputs are endpoint/trailing estimates of
% x(t) and v(t), together with x-v plots for each dataset.

    assert(nargin >= 1 && ~isempty(dataFile),'Provide the path to the analysis panel.');
    if nargin < 2
        userOpts = struct();
    end

    opts = default_options_local();
    opts = merge_structs_local(opts, userOpts);

    if ~exist(dataFile, 'file') && exist(fullfile('Data', dataFile), 'file')
        dataFile = fullfile('Data', dataFile);
    end
    if ~exist(dataFile, 'file')
        error('Data file not found: %s', dataFile);
    end
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end
    figDir = fullfile(opts.outputDir, 'figures');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end

    rng(opts.seed, 'twister');
    raw = readtable(dataFile);
    data = prepare_data_local(raw, opts);

    nD = numel(opts.datasets);
    posRows = cell(nD, 1);
    velRows = cell(nD, 1);
    detail = struct();

    fprintf('\n============================================================\n');
    fprintf('Velocity-RW robustness state-space model\n');
    fprintf('============================================================\n');
    fprintf('Datasets: %s\n', strjoin(opts.datasets, ', '));
    fprintf('Model: dx = v dt, dv = sigma dW\n');

    for j = 1:nD
        datasetName = opts.datasets{j};

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Dataset: %s\n', datasetName);
        fprintf('------------------------------------------------------------\n');

        sample = build_sample_local(data, datasetName, opts);
        ssOpts = state_space_options_for_dataset_local(opts.stateSpace, datasetName);
        est = estimate_model_local(sample, ssOpts, []);
        states = est.smoother.alphaSmooth(1:2, :);
        statsRaw = make_estimands_local(states, opts.trailingMonths);
        stats = scale_estimands_local(statsRaw);

        posRows{j} = make_table_row_local(datasetName, stats.x);
        velRows{j} = make_table_row_local(datasetName, stats.v);

        detail.(datasetName) = struct();
        detail.(datasetName).sample = sample.meta;
        detail.(datasetName).estimate = est;
        detail.(datasetName).estimandsRaw = statsRaw;
        detail.(datasetName).estimandsScaled = stats;

        figFile = save_xv_plot_local(sample, est, stats, opts, figDir);
        detail.(datasetName).figure = figFile;

        fprintf('  x(T)      = %10.4f\n', stats.x.endpoint);
        fprintf('  xbar(60)  = %10.4f\n', stats.x.trailing60);
        fprintf('  v(T)      = %10.4f degC/decade\n', stats.v.endpoint);
        fprintf('  vbar(60)  = %10.4f degC/decade\n', stats.v.trailing60);
        fprintf('  figure    = %s\n', figFile);
    end

    posTable = struct2table(vertcat(posRows{:}));
    velTable = struct2table(vertcat(velRows{:}));

    writetable(posTable, fullfile(opts.outputDir, 'XVRobustness_Table1_Position.csv'));
    writetable(velTable, fullfile(opts.outputDir, 'XVRobustness_Table2_Velocity.csv'));
    save(fullfile(opts.outputDir, 'XVRobustness_Output.mat'), 'detail', 'opts', '-v7.3');

    out = struct();
    out.options = opts;
    out.positionTable = posTable;
    out.velocityTable = velTable;
    out.detail = detail;

    fprintf('\nSaved velocity-RW robustness outputs to:\n  %s\n', opts.outputDir);
end

%% ========================================================================
function opts = default_options_local()
    opts = struct();
    opts.seed = 24680;
    opts.outputDir = fullfile(pwd, 'results_velocity_rw_robustness');
    opts.datasets = {'NASA','NOAA','HadCRU','Berk','ERA5'};
    opts.rebaseData = true;
    opts.trailingMonths = 60;

    opts.covariates = struct();
    opts.covariates.names = {'NINO34','volc2','ssn'};
    opts.covariates.standardize = true;
    opts.covariates.maxLags = [12 12 12];
    opts.covariates.almonDegrees = [2 1 1];

    opts.stateSpace = struct();
    opts.stateSpace.dt = 1/12;
    opts.stateSpace.rhoBound = 0.995;
    opts.stateSpace.maBound = 0.995;
    opts.stateSpace.useMonthlyVariance = true;
    opts.stateSpace.varianceBreakYears = [];
    opts.stateSpace.jitter = 1e-8;
    opts.stateSpace.minVariance = 1e-10;
    opts.stateSpace.qLowerBound = 0;
    opts.stateSpace.fixedQ = [];
    opts.stateSpace.fixedQByDataset = struct();
    opts.stateSpace.optimizerMaxIter = 600;
    opts.stateSpace.optimizerTol = 1e-6;

    opts.figure = struct();
    opts.figure.visible = 'off';
    opts.figure.format = 'png';
    opts.figure.width = 950;
    opts.figure.height = 680;
    opts.figure.plotStart = -Inf;
end

%% ========================================================================
function ssOpts = state_space_options_for_dataset_local(baseOpts, datasetName)
    ssOpts = baseOpts;
    if isfield(baseOpts, 'fixedQByDataset') && isstruct(baseOpts.fixedQByDataset) && ...
            isfield(baseOpts.fixedQByDataset, datasetName)
        ssOpts.fixedQ = baseOpts.fixedQByDataset.(datasetName);
    end
end

%% ========================================================================
function data = prepare_data_local(raw, opts)
    requiredCols = [{'t'}, opts.datasets, opts.covariates.names];
    for i = 1:numel(requiredCols)
        if ~ismember(requiredCols{i}, raw.Properties.VariableNames)
            error('Missing required column in CSV: %s', requiredCols{i});
        end
    end

    data = struct();
    data.t = coerce_numeric_local(raw.t);
    data.n = numel(data.t);
    data.calendarYear = floor(data.t + 1e-9);
    data.month = mod((0:data.n-1)', 12) + 1;

    for i = 1:numel(opts.datasets)
        nm = opts.datasets{i};
        data.temp.(nm) = coerce_numeric_local(raw.(nm));
    end

    if opts.rebaseData
        data.temp = rebase_temperature_local(data.temp, data.t);
    end

    for i = 1:numel(opts.covariates.names)
        nm = opts.covariates.names{i};
        x = coerce_numeric_local(raw.(nm));
        if opts.covariates.standardize
            mu = nanmean_local(x);
            sd = nanstd_local(x);
            if ~(isfinite(sd) && sd > 0)
                error('Covariate %s has zero or undefined standard deviation.', nm);
            end
            x = (x - mu) ./ sd;
        end
        data.covariates.(nm) = x;
    end
end

function temp = rebase_temperature_local(temp, t)
    targetMean = 0.88;
    ref = t >= 1991 & t < 2021;
    names = fieldnames(temp);
    for i = 1:numel(names)
        nm = names{i};
        mu = nanmean_local(temp.(nm)(ref));
        if ~isfinite(mu)
            error('Cannot rebase %s: no finite values over 1991--2020.', nm);
        end
        temp.(nm) = temp.(nm) + (targetMean - mu);
    end
end

%% ========================================================================
function sample = build_sample_local(data, datasetName, opts)
    almon = build_almon_basis_local(data.covariates, opts.covariates.names, ...
        opts.covariates.maxLags, opts.covariates.almonDegrees);

    yFull = data.temp.(datasetName);
    validW = all(isfinite(almon.W), 2);
    obsY = ~isnan(yFull);
    idx = find(validW & obsY);
    if isempty(idx)
        error('No valid observations for %s after building lag basis.', datasetName);
    end

    use = false(size(yFull));
    use(idx(1):idx(end)) = validW(idx(1):idx(end));

    sample = struct();
    sample.name = datasetName;
    sample.y = yFull(use);
    sample.W = almon.W(use, :);
    sample.WNames = almon.WNames;
    sample.t = data.t(use);
    sample.month = data.month(use);
    sample.era = make_era_index_local(data.calendarYear(use), opts.stateSpace.varianceBreakYears);
    sample.yObserved = ~isnan(sample.y);
    sample.T = numel(sample.y);

    sample.meta = struct();
    sample.meta.dataset = datasetName;
    sample.meta.startDate = sample.t(1);
    sample.meta.endDate = sample.t(end);
    sample.meta.T = sample.T;
    sample.meta.observedY = sum(sample.yObserved);
    sample.meta.internalMissingY = sum(~sample.yObserved);
    sample.meta.WNames = sample.WNames;
end

function almon = build_almon_basis_local(covStruct, covNames, maxLags, degrees)
    nCov = numel(covNames);
    n = numel(covStruct.(covNames{1}));
    K = sum(degrees + 1);
    W = nan(n, K);
    WNames = cell(1, K);
    k = 0;

    for i = 1:nCov
        x = covStruct.(covNames{i});
        L = maxLags(i);
        deg = degrees(i);
        lagGrid = (0:L)';
        lagScaled = lagGrid;
        if max(lagGrid) > 0
            lagScaled = lagGrid / max(lagGrid);
        end

        for p = 0:deg
            col = nan(n, 1);
            for t = (L+1):n
                vals = x(t-lagGrid);
                if all(isfinite(vals))
                    col(t) = sum((lagScaled.^p) .* vals);
                end
            end
            k = k + 1;
            W(:, k) = col;
            WNames{k} = sprintf('%s_polydeg%d', covNames{i}, p);
        end
    end

    almon = struct('W', W, 'WNames', {WNames});
end

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
function est = estimate_model_local(sample, opts, thetaStart)
    setup = build_setup_local(sample, opts);
    if nargin < 3 || isempty(thetaStart)
        thetaStart = initial_theta_local(sample, setup);
    end

    [thetaHat, logLik, exitflag, optOut] = optimize_local(thetaStart, sample, setup, opts);
    [~, filt, pred, params] = negloglik_local(thetaHat, sample, setup, opts, false);
    F = state_transition_local(params, setup);
    sm = rts_smoother_local(filt, pred, F);

    est = struct();
    est.thetaHat = thetaHat;
    est.logLik = logLik;
    est.exitflag = exitflag;
    est.optimizerOutput = optOut;
    est.filter = filt;
    est.pred = pred;
    est.smoother = sm;
    est.paramsHat = params;
    est.setup = setup;
end

function setup = build_setup_local(sample, opts)
    setup = struct();
    setup.dt = opts.dt;
    setup.T = sample.T;
    setup.K = size(sample.W, 2);
    setup.month = sample.month(:);
    setup.era = sample.era(:);
    setup.nEra = max(sample.era(:));
    setup.rhoBound = opts.rhoBound;
    setup.maBound = opts.maBound;
    setup.useMonthlyVariance = opts.useMonthlyVariance;
    setup.jitter = opts.jitter;
    setup.minVariance = opts.minVariance;
    setup.qLowerBound = opts.qLowerBound;
    if isfield(opts, 'fixedQ')
        setup.fixedQ = opts.fixedQ;
    else
        setup.fixedQ = [];
    end

    dt = setup.dt;
    setup.SigmaXV = [dt^3/3, dt^2/2; ...
                     dt^2/2, dt];

    setup.numGammaFree = 0;
    if setup.useMonthlyVariance
        setup.numGammaFree = 11;
    end
    setup.numDeltaFree = max(setup.nEra - 1, 0);
    setup.stateDim = 4;
    setup.H = [1, 0, 1, 0];
end

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

    thetaPieces = {beta0(:); ...
        atanh(rho0 / setup.rhoBound); ...
        atanh(0.1 / setup.maBound)};

    if isempty(setup.fixedQ)
        q0 = max(1e-8, 1e-4 * nanvar_local(u));
        thetaPieces{end+1} = log(max(q0, 1e-10)); %#ok<AGROW>
    end

    sig2u0 = max(1e-6, 0.25 * nanvar_local(u));
    thetaPieces = [thetaPieces; ...
        {log(max(sig2u0, 1e-10)); ...
         zeros(setup.numGammaFree, 1); ...
         zeros(setup.numDeltaFree, 1)}];
    theta0 = vertcat(thetaPieces{:});
end

function [thetaHat, logLik, exitflag, out] = optimize_local(theta0, sample, setup, opts)
    obj = @(th) negloglik_local(th, sample, setup, opts, true);

    if exist('fminunc', 'file') == 2
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
    firstObs = find(sample.yObserved, 1, 'first');
    if isempty(firstObs)
        error('Sample contains no observed y values.');
    end

    aPrev = zeros(m, 1);
    aPrev(1) = y(firstObs) - sample.W(firstObs, :) * params.beta;
    PPrev = 100 * eye(m);
    PPrev(3,3) = 10;
    PPrev(4,4) = 10;
    nll = 0;

    if ~onlyNll
        filt = struct();
        filt.aPred = nan(m, T);
        filt.PPred = nan(m, m, T);
        filt.aFilt = nan(m, T);
        filt.PFilt = nan(m, m, T);
        filt.v = nan(T, 1);
        filt.F = nan(T, 1);
        filt.vStd = nan(T, 1);
        pred = struct('a0', aPrev, 'P0', PPrev);
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
            if isfinite(vt) && isfinite(Ft) && Ft > 0
                filt.vStd(t) = vt / sqrt(Ft);
            end
        end
    end

    if ~isfinite(nll)
        nll = 1e12;
    end
end

function params = unpack_params_local(thetaRaw, setup)
    k = 0;
    K = setup.K;
    params = struct();
    params.beta = thetaRaw(1:K);
    k = K;

    params.rho = setup.rhoBound * tanh(thetaRaw(k+1));
    k = k + 1;
    params.ma = setup.maBound * tanh(thetaRaw(k+1));
    k = k + 1;

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

    if setup.numGammaFree > 0
        gammaFree = thetaRaw(k+1:k+setup.numGammaFree);
        k = k + setup.numGammaFree;
        params.gamma = [gammaFree(:); -sum(gammaFree)];
    else
        params.gamma = zeros(12, 1);
    end

    if setup.numDeltaFree > 0
        deltaFree = thetaRaw(k+1:k+setup.numDeltaFree);
        params.delta = [0; deltaFree(:)];
    else
        params.delta = 0;
    end
end

function F = state_transition_local(params, setup)
    dt = setup.dt;
    F = [1, dt, 0,          0; ...
         0, 1,  0,          0; ...
         0, 0,  params.rho, params.ma; ...
         0, 0,  0,          0];
end

function Qt = state_noise_local(params, setup, t)
    Qxv = params.q * setup.SigmaXV;
    scale = 1;
    if setup.useMonthlyVariance
        scale = scale * exp(params.gamma(setup.month(t)));
    end
    if numel(params.delta) >= setup.era(t)
        scale = scale * exp(params.delta(setup.era(t)));
    end
    Rt = max(params.sig2u * scale, setup.minVariance);
    Qu = Rt * [1 1; 1 1];
    Qt = blkdiag(Qxv, Qu);
    Qt = symmetrize_local(Qt);
end

function sm = rts_smoother_local(filt, pred, F)
    %#ok<INUSD>
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

    sm = struct('alphaSmooth', alphaSmooth, 'PSmooth', PSmooth);
end

%% ========================================================================
function estimands = make_estimands_local(alphaPath, trailingMonths)
    x = alphaPath(1, :);
    v = alphaPath(2, :);
    L = min(trailingMonths, size(alphaPath, 2));
    idx = (size(alphaPath, 2) - L + 1):size(alphaPath, 2);
    estimands = struct();
    estimands.x_T = x(end);
    estimands.x_bar = mean(x(idx));
    estimands.v_T = v(end);
    estimands.v_bar = mean(v(idx));
end

function scaled = scale_estimands_local(raw)
    scaled = struct();
    scaled.x = struct('endpoint', raw.x_T, 'trailing60', raw.x_bar);
    scaled.v = struct('endpoint', 10 * raw.v_T, 'trailing60', 10 * raw.v_bar);
end

function row = make_table_row_local(datasetName, s)
    row = struct();
    row.Dataset = string(datasetName);
    row.Endpoint_Est = s.endpoint;
    row.Trailing60_Est = s.trailing60;
end

%% ========================================================================
function figFile = save_xv_plot_local(sample, est, stats, opts, figDir)
    t = sample.t(:);
    alpha = est.smoother.alphaSmooth;
    x = alpha(1, :)';
    v = 10 * alpha(2, :)';
    y = sample.y(:);
    obs = sample.yObserved(:);

    keep = t >= opts.figure.plotStart;
    fig = figure('Visible', opts.figure.visible, 'Color', 'w', ...
        'Units', 'pixels', 'Position', [100 100 opts.figure.width opts.figure.height]);

    ax1 = subplot(2,1,1);
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
    plot(ax1, t(keep & obs), y(keep & obs), '-', 'Color', [0.72 0.78 1.00], 'LineWidth', 0.8);
    plot(ax1, t(keep), x(keep), 'k-', 'LineWidth', 1.8);
    ylabel(ax1, '$^\circ$C', 'Interpreter', 'latex');
    title(ax1, sprintf('%s: velocity-RW level $x(t)$', sample.name), 'Interpreter', 'latex');
    legend(ax1, {'Observed','Smoothed XV'}, 'Location', 'northwest', 'Box', 'off');

    ax2 = subplot(2,1,2);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
    plot(ax2, t(keep), v(keep), 'k-', 'LineWidth', 1.8);
    yline(ax2, 0, 'k:', 'HandleVisibility', 'off');
    ylabel(ax2, '$^\circ$C decade$^{-1}$', 'Interpreter', 'latex');
    xlabel(ax2, 'Year', 'Interpreter', 'latex');
    title(ax2, sprintf('%s: velocity-RW warming rate $v(t)$', sample.name), 'Interpreter', 'latex');
    text(ax2, 0.02, 0.92, sprintf('$v(T)=%.3f$', stats.v.endpoint), ...
        'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', 10);

    linkaxes([ax1 ax2], 'x');
    xlim(ax1, [min(t(keep)) max(t(keep))]);
    set([ax1 ax2], 'TickLabelInterpreter', 'latex', 'FontSize', 10);

    figFile = fullfile(figDir, sprintf('XVRobustness_XV_%s.%s', sample.name, lower(opts.figure.format)));
    save_figure_local(fig, figFile, opts.figure.format);
    close(fig);
end

function save_figure_local(fig, fileName, fmt)
    switch lower(fmt)
        case 'png'
            print(fig, fileName, '-dpng', '-r300');
        case 'pdf'
            print(fig, fileName, '-dpdf', '-r300');
        otherwise
            saveas(fig, fileName);
    end
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

function A = symmetrize_local(A)
    A = 0.5 * (A + A');
end

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

function m = nanmean_local(x)
    x = x(isfinite(x));
    if isempty(x), m = NaN; else, m = mean(x); end
end

function s = nanstd_local(x)
    x = x(isfinite(x));
    if numel(x) < 2, s = NaN; else, s = std(x); end
end

function v = nanvar_local(x)
    x = x(isfinite(x));
    if isempty(x), v = 1; else, v = var(x); end
end

function out = merge_structs_local(a, b)
    out = a;
    if isempty(b)
        return;
    end
    f = fieldnames(b);
    for i = 1:numel(f)
        fn = f{i};
        if isstruct(b.(fn)) && isfield(out, fn) && isstruct(out.(fn))
            out.(fn) = merge_structs_local(out.(fn), b.(fn));
        else
            out.(fn) = b.(fn);
        end
    end
end
