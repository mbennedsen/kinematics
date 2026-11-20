
function out = estimate_monthly_state_space(dataFile, userOpts)
% Fit the monthly kinematic state-space model with optional covariates.
% Point estimates use the Kalman smoother. Bootstrap options govern parameter
% refits, path intervals and per-draw checkpoints for the coverage experiment.

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
    fprintf('Estimating monthly state-space estimates\n');
    fprintf('============================================================\n');
    fprintf('Datasets: %s\n', strjoin(opts.datasets, ', '));
    fprintf('Bootstrap replications: %d\n', opts.bootstrap.B);
    fprintf('Include parameter uncertainty: %s\n', tf_to_yesno_local(opts.bootstrap.includeParameterUncertainty));
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

        fprintf('start t = %.4f\n', sample.t(1));
        fprintf('end   t = %.4f\n', sample.t(end));
        fprintf('T = %d, observed y = %d, internal missing y = %d\n', ...
            sample.T, sum(sample.yObserved), sum(~sample.yObserved));
        fprintf('Using maxLags      = [%s]\n', num2str(opts.covariates.maxLags));
        fprintf('Using almonDegrees = [%s]\n', num2str(opts.covariates.almonDegrees));

        ssOpts = state_space_options_for_dataset_local(opts.stateSpace, datasetName);

        thetaStart = [];
        if isfield(ssOpts, 'thetaStartByDataset') && isstruct(ssOpts.thetaStartByDataset) && ...
                isfield(ssOpts.thetaStartByDataset, datasetName)
            thetaStart = ssOpts.thetaStartByDataset.(datasetName);
        end

        if isfield(opts.bootstrap,'checkpointDir') && ~isempty(opts.bootstrap.checkpointDir)
            pf=fullfile(opts.bootstrap.checkpointDir,'PointFit.mat');
            key=opts.bootstrap.checkpointKey;
            if isfile(pf)
                L=load(pf,'key','pointEst');assert(strcmp(key,L.key),'Point checkpoint mismatch.');pointEst=L.pointEst;
            else
                pointEst=estimate_individual_model_local(sample,ssOpts,thetaStart);
                assert(pointEst.exitflag>0,'Outer point fit did not converge.');
                tmp=[tempname(opts.bootstrap.checkpointDir) '.mat'];save(tmp,'key','pointEst','-v7');movefile(tmp,pf,'f');
            end
        else
            pointEst = estimate_individual_model_local(sample, ssOpts, thetaStart);
        end
        pointStates = pointEst.smoother.alphaSmooth(1:3, :);
        pointStatsRaw = make_estimands_local(pointStates, opts.trailingMonths);
        pointStats = scale_estimands_local(pointStatsRaw);
        if average_acceleration_enabled_local(opts)
            pointAverageAcceleration = make_average_acceleration_detail_local( ...
                pointEst.smoother.alphaSmooth, sample.t, opts.averageAcceleration);
        else
            pointAverageAcceleration = struct();
        end

        optsDataset = opts;
        optsDataset.stateSpace = ssOpts;
        boot = run_bootstrap_local(sample, pointEst, optsDataset);

        posRows{j} = make_table_row_local(datasetName, pointStats.x, boot.position);
        velRows{j} = make_table_row_local(datasetName, pointStats.v, boot.velocity);
        accRows{j} = make_table_row_local(datasetName, pointStats.a, boot.acceleration);

        detail.(datasetName) = struct();
        detail.(datasetName).sample = sample.meta;
        % Retain the estimation sample as well as its compact metadata.  This
        % makes fitted-model simulation checks reproducible without changing
        % the estimator or bootstrap calculations.
        detail.(datasetName).sampleData = sample;
        detail.(datasetName).pointEstimate = pointEst;
        detail.(datasetName).pointEstimandsRaw = pointStatsRaw;
        detail.(datasetName).pointEstimandsScaled = pointStats;
        detail.(datasetName).pointAverageAcceleration = pointAverageAcceleration;
        detail.(datasetName).bootstrap = boot;

        fprintf('  x(T)      = %10.4f\n', pointStats.x.endpoint);
        fprintf('  xbar(60)  = %10.4f\n', pointStats.x.trailing60);
        fprintf('  v(T)      = %10.4f degC/decade\n', pointStats.v.endpoint);
        fprintf('  vbar(60)  = %10.4f degC/decade\n', pointStats.v.trailing60);
        fprintf('  a(T)      = %10.4f degC/decade^2\n', pointStats.a.endpoint);
        fprintf('  abar(60)  = %10.4f degC/decade^2\n', pointStats.a.trailing60);
        if average_acceleration_enabled_local(opts)
            fprintf('  abar(%g,T) = %10.4f degC/decade^2\n', ...
                opts.averageAcceleration.T0, pointAverageAcceleration.estimate);
        end

        if opts.bootstrap.B > 0
            ciPct = round(100 * (1 - opts.bootstrap.ciAlpha));
            fprintf('  %d%% CI x(T):     [%10.4f, %10.4f]\n', ciPct, boot.position.endpoint_ci_low, boot.position.endpoint_ci_high);
            fprintf('  %d%% CI v(T):     [%10.4f, %10.4f]\n', ciPct, boot.velocity.endpoint_ci_low, boot.velocity.endpoint_ci_high);
            fprintf('  %d%% CI a(T):     [%10.4f, %10.4f]\n', ciPct, boot.acceleration.endpoint_ci_low, boot.acceleration.endpoint_ci_high);
            fprintf('  Bootstrap SE x(T), v(T), a(T): [%10.4f, %10.4f, %10.4f]\n', ...
                boot.position.endpoint_se, boot.velocity.endpoint_se, boot.acceleration.endpoint_se);
        end
    end

    posTable = struct2table(vertcat(posRows{:}));
    velTable = struct2table(vertcat(velRows{:}));
    accTable = struct2table(vertcat(accRows{:}));

    writetable(posTable, fullfile(opts.outputDir, 'SS_Table1_Position.csv'));
    writetable(velTable, fullfile(opts.outputDir, 'SS_Table2_Velocity.csv'));
    writetable(accTable, fullfile(opts.outputDir, 'SS_Table3_Acceleration.csv'));
    save(fullfile(opts.outputDir, 'SS_Table_Bootstrap_Output.mat'), 'detail', 'opts', '-v7.3');

    out = struct();
    out.options = opts;
    out.positionTable = posTable;
    out.velocityTable = velTable;
    out.accelerationTable = accTable;
    out.detail = detail;

    fprintf('\nSaved SS table outputs to:\n  %s\n', opts.outputDir);
end

%% ========================================================================
function opts = default_options_local()

    opts = struct();
    opts.seed = 12345;
    opts.outputDir = fullfile(pwd, 'results_monthly_state_space');
    opts.datasets = {'NASA','NOAA','HadCRU','Berk','ERA5'};
    opts.rebaseData = true;
    opts.trailingMonths = 60;
    % Optional supplementary estimand used by the coverage-study extension.
    % It is deliberately not written to Tables 1--3.
    opts.averageAcceleration = struct();
    opts.averageAcceleration.enabled = false;
    opts.averageAcceleration.T0 = 2015;

    opts.covariates = struct();
    opts.covariates.names = {'NINO34','volc2','ssn'};
    opts.covariates.standardize = true;
    opts.covariates.maxLags = [12 12 12];
    opts.covariates.almonDegrees = [2 1 1];
    opts.covariates.includeMonthlyDummies = false;

    opts.stateSpace = struct();
    opts.stateSpace.dt = 1/12;
    opts.stateSpace.rhoBound = 0.995;
    opts.stateSpace.maBound = 0.995;
    opts.stateSpace.arOrder = 1;
    opts.stateSpace.maOrder = 1;
    opts.stateSpace.enforceARStationarity = true;
    opts.stateSpace.enforceMAInvertibility = false;
    opts.stateSpace.useMonthlyVariance = true;
    opts.stateSpace.varianceBreakYears = [];
    opts.stateSpace.jitter = 1e-8;
    opts.stateSpace.minVariance = 1e-10;
    opts.stateSpace.qLowerBound = 0;
    opts.stateSpace.fixedQ = [];
    opts.stateSpace.fixedQByDataset = struct();
    opts.stateSpace.optimizerMaxIter = 600;
    opts.stateSpace.optimizerTol = 1e-6;
    opts.stateSpace.verbose = false;
    opts.stateSpace.computeHessian = false;

    opts.bootstrap = struct();
    opts.bootstrap.B = 499;
    % Enable full-path error bands when requested.
    opts.bootstrap.pathBands = false;
    opts.bootstrap.method = 'stateSimulation';
    % 'anchoredDraw' (default): parameter draws are re-smoothed on y.
    % 'bootstrapError': simulate a known latent path, re-estimate on y*,
    % and use the bootstrap distribution of estimator minus true state.
    opts.bootstrap.inferenceMode = 'anchoredDraw';
    % 'recursive' implements the standardized-innovation recursion in
    % Rodríguez and Ruiz (2012, Eq. 12).  'fixedPrediction' is retained for
    % the fixed-prediction innovation bootstrap.
    opts.bootstrap.innovationReconstruction = 'fixedPrediction';
    opts.bootstrap.includeParameterUncertainty = true;
    opts.bootstrap.drawInitialState = true;
    % Initialization for a forward parametric state bootstrap.  The default
    % reproduces the historical behaviour.  'fixedInitial' fixes alpha(1) at
    % its fitted smoothed mean; 'conditionalPrefix' holds observations through
    % anchorYear fixed, draws alpha(anchor) from its filtered distribution,
    % and simulates every later state and observation innovation forward.
    opts.bootstrap.parametricInitialization = 'smoothedInitial';
    opts.bootstrap.anchorYear = 1970;
    opts.bootstrap.startAtOriginalTheta = true;
    opts.bootstrap.ciAlpha = 0.05;
    opts.bootstrap.maxRefitTries = 2;
    opts.bootstrap.useParallel = false;
    opts.bootstrap.nWorkers = [];

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
function boot = run_bootstrap_local(sample, pointEst, opts)

    B = opts.bootstrap.B;
    useAverageAcceleration = average_acceleration_enabled_local(opts);
    usePathBands = opts.bootstrap.pathBands;
    if usePathBands && (~is_bootstrap_error_mode_local(opts) || ...
            ~opts.bootstrap.includeParameterUncertainty)
        error('Full path bands require bootstrapError with parameter re-estimation.');
    end
    if B <= 0
        boot = summarize_bootstrap_local([], [], [], opts.bootstrap.ciAlpha);
        if useAverageAcceleration
            boot.averageAcceleration = empty_average_acceleration_bootstrap_local(opts.averageAcceleration);
        end
        return;
    end

    if (is_bootstrap_error_mode_local(opts) || is_conditional_pmse_mode_local(opts)) && ...
            ~opts.bootstrap.includeParameterUncertainty
        error('bootstrapError and conditionalPMSE inference require parameter re-estimation.');
    end
    if useAverageAcceleration && ~is_bootstrap_error_mode_local(opts)
        error(['averageAcceleration is currently implemented for the local ', ...
            'bootstrap-error procedure only.']);
    end

    fprintf('  Bootstrapping (%d reps', B);
    if opts.bootstrap.includeParameterUncertainty
        fprintf(', re-estimating parameters');
        if is_bootstrap_error_mode_local(opts)
            fprintf('; estimator-error pivots');
        elseif is_conditional_pmse_mode_local(opts)
            fprintf('; conditional PMSE');
        end
    elseif is_kalman_innovation_bootstrap_local(opts)
        fprintf(', fixed parameters via Kalman innovation bootstrap');
    else
        fprintf(', fixed parameters via simulation smoother');
    end
    if opts.bootstrap.useParallel
        fprintf(', parfor');
    end
    fprintf(') ...\n');

    if ~opts.bootstrap.includeParameterUncertainty
        if is_kalman_innovation_bootstrap_local(opts)
            [repsX, repsV, repsA] = run_fixed_parameter_innovation_bootstrap_local(sample, pointEst, opts);
        else
            [repsX, repsV, repsA] = run_conditional_uncertainty_local(thetaHat_from_pointest_local(pointEst), sample, pointEst.setup, opts);
        end
        boot = summarize_bootstrap_local(repsX, repsV, repsA, opts.bootstrap.ciAlpha);
        return;
    end

    if is_conditional_pmse_mode_local(opts)
        [means, withinVariances] = run_conditional_pmse_bootstrap_local(sample, pointEst, opts);
        pointScaled = scale_estimands_local(make_estimands_local( ...
            pointEst.smoother.alphaSmooth(1:3, :), opts.trailingMonths));
        boot = summarize_conditional_pmse_local(means, withinVariances, pointScaled, ...
            opts.bootstrap.ciAlpha);
        return;
    end

    repsX = nan(B, 2);
    repsV = nan(B, 2);
    repsA = nan(B, 2);
    repsAverageAcceleration = nan(B, 1);
    pathErrors = cell(B, 1);

    if opts.bootstrap.useParallel
        parfor b = 1:B
            [xRow, vRow, aRow, averageAccelerationRow, pathError] = ...
                run_saved_rep_local(b, sample, pointEst, opts);
            pathErrors{b} = pathError;
            repsX(b, :) = xRow;
            repsV(b, :) = vRow;
            repsA(b, :) = aRow;
            repsAverageAcceleration(b) = averageAccelerationRow;

            fprintf('    %d / %d\n', b, B);
        end
    else
        for b = 1:B
            %if mod(b, max(1, round(B/10))) == 0 || b == 1 || b == B
                fprintf('    %d / %d\n', b, B);
            %end

            [xRow, vRow, aRow, averageAccelerationRow, pathError] = ...
                run_saved_rep_local(b, sample, pointEst, opts);
            pathErrors{b} = pathError;
            repsX(b, :) = xRow;
            repsV(b, :) = vRow;
            repsA(b, :) = aRow;
            repsAverageAcceleration(b) = averageAccelerationRow;
        end
    end

    repsX = repsX(all(isfinite(repsX),2), :);
    repsV = repsV(all(isfinite(repsV),2), :);
    repsA = repsA(all(isfinite(repsA),2), :);

    if is_bootstrap_error_mode_local(opts)
        pointRaw = make_estimands_local(pointEst.smoother.alphaSmooth(1:3, :), opts.trailingMonths);
        pointScaled = scale_estimands_local(pointRaw);
        boot = summarize_bootstrap_error_local(repsX, repsV, repsA, pointScaled, opts.bootstrap.ciAlpha);
        if usePathBands
            boot.pathBands = bootstrap_path_intervals( ...
                pointEst.smoother.alphaSmooth, cat(3, pathErrors{:}), opts.bootstrap.ciAlpha);
        end
        if useAverageAcceleration
            pointAverageAcceleration = average_acceleration_from_state_path_local( ...
                pointEst.smoother.alphaSmooth, sample.t, opts.averageAcceleration);
            boot.averageAcceleration = summarize_bootstrap_error_scalar_local( ...
                repsAverageAcceleration, pointAverageAcceleration, opts.bootstrap.ciAlpha, ...
                opts.averageAcceleration, sample.t);
        end
    else
        boot = summarize_bootstrap_local(repsX, repsV, repsA, opts.bootstrap.ciAlpha);
    end
end

%% ========================================================================
function tf = is_kalman_innovation_bootstrap_local(opts)

    method = 'stateSimulation';
    if isfield(opts.bootstrap, 'method') && ~isempty(opts.bootstrap.method)
        method = char(opts.bootstrap.method);
    end
    tf = any(strcmpi(method, {'kalmanInnovation', 'innovation', 'stofferWall', 'stoffer-wall'}));
end

%% ========================================================================
function tf = is_bootstrap_error_mode_local(opts)

    mode = 'anchoredDraw';
    if isfield(opts.bootstrap, 'inferenceMode') && ~isempty(opts.bootstrap.inferenceMode)
        mode = char(opts.bootstrap.inferenceMode);
    end
    tf = any(strcmpi(mode, {'bootstrapError', 'error', 'basicError'}));
end

%% ========================================================================
function tf = is_conditional_pmse_mode_local(opts)

    mode = 'anchoredDraw';
    if isfield(opts.bootstrap, 'inferenceMode') && ~isempty(opts.bootstrap.inferenceMode)
        mode = char(opts.bootstrap.inferenceMode);
    end
    tf = any(strcmpi(mode, {'conditionalPMSE', 'conditionalPMSEnormal', 'rodriguezRuiz'}));
end

%% ========================================================================
function [repsX, repsV, repsA] = run_fixed_parameter_innovation_bootstrap_local(sample, pointEst, opts)

    B = opts.bootstrap.B;
    trailingMonths = opts.trailingMonths;
    repsX = nan(B, 2);
    repsV = nan(B, 2);
    repsA = nan(B, 2);
    thetaHat = thetaHat_from_pointest_local(pointEst);

    if opts.bootstrap.useParallel
        parfor b = 1:B
            bsSample = simulate_kalman_innovation_bootstrap_sample_local(sample, pointEst, opts);
            condEst = evaluate_given_theta_local(bsSample, opts.stateSpace, thetaHat);
            bsStats = scale_estimands_local(make_estimands_local(condEst.smoother.alphaSmooth(1:3, :), trailingMonths));
            repsX(b, :) = [bsStats.x.endpoint, bsStats.x.trailing60];
            repsV(b, :) = [bsStats.v.endpoint, bsStats.v.trailing60];
            repsA(b, :) = [bsStats.a.endpoint, bsStats.a.trailing60];
        end
    else
        for b = 1:B
            if mod(b, max(1, round(B/10))) == 0 || b == 1 || b == B
                fprintf('    %d / %d\n', b, B);
            end
            bsSample = simulate_kalman_innovation_bootstrap_sample_local(sample, pointEst, opts);
            condEst = evaluate_given_theta_local(bsSample, opts.stateSpace, thetaHat);
            bsStats = scale_estimands_local(make_estimands_local(condEst.smoother.alphaSmooth(1:3, :), trailingMonths));
            repsX(b, :) = [bsStats.x.endpoint, bsStats.x.trailing60];
            repsV(b, :) = [bsStats.v.endpoint, bsStats.v.trailing60];
            repsA(b, :) = [bsStats.a.endpoint, bsStats.a.trailing60];
        end
    end
end


%% ========================================================================
function [repsX, repsV, repsA] = run_conditional_uncertainty_local(thetaHat, sample, setup, opts)

    B = opts.bootstrap.B;
    trailingMonths = opts.trailingMonths;

    [~, filt, pred, params] = negloglik_local(thetaHat, sample, setup, opts.stateSpace, false);
    F = state_transition_local(params, setup);
    sm = rts_smoother_local(filt, pred, F);

    repsX = nan(B, 2);
    repsV = nan(B, 2);
    repsA = nan(B, 2);

    if opts.bootstrap.useParallel
        parfor b = 1:B
            alphaDraw = simulation_smoother_draw_local(filt, sm, F);
            bsStats = scale_estimands_local(make_estimands_local(alphaDraw(1:3, :), trailingMonths));
            repsX(b, :) = [bsStats.x.endpoint, bsStats.x.trailing60];
            repsV(b, :) = [bsStats.v.endpoint, bsStats.v.trailing60];
            repsA(b, :) = [bsStats.a.endpoint, bsStats.a.trailing60];
        end
    else
        for b = 1:B
            if mod(b, max(1, round(B/10))) == 0 || b == 1 || b == B
                fprintf('    %d / %d\n', b, B);
            end

            alphaDraw = simulation_smoother_draw_local(filt, sm, F);
            bsStats = scale_estimands_local(make_estimands_local(alphaDraw(1:3, :), trailingMonths));

            repsX(b, :) = [bsStats.x.endpoint, bsStats.x.trailing60];
            repsV(b, :) = [bsStats.v.endpoint, bsStats.v.trailing60];
            repsA(b, :) = [bsStats.a.endpoint, bsStats.a.trailing60];
        end
    end
end


%% ========================================================================
function [bsSample, alphaTrue] = simulate_bootstrap_sample_local(sample, pointEst, opts)

    params = pointEst.paramsHat;
    setup = pointEst.setup;
    F = state_transition_local(params, setup);
    H = setup.H(:)';
    T = sample.T;
    m = setup.stateDim;

    alpha = nan(m, T);
    y = nan(T, 1);

    initialization = 'smoothedInitial';
    if isfield(opts.bootstrap, 'parametricInitialization') && ...
            ~isempty(opts.bootstrap.parametricInitialization)
        initialization = char(opts.bootstrap.parametricInitialization);
    end

    switch lower(initialization)
        case {'smoothedinitial', 'smoothed'}
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
            simulationStart = 2;

        case {'fixedinitial', 'fixed1880'}
            % Conditional on the fitted state at the first observed month.
            alpha(:, 1) = pointEst.smoother.alphaSmooth(:, 1);
            if sample.yObserved(1)
                y(1) = H * alpha(:, 1) + sample.W(1, :) * params.beta;
            end
            simulationStart = 2;

        case {'conditionalprefix', 'conditional1970', 'prefix'}
            anchor = bootstrap_anchor_index_local(sample, opts.bootstrap);
            % Preserve the realised pre-anchor record and draw only from the
            % filtered state distribution given that record.  Filtering, not
            % smoothing, avoids conditioning the anchor state on future data.
            y(1:anchor) = sample.y(1:anchor);
            alpha(:, anchor) = mvnrnd_local(pointEst.filter.aFilt(:, anchor), ...
                regularize_cov_local(pointEst.filter.PFilt(:, :, anchor)));
            simulationStart = anchor + 1;

        otherwise
            error('Unknown parametricInitialization: %s', initialization);
    end

    for t = simulationStart:T
        shock = draw_state_noise_local(params, setup, t);
        alpha(:, t) = F * alpha(:, t-1) + shock;
        if sample.yObserved(t)
            y(t) = H * alpha(:, t) + sample.W(t, :) * params.beta;
        end
    end

    bsSample = sample;
    bsSample.y = y;
    bsSample.yObserved = ~isnan(y);
    alphaTrue = alpha;
end

%% ========================================================================
function anchor = bootstrap_anchor_index_local(sample, bootstrapOpts)

    anchorYear = 1970;
    if isfield(bootstrapOpts, 'anchorYear') && ~isempty(bootstrapOpts.anchorYear)
        anchorYear = bootstrapOpts.anchorYear;
    end
    validateattributes(anchorYear, {'numeric'}, {'scalar', 'finite'});
    anchor = find(sample.t >= anchorYear, 1, 'first');
    if isempty(anchor) || anchor >= sample.T
        error('bootstrap.anchorYear must lie inside the estimation sample before its endpoint.');
    end
    if ~sample.yObserved(anchor)
        error('The selected bootstrap anchor month has a missing observation.');
    end
end

%% ========================================================================
function [bsSample, alphaTrue] = simulate_bootstrap_sample_dispatch_local(sample, pointEst, opts)

    method = 'stateSimulation';
    if isfield(opts.bootstrap, 'method') && ~isempty(opts.bootstrap.method)
        method = char(opts.bootstrap.method);
    end

    switch lower(method)
        case {'statesimulation', 'parametric', 'parametricstate'}
            [bsSample, alphaTrue] = simulate_bootstrap_sample_local(sample, pointEst, opts);
        case {'posteriorresidual', 'localposterior', 'conditionalposterior'}
            [bsSample, alphaTrue] = simulate_posterior_residual_bootstrap_sample_local(sample, pointEst, opts);
        case {'kalmaninnovation', 'innovation', 'stofferwall', 'stoffer-wall'}
            bsSample = simulate_kalman_innovation_bootstrap_sample_local(sample, pointEst, opts);
            alphaTrue = [];
        otherwise
            error('Unknown bootstrap method: %s', method);
    end
end

%% ========================================================================
function [bsSample, alphaTrue] = simulate_posterior_residual_bootstrap_sample_local(sample, pointEst, opts)
% Conditional/local pseudo-sample generator.  The x-v-a trajectory is a
% simulation-smoother draw given the observed record; a fresh ARMA residual
% path is then generated around it.  This avoids long-horizon unconditional
% integrated-RW trajectories while retaining a known latent truth for an
% estimator-error bootstrap.

    %#ok<INUSD> opts retained for a common simulator interface
    params = pointEst.paramsHat;
    setup = pointEst.setup;
    F = state_transition_local(params, setup);
    H = setup.H(:)';
    T = sample.T;

    alphaTrue = simulation_smoother_draw_local(pointEst.filter, pointEst.smoother, F);
    residIdx = 4:setup.stateDim;
    for t = 2:T
        shock = draw_state_noise_local(params, setup, t);
        alphaTrue(residIdx, t) = F(residIdx, residIdx) * alphaTrue(residIdx, t-1) + shock(residIdx);
    end

    y = nan(T, 1);
    for t = 1:T
        if sample.yObserved(t)
            y(t) = H * alphaTrue(:, t) + sample.W(t, :) * params.beta;
        end
    end

    bsSample = sample;
    bsSample.y = y;
    bsSample.yObserved = ~isnan(y);
end

%% ========================================================================
function bsSample = simulate_kalman_innovation_bootstrap_sample_local(sample, pointEst, opts)
% Resample standardized one-step-ahead Kalman innovations and reconstruct a
% bootstrap observation sequence through the fitted Kalman recursions.

    params = pointEst.paramsHat;
    setup = pointEst.setup;
    H = setup.H(:)';
    T = sample.T;
    reconstruction = 'fixedPrediction';
    if nargin >= 3 && isfield(opts.bootstrap, 'innovationReconstruction') && ...
            ~isempty(opts.bootstrap.innovationReconstruction)
        reconstruction = char(opts.bootstrap.innovationReconstruction);
    end

    zPool = pointEst.filter.vStd(sample.yObserved & isfinite(pointEst.filter.vStd));
    zPool = zPool(:);
    if isempty(zPool)
        error('Cannot run Kalman innovation bootstrap: no finite standardized innovations.');
    end
    zPool = zPool - mean(zPool);
    sdZ = std(zPool);
    if isfinite(sdZ) && sdZ > 0
        zPool = zPool ./ sdZ;
    end

    yStar = nan(T, 1);

    if any(strcmpi(reconstruction, {'fixedPrediction', 'fixedPred', 'prediction'}))
        for t = 1:T
            if sample.yObserved(t)
                Ft = max(pointEst.filter.F(t), setup.jitter);
                muObs = H * pointEst.filter.aPred(:, t) + sample.W(t, :) * params.beta;
                zStar = zPool(randi(numel(zPool)));
                yStar(t) = muObs + sqrt(Ft) * zStar;
            end
        end
        bsSample = sample;
        bsSample.y = yStar;
        bsSample.yObserved = sample.yObserved & isfinite(yStar);
        return;
    end

    if ~any(strcmpi(reconstruction, {'recursive', 'kalmanRecursive'}))
        error('Unknown innovation reconstruction mode: %s', reconstruction);
    end

    Fmat = state_transition_local(params, setup);
    aPrev = pointEst.pred.a0;
    PPrev = pointEst.pred.P0;

    for t = 1:T
        if t == 1
            aPred = aPrev;
            PPred = PPrev;
        else
            Qt = state_noise_local(params, setup, t);
            aPred = Fmat * aPrev;
            PPred = Fmat * PPrev * Fmat' + Qt;
        end
        PPred = regularize_cov_local(PPred);

        if sample.yObserved(t)
            muObs = H * aPred + sample.W(t, :) * params.beta;
            Ft = H * PPred * H' + observation_noise_variance_local(params, setup, t) + setup.jitter;
            Ft = max(Ft, setup.jitter);
            Kt = (PPred * H') / Ft;

            zStar = zPool(randi(numel(zPool)));
            vStar = sqrt(Ft) * zStar;
            yStar(t) = muObs + vStar;

            aFilt = aPred + Kt * vStar;
            PFilt = PPred - Kt * Ft * Kt';
            PFilt = regularize_cov_local(PFilt);
        else
            aFilt = aPred;
            PFilt = PPred;
        end

        aPrev = aFilt;
        PPrev = PFilt;
    end

    bsSample = sample;
    bsSample.y = yStar;
    bsSample.yObserved = sample.yObserved & isfinite(yStar);
end

%% ========================================================================
function shock = draw_state_noise_local(params, setup, t)

    qxva = params.q * setup.SigmaXVA;
    xvaShock = mvnrnd_local(zeros(3,1), qxva);

    scale = 1;
    if setup.useMonthlyVariance
        scale = scale * exp(params.gamma(setup.month(t)));
    end
    if numel(params.delta) >= setup.era(t)
        scale = scale * exp(params.delta(setup.era(t)));
    end
    Rt = max(params.sig2u * scale, setup.minVariance);
    uShock = sqrt(Rt) * randn;
    residShock = zeros(setup.residDim, 1);
    if setup.residDim > 0
        residShock(1) = uShock;
    end
    if setup.maOrder > 0
        residShock(setup.arOrder + 1) = uShock;
    end

    shock = [xvaShock; residShock];
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
function [xRow, vRow, aRow, averageAccelerationRow, pathError] = run_parameter_bootstrap_rep_local(sample, pointEst, opts)

    xRow = [NaN, NaN];
    vRow = [NaN, NaN];
    aRow = [NaN, NaN];
    averageAccelerationRow = NaN;
    pathError = [];
    if opts.bootstrap.pathBands
        pathError = nan(sample.T, 4);
    end

    success = false;
    thetaStart = pointEst.thetaHat;
    alphaTrue = [];

    for tryNum = 1:opts.bootstrap.maxRefitTries
        [bsSample, alphaTrue] = simulate_bootstrap_sample_dispatch_local(sample, pointEst, opts);

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
        if is_bootstrap_error_mode_local(opts)
            if isempty(alphaTrue)
                error('bootstrapError inference requires stateSimulation pseudo-samples.');
            end
            estimateStats = scale_estimands_local(make_estimands_local( ...
                bsEst.smoother.alphaSmooth(1:3, :), opts.trailingMonths));
            trueStats = scale_estimands_local(make_estimands_local( ...
                alphaTrue(1:3, :), opts.trailingMonths));

            xRow = [estimateStats.x.endpoint - trueStats.x.endpoint, ...
                    estimateStats.x.trailing60 - trueStats.x.trailing60];
            vRow = [estimateStats.v.endpoint - trueStats.v.endpoint, ...
                    estimateStats.v.trailing60 - trueStats.v.trailing60];
            aRow = [estimateStats.a.endpoint - trueStats.a.endpoint, ...
                    estimateStats.a.trailing60 - trueStats.a.trailing60];
            if opts.bootstrap.pathBands
                pathError = kinematic_path_values(bsEst.smoother.alphaSmooth) - ...
                    kinematic_path_values(alphaTrue);
            end
            if average_acceleration_enabled_local(opts)
                estimateAverageAcceleration = average_acceleration_from_state_path_local( ...
                    bsEst.smoother.alphaSmooth, sample.t, opts.averageAcceleration);
                trueAverageAcceleration = average_acceleration_from_state_path_local( ...
                    alphaTrue, sample.t, opts.averageAcceleration);
                averageAccelerationRow = estimateAverageAcceleration - trueAverageAcceleration;
            end
            return;
        end

        condEst = evaluate_given_theta_local(sample, opts.stateSpace, thetaStar);
        F = state_transition_local(condEst.paramsHat, condEst.setup);
        alphaDraw = simulation_smoother_draw_local(condEst.filter, condEst.smoother, F);
        bsStats = scale_estimands_local(make_estimands_local(alphaDraw(1:3, :), opts.trailingMonths));

        xRow = [bsStats.x.endpoint, bsStats.x.trailing60];
        vRow = [bsStats.v.endpoint, bsStats.v.trailing60];
        aRow = [bsStats.a.endpoint, bsStats.a.trailing60];
    catch
        % leave NaN rows
    end
end

%% ========================================================================
function [means, withinVariances] = run_conditional_pmse_bootstrap_local(sample, pointEst, opts)
% Conditional-PMSE bootstrap of Rodríguez and Ruiz (2012).  Each
% pseudo-sample is used only to obtain thetaStar; the Kalman smoother is then
% evaluated on the original sample.  The uncertainty decomposition is
%   E_{thetaStar}[Var(psi | y, thetaStar)]
%     + Var_{thetaStar}[E(psi | y, thetaStar)].

    B = opts.bootstrap.B;
    means = nan(B, 6);
    withinVariances = nan(B, 6);

    if opts.bootstrap.useParallel
        parfor b = 1:B
            [means(b, :), withinVariances(b, :)] = ...
                run_conditional_pmse_rep_local(sample, pointEst, opts);
            fprintf('    %d / %d\n', b, B);
        end
    else
        for b = 1:B
            fprintf('    %d / %d\n', b, B);
            [means(b, :), withinVariances(b, :)] = ...
                run_conditional_pmse_rep_local(sample, pointEst, opts);
        end
    end

    ok = all(isfinite(means), 2) & all(isfinite(withinVariances), 2);
    means = means(ok, :);
    withinVariances = withinVariances(ok, :);
end

%% ========================================================================
function [meanRow, withinRow] = run_conditional_pmse_rep_local(sample, pointEst, opts)

    meanRow = nan(1, 6);
    withinRow = nan(1, 6);
    success = false;
    thetaStart = pointEst.thetaHat;

    for tryNum = 1:opts.bootstrap.maxRefitTries
        [bsSample, ~] = simulate_bootstrap_sample_dispatch_local(sample, pointEst, opts);
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
        meanRow = estimands_row_local(condEst.smoother.alphaSmooth, opts.trailingMonths);
        withinRow = conditional_estimand_variances_local(condEst, opts.trailingMonths);
    catch ME
        warning('Conditional-PMSE bootstrap replicate failed: %s', ME.message);
        % Leave this replicate unavailable.
    end
end

%% ========================================================================
function values = estimands_row_local(alphaPath, trailingMonths)

    stats = scale_estimands_local(make_estimands_local(alphaPath(1:3, :), trailingMonths));
    values = [stats.x.endpoint, stats.x.trailing60, ...
              stats.v.endpoint, stats.v.trailing60, ...
              stats.a.endpoint, stats.a.trailing60];
end

%% ========================================================================
function variances = conditional_estimand_variances_local(condEst, trailingMonths)
% Exact conditional variances for endpoint and trailing-window means under
% the Gaussian smoothing distribution.  The off-diagonal smoothed state
% covariances are propagated with the RTS gains, avoiding an additional
% simulation-smoother Monte Carlo layer.

    filt = condEst.filter;
    sm = condEst.smoother;
    F = state_transition_local(condEst.paramsHat, condEst.setup);
    T = size(sm.alphaSmooth, 2);
    L = min(trailingMonths, T);
    first = T - L + 1;
    m = size(sm.alphaSmooth, 1);

    gains = cell(T - 1, 1);
    for t = first:T-1
        Ptt = regularize_cov_local(filt.PFilt(:, :, t));
        Ptpred = regularize_cov_local(filt.PPred(:, :, t + 1));
        gains{t} = (Ptt * F') / Ptpred;
    end

    sumCov = zeros(m, m);
    for s = first:T
        Cts = regularize_cov_local(sm.PSmooth(:, :, s));
        sumCov = sumCov + Cts;
        for t = s-1:-1:first
            Cts = gains{t} * Cts;
            sumCov = sumCov + Cts + Cts';
        end
    end

    endpointVar = diag(sm.PSmooth(1:3, 1:3, T));
    trailingVar = diag(sumCov(1:3, 1:3)) / (L^2);
    scale = [1; 10; 100];
    endpointVar = (scale .^ 2) .* endpointVar;
    trailingVar = (scale .^ 2) .* trailingVar;
    variances = [endpointVar(1), trailingVar(1), ...
                 endpointVar(2), trailingVar(2), ...
                 endpointVar(3), trailingVar(3)];
    variances = max(real(variances), 0);
end

%% ========================================================================
function summary = summarize_conditional_pmse_local(means, withinVariances, pointStats, alpha)

    if isempty(means)
        summary = summarize_bootstrap_local([], [], [], alpha);
        summary.inferenceMode = 'conditionalPMSE';
        summary.conditionalPMSE = struct('validParameterDraws', 0, ...
            'withinVariances', [], 'parameterVariances', [], 'pmseVariances', []);
        return;
    end

    if size(means, 1) > 1
        parameterVariances = var(means, 0, 1);
    else
        parameterVariances = zeros(1, size(means, 2));
    end
    withinVariancesMean = mean(withinVariances, 1);
    pmseVariances = max(withinVariancesMean + parameterVariances, 0);
    pointRow = [pointStats.x.endpoint, pointStats.x.trailing60, ...
                pointStats.v.endpoint, pointStats.v.trailing60, ...
                pointStats.a.endpoint, pointStats.a.trailing60];
    z = sqrt(2) * erfcinv(alpha);

    summary = struct();
    summary.position = summarize_conditional_pmse_two_stats_local(pointRow(1:2), pmseVariances(1:2), z);
    summary.velocity = summarize_conditional_pmse_two_stats_local(pointRow(3:4), pmseVariances(3:4), z);
    summary.acceleration = summarize_conditional_pmse_two_stats_local(pointRow(5:6), pmseVariances(5:6), z);
    summary.inferenceMode = 'conditionalPMSE';
    summary.replicates = struct();
    summary.replicates.position = means(:, 1:2);
    summary.replicates.velocity = means(:, 3:4);
    summary.replicates.acceleration = means(:, 5:6);
    summary.conditionalPMSE = struct();
    summary.conditionalPMSE.validParameterDraws = size(means, 1);
    summary.conditionalPMSE.conditionalMeans = means;
    summary.conditionalPMSE.withinVariances = withinVariances;
    summary.conditionalPMSE.meanWithinVariances = withinVariancesMean;
    summary.conditionalPMSE.parameterVariances = parameterVariances;
    summary.conditionalPMSE.pmseVariances = pmseVariances;
end

%% ========================================================================
function s = summarize_conditional_pmse_two_stats_local(pointValues, variances, z)

    se = sqrt(max(variances, 0));
    s = struct();
    s.endpoint_ci_low = pointValues(1) - z * se(1);
    s.endpoint_ci_high = pointValues(1) + z * se(1);
    s.endpoint_se = se(1);
    s.endpoint_p = NaN;
    s.trailing60_ci_low = pointValues(2) - z * se(2);
    s.trailing60_ci_high = pointValues(2) + z * se(2);
    s.trailing60_se = se(2);
    s.trailing60_p = NaN;
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
function summary = summarize_bootstrap_error_local(errX, errV, errA, pointStats, alpha)
% Basic parametric-bootstrap intervals from estimator-minus-truth errors.
% The stored replicates remain errors, not pseudo-estimates.

    summary = struct();
    summary.position = summarize_bootstrap_error_two_stats_local(errX, pointStats.x, alpha);
    summary.velocity = summarize_bootstrap_error_two_stats_local(errV, pointStats.v, alpha);
    summary.acceleration = summarize_bootstrap_error_two_stats_local(errA, pointStats.a, alpha);
    summary.inferenceMode = 'bootstrapError';

    summary.replicates = struct();
    summary.replicates.position = errX;
    summary.replicates.velocity = errV;
    summary.replicates.acceleration = errA;
end

%% ========================================================================
function tf = average_acceleration_enabled_local(opts)

    tf = isfield(opts, 'averageAcceleration') && isstruct(opts.averageAcceleration) && ...
        isfield(opts.averageAcceleration, 'enabled') && ...
        logical(opts.averageAcceleration.enabled);
end

%% ========================================================================
function detail = make_average_acceleration_detail_local(alphaPath, t, averageOpts)

    [estimate, index0, actualT0, elapsedDecades] = ...
        average_acceleration_from_state_path_local(alphaPath, t, averageOpts);
    detail = struct('estimate', estimate, 'T0', averageOpts.T0, ...
        'index0', index0, 'actualT0', actualT0, 'endTime', t(end), ...
        'elapsedDecades', elapsedDecades);
end

%% ========================================================================
function [estimate, index0, actualT0, elapsedDecades] = ...
        average_acceleration_from_state_path_local(alphaPath, t, averageOpts)
% Average acceleration between the first sample date at or after T0 and T:
% [v(T)-v(T0)] / ((T-T0)/10), in deg C per decade squared.

    if ~isfield(averageOpts, 'T0') || ~isscalar(averageOpts.T0) || ...
            ~isfinite(averageOpts.T0)
        error('averageAcceleration.T0 must be a finite scalar year.');
    end
    t = t(:);
    index0 = find(t >= averageOpts.T0, 1, 'first');
    if isempty(index0) || index0 >= numel(t)
        error('averageAcceleration.T0 must lie inside the estimation sample before T.');
    end
    actualT0 = t(index0);
    elapsedDecades = (t(end) - actualT0) / 10;
    if ~isfinite(elapsedDecades) || elapsedDecades <= 0
        error('Average-acceleration elapsed time is not positive.');
    end
    v = 10 * alphaPath(2, :);
    estimate = (v(end) - v(index0)) / elapsedDecades;
end

%% ========================================================================
function s = empty_average_acceleration_bootstrap_local(averageOpts)

    s = struct('estimate', NaN, 'endpoint_ci_low', NaN, ...
        'endpoint_ci_high', NaN, 'endpoint_se', NaN, 'endpoint_p', NaN, ...
        'replicates', zeros(0, 1), 'T0', averageOpts.T0, ...
        'actualT0', NaN, 'endTime', NaN, 'elapsedDecades', NaN);
end

%% ========================================================================
function s = summarize_bootstrap_error_scalar_local(errors, pointEstimate, alpha, averageOpts, t)

    errors = errors(isfinite(errors));
    if isempty(errors)
        s = empty_average_acceleration_bootstrap_local(averageOpts);
        s.estimate = pointEstimate;
        return;
    end
    [~, ~, actualT0, elapsedDecades] = average_acceleration_from_state_path_local( ...
        zeros(3, numel(t)), t, averageOpts);
    lo = 100 * (alpha / 2);
    hi = 100 * (1 - alpha / 2);
    s = struct();
    s.estimate = pointEstimate;
    s.endpoint_ci_low = pointEstimate - prctile(errors, hi);
    s.endpoint_ci_high = pointEstimate - prctile(errors, lo);
    s.endpoint_se = std(errors, 0);
    s.endpoint_p = NaN;
    s.replicates = errors;
    s.T0 = averageOpts.T0;
    s.actualT0 = actualT0;
    s.endTime = t(end);
    s.elapsedDecades = elapsedDecades;
end

%% ========================================================================
function s = summarize_bootstrap_error_two_stats_local(errors, pointObj, alpha)

    errors = errors(all(isfinite(errors), 2), :);
    if isempty(errors)
        s = summarize_two_stats_local([], alpha);
        return;
    end

    lo = 100 * (alpha / 2);
    hi = 100 * (1 - alpha / 2);
    endpointErrors = errors(:, 1);
    trailingErrors = errors(:, 2);

    s = struct();
    s.endpoint_ci_low = pointObj.endpoint - prctile(endpointErrors, hi);
    s.endpoint_ci_high = pointObj.endpoint - prctile(endpointErrors, lo);
    s.endpoint_se = std(endpointErrors, 0);
    % A bootstrap error distribution alone does not define the posterior-style
    % one-sided state probabilities reported for anchored draws.
    s.endpoint_p = NaN;
    s.trailing60_ci_low = pointObj.trailing60 - prctile(trailingErrors, hi);
    s.trailing60_ci_high = pointObj.trailing60 - prctile(trailingErrors, lo);
    s.trailing60_se = std(trailingErrors, 0);
    s.trailing60_p = NaN;
end

%% ========================================================================
function s = summarize_two_stats_local(X, alpha)

    if isempty(X)
        s = struct();
        s.endpoint_ci_low = NaN;
        s.endpoint_ci_high = NaN;
        s.endpoint_se = NaN;
        s.endpoint_p = NaN;
        s.trailing60_ci_low = NaN;
        s.trailing60_ci_high = NaN;
        s.trailing60_se = NaN;
        s.trailing60_p = NaN;
        return;
    end

    X = X(all(isfinite(X), 2), :);
    if isempty(X)
        s = struct();
        s.endpoint_ci_low = NaN;
        s.endpoint_ci_high = NaN;
        s.endpoint_se = NaN;
        s.endpoint_p = NaN;
        s.trailing60_ci_low = NaN;
        s.trailing60_ci_high = NaN;
        s.trailing60_se = NaN;
        s.trailing60_p = NaN;
        return;
    end

    lo = 100 * (alpha / 2);
    hi = 100 * (1 - alpha / 2);

    endpoint = X(:, 1);
    trailing = X(:, 2);

    s = struct();
    s.endpoint_ci_low = prctile(endpoint, lo);
    s.endpoint_ci_high = prctile(endpoint, hi);
    s.endpoint_se = std(endpoint, 0);
    s.endpoint_p = empirical_one_sided_p_local(endpoint, 0);
    s.trailing60_ci_low = prctile(trailing, lo);
    s.trailing60_ci_high = prctile(trailing, hi);
    s.trailing60_se = std(trailing, 0);
    s.trailing60_p = empirical_one_sided_p_local(trailing, 0);
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
    row.Endpoint_SE = bootObj.endpoint_se;
    row.Endpoint_CI_Low = bootObj.endpoint_ci_low;
    row.Endpoint_CI_High = bootObj.endpoint_ci_high;
    row.Endpoint_P = bootObj.endpoint_p;
    row.Trailing60_Est = pointObj.trailing60;
    row.Trailing60_SE = bootObj.trailing60_se;
    row.Trailing60_CI_Low = bootObj.trailing60_ci_low;
    row.Trailing60_CI_High = bootObj.trailing60_ci_high;
    row.Trailing60_P = bootObj.trailing60_p;
end

%% ========================================================================
function scaled = scale_estimands_local(raw)

    scaled = struct();

    scaled.x = struct();
    scaled.x.endpoint = raw.x_T;
    scaled.x.trailing60 = raw.x_bar;

    scaled.v = struct();
    scaled.v.endpoint = 10 * raw.v_T;
    scaled.v.trailing60 = 10 * raw.v_bar;

    scaled.a = struct();
    scaled.a.endpoint = 100 * raw.a_T;
    scaled.a.trailing60 = 100 * raw.a_bar;
end

%% ========================================================================
function estimands = make_estimands_local(alphaPath, trailingMonths)

    x = alphaPath(1, :);
    v = alphaPath(2, :);
    a = alphaPath(3, :);

    L = min(trailingMonths, size(alphaPath, 2));
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

    data.covariates = struct();
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
function sample = build_sample_local(data, datasetName, opts)

    if isempty(opts.covariates.names)
        almon = struct('W', zeros(numel(data.t),0), 'WNames', {{}}, 'coeffInfo', {{}});
    else
        almon = build_almon_basis_local(data.covariates, opts.covariates.names, ...
            opts.covariates.maxLags, opts.covariates.almonDegrees);
    end

    yFull = data.temp.(datasetName);
    validW = all(isfinite(almon.W), 2);
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
    sample.W = almon.W(use, :);
    sample.WNames = almon.WNames;
    sample.t = data.t(use);
    sample.month = data.month(use);
    if isfield(opts.covariates, 'includeMonthlyDummies') && opts.covariates.includeMonthlyDummies
        [sample.W, sample.WNames] = append_monthly_mean_dummies_local( ...
            sample.W, sample.WNames, sample.month);
    end
    sample.era = make_era_index_local(data.calendarYear(use), opts.stateSpace.varianceBreakYears);
    sample.yObserved = ~isnan(sample.y);
    sample.T = numel(sample.y);
    sample.lagInfo = almon;

    sample.meta = struct();
    sample.meta.dataset = datasetName;
    sample.meta.startDate = sample.t(1);
    sample.meta.endDate = sample.t(end);
    sample.meta.T = sample.T;
    sample.meta.observedY = sum(sample.yObserved);
    sample.meta.internalMissingY = sum(~sample.yObserved);
    sample.meta.WNames = sample.WNames;
end

%% ========================================================================
function almon = build_almon_basis_local(covStruct, covNames, maxLags, degrees)

    nCov = numel(covNames);
    n = numel(covStruct.(covNames{1}));
    K = sum(degrees + 1);

    W = nan(n, K);
    WNames = cell(1, K);
    coeffInfo = cell(1, nCov);

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

        coeffInfo{i} = struct();
        coeffInfo{i}.covariate = covNames{i};
        coeffInfo{i}.maxLag = L;
        coeffInfo{i}.degree = deg;
        coeffInfo{i}.lagGrid = lagGrid;

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

    almon = struct();
    almon.W = W;
    almon.WNames = WNames;
    almon.coeffInfo = coeffInfo;
end

%% ========================================================================
function [W, WNames] = append_monthly_mean_dummies_local(W, WNames, month)
% Add 11 seasonal mean dummies to the observation equation. December is the
% omitted reference month, so coefficients are relative to December.

    month = month(:);
    D = zeros(numel(month), 11);
    for m = 1:11
        D(:, m) = double(month == m);
    end

    W = [W, D];
    dummyNames = arrayfun(@(m) sprintf('month_%02d_vs_Dec', m), 1:11, ...
        'UniformOutput', false);
    WNames = [WNames, dummyNames];
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

    [thetaHat, logLik, exitflag, out, hessian] = optimize_local(thetaStart, sample, setup, opts);
    [~, filt, pred, params] = negloglik_local(thetaHat, sample, setup, opts, false);
    F = state_transition_local(params, setup);
    sm = rts_smoother_local(filt, pred, F);

    est = struct();
    est.thetaHat = thetaHat;
    est.logLik = logLik;
    est.exitflag = exitflag;
    est.optimizerOutput = out;
    est.hessian = hessian;
    est.filter = filt;
    est.pred = pred;
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
    setup.month = sample.month(:);
    setup.era = sample.era(:);
    setup.nEra = max(sample.era(:));
    setup.rhoBound = opts.rhoBound;
    setup.maBound = opts.maBound;
    setup.arOrder = getfield_default_local(opts, 'arOrder', 1);
    setup.maOrder = getfield_default_local(opts, 'maOrder', 1);
    setup.enforceARStationarity = getfield_default_local(opts, 'enforceARStationarity', true);
    setup.enforceMAInvertibility = getfield_default_local(opts, 'enforceMAInvertibility', false);
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
    setup.SigmaXVA = [dt^5/20, dt^4/8, dt^3/6; ...
                      dt^4/8,  dt^3/3, dt^2/2; ...
                      dt^3/6,  dt^2/2, dt      ];

    setup.numGammaFree = 0;
    if setup.useMonthlyVariance
        setup.numGammaFree = 11;
    end
    setup.numDeltaFree = max(setup.nEra - 1, 0);

    setup.residDim = setup.arOrder + setup.maOrder;
    setup.stateDim = 3 + setup.residDim;
    setup.H = zeros(1, setup.stateDim);
    setup.H(1) = 1;
    if setup.residDim > 0
        setup.H(4) = 1;
    end
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
    ar0 = zeros(setup.arOrder, 1);
    if setup.arOrder >= 1
        ar0(1) = rho0;
    end
    ma0 = zeros(setup.maOrder, 1);
    if setup.maOrder >= 1
        ma0(1) = 0.1;
    end

    sig2u0 = max(1e-6, 0.25 * nanvar_local(u));

    gamma0 = zeros(setup.numGammaFree, 1);
    delta0 = zeros(setup.numDeltaFree, 1);

    thetaPieces = {beta0(:); ...
                   atanh(ar0(:) / setup.rhoBound); ...
                   atanh(ma0(:) / setup.maBound)};

    if isempty(setup.fixedQ)
        q0 = max(1e-8, 1e-4 * nanvar_local(u));
        thetaPieces{end+1} = log(max(q0, 1e-10)); %#ok<AGROW>
    end

    thetaPieces = [thetaPieces; {log(max(sig2u0, 1e-10)); gamma0; delta0}];
    theta0 = vertcat(thetaPieces{:});
end

%% ========================================================================
function [thetaHat, logLik, exitflag, out, hessian] = optimize_local(theta0, sample, setup, opts)

    obj = @(th) negloglik_local(th, sample, setup, opts, true);
    f0 = obj(theta0);

    hasFminunc = exist('fminunc', 'file') == 2;
    if hasFminunc
        optimOpts = optimoptions('fminunc', ...
            'Algorithm', 'quasi-newton', ...
            'Display', 'off', ...
            'MaxIterations', opts.optimizerMaxIter, ...
            'MaxFunctionEvaluations', 50000, ...
            'OptimalityTolerance', opts.optimizerTol, ...
            'StepTolerance', 1e-10);
        if isfield(opts, 'computeHessian') && opts.computeHessian
            [thetaHat, fval, exitflag, out, ~, hessian] = fminunc(obj, theta0, optimOpts);
        else
            [thetaHat, fval, exitflag, out] = fminunc(obj, theta0, optimOpts);
            hessian = NaN(numel(theta0));
        end
    else
        optimOpts = optimset('Display', 'off', ...
            'MaxIter', opts.optimizerMaxIter, ...
            'MaxFunEvals', 50000, ...
            'TolX', opts.optimizerTol, ...
            'TolFun', opts.optimizerTol);
        [thetaHat, fval, exitflag, out] = fminsearch(obj, theta0, optimOpts);
        hessian = NaN(numel(theta0));
    end

    if ~isfinite(fval) || (isfinite(f0) && fval > f0)
        thetaHat = theta0;
        fval = f0;
        exitflag = 0;
        if isstruct(out)
            out.retainedStart = true;
        end
        hessian = NaN(numel(theta0));
    elseif isstruct(out)
        out.retainedStart = false;
    end

    logLik = -fval;
end

%% ========================================================================
function [nll, filt, pred, params] = negloglik_local(thetaRaw, sample, setup, opts, onlyNll)

    if nargin < 5
        onlyNll = false;
    end

    params = unpack_params_local(thetaRaw, setup);
    if ~arma_is_admissible_local(params, setup)
        nll = 1e12;
        filt = [];
        pred = [];
        return;
    end
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
    if m >= 4
        PPrev(4:end, 4:end) = 10 * eye(m - 3);
    end

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
            Ft = H * PPred * H' + observation_noise_variance_local(params, setup, t) + setup.jitter;
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

    if setup.arOrder > 0
        params.ar = setup.rhoBound * tanh(thetaRaw(k+1:k+setup.arOrder));
        k = k + setup.arOrder;
    else
        params.ar = zeros(0, 1);
    end
    if setup.maOrder > 0
        params.maCoeffs = setup.maBound * tanh(thetaRaw(k+1:k+setup.maOrder));
        k = k + setup.maOrder;
    else
        params.maCoeffs = zeros(0, 1);
    end
    params.rho = first_or_nan_local(params.ar);
    params.ma = first_or_nan_local(params.maCoeffs);

    if isfield(setup, 'fixedQ') && ~isempty(setup.fixedQ)
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
    m = setup.stateDim;
    F = zeros(m, m);
    F(1:3, 1:3) = [1, dt, 0.5*dt^2; ...
                   0, 1,  dt; ...
                   0, 0,  1];

    r0 = 4;
    p = setup.arOrder;
    q = setup.maOrder;

    if p > 0
        F(r0, r0:(r0+p-1)) = params.ar(:)';
    end
    if q > 0
        e0 = r0 + p;
        F(r0, e0:(e0+q-1)) = params.maCoeffs(:)';
    end

    for i = 2:p
        F(r0+i-1, r0+i-2) = 1;
    end
    for j = 2:q
        F(r0+p+j-1, r0+p+j-2) = 1;
    end
end

%% ========================================================================
function Qt = state_noise_local(params, setup, t)

    Qxva = params.q * setup.SigmaXVA;

    scale = 1;
    if setup.useMonthlyVariance
        scale = scale * exp(params.gamma(setup.month(t)));
    end
    if numel(params.delta) >= setup.era(t)
        scale = scale * exp(params.delta(setup.era(t)));
    end

    Rt = max(params.sig2u * scale, setup.minVariance);
    b = zeros(setup.residDim, 1);
    if setup.residDim > 0
        b(1) = 1;
    end
    if setup.maOrder > 0
        b(setup.arOrder + 1) = 1;
    end
    Qu = Rt * (b * b');

    Qt = blkdiag(Qxva, Qu);
    Qt = symmetrize_local(Qt);
end

%% ========================================================================
function R = observation_noise_variance_local(params, setup, t)

    if setup.residDim > 0
        R = 0;
        return;
    end

    scale = 1;
    if setup.useMonthlyVariance
        scale = scale * exp(params.gamma(setup.month(t)));
    end
    if numel(params.delta) >= setup.era(t)
        scale = scale * exp(params.delta(setup.era(t)));
    end
    R = max(params.sig2u * scale, setup.minVariance);
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
function val = getfield_default_local(S, fieldName, defaultValue)
    if isstruct(S) && isfield(S, fieldName)
        val = S.(fieldName);
    else
        val = defaultValue;
    end
end

%% ========================================================================
function x = first_or_nan_local(v)
    if isempty(v)
        x = NaN;
    else
        x = v(1);
    end
end

%% ========================================================================
function ok = arma_is_admissible_local(params, setup)
    ok = true;
    if setup.enforceARStationarity
        ok = ok && ar_is_stationary_local(params.ar);
    end
    if setup.enforceMAInvertibility
        ok = ok && ma_is_invertible_local(params.maCoeffs);
    end
end

%% ========================================================================
function ok = ar_is_stationary_local(ar)
    ar = ar(:)';
    p = numel(ar);
    if p == 0
        ok = true;
        return;
    end
    if any(~isfinite(ar))
        ok = false;
        return;
    end
    if p == 1
        ok = abs(ar(1)) < 0.999;
        return;
    end
    C = [ar; eye(p-1), zeros(p-1, 1)];
    ok = max(abs(eig(C))) < 0.999;
end

%% ========================================================================
function ok = ma_is_invertible_local(ma)
    ma = ma(:)';
    q = numel(ma);
    if q == 0
        ok = true;
        return;
    end
    if any(~isfinite(ma))
        ok = false;
        return;
    end
    if q == 1
        ok = abs(ma(1)) < 0.999;
        return;
    end
    polyCoef = [fliplr(ma), 1];
    lead = find(abs(polyCoef) > 1e-12, 1, 'first');
    if isempty(lead)
        ok = true;
        return;
    end
    rootsVal = roots(polyCoef(lead:end));
    ok = isempty(rootsVal) || min(abs(rootsVal)) > 1.001;
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

function [xRow,vRow,aRow,averageAccelerationRow,pathError]=run_saved_rep_local(b,sample,pointEst,opts)
% Scheduling-independent RNG and atomic per-draw restart, no estimator changes.
if ~isfield(opts.bootstrap,'checkpointDir') || isempty(opts.bootstrap.checkpointDir)
    [xRow,vRow,aRow,averageAccelerationRow,pathError]=run_parameter_bootstrap_rep_local(sample,pointEst,opts);
    return
end
seed=opts.seed+b; key=opts.bootstrap.checkpointKey;
file=fullfile(opts.bootstrap.checkpointDir,sprintf('draw_%05d.mat',b));
if isfile(file)
 L=load(file);assert(strcmp(L.key,key)&&L.seed==seed,'Bootstrap checkpoint mismatch.');
 xRow=L.xRow;vRow=L.vRow;aRow=L.aRow;averageAccelerationRow=L.averageAccelerationRow;pathError=L.pathError;return
end
rng(seed,'twister');
[xRow,vRow,aRow,averageAccelerationRow,pathError]=run_parameter_bootstrap_rep_local(sample,pointEst,opts);
tmp=[tempname(opts.bootstrap.checkpointDir) '.mat'];
save(tmp,'key','seed','xRow','vRow','aRow','averageAccelerationRow','pathError');movefile(tmp,file,'f');
end
