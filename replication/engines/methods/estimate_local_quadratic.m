
function out = estimate_local_quadratic(dataFile, userOpts)
% ESTIMATE_LOCAL_QUADRATIC: trend-aware partially linear adjustment.
% Non-parametric benchmark for the temperature-kinematics paper.
%
% The partially linear estimator uses the following specification.
%   - local polynomial regression of degree 2
%   - partial out smooth time trends from BOTH y and Z before estimating beta
%   - re-estimate beta in every bootstrap draw (fixed-design linear map)
%   - endpoint and 5-year trailing-average estimates of x(t), v(t), a(t)
%   - trend-preserving moving-block bootstrap for uncertainty
%   - direct 90% percentile intervals (5th and 95th bootstrap percentiles)
%
% The script writes CSV files with the information needed for Tables 1-3:
%   Kernel_Table1_Position.csv
%   Kernel_Table2_Velocity.csv
%   Kernel_Table3_Acceleration.csv
%   Kernel_AverageAccelerationSince2015.csv

    assert(nargin >= 1 && ~isempty(dataFile),'Provide the path to the analysis panel.');
    if nargin < 2
        userOpts = struct();
    end

    opts = default_options_local();
    opts = merge_structs_local(opts, userOpts);
    validateattributes(opts.bootstrap.ciAlpha, {'numeric'}, ...
        {'scalar','real','finite','>',0,'<',1});

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
    data = prepare_base_data_local(raw, opts);

    datasets = opts.datasets(:)';
    nD = numel(datasets);

    posRows = cell(nD, 1);
    velRows = cell(nD, 1);
    accRows = cell(nD, 1);
    detail = struct();

    fprintf('\n============================================================\n');
    fprintf('Estimating kernel benchmark model\n');
    fprintf('============================================================\n');
    fprintf('Bandwidth h (years): %.3f\n', opts.kernel.bandwidthYears);
    fprintf('Kernel: %s\n', opts.kernel.name);
    fprintf('Bootstrap replications: %d\n', opts.bootstrap.B);
    fprintf('Bootstrap block length (observations): %d\n', opts.bootstrap.blockLengthMonths);
    fprintf('Bootstrap CI level: %.1f%%\n', 100 * (1 - opts.bootstrap.ciAlpha));
    fprintf('Use parfor: %s\n', tf_to_yesno_local(opts.bootstrap.useParallel));

    for i = 1:nD
        datasetName = datasets{i};
        fprintf('\n------------------------------------------------------------\n');
        fprintf('Dataset: %s\n', datasetName);
        fprintf('------------------------------------------------------------\n');

        sample = build_sample_local(data, datasetName, opts);

        est = fit_kernel_model_local(sample, opts);

        fprintf('  x(T)      = %10.4f\n', est.position.endpoint);
        fprintf('  xbar(60)  = %10.4f\n', est.position.trailing60);
        fprintf('  v(T)      = %10.4f degC/decade\n', est.velocity.endpoint);
        fprintf('  vbar(60)  = %10.4f degC/decade\n', est.velocity.trailing60);
        fprintf('  a(T)      = %10.4f degC/decade^2\n', est.acceleration.endpoint);
        fprintf('  abar(60)  = %10.4f degC/decade^2\n', est.acceleration.trailing60);

        posRows{i} = make_row_local(datasetName, est.position, opts.bootstrap.ciAlpha);
        velRows{i} = make_row_local(datasetName, est.velocity, opts.bootstrap.ciAlpha);
        accRows{i} = make_row_local(datasetName, est.acceleration, opts.bootstrap.ciAlpha);

        detail.(datasetName) = est;
    end

    posTable = struct2table(vertcat(posRows{:}));
    velTable = struct2table(vertcat(velRows{:}));
    accTable = struct2table(vertcat(accRows{:}));

    writetable(posTable, fullfile(opts.outputDir, 'Kernel_Table1_Position.csv'));
    writetable(velTable, fullfile(opts.outputDir, 'Kernel_Table2_Velocity.csv'));
    writetable(accTable, fullfile(opts.outputDir, 'Kernel_Table3_Acceleration.csv'));
    avgRows = cell(nD,1);
    for i=1:nD
        E=detail.(datasets{i});
        avgRows{i}=table(string(datasets{i}),E.actualT0,E.sample.t(end), ...
            E.averageAcceleration.estimate,E.averageAcceleration.ci_low, ...
            E.averageAcceleration.ci_high,1-opts.bootstrap.ciAlpha, ...
            'VariableNames',{'Dataset','ActualT0','EndTime','Estimate','CI_Low','CI_High','CILevel'});
    end
    averageTable=vertcat(avgRows{:});
    writetable(averageTable,fullfile(opts.outputDir,'Kernel_AverageAccelerationSince2015.csv'));
    save(fullfile(opts.outputDir, 'Kernel_Benchmark_Output.mat'), 'detail', 'opts', '-v7.3');

    out = struct();
    out.options = opts;
    out.positionTable = posTable;
    out.velocityTable = velTable;
    out.accelerationTable = accTable;
    out.detail = detail;
    out.averageAccelerationTable = averageTable;

    fprintf('\nSaved kernel benchmark output to:\n  %s\n', opts.outputDir);
end

%% ========================================================================
function opts = default_options_local()

    opts = struct();
    opts.seed = 12345;
    opts.outputDir = fullfile(pwd, 'results_local_quadratic');
    opts.datasets = {'NASA','NOAA','HadCRU','Berk','ERA5'};
    opts.rebaseData = true;

    opts.covariates = struct();
    opts.covariates.names = {'NINO34','volc2','ssn'};
    opts.covariates.standardize = true;
    opts.covariates.maxLags = [12 12 12];
    opts.covariates.almonDegrees = [2 1 1];

    opts.kernel = struct();
    opts.kernel.bandwidthYears = 25;
    opts.kernel.name = 'tricube';
    opts.kernel.degree = 2;
    opts.kernel.trailingMonths = 60;
    opts.kernel.averageStartYear = 2015;

    opts.bootstrap = struct();
    opts.bootstrap.B = 2000;
    opts.bootstrap.blockLengthMonths = 24;
    opts.bootstrap.centerResiduals = true;
    opts.bootstrap.ciAlpha = 0.10;
    opts.bootstrap.useParallel = false;
    opts.bootstrap.nWorkers = [];
    opts.bootstrap.verboseEvery = 250;
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
function sample = build_sample_local(data, datasetName, opts)

    y = data.temp.(datasetName);
    [Z, zNames] = build_almon_design_local(data.covariates, opts.covariates.names, ...
        opts.covariates.maxLags, opts.covariates.almonDegrees);

    valid = isfinite(y) & all(isfinite(Z), 2);
    idx = find(valid);

    if isempty(idx)
        error('No valid observations for dataset %s.', datasetName);
    end

    use = idx(1):idx(end);

    sample = struct();
    sample.name = datasetName;
    sample.t = data.t(use);
    sample.y = y(use);
    sample.Z = Z(use, :);
    sample.ZNames = zNames;
    sample.n = numel(sample.y);

    tail = opts.kernel.trailingMonths;
    if sample.n < tail
        error('Sample too short for trailing %d months in dataset %s.', tail, datasetName);
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

function est = fit_kernel_model_local(sample, opts)

    tail = opts.kernel.trailingMonths;
    targetIdx = (sample.n-tail+1):sample.n;
    startIdx=find(sample.t>=opts.kernel.averageStartYear,1);
    assert(~isempty(startIdx) && startIdx<sample.n,'Average acceleration window is empty.');
    targets=[sample.t(targetIdx);sample.t(startIdx)];
    projFull=precompute_local_poly_projections_local(sample.t,sample.t, ...
        opts.kernel.bandwidthYears,opts.kernel.name,opts.kernel.degree);
    S=zeros(sample.n);
    for j=1:sample.n, S(j,:)=projFull{j}.A(1,:); end
    % beta = argmin ||(I-S)(y-Z*beta)||^2. Columns are scaled for stability.
    Ztilde=sample.Z-S*sample.Z;
    scale=sqrt(sum(Ztilde.^2,1));
    assert(all(scale>0),'Covariates are unidentified after smoothing over time.');
    Zscaled=Ztilde./scale;
    assert(rank(Zscaled)==size(Zscaled,2),'Residualized covariates are rank deficient.');
    [Q,R]=qr(Zscaled,0);
    betaMap=(R\(Q'-Q'*S))./scale';
    betaHat=betaMap*sample.y;
    yAdj=sample.y-sample.Z*betaHat;
    xFull=S*yAdj;
    proj=precompute_local_poly_projections_local(sample.t,targets, ...
        opts.kernel.bandwidthYears,opts.kernel.name,opts.kernel.degree);
    m=numel(targets); L=zeros(3*m,sample.n);
    for j=1:m
        A=proj{j}.A;
        A=A-(A*sample.Z)*betaMap;
        L(j,:)=A(1,:); L(m+j,:)=10*A(2,:); L(2*m+j,:)=200*A(3,:);
    end
    theta=L*sample.y;
    xTail=theta(1:tail); vTail=theta(m+(1:tail)); aTail=theta(2*m+(1:tail));
    point=build_estimands_local(xTail,vTail,aTail);
    dt=sample.t(end)-sample.t(startIdx);
    % v is degC/decade; multiply its change/year by 10 for degC/decade^2.
    avgL=10*(L(m+tail,:)-L(2*m,:))/dt;
    statsL=[L(tail,:);mean(L(1:tail,:),1); ...
        L(m+tail,:);mean(L(m+(1:tail),:),1); ...
        L(2*m+tail,:);mean(L(2*m+(1:tail),:),1);avgL];
    fitted=xFull+sample.Z*betaHat;
    resid=sample.y-fitted;
    if opts.bootstrap.centerResiduals, resid=resid-mean(resid); end
    bootStats=nan(opts.bootstrap.B,7);
    % statsL incorporates betaMap: applying it to yStar refits BOTH stages.
    if opts.bootstrap.useParallel
        parfor b=1:opts.bootstrap.B
            bootStats(b,:)=(statsL*(fitted+moving_block_bootstrap_local(resid,opts.bootstrap.blockLengthMonths)))';
        end
    else
        for b=1:opts.bootstrap.B
            bootStats(b,:)=(statsL*(fitted+moving_block_bootstrap_local(resid,opts.bootstrap.blockLengthMonths)))';
            if opts.bootstrap.verboseEvery>0 && mod(b,opts.bootstrap.verboseEvery)==0
                fprintf('  Bootstrap %d / %d\n',b,opts.bootstrap.B);
            end
        end
    end
    est=struct('sample',sample,'betaHat',betaHat,'yAdj',yAdj,'xFull',xFull, ...
        'tailTargets',sample.t(targetIdx),'xTail',xTail,'vTail',vTail,'aTail',aTail, ...
        'bootstrapStats',bootStats(:,1:6),'actualT0',sample.t(startIdx));
    est.adjustment='partially linear: residualize both y and Z on smooth time';
    est.diagnostics.residualizedDesignCondition=cond(Zscaled);
    est.diagnostics.betaNormalEquation=norm(Ztilde'*(yAdj-S*yAdj),inf);
    est.position=summarize_bootstrap_estimand_local(point.position,bootStats(:,1),bootStats(:,2),opts.bootstrap.ciAlpha);
    est.velocity=summarize_bootstrap_estimand_local(point.velocity,bootStats(:,3),bootStats(:,4),opts.bootstrap.ciAlpha);
    est.acceleration=summarize_bootstrap_estimand_local(point.acceleration,bootStats(:,5),bootStats(:,6),opts.bootstrap.ciAlpha);
    est.averageAcceleration=struct('estimate',avgL*sample.y, ...
        'ci_low',quantile_local(bootStats(:,7),opts.bootstrap.ciAlpha/2), ...
        'ci_high',quantile_local(bootStats(:,7),1-opts.bootstrap.ciAlpha/2), ...
        'bootstrap',bootStats(:,7));
    % Keep small maps for independent validation and alternative contrasts.
    est.betaMap=betaMap; est.statsMap=statsL;
end

%% ========================================================================
function proj = precompute_local_poly_projections_local(t, targets, h, kernelName, degree)

    if degree ~= 2
        error('This implementation currently expects local polynomial degree = 2.');
    end

    n = numel(t);
    m = numel(targets);
    proj = cell(m, 1);

    for j = 1:m
        tau = targets(j);
        z = t - tau;
        u = z ./ h;
        w = kernel_weights_local(u, kernelName);

        use=w>0; sw=sqrt(w(use));
        X=[ones(sum(use),1),u(use),u(use).^2];
        Xw=X.*sw;
        assert(rank(Xw)==3,'Insufficient local support for quadratic fit.');
        [Q,R]=qr(Xw,0);
        A=zeros(3,n);
        A(:,use)=(R\Q').*sw';
        A(2,:)=A(2,:)/h; A(3,:)=A(3,:)/h^2;

        pj = struct();
        pj.target = tau;
        pj.A = A;
        proj{j} = pj;
    end
end

%% ========================================================================
function w = kernel_weights_local(u, kernelName)

    switch lower(strtrim(kernelName))
        case 'gaussian'
            w = exp(-0.5 * u.^2);
        case 'epanechnikov'
            w = 0.75 * max(1 - u.^2, 0);
        case 'tricube'
            a = abs(u);
            w = (a < 1) .* (1 - a.^3).^3;
        otherwise
            error('Unknown kernel: %s', kernelName);
    end
end

%% ========================================================================
function [xTail, vTail, aTail] = apply_local_poly_projections_local(proj, y)

    m = numel(proj);
    xTail = nan(m, 1);
    vTail = nan(m, 1);
    aTail = nan(m, 1);

    for j = 1:m
        b = proj{j}.A * y;
        xTail(j) = b(1);
        vTail(j) = b(2) * 10;
        aTail(j) = 2 * b(3) * 100;
    end
end

%% ========================================================================
function xFit = fitted_local_from_targets_local(t, targets, xTail)

    n = numel(t);
    xFit = nan(n, 1);

    for i = 1:n
        ti = t(i);
        if ti <= targets(1)
            xFit(i) = xTail(1);
        elseif ti >= targets(end)
            xFit(i) = xTail(end);
        else
            j = find(targets <= ti, 1, 'last');
            t0 = targets(j);
            t1 = targets(j + 1);
            x0 = xTail(j);
            x1 = xTail(j + 1);
            lam = (ti - t0) / (t1 - t0);
            xFit(i) = (1 - lam) * x0 + lam * x1;
        end
    end
end

%% ========================================================================
function stats = build_estimands_local(xTail, vTail, aTail)

    stats = struct();

    stats.position = struct();
    stats.position.endpoint = xTail(end);
    stats.position.trailing60 = mean(xTail);

    stats.velocity = struct();
    stats.velocity.endpoint = vTail(end);
    stats.velocity.trailing60 = mean(vTail);

    stats.acceleration = struct();
    stats.acceleration.endpoint = aTail(end);
    stats.acceleration.trailing60 = mean(aTail);
end

%% ========================================================================
function s = summarize_bootstrap_estimand_local(pointStruct, endpointBoot, trailingBoot, alpha)

    s = struct();
    s.endpoint = pointStruct.endpoint;
    s.trailing60 = pointStruct.trailing60;

    s.endpoint_ci_low = quantile_local(endpointBoot, alpha / 2);
    s.endpoint_ci_high = quantile_local(endpointBoot, 1 - alpha / 2);
    s.endpoint_p = (1 + sum(endpointBoot <= 0)) / (numel(endpointBoot) + 1);

    s.trailing60_ci_low = quantile_local(trailingBoot, alpha / 2);
    s.trailing60_ci_high = quantile_local(trailingBoot, 1 - alpha / 2);
    s.trailing60_p = (1 + sum(trailingBoot <= 0)) / (numel(trailingBoot) + 1);
end

%% ========================================================================
function row = make_row_local(datasetName, s, ciAlpha)
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
    row.CILevel = 1 - ciAlpha;
end

%% ========================================================================
function xStar = moving_block_bootstrap_local(x, blockLen)

    x = x(:);
    n = numel(x);
    starts = 1:(n - blockLen + 1);
    nBlocks = ceil(n / blockLen);

    draw = zeros(nBlocks, 1);
    for i = 1:nBlocks
        draw(i) = starts(randi(numel(starts)));
    end

    xStar = nan(nBlocks * blockLen, 1);
    k = 0;
    for i = 1:nBlocks
        ii = draw(i):(draw(i) + blockLen - 1);
        xStar((k+1):(k+blockLen)) = x(ii);
        k = k + blockLen;
    end
    xStar = xStar(1:n);
end

%% ========================================================================
function q = quantile_local(x, p)

    x = x(isfinite(x));
    if isempty(x)
        q = NaN;
        return;
    end

    x = sort(x);
    n = numel(x);
    pos = 1 + (n - 1) * p;
    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        q = x(lo);
    else
        w = pos - lo;
        q = (1 - w) * x(lo) + w * x(hi);
    end
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
function A = regularize_matrix_local(A)
    A = 0.5 * (A + A');
    ridge = 1e-10 * max(1, trace(abs(A)));
    A = A + ridge * eye(size(A));
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
