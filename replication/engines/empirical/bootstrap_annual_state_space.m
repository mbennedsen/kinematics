function out=bootstrap_annual_state_space(action,varargin)
% Uses a frozen snapshot of the validated annual likelihood/profile helpers.
switch action
    case 'record', out=make_record(varargin{:});
    case 'run', out=run_record(varargin{:});
    otherwise, error('Unknown action');
end
end

function record=make_record(input,opts,est,dataset,seed,inputSHA)
sample=build_sample_local(prepare_data_local(readtable(input),opts),dataset,opts);
assert(~opts.rebaseData && sample.t(end)==2025 && sample.t(end-10)==2015);
assert(est.paramsHat.q==opts.stateSpace.fixedQ && est.setup.fixedQ==opts.stateSpace.fixedQ);
[nll,filt,pred,params]=negloglik_local(est.thetaHat,sample,est.setup,opts.stateSpace,false);
sm=rts_smoother_local(filt,pred,state_transition_local(params,est.setup));
assert(abs(nll+est.logLik)<1e-5 && max(abs(sm.alphaSmooth(:)-est.smoother.alphaSmooth(:)))<1e-7);
record=struct('dataset',dataset,'sample',sample,'pointEst',est,'stateOpts',opts.stateSpace, ...
 'sigma',sqrt(opts.stateSpace.fixedQ),'seed',seed,'source',inputSHA);
end

function out=run_record(record,B,nWorkers,outputDir,batchSize)
ds=record.dataset; folder=outputDir;
if ~exist(folder,'dir'),mkdir(folder);end
config=struct('version',1,'B',B,'ciAlpha',.10,'method','posteriorResidual-bootstrapError-fixedSigma');
checkpoint=fullfile(folder,'checkpoint.mat'); T=record.sample.T;
if exist(checkpoint,'file')
    c=load(checkpoint);
    assert(isequaln(c.config,config) && isequaln(c.record,record), ...
        'Checkpoint inputs/settings differ. Use a new output directory.');
    errors=c.errors; pathErrors=c.pathErrors; done=c.done; messages=c.messages;
else
    errors=nan(B,7);pathErrors=nan(T,3,B);done=false(B,1);messages=strings(B,1);
end
pending=find(~done);
for first=1:batchSize:numel(pending)
    ids=pending(first:min(first+batchSize-1,numel(pending)));
    E=nan(numel(ids),7);P=nan(T,3,numel(ids));M=strings(numel(ids),1);
    if nWorkers>0
        parfor j=1:numel(ids)
            [E(j,:),P(:,:,j),M(j)]=one_draw(record,ids(j));
        end
    else
        for j=1:numel(ids)
            [E(j,:),P(:,:,j),M(j)]=one_draw(record,ids(j));
        end
    end
    errors(ids,:)=E;pathErrors(:,:,ids)=P;messages(ids)=M;done(ids)=true;
    save(fullfile(folder,'checkpoint.tmp.mat'),'config','record','errors','pathErrors','done','messages','-v7.3');
    movefile(fullfile(folder,'checkpoint.tmp.mat'),checkpoint,'f');
    fprintf('%s: completed %d/%d; successful %d\n',ds,sum(done),B,sum(all(isfinite(errors),2)));
end
valid=all(isfinite(errors),2) & reshape(all(all(isfinite(pathErrors),1),2),[],1);
failures=table((1:B)',valid,messages,'VariableNames',{'Replication','Success','Message'});
writetable(failures,fullfile(folder,'BootstrapDiagnostics.csv'));
assert(sum(valid)>=ceil(.95*B),'Fewer than 95%% successful draws; inspect BootstrapDiagnostics.csv.');
point=stats(record.pointEst.smoother.alphaSmooth,record.sample.t);
ci=point-prctile(errors(valid,:),[95 5],1);se=std(errors(valid,:),0,1);
labels=["Level";"Rate";"Acceleration";"Trailing5Level";"Trailing5Rate";"Trailing5Acceleration";"AverageAcceleration2015"];
units=["deg C";"deg C/decade";"deg C/decade^2";"deg C";"deg C/decade";"deg C/decade^2";"deg C/decade^2"];
summary=table(repmat(string(ds),7,1),labels,units,point',se',ci(1,:)',ci(2,:)', ...
    repmat(record.sigma,7,1),repmat(sum(valid),7,1),repmat(B,7,1), ...
    'VariableNames',{'Dataset','Estimand','Units','Estimate','SE','CI_Low','CI_High','FixedSigma','SuccessfulDraws','RequestedDraws'});
writetable(summary,fullfile(folder,'AnnualSS_Intervals.csv'));
pointPath=record.pointEst.smoother.alphaSmooth(1:3,:)'.*[1 10 100];
quant=prctile(pathErrors(:,:,valid),[5 95],3);
low=pointPath-quant(:,:,2);high=pointPath-quant(:,:,1);
paths=table(record.sample.t,pointPath(:,1),low(:,1),high(:,1),pointPath(:,2),low(:,2),high(:,2), ...
    pointPath(:,3),low(:,3),high(:,3),'VariableNames',{'Year','Level','Level_Low','Level_High','Rate','Rate_Low','Rate_High','Acceleration','Acceleration_Low','Acceleration_High'});
writetable(paths,fullfile(folder,'AnnualSS_PathBands.csv'));
% Table-9-compatible endpoint/trailing estimates; all intervals are 90%.
names={'SS_Annual_Table1_Position.csv','SS_Annual_Table2_Velocity.csv','SS_Annual_Table3_Acceleration.csv'};
for k=1:3
    tab=table(string(ds),point(k),se(k),ci(1,k),ci(2,k),point(k+3),se(k+3),ci(1,k+3),ci(2,k+3), ...
        'VariableNames',{'Dataset','Endpoint_Est','Endpoint_SE','Endpoint_CI_Low','Endpoint_CI_High','TrailingN_Est','TrailingN_SE','TrailingN_CI_Low','TrailingN_CI_High'});
    writetable(tab,fullfile(folder,names{k}));
end
out=struct('summary',summary,'paths',paths,'config',config,'fixedSigma',record.sigma,'conditionalOnMonthlySigma',true);
save(fullfile(folder,'AnnualSS_Bootstrap_Output.mat'),'out','record','errors','pathErrors','valid','-v7.3');
disp(summary);
end

function [e,p,message]=one_draw(record,b)
T=record.sample.T;e=nan(1,7);p=nan(T,3);message="";
rng(record.seed+b,'twister');
try
    sample=record.sample;point=record.pointEst;setup=point.setup;opts=record.stateOpts;
    assert(setup.fixedQ==opts.fixedQ && point.paramsHat.q==opts.fixedQ);
    F=state_transition_local(point.paramsHat,setup);H=setup.H(:)';
    truth=simulation_smoother_draw_local(point.filter,point.smoother,F);
    for t=2:T
        shock=draw_state_noise_local(point.paramsHat,setup,t);
        truth(4:end,t)=F(4:end,4:end)*truth(4:end,t-1)+shock(4:end);
    end
    sample.y=nan(T,1);
    for t=1:T
        if sample.yObserved(t), sample.y(t)=H*truth(:,t)+sample.W(t,:)*point.paramsHat.beta;end
    end
    if setup.estimateRho
        starts=[atanh(point.paramsHat.rho/setup.rhoBound),log(point.paramsHat.sig2u);0,log(point.paramsHat.sig2u)];
    else
        starts=[log(point.paramsHat.sig2u);log(point.paramsHat.sig2u)+1];
    end
    best=Inf;bestZ=[];
    searchOpts=optimset('Display','off','MaxIter',1200,'MaxFunEvals',4000,'TolX',2e-7,'TolFun',1e-9);
    % Retry optimization on the SAME pseudo-sample; never replace a failed draw.
    for j=1:size(starts,1)
        [z,nll,flag]=fminsearch(@(z)profile_objective(z,opts.fixedQ,sample,setup,opts),starts(j,:)',searchOpts);
        if flag>0 && isfinite(nll) && nll<best, best=nll;bestZ=z;end
    end
    assert(~isempty(bestZ),'Both nuisance optimizations failed');
    [~,theta]=profile_objective(bestZ,opts.fixedQ,sample,setup,opts);
    [check,filt,pred,params]=negloglik_local(theta,sample,setup,opts,false);
    assert(params.q==opts.fixedQ && abs(check-best)<1e-6);
    sm=rts_smoother_local(filt,pred,state_transition_local(params,setup));
    e=stats(sm.alphaSmooth,sample.t)-stats(truth,sample.t);
    p=(sm.alphaSmooth(1:3,:)'-truth(1:3,:)').*[1 10 100];
    assert(all(isfinite(e)) && all(isfinite(p(:))));
catch err
    e=nan(1,7);p=nan(T,3);message=string(err.message);
end
end

function v=stats(alpha,t)
scaled=alpha(1:3,:)'.*[1 10 100];
i0=find(t>=2015,1);assert(~isempty(i0) && t(i0)==2015 && t(end)>2015);
av=(scaled(end,2)-scaled(i0,2))/((t(end)-t(i0))/10);
v=[scaled(end,:),mean(scaled(max(1,end-4):end,:),1),av];
end

function [nll,theta]=profile_objective(z,q,sample,setup,opts)
% For fixed covariance parameters, innovations are linear in beta. GLS
% concentrates beta exactly, including beta-dependent initial state mean.
if any(~isfinite(z)) || any(abs(z)>35) || ~isfinite(q) || q>1 || q<0
    nll=1e10+sum(z.^2); theta=[]; return;
end
setup.fixedQ=q;
theta=[zeros(setup.K,1);z(:)];
[~,filt,~,params]=negloglik_local(theta,sample,setup,opts,false);
F=state_transition_local(params,setup); H=setup.H(:)';
B=zeros(setup.stateDim,setup.K); B(1,:)=-sample.W(1,:);
X=zeros(sample.T,setup.K);
for t=1:sample.T
    if t>1, B=F*B; end
    X(t,:)=sample.W(t,:)+H*B;
    if sample.yObserved(t), B=B-filt.K(:,t)*X(t,:); end
end
use=sample.yObserved;
sd=sqrt(filt.F(use));
Xw=X(use,:)./sd; yw=filt.v(use)./sd;
beta=Xw\yw;
res=yw-Xw*beta;
nll=.5*sum(log(2*pi)+log(filt.F(use))+res.^2);
theta(1:setup.K)=beta;
end

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
