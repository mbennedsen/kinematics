function out = write_ss_parameter_appendix_tables(inputMatFile, userOpts)
% WRITE_SS_PARAMETER_APPENDIX_TABLES
% Create appendix-ready parameter tables for the individual BM state-space
% models from saved estimation output.
%
% Outputs:
%   SS_ParameterAppendix_Dynamics.csv/.tex
%   SS_ParameterAppendix_BetaBasis.csv
%   SS_ParameterAppendix_ImpliedLagWeights.csv
%   SS_ParameterAppendix_MonthlyVarianceScales.csv
%   SS_ParameterAppendix_Tables.tex
%
% inputMatFile contains the combined monthly point fits.

assert(nargin >= 1 && ~isempty(inputMatFile),'Provide the combined monthly point fits.');
if nargin < 2 || isempty(userOpts)
    userOpts = struct();
end

opts = default_options_local();
opts = merge_structs_local(opts, userOpts);

if exist(inputMatFile, 'file') ~= 2
    error('Input MAT file not found: %s', inputMatFile);
end
if ~exist(opts.outputDir, 'dir')
    mkdir(opts.outputDir);
end

S = load(inputMatFile, 'detail', 'opts');
detail = S.detail;
modelOpts = S.opts;
datasets = fieldnames(detail);

dynTable = make_dynamics_table_local(detail, datasets);
betaTable = make_beta_basis_table_local(detail, datasets);
lagTable = make_lag_weight_table_local(detail, datasets, modelOpts);
monthTable = make_monthly_variance_table_local(detail, datasets);

files = struct();
files.dynamicsCsv = fullfile(opts.outputDir, 'SS_ParameterAppendix_Dynamics.csv');
files.betaCsv = fullfile(opts.outputDir, 'SS_ParameterAppendix_BetaBasis.csv');
files.lagCsv = fullfile(opts.outputDir, 'SS_ParameterAppendix_ImpliedLagWeights.csv');
files.monthCsv = fullfile(opts.outputDir, 'SS_ParameterAppendix_MonthlyVarianceScales.csv');
files.tex = fullfile(opts.outputDir, 'SS_ParameterAppendix_Tables.tex');

writetable(dynTable, files.dynamicsCsv);
writetable(betaTable, files.betaCsv);
writetable(lagTable, files.lagCsv);
writetable(monthTable, files.monthCsv);
write_latex_local(files.tex, dynTable, betaTable, lagTable, monthTable, opts);

out = struct();
out.files = files;
out.dynamicsTable = dynTable;
out.betaTable = betaTable;
out.lagWeightTable = lagTable;
out.monthlyVarianceTable = monthTable;

fprintf('Saved SS parameter appendix tables to:\n  %s\n', opts.outputDir);
fprintf('Main LaTeX file:\n  %s\n', files.tex);
end

%% ========================================================================
function opts = default_options_local()
opts = struct();
opts.outputDir = fullfile(pwd, 'results_ss_parameter_appendix');
opts.captionPrefix = '';
opts.decimals = 4;
opts.smallDecimals = 3;
opts.includeLagWeightTable = true;
end

%% ========================================================================
function T = make_dynamics_table_local(detail, datasets)

rows = {};
for i = 1:numel(datasets)
    ds = datasets{i};
    est = detail.(ds).pointEstimate;
    p = est.paramsHat;
    setup = est.setup;
    q = p.q;
    dt = setup.dt;
    accelMonthlySd = 100 * sqrt(max(q * dt, 0));
    residBaseSd = sqrt(max(p.sig2u, 0));
    logLik = est.logLik;
    nObs = sum(isfinite(est.filter.vStd));
    nParam = numel(est.thetaHat);
    bic = -2 * logLik + nParam * log(max(nObs, 1));
    rows(end+1, :) = {string(ds), q, sqrt(max(q, 0)), accelMonthlySd, ...
        p.rho, p.ma, p.sig2u, residBaseSd, logLik, nParam, nObs, bic}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', { ...
    'Dataset','q','sqrt_q','MonthlyAccelerationShockSD_decade2', ...
    'ARMA_AR1','ARMA_MA1','BaseResidualVariance','BaseResidualSD', ...
    'LogLik','NumParameters','NObs','BIC'});
end

%% ========================================================================
function T = make_beta_basis_table_local(detail, datasets)

rows = {};
firstDs = datasets{1};
names = string(detail.(firstDs).sample.WNames(:));
for k = 1:numel(names)
    row = cell(1, 2 + numel(datasets));
    row{1} = names(k);
    row{2} = pretty_beta_name_local(names(k));
    for i = 1:numel(datasets)
        ds = datasets{i};
        beta = detail.(ds).pointEstimate.paramsHat.beta(:);
        row{2+i} = beta(k);
    end
    rows(end+1, :) = row; %#ok<AGROW>
end

varNames = [{'Parameter','Label'}, matlab.lang.makeValidName(datasets(:)')];
T = cell2table(rows, 'VariableNames', varNames);
end

%% ========================================================================
function T = make_lag_weight_table_local(detail, datasets, modelOpts)

covNames = string(modelOpts.covariates.names(:));
maxLags = modelOpts.covariates.maxLags(:);
degrees = modelOpts.covariates.almonDegrees(:);

rows = {};
for c = 1:numel(covNames)
    L = maxLags(c);
    deg = degrees(c);
    basisNames = strings(deg + 1, 1);
    for p = 0:deg
        basisNames(p+1) = sprintf('%s_polydeg%d', covNames(c), p);
    end

    for ell = 0:L
        row = cell(1, 3 + numel(datasets));
        row{1} = covNames(c);
        row{2} = ell;
        row{3} = sprintf('%s lag %d', covNames(c), ell);
        powers = (ell / max(L, 1)) .^ (0:deg);
        for i = 1:numel(datasets)
            ds = datasets{i};
            names = string(detail.(ds).sample.WNames(:));
            beta = detail.(ds).pointEstimate.paramsHat.beta(:);
            b = nan(deg + 1, 1);
            for p = 1:(deg + 1)
                idx = find(names == basisNames(p), 1);
                if ~isempty(idx)
                    b(p) = beta(idx);
                end
            end
            row{3+i} = sum(powers(:) .* b(:), 'omitnan');
        end
        rows(end+1, :) = row; %#ok<AGROW>
    end
end

varNames = [{'Covariate','LagMonths','Label'}, matlab.lang.makeValidName(datasets(:)')];
T = cell2table(rows, 'VariableNames', varNames);
end

%% ========================================================================
function T = make_monthly_variance_table_local(detail, datasets)

monthNames = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]';
rows = {};
for m = 1:12
    row = cell(1, 2 + 3 * numel(datasets));
    row{1} = m;
    row{2} = monthNames(m);
    k = 3;
    for i = 1:numel(datasets)
        ds = datasets{i};
        p = detail.(ds).pointEstimate.paramsHat;
        scale = exp(p.gamma(m));
        variance = p.sig2u * scale;
        sd = sqrt(max(p.sig2u * scale, 0));
        row{k} = scale;
        row{k+1} = variance;
        row{k+2} = sd;
        k = k + 3;
    end
    rows(end+1, :) = row; %#ok<AGROW>
end

varNames = {'Month','MonthName'};
for i = 1:numel(datasets)
    safe = matlab.lang.makeValidName(datasets{i});
    varNames{end+1} = [safe '_VarianceScale']; %#ok<AGROW>
    varNames{end+1} = [safe '_ResidualVariance']; %#ok<AGROW>
    varNames{end+1} = [safe '_ResidualSD']; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', varNames);
end

%% ========================================================================
function write_latex_local(texFile, dynT, betaT, lagT, monthT, opts)

fid = fopen(texFile, 'w');
if fid < 0
    error('Could not open LaTeX output file: %s', texFile);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, '%% Auto-generated by write_ss_parameter_appendix_tables.m\n');
fprintf(fid, '%% Requires booktabs; lag-weight table uses longtable.\n\n');

write_dynamics_latex_local(fid, dynT);
fprintf(fid, '\n\n');
write_beta_latex_local(fid, betaT, opts);
fprintf(fid, '\n\n');
write_month_scale_latex_local(fid, monthT, opts);
fprintf(fid, '\n\n');
write_month_variance_latex_local(fid, monthT);
if opts.includeLagWeightTable
    fprintf(fid, '\n\n');
    write_lag_latex_local(fid, lagT, opts);
end
end

%% ========================================================================
function write_dynamics_latex_local(fid, T)

fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n\\scriptsize\n\\setlength{\\tabcolsep}{4pt}\n');
fprintf(fid, '\\caption{\\textbf{State-space dynamic and variance parameter estimates.} ');
fprintf(fid, 'The parameter $q$ scales the continuous-time integrated-Brownian process covariance for acceleration; ');
fprintf(fid, '$100\\sqrt{q\\Delta t}$ reports the monthly acceleration shock standard deviation in $^\\circ$C decade$^{-2}$. ');
fprintf(fid, '$\\sigma_u^2$ is the baseline residual innovation variance before monthly variance scaling.}\\label{tab:ss_param_dynamic}\n');
fprintf(fid, '\\begin{tabular}{lcccccc}\n\\toprule\n');
fprintf(fid, 'Dataset & $q$ & $100\\sqrt{q\\Delta t}$ & AR(1) & MA(1) & $\\sigma_u^2$ & $\\sigma_u$ \\\\\n');
fprintf(fid, '\\midrule\n');
for i = 1:height(T)
    fprintf(fid, '%s & %s & %.4f & %.3f & %.3f & %s & %.4f \\\\\n', ...
        tex_escape_local(T.Dataset(i)), sci_local(T.q(i)), ...
        T.MonthlyAccelerationShockSD_decade2(i), T.ARMA_AR1(i), T.ARMA_MA1(i), ...
        sci_local(T.BaseResidualVariance(i)), T.BaseResidualSD(i));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

%% ========================================================================
function write_beta_latex_local(fid, T, opts)

dsNames = string(T.Properties.VariableNames(3:end));
fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n\\scriptsize\n\\setlength{\\tabcolsep}{4pt}\n');
fprintf(fid, '\\caption{\\textbf{Covariate coefficient estimates.} Coefficients are for the Almon distributed-lag basis used in the state-space observation equation. Covariates are standardized before constructing lag bases.}\\label{tab:ss_param_beta}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('c', 1, numel(dsNames)));
fprintf(fid, 'Coefficient');
for j = 1:numel(dsNames)
    fprintf(fid, ' & %s', tex_escape_local(dsNames(j)));
end
fprintf(fid, ' \\\\\n\\midrule\n');
for i = 1:height(T)
    fprintf(fid, '%s', tex_escape_local(T.Label(i)));
    for j = 1:numel(dsNames)
        fprintf(fid, ' & %s', fmt_num_local(T.(dsNames(j))(i), opts.decimals));
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

%% ========================================================================
function write_month_scale_latex_local(fid, T, opts)

ds = extract_dataset_prefixes_local(T.Properties.VariableNames(3:end));
fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n\\tiny\n\\setlength{\\tabcolsep}{3pt}\n');
fprintf(fid, '\\caption{\\textbf{Monthly residual variance scales.} Entries report $\\exp(\\gamma_m)$, the multiplicative variance scale for month $m$ relative to the baseline residual innovation variance.}\\label{tab:ss_param_monthly_variance}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('c', 1, numel(ds)));
fprintf(fid, 'Month');
for j = 1:numel(ds)
    fprintf(fid, ' & %s', tex_escape_local(ds(j)));
end
fprintf(fid, ' \\\\\n\\midrule\n');
for i = 1:height(T)
    fprintf(fid, '%s', tex_escape_local(T.MonthName(i)));
    for j = 1:numel(ds)
        col = matlab.lang.makeValidName(ds(j) + "_VarianceScale");
        fprintf(fid, ' & %s', fmt_num_local(T.(col)(i), opts.smallDecimals));
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

%% ========================================================================
function write_month_variance_latex_local(fid, T)

ds = extract_dataset_prefixes_local(T.Properties.VariableNames(3:end));
fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n\\tiny\n\\setlength{\\tabcolsep}{3pt}\n');
fprintf(fid, '\\caption{\\textbf{Monthly residual innovation variances.} Entries report $\\sigma_u^2\\exp(\\gamma_m)$ for each month.}\\label{tab:ss_param_monthly_residual_variance}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('c', 1, numel(ds)));
fprintf(fid, 'Month');
for j = 1:numel(ds)
    fprintf(fid, ' & %s', tex_escape_local(ds(j)));
end
fprintf(fid, ' \\\\\n\\midrule\n');
for i = 1:height(T)
    fprintf(fid, '%s', tex_escape_local(T.MonthName(i)));
    for j = 1:numel(ds)
        col = matlab.lang.makeValidName(ds(j) + "_ResidualVariance");
        fprintf(fid, ' & %s', sci_or_fixed_local(T.(col)(i)));
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

%% ========================================================================
function write_lag_latex_local(fid, T, opts)

dsNames = string(T.Properties.VariableNames(4:end));
covariates = unique(T.Covariate, 'stable');
fprintf(fid, '\\begin{longtable}{l%s}\n', repmat('c', 1, numel(dsNames)));
fprintf(fid, '\\caption{\\textbf{Implied distributed-lag covariate weights.} Entries are the implied coefficient at each monthly lag from the estimated Almon basis coefficients.}\\label{tab:ss_param_lag_weights}\\\\\n');
fprintf(fid, '\\toprule\nCovariate/lag');
for j = 1:numel(dsNames)
    fprintf(fid, ' & %s', tex_escape_local(dsNames(j)));
end
fprintf(fid, ' \\\\\n\\midrule\n\\endfirsthead\n');
fprintf(fid, '\\toprule\nCovariate/lag');
for j = 1:numel(dsNames)
    fprintf(fid, ' & %s', tex_escape_local(dsNames(j)));
end
fprintf(fid, ' \\\\\n\\midrule\n\\endhead\n');
for c = 1:numel(covariates)
    C = T(T.Covariate == covariates(c), :);
    fprintf(fid, '\\multicolumn{%d}{l}{\\textbf{%s}}\\\\\n', 1 + numel(dsNames), tex_escape_local(covariates(c)));
    for i = 1:height(C)
        fprintf(fid, '%d', C.LagMonths(i));
        for j = 1:numel(dsNames)
            fprintf(fid, ' & %s', fmt_num_local(C.(dsNames(j))(i), opts.decimals));
        end
        fprintf(fid, ' \\\\\n');
    end
    if c < numel(covariates)
        fprintf(fid, '\\addlinespace[2pt]\n');
    end
end
fprintf(fid, '\\bottomrule\n\\end{longtable}\n');
end

%% ========================================================================
function label = pretty_beta_name_local(name)
name = string(name);
parts = split(name, "_polydeg");
if numel(parts) == 2
    label = sprintf('%s Almon degree %s', parts(1), parts(2));
else
    label = name;
end
label = string(label);
end

%% ========================================================================
function out = extract_dataset_prefixes_local(varNames)
out = strings(0, 1);
for i = 1:numel(varNames)
    nm = string(varNames{i});
    suffix = "_VarianceScale";
    if endsWith(nm, suffix)
        out(end+1, 1) = extractBefore(nm, strlength(nm) - strlength(suffix) + 1); %#ok<AGROW>
    end
end
end

%% ========================================================================
function s = tex_escape_local(x)
x = char(string(x));
x = strrep(x, '_', '\\_');
s = x;
end

%% ========================================================================
function s = fmt_num_local(x, nd)
if ~isfinite(x)
    s = '--';
elseif abs(x) > 0 && abs(x) < 10^(-nd)
    s = sci_local(x);
else
    s = sprintf(['%0.' num2str(nd) 'f'], x);
end
end

%% ========================================================================
function s = sci_local(x)
if ~isfinite(x)
    s = '--';
elseif x == 0
    s = '0';
else
    expo = floor(log10(abs(x)));
    mant = x / 10^expo;
    s = sprintf('$%.3f\\times 10^{%d}$', mant, expo);
end
end

%% ========================================================================
function s = sci_or_fixed_local(x)
if ~isfinite(x)
    s = '--';
elseif abs(x) < 0.001
    s = sci_local(x);
else
    s = sprintf('%.4f', x);
end
end

%% ========================================================================
function out = merge_structs_local(a, b)
out = a;
if isempty(b)
    return;
end
fields = fieldnames(b);
for i = 1:numel(fields)
    f = fields{i};
    if isstruct(b.(f)) && isfield(out, f) && isstruct(out.(f))
        out.(f) = merge_structs_local(out.(f), b.(f));
    else
        out.(f) = b.(f);
    end
end
end
