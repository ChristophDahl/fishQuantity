function [M_small, M_large, T_small, T_large] = fit_laterality_condition_models(resultsfolder, varargin)
%FIT_LATERALITY_CONDITION_MODELS  Merge laterality summaries with stimulus metrics and fit LMEs.
%
% This function:
%   1) loads laterality condition summaries (small/large) from <resultsfolder>
%   2) loads stimulus metrics (Gs/Gl) from stimuli_metrics_small.mat / _large.mat
%   3) constructs Δ and Σ covariates (Coverage; and NNmean if available)
%   4) fits separate linear models for SMALL and LARGE conditions
%
% Output:
%   M_small, M_large : fitted linear models (fitlm)
%   T_small, T_large : merged tables used for fitting (saved to CSV as side-effect)
%


ip = inputParser;
addParameter(ip,'FileSuffix','_raw');
addParameter(ip,'LateralityField','LateralityIdx');
parse(ip,varargin{:});
suffix  = ip.Results.FileSuffix;
latField= ip.Results.LateralityField;

% ---------- 1) Load condition-level laterality summaries ----------
Bs = readtable(fullfile(resultsfolder, ['laterality_conditions_small' suffix '.csv']));
Bl = readtable(fullfile(resultsfolder, ['laterality_conditions_large' suffix '.csv']));

% Parse 'Contrast' like '1 vs 16' -> N1,N2
[Bs.N1, Bs.N2] = parse_contrast(Bs.Contrast);
[Bl.N1, Bl.N2] = parse_contrast(Bl.Contrast);

% Ensure Log2Ratio exists
if ~ismember('Log2Ratio', Bs.Properties.VariableNames), Bs.Log2Ratio = log2(Bs.Ratio); end
if ~ismember('Log2Ratio', Bl.Properties.VariableNames), Bl.Log2Ratio = log2(Bl.Ratio); end

% Pick laterality DV
assert(ismember(latField, Bs.Properties.VariableNames), 'Laterality field not found in small: %s', latField);
assert(ismember(latField, Bl.Properties.VariableNames), 'Laterality field not found in large: %s', latField);
Bs.Laterality = Bs.(latField);
Bl.Laterality = Bl.(latField);

% ---------- 2) Load stimulus metrics & aggregate by (Nobj,Side if available) ----------
S_s = load(fullfile(resultsfolder,'stimuli_metrics_small.mat'),'Gs'); S_s = S_s.Gs;
S_l = load(fullfile(resultsfolder,'stimuli_metrics_large.mat'),'Gl'); S_l = S_l.Gl;

Smean_s = side_means(S_s);
Smean_l = side_means(S_l);

% ---------- 3) Attach Δ/Σ coverage (and ΔNN if available) ----------
T_small = attach_covariates(Bs, Smean_s);
T_large = attach_covariates(Bl, Smean_l);

% ---------- 4) Fit per-size linear models ----------
form = 'Laterality ~ Log2Ratio + DeltaCoverage + SumCoverage';
form_s = form; form_l = form;
if ismember('DeltaNN', T_small.Properties.VariableNames), form_s = [form_s ' + DeltaNN']; end
if ismember('DeltaNN', T_large.Properties.VariableNames), form_l = [form_l ' + DeltaNN']; end

M_small = fitlm(T_small, form_s);
M_large = fitlm(T_large, form_l);

fprintf('\nSMALL (DV=%s, %s):\n', latField, erase(suffix,'_'));  disp(M_small.Coefficients(:,{'Estimate','SE','pValue'}));
fprintf('\nLARGE (DV=%s, %s):\n', latField, erase(suffix,'_'));  disp(M_large.Coefficients(:,{'Estimate','SE','pValue'}));

% Optional: save merged tables
writetable(T_small, fullfile(resultsfolder, ['merged_condition_table_small' suffix '.csv']));
writetable(T_large, fullfile(resultsfolder, ['merged_condition_table_large' suffix '.csv']));
end

% ================= helpers =================

function [n1, n2] = parse_contrast(C)
% C like "1 vs 16"
pat = '^\s*(\d+)\s*vs\s*(\d+)\s*$';
tokens = regexp(string(C), pat, 'tokens', 'once');
n1 = nan(numel(tokens),1); n2 = nan(numel(tokens),1);
for i=1:numel(tokens)
    if ~isempty(tokens{i})
        n1(i) = str2double(tokens{i}{1});
        n2(i) = str2double(tokens{i}{2});
    end
end
end

function Smean = side_means(S)
% Robustly compute means by Nobj and (if available) Side.
% 1) Identify/derive Side. 2) pick coverage-like column. 3) groupsummary.

% --- (a) ensure Side exists if we can infer it ---
if ~ismember('Side', S.Properties.VariableNames)
    % look for a filename-like column to infer L/R
    fncol = intersect(S.Properties.VariableNames, {'Fn','File','Filename','FileName','Path','Image','ImagePath'}, 'stable');
    if ~isempty(fncol)
        fstr = string(S.(fncol{1}));
        tok  = regexp(fstr, '_([LR])\.[A-Za-z0-9]+$','tokens','once');
        Side = strings(numel(fstr),1);
        for i=1:numel(fstr)
            if ~isempty(tok{i}), Side(i) = tok{i}{1}; else, Side(i) = ""; end
        end
        % only set if we actually found usable tags
        if any(Side~="")
            S.Side = categorical(Side);
        end
    end
end
if ~ismember('Side', S.Properties.VariableNames)
    % still no side info: create a missing Side column (we will fallback later)
    S.Side = categorical(repmat(missing,height(S),1));
end

% --- (b) pick a coverage-like column (or accept pre-aggregated mean_Coverage) ---
covCands = {'mean_Coverage','Coverage','coverage','Cover','AreaFrac','AreaFraction'};
covCol = '';
for k=1:numel(covCands)
    if ismember(covCands{k}, S.Properties.VariableNames)
        covCol = covCands{k}; break
    end
end
assert(~isempty(covCol), 'No coverage-like column found in stimuli metrics.');

% (optional) NN
hasNN = ismember('NNmean', S.Properties.VariableNames);

% --- (c) If already pre-aggregated (mean_Coverage), just pass through (ensure naming) ---
if strcmp(covCol,'mean_Coverage')
    Smean = S;
    if hasNN && ~ismember('mean_NNmean', Smean.Properties.VariableNames) && ismember('NNmean',Smean.Properties.VariableNames)
        Smean.mean_NNmean = Smean.NNmean;
    end
    return
end

% --- (d) Aggregate by available keys ---
grouping = {'Nobj'};
if any(~ismissing(S.Side)), grouping = {'Nobj','Side'}; end

aggVars = {covCol};
if hasNN, aggVars{end+1} = 'NNmean'; end

Smean = groupsummary(S, grouping, 'mean', aggVars);

% unify names
Smean.Properties.VariableNames = strrep(Smean.Properties.VariableNames, ['mean_' covCol], 'mean_Coverage');
end

function T = attach_covariates(B, Smean)
% Build Delta/Sum coverage (and DeltaNN if present) for each contrast row in B.

% helpers
hasSide = ismember('Side', Smean.Properties.VariableNames);
hasNNm  = ismember('mean_NNmean', Smean.Properties.VariableNames);

getMean = @(N,side,field) local_get_mean(Smean,N,side,field,hasSide);
fallback= @(N,field) mean(Smean.(sprintf('mean_%s',field))(Smean.Nobj==N),'omitnan');

DeltaCoverage = nan(height(B),1);
SumCoverage   = nan(height(B),1);
DeltaNN       = nan(height(B),1);

for i=1:height(B)
    nL = B.N1(i); nR = B.N2(i);
    cL = getMean(nL,'L','Coverage'); if isnan(cL), cL = fallback(nL,'Coverage'); end
    cR = getMean(nR,'R','Coverage'); if isnan(cR), cR = fallback(nR,'Coverage'); end
    DeltaCoverage(i) = cL - cR;
    SumCoverage(i)   = cL + cR;

    if hasNNm
        nnL = getMean(nL,'L','NNmean'); if isnan(nnL), nnL = fallback(nL,'NNmean'); end
        nnR = getMean(nR,'R','NNmean'); if isnan(nnR), nnR = fallback(nR,'NNmean'); end
        DeltaNN(i) = nnL - nnR;
    end
end

T = B(:, {'Contrast','Size','N1','N2','Ratio','Log2Ratio','Laterality'});
T.DeltaCoverage = DeltaCoverage;
T.SumCoverage   = SumCoverage;
if any(~isnan(DeltaNN)), T.DeltaNN = DeltaNN; end
end

function m = local_get_mean(Smean,N,side,field,hasSide)
vname = sprintf('mean_%s',field);
if hasSide
    idx = (Smean.Nobj==N) & (Smean.Side==categorical({side}));
else
    idx = (Smean.Nobj==N);
end
vals = Smean.(vname)(idx);
if isempty(vals), m = NaN; else, m = vals(1); end
end