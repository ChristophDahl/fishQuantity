function out = interactionFishQuantity(resultsfolder)
%INTERACTIONFISHQUANTITY  Size × log2_ratio interaction model for laterality index.
%
%   out = interactionFishQuantity(resultsfolder)
%
% Purpose
%   Load small/large trial tables, apply strict QC (drop NaN/Inf/undefined),
%   and fit an LME with a Size × log2_ratio interaction on lat_idx.
%
% Required input files (in resultsfolder)
%   - Tsmall_trials.mat (table variable: Tsmall)
%   - Tlarge_trials.mat (table variable: Tlarge)
%
% Outputs (written to resultsfolder)
%   - Interaction_LME_summary.txt
%   - Interaction_LME_results.mat
%   - Interaction_badRows.csv (if any rows are dropped)
%
% Example
%   out = interactionFishQuantity('I:\fishQuantity\results');

    if nargin < 1 || isempty(resultsfolder)
        resultsfolder = 'I:\fishQuantity\results';
    end

    smallFile = fullfile(resultsfolder, 'Tsmall_trials.mat');
    largeFile = fullfile(resultsfolder, 'Tlarge_trials.mat');

    outMat  = fullfile(resultsfolder, 'Interaction_LME_results.mat');
    outTxt  = fullfile(resultsfolder, 'Interaction_LME_summary.txt');
    outBad  = fullfile(resultsfolder, 'Interaction_badRows.csv');

    % If you want to INCLUDE 0-vs-n contrasts inside the log2-ratio analysis,
    % set includeZeroContrasts=true and choose a pseudo-count (e.g., 0.5 or 1).
    % Otherwise, those trials will be dropped (recommended).
    includeZeroContrasts = false;
    pseudo = 0.5;

    %% ---- LOAD TABLES ----
    assert(isfile(smallFile), 'Missing file: %s', smallFile);
    assert(isfile(largeFile), 'Missing file: %s', largeFile);

    Tsmall = loadTableFromMat(smallFile, "Tsmall");
    Tlarge = loadTableFromMat(largeFile, "Tlarge");

    fprintf('Loaded: Tsmall=%d rows, Tlarge=%d rows\n\n', height(Tsmall), height(Tlarge));

    %% ---- STANDARDIZE / VALIDATE ----
    Tsmall = standardizeTrialTable(Tsmall, "small");
    Tlarge = standardizeTrialTable(Tlarge, "large");

    Tall = [Tsmall; Tlarge];

    Tall.subject = categorical(Tall.subject);
    Tall.size    = categorical(Tall.size);
    Tall.size    = reordercats(Tall.size, {'small','large'});  % ensure "small" reference

    %% ---- (RE)COMPUTE log2_ratio DEFENSIVELY ----
    if includeZeroContrasts
        Tall.log2_ratio = log2((Tall.int2 + pseudo) ./ (Tall.int1 + pseudo));
    else
        Tall.log2_ratio = NaN(height(Tall), 1);
        validRatio = isfinite(Tall.int1) & isfinite(Tall.int2) & (Tall.int1 > 0) & (Tall.int2 > 0);
        Tall.log2_ratio(validRatio) = log2(Tall.int2(validRatio) ./ Tall.int1(validRatio));
    end

    %% ---- QC + HARD CLEAN ----
    n0 = height(Tall);

    nan_lr  = isnan(Tall.log2_ratio);
    inf_lr  = isinf(Tall.log2_ratio);
    nan_li  = isnan(Tall.lat_idx);
    inf_li  = isinf(Tall.lat_idx);
    undef_s = isundefined(Tall.subject);
    undef_z = isundefined(Tall.size);

    nonpos  = (Tall.int1 <= 0) | (Tall.int2 <= 0) | ~isfinite(Tall.int1) | ~isfinite(Tall.int2);

    fprintf('QC BEFORE DROP\n');
    fprintf('  NaN log2_ratio: %d\n', sum(nan_lr));
    fprintf('  Inf log2_ratio: %d\n', sum(inf_lr));
    fprintf('  NaN lat_idx:    %d\n', sum(nan_li));
    fprintf('  Inf lat_idx:    %d\n', sum(inf_li));
    fprintf('  undef subject:  %d\n', sum(undef_s));
    fprintf('  undef size:     %d\n', sum(undef_z));
    fprintf('  nonpositive/invalid int1/int2: %d\n', sum(nonpos));

    bad = nan_lr | inf_lr | nan_li | inf_li | undef_s | undef_z;

    TallBad = table();
    if any(bad)
        TallBad = Tall(bad,:);
        writetable(TallBad, outBad);
        fprintf('Dropping %d/%d bad rows. Saved dropped rows to: %s\n', sum(bad), n0, outBad);
        disp(TallBad(1:min(10,height(TallBad)), {'subject','size','cond_key','int1','int2','log2_ratio','lat_idx'}));
        Tall(bad,:) = [];
    else
        fprintf('No bad rows detected.\n');
    end

    fprintf('\nQC AFTER DROP: Tall=%d rows (small=%d, large=%d), subjects=%d\n\n', ...
        height(Tall), sum(Tall.size=='small'), sum(Tall.size=='large'), numel(categories(Tall.subject)));

    if sum(Tall.size=='small')==0 || sum(Tall.size=='large')==0
        error('After cleaning, one size level disappeared. Cannot fit interaction model.');
    end

    %% ---- FIT MODELS (ML for LRT) ----
    % Try random slope first; fallback to random intercept.
    form_int_rs  = 'lat_idx ~ log2_ratio*size + (1 + log2_ratio|subject)';
    form_main_rs = 'lat_idx ~ log2_ratio + size + (1 + log2_ratio|subject)';

    form_int_ri  = 'lat_idx ~ log2_ratio*size + (1|subject)';
    form_main_ri = 'lat_idx ~ log2_ratio + size + (1|subject)';

    useRandomSlope = true;

    try
        mMain = fitlme(Tall, form_main_rs, 'FitMethod','ML');
        mInt  = fitlme(Tall, form_int_rs,  'FitMethod','ML');
    catch ME
        fprintf('Random-slope model failed (%s). Falling back to random-intercept only.\n\n', ME.message);
        useRandomSlope = false;
        mMain = fitlme(Tall, form_main_ri, 'FitMethod','ML');
        mInt  = fitlme(Tall, form_int_ri,  'FitMethod','ML');
    end

    cmp  = compare(mMain, mInt);     % LRT for interaction
    coef = mInt.Coefficients;

    %% ---- EXTRACT FIXED EFFECTS + DERIVED SLOPES ----
    names = string(coef.Name);

    idx_log2 = find(names == "log2_ratio", 1);

    % Robust interaction detection:
    % Works for "size_large:log2_ratio" and "log2_ratio:size_large"
    idx_int = find(contains(names, "log2_ratio") & contains(names, "size") & contains(names, ":"), 1);

    if isempty(idx_log2) || isempty(idx_int)
        error('Could not find log2_ratio or interaction term in fixed effects. Names were: %s', strjoin(names, ', '));
    end

    beta = coef.Estimate;
    V    = mInt.CoefficientCovariance;

    slope_small = beta(idx_log2);
    se_small    = sqrt(V(idx_log2, idx_log2));
    ci_small    = [slope_small - 1.96*se_small, slope_small + 1.96*se_small];

    slope_large = beta(idx_log2) + beta(idx_int);
    se_large    = sqrt(V(idx_log2, idx_log2) + V(idx_int, idx_int) + 2*V(idx_log2, idx_int));
    ci_large    = [slope_large - 1.96*se_large, slope_large + 1.96*se_large];

    %% ---- WRITE SUMMARY ----
    fid = fopen(outTxt, 'w');
    assert(fid>0, 'Could not open output text file: %s', outTxt);

    fprintf(fid, '=== Size x log2_ratio interaction LME (trial-level lat_idx) ===\n');
    fprintf(fid, 'resultsfolder: %s\n', resultsfolder);
    fprintf(fid, 'includeZeroContrasts: %d (pseudo=%.3f)\n\n', includeZeroContrasts, pseudo);

    fprintf(fid, 'Data after cleaning: %d trials; subjects=%d; small=%d; large=%d\n', ...
        height(Tall), numel(categories(Tall.subject)), sum(Tall.size=='small'), sum(Tall.size=='large'));
    fprintf(fid, 'Random effects: %s\n\n', ternary(useRandomSlope, '(1 + log2_ratio | subject)', '(1 | subject)'));

    fprintf(fid, '--- Model comparison (no interaction vs interaction; ML) ---\n');
    fprintf(fid, '%s\n\n', evalc('disp(cmp)'));

    fprintf(fid, '--- Fixed effects (interaction model; ML) ---\n');
    fprintf(fid, '%s\n\n', evalc('disp(coef)'));

    fprintf(fid, 'Interaction term used: %s\n', coef.Name{idx_int});
    fprintf(fid, 'Interaction beta=%.6f, SE=%.6f, t=%.3f, df=%.1f, p=%.6g\n\n', ...
        coef.Estimate(idx_int), coef.SE(idx_int), coef.tStat(idx_int), coef.DF(idx_int), coef.pValue(idx_int));

    fprintf(fid, '--- Derived slopes (dLI/dlog2_ratio) ---\n');
    fprintf(fid, 'Small: %.6f (SE %.6f), 95%% CI [%.6f, %.6f]\n', slope_small, se_small, ci_small(1), ci_small(2));
    fprintf(fid, 'Large: %.6f (SE %.6f), 95%% CI [%.6f, %.6f]\n', slope_large, se_large, ci_large(1), ci_large(2));

    fclose(fid);

    %% ---- SAVE RESULTS ----
    save(outMat, ...
        'Tall','TallBad','mMain','mInt','cmp','coef', ...
        'slope_small','se_small','ci_small','slope_large','se_large','ci_large', ...
        'useRandomSlope','includeZeroContrasts','pseudo');

    fprintf('Saved summary: %s\n', outTxt);
    fprintf('Saved results: %s\n', outMat);
    if any(bad)
        fprintf('Saved dropped rows: %s\n', outBad);
    end

    %% ---- RETURN STRUCT ----
    out = struct();
    out.Tall = Tall;
    out.TallBad = TallBad;
    out.mMain = mMain;
    out.mInt  = mInt;
    out.compare = cmp;
    out.coef = coef;
    out.slope_small = slope_small;
    out.se_small = se_small;
    out.ci_small = ci_small;
    out.slope_large = slope_large;
    out.se_large = se_large;
    out.ci_large = ci_large;
    out.useRandomSlope = useRandomSlope;
    out.includeZeroContrasts = includeZeroContrasts;
    out.pseudo = pseudo;
    out.outTxt = outTxt;
    out.outMat = outMat;
    out.outBad = outBad;

end

%% ========================================================================
% Local helper functions
% ========================================================================

function T = loadTableFromMat(matPath, preferredVar)
    S = load(matPath);

    if isfield(S, preferredVar)
        T = S.(preferredVar);
        if ~istable(T)
            error('%s exists in %s but is not a table.', preferredVar, matPath);
        end
        return;
    end

    % Fallback: first table in file
    fn = fieldnames(S);
    for i = 1:numel(fn)
        if istable(S.(fn{i}))
            T = S.(fn{i});
            warning('Variable "%s" not found in %s. Using table "%s" instead.', preferredVar, matPath, fn{i});
            return;
        end
    end

    error('No table found in %s.', matPath);
end

function T = standardizeTrialTable(T, defaultSize)
    if ~istable(T)
        error('Input is not a table.');
    end

    % subject
    if ~ismember('subject', T.Properties.VariableNames)
        error('Table is missing required variable: subject');
    end

    % size
    if ~ismember('size', T.Properties.VariableNames)
        if ismember('cond_key', T.Properties.VariableNames)
            s = string(T.cond_key);
            inferred = repmat(defaultSize, height(T), 1);
            inferred(endsWith(s,'_s')) = "small";
            inferred(endsWith(s,'_l')) = "large";
            T.size = categorical(inferred);
        else
            T.size = categorical(repmat({defaultSize}, height(T), 1));
        end
    end

    % cond_key (useful but not strictly required)
    if ~ismember('cond_key', T.Properties.VariableNames)
        T.cond_key = categorical(repmat({""}, height(T), 1));
    else
        % keep as categorical for compactness
        T.cond_key = categorical(T.cond_key);
    end

    % int1, int2
    if ~ismember('int1', T.Properties.VariableNames) || ~ismember('int2', T.Properties.VariableNames)
        error('Table must contain int1 and int2.');
    end
    T.int1 = double(T.int1);
    T.int2 = double(T.int2);

    % lat_idx
    if ~ismember('lat_idx', T.Properties.VariableNames)
        error('Table is missing required variable: lat_idx');
    end
    T.lat_idx = double(T.lat_idx);

    % log2_ratio exists (will be recomputed anyway)
    if ~ismember('log2_ratio', T.Properties.VariableNames)
        T.log2_ratio = NaN(height(T), 1);
    else
        T.log2_ratio = double(T.log2_ratio);
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end