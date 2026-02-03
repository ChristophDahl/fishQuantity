%% runFishQuantityPipeline.m

clc; close all;

fprintf('\n=== runFishQuantityPipeline: START ===\n');

%% ------------------------------------------------------------------------
% 1) SMALL fish: descriptive profiles + overview/export + laterality (LME)
% -------------------------------------------------------------------------
fprintf('\n--- SMALL fish: stage 1 profiles ---\n');
stage1_smallFish_profiles_AB

fprintf('\n--- SMALL fish: stage 1 overview/export ---\n');
stage1_smallFish_laterality_export

fprintf('\n--- SMALL fish: stage 2 laterality (LME) ---\n');
stage2_smallFish_laterality_LME

% OTS vs ANS (expects all_trials to exist after stage2 script)
fprintf('\n--- SMALL fish: OTS vs ANS ---\n');
resS = fitOTSvsANSModels(all_trials); %#ok<NASGU>

% Save the current figure (same behaviour as your earlier main)
exportFigureDual(figurefolder, 'SmallFishComparisons');

%% ------------------------------------------------------------------------
% 2) LARGE fish: descriptive profiles + laterality (LME)
% -------------------------------------------------------------------------
fprintf('\n--- LARGE fish: stage 1 profiles ---\n');
stage1_largeFish_profiles_AB

fprintf('\n--- LARGE fish: stage 2 laterality (LME) ---\n');
stage2_largeFish_laterality_LME

fprintf('\n--- LARGE fish: OTS vs ANS ---\n');
resL = fitOTSvsANSModels(all_trials, 'IncludeSizes', 'l'); %#ok<NASGU>

exportFigureDual(figurefolder, 'LargeFishComparisons');

%% ------------------------------------------------------------------------
% 3) Single contrasts vs 0
% -------------------------------------------------------------------------
fprintf('\n--- Single contrasts vs 0 (small) ---\n');
stage3_singleContrast_smallFish

fprintf('\n--- Single contrasts vs 0 (large) ---\n');
stage3_singleContrast_largeFish

exportFigureDual(figurefolder, 'SingleComparisons');

%% ------------------------------------------------------------------------
% 4) Stimulus controls: image-derived metrics + condition models
% -------------------------------------------------------------------------
fprintf('\n--- Stimulus controls: page metrics ---\n');
inDir = "I:\fishQuantity\stimuli\pages_halves\cropped\vtrim_manual";

% SMALL fish pages
T = analyze_fish_pages_range(inDir, 2, 16, ...
    'Sides',["L","R"], ...
    'Polarity',"dark",'Sensitivity',0.40,'Sigma',0.7, ...
    'OpenR',1,'CloseR',0, ...
    'MinAreaAbs',30,'BandLow',0.40,'BandHigh',2.50, ...
    'BigFrac',0.6,'DilateFrac',0.6,'MaxLinkFactor',1.2, ...
    'FillNNWithZero',true);
disp(T(:, {'Fn','Nobj','Density','Coverage','CentVar','CentStdX','CentStdY','NNmean'}));
T_small = T; %#ok<NASGU>
clear T;

% LARGE fish pages
T = analyze_fish_pages_range(inDir, 17, 31, ...
    'Sides',["L","R"], ...
    'Polarity',"dark",'Sensitivity',0.40,'Sigma',0.1, ...
    'OpenR',0,'CloseR',1, ...
    'MinAreaAbs',30,'BandLow',0.40,'BandHigh',2.50, ...
    'BigFrac',0.6,'DilateFrac',0.6,'MaxLinkFactor',1.2, ...
    'FillNNWithZero',true);
disp(T(:, {'Fn','Nobj','Density','Coverage','CentVar','CentStdX','CentStdY','NNmean'}));
T_large = T; %#ok<NASGU>
clear T;

fprintf('\n--- Stimulus controls: summarize small vs large ---\n');
[Gs, Gl] = summarize_small_large_runs(inDir);

save(fullfile(resultsfolder,'stimuli_metrics_small.mat'),'Gs');
save(fullfile(resultsfolder,'stimuli_metrics_large.mat'),'Gl');

fprintf('\n--- Stimulus controls: fit laterality/condition models ---\n');
[M_small, M_large, T_small, T_large] = fit_laterality_condition_models(resultsfolder, ...
    'FileSuffix','_raw', 'LateralityField','LateralityIdx'); %#ok<ASGLU>

%% ------------------------------------------------------------------------
% 5) Example stimulus figure + interaction model
% -------------------------------------------------------------------------
fprintf('\n--- Stimulus example figure ---\n');
imgSmall = 'I:\fishQuantity\stimuli\pages_halves\cropped\vtrim_manual\page_0006_R.png';
imgLarge = 'I:\fishQuantity\stimuli\pages_halves\cropped\vtrim_manual\page_0012_L.png';

outPng = fullfile(figurefolder, 'stimulus_examples.png');
make_stimulus_pair_figure(outPng, imgSmall, imgLarge);
fprintf('Saved: %s\n', outPng);

fprintf('\n--- Interaction model ---\n');
out = interactionFishQuantity;

idx = find(contains(string(out.coef.Name),"log2_ratio") & ...
           contains(string(out.coef.Name),"size") & ...
           contains(string(out.coef.Name),":"));
out.coef(idx,:)

out.compare

fprintf('\n=== runFishQuantityPipeline: DONE ===\n');

%% ------------------------------------------------------------------------
% Local helper: export current figure
% -------------------------------------------------------------------------
function exportFigureDual(figDir, baseName)
% Save current figure as vector PDF + 600 dpi PNG.
    exportgraphics(gcf, fullfile(figDir, baseName + ".pdf"), ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'none');
    exportgraphics(gcf, fullfile(figDir, baseName + ".png"), ...
        'Resolution', 600);
end
