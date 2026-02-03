%% stage3_singleContrast_largeFish.m
% Single-contrast analysis: 0 vs N compartment occupancy profiles (large stimuli)
%
% Expected inputs/files:
% - resultsfolder contains: roi_rect.mat and QF*_all_trials.mat
% - extract_profiles (local function below) expects each subj_data trial to
%   include fields used in the original pipeline (e.g., trials{k}.prop_per_bin)
%
% Output:
% - Creates figures and workspace variables exactly as the original script.


clear all
clc


%% ------------------------------------------------------------------------
% 1) Paths, loading, settings
% -------------------------------------------------------------------------
% --- Paths, loading, settings ---
filefolder = 'I:\fishQuantity\files';
programfolder = 'I:\fishQuantity\programs';
resultsfolder = 'I:\fishQuantity\results';
figurefolder = 'I:\fishQuantity\figures';
addpath(genpath('I:\fishQuantity'))
cd(resultsfolder);

load(fullfile(resultsfolder, 'roi_rect.mat'), 'roi_rect');
files = dir(fullfile(resultsfolder, 'QF*_all_trials.mat'));
subject_ids = cell(1, numel(files));
all_subjects = cell(1, numel(files));
for s = 1:numel(files)
    subj_id = regexp(files(s).name, '(QF\d{2})_all_trials.mat', 'tokens', 'once');
    subject_ids{s} = subj_id{1};
    tmp = load(fullfile(resultsfolder, files(s).name));
    all_subjects{s} = tmp.subj_data;  % subj_data was saved as variable
end

roi_width = roi_rect(3);
n_bins = 7; % Number of compartments

%% ------------------------------------------------------------------------
% 2) Define contrasts: only 0 vs N (size-coded)
% -------------------------------------------------------------------------
% ***** MODIFIED SECTION: Only 0 vs N ("s") contrasts *****
group_sizes = [1 2 4 8 16];
desired_order = arrayfun(@(n) sprintf('0_%d', n), group_sizes, 'UniformOutput', false);
labels_large = strcat(desired_order, '_l');
n = length(desired_order);

%% ------------------------------------------------------------------------
% 3) Extract trial-wise profiles and compute mean±SEM (raw + bias-corrected)
% -------------------------------------------------------------------------
% --- Store all results in a struct: results_raw and results_corr ---
% results_raw = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_small, false);
% results_corr = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_small, true);
[results_raw, all_laterality] = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_large, false);
results_corr = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_large, true);

means_raw = results_raw.means;
sems_raw = results_raw.sems;
xlabels_raw = labels_large; % order is preserved

means_corr = results_corr.means;
sems_corr = results_corr.sems;
xlabels_corr = labels_large;

%% ------------------------------------------------------------------------
% 4) Curve fitting (linear + quadratic) on mean profiles
% -------------------------------------------------------------------------
% ---- Curve fitting for each 0 vs N contrast (mean profile, RAW) ----
curvefit_results = struct();
fprintf('\nCurve fits for 0 vs N (large, raw means):\n');
fprintf('N\tSlope\tCurvature\tR2\n');

for i = 1:numel(labels_large)
    mean_profile = means_raw(i, :);  % mean compartment occupancy
    x = (1:n_bins)';
    y = mean_profile(:);
    p1 = polyfit(x, y, 1);    % Linear fit
    p2 = polyfit(x, y, 2);    % Quadratic fit
    y1 = polyval(p1, x);
    y2 = polyval(p2, x);
    R2_1 = 1 - sum((y - y1).^2)/sum((y-mean(y)).^2);
    R2_2 = 1 - sum((y - y2).^2)/sum((y-mean(y)).^2);
    % Save linear fit
    curvefit_results(i).label = labels_large{i};
    curvefit_results(i).N = group_sizes(i);
    curvefit_results(i).Slope = p1(1);
    curvefit_results(i).Curvature = p2(1);
    curvefit_results(i).R2 = R2_2; % You can also report R2_1 if preferred
    fprintf('%2d\t% .5f\t% .5f\t%.4f\n', group_sizes(i), p1(1), p2(1), R2_2);
end

% (Optional) Save as a table:
T_curvefits = struct2table(curvefit_results);
writetable(T_curvefits, fullfile(resultsfolder, 'curvefits_zero_vs_N_small.csv'));

% ---- Curve fitting for each 0 vs N contrast (mean profile, CORRECTED) ----
curvefit_corr_results = struct();
fprintf('\nCurve fits for 0 vs N (small, bias-corrected means):\n');
fprintf('N\tSlope\tCurvature\tR2\n');

for i = 1:numel(labels_large)
    mean_profile = means_corr(i, :);  % bias-corrected mean
    x = (1:n_bins)';
    y = mean_profile(:);
    p1 = polyfit(x, y, 1);
    p2 = polyfit(x, y, 2);
    y1 = polyval(p1, x);
    y2 = polyval(p2, x);
    R2_1 = 1 - sum((y - y1).^2)/sum((y-mean(y)).^2);
    R2_2 = 1 - sum((y - y2).^2)/sum((y-mean(y)).^2);
    curvefit_corr_results(i).label = labels_large{i};
    curvefit_corr_results(i).N = group_sizes(i);
    curvefit_corr_results(i).Slope = p1(1);
    curvefit_corr_results(i).Curvature = p2(1);
    curvefit_corr_results(i).R2 = R2_2;
%     fprintf('%2d\t% .5f\t%  .5f\t%.4f\n', group_sizes(i), p1(1), p2(1), R2_2);
    fprintf('%2d\t% .5f\t% .5f\t%.4f\n', group_sizes(i), p1(1), p2(1), R2_2);

end

T_curvefits_corr = struct2table(curvefit_corr_results);
writetable(T_curvefits_corr, fullfile(resultsfolder, 'curvefits_zero_vs_N_large_corr.csv'));

% save(fullfile(resultsfolder, 'laterality_zero_vs_large.mat'), 'all_laterality', 'labels_large');
save(fullfile(resultsfolder, 'laterality_zero_vs_large.mat'), ...
     'all_laterality', ...
     'labels_large', ...
     'curvefit_results', ...
     'curvefit_corr_results', ...
     'T_curvefits', ...
     'T_curvefits_corr');
 


row_names = { ...
    labels_large ...
    };
n_rows = numel(row_names);
n_cols = max(cellfun(@numel, row_names)); % Only other cols, diagonal-corrected is extra

% figure('Position',[100 100 1400 1100]);
% plot_idx = 1;
for r = 1:n_rows
    n_this_row = numel(row_names{r});
    % --- LEFT-most column: corrected diagonal profile ---
    cond_diag = row_names{r}{1};
    idx_corr = find(strcmp(xlabels_corr, cond_diag));
%     subplot(n_rows, n_cols+1, (r-1)*(n_cols+1) + 1);
%     if ~isempty(idx_corr)
%         meanx = means_corr(idx_corr,:);
%         errs = sems_corr(idx_corr,:);
%         plot(1:n_bins, meanx, '.', 'Markersize', 10, ...
%             'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1]);
%         hold on
%         for k = 1:length(errs)
%             plot([k k], [meanx(k) - errs(k), meanx(k) + errs(k)], '-', 'Color', [0 0 1])
%         end
%         % Fit, blue dashed
%         x = (1:n_bins)';
%         y = meanx(:);
%         p1 = polyfit(x, y, 1); y1 = polyval(p1, x);
%         p2 = polyfit(x, y, 2); y2 = polyval(p2, x);
%         R2_1 = 1 - sum((y - y1).^2)/sum((y-mean(y)).^2);
%         R2_2 = 1 - sum((y - y2).^2)/sum((y-mean(y)).^2);
%         if R2_2 > R2_1, best_fit = y2; else, best_fit = y1; end
%         plot(x, best_fit, 'b--', 'LineWidth', 1.3);
% 
%         title(contrast_title(cond_diag, true), 'FontSize', 10, 'Interpreter','tex');
%         if r == n_rows
%             xlabel('Compartment'); 
%         end
%         ylabel('Proportion');
%         set(gca, ...
%             'Box', 'off', ...
%             'TickDir', 'out', ...
%             'TickLength', [.01 .01], ...
%             'XColor', [0 0 0], ...
%             'YColor', [0 0 0], ...
%             'XTick', 1:n_bins, ...
%             'YTick', 0:.1:.4, ...
%             'LineWidth', .75);
%         axis square
%         axis([0 n_bins+1 0 .4])
%         hold off
%     else
%         axis off
%     end
%     if r == 1
%         text(-2, .5, 'A', 'Fontsize',14, 'FontWeight', 'bold')
%     end
    
    % --- All other columns: raw version ---
    for c = 1:n_this_row
        cond = row_names{r}{c};
        idx_raw = find(strcmp(xlabels_raw, cond));
        subplot(n_rows*2, n_cols + 2, (r-1)*(n_cols) + c + n_cols + 4);
        if ~isempty(idx_raw)
            meanx = means_raw(idx_raw,:);
            errs = sems_raw(idx_raw,:);
            plot(1:n_bins, meanx, '.', 'Markersize', 10, ...
                'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
            hold on
            for k = 1:length(errs)
                plot([k k], [meanx(k) - errs(k), meanx(k) + errs(k)], '-', 'Color', [0 0 0])
            end
            % Red dashed fit
            x = (1:n_bins)';
            y = meanx(:);
            p1 = polyfit(x, y, 1); y1 = polyval(p1, x);
            p2 = polyfit(x, y, 2); y2 = polyval(p2, x);
            R2_1 = 1 - sum((y - y1).^2)/sum((y-mean(y)).^2);
            R2_2 = 1 - sum((y - y2).^2)/sum((y-mean(y)).^2);
            if R2_2 > R2_1, best_fit = y2; else, best_fit = y1; end
            plot(x, best_fit, 'r--', 'LineWidth', 1.3);

            title(contrast_title(cond), 'FontSize', 10, 'Interpreter','tex');
            if c == n_this_row
                xlabel('Compartment');
            end
            if c == 1
                ylabel('Proportion');
            end
            set(gca, ...
                'Box', 'off', ...
                'TickDir', 'out', ...
                'TickLength', [.01 .01], ...
                'XColor', [0 0 0], ...
                'YColor', [0 0 0], ...
                'XTick', 1:n_bins, ...
                'YTick', 0:.1:.4, ...
                'LineWidth', .75);
            axis square
            axis([0 n_bins+1 0 .4])
            hold off
        else
            axis off
        end
        if c == 1 & r == 1
            text(-2, .5, 'F', 'Fontsize',14, 'FontWeight', 'bold')
        end
    end
end


Nvals = [1 2 4 8 16];
laterality_all = [];
group_N = [];

for i = 1:numel(Nvals)
    cond_key = sprintf('0_%d_l', Nvals(i));
    if isKey(all_laterality, cond_key)
        vals = all_laterality(cond_key);
        laterality_all = [laterality_all; vals(:)];
        group_N = [group_N; repmat(Nvals(i), numel(vals), 1)];
    end
end

% -- Combine into a table if you want --
T = table(group_N, laterality_all, 'VariableNames', {'N','LateralityIndex'});

% % --- For each N, test if LateralityIndex > 0 (rightward preference, 1-sided) ---
% fprintf('\n0 vs N (small): Test if laterality index > 0 for each N\n');
% for i = 1:numel(Nvals)
%     cond_key = sprintf('0_%d_s', Nvals(i));
%     if isKey(all_laterality, cond_key)
%         vals = all_laterality(cond_key);
%         [p,~,stats] = signrank(vals, 0, 'tail', 'right');
%         fprintf('N=%2d: median=%.3f, n=%2d, p(one-sided)=%.4f\n', ...
%             Nvals(i), median(vals), numel(vals), p);
%     end
% end

% --- (Optional) Test overall effect across all N (collapsed) ---
[p_overall,~,stats_overall] = signrank(laterality_all, 0, 'tail', 'right');
fprintf('\nOverall 0 vs N (small): median=%.3f, n=%d, p(one-sided)=%.4f\n', ...
    median(laterality_all), numel(laterality_all), p_overall);


pvals = zeros(numel(Nvals),1);

for i = 1:numel(Nvals)
    cond_key = sprintf('0_%d_l', Nvals(i));
    if isKey(all_laterality, cond_key)
        vals = all_laterality(cond_key);
        [p,~,~] = signrank(vals, 0, 'tail', 'right');
        pvals(i) = p;
        fprintf('N=%2d: median=%.3f, n=%2d, p(one-sided)=%.4f\n', ...
            Nvals(i), median(vals), numel(vals), p);
    end
end

% Bonferroni
alpha = 0.05;
alpha_bonf = alpha / numel(Nvals);

fprintf('\nBonferroni corrected significance level: %.4f\n', alpha_bonf);
for i = 1:numel(Nvals)
    sig = pvals(i) < alpha_bonf;
    fprintf('N=%2d: p=%.4f %s\n', Nvals(i), pvals(i), ternary(sig,'*',''));
end

filefolder = 'I:\fishQuantity\stimuli\pages_halves\cropped'
% --- Stimulus bitmaps (adjust filenames/paths to yours) ---
stimfolder = fullfile(filefolder);     % e.g., I:\fishQuantity\files\stimuli
img1 = imread(fullfile(stimfolder, 'page_0024_R.png'));
img2 = imread(fullfile(stimfolder, 'page_0018_L.png'));

% If the images are indexed/grayscale, imagesc works; keep aspect; no axes
subplot(2,7,8);  % ← first empty slot
imshow(img2(920:1170, 570:1150)); axis image off;            % or: imshow(img1) if you prefer
title('Predator-sized fish', 'FontSize', 9);   % optional; you said you’ll relabel later
text(-110, -120, 'D', 'Fontsize',14, 'FontWeight','bold');
text(5, 290, 'Stimulus dimensions:', 'FontSize', 9);
text(5, 360, 'Length: 47mm', 'FontSize', 9);
text(5, 430, 'Width: 18mm', 'FontSize', 9);


subplot(2,7,9);  % ← second empty slot
imshow(img1); axis image off;
hold on
plot([0 size(img1,2)],[0 0],'k-')
plot([0 size(img1,2)],[size(img1,1) size(img1,1)],'k-')
plot([size(img1,2) size(img1,2)],[0 size(img1,1)],'k-')
plot([0 0],[0 size(img1,1)],'k-')
text(-190, size(img1,1)+150, 'Example stimulus: 8 fish', 'FontSize', 9);
text(-220, -210, 'E', 'Fontsize',14, 'FontWeight','bold');

hold off


hold off

function out = ternary(cond, tval, fval)
    if cond, out = tval; else, out = fval; end
end

function [out, all_laterality] = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_large, do_bias_correction)
    n = length(labels_large);
    all_conditions = containers.Map('KeyType','char','ValueType','any');
    all_laterality = containers.Map('KeyType','char','ValueType','any');
    roi_width = roi_rect(3);

    for s = 1:length(all_subjects)
        trials = all_subjects{s};
        subj_core = subject_ids{s};

        load(fullfile(resultsfolder, [subj_core '_detections.mat']), 'subject_detections');
        for tIdx = 1:length(trials)
            fish_traj = subject_detections{tIdx};
            if isempty(fish_traj), continue; end
            traj_x = fish_traj(:,2);
            int1 = trials{tIdx,4};
            int2 = trials{tIdx,5};
            size_letter = trials{tIdx,6};

            x_edges = linspace(1, roi_width+1, n_bins+1);
            compartment_assignment = discretize(traj_x, x_edges);
            prop_per_bin = zeros(1, n_bins);
            for b = 1:n_bins
                prop_per_bin(b) = sum(compartment_assignment == b) / numel(compartment_assignment);
            end

            do_swap = int1 > int2;
            if do_swap
                tmp = int1; int1 = int2; int2 = tmp;
                traj_x = roi_rect(3) - traj_x + 1;
                compartment_assignment = n_bins + 1 - compartment_assignment;
                prop_per_bin = fliplr(prop_per_bin);
            end
            cond_key = sprintf('%d_%d_%s', int1, int2, size_letter);
            lat_idx = mean(prop_per_bin(end-1:end)) - mean(prop_per_bin(1:2));

            result.prop_per_bin = prop_per_bin;
            result.int1 = int1;
            result.int2 = int2;
            result.size_letter = size_letter;
            result.lat_idx = lat_idx;

            if ~isKey(all_conditions, cond_key)
                all_conditions(cond_key) = {result};
            else
                tmp = all_conditions(cond_key);
                tmp{end+1} = result;
                all_conditions(cond_key) = tmp;
            end
            if ~isKey(all_laterality, cond_key)
                all_laterality(cond_key) = lat_idx;
            else
                tmp = all_laterality(cond_key);   % get the current array
                tmp(end+1) = lat_idx;             % append
                all_laterality(cond_key) = tmp;   % set back to the map
            end
        end
    end

    % ---- Collate by condition and bias correct if needed ----
    all_bins = cell(n, n_bins);
    for i = 1:n
        cond = labels_large{i};
        if ~isKey(all_conditions, cond), continue; end
        trials = all_conditions(cond);
        trialwise_profiles = [];
        for k = 1:length(trials)
            ppb = trials{k}.prop_per_bin;
            split_cond = split(cond, '_');
            n1 = str2double(split_cond{1});
            n2 = str2double(split_cond{2});
            is_equal = n1 == n2;
            if do_bias_correction && is_equal
                sum_left = sum(ppb(1:2));
                sum_right = sum(ppb(end-1:end));
                if sum_right < sum_left
                    ppb = fliplr(ppb);
                end
            end
            trialwise_profiles = [trialwise_profiles; ppb];
        end
        for b = 1:n_bins
            if ~isempty(trialwise_profiles)
                all_bins{i, b} = trialwise_profiles(:, b);
            end
        end
    end
    means_bins = nan(n, n_bins);
    sems_bins  = nan(n, n_bins);
    for i = 1:n
        for b = 1:n_bins
            vals = all_bins{i, b};
            if ~isempty(vals)
                means_bins(i, b) = mean(vals);
                sems_bins(i, b)  = std(vals)/sqrt(numel(vals));
            end
        end
    end
    out.means = means_bins;
    out.sems = sems_bins;
    out.all_bins = all_bins;
    out.all_conditions = all_conditions; 
end


function t = contrast_title(label, corrected)
    parts = split(label, '_');
    t = sprintf('%s vs %s', parts{1}, parts{2});
    if nargin > 1 && corrected
        t = [t, ' (corrected)'];
    end
end
