%% analysisFishQuantity4.m
% Small-fish overview + ratio/laterality summaries (Stage-1 descriptive + export)
% -------------------------------------------------------------------------
%
% Prerequisite:
%   This script expects `all_subjects` and `subject_ids` in the base workspace
%   (e.g., created by analysisFishQuantity7b.m).
%
% clc

do_bias_correction = false; % << Set to true for trial-wise bias correction

filefolder = 'I:\fishQuantity\files';
programfolder = 'I:\fishQuantity\programs';
resultsfolder = 'I:\fishQuantity\results';
figurefolder = 'I:\fishQuantity\figures';
addpath(genpath('I:\fishQuantity'))

load(fullfile(resultsfolder, 'roi_rect.mat'), 'roi_rect');
% subject_files = dir(fullfile(resultsfolder, '*_sequences.csv'));

% all_subjects = cell(1, numel(subject_files));
% subject_ids  = cell(1, numel(subject_files));

% for s = 1:numel(subject_files)
%     T = readtable(fullfile(resultsfolder, subject_files(s).name));
%     all_subjects{s} = table2cell(T);
%     [~, base, ~] = fileparts(subject_files(s).name);
%     subject_ids{s} = regexprep(base, '_sequences$', '');
% end
% 
% imageon = 0; % Enable video writing
% output_folder = 'output_videos';
% if imageon && ~exist(output_folder, 'dir')
%     mkdir(output_folder);
% end

all_conditions = containers.Map('KeyType','char','ValueType','any');

roi_width = roi_rect(3);
n_bins = 7; % Number of compartments

% --- keep only subjects with the full 30-trial protocol ---
nTrialsPerSubj = cellfun(@(x) size(x,1), all_subjects);   % e.g., 30 or 10
idxFull = find(nTrialsPerSubj == 30);

if numel(idxFull) < 12
    error('Expected at least 12 subjects with 30 trials, found %d.', numel(idxFull));
end

% If you truly want "the first 12" full-protocol subjects:
idxUse = idxFull(1:12);

% --- main loop over ONLY those subjects ---
for ii = 1:numel(idxUse)
    s = idxUse(ii);
    trials = all_subjects{s};
    subj_core = subject_ids{s};
%     load(fullfile(resultsfolder, [subj_core '_detections.mat']), 'subject_detections');

%     tokens = regexp(subject_files(s).name, '^(QF\d+)_', 'tokens', 'once');
%     if isempty(tokens)
%         error('Could not extract core subject ID from file: %s', subject_files(s).name);
%     end
%     subj_core = tokens{1};
    load(fullfile(resultsfolder, [subj_core '_detections.mat']), 'subject_detections');
    for tIdx = 1:length(trials)
        fish_traj = subject_detections{tIdx};
        if isempty(fish_traj), continue; end

        traj_x = fish_traj(:,2);
        traj_y = fish_traj(:,3);

        int1 = trials{tIdx,4};
        int2 = trials{tIdx,5};
        size_letter = trials{tIdx,6};

        dist_left = traj_x - 1;
        dist_right = roi_width - traj_x;
        [~, side_assignment] = min([dist_left, dist_right], [], 2);

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
            sa_new = side_assignment;
            sa_new(side_assignment == 1) = 2;
            sa_new(side_assignment == 2) = 1;
            side_assignment = sa_new;
        end

        cond_key = sprintf('%d_%d_%s', int1, int2, size_letter);
        time_left = mean(side_assignment==1);
        time_right = mean(side_assignment==2);

        result.traj_x = traj_x;
        result.traj_y = traj_y;
        result.time_left = time_left;
        result.time_right = time_right;
        result.side_assignment = side_assignment;
        result.compartment_assignment = compartment_assignment;
        result.prop_per_bin = prop_per_bin;
        result.tIdx = tIdx;
        result.subject = subject_ids{s};
        result.int1 = int1;
        result.int2 = int2;
        result.size_letter = size_letter;

        if ~isKey(all_conditions, cond_key)
            all_conditions(cond_key) = {result};
        else
            tmp = all_conditions(cond_key);
            tmp{end+1} = result;
            all_conditions(cond_key) = tmp;
        end
    end
end

%% ===== EXPORT trial-level table for SMALL fish (for later interaction model) =====
% This uses the raw per-trial prop_per_bin stored in all_conditions.
% IMPORTANT: no equal-number "bias correction" here (inferential dataset).

all_keys = keys(all_conditions);

% Keep only small-fish conditions (suffix "_s")
small_keys = all_keys(endsWith(all_keys, '_s'));

% Preallocate (simple grow is fine at this scale, but we keep it tidy)
subjC   = {};
condC   = {};
int1V   = [];
int2V   = [];
log2V   = [];
latV    = [];

for k = 1:numel(small_keys)
    key   = small_keys{k};
    trials_cell = all_conditions(key);   % cell array of result structs for this condition

    for j = 1:numel(trials_cell)
        r = trials_cell{j};

        % Defensive checks
        if ~isfield(r,'prop_per_bin') || isempty(r.prop_per_bin), continue; end
        if ~isfield(r,'subject') || isempty(r.subject), continue; end

        ppb = r.prop_per_bin;

        % Laterality index (raw): mean rightmost 2 bins - mean leftmost 2 bins
        lat = mean(ppb(end-1:end)) - mean(ppb(1:2));

        % Predictors
        i1 = r.int1;
        i2 = r.int2;
        lr = log2(i2 / i1);

        % Store
        subjC{end+1,1} = r.subject; %#ok<SAGROW>
        condC{end+1,1} = key;       %#ok<SAGROW>
        int1V(end+1,1) = i1;        %#ok<SAGROW>
        int2V(end+1,1) = i2;        %#ok<SAGROW>
        log2V(end+1,1) = lr;        %#ok<SAGROW>
        latV(end+1,1)  = lat;       %#ok<SAGROW>
    end
end

Tsmall = table( ...
    categorical(subjC), categorical(repmat({'small'}, numel(latV), 1)), ...
    categorical(condC), int1V, int2V, log2V, latV, ...
    'VariableNames', {'subject','size','cond_key','int1','int2','log2_ratio','lat_idx'});

% Save for later (next step: create Tlarge similarly, then concatenate)
save(fullfile(resultsfolder, 'Tsmall_trials.mat'), 'Tsmall');
writetable(Tsmall, fullfile(resultsfolder, 'Tsmall_trials.csv'));

fprintf('Saved Tsmall with %d trials across %d subjects.\n', height(Tsmall), numel(unique(Tsmall.subject)));


all_keys = keys(all_conditions);
n_cond = length(all_keys);

desired_order = { ...
    '1_1',  '1_2',  '1_4',  '1_8',  '1_16', ...
    '2_2',  '2_4',  '2_8',  '2_16', ...
    '4_4',  '4_8',  '4_16', ...
    '8_8',  '8_16', ...
    '16_16'};

labels_small = strcat(desired_order, '_s');
n = length(desired_order);

% -- Trial-wise bias-corrected aggregation --
all_bins_s = cell(n, n_bins);
laterality_trialwise = cell(n,1); % for trialwise laterality
for i = 1:n
    idx = find(strcmp(all_keys, labels_small{i}));
    if isempty(idx), continue; end
    trials = all_conditions(labels_small{i});
    trialwise_profiles = [];
    trialwise_laterality = [];
    for k = 1:length(trials)
        ppb = trials{k}.prop_per_bin;
        split_cond = split(labels_small{i}, '_');
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
        % Trialwise laterality index (mean rightmost 2 bins - mean leftmost 2 bins)
        trialwise_laterality = [trialwise_laterality; mean(ppb(end-1:end)) - mean(ppb(1:2))];
    end
    for b = 1:n_bins
        if ~isempty(trialwise_profiles)
            all_bins_s{i, b} = trialwise_profiles(:, b);
        end
    end
    if ~isempty(trialwise_laterality)
        laterality_trialwise{i} = trialwise_laterality;
    end
end

means_bins_s = nan(n, n_bins);
sems_bins_s  = nan(n, n_bins);
for i = 1:n
    for b = 1:n_bins
        vals = all_bins_s{i, b};
        if ~isempty(vals)
            means_bins_s(i, b) = mean(vals);
            sems_bins_s(i, b)  = std(vals)/sqrt(numel(vals));
        end
    end
end
valid_idx_s = ~all(isnan(means_bins_s), 2);
means_s = means_bins_s(valid_idx_s,:);
sems_s  = sems_bins_s(valid_idx_s,:);
xlabels_s = labels_small(valid_idx_s);

row_names = { ...
    {'1_1_s','1_2_s','1_4_s','1_8_s','1_16_s'}, ...
    {'2_2_s','2_4_s','2_8_s','2_16_s'}, ...
    {'4_4_s','4_8_s','4_16_s'}, ...
    {'8_8_s','8_16_s'}, ...
    {'16_16_s'} ...
    };

n_rows = numel(row_names);
n_cols = max(cellfun(@numel, row_names));

if do_bias_correction
    suffix = '_biascorr';
else
    suffix = '_raw';
end

n_profiles = numel(desired_order);
ratios = nan(n_profiles,1);
slopes = nan(n_profiles,1);
curvatures = nan(n_profiles,1);
laterality = nan(n_profiles,1);
labels = strings(n_profiles,1);
profile_idx = 1;

% figure('Position',[100 200 1100 1100])
plot_idx = 1;
for r = 1:n_rows
    n_this_row = numel(row_names{r});
    for c = 1:n_this_row
%         subplot(n_rows, n_cols, (r-1)*n_cols + c);
        cond = row_names{r}{c};
        idx = find(strcmp(xlabels_s, cond));
%         if isempty(idx), axis off; continue; end

        meanx = means_s(idx,:);
        errs  = sems_s(idx,:);

        split_cond = split(cond, '_');
        n1 = str2double(split_cond{1});
        n2 = str2double(split_cond{2});
        is_equal = n1 == n2;
%         if is_equal && do_bias_correction
%             plot_color = [0 0 1];
%         else
%             plot_color = [0 0 0];
%         end

%         plot(1:n_bins, meanx, '.', 'Markersize', 8, ...
%             'MarkerEdgeColor', plot_color, 'MarkerFaceColor', plot_color);
%         hold on
%         for k = 1:length(errs)
%             plot([k k],[meanx(k) - errs(k) meanx(k) + errs(k)], '-', 'Color', plot_color)
%         end

        % ----- Polynomial Fitting -----
        x = (1:n_bins)';
        y = meanx(:);

        p1 = polyfit(x, y, 1);
        y1 = polyval(p1, x);
        R2_1 = 1 - sum((y - y1).^2)/sum((y-mean(y)).^2);

        p2 = polyfit(x, y, 2);
        y2 = polyval(p2, x);
        R2_2 = 1 - sum((y - y2).^2)/sum((y-mean(y)).^2);

        if R2_2 > R2_1
            best_fit = y2; best_type = 'quad'; best_R2 = R2_2;
        else
            best_fit = y1; best_type = 'lin'; best_R2 = R2_1;
        end

%         plot(x, best_fit, 'r--', 'LineWidth', 1.3);

        % -- Title: "1 vs 16" style, no coefficients --
        contrast_str = sprintf('%s vs %s', split_cond{1}, split_cond{2});
%         title(contrast_str, 'FontSize', 11, 'Interpreter','tex');

        % --- Store ratio, slope, curvature, laterality ---
        ratios(profile_idx) = min(n1,n2)/max(n1,n2);
        labels(profile_idx) = contrast_str;
        slopes(profile_idx) = p1(1);         % always store linear slope
        curvatures(profile_idx) = p2(1);     % always store quadratic coefficient
        % laterality (use *mean* of trialwise laterality indices)
        idx_lat = find(strcmp(labels_small, cond));
        if ~isempty(laterality_trialwise{idx_lat})
            laterality(profile_idx) = mean(laterality_trialwise{idx_lat});
        end
        profile_idx = profile_idx + 1;

        % --- Formatting ---
%         if ismember(plot_idx,[5,9,12,14,15])
%             xlabel('Compartment');
%         end
%         if ismember(plot_idx,[1,6,10,13,15])
%             ylabel('Proportion');
%         end
%         set(gca, 'XTick', 1:n_bins);
%         set(gca, ...
%           'Box', 'off', ...
%           'TickDir', 'out', ...
%           'TickLength', [.01 .01], ...
%           'XColor', [0 0 0], ...
%           'YColor', [0 0 0], ...
%           'XTick', 1:n_bins, ...
%           'YTick', 0:.2:.4, ...
%           'LineWidth', .75);
%         axis square
%         axis([0 n_bins+1 0 .4])
%         plot_idx = plot_idx + 1;
    end
%     for c = n_this_row+1:n_cols
%         subplot(n_rows, n_cols, (r-1)*n_cols + c);
%         axis off;
%     end
end

% sgtitle(sprintf('Small fish: time per compartment, grouped by left value [%s]', ...
%     ternary(do_bias_correction, 'bias corrected', 'raw')));

% --- Save the figure ---
fig_name = sprintf('SmallFish_compartmentProfiles%s.png', suffix);
save_path = fullfile(figurefolder, fig_name);
saveas(gcf, save_path);

% --- OUTPUT: print or display as table ---
keep_idx = ~isnan(ratios); % only keep contrasts with data
fprintf('\nContrast\t\tRatio\t\tSlope\t\tCurvature\tLateralityIdx\n');
for i = 1:profile_idx-1
    if keep_idx(i)
        fprintf('%-10s\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n', ...
            labels(i), ratios(i), slopes(i), curvatures(i), laterality(i));
    end
end

results_table = table(labels(keep_idx), ratios(keep_idx), slopes(keep_idx), curvatures(keep_idx), laterality(keep_idx), ...
    'VariableNames', {'Contrast', 'Ratio', 'Slope', 'Curvature', 'LateralityIdx'});
disp(results_table);

% --- Save condition-level laterality summary ---
size_tag = 'small';  % <<< change to 'large' in the large-fish script
if ~exist('suffix','var') || isempty(suffix)
    if exist('do_bias_correction','var') && do_bias_correction, suffix = '_biascorr'; else, suffix = '_raw'; end
end

% Add convenience columns
results_table.Log2Ratio     = log2(results_table.Ratio);
results_table.Size          = repmat(string(size_tag), height(results_table), 1);
results_table.BiasCorrected = repmat(logical(exist('do_bias_correction','var') && do_bias_correction), height(results_table), 1);

% File paths
csv_path = fullfile(resultsfolder, sprintf('laterality_conditions_%s%s.csv', size_tag, suffix));
mat_path = fullfile(resultsfolder, sprintf('laterality_conditions_%s%s.mat',  size_tag, suffix));

% Save CSV
writetable(results_table, csv_path);

% Save MAT (with a tidy struct)
laterality_summary = struct( ...
    'Size',           size_tag, ...
    'BiasCorrected',  exist('do_bias_correction','var') && do_bias_correction, ...
    'Contrast',       results_table.Contrast, ...
    'Ratio',          results_table.Ratio, ...
    'Log2Ratio',      results_table.Log2Ratio, ...
    'Slope',          results_table.Slope, ...
    'Curvature',      results_table.Curvature, ...
    'LateralityIdx',  results_table.LateralityIdx, ...
    'n_bins',         n_bins, ...
    'roi_width',      roi_width ...
);
save(mat_path, 'laterality_summary');

fprintf('Saved condition-level laterality to:\n  %s\n  %s\n', csv_path, mat_path);


valid = ~isnan(results_table.Ratio) & ~isnan(results_table.Slope);

log_ratio = log2(results_table.Ratio(valid));
slope = results_table.Slope(valid);
curv = results_table.Curvature(valid);
lat_idx = results_table.LateralityIdx(valid);

lm_slope = fitlm(log_ratio, slope, 'linear');
fprintf('\nSlope vs. log2(Ratio):\n');
disp(lm_slope);

lm_curv = fitlm(log_ratio, curv, 'linear');
fprintf('\nCurvature vs. log2(Ratio):\n');
disp(lm_curv);

lm_lat = fitlm(log_ratio, lat_idx, 'linear');
fprintf('\nLaterality Index vs. log2(Ratio):\n');
disp(lm_lat);

yfit_lat = predict(lm_lat, log_ratio);
residuals_lat = lat_idx - yfit_lat;


if ~do_bias_correction
%     % --- Plot ---
%     % figure;
%     subplot(5, 6, 24)
%     hold off
%     plot(log_ratio, slope, 'k.', 'markersize',8);
%     hold on;
%     plot(log_ratio, predict(lm_slope, log_ratio), 'r-', 'LineWidth', 1);
%     xlabel('log_2(Ratio)');
%     ylabel('Slope');
%     % title('Slope vs. log_2(Ratio)');
%     set(gca, ...
%       'Box', 'off', ...
%       'TickDir', 'out', ...
%       'TickLength', [.01 .01], ...
%       'XColor', [0 0 0], ...
%       'YColor', [0 0 0], ...
%       'XTick', -4:1:0, ...
%       'YTick', -0.02:.02:.04, ...
%       'LineWidth', .75);
%     axis square
%     axis([-4.5 .5 -0.03 .045])
%     text(-5.5, 0.053, 'E', 'Fontsize',14, 'FontWeight', 'bold')
% 
% %     subplot(5, 6, 23)
% %     hold off
% %     plot(log_ratio, curv, 'k.', 'markersize',8);
% %     hold on;
% %     plot(log_ratio, predict(lm_curv, log_ratio), 'r-', 'LineWidth', 1);
% % %     xlabel('log_2(Ratio)');
% %     ylabel('Curvature');
% %     % title('Curvature vs. log_2(Ratio)');
% %     % grid on
% %     set(gca, ...
% %       'Box', 'off', ...
% %       'TickDir', 'out', ...
% %       'TickLength', [.01 .01], ...
% %       'XColor', [0 0 0], ...
% %       'YColor', [0 0 0], ...
% %       'XTick', -4:1:0, ...
% %       'YTick', [-15:5: 10]*10^-3, ...
% %       'LineWidth', .75);
% %     axis square
% %     axis([-4.5 .5 [-15 7]*10^-3])
% %     text(-5.5, 9.5*10^-3, 'E', 'Fontsize',14, 'FontWeight', 'bold')


%     subplot(5, 6, 23)
%     hold off
%     plot(log_ratio, lat_idx, 'k.', 'markersize',8);
%     hold on;
%     plot(log_ratio, predict(lm_lat, log_ratio), 'r-', 'LineWidth', 1);
%     xlabel('log_2(Ratio)');
%     ylabel('Laterality Index');
%     % title('Laterality Index vs. log_2(Ratio)');
%     set(gca, ...
%       'Box', 'off', ...
%       'TickDir', 'out', ...
%       'TickLength', [.01 .01], ...
%       'XColor', [0 0 0], ...
%       'YColor', [0 0 0], ...
%       'XTick', -4:1:0, ...
%       'YTick', [-.1:.1:.4], ...
%       'LineWidth', .75);
%     axis square
%     axis([-4.5 .5 -.15 .25])
%     text(-5.5, .285, 'D', 'Fontsize',14, 'FontWeight', 'bold')
% ---------- Panel D: LI vs log2(Right/Left) with LME population fit ----------
% condition means for compact scatter

   
%     %% -------- Panel F: Equal-number trials (per-trial LI) --------
%     % Uses: laterality_trialwise, labels_small
%     equal_labs = {'1_1_s','2_2_s','4_4_s','8_8_s','16_16_s'};
%     Ns         = [1 2 4 8 16];
% 
%     % Collect per-trial LIs and build jittered x-positions
%     x = []; y = []; grp = [];
%     rng(7);                               % fixed jitter for reproducibility
%     for i = 1:numel(Ns)
%         li = [];                          % default empty
%         li_idx = find(strcmp(labels_small, equal_labs{i}), 1);
%         if ~isempty(li_idx) && ~isempty(laterality_trialwise{li_idx})
%             li = laterality_trialwise{li_idx};
%             n  = numel(li);
%             x  = [x; Ns(i) + (rand(n,1)-0.5)*0.6];  %#ok<AGROW>  % jitter width
%             y  = [y; li(:)];                         %#ok<AGROW>
%             grp= [grp; repmat(Ns(i), n, 1)];        %#ok<AGROW>
%         end
%     end
% 
%     subplot(5,6,28); cla; hold on
%     % points
%     plot(x, y, 'k.', 'markersize', 8);
% 
%     % mean ± 95% CI per N
%     ux = unique(grp);
%     for ii = 1:numel(ux)
%         idx = grp==ux(ii);
%         m   = mean(y(idx), 'omitnan');
%         s   = std(y(idx),  'omitnan');
%         n   = sum(idx);
%         ci  = 1.96 * s / max(sqrt(n),1);          % 95% CI
%         plot([ux(ii) ux(ii)], [m-ci m+ci], 'r-', 'LineWidth', 1.3);
%         plot(ux(ii), m, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 4);
%     end
% 
%     yline(0, ':k');
%     xlim([0 17]); set(gca,'XTick',Ns);
%     xlabel('# fish (per side)'); ylabel('Laterality Index');
%     title('Equal-number trials');
%     set(gca, ...
%       'Box','off','TickDir','out','TickLength',[.01 .01], ...
%       'XColor',[0 0 0],'YColor',[0 0 0], ...
%       'LineWidth', .75);
%     axis square
%     text(0-1.5, max(get(gca,'YLim'))*1.1, 'F', 'Fontsize',14, 'FontWeight','bold');


end


% figure('Position',[200,200,1000,400]); % Wide layout for more panels
% 
% % Panel A: Residuals vs Fitted
% subplot(1,2,1);
% scatter(yfit_lat, residuals_lat, 60, 'filled');
% xlabel('Fitted Laterality Index');
% ylabel('Residual');
% title('A. Residuals vs. Fitted');
% yline(0, '--k');
% set(gca, 'FontSize', 13);
% 
% % Panel B: Residuals vs log2(Ratio)
% subplot(1,2,2);
% scatter(log_ratio, residuals_lat, 60, 'filled');
% xlabel('log_2(Ratio)');
% ylabel('Residual');
% title('B. Residuals vs. log_2(Ratio)');
% yline(0, '--k');
% set(gca, 'FontSize', 13);
% 
% % You can add more panels below as needed
% % e.g., subplot(2,3,3), etc.
% 
% sgtitle('Supplementary Figure S1: Model Diagnostics and Exploratory Analyses');
% 
% saveas(gcf, fullfile(figurefolder, 'SupplementaryFigure_S1.png'));

% --- Helper ternary function ---
function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end
