%% analysisFishQuantity7b.m
% Small-fish compartment profiles (Stage-1 descriptive)
% -------------------------------------------------------------------------
%
%% Small-fish compartment profiles (A: corrected equal-number; B: raw all contrasts)
% Stage-1 (descriptive). No inferential stats here.
% Outputs:
%   - Figure: A/B panel grid
%   - CSV: fit_parameters_small.csv (QC/SI only)
%   - MAT: profiles_small.mat   (means, sems, Ns, labels) for downstream use

clear; clc; close all;

% --- Paths, loading, settings ---
filefolder    = 'I:\fishQuantity\files';
programfolder = 'I:\fishQuantity\programs';
resultsfolder = 'I:\fishQuantity\results';
figurefolder  = 'I:\fishQuantity\figures';
addpath(genpath('I:\fishQuantity'))
cd(resultsfolder);

S = load(fullfile(resultsfolder, 'roi_rect.mat'), 'roi_rect');
roi_rect  = S.roi_rect;
roi_width = roi_rect(3);

n_bins = 7; % number of x-compartments (fixed across panels)

% Canonical small-fish contrast order (left_right), then tag with '_s'
desired_order = { ...
    '1_1',  '1_2',  '1_4',  '1_8',  '1_16', ...
    '2_2',  '2_4',  '2_8',  '2_16', ...
    '4_4',  '4_8',  '4_16', ...
    '8_8',  '8_16', ...
    '16_16'};
labels_small = strcat(desired_order, '_s');
n = numel(desired_order);

% Subject trial files already exported earlier as QF##_all_trials.mat
files         = dir(fullfile(resultsfolder, 'QF*_all_trials.mat'));
subject_ids   = cell(1, numel(files));
all_subjects  = cell(1, numel(files));
for s = 1:numel(files)
    subj_id          = regexp(files(s).name, '(QF\d{2})_all_trials.mat', 'tokens', 'once');
    subject_ids{s}   = subj_id{1};
    tmp              = load(fullfile(resultsfolder, files(s).name));  % contains subj_data
    all_subjects{s}  = tmp.subj_data;
end

% --- Aggregate profiles: raw and equal-number corrected (for diagonal only) ---
results_raw  = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_small, false);
results_corr = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_small, true);

means_raw   = results_raw.means;
sems_raw    = results_raw.sems;
Ns_raw      = results_raw.N;           % trials per contrast
xlabels_raw = labels_small;            % preserved order

means_corr   = results_corr.means;
sems_corr    = results_corr.sems;
xlabels_corr = labels_small;

% --- Panel layout (A=left-most corrected equal-number; B=others raw) ---
row_names = { ...
    {'1_1_s','1_2_s','1_4_s','1_8_s','1_16_s'}, ...
    {'2_2_s','2_4_s','2_8_s','2_16_s'}, ...
    {'4_4_s','4_8_s','4_16_s'}, ...
    {'8_8_s','8_16_s'}, ...
    {'16_16_s'} ...
    };
n_rows = numel(row_names);
n_cols = max(cellfun(@numel, row_names)); % raw columns; +1 for corrected left column

hFig = figure('Position',[100 100 1400 1100]);

for r = 1:n_rows
    n_this_row = numel(row_names{r});

    % ---- LEFT-most: corrected diagonal (panel A instances per row) ----
    cond_diag = row_names{r}{1};
    idx_corr  = find(strcmp(xlabels_corr, cond_diag));
    subplot(n_rows, n_cols+1, (r-1)*(n_cols+1) + 1);

    if ~isempty(idx_corr)
        meanx = means_corr(idx_corr,:);
        errs  = sems_corr(idx_corr,:);

        plot(1:n_bins, meanx, '.', 'MarkerSize',10, ...
            'MarkerEdgeColor',[0 0 1], 'MarkerFaceColor',[0 0 1]); hold on
        for k = 1:numel(errs)
            plot([k k], [meanx(k)-errs(k), meanx(k)+errs(k)], '-', 'Color',[0 0 1]);
        end

        % Blue dashed: best of linear/quadratic fit (descriptive)
        x = (1:n_bins)'; y = meanx(:);
        p1 = polyfit(x,y,1);  y1 = polyval(p1,x);
        p2 = polyfit(x,y,2);  y2 = polyval(p2,x);
        R2_1 = 1 - sum((y-y1).^2)/sum((y-mean(y)).^2);
        R2_2 = 1 - sum((y-y2).^2)/sum((y-mean(y)).^2);
        plot(x, (R2_2>R2_1).*y2 + (R2_2<=R2_1).*y1, 'b--','LineWidth',1.3);

        title(contrast_title(cond_diag, true), 'FontSize',10, 'Interpreter','tex');
        if r == n_rows, xlabel('Compartment'); end
        ylabel('Proportion');
        set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
            'XColor',[0 0 0],'YColor',[0 0 0], ...
            'XTick',1:n_bins,'YTick',0:.1:.4,'LineWidth',.75);
        axis square; axis([0 n_bins+1 0 .4]); hold off
    else
        axis off
    end
    if r == 1
        text(-2, .5, 'A', 'FontSize',14,'FontWeight','bold');
    end

    % ---- Remaining columns in row: raw profiles (panel B instances) ----
    for c = 1:n_this_row
        cond    = row_names{r}{c};
        idx_raw = find(strcmp(xlabels_raw, cond));
        subplot(n_rows, n_cols+1, (r-1)*(n_cols+1) + c + 1);

        if ~isempty(idx_raw)
            meanx = means_raw(idx_raw,:);
            errs  = sems_raw(idx_raw,:);

            plot(1:n_bins, meanx, '.', 'MarkerSize',10, ...
                'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0]); hold on
            for k = 1:numel(errs)
                plot([k k], [meanx(k)-errs(k), meanx(k)+errs(k)], '-', 'Color',[0 0 0]);
            end

            % Red dashed: best of linear/quadratic fit (descriptive)
            x = (1:n_bins)'; y = meanx(:);
            p1 = polyfit(x,y,1);  y1 = polyval(p1,x);
            p2 = polyfit(x,y,2);  y2 = polyval(p2,x);
            R2_1 = 1 - sum((y-y1).^2)/sum((y-mean(y)).^2);
            R2_2 = 1 - sum((y-y2).^2)/sum((y-mean(y)).^2);
            plot(x, (R2_2>R2_1).*y2 + (R2_2<=R2_1).*y1, 'r--','LineWidth',1.3);

            title(contrast_title(cond), 'FontSize',10, 'Interpreter','tex');
            if c == n_this_row, xlabel('Compartment'); end
            if c == 1,         ylabel('Proportion');   end
            set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
                'XColor',[0 0 0],'YColor',[0 0 0], ...
                'XTick',1:n_bins,'YTick',0:.1:.4,'LineWidth',.75);
            axis square; axis([0 n_bins+1 0 .4]); hold off
        else
            axis off
        end

        if c == 1 && r == 1
            text(-2, .5, 'B', 'FontSize',14,'FontWeight','bold');
        end
    end
end

% --- Save figure (main panels only; inferential lines come in Stage-2 figure) ---
if ~exist(figurefolder,'dir'); mkdir(figurefolder); end
saveas(hFig, fullfile(figurefolder, 'SmallFish_compartmentProfiles_AB.png'));

%% QC/SI table: simple linear/quadratic fits per contrast (raw & corrected)
% Note: retained for supplemental rigor; main inference will use per-trial LME (Stage-2).
fit_table = table('Size',[n,8], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Contrast','Slope_Raw','Curvature_Raw','R2_Raw', ...
                      'Slope_Corr','Curvature_Corr','R2_Corr','N'});

x = (1:n_bins)';

for i = 1:n
    fit_table.Contrast(i) = string(xlabels_raw{i});

    % RAW
    y_raw = means_raw(i,:);
    if ~any(isnan(y_raw))
        p1 = polyfit(x, y_raw(:), 1);  y1 = polyval(p1, x);
        p2 = polyfit(x, y_raw(:), 2);  y2 = polyval(p2, x);
        R2_1 = 1 - sum((y_raw(:)-y1).^2) / sum((y_raw(:)-mean(y_raw(:))).^2);
        R2_2 = 1 - sum((y_raw(:)-y2).^2) / sum((y_raw(:)-mean(y_raw(:))).^2);
        fit_table.Slope_Raw(i)     = p1(1);
        fit_table.Curvature_Raw(i) = p2(1);
        fit_table.R2_Raw(i)        = max(R2_1, R2_2);
    else
        fit_table.Slope_Raw(i)     = NaN;
        fit_table.Curvature_Raw(i) = NaN;
        fit_table.R2_Raw(i)        = NaN;
    end

    % CORRECTED
    y_corr = means_corr(i,:);
    if ~any(isnan(y_corr))
        p1c = polyfit(x, y_corr(:), 1);  y1c = polyval(p1c, x);
        p2c = polyfit(x, y_corr(:), 2);  y2c = polyval(p2c, x);
        R2_1c = 1 - sum((y_corr(:)-y1c).^2) / sum((y_corr(:)-mean(y_corr(:))).^2);
        R2_2c = 1 - sum((y_corr(:)-y2c).^2) / sum((y_corr(:)-mean(y_corr(:))).^2);
        fit_table.Slope_Corr(i)     = p1c(1);
        fit_table.Curvature_Corr(i) = p2c(1);
        fit_table.R2_Corr(i)        = max(R2_1c, R2_2c);
    else
        fit_table.Slope_Corr(i)     = NaN;
        fit_table.Curvature_Corr(i) = NaN;
        fit_table.R2_Corr(i)        = NaN;
    end

    % N = trials per contrast (from aggregator)
    fit_table.N(i) = Ns_raw(i);
end

writetable(fit_table, fullfile(resultsfolder, 'fit_parameters_small.csv'));
disp(fit_table);

% Keep a tidy MAT with means/SEMs/Ns/labels for downstream use
profiles_small = struct( ...
    'labels',        labels_small, ...
    'means_raw',     means_raw, ...
    'sems_raw',      sems_raw, ...
    'means_corr',    means_corr, ...
    'sems_corr',     sems_corr, ...
    'N',             Ns_raw, ...
    'n_bins',        n_bins, ...
    'roi_width',     roi_width);
save(fullfile(resultsfolder,'profiles_small.mat'), 'profiles_small');

% Also keep a reference to the original per-subject trial tables if needed later
trials_small = all_subjects;   % rename for clarity
save(fullfile(resultsfolder,'behavior_trials_small.mat'),'trials_small');

%% ----------------------- helpers -----------------------
function out = extract_profiles(resultsfolder, all_subjects, subject_ids, roi_rect, n_bins, labels_small, do_bias_correction)
    % Aggregate per-trial x-occupancy into per-contrast means/SEMs.
    % If do_bias_correction=true, apply equal-number flipping on a per-trial basis.

    n         = numel(labels_small);
    roi_width = roi_rect(3);

    % Collect trials keyed by contrast "N1_N2_s/l"
    all_conditions = containers.Map('KeyType','char','ValueType','any');

    for s = 1:numel(all_subjects)
        trials    = all_subjects{s};
        subj_core = subject_ids{s};

        % Load detections for this subject
        S = load(fullfile(resultsfolder, [subj_core '_detections.mat']), 'subject_detections');
        if ~isfield(S, 'subject_detections') || isempty(S.subject_detections)
            warning('No detections for %s. Skipping subject.', subj_core);
            continue;
        end
        subject_detections = S.subject_detections;

        % --- robust trial count (rows), not numel ---
        if istable(trials)
            nTrials = height(trials);
        elseif iscell(trials)
            nTrials = size(trials, 1);   % rows = trials (expect m x k cell)
        else
            nTrials = size(trials, 1);
        end

        nDet  = numel(subject_detections);
        nUse  = min(nTrials, nDet);

        if nTrials ~= nDet
            warning('Trial/detection mismatch %s: trials(rows)=%d, det=%d. Using %d paired entries.', ...
                subj_core, nTrials, nDet, nUse);
        end

        for tIdx = 1:nUse
            % Guard: some subjects may have fewer trial rows than nUse (defensive)
            if tIdx > nTrials
                warning('Skipping tIdx=%d for %s: exceeds trial rows (%d).', tIdx, subj_core, nTrials);
                continue;
            end

            fish_traj = subject_detections{tIdx};
            if isempty(fish_traj) || size(fish_traj,2) < 2
                continue;
            end

            traj_x = fish_traj(:,2);

            % --- read left/right numbers and size tag (table or cell) ---
            if istable(trials)
                int1        = trials{tIdx, 4};
                int2        = trials{tIdx, 5};
                size_letter = trials{tIdx, 6};
            else % cell array
                int1        = trials{tIdx, 4};
                int2        = trials{tIdx, 5};
                size_letter = trials{tIdx, 6};
            end

            % Discretize x into compartments
            x_edges = linspace(1, roi_width+1, n_bins+1);
            compartment_assignment = discretize(traj_x, x_edges);
            ppb = zeros(1, n_bins);
            for b = 1:n_bins
                ppb(b) = sum(compartment_assignment == b) / numel(compartment_assignment);
            end

            % Normalise: ensure "more" is on the right
            if int1 > int2
                tmp = int1; int1 = int2; int2 = tmp;
                ppb = fliplr(ppb);
            end

            % Equal-number bias correction (trial-wise)
            if do_bias_correction && (int1 == int2)
                sum_left  = sum(ppb(1:2));
                sum_right = sum(ppb(end-1:end));
                if sum_right < sum_left
                    ppb = fliplr(ppb);
                end
            end

            cond_key = sprintf('%d_%d_%s', int1, int2, size_letter);
            rec = struct('prop_per_bin', ppb, 'int1', int1, 'int2', int2, 'size_letter', size_letter);

            if ~isKey(all_conditions, cond_key)
                all_conditions(cond_key) = {rec};
            else
                tmp = all_conditions(cond_key);
                tmp{end+1} = rec;
                all_conditions(cond_key) = tmp;
            end
        end
    end

    % Collate by requested label order
    all_bins = cell(n, n_bins);
    Ns       = zeros(n,1);

    for i = 1:n
        cond = labels_small{i};
        if ~isKey(all_conditions, cond), continue; end

        trials_list = all_conditions(cond);
        Ns(i) = numel(trials_list);

        trialwise_profiles = zeros(Ns(i), n_bins);
        for k = 1:Ns(i)
            trialwise_profiles(k,:) = trials_list{k}.prop_per_bin;
        end
        for b = 1:n_bins
            if ~isempty(trialwise_profiles)
                all_bins{i,b} = trialwise_profiles(:,b);
            end
        end
    end

    % Means/SEMs across trials
    means_bins = nan(n, n_bins);
    sems_bins  = nan(n, n_bins);
    for i = 1:n
        for b = 1:n_bins
            vals = all_bins{i,b};
            if ~isempty(vals)
                means_bins(i,b) = mean(vals, 'omitnan');
                sems_bins(i,b)  = std(vals,  'omitnan') / sqrt(numel(vals));
            end
        end
    end

    out.means    = means_bins;
    out.sems     = sems_bins;
    out.all_bins = all_bins;
    out.N        = Ns;
end



function t = contrast_title(label, corrected)
    parts = split(label, '_');
    t = sprintf('%s vs %s', parts{1}, parts{2});
    if nargin > 1 && corrected
        t = [t, ' (corrected)'];
    end
end
