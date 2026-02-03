%% lateralityIndexSmallFishQuantity2.m
% Compute per-trial laterality index (LI) for SMALL fish trials and fit a per-trial
% linear mixed-effects model:
%   lat_idx ~ 1 + log2_ratio + (1|subject)
%
% Prerequisites / inputs:
%   - resultsfolder contains:
%       * roi_rect.mat
%       * QF##_all_trials.mat  (each contains subj_data)
%       * QF##_detections.mat  (each contains subject_detections)
%
% Outputs (written to resultsfolder):
%   - laterality_trials_small.mat        (struct array all_trials)
%   - fixed_effects_small.csv            (fixed-effect table for manuscript)
%   - lme_small.mat                      (fitted model + R2)
%   - pred_line_small.csv                (population prediction line)
%
% Laterality analysis (SMALL fish): per-trial LI + mixed-effects model + panels C/D (+ optional F)

% clear all
% clc

% IMPORTANT: do NOT bias-correct equal-number trials for the per-trial LME
do_bias_correction_equalN = false;

filefolder    = 'I:\fishQuantity\files';
programfolder = 'I:\fishQuantity\programs';
resultsfolder = 'I:\fishQuantity\results';
figurefolder  = 'I:\fishQuantity\figures';
addpath(genpath('I:\fishQuantity'))

S = load(fullfile(resultsfolder, 'roi_rect.mat'), 'roi_rect');
roi_rect  = S.roi_rect;
roi_width = roi_rect(3);

% Ensure subjects are loaded (if not already in workspace)
if ~exist('all_subjects','var') || ~exist('subject_ids','var')
    files = dir(fullfile(resultsfolder, 'QF*_all_trials.mat'));
    subject_ids  = cell(1, numel(files));
    all_subjects = cell(1, numel(files));
    for s = 1:numel(files)
        subj_id            = regexp(files(s).name, '(QF\d{2})_all_trials.mat', 'tokens', 'once');
        subject_ids{s}     = subj_id{1};
        tmp                = load(fullfile(resultsfolder, files(s).name));  % contains subj_data
        all_subjects{s}    = tmp.subj_data;
    end
end

imageon = 0;  % Enable video writing if you need it
output_folder = 'output_videos';
if imageon && ~exist(output_folder, 'dir'); mkdir(output_folder); end

all_trials = [];               % per-trial collection
n_bins     = 7;                % number of compartments

%% --------- Build per-trial LI (small fish only), robust trial↔detection pairing ----------
for s = 1:length(all_subjects)
    trials    = all_subjects{s};
    subj_core = subject_ids{s};

    % Load detections and pair robustly
    D = load(fullfile(resultsfolder, [subj_core '_detections.mat']), 'subject_detections');
    if ~isfield(D,'subject_detections') || isempty(D.subject_detections)
        warning('No detections for %s; skipping.', subj_core);
        continue;
    end
    subject_detections = D.subject_detections;

    if istable(trials), nTrials = height(trials); else, nTrials = size(trials,1); end
    nDet = numel(subject_detections);
    nUse = min(nTrials, nDet);
    if nTrials ~= nDet
        warning('Trial/detection mismatch %s: trials(rows)=%d, det=%d. Using %d paired.', ...
            subj_core, nTrials, nDet, nUse);
    end

    for tIdx = 1:nUse
        fish_traj = subject_detections{tIdx};
        if isempty(fish_traj) || size(fish_traj,2) < 2, continue; end
        traj_x = fish_traj(:,2);

        % Read trial meta (left/right numbers, size)
        if istable(trials)
            int1        = trials{tIdx,4}; % left
            int2        = trials{tIdx,5}; % right
            size_letter = trials{tIdx,6};
        else
            int1        = trials{tIdx,4};
            int2        = trials{tIdx,5};
            size_letter = trials{tIdx,6};
        end

        % Small-fish only in this script
        if ~strcmp(size_letter,'s'), continue; end

        % 7-bin discretization
        x_edges = linspace(1, roi_width+1, n_bins+1);
        compartment_assignment = discretize(traj_x, x_edges);
        prop_per_bin = zeros(1, n_bins);
        for b = 1:n_bins
            prop_per_bin(b) = sum(compartment_assignment == b) / numel(compartment_assignment);
        end

        % Normalize so "more is on the right"
        if int1 > int2
            tmp = int1; int1 = int2; int2 = tmp;
            prop_per_bin = fliplr(prop_per_bin);
        end

        % DO NOT flip equal-number trials for the LME (keeps LI centered near 0)
        if do_bias_correction_equalN && (int1 == int2)
            sum_left  = sum(prop_per_bin(1:2));
            sum_right = sum(prop_per_bin(end-1:end));
            if sum_right < sum_left
                prop_per_bin = fliplr(prop_per_bin);
            end
        end

        % Laterality index: mean(rightmost 2) - mean(leftmost 2)
        lat_idx = mean(prop_per_bin(end-1:end)) - mean(prop_per_bin(1:2));

        % Store per-trial record (include subject for LME)
        all_trials = [all_trials; struct( ...
            'subject',      subj_core, ...
            'int1',         int1, ...
            'int2',         int2, ...
            'size_letter',  size_letter, ...
            'lat_idx',      lat_idx, ...
            'prop_per_bin', prop_per_bin )]; %#ok<AGROW>
    end
end

% Save per-trial data for reproducibility (and later combining with large-fish)
if ~exist(resultsfolder, 'dir'); mkdir(resultsfolder); end
save(fullfile(resultsfolder, 'laterality_trials_small.mat'), 'all_trials');

%% ---------- Panel C (descriptive): aggregate LI by absolute N on right ----------
abs_numbers   = [1 2 4 8 16];
lat_by_number = cell(numel(abs_numbers),1);

for i = 1:numel(all_trials)
    t = all_trials(i);
    idx = find(abs_numbers == t.int2, 1);
    if ~isempty(idx)
        lat_by_number{idx}(end+1) = t.lat_idx; %#ok<SAGROW>
    end
end

mean_lat = nan(1, numel(abs_numbers));
sem_lat  = nan(1, numel(abs_numbers));
for i = 1:numel(abs_numbers)
    vals = lat_by_number{i};
    if ~isempty(vals)
        mean_lat(i) = mean(vals);
        sem_lat(i)  = std(vals)/sqrt(numel(vals));
    end
end

subplot(5,6,23); cla; hold on
xa = abs_numbers; meanx = mean_lat; stderrorx = sem_lat;
plot(xa, meanx, '-', 'LineWidth', 0.75, 'Color', [0 0 0]);
plot(xa, meanx-stderrorx, '--', 'LineWidth', 0.75, 'Color', [0 0 0]);
plot(xa, meanx+stderrorx, '--', 'LineWidth', 0.75, 'Color', [0 0 0]);
xlabel({'# fish (right compartment)'});
ylabel('Laterality Index');
set(gca, 'Box','off','TickDir','out','TickLength',[.01 .01], ...
    'XColor',[0 0 0],'YColor',[0 0 0], 'XTick',[1,2,4,8,16], ...
    'YTick', -.2:.1:.2, 'LineWidth', .75);
axis square; axis([0.5 16.5 -.2 .2]);
text(-4.5, .20, 'E', 'Fontsize',14, 'FontWeight', 'bold');

%% ---------- Primary analysis: per-trial mixed-effects model (small fish) ----------
T = struct2table(all_trials);
T = T(strcmp(T.size_letter,'s') & T.int1>0 & T.int2>0, :);

% Predictor: log2(N_right/N_left) (0..4). For plotting we'll flip to -4..0 later.
T.log2_ratio = log2(T.int2 ./ T.int1);
T.subject    = categorical(T.subject);

% Mixed-effects model: LI ~ log2_ratio + (1|subject)
lme_small = fitlme(T, 'lat_idx ~ 1 + log2_ratio + (1|subject)');

% Fixed effects summary
fe  = lme_small.Coefficients;
ci  = coefCI(lme_small);
row = strcmp(fe.Name, 'log2_ratio');
beta     = fe.Estimate(row);
beta_se  = fe.SE(row);
beta_p   = fe.pValue(row);
beta_ciL = ci(row,1);
beta_ciU = ci(row,2);

% R^2 (try built-in; fallback to Nakagawa pseudo-R2)
try
    R2m = lme_small.Rsquared.Marginal;
    R2c = lme_small.Rsquared.Conditional;
catch
    [R2m, R2c] = pseudoR2_Nakagawa(lme_small, T);
end

% Console report
disp(lme_small);
fprintf('β per doubling (log2_ratio): %.4f (SE=%.4f), 95%% CI [%.4f, %.4f], p=%.4g\n', ...
    beta, beta_se, beta_ciL, beta_ciU, beta_p);
fprintf('R^2_marginal=%.3f, R^2_conditional=%.3f\n', R2m, R2c);

% Save stats for manuscript
fixed_effects_small = table(fe.Name, fe.Estimate, fe.SE, fe.DF, fe.tStat, fe.pValue, ...
    'VariableNames', {'Term','Estimate','SE','DF','t','p'});
fixed_effects_small.CI_low  = ci(:,1);
fixed_effects_small.CI_high = ci(:,2);
writetable(fixed_effects_small, fullfile(resultsfolder,'fixed_effects_small.csv'));
save(fullfile(resultsfolder,'lme_small.mat'), 'lme_small','fixed_effects_small','R2m','R2c');

% Predicted population line for plotting (pass BOTH vars to avoid predict() error)
xx_pos = (0:4)';                                      % 0..4 = log2(N_right/N_left)
subj_ref = repmat(T.subject(1), numel(xx_pos), 1);    % any valid subject (ignored when Conditional=false)
newT = table(xx_pos, subj_ref, 'VariableNames', {'log2_ratio','subject'});
yhat = predict(lme_small, newT, 'Conditional', false);

% Flip x to -4..0 for display (to match your figure style)
xx_plot = -xx_pos; % -4..0 = log2(min/max)
pred_small = table(xx_plot, yhat, 'VariableNames', {'log2ratio_plot','LI_hat'});
writetable(pred_small, fullfile(resultsfolder,'pred_line_small.csv'));

%% ---------- Panel D: LI vs log2(Ratio) with LME population fit ----------
% Compact scatter: condition means (still, inference comes from per-trial LME)
% conds = unique([T.int1 T.int2], 'rows');
% x_scatter = []; y_scatter = [];
% for i = 1:size(conds,1)
%     n1 = conds(i,1); n2 = conds(i,2);
%     mask = (T.int1==n1 & T.int2==n2);
%     if ~any(mask), continue; end
%     y_scatter(end+1,1) = mean(T.lat_idx(mask));    %#ok<SAGROW>
%     x_scatter(end+1,1) = log2(n1 / n2);            %#ok<SAGROW>  % negative values (−4..0)
% end
% 
% subplot(5,6,23); cla; hold on
% plot(x_scatter, y_scatter, 'k.', 'markersize', 8);   % condition means
% plot(xx_plot, yhat, 'r-', 'LineWidth', 1.0);         % LME population fit
% xlabel('log_2(Ratio)'); ylabel('Laterality Index');
% set(gca, 'Box','off','TickDir','out','TickLength',[.01 .01], ...
%     'XColor',[0 0 0],'YColor',[0 0 0], ...
%     'XTick', -4:1:0, 'YTick', -.1:.1:.4, 'LineWidth', .75);
% yline(0, ':k');
% axis square; axis([-4.5 .5 -.15 .25]);
% text(-5.5, .285, 'D', 'Fontsize',14, 'FontWeight', 'bold');

% ---------- Panel D: LI vs log2(Right/Left) with LME population fit ----------
% condition means for compact scatter
conds = unique([T.int1 T.int2], 'rows');
x_scatter = []; y_scatter = [];
for i = 1:size(conds,1)
    n1 = conds(i,1); n2 = conds(i,2);
    mask = (T.int1==n1 & T.int2==n2);
    if ~any(mask), continue; end
    y_scatter(end+1,1) = mean(T.lat_idx(mask));             %#ok<AGROW>
    x_scatter(end+1,1) = log2(n2 / n1);                     %#ok<AGROW> % 0..4
end

subplot(5,6,30); cla; hold on
plot(x_scatter, y_scatter, 'k.', 'markersize', 8);          % condition means
plot(xx_pos, yhat, 'r-', 'LineWidth', 1.0);                 % LME population fit (xx_pos = 0:4)

xlabel('log_2(Right/Left)','Interpreter','tex');
ylabel('Laterality Index');
set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
    'XColor',[0 0 0],'YColor',[0 0 0], ...
    'XTick',0:4,'YTick',-.1:.1:.4,'LineWidth',.75);
yline(0, ':k');
axis square; axis([-0.2 4.2 -.15 .25]);
text(-0.75, .285, 'H', 'Fontsize',14, 'FontWeight','bold');


%% ---------- Panel F: Equal-number trials (flip per trial; descriptive only) ----------
equalNs = [1 2 4 8 16];
xF = []; yF = []; grpF = [];
rng(7);  % reproducible jitter

for i = 1:numel(equalNs)
    N = equalNs(i);
    % equal numbers, small fish only
    mask = [all_trials.int1] == N & [all_trials.int2] == N ...
         & strcmp({all_trials.size_letter}, 's');
    if ~any(mask), continue; end
    ati = all_trials(mask);

    for k = 1:numel(ati)
        ppb = ati(k).prop_per_bin;
        l_end = mean(ppb(1:2));
        r_end = mean(ppb(end-1:end));
        li = r_end - l_end;

        % Flip ONLY equal-N trials so higher end occupancy is positive
        if r_end < l_end
            li = -li;
        end

        xF   = [xF; N + (rand(1)-0.5)*0.6]; %#ok<AGROW>
        yF   = [yF; li];                    %#ok<AGROW>
        grpF = [grpF; N];                   %#ok<AGROW>
    end
end

subplot(5,6,28); cla; hold on
plot(xF, yF, 'k.', 'markersize', 7);

% mean ± 95% CI per N
uN = unique(grpF);
for ii = 1:numel(uN)
    idx = grpF==uN(ii);
    m   = mean(yF(idx), 'omitnan');
    s   = std(yF(idx),  'omitnan');
    n   = sum(idx);
    ci  = 1.96 * s / max(sqrt(n),1);
    plot([uN(ii) uN(ii)], [m-ci m+ci], 'r-', 'LineWidth', 1.0);
    plot(uN(ii), m, 'r.', 'MarkerSize', 12);
end

yline(0, ':k');
xlim([0 17]); set(gca,'XTick',equalNs);
ylim([-0.75 0.75]);                       % wide scale for raw dispersion
xlabel('# fish per side');
ylabel('Laterality Index');              % mean ± 95% CI
set(gca,'Box','off','TickDir','out','LineWidth', .75,'TickLength',[.01 .01], ...
  'XColor',[0 0 0], 'YColor',[0 0 0], ...
  'YTick', -.5:.25:.5, 'LineWidth', .75);
axis square;
text(-3, 0.85, 'F', 'Fontsize',14, 'FontWeight','bold');

%% ---------- Panel G: Unequal-number trials (UNFLIPPED; grouped by log2 ratio) ----------
xG = []; yG = []; grpG = [];
rng(8);  % independent jitter for this panel

% select unequal-number trials, small fish only
maskUnequal = strcmp({all_trials.size_letter}, 's') & [all_trials.int1] > 0 & [all_trials.int2] > 0 ...
            & ([all_trials.int1] ~= [all_trials.int2]);

ati = all_trials(maskUnequal);

for k = 1:numel(ati)
    ppb = ati(k).prop_per_bin;
    li  = mean(ppb(end-1:end)) - mean(ppb(1:2));  % UNFLIPPED
    % bin by log2 ratio (should be 1..4 after normalization)
    bin = round(log2(ati(k).int2 / ati(k).int1));
    xG   = [xG; bin + (rand(1)-0.5)*0.35];  %#ok<AGROW>  % light jitter within each bin
    yG   = [yG; li];                         %#ok<AGROW>
    grpG = [grpG; bin];                      %#ok<AGROW>
end

subplot(5,6,29); cla; hold on
plot(xG, yG, 'k.', 'markersize', 7);

% mean ± 95% CI per log2 ratio bin
ub = unique(grpG);
for ii = 1:numel(ub)
    idx = grpG==ub(ii);
    m   = mean(yG(idx), 'omitnan');
    s   = std(yG(idx),  'omitnan');
    n   = sum(idx);
    ci  = 1.96 * s / max(sqrt(n),1);
    plot([ub(ii) ub(ii)], [m-ci m+ci], 'r-', 'LineWidth', 1.0);
    plot(ub(ii), m, 'r.', 'MarkerSize', 12);
end

yline(0, ':k');
xlim([0.5 4.5]); set(gca,'XTick',1:4);
ylim([-0.75 0.75]);                       % same wide scale for comparability with F
xlabel('log_2(Right/Left)');            % 1..4 (2×, 4×, 8×, 16×)
ylabel('Laterality Index');             % mean ± 95% CI
set(gca,'Box','off','TickDir','out','LineWidth', .75,'TickLength',[.01 .01], ...
  'XColor',[0 0 0], 'YColor',[0 0 0], ...
  'YTick', -.5:.25:.5, 'LineWidth', .75);axis square;
text(-.2, 0.85, 'G', 'Fontsize',14, 'FontWeight','bold');



% filefolder = 'I:\fishQuantity\stimuli\pages_halves\cropped'
% % --- Stimulus bitmaps (adjust filenames/paths to yours) ---
% stimfolder = fullfile(filefolder);     % e.g., I:\fishQuantity\files\stimuli
% img1 = imread(fullfile(stimfolder, 'page_0015_L.png'));
% img2 = imread(fullfile(stimfolder, 'page_0002_R.png'));
% 
% % If the images are indexed/grayscale, imagesc works; keep aspect; no axes
% subplot(10,6,30);  % ← first empty slot
% imshow(img2(890:990, 740:1050)); axis image off;            % or: imshow(img1) if you prefer
% title('Conspecific-sized fish', 'FontSize', 9);   % optional; you said you’ll relabel later
% text(-40, -10, 'C', 'Fontsize',14, 'FontWeight','bold');
% text(5, 120, 'Stimulus dimensions:', 'FontSize', 9);
% text(5, 150, 'Length: 25mm', 'FontSize', 9);
% text(5, 180, 'Width: 6mm', 'FontSize', 9);


% subplot(10,6,[36,42,48]);  % ← second empty slot
% imshow(img1); axis image off;
% hold on
% plot([0 size(img1,2)],[0 0],'k-')
% plot([0 size(img1,2)],[size(img1,1) size(img1,1)],'k-')
% plot([size(img1,2) size(img1,2)],[0 size(img1,1)],'k-')
% plot([0 0],[0 size(img1,1)],'k-')
% text(80, size(img1,1)+100, 'Example stimulus: 8 fish', 'FontSize', 9);
% text(-200, -150, 'D', 'Fontsize',14, 'FontWeight','bold');
% 
% hold off

%% ----------------- helpers -----------------
function [R2m, R2c] = pseudoR2_Nakagawa(lme, T)
% Pseudo R² (Nakagawa & Schielzeth) for linear mixed models
% Marginal:    Var(Fixed) / (Var(Fixed)+Var(Random)+Var(Resid))
% Conditional: (Var(Fixed)+Var(Random)) / total
y   = T.lat_idx;
yF  = predict(lme, T, 'Conditional', false); % fixed effects only
yC  = predict(lme, T, 'Conditional', true);  % fixed + random
rfx = yC - yF;                               % random-effects contribution
res = y  - yC;                               % residuals
varF = var(yF,  1);
varR = var(rfx, 1);
varE = var(res, 1);
den  = varF + varR + varE + eps;
R2m  = varF / den;
R2c  = (varF + varR) / den;
end
