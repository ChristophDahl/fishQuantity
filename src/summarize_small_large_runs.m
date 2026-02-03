function [Gs, Gl] = summarize_small_large_runs(inDir, smallFile, largeFile)
%SUMMARIZE_SMALL_LARGE_RUNS  Aggregate stimulus-metric runs for SMALL vs LARGE sets.
%
% This function reads the per-image/per-page stimulus-metrics CSVs produced by
% your stimulus pipeline and aggregates them by numerosity (Nobj) separately
% for:
%   - SMALL run: pages 0002–0016
%   - LARGE run: pages 0017–0031
%
% Output:
%   Gs, Gl : tables grouped by Nobj with GroupCount and mean_/std_ columns for
%            metrics present in the input CSVs (e.g., Coverage, Density, NNmean).
%
% Side effects (unchanged):
%   Writes the aggregated summaries into:
%     <inDir>/summary_SMALL_by_Nobj.csv
%     <inDir>/summary_LARGE_by_Nobj.csv
%
% Usage:
%   inDir = "I:\fishQuantity\stimuli\pages_halves\cropped\vtrim_manual";
%   [Gs, Gl] = summarize_small_large_runs(inDir);
%
% Optional:
%   [Gs, Gl] = summarize_small_large_runs(inDir, ...
%       "run_0002_0016_results.csv", "run_0017_0031_results.csv");
%


inDir = string(inDir);

% ---- locate files (robust patterns unless explicitly provided) ----
if nargin < 2 || strlength(string(smallFile)) == 0
    smallPath = locate_file(inDir, ...
        ["run_0002_0016_results*.csv", "*small*0002*0016*.csv", "*0002*0016*.csv"]);
else
    smallPath = fullfile(inDir, string(smallFile));
end

if nargin < 3 || strlength(string(largeFile)) == 0
    largePath = locate_file(inDir, ...
        ["run_0017_0031_results*.csv", "*large*0017*0031*.csv", "*0017*0031*.csv"]);
else
    largePath = fullfile(inDir, string(largeFile));
end

assert(~(smallPath == ""), "Small run file not found in %s.", inDir);
assert(~(largePath == ""), "Large run file not found in %s.", inDir);

% ---- read tables ----
Tsmall = readtable(smallPath);
Tlarge = readtable(largePath);

% ---- summarize each run separately ----
Gs = summarize_by_Nobj(Tsmall);
Gl = summarize_by_Nobj(Tlarge);

% ---- save CSVs ----
outSmall = fullfile(inDir, "summary_SMALL_by_Nobj.csv");
outLarge = fullfile(inDir, "summary_LARGE_by_Nobj.csv");
writetable(Gs, outSmall);
writetable(Gl, outLarge);

fprintf("Saved:\n  %s\n  %s\n", outSmall, outLarge);

end

% ================= helpers (unchanged computations) =======================
function p = locate_file(inDir, patterns)
p = "";
for pat = string(patterns)
    dd = dir(fullfile(inDir, pat));
    if ~isempty(dd)
        [~, ix] = max([dd.datenum]);  % newest match
        p = fullfile(dd(ix).folder, dd(ix).name);
        return;
    end
end
end

function G = summarize_by_Nobj(T)
% Ensure we have Nobj; fallback to KeptN if needed
if ~ismember("Nobj", T.Properties.VariableNames)
    if ismember("KeptN", T.Properties.VariableNames)
        T.Nobj = T.KeptN;
    else
        error("Neither Nobj nor KeptN found in the table.");
    end
end
T.Nobj = double(T.Nobj);

% Metrics to summarize (use those present)
preferred = ["Density","Coverage","CentVar","CentStdX","CentStdY","NNmean", ...
             "AreaTot","AreaMean","AreaStd","PerimTot", ...
             "HullArea","HullPerim","PairMean", ...
             "MeanEccentricity","MeanSolidity","MeanCompactness"];
vars = intersect(preferred, string(T.Properties.VariableNames));

% Manual grouping (version-proof)
[Nuniq, ~, g] = unique(T.Nobj);
K = numel(Nuniq);
counts = accumarray(g, 1);

% Preallocate
G = table(Nuniq, counts, "VariableNames", {"Nobj","GroupCount"});

for v = vars
    x = T.(v);
    m = nan(K, 1);
    s = nan(K, 1);
    for i = 1:K
        xi = x(g == i);
        m(i) = mean(xi, "omitnan");
        s(i) = std(xi,  "omitnan");
    end
    G.("mean_" + v) = m;
    G.("std_"  + v) = s;
end

G = sortrows(G, "Nobj");
end
