function modelResults = fitOTSvsANSModels(all_trials, varargin)
% fitOTSvsANSModels  Alias wrapper around runOTSvsANS (kept for backwards compatibility).
%
% This wrapper exists solely for repository readability. It does not change
% any computation; it forwards all inputs to runOTSvsANS unchanged.
modelResults = runOTSvsANS(all_trials, varargin{:});
end
