function result = FigureS5(mode, options)
% FIGURES5 NOAA full-sample paths.
% FigureS5('saved') reproduces the paper from outputs/saved (default).
% FigureS5('reestimate') estimates again using the frozen data and specifications.
% Optional options.Workers (default 1); options.RunID reuses a new-estimation run.
% Re-estimation never loads outputs/saved. See README.md for cost and outputs.
if nargin < 1, mode = 'saved'; end
if nargin < 2, options = struct(); end
result = replicate_item('FigureS5', mode, options);
end
