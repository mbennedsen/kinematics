function result = TableS3(mode, options)
% TABLES3 Conditional continuations to 2 degrees C.
% TableS3('saved') reproduces the paper from outputs/saved (default).
% TableS3('reestimate') estimates again using the frozen data and specifications.
% Optional options.Workers (default 1); options.RunID reuses a new-estimation run.
% Re-estimation never loads outputs/saved. See README.md for cost and outputs.
if nargin < 1, mode = 'saved'; end
if nargin < 2, options = struct(); end
result = replicate_item('TableS3', mode, options);
end
