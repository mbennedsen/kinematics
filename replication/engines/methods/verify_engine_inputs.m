function snapshot = verify_engine_inputs(study)
% Check the estimator code and frozen data before fitting or resuming a job.
work = fileparts(study);
manifest = fullfile(work,'SHA256SUMS.tsv');
fid = fopen(manifest,'r');
assert(fid >= 0,'Missing estimation manifest: %s',manifest);
guard = onCleanup(@() fclose(fid));
rows = textscan(fid,'%s%s','Delimiter','\t');
for j = 1:numel(rows{1})
    relative = rows{2}{j};
    assert(strcmp(file_sha256(fullfile(work,relative)),rows{1}{j}), ...
        'Estimation input or code changed: %s. Use a fresh RunID.',relative);
end
snapshot = file_sha256(manifest);
end
