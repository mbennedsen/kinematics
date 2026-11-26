function verify_snapshot
% VERIFY_SNAPSHOT Check the frozen data, saved results and supplied MATLAB code.
root=fileparts(fileparts(mfilename('fullpath')));
old=path;cleanup=onCleanup(@()path(old));
addpath(fullfile(root,'replication','engines','methods'),'-begin');
T=readtable(fullfile(root,'data','SNAPSHOT.csv'),'TextType','string');
for k=1:height(T)
 file=fullfile(root,T.RelativePath(k));assert(isfile(file),'Missing file: %s',file);
 info=dir(file);assert(info.bytes==T.Bytes(k),'Size mismatch: %s',file);
 hash=file_sha256(file);
 assert(strcmp(hash,T.SHA256(k)),'replication:SnapshotHashMismatch', ...
     'Changed snapshot file: %s',file);
end
fprintf('SNAPSHOT PASS: %d data, result and code files verified.\n',height(T));
end
