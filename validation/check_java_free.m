function check_java_free
% CHECK_JAVA_FREE Check hashing, fresh estimation, resume and corruption guards.
% Start MATLAB with -nojvm, then addpath('validation'); check_java_free.
assert(~usejava('jvm'),'Run this check in MATLAB with -nojvm.');
root=fileparts(fileparts(mfilename('fullpath')));
old=path;restore=onCleanup(@()path(old));
addpath(fullfile(root,'replication'),fullfile(root,'validation'),'-begin');
check_hashing;
started=tic;verify_snapshot;
fprintf('SNAPSHOT TIMING: %.2f s.\n',toc(started));
check_snapshot_rejection(root);

% Stage through the public interface, using the same data and source manifests.
opt=struct('Workers',1,'RunID',['java_free_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
TableS4('reestimate',opt);
work=fullfile(root,'outputs','reestimated',opt.RunID,'work');
addpath(fullfile(work,'empirical'),fullfile(work,'methods'),'-begin');
first=run_empirical(5,3,1,{'quadratic','kernel'});
saved=readtable(fullfile(root,'outputs','saved','empirical','Empirical_Results.csv'),'TextType','string');
largest=0;
for k=1:height(first)
    z=saved(saved.Dataset==first.Dataset(k) & saved.Method==first.Method(k) & ...
        saved.Object==first.Object(k) & saved.Specification=="monthly",:);
    assert(height(z)==1);
    difference=abs(z.Estimate-first.Estimate(k));largest=max(largest,difference);
    assert(difference<1e-7,'Fresh point estimate does not match the archive.');
end
caches=dir(fullfile(work,'empirical','results','monthly','ERA5','*_B3','CompletedFit.mat'));
assert(numel(caches)==2);
before=cell(2,1);
for k=1:2,before{k}=file_sha256(fullfile(caches(k).folder,caches(k).name));end
% Re-enter the same RunID before resuming the fits.
TableS4('reestimate',opt);
second=run_empirical(5,3,1,{'quadratic','kernel'});
assert(isequaln(first,second),'Resumed summaries differ.');
for k=1:2
    assert(strcmp(before{k},file_sha256(fullfile(caches(k).folder,caches(k).name))), ...
        'Resume rewrote a completed fit.');
end
fprintf('ESTIMATION/RESUME PASS: 8 point estimates (max error %.3g), unchanged caches.\n',largest);

% Mutations occur only in this disposable run's staged copies.
check_engine_rejection(fullfile(work,'data','Monthly_HadISST.csv'),work);
check_engine_rejection(fullfile(work,'methods','estimate_quadratic_benchmark.m'),work);
verify_engine_inputs(fullfile(work,'empirical'));
fprintf('ENGINE CORRUPTION PASS: same-size data and source changes rejected.\n');
fprintf('JAVA-FREE PASS: %s; MATLAB %s; JVM=%d.\n',opt.RunID,version,usejava('jvm'));
end

function check_snapshot_rejection(root)
% Exercise the actual snapshot verifier in a small independent fixture.
fixture=tempname;mkdir(fixture);
old=path;cleanup=onCleanup(@()cleanup_fixture(fixture,old));
rep=fullfile(fixture,'replication');mkdir(rep);mkdir(fullfile(rep,'engines'));
mkdir(fullfile(fixture,'data'));
copyfile(fullfile(root,'replication','verify_snapshot.m'),rep);
copyfile(fullfile(root,'replication','engines','methods'),fullfile(rep,'engines','methods'));
addpath(fullfile(rep,'engines','methods'),'-begin');
file=fullfile(fixture,'fixture.bin');write_bytes(file,uint8([0 1 127 128 255]));
T=table("fixture.bin",5,string(file_sha256(file)), ...
    'VariableNames',{'RelativePath','Bytes','SHA256'});
writetable(T,fullfile(fixture,'data','SNAPSHOT.csv'));
addpath(rep,'-begin');clear verify_snapshot
verify_snapshot;
write_bytes(file,uint8([0 1 126 128 255]));
try
    verify_snapshot;
    error('test:ExpectedFailure','Changed snapshot was accepted.');
catch err
    assert(strcmp(err.identifier,'replication:SnapshotHashMismatch'),'%s',err.message);
end
fprintf('SNAPSHOT CORRUPTION PASS: same-size file change rejected.\n');
end

function cleanup_fixture(folder,oldpath)
path(oldpath);
rmdir(folder,'s');
end

function check_engine_rejection(file,work)
fid=fopen(file,'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
restore=onCleanup(@()write_bytes(file,bytes));
altered=bytes;altered(1)=bitxor(altered(1),uint8(1));write_bytes(file,altered);
try
    verify_engine_inputs(fullfile(work,'empirical'));
    error('test:ExpectedFailure','Changed engine input was accepted.');
catch err
    assert(contains(err.message,'Estimation input or code changed:'),'%s',err.message);
end
end

function write_bytes(file,bytes)
fid=fopen(file,'wb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
assert(fwrite(fid,bytes,'uint8')==numel(bytes));
end
