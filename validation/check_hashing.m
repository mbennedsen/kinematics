function check_hashing
% CHECK_HASHING Verify raw-byte SHA-256, including padding/chunk boundaries.
% Run after addpath('validation'); no Java or external tools are used.
root=fileparts(fileparts(mfilename('fullpath')));
old=path;restore=onCleanup(@()path(old));
addpath(fullfile(root,'replication','engines','methods'),'-begin');
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s'));
file=fullfile(folder,'binary ''quoted'' space æøå.bin');
vectors=readtable(fullfile(root,'validation','sha256_vectors.csv'),'TextType','string');
started=tic;
for k=1:height(vectors)
    switch vectors.Kind(k)
        case "empty",bytes=zeros(0,1,'uint8');
        case "abc",bytes=uint8('abc');
        case "nist56",bytes=uint8('abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq');
        case "nist112"
            bytes=uint8(['abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn' ...
                'ohijklmnopijklmnopqjklmnopqrklmnopqrslmnopqrsmnopqrstnopqrstu']);
        case "million_a",bytes=repmat(uint8('a'),1,1000000);
        case "binary",bytes=uint8(mod(37*(0:vectors.Bytes(k)-1)+11,256));
        otherwise,error('Unknown test vector.');
    end
    assert(numel(bytes)==vectors.Bytes(k));
    write_bytes(file,bytes);
    assert(strcmp(file_sha256(file),vectors.SHA256(k)), ...
        'SHA-256 mismatch for vector %d (%s, %d bytes).',k,vectors.Kind(k),numel(bytes));
end
% Equal-length changes must change the digest; file size alone cannot suffice.
before=file_sha256(file);bytes(1)=bitxor(bytes(1),uint8(1));write_bytes(file,bytes);
assert(~strcmp(file_sha256(file),before),'Same-length corruption was missed.');
try
    file_sha256(fullfile(folder,'missing.bin'));
    error('test:ExpectedFailure','Missing file was accepted.');
catch err
    assert(strcmp(err.identifier,'replication:FileOpen'),'%s',err.message);
end
fprintf('HASHING PASS: %d vectors, same-size corruption and missing file (%.2f s).\n', ...
    height(vectors),toc(started));
end

function write_bytes(file,bytes)
fid=fopen(file,'wb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
assert(fwrite(fid,bytes,'uint8')==numel(bytes));
end
