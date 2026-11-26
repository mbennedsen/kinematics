function digest=file_sha256(file)
% FILE_SHA256 Hash raw file bytes using native MATLAB and bounded memory.
% The private SHA-256 core is from DataHash (see private/DataHash_NOTICE.md).
f=fopen(file,'rb');
assert(f>=0,'replication:FileOpen','Cannot open %s',file);
cleanup=onCleanup(@() fclose(f));
context=sha256Initialize();
tail=zeros(0,1,'uint8');
while true
    bytes=fread(f,1024*1024,'*uint8');
    if isempty(bytes),break;end
    context.totalLength=context.totalLength+uint64(numel(bytes));
    bytes=[tail;bytes];
    count=64*floor(numel(bytes)/64);
    if count>0
        context.state=sha256Compress(context.state,reshape(bytes(1:count),64,[]));
    end
    tail=bytes(count+1:end);
end
[message,number]=ferror(f);
assert(number==0 || feof(f),'replication:FileRead','Cannot read %s: %s',file,message);
value=sha256Finalize(context.state,tail,context.totalLength);
digest=lower(reshape(dec2hex(value,2).',1,[]));
end
