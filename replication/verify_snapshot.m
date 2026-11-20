function verify_snapshot
% VERIFY_SNAPSHOT Check the frozen data, saved results and supplied MATLAB code.
root=fileparts(fileparts(mfilename('fullpath')));
T=readtable(fullfile(root,'data','SNAPSHOT.csv'),'TextType','string');
for k=1:height(T)
 file=fullfile(root,T.RelativePath(k));assert(isfile(file),'Missing file: %s',file);
 f=fopen(file,'rb');assert(f>0);c=onCleanup(@()fclose(f));bytes=fread(f,Inf,'*uint8');clear c
 assert(numel(bytes)==T.Bytes(k),'Size mismatch: %s',file);
 md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
 hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
 assert(strcmp(hash,T.SHA256(k)),'Changed snapshot file: %s',file);
end
fprintf('SNAPSHOT PASS: %d data, result and code files verified.\n',height(T));
end
