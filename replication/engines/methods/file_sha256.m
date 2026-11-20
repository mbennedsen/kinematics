function digest=file_sha256(file)
% Hash raw bytes to reject altered saved outer samples on resume.
f=fopen(file,'rb'); assert(f>=0,'Cannot open %s',file); cleanup=onCleanup(@() fclose(f));
bytes=fread(f,Inf,'*uint8'); md=java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes); value=typecast(md.digest(),'uint8');
digest=lower(reshape(dec2hex(value,2)',1,[]));
end
