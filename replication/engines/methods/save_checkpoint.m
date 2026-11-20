function save_checkpoint(file,R)
folder=fileparts(file);if ~isfolder(folder),mkdir(folder);end
tmp=[tempname(folder) '.mat'];save(tmp,'R','-v7');movefile(tmp,file,'f');
end
