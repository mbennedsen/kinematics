function merge_rows(file,T,keys)
if isfile(file)
 old=readtable(file,'TextType','string');remove=ismember(old(:,keys),T(:,keys));
 T=[old(~remove,:);T];
end
writetable(T,file);
end
