function retry_failed_draws(folder,key,B)
% Preserve unsuccessful restart draws separately so a retry can recompute them.
for b=1:B
 file=fullfile(folder,sprintf('draw_%05d.mat',b));if ~isfile(file),continue;end
 L=load(file);assert(strcmp(L.key,key),'Bootstrap checkpoint mismatch.');
 values=[L.xRow(1),L.vRow(1),L.aRow(1),L.averageAccelerationRow];
 if ~all(isfinite(values))
  archive=fullfile(folder,'failed_draws');if ~isfolder(archive),mkdir(archive);end
  [~,name]=fileparts(file);backup=[tempname(archive) '_' name '.mat'];movefile(file,backup);
 end
end
end
