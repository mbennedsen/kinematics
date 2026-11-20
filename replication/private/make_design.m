function make_design(dest)
% Regenerate deterministic latent paths and their truths from frozen calibration.
[~,~,C]=prepare_accuracy;
out=fullfile(dest,'simulations','accuracy');if ~isfolder(out),mkdir(out);end
truth=table();
for k=1:numel(C.selected)
 d=C.selected(k);s=simulate_temperature_record(C,d,1);a=s.alpha;
 T=array2table([C.sample.t(:),a(1,:)',10*a(2,:)',100*a(3,:)'],'VariableNames',{'t','x','v','a'});
 writetable(T,fullfile(out,['DGP_' C.scenarios{d} '.csv']));
 r=table(string(C.scenarios{d}),s.truth(1),s.truth(2),s.truth(3),s.truth(4), ...
  'VariableNames',{'DGP','x_T','v_T','a_T','average_a_past_decade'});truth=[truth;r]; %#ok<AGROW>
end
writetable(truth,fullfile(out,'DGPTruth.csv'));
end
