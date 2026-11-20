function [raw,truth,info]=simulate_smooth_record(C,rep)
% Reuse the exact accuracy generator. The third selected calibration is the smoothed NASA path.
d=C.selected(3);assert(d==4 && strcmp(C.scenarios{d},'NASA_smoothed'));
sim=simulate_temperature_record(C,d,rep);
raw=C.raw;raw.NASA(:)=NaN;raw.NASA(C.sampleIndex)=sim.y;
truth=sim.truth;
info=struct('alpha',sim.alpha,'noise',sim.noise,'rebaseShift',0, ...
 'seed',C.cfg.seed+1000000+rep,'manuscriptDGP',3,'calibrationIndex',d, ...
 'dgpMode','fixedIntegratedNASA_DGP3','averageStart',C.actualT0,'endTime',C.endTime);
end
