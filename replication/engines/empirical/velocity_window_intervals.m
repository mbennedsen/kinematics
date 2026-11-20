function tableOut=velocity_window_intervals(E,name,folder)
% All windows use paired errors from the SAME bootstrap paths; no refitting.
t=E.sampleData.t;v=10*E.pointEstimate.smoother.alphaSmooth(2,:)';
B=E.bootstrap.pathBands;e=B.velocityErrors;valid=B.validPathDraws;
starts=[t(end-120),2010,2012,2014,2015,2016,2018,2020];rows=table();
for j=1:numel(starts)
 i=find(t>=starts(j)-1e-9,1);elapsed=(t(end)-t(i))/10;
 point=(v(end)-v(i))/elapsed;errors=(e(end,valid)-e(i,valid))/elapsed;
 ci=point-prctile(errors,[95 5]);
 if j==1,label="past_decade";assert(abs(elapsed-1)<1e-8);assert(max(abs(ci-[B.low(end,4),B.high(end,4)]))<1e-10);
 else,label="since_"+string(starts(j));end
 rows=[rows;table(string(name),label,t(i),t(end),point,ci(1),ci(2),sum(valid), ...
 'VariableNames',{'Dataset','Window','ActualT0','EndTime','Estimate','CILow','CIHigh','ValidDraws'})]; %#ok<AGROW>
end
tableOut=rows;writetable(rows,fullfile(folder,'AverageAcceleration_Windows.csv'));
end
