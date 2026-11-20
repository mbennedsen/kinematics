function validate_smooth_design
[root,C,S]=prepare_smooth_coverage('smoke');
for na=[1 12 20 32]
 ids=[];for a=1:na,ids=[ids a:na:400];end %#ok<AGROW>
 assert(isequal(sort(ids),1:400),'Allocation has gaps or overlaps.');
end
[a,t,i]=simulate_smooth_record(C,1);[b,u,j]=simulate_smooth_record(C,1);[~,v,k]=simulate_smooth_record(C,400);
assert(isequaln(a,b)&&isequaln(t,u)&&isequaln(i,j));
assert(isequaln(t,v)&&isequaln(i.alpha,k.alpha)&&~isequaln(i.noise,k.noise));
x=simulate_temperature_record(C,C.selected(3),1);
assert(isequaln(i.alpha,x.alpha)&&isequaln(a.NASA(C.sampleIndex),x.y)&&isequaln(t,x.truth));
assert(i.rebaseShift==0 && max(abs(a.NASA(C.sampleIndex)-i.alpha(1,:)'-C.sample.W*C.est.paramsHat.beta-i.noise))<1e-12);
truthTable=readtable(fullfile(fileparts(root),'data','Simulation','DGPTruth.csv'));row=strcmp(string(truthTable.DGP),'NASA_smoothed');
assert(sum(row)==1 && max(abs(table2array(truthTable(row,2:5))-t))<1e-10);
dt=1/12;z=i.alpha;da=diff(z(3,:));
assert(max(abs(diff(z(2,:))-dt*(z(3,1:end-1)+da/2)))<1e-14);
assert(max(abs(diff(z(1,:))-dt*z(2,1:end-1)-dt^2*(z(3,1:end-1)/2+da/6)))<1e-14);
assert(abs(t(4)-100*(z(2,end)-z(2,C.index0))/(C.endTime-C.actualT0))<1e-12);
P=table(C.sample.t,z(1,:)',10*z(2,:)',100*z(3,:)','VariableNames',{'t','x','v','a'});
writetable(P,fullfile(root,'validation','FixedDGP3Path.csv'));
writetable(array2table(t,'VariableNames',{'x_T','v_T','a_T','average_a_past_decade'}),fullfile(root,'validation','FixedDGP3Truth.csv'));
fprintf('PASS: exact accuracy DGP 3, fixed integrated path, no rebasing, independent noise, paired average, and gap-free worker-independent allocation.\n');
fprintf('Truth: x=%.12f v=%.12f a=%.12f average a=%.12f.\n',t);
end
