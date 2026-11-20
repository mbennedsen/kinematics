function collect_windows(rob,emp,dest)
R=table();names={'NASA','NOAA','HadCRU','Berk','ERA5'};methods={'ss','kernel','quadratic','changepoint'};
for d=1:5
 name=names{d};s=fullfile(emp,'results','monthly',name);
 W=readtable(fullfile(s,'ss_B1999','AverageAcceleration_Windows.csv'),'TextType','string');
 K=readtable(fullfile(rob,['AlternativeWindows_Kernel_' name '.csv']),'TextType','string');
 Q=readtable(fullfile(s,'quadratic_B1999','Summary.csv'),'TextType','string');Q=Q(Q.Object=="a",:);
 C=readtable(fullfile(s,'changepoint_B1999','AverageAcceleration_Windows.csv'),'TextType','string');
 for year=[2010 2012 2014 2015 2016 2018 2020]
  for m=1:4
   switch m
    case 1,z=W(W.Window=="since_"+year,:);
    case 2,z=K(K.RequestedT0==year,:);
    case 3,z=Q;
    case 4,z=C(C.Window=="since_"+year,:);
   end
   assert(height(z)==1);
   Z=table(string(name),string(methods{m}),"average_a",year,z.Estimate,z.CILow,z.CIHigh, ...
    'VariableNames',{'Dataset','Method','Object','Year','Estimate','CILow','CIHigh'});R=[R;Z]; %#ok<AGROW>
  end
 end
end
writetable(R,fullfile(dest,'empirical','reporting','AlternativeWindows_AllMethods.csv'));
end
