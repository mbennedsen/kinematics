function f=paper_figure(item,source,data,visible)
rd=@(p)readtable(p,'TextType','string');e=fullfile(source,'empirical');
colors=[23 104 172;213 94 0;179 138 0;122 31 162]/255;
methods={'ss','kernel','quadratic','changepoint'};labels={'Semiparametric','Nonparametric','Quadratic','Changepoint'};
f=figure('Visible',visible,'Color','w','Name',item,'Position',[50 50 1100 650]);
set(f,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',11);
if strcmp(item,'figure1')||ismember(item,{'FigureS3','FigureS4','FigureS5','FigureS6','FigureS7','FigureS8'})
 recent=strcmp(item,'figure1');name='Mean';
 if ~recent,ds={'Mean','NASA','NOAA','HadCRU','Berk','ERA5'};name=ds{str2double(extractAfter(item,'FigureS'))-2};end
 z=rd(fullfile(e,'validated_monthly',[name '_ss_paths.csv']));
 if recent,set(f,'Position',[50 50 1250 430]);tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');else,set(f,'Position',[50 50 850 900]);tl=tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');end
 objs={'x','v','a'};units={'Temperature anomaly (\circC)','Rate (\circC decade^{-1})','Acceleration (\circC decade^{-2})'};titles={'Warming level','Warming rate','Warming acceleration'};limits={[.35 1.85],[0 .75],[-.15 .55]};
 handles=gobjects(0);legendLabels={};
 for j=1:3
  ax=nexttile(tl);hold(ax,'on');o=objs{j};
  if j==1,h=plot(ax,z.t,z.y,'Color',[.68 .68 .68],'LineWidth',.6);handles(end+1)=h;legendLabels{end+1}='Observed';end
  fill(ax,[z.t;flipud(z.t)],[z.([o 'Low']);flipud(z.([o 'High']))],colors(1,:),'FaceAlpha',.14,'EdgeColor','none','HandleVisibility','off');
  for m=1:4
   if j==3&&m==4,continue;end
   if m==1,p=z;else,p=rd(fullfile(e,'reporting',[name '_' methods{m} '_paths.csv']));end
   style='-';if m==4,style='--';end
   if recent && strcmp(o,'v') && strcmp(methods{m},'changepoint')
    h=stairs(ax,p.t,p.v,'Color',colors(m,:),'LineWidth',1.5,'LineStyle',style);
   else
    h=plot(ax,p.t,p.(o),'Color',colors(m,:),'LineWidth',1.5,'LineStyle',style);
   end
   if j==1,handles(end+1)=h;legendLabels{end+1}=labels{m};end
  end
  ylabel(ax,units{j});grid(ax,'on');ax.GridAlpha=.16;box(ax,'off');
  if recent,xlim(ax,[2000 2026.8]);ylim(ax,limits{j});xticks(ax,[2000 2010 2020 2026]);xlabel(ax,'Year');title(ax,[char(96+j) '  ' titles{j}]);
  else,xlim(ax,[z.t(1) 2026.8]);xline(ax,z.t(end)-10,':','Color',[.5 .5 .5]);title(ax,[char(96+j) '  ' titles{j}]);end
  if j==3,yline(ax,0,'Color',[.5 .5 .5],'HandleVisibility','off');end
 end
 lg=legend(handles,legendLabels,'Orientation','horizontal','Box','off');lg.Layout.Tile='north';
 if ~recent,xlabel(ax,'Year');title(tl,name,'Interpreter','none');end
elseif strcmp(item,'figure2')
 plot_endpoints(f,rd(fullfile(e,'Empirical_Results.csv')),methods,labels);
elseif ismember(item,{'figure3','FigureS9'})
 plot_continuations(f,e,colors,strcmp(item,'FigureS9'));
elseif strcmp(item,'FigureS1')
 r=rd(fullfile(data,'Monthly_HadISST.csv'));ds={'NASA','NOAA','HadCRU','Berk','ERA5'};ax=axes(f);hold(ax,'on');
 for k=1:5,z=r.(ds{k});b=mean(z(r.t>=1991&r.t<2021),'omitnan');plot(ax,r.t,z-b+.88,'LineWidth',.7,'DisplayName',ds{k});end
 xlim(ax,[1880 2027]);xlabel(ax,'Year');ylabel(ax,'Temperature anomaly (\circC)');legend(ax,'Location','northwest','Orientation','horizontal');grid(ax,'on');box(ax,'off');
elseif strcmp(item,'FigureS2')
 r=rd(fullfile(data,'Monthly_HadISST.csv'));tl=tiledlayout(f,3,1,'TileSpacing','compact');cols={'NINO34','volc2','ssn'};names={'ENSO','Volcanic aerosols','Sunspots'};
 for k=1:3,ax=nexttile(tl);z=r.(cols{k});plot(ax,r.t,(z-mean(z,'omitnan'))/std(z,'omitnan'),'Color',colors(1,:));xlim(ax,[1880 2027]);ylabel(ax,names{k});grid(ax,'on');box(ax,'off');end
 xlabel(ax,'Year');
elseif strcmp(item,'FigureS10')
 tl=tiledlayout(f,2,2,'TileSpacing','compact');ds={'zero','ramp','NASA_smoothed','NASA_quadratic'};names={'Zero acceleration','Increasing acceleration','NASA-smoothed acceleration','Constant positive acceleration'};
 for k=1:4
  z=rd(fullfile(source,'simulations','accuracy',['DGP_' ds{k} '.csv']));ax=nexttile(tl);hold(ax,'on');plot(ax,z.t,z.a,'k','LineWidth',1.7);yline(ax,0,':');xline(ax,z.t(end)-10,'--');
  xlim(ax,[1881 z.t(end)]);ylim(ax,[-.25 .32]);xticks(ax,[1881 1920 1960 2000 2026]);ylabel(ax,'\circC decade^{-2}');title(ax,[char(96+k) '  ' names{k}]);grid(ax,'on');box(ax,'off');if k>2,xlabel(ax,'Year');end
 end
else,error('Unknown figure: %s',item);
end
end
