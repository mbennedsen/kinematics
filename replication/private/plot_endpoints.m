function plot_endpoints(f,fit,methods,labels)
% Group records by estimator as in main Figure 2.
fit=fit(fit.Specification=="monthly",:);
records=["NASA","NOAA","HadCRU","Berk","ERA5","Mean"];
names={'NASA','NOAA','HadCRUT5','Berkeley Earth','ERA5','Arithmetic mean'};
colors=[68 119 170;238 119 51;34 136 51;204 102 119;102 170 170;34 34 34]/255;
shades=[234 242 250;252 240 232;247 243 223;242 234 247]/255;
markers={'o','s','^','d','v','p'};
objects=["x","v","a","average_a"];
titles={'a   Warming level','b   Warming rate','c   Warming acceleration',{'d   Average acceleration','over past decade'}};
units={'°C','°C decade^{-1}','°C decade^{-2}','°C decade^{-2}'};
ranges=[1.25 1.60;-.04 .80;-.06 .68;-.04 .50];
set(f,'Position',[40 40 1350 800]); layout=tiledlayout(f,1,4,'TileSpacing','compact','Padding','compact');
for k=1:4
    ax=nexttile(layout);hold(ax,'on');lim=ranges(k,:);
    for m=1:4
        top=(m-1)*7.8;
        patch(ax,lim([1 2 2 1]),top+[-.4 -.4 7 7],shades(m,:),'EdgeColor','none');
        text(ax,lim(2)-.025*diff(lim),top+.2,labels{m},'HorizontalAlignment','right','FontWeight','bold');
        for j=1:6
            z=fit(fit.Method==methods{m} & fit.Dataset==records(j) & fit.Object==objects(k),:);
            assert(height(z)==1); y=top+.5+j;
            if isfinite(z.CILow) && isfinite(z.CIHigh)
                plot(ax,[z.CILow z.CIHigh],[y y],'-','Color',colors(j,:));
                plot(ax,[z.CILow z.CIHigh],[y y],'|','Color',colors(j,:));
            end
            plot(ax,z.Estimate,y,markers{j},'Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',5);
        end
    end
    ax.YDir='reverse';ylim(ax,[-.6 30.7]);xlim(ax,lim);
    positions=(1.5:6.5)' + (0:3)*7.8;yticks(ax,positions(:)');
    if k==1,yticklabels(ax,repmat(names,1,4));else,yticklabels(ax,{});end
    if k==1,xline(ax,1.5,'--');else,xline(ax,0,'--');end
    title(ax,titles{k});xlabel(ax,units{k});box(ax,'off');
end
end
