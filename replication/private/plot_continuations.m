function plot_continuations(f,empirical,colors,supplement)
% Use the same full-precision paths, crossing formula and endpoints as Tables S2/S3.
methods = ["ss","kernel","quadratic","changepoint"];
names = ["Semiparametric","Nonparametric","Quadratic","Changepoint"];
half = continuation_results(empirical,1.5);
two = continuation_results(empirical,2);
set(f,'Position',[50 50 990 490]);
if supplement
    main = axes(f,'Position',[.08 .17 .54 .66]);
    zoom = axes(f,'Position',[.70 .17 .27 .66]);
else
    main = axes(f,'Position',[.09 .14 .88 .65]);
    backdrop = axes(f,'Position',[.60 .155 .37 .33],'Color','white','Box','on');
    set(backdrop,'XTick',[],'YTick',[],'XColor',[.75 .75 .75],'YColor',[.75 .75 .75]);
    zoom = axes(f,'Position',[.65 .22 .27 .20],'Color','white');
end
hold(main,'on'); hold(zoom,'on');
xlim(main,[2014.5 2048.2]); ylim(main,[.93 2.14]);
xlim(zoom,[2026 2031.65]); ylim(zoom,[1.431 1.569]);
if supplement, legendPos = [.08 .89 .89 .08]; else, legendPos = [.09 .90 .87 .07]; end
for ax = [main zoom]
    grid(ax,'on'); box(ax,'off'); ax.GridAlpha=.15;
end
ylabel(main,'Warming above pre-industrial (°C)'); xlabel(main,'Year');
title(main,'a   Continuations to 2 °C'); title(zoom,'b   Crossings of 1.5 °C');
yline(main,2,'Color',[.5 .5 .5]); xline(main,two.Endpoint(1),':');
rectangle(main,'Position',[2026 1.431 5.65 .138],'EdgeColor',[.55 .55 .55]);
yline(zoom,1.5,'Color',[.5 .5 .5]); xticks(zoom,[2026 2028 2030]); yticks(zoom,[1.45 1.5 1.55]);
handles = gobjects(1,4);
dy = [.065 .035 .057 .057]; dx = [.2 -.55 0 0];
dyZoom = [.022 -.019 .021 -.019]; dxZoom = [-.08 -.09 0 .22];
for j = 1:4
    method = methods(j);
    if method=="ss", file=fullfile(empirical,'validated_monthly','Mean_ss_paths.csv');
    else, file=fullfile(empirical,'reporting',"Mean_"+method+"_paths.csv"); end
    p = readtable(file);
    selected=two.Method==method & ~endsWith(two.Continuation,"constant rate");
    model=two(selected,:); small=half(selected,:);
    assert(height(model)==1 && max(abs([p.t(end),p.x(end),p.v(end),p.a(end)]- ...
        [model.Endpoint,model.Level,model.Rate,model.Acceleration]))<1e-7,'Endpoint/path mismatch.');
    col=colors(j,:);
    for ax=[main zoom]
        h=plot(ax,p.t,p.x,'Color',col,'LineWidth',1.6);
        if ax==main,handles(j)=h;end
        curve(ax,model,'--',col);
        if supplement && j<=2
            flat=two(two.Method==method & endsWith(two.Continuation,"constant rate"),:);
            curve(ax,flat,':',col);
            if ax==main
                plot(ax,flat.CrossingYear,2,'o','Color',col,'MarkerFaceColor','w');
                text(ax,flat.CrossingYear,2-.055,sprintf('%.1f',flat.CrossingYear),'Color',col,'HorizontalAlignment','center');
            end
        end
    end
    for ax=[main zoom],plot(ax,small.CrossingYear,1.5,'o','MarkerSize',4,'Color',col,'MarkerFaceColor',col);end
    plot(main,model.CrossingYear,2,'o','MarkerSize',5,'Color',col,'MarkerFaceColor',col);
    text(main,model.CrossingYear+dx(j),2+dy(j),sprintf('%.1f',model.CrossingYear),'Color',col,'HorizontalAlignment','center');
    text(zoom,small.CrossingYear+dxZoom(j),1.5+dyZoom(j),sprintf('%.1f',small.CrossingYear),'Color',col,'FontSize',8,'HorizontalAlignment','center');
end
leg=legend(main,handles,cellstr(names),'Orientation','horizontal','Box','off');
leg.Units='normalized'; leg.Position=legendPos;
end

function curve(ax,model,lineStyle,col)
years=linspace(model.Endpoint,model.CrossingYear,300);
h=(years-model.Endpoint)/10;
levels=model.Level+model.Rate*h+0.5*model.Acceleration*h.^2;
plot(ax,years,levels,'LineStyle',lineStyle,'Color',col,'LineWidth',1.5);
end
