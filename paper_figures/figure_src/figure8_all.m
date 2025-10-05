%Figure 8

close all
clear
%%

transect = {'AX32','AX97','PX30','PX36','PX40'}
legs = {'a)','b)','c)','d)','e)'}
currs = {'Gulf Stream','Brazil Current','East Australian Current',...
    'Antarctic Circumpolar Current','Kuroshio Current'}
fig3 = myFigSize(3,10,6); clf % papersize 

for t = 1:length(transect)
    file = [transect{t},'/',transect{t},'_gridded_new.nc']
    if ~exist(file,'file')
        file = [transect{t},'/',transect{t},'_gridded.nc']
    end

    GStransp_alt = ncread(file,'altimetric_ssh_gradient');
    GStransp = ncread(file,'geostrophic_transport');
    time_avg = ncread(file,'time')+datetime(1950,1,1);

    %Calculate Correlation between altimetry and XBT

    inan = ~isnan(GStransp_alt+GStransp);%inan=inan(1:end-1);
    C = corr(GStransp_alt(inan),GStransp(inan))

    %%
    subplot(3,2,t)
    plot(time_avg(inan),GStransp(inan)/1e6,'color',[0 .45 .75],'LineWidth',1);
    hold on
    plot(time_avg(inan),GStransp(inan)/1e6,'.','color',[0 .45 .75],'markersize',6);
    ylabel('Transport (Sv)')
    %set(gca,'position',[0.1 0.1 0.8 0.8])
    h1=gca;h2=axes('position',get(h1,'position'));
    plot(time_avg(inan),GStransp_alt(inan),'color','r','LineWidth',1);hold on
    ylabel('cross-front \DeltaSSH (m)')
    set(h2,'yaxislocation','right','xaxislocation','top','color','none')
    set(h2,'xcolor','k','ycolor','r','box','off','xticklabel',[])
    set(h2,'xlim',[time_avg(1)-30 time_avg(end)+30])
    set(h1,'yaxislocation','left','xaxislocation','bottom','color','none')
    set(h1,'xcolor',[.3 .3 .3],'ycolor',[0 0.45 0.75],'box','off')
    set(h1,'xlim',[time_avg(1)-30 time_avg(end)+30])
    y_lim = get(h2,'ylim');
    textc = sprintf('%s %s %s',legs{t},transect{t},currs{t})
    text(time_avg(1)-30,y_lim(2)+diff(y_lim)/10,textc)
    text(time_avg(1)+90,y_lim(1)+diff(y_lim)/10,sprintf('R = %4.2f',C))
end
        exportgraphics(gcf, 'Figure8.png', 'Resolution', 400);



