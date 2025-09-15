%addpath(genpath('/datasets/work/soop-xbt/work/cowley/matlab_localtools/toolbox/local/gsw_toolbox/'))
%addpath /datasets/work/soop-xbt/work/UOT/programs/MATLAB/seawater_ver3_3/
%Output: GStransp, GSspd_z, GStransp_z; 
% Extra Output: GStransp_alt, GSlat_alt, GSspd_alt
close all
clear

addpath([pwd '/helper_functions']);
addpath([pwd '/helper_functions/colorbar_yyw/']);
addpath([pwd '/data']);
addpath([pwd '/results']);
addpath([pwd '/m_map']);

%%
load 'output_p40.mat';
vel_ref = permute(vel_ref,[2 1 3]);

% load ssh for the reference section (matching the XBT grid)
%load /phodnet/share/dong/gulfstream/oleander/
load ssh4px40.mat %%%% ssh_ref_px30

%%
disp('transp calculated')
calc_transp

%save Output.mat -append GStransp_alt GSlat_alt GSspd_alt GSlat GSlon GStransp

%Calculate Correlation between altimetry and XBT
inan = ~isnan(GStransp_alt+GStransp);
C = corr(GStransp_alt(inan),GStransp(inan))

%%
fig3 = myFigSize(3,7,3); clf % papersize (10,10)

% time_avg is a bit weird timeYMD is perfect to convert that
time_avg = datenum(timeYMD);

plot(time_avg,GStransp/1e6,'color',[0 .45 .75],'LineWidth',1);
hold on
plot(time_avg,GStransp/1e6,'.','color',[0 .45 .75],'markersize',6);
ylabel('Transport (Sv)')
set(gca,'position',[0.1 0.1 0.8 0.8])
h1=gca;h2=axes('position',get(h1,'position'));
plot(time_avg,GStransp_alt,'color','r','LineWidth',1);hold on; 
ylabel('cross-front \DeltaSSH (m)')
set(h2,'yaxislocation','right','xaxislocation','top','color','none')
set(h2,'xcolor','k','ycolor','r','box','off','xticklabel',[])
%set(h2,'ylim',[.7 1.4], ...
    set(h2,'xlim',[time_avg(1)-30 time_avg(end)+30])
set(h2,'xtick',datenum(2012:2:2024,1,1),'xticklabel',[])
set(h1,'yaxislocation','left','xaxislocation','bottom','color','none')
set(h1,'xcolor','k','ycolor',[0 0.45 0.75],'box','off')
set(h1,'xtick',datenum(2012:2:2024,1,1),'xticklabel',2012:2:2024)
%set(h1,'ylim',[5 55], ...
    set(h1,'xlim',[time_avg(1)-30 time_avg(end)+30])
%text(double(time_avg(35)),-1.4,'XBT','Color',[0 .45 .75])
%text(double(time_avg(35)),-1.5,'Altimeter','Color','r')
%text(double(time_avg(1)-30),-0.1,'Kuroshio transport')
text(double(time_avg(1)+90),1.4,sprintf('R = %4.2f',C))

print -dpng -r600 figure2.PX40_GS.png


