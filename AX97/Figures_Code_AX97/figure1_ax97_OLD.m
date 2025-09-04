
clear all, close all, clc
addpath('/data2/tayannepf/toolbox_movar/scripts_matlab/m_map')
addpath('/data2/tayannepf/toolbox_movar/scripts_matlab/colorbar_yyw/' )
% addpath('/home2/tayannepf/matlab/gsw/')

clr_T = jet; clr_T(1,:) = [.7 .7 .7];
clr_S = colormap_yyw('MPL_rainbow'); clr_S(1,:) = [.7 .7 .7];
clr_V = colormap_yyw('GMT_polar'); 

%%  set colormaps
disp('Plotting section')
%addpath /phodnet/share/mgoes/matlab/mapcolor/ 

clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
% clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
% clr_S = cmocean('haline'); clr_S(1,:) = [.7 .7 .7];
% clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
xyclr = [.4 .4 .4];
lon_lim = [-42 -28];                                   %Define lon lim for map
lat_lim = [-28 -16];                                     %Define lat lim for map

%%  ======= Figures =======
% load topography
x=ncread('ETOPO1_Ice_g_gmt4.grd','x');
y=ncread('ETOPO1_Ice_g_gmt4.grd','y');
z=ncread('ETOPO1_Ice_g_gmt4.grd','z');
kx=find(x>=lon_lim(1)-0.25 & x<=lon_lim(2)+0.25);
ky=find(y>=lat_lim(1)-0.25 & y<=lat_lim(2)+0.25);
lonb=x(kx); latb=y(ky); depth_btm=z(kx,ky)';
depth_btm(depth_btm>0)=NaN;
clear kx ky x y z

%%
load /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/output_ax97.mat
temp_ref=ncread('/data2/Movar_analises/nc_final/transporte_geostrofico.nc','CT');
salt_ref  =ncread('/data2/Movar_analises/nc_final/transporte_geostrofico.nc','SA');
vel_ref  =ncread('/data2/Movar_analises/nc_final/transporte_geostrofico.nc','v_geo');

label_lat = {'42\circW','39\circW','36\circW','33\circW','30\circW'};  %Define horizontal Label ticks

%% Retirando cruzeiros ruins

ind_rm=[24, 28, 42, 76];

lon_xbt(ind_rm)=[];
lat_xbt(ind_rm)=[];



fig1 = myFigSize(1,10,6); clf % papersize (10,6)
subplot(2,3,1);hold off
m_proj('miller','lon',lon_lim,'lat',lat_lim);
[CS,CH]=m_contourf(lonb,latb,depth_btm,[-7000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
colormap(gca,clr_topo(50:end,:))
m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on
for i = 1:length(lon_xbt)
    m_plot(lon_xbt{i},lat_xbt{i},'.r','markersize',4)
end
m_plot(mean(lon_ref,2),lat_ref,'k','linewidth',2)
m_grid('tickdir','in','color',xyclr)
set(gca,'Position',[.01 .56 .34 .36])
ax=m_contfbar(0.31,[.57 .94],CS,CH,'axfrac',.025);
% set(ax,'colormap',clr_topo(50:end,:))
set(ax,'tickdir','in','fontsize',8,'ycolor',xyclr,'xcolor',xyclr)
title(ax,{'m'},'fontweight','normal','fontsize',8);
m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,'(a) AX97')

% histogram # of cruises (month)
subplot(2,3,2); hold off
% [~,m] = datevec(time_avg);
m=meses_ax';
h = histogram(m,0.5:1:12.5); hold on
set(h,'EdgeColor',[0.1294 0.1294 0.1294])
set(gca,'XLim',[0.5 12.5],'YLim',[0 20])
set(gca,'Position',[0.38 0.56 0.23 0.36])
text(0.5,21,'(b) Monthly distribution of cruises')
ylabel('number of cruises')



% mean temperature
subplot(2,3,4); hold off
mtemp = mean(temp_ref,1,'omitnan');
tmp = squeeze(mtemp);
tmp(~isnan(tmp)) = 26;
tmp(isnan(tmp)) = -5;
contourf(lon_ref,-depth_ref,tmp(1:end-1,:)',-6:32:26); hold on
[CS,CH] = contourf(lon_ref,-depth_ref,squeeze(mtemp(1,1:end-1,:))',2:1:26,'LineStyle','none'); hold on
colormap(gca,clr_T);
caxis([0 26])
plot(lon_xbt{82},0,'.k')
ax=m_contfbar(0.82,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[2 26],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
% set(ax,'colormap',clr_T)
title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.06 0.1 0.23 0.35])
ylabel('Depth (m)')
set(gca,'xtick',-42:3:-28,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lon_ref(1),50,'(d) Mean temperature')

xlim([min(lon_ref) max(lon_ref)])
% set(gca,'color',[.7 .7 .7])

% mean salinity
msalt = mean(salt_ref,1,'omitnan');
msalt(msalt<32.5) = 32.5;
tmp = squeeze(msalt);
tmp(~isnan(tmp)) = 37;
tmp(isnan(tmp)) = 27;
subplot(2,3,5); hold off
% set(gca,'color',[.7 .7 .7])

contourf(lon_ref,-depth_ref,tmp(1:end-1,:)',22:15:37); hold on
[CS,CH] = contourf(lon_ref,-depth_ref,squeeze(msalt(1,1:end-1,:))',32:0.2:37,'LineStyle','none'); hold on
colormap(gca,clr_S);
caxis([30 37])
ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[32 37],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
% set(ax,'colormap',clr_S)
title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.38 0.1 0.23 0.35])
set(gca,'xtick',-42:3:-28,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lon_ref(1),50,'(e) Mean salinity')
xlim([min(lon_ref) max(lon_ref)])


% mean velocity
mvel = nanmean(vel_ref,1);
tmp = squeeze(mvel);
tmp(~isnan(tmp)) = 0.35;
tmp(isnan(tmp)) = -0.35;
subplot(2,3,6); hold off
contourf(lon_ref,-depth_ref,tmp',-0.7:1.4:0.7); hold on
colormap(gca,clr_V2);
[CS,CH] = contourf(lon_ref,-depth_ref,squeeze(mvel)',-0.15:0.01:0.15,'LineStyle','none');
contour(lon_ref,-depth_ref,squeeze(mvel)',[0 0],'color',[.7 .7 .7]);
caxis([-0.15 0.15])
ax=m_contfbar(1.2,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
% ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');

set(ax,'tickdir','in','ylim',[-0.15 0.15],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',-0.15:0.05:0.15,'yticklabel',-0.15:0.05:0.15)
% set(ax,'colormap',clr_V2)
title(ax,'m/s','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.7 0.1 0.23 0.35])
set(gca,'xtick',-42:3:-28,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lon_ref(1),50,'(f) Mean velocity')
xlim([min(lon_ref) max(lon_ref)])

%PIE
%CREATE X,Y AXES
% timeYMD = datestr(time_avg','yyyymmdd');
% data1 = str2num(timeYMD(:,5:6))'; %Month numbers
% data2 = str2num(timeYMD(:,1:4))'; %year numbers
data1 = meses_ax; %Month numbers
data2 = anos_ax; %year numbers


%GET DATA RANGE
axis1 = [min(data1):max(data1)]';
axis2 = [min(data2):max(data2)]';

%create month axis if wanted
leg1 = {'jan','feb','mar','apr','may','jun','jul','ago','sep','oct','nov','dec'};

%Add N number of spaces in the middle
N = 3;
C = npf; %number of profiles/transect

%MAKE PIE PLOTS
%1: month, year
%clf
[g,col] = make_pie_subplot_tay(axis2,axis1,[],leg1,data2,data1,C,N,[.62 .51 .345 .415]);
title(col,'# prof')
print -dpng -r400 /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/figure1_ax97_2.png

%%
%Redo auxiliary month dates
[y_xbt,m_xbt,d_xbt] = datevec(time_avg);
Yday = datenum(0,m_xbt,d_xbt);
for i = 1:12
   DayM(i) = mean(datenum(0,i,1):datenum(0,i+1,1));
end


label_GSpos = {'36.5\circN','37\circN','37.5\circN','38\circN','38.5\circN'};
fig2 = myFigSize(2,10,6); clf

subplot(2,2,1); hold off
plot(Yday,GSlat,'.k'); hold on
plot(DayM,mGSlat(:,1),'Color',[0 .45 .75],'LineWidth',1.5)
set(gca,'xlim',[0 366],'ylim',[36.5 38.7],'ytick',36.5:0.5:38.5,'YTickLabel',label_GSpos)
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'position',[0.1 0.56 0.36 0.35])
text(0,38.88,'(a) GS position (seasonal)')

subplot(2,2,2); hold off
plot(Yday,GStransp/1e6,'.k'); hold on
plot(DayM,mGStransp(:,1)/1e6,'Color',[0 .45 .75],'LineWidth',1.5)
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'xlim',[0 366],'ylim',[26 44])
set(gca,'position',[0.56 0.56 0.36 0.35])
text(0,45.5,'(b) GS transport (seasonal)')

subplot(2,2,3);hold off
plot(time_avg,GSlat,'color',[0 .45 .75],'LineWidth',1);hold on
set(gca,'ylim',[36.5 38.7],'xlim',[time_avg(1)-30 time_avg(end)+30])
set(gca,'xtick',datenum(2010:2:2024,1,1),'xticklabel',2010:2:2024,'XTickLabelRotation',0)
set(gca,'ytick',36.5:0.5:38.5,'YTickLabel',label_GSpos)
text(time_avg(1)-30,38.88,'(c) GS position')
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'position',[0.1 0.1 0.36 0.35])

subplot(2,2,4);hold off
plot(time_avg,GStransp/1e6,'color',[0 .45 .75],'LineWidth',1);hold on
ylabel('Transport (Sv)')
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'ylim',[26 44],'xlim',[time_avg(1)-30 time_avg(end)+30])
set(gca,'yaxislocation','left','xaxislocation','bottom')%,'color','none')
set(gca,'xtick',datenum(2010:2:2024,1,1),'xticklabel',2010:2:2024,'XTickLabelRotation',0)
set(gca,'position',[0.56 0.1 0.36 0.35])
text(time_avg(1)-30,1.774,'(d) GS transport')

print -dpng -r400 figure2_ax32.png


%% Hovmoller of dynamic height
fig3 = myFigSize(3,5,8); clf
label_lat = {'34\circN','','36\circN','','38\circN',''};
nt = length(time_avg);
time_label=str2num(datestr(time_avg,'yyyy'));
%THIS USES THE REAL AXIS
%contourf(lat_ref,time_avg,dh0_ref',[0.6:.1:1.8]); hold on
%set(gca,'ytick',datenum(2010:2:2024,1,1),'yticklabel',2010:2:2024,'XTickLabelRotation',0)

%THIS USES A FAKE AXIS
contourf(lat_ref,1:nt,dh0_ref',[0.6:.1:1.8],'edgecolor','none'); hold on
set(gca,'ytick',10:10:nt,'yticklabel',time_label(10:10:nt))
%
colormap(gca,clr_V2);
xlim([33 39])
set(gca,'xtick',34:1:39,'xticklabel',label_lat,'XTickLabelRotation',0)
set(gca,'position',[0.13 0.2 0.7 0.7])
ax = colorbar('position',[.85 .2 .03 .7]);
set(ax,'tickdir','in','fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',0.6:0.2:1.8,'yticklabel',0.6:0.2:1.8)
set(ax,'colormap',clr_V2)
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
title(ax,'m','FontSize',10,'color',xyclr,'FontWeight','normal')
title('Surface Dynamic Height ')
print -dpng -r400 figure3_ax32.png
