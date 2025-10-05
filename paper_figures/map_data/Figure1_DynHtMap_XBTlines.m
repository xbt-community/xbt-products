% program to map the Argo dynamic height for globe and overlay XBT
% transects
% load('/Users/jsprintall/Documents/JanetsDocuments/indianocean/papers-ours/Ummenhofer&Hood_InterbasinConnections/Figures/crameri/CrameriColourMaps4.0.mat')
addpath /phodnet/share/mgoes/matlab/mapcolor/colors/
addpath /phodnet/share/mgoes/matlab/mapcolor/
addpath /home/mgoes/matlab/mystuff/
addpath /phodnet/share/mgoes/matlab/m_map/
clear all 
close all

%mainpath = '/Users/jsprintall/Documents/JanetsDocuments/indianocean/papers-ours/Ummenhofer&Hood_InterbasinConnections/Figures/';
mainpath = '/WSmounts/phoddat/share/mgoes/KANDAGA/XBT/AX32/'
% % load([mainpath 'gshhs_highres_shoreline.mat'])
% % gshhs_lat = lat;
% % gshhs_long = lon;

%FROM ETOPO1
topo_long=ncread([mainpath 'ETOPO1_Ice_g_gmt4.grd'],'x');
topo_lat=ncread([mainpath 'ETOPO1_Ice_g_gmt4.grd'],'y');
topo_bath=ncread([mainpath 'ETOPO1_Ice_g_gmt4.grd'],'z');
topo_long = topo_long(1:5:end);
topo_lat = topo_lat(1:5:end);
topo_bath = topo_bath(1:5:end,1:5:end);

% kx=find(x>=lon_lim(1)-0.25 & x<=lon_lim(2)+0.25);
% ky=find(y>=lat_lim(1)-0.25 & y<=lat_lim(2)+0.25);
% lonb=x(kx); latb=y(ky); depth_btm=z(kx,ky)';
% depth_btm(depth_btm>0)=NaN;
% clear kx ky x y z

% addpath Data\Bathy
%topo_long = ncread([mainpath 'etopo6.nc'],'LON6');
%topo_lat = ncread([mainpath 'etopo6.nc'],'LAT6101_1700');
%topo_bath = ncread([mainpath 'etopo6.nc'],'BATH6');

%Shift topo data to be between 0-360
idx_above = find(topo_long>=360);
idx_below = find(topo_long<360);
topo_long2 = [topo_long(idx_above)-360; topo_long(idx_below)];
topo_bath2 = [topo_bath(idx_above,:); topo_bath(idx_below,:)];
clear topo_bath topo_long

load(['argomeanDH_02000db.mat']);
% % % %convert longitude to be in the same range as the topo world map data
% % % tsubxlon_i(find(tsubxlon_i < min(topo_long2))) = 360+tsubxlon_i(find(tsubxlon_i < min(topo_long2)));
% % % argo_long =  tsubxlon_i;
% % % argo_lat =  tsubxlat_i;

% now wrap around
%
Dtopo_long = [topo_long2; topo_long2+360];
Dtopo_bath = [topo_bath2; topo_bath2];
Dargo_LO = [argo_LO; argo_LO+360];
Dargo_DH = [argoDH02000; argoDH02000];
clear topo_bath2 topo_long2

lon1 = 10;
      lon2 = 400;
      lat1 = -65;
      lat2 = 70;

%%
 fig3 = myFigSize(5,9,4); clf % papersize (10,10)
%     figure(1);
    set(gcf,'DefaultAxesFontSize',12)
%    set(gcf,'DefaultAxesFontWeight','bold');
    set(gcf,'DefaultAxesLineWidth',2.);
    set(gcf,'DefaultTextFontSize',12)
    set(gcf,'DefaultLineLineWidth',2);


%             f1.WindowState = 'maximized';
hold on
% hold on


m_proj('mercator','long',[lon1 lon2],'lat',[lat1 lat2]);
%m_contour(topo_long2,topo_lat,topo_bath2',[-1000 -1000],'LineColor',rgb('grey'),'LineWidth',1)

hold on

%       xlim([0 160])
% 
%       ylim([-55 30])
      
%       daspect([1 1 1])
      xlabel('Longitude')
      ylabel('Latitude')
m_contour(Dtopo_long,topo_lat,Dtopo_bath',[0 0],'LineColor',rgb('black'),'LineWidth',1)
 m_contour(Dtopo_long,topo_lat,Dtopo_bath',[0 10000]);
 %colormap('gray');

% edgecolor removes contour lines 
m_contourf(Dargo_LO,argo_LA,Dargo_DH',[5:1:35],'LineColor','k','edgecolor','none')
m_contour(Dargo_LO,argo_LA,Dargo_DH',[5:5:35],'LineColor','k')
caxis([5 35])
      colormap(flipud(crameri('roma')))
      colorbar('eastoutside');
ylabel(colorbar,'Dynamic Height [m^2 s^{-2}]');
hold on

m_contour(Dtopo_long,topo_lat,Dtopo_bath',[-1000 -1000],'LineColor',rgb('black'),'LineWidth',0.25)

% now plot on the XBT transects
T1 = importdata('p05stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1, 0.25,0.75])
T1 = importdata('po9stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1, 0.25,0.75])
T1 = importdata('p13stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('p28stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('p31stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('s37stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('i21stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('p37stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('p22stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('p34stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',0.8,'Color',[1, 0.25, 0.75])

load ax_lines.mat
oi = load('../../AX25/Lat_lon_time_ax25.mat')
lon_ax25 = oi.LON(:,5);% nanmean(oi.LON,2);
lat_ax25 = oi.LAT(:,5);%nanmean(oi.LAT,2);

m_plot(360+lon_ax1,lat_ax1,'Linewidth',0.8,'Color',[1,0.25,0.75])
m_plot(360+lon_ax10,lat_ax10,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax18,lat_ax18,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax1a,lat_ax1a,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax1b,lat_ax1b,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax1c,lat_ax1c,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax2,lat_ax2,'Linewidth',0.8,'Color',[1,0.25, 0.75])
% m_plot(360+lon_ax20,lat_ax20,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax25,lat_ax25,'Linewidth',0.8,'Color',[1,0.25, 0.75])
% m_plot(360+lon_ax4,lat_ax4,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax7,lat_ax7,'Linewidth',0.8,'Color',[1,0.25, 0.75])
% m_plot(360+lon_ax90,lat_ax90,'Linewidth',0.8,'Color',[1,0.25, 0.75])
% m_plot(360+lon_axcs,lat_axcs,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_axwbts,lat_axwbts,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_ax8,lat_ax8,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(lon_mx1,lat_mx1,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(lon_mx1a,lat_mx1a,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(360+lon_mx2,lat_mx2,'Linewidth',0.8,'Color',[1,0.25, 0.75])
m_plot(lon_mx4,lat_mx4,'Linewidth',0.8,'Color',[1,0.25, 0.75])

% BOM lines
T1 = importdata('reference_lines_PX02.csv');
m_plot(T1.data(:,2), T1.data(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('reference_lines_IX22-PX11.csv');
m_plot(T1.data(:,2), T1.data(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])
T1 = importdata('reference_lines_IX01.csv');
m_plot(T1.data(:,2), T1.data(:,1),'Linewidth',0.8,'Color',[1,0.25 0.75])



% transects described in paper
T1 = importdata('p40stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',1.5,'Color','r')
T1 = importdata('p31stnpos.fer');
m_plot(T1(:,2),T1(:,1),'Linewidth',1.5,'Color','r')

load reference_transect_AX97.mat
m_plot(360+lon_ax97,lat_ax97,'Linewidth',1.5,'Color','r')


load Reference_Section_PX36.mat
lon_px36 = lon;
lat_px36 = lat;
m_plot(lon_px36,lat_px36,'Linewidth',1.5,'Color','r')
m_plot(360+lon_ax32,lat_ax32,'Linewidth',1.5,'Color','r')


m_grid('box','fancy','tickdir','in');

set(gca,'plotboxaspectratio',[abs(lon2-lon1)/abs(lat2-lat1) 1 1]);

%print -djpeg Figure1_Map_LandOutlineOnly.jpg

print -dpng -r400 figure1_XBTMap.png

exportgraphics(gcf, 'figure1_XBTmap_new.png', 'Resolution', 400);

