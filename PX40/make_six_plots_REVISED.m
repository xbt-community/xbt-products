%% New PLot for netCDF files
%
%
% BY Jared Brzenski
%
% Requires the GSW Toolbox, available:
% https://www.teos-10.org/software.htm
%
% ALSO requires the M_Map toolbox, available:
% https://www-old.eoas.ubc.ca/~rich/map.html
%
% Add the folder HERE, in this working path, and it will be added to the
% path by this script.
%
% XBT Data 2025 Scripps Institute of Oceanography
%
% Data should be in the format ( per GSW )
% DEPTH x LAT?LON x Time
%
% Should only have to change 'filename' below.
close all; clear all; clc;

netcdf_filename = 'data/p40tem1211_2312.nc';
outmatfilename = 'data/output_p40.mat';
save_filename = 'results/PX40_six-plots_400dpi.png';
gridfile = 'data/ETOPO1_Ice_g_gmt4.grd';

BLim = [19.8 34.9];
%BLim = [36.5 38.5];

%% Check if GSW Tools is available
exists = isfolder('GSW_Tools');
if exists == 0
    disp('Need GSW tools to continue.');
    disp('Please download and place folder here as GSW_Tools');
    return 
end

% Add folder and all subfolders
addpath(genpath('GSW_Tools'));

%% Check if m_map is available
exists = isfolder('m_map');
if exists == 0
    disp('Need m_map tools to plot.');
    disp('Will run, but will not print plots :(');
end

if exists == 1
    addpath(genpath('m_map'));
end

%% Check if Helper_functions is available
exists = isfolder('helper_functions');
if exists == 0
    disp('Need helper functions for plot coloring');
    disp('Will not display properly :(');
    return
end

if exists == 1
    addpath(genpath('helper_functions'));
end
%% Check if Salinity Calculation function ( Thacker and Goes ) is available
exists = isfolder('Calc_sal');
if exists == 0
    disp('Need Calc_sal_Thacker_Goes_EmDr_Stom function for salinity.');
    disp('Will not calculate salinity :(');
    return
end

if exists == 1
    addpath(genpath('Calc_sal'));
end

%% Load netCDF Data
loadNetCDFData( netcdf_filename, outmatfilename);

%% =======================Marlos Salinity===========
disp('Calc Salinity (linear)')
calculateSalinity( outmatfilename );

%% Compute Velocity with GSW Toolbox
% GSW Wants (DEPTH x LAT/LON x Time).
% Adjust as needed

disp('Calc Velocity profiles')
calculateVelocityProfiles( outmatfilename );

%% Calculate transport - 
disp('Calc transport')
calculateTransport( outmatfilename );

%% Calculate Seasonal Cycle - removed for now as not needed
% % disp('Cal seasonal cycle')
% % calculateSeason( outmatfilename );

%% Plotting
%%  set colormaps
disp('Plotting section')

clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
% clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
% clr_S = cmocean('haline'); clr_S(1,:) = [.7 .7 .7];
% clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
addpath([pwd '/helper_functions/colorbar_yyw/']);
clr_T = jet;
clr_S = colormap_yyw('MPL_rainbow');
clr_V = colormap_yyw('GMT_polar');

xyclr = [.4 .4 .4];
lon_lim = [138 165];                       % Define lon lim for map
lat_lim = [16 40];                                     % Define lat lim for map

%%  ======= Figures =======
% load topography
x=ncread(gridfile,'x');
y=ncread(gridfile,'y');
z=ncread(gridfile,'z');

x(x<0) = x(x<0) + 360;
[x,idx] = sort(x);
z = z(idx,:);

kx=find(x>=lon_lim(1)-0.25 & x<=lon_lim(2)+0.25);
ky=find(y>=lat_lim(1)-0.25 & y<=lat_lim(2)+0.25);
lonb=x(kx); latb=y(ky); depth_btm=z(kx,ky)';

depth_btm(depth_btm>0)=NaN;
clear kx ky x y z

% +++++++++++++++++++++++++++++++++++++++++++
% Load data
load(outmatfilename);
% +++++++++++++++++++++++++++++++++++++++++++
npf = length(time);

label_lat = {'22\circN','','26\circN','','30\circN','','34\circN'};  %Define horizontal Label ticks
label_lon = {'140\circE', '', '150\circE', '', '160\circE', '', '170\circE'};



%% Begin Plotting
%

fig1 = myFigSize(1,10,6); clf % papersize (10,6)

tiledlayout(2,3);
%=========================================================================
nexttile(1);hold off
m_proj('miller', 'lon', lon_lim, 'lat', lat_lim);
[CS,CH]=m_contourf(lonb,latb,depth_btm,[-11000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
colormap(gca,clr_topo(50:end,:))
m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on
for i = 1:length(time)
    m_plot(lon_ref(:,i),lat_ref(:,i),'.r','markersize',4); hold on;
end
m_plot(lon_ref(:,1),mean(lat_ref,2),'k','linewidth',1)
m_grid('tickdir','in','color',xyclr)
%set(gca,'Position',[.06 .56 .28 .46])
axis square;

ax=m_contfbar(0.2,[.05 .4],CS,CH,'axfrac',.025);
set(ax,'colormap',clr_topo(50:end,:))
set(ax,'tickdir','in','fontsize',8,'ycolor',xyclr,'xcolor',xyclr)
title(ax,{'m'},'fontweight','normal','fontsize',8);
m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,'(a) AX32')


%=========================================================================
% histogram # of cruises (month)
nexttile(2); hold off
[~,m] = datevec(timeYMD);
h = histogram(m,0.5:1:12.5); hold on
set(h,'EdgeColor',[0.1294 0.1294 0.1294])
set(gca,'XLim',[0.5 12.5],'YLim',[0 8])
set(gca,'Position',[0.38 0.56 0.23 0.36])
text(0.5,8.8,'(b) Monthly distribution of cruises')
ylabel('number of cruises')

%=========================================================================
% mean temperature
nexttile(4); hold off

% some bad values at beginning
mtemp = mean(temp_ref,3);
tmp = mtemp;
% Not NaNs, just really low < -1000!
% tmp(tmp>-5) = 26;
% tmp(tmp<=-5) = 0;
%tmp(~isnan(tmp)) = 26;
%tmp(isnan(tmp)) = -5;

lon_depth = repmat(squeeze(lon_ref(:,1)),1, length(depth_ref));
depth_plot = repmat(depth_ref, 1, length(lon))';

contourf(lon_depth,-depth_plot,tmp,-6:32:26); hold on
%contourf(lat_ref,-depth_ref,tmp,-6:32:26); hold on
colormap(gca,clr_T);
%[CS,CH] = contourf(lat_ref(:,1),-depth_ref,mtemp,2:1:26,'LineStyle','none');
[CS,CH] = contourf(lon_depth,-depth_plot,mtemp,0:1:26,'LineStyle','none');

caxis([0 26])
plot(lon(:,1),0,'.k')

ax=m_contfbar(1.05,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[2 26],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_T)
title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.06 0.1 0.23 0.35])
ylabel('Depth (m)')
set(gca,'xtick',lon_lim(1):5:lon_lim(2),'xticklabel',label_lon,'XTickLabelRotation',0)
tlab = text(double(lon_ref(1,1)),50,'(d) Mean temperature');


%=========================================================================
% mean salinity
msalt = mean(salt_ref,3,'omitnan');
msalt(msalt<32.5) = 32.5;
tmp = msalt;
tmp(~isnan(tmp)) = 37;
tmp(isnan(tmp)) = 27;
nexttile(5); hold off
contourf(lon_depth,-depth_plot,tmp,22:15:36); hold on
colormap(gca,clr_S);
[CS,CH] = contourf(lon_depth,-depth_plot,msalt,32:0.2:36,'LineStyle','none');
caxis([32 36])
ax=m_contfbar(1.05,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[32 36],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_S)
title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.38 0.1 0.23 0.35])
set(gca,'xtick',lon_lim(1):10:lon_lim(2),'xticklabel',label_lon,'XTickLabelRotation',0)
text(double(lon_ref(1,1)),50,'(e) Mean salinity')


%=========================================================================
% mean velocity
mvel = mean(vel_ref,3,'omitnan');
tmp = mvel;
tmp(~isnan(tmp)) = 0.65;
tmp(isnan(tmp)) = -0.65;
nexttile(6); hold off
% lat_depth_v = lat_depth(2:end,:);
depth_plot_v = depth_plot(2:end,:);
lon_depth_v = lon_depth(2:end,:);

contourf(lon_depth_v,-depth_plot_v,tmp,-0.7:1.4:0.7); hold on
colormap(gca,clr_V2);
[CS,CH] = contourf(lon_depth_v,-depth_plot_v,mvel,-0.6:0.05:0.6,'LineStyle','none');
contour(lon_depth_v,-depth_plot_v,mvel,[0 0],'color',[.7 .7 .7]);
caxis([-0.65 0.65])
ax=m_contfbar(1.05,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[-0.6 0.6],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',-0.6:0.2:0.6,'yticklabel',-0.6:0.2:0.6)
set(ax,'colormap',clr_V2)
title(ax,'  m/s','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.7 0.1 0.23 0.35])
set(gca,'xtick',lon_lim(1):10:lon_lim(2),'xticklabel',label_lon,'XTickLabelRotation',0)
text(double(lon_ref(1,1)),50,'(f) Mean velocity')

%% Pie Chart
 ax_pie = nexttile(3);
% %n_vals = sum(~isnan(OG_lat), 1);
drops = [ 190 187 179 173 171 196 179 194 187 185 ...
          190 207 192 193 192 190 188 187 189 192 ...
          190 164 187 188 193 191 184 185 188 184 ...
          190 189 174 171 182 189 188 198 196 193 ...
          198 ];
plotPieChart( ax_pie, timeYMD, drops);







% 
% OR, If already have a nice plot!
% img = imread('pie_chart.png');
% imshow(img, 'Parent', ax_pie, 'InitialMagnification',150);
% axis(ax_pie, 'image');
