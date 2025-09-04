%% quick_plot_6_graphs
% This is to be used if you already have the outmatfilename done from
% a previous run of the 'make_six_plots.m'. It skips all of the
% recalculation, and just does the plotting, useful for fine tuning the
% plot, not the values plotted.
%
% BY Jared Brzenski
% 2025-08-01
% 
% 
close all; clear all; clc;

% netCDF Filename is used for historical purpose, it does not actually
% download this data.

addpath([pwd '/helper_functions']);
addpath([pwd '/helper_functions/colorbar_yyw/']);
addpath([pwd '/data']);
addpath([pwd '/results']);
addpath([pwd '/m_map']);

netcdf_filename = 'data/p40tem1211_2312.nc';

% This si the filename of the processed data. see 'make_six_plots.m'
outmatfilename = 'data/output_p40.mat';

% The save name of the generated image.
save_filename = 'results/PX40_six-plots_400dpi.png';

% gridfile for the topography
gridfile = 'data/ETOPO1_Ice_g_gmt4.grd';

%% Plotting
%  set colormaps
disp('Plotting section')

clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
% clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
% clr_S = cmocean('haline');  clr_S(1,:) = [.7 .7 .7];
% clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
addpath([pwd '/helper_functions/colorbar_yyw/']);
clr_T = jet;
clr_S = colormap_yyw('MPL_rainbow');
clr_V = colormap_yyw('GMT_polar');
xyclr = [.4 .4 .4];
lon_lim = [139.5 165];                       % Define lon lim for map
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
T1 = importdata('p40stnpos.txt');
T1 = T1(1:104,:); % only to 165

% +++++++++++++++++++++++++++++++++++++++++++
npf = length(time);

label_lat = {'22\circN','','26\circN','','30\circN','','34\circN'};  %Define horizontal Label ticks
label_lon = {'140\circE', '145\circE', '150\circE', '155\circE', '160\circE','165\circE'};

%% Begin Plotting
fig1 = myFigSize(1,10,6); clf % papersize (10,6)

tiledlayout(2,3);
%%
%=========================================================================
nexttile(1);hold off
m_proj('miller', 'lon', lon_lim, 'lat', lat_lim);
[CS,CH]=m_contourf(lonb,latb,depth_btm,[-11000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
colormap(gca,clr_topo(50:end,:))
m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on
for i = 1:length(time)
    m_plot(lon_ref(1:2:end,i),lat_ref(1:2:end,i),'.r','markersize',4); hold on;
end
m_plot(lon_ref(:,1),mean(lat_ref,2),'k','linewidth',1)
m_grid('tickdir','out','color',xyclr)
% orig removed next line
set(gca,'Position',[.05 .56 .24 .36])
axis normal;
%orig: axis square;

ax=m_contfbar(0.2,[.05 .4],CS,CH,'axfrac',.025);
set(ax,'colormap',clr_topo(50:end,:))
set(ax,'tickdir','in','fontsize',8,'ycolor',xyclr,'xcolor',xyclr)
title(ax,{'m'},'fontweight','normal','fontsize',8);
m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,'(a) PX40')


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
tmp = temp_ref;
% % % Not NaNs, just really low < -1000!
% tmp(tmp>-5) = 26;
tmp(tmp<=-5) = NaN;
temp_ref = tmp;
nan_in_temp = isnan(tmp);
mtemp = nanmean(temp_ref,3);
tmp = mtemp;

lon_depth = repmat(squeeze(lon_ref(:,1)),1, length(depth_ref));
depth_plot = repmat(depth_ref, 1, length(lon))';

contourf(lon_depth,-depth_plot,tmp,0:32:26); hold on
colormap(gca,clr_T);
[CS,CH] = contourf(lon_depth,-depth_plot,mtemp,0:1:26,'LineStyle','none');
set(gca,'tickdir','out');
caxis([0 26])
set(gca,'xtick',[140:5:165]);
get(gca,'xtick')
set(gca,'xticklabel',label_lon);

plot(T1(:,2),0,'.k')
ax=m_contfbar(1.05,[.09 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[2 26],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_T)
title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.06 0.1 0.23 0.35])


ylabel('Depth (m)')
%set(gca,'xtick',lon_lim(1):10:lon_lim(2),'xticklabel',label_lon,'XTickLabelRotation',0)
tlab = text(double(lon_ref(1,1)),50,'(d) Mean temperature');

%=========================================================================
% mean salinity

saltref(nan_in_temp) = NaN;
msalt = nanmean(salt_ref,3);

% % msalt(msalt<32.5) = 32.5;
tmp = msalt;
% % tmp(~isnan(tmp)) = 37;
% % tmp(isnan(tmp)) = 27;

nexttile(5); hold off
contourf(lon_depth,-depth_plot,tmp,33.8:15:35.5); hold on
colormap(gca,clr_S);
[CS,CH] = contourf(lon_depth,-depth_plot,msalt,33.8:0.1:35.5,'LineStyle','none');
set(gca,'tickdir','out');
caxis([33.8 35.5])
% get(gca,'xtick')
% set(gca,'xtick',[140:5:165]);
% get(gca,'xtick')
set(gca,'xticklabel',label_lon);


ax=m_contfbar(1.05,[.04 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[33.8 35.5],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_S)
title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.38 0.1 0.23 0.35])
%set(gca,'xtick',lon_lim(1):10:lon_lim(2),'xticklabel',label_lon,'XTickLabelRotation',0)
text(double(lon_ref(1,1)),50,'(e) Mean salinity')


%=========================================================================
% Calculate mean velocity
mvel = mean(vel_ref,3,'omitnan');

tmp = mvel;
tmp(~isnan(tmp)) = 0.65;
tmp(isnan(tmp)) = -0.65;

nexttile(6); hold off
depth_plot_v = depth_plot(2:end,:);
lon_depth_v = lon_depth(2:end,:);

contourf(lon_depth_v,-depth_plot_v,tmp,-0.7:1.4:0.7); hold on
colormap(gca,clr_V2);
[CS,CH] = contourf(lon_depth_v,-depth_plot_v,mvel,-0.6:0.05:0.6,'LineStyle','none');
contour(lon_depth_v,-depth_plot_v,mvel,[0 0],'color',[.7 .7 .7]);
set(gca,'tickdir','out');
caxis([-0.6 0.6])
get(gca,'xtick')
set(gca,'xtick',[140:5:165]);
get(gca,'xtick')
set(gca,'xticklabel',label_lon);


ax=m_contfbar(1.05,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[-0.6 0.6],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',-0.6:0.2:0.6,'yticklabel',-0.6:0.2:0.6)
set(ax,'colormap',clr_V2)
title(ax,'  m/s','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.7 0.1 0.23 0.35])
%set(gca,'xtick',lon_lim(1):10:lon_lim(2),'xticklabel',label_lon,'XTickLabelRotation',0)
text(double(lon_ref(1,1)),50,'(f) Mean velocity')

%% Pie Chart
% ax_pie = nexttile(3);
% ax_pie.Color = 'none';
%n_vals = sum(~isnan(OG_lat), 1);
drops = [ 190 187 179 173 171 196 179 194 187 185 ...
          190 207 192 193 192 190 188 187 189 192 ...
          190 164 187 188 193 191 184 185 188 184 ...
          190 189 174 171 182 189 188 198 196 193 ...
          198 ]';
      
      
      
      % All below is new from Marlos
timeYMD2 = datestr((datenum(timeYMD))','yyyymmdd');
data1 = str2num(timeYMD2(:,5:6))'; %Month numbers
data2 = str2num(timeYMD2(:,1:4))'; %year numbers

%GET DATA RANGE
axis1 = [min(data1):max(data1)]';
axis2 = [min(data2):max(data2)+1]';

%create month axis if wanted
leg1 = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};

%Add N number of spaces in the middle
N = 3;
C = drops; %number of profiles/transect
%axis2(end+1) = axis2(1)-1;
%axis2 = [axis2(end); axis2(1:end-1)];
leg2 = cell(1,length(axis2));
leg2{end} = '';
for ii=1:length(axis2)-1
    aux = num2str(axis2(ii));
    leg2{ii} = aux(end-1:end);
end
%MAKE PIE PLOTS
%1: month, year
%clf
[g,col] = make_pie_subplot2(axis2,axis1,leg2,leg1,data2,data1,C,N,[.62 .51 .345 .415]);
title(col,'# prof')

      
      
% % % from Jared original 
%%%%%plotPieChart( ax_pie, timeYMD, drops);
% % % ax_pie.Color='none';
% % % box(ax_pie, 'off');
% % % axis(ax_pie, 'off');

%% Save File
%print( gcf, save_filename, "-dpng", "-r400" );

print -dpng -r400 figure1_px40.png

