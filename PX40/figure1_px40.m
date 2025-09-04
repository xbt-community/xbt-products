%%  set colormaps

clear

disp('Plotting section')
addpath /phodnet/share/mgoes/matlab/m_map/ 
addpath ../Colormap/colorbar_yyw/
addpath ../Colormap/
clr_T = jet;
clr_S = colormap_yyw('MPL_rainbow');
clr_V = colormap_yyw('GMT_polar');

clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
%clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
%clr_S = cmocean('haline'); clr_S(1,:) = [.7 .7 .7];
%clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
xyclr = [.4 .4 .4];
lon_lim = [140 165];                                   %Define lon lim for map
lat_lim = [15 40];                                     %Define lat lim for map

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
load output_p40.mat
npf = [ 190 187 179 173 171 196 179 194 187 185 ...
          190 207 192 193 192 190 188 187 189 192 ...
          190 164 187 188 193 191 184 185 188 184 ...
          190 189 174 171 182 189 188 198 196 193 ...
          198 ]';
      

tick_lat = 145:10:165;
label_lat = strcat(num2str(tick_lat'),'\circE');

%label_lat = {'34\circN','','36\circN','','38\circN','','40\circN'};  %Define horizontal Label ticks

fig1 = myFigSize(1,10,6); clf % papersize (10,6)
subplot('position',[-0.05  0.56  0.52  0.36])
subplot(2,3,1);hold off
m_proj('miller','lon',lon_lim,'lat',lat_lim);
[CS,CH]=m_contourf(lonb,latb,depth_btm,[-7000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
colormap(gca,clr_topo(50:end,:))
m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on
for i = 1:length(npf)
%    m_plot(lon_xbt{i},lat_xbt{i},'.r','markersize',4)
    m_plot(lon_ref(1:2:end,i),lat_ref(1:2:end,i),'.r','markersize',4); hold on;
end
m_plot(lon_ref(:,1),mean(lat_ref,2),'k','linewidth',1)


%m_plot(mean(lon_ref,2),lat_ref,'k','linewidth',1)
m_grid('tickdir','in','color',xyclr,'xtick',tick_lat)
set(gca,'Position',[.01 .56 .34 .36])
ax=m_contfbar(0.32,[.12 .49],CS,CH,'axfrac',.025);
set(ax,'colormap',clr_topo(50:end,:))
set(ax,'tickdir','in','fontsize',8,'ycolor',xyclr,'xcolor',xyclr)
title(ax,{'m'},'fontweight','normal','fontsize',8);
m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,'(a) PX40')
%axis normal

tick_lat = 145:5:165;
label_lat = strcat(num2str(tick_lat'),'\circE');

% histogram # of cruises (month)
subplot(2,3,2); hold off
[~,m] = datevec(timeYMD);%time_avg);
h = histogram(m,0.5:1:12.5); hold on
set(h,'EdgeColor',[0.1294 0.1294 0.1294])
set(gca,'XLim',[0.5 12.5],'YLim',[0 12])
set(gca,'Position',[0.38 0.56 0.23 0.36])
%text(0.5,12.8,'(b) Monthly distribution of cruises')
y_lim = ylim;
text(0.5,y_lim(2)+diff(y_lim)/15,'(b) Monthly distribution of cruises')
ylabel('number of cruises')

% mean temperature
subplot(2,3,4); hold off
mtemp = squeeze(mean(temp_ref,3,'omitmissing'))';
tmp = mtemp;
tmp(~isnan(tmp)) = 26;
tmp(isnan(tmp)) = -5;
[CS,CH] = contourf(lon_ref(:,1),-depth_ref,mtemp,2:1:26,'LineStyle','none');
clim([0 26])
colormap(gca,clr_T);
set(gca,'color',[.7 .7 .7])%,'xdir','reverse')
hold on
plot(lon(:,10),0,'.k')
ax=m_contfbar(0.82,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[2 26],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_T)
title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.06 0.1 0.23 0.35])
ylabel('Depth (m)')
set(gca,'xtick',tick_lat,'xticklabel',label_lat,'XTickLabelRotation',0)
%set(gca,'xtick',34:1:40,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lon_ref(1),50,'(d) Mean temperature')

% mean salinity
msalt = mean(salt_ref,3,'omitmissing')';
msalt(msalt<32.5) = 32.5;
tmp = msalt;
tmp(~isnan(tmp)) = 37;
tmp(isnan(tmp)) = 27;
subplot(2,3,5); hold off
%contourf(lat_ref,-depth_ref,tmp,22:15:37); hold on
%colormap(gca,clr_S);
[CS,CH] = contourf(lon_ref(:,1),-depth_ref,msalt,34:0.1:35,'LineStyle','none');
clim([34 35])
set(gca,'color',[.7 .7 .7])%,'xdir','reverse')
colormap(gca,clr_S);
ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[34 35],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_S)
title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.38 0.1 0.23 0.35])
set(gca,'xtick',tick_lat,'xticklabel',label_lat,'XTickLabelRotation',0)
%set(gca,'xtick',34:1:40,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lon_ref(1),50,'(e) Mean salinity')

% mean velocity
mvel = mean(vel_ref,3,'omitmissing')';
tmp = mvel;
tmp(~isnan(tmp)) = 0.65;
tmp(isnan(tmp)) = -0.65;
subplot(2,3,6); hold off
%contourf(lat_vel,-depth_ref,tmp,-0.7:1.4:0.7); hold on
%colormap(gca,clr_V2);
[CS,CH] = contourf(lon_vel(:,10),-depth_ref,mvel,-0.6:0.05:0.6,'LineStyle','none');
hold on
contour(lon_vel(:,10),-depth_ref,mvel,[0 0],'color',[.7 .7 .7]);
clim([-0.65 0.65])
set(gca,'color',[.7 .7 .7])%,'xdir','reverse')
colormap(gca,clr_V2);
ax=m_contfbar(1.2,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[-0.6 0.6],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',-0.6:0.2:0.6,'yticklabel',-0.6:0.2:0.6)
set(ax,'colormap',clr_V2)
title(ax,'m/s','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.7 0.1 0.23 0.35])
%            set(gca,'xdir','reverse')
set(gca,'xtick',tick_lat,'xticklabel',label_lat,'XTickLabelRotation',0)
%set(gca,'xtick',34:1:40,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lon_ref(1),50,'(f) Mean velocity')

%PIE
%CREATE X,Y AXES
timeYMD1 = datestr(timeYMD','yyyymmdd');
data1 = str2num(timeYMD1(:,5:6))'; %Month numbers
data2 = str2num(timeYMD1(:,1:4))'; %year numbers

%GET DATA RANGE
axis1 = [min(data1):max(data1)]';
axis2 = [min(data2):max(data2)+1]';

%create month axis if wanted
leg1 = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};

%Add N number of spaces in the middle
N = 3;
C = npf; %number of profiles/transect
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
%[g,col] = make_pie_subplot2(axis2,axis1,leg2,leg1,data2,data1,C,N,[.62 .55 .345 .415]);
title(col,'# prof')
%text(1,13,'(c) Transects')
text(-3.5,-.5,'(c)')

%print -dpng -r400 figure1_ax32.png
exportgraphics(gcf, 'figure1_ax40.png', 'Resolution', 400);