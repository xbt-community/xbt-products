%%
close all
clear

%Define parameters for each section
direc = {'AX32','AX97','PX30','PX36','PX40'};
letters={'a) ','b) ','c) ','d) ','e) '};
vlim = [.8 .2 .7 .2 .8];         %velocity maximum for colorbar
labE = {'^oN','^oW','^oE','^oS','^oE'};
%stp = [0.05 0.01 0.05 0.05 .05];
lstp = [5,4,10,10,10];           %xtick step for map
lalim = [ 31.5000   41.5000;...
    -24.4400 -17.8500;...
    -33.9030 -14.9530;...
    -73.8300 -44.0900;...
     21.3750  37.1980];         %latitude limits for T/S
lolim = [-74.9470  -62.9720;...
   -41.9000  -28.7000;...
    152.1300  179.3300;...
    165.0270  185.9820;...
    139.00    165.00];           %longitude limits for T/S
blim = [32.5500  40.4500;...
    -40.8167  -38;...
     153.1800  160.00;...
    -71.8300  -46.0900;...
     139.00    150.00];          %main axis limits for velocity


%%  set colormaps
disp('Plotting section')
addpath ../Colormap/
addpath /phodnet/share/mgoes/matlab/m_map
addpath ../Colormap/colorbar_yyw/
clr_T = jet;
clr_S = colormap_yyw('MPL_rainbow');
clr_V = colormap_yyw('GMT_polar');

clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256));
clr_V2 = interp1(1:length(clr_V2),clr_V2,linspace(1,length(clr_V2),14));
xyclr = [.4 .4 .4];

% load topography
x=ncread('ETOPO1_Ice_g_gmt4.grd','x');
y=ncread('ETOPO1_Ice_g_gmt4.grd','y');
z=ncread('ETOPO1_Ice_g_gmt4.grd','z');


for ii = 1:length(direc)
    file = [direc{ii},'/',direc{ii},'_gridded_new.nc'];
    if ~exist(file,'file')
        file = [direc{ii},'/',direc{ii},'_gridded.nc'];
    end

 %LOAD DATA
    vel_ref = ncread(file,'velocity');
    temp_ref = ncread(file,'temperature');
    salt_ref = ncread(file,'reference_salinity');
    time_avg = ncread(file,'time') + datenum(1950,1,1);
    lon_xbt = ncread(file,'longitude');
    lat_xbt = ncread(file,'latitude');
    depth_ref = ncread(file,'depth');

  %Load number of profiles
  load npf_goxbt.mat
  eval(['npf=npf_',direc{ii},';'])

  %load reference sample of transects
  load latlon_goxbt.mat
  eval(['lon_ref = lon_',direc{ii},';'])
  eval(['lat_ref = lat_',direc{ii},';'])

    [nz_ref,nx_ref,nt] = size(temp_ref);

%    mvel = median(vel_ref,3,'omitmissing');
%    mask = double(~isnan(mvel)); mask(mask==0)=nan;

    %Define map boundaries
    Blim = blim(ii,:);
    lon_lim = lolim(ii,:);   %Define lon lim for map
    lat_lim = lalim(ii,:);   %Define lat lim for map
   
    %  ======= Define map bathymetry =======
    x2 = x;
    z2 = z;
    if ii==4
        x2(x2<0)=x2(x2<0)+360;
        [x2,ilon] = sort(x2);
        z2 = z2(ilon,:);
    end
    kx = find(x2>=lon_lim(1)-0.25 & x2<=lon_lim(2)+0.25);
    ky = find(y>=lat_lim(1)-0.25 & y<=lat_lim(2)+0.25);
    lonb=x2(kx); latb=y(ky); depth_btm=z2(kx,ky)';
    depth_btm(depth_btm>0)=NaN;
    clear kx ky 

    %  ========== Define main axis==========
    if any(ii==[2 3 5])
        ax_ref = lon_xbt;
        lon_xbt = lon_xbt(:,ones(1,nt));
     %   lat_xbt(isnan(squeeze(temp_ref(1,:,:))))=nan;
        ax_sample = lon_ref;
   else
        ax_ref = lat_xbt;
        lat_xbt = lat_xbt(:,ones(1,nt));
     %   lon_xbt(isnan(squeeze(temp_ref(1,:,:))))=nan;
        ax_sample = lat_ref;
    end

    %%   ============ Make Figure ====================
    fig1 = myFigSize(ii,10,6); clf % papersize (10,6)

    %map of profiles
    subplot(2,3,1);hold off
    m_proj('miller','lon',lon_lim,'lat',lat_lim);
    [CS,CH]=m_contourf(lonb,latb,depth_btm,[-7000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
    colormap(gca,clr_topo(50:end,:))
    m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on

    for i = 1:length(npf)
        m_plot(lon_xbt(:,i),lat_xbt(:,i),'.r','markersize',4)
    end

    m_plot(mean(lon_xbt,2,'omitmissing'),mean(lat_xbt,2,'omitmissing'),'k','linewidth',1)
    m_grid('tickdir','in','color',xyclr,'xtick',round(lon_lim(1)):lstp(ii):round(lon_lim(2)))
    if ii==4
        set(gca,'Position',[.05 .56 .24 .36])
        ax=m_contfbar(0.05,[.57 .94],CS,CH,'axfrac',.025);
    elseif ii==1 || ii==5
        set(gca,'Position',[.05 .56 .24 .36])
        ax=m_contfbar(0.17,[.10 .47],CS,CH,'axfrac',.025);
        axis normal
    else 
        set(gca,'Position',[.05 .56 .24 .36])
        ax=m_contfbar(0.17,[.57 .94],CS,CH,'axfrac',.025);
        axis normal
    end
    set(ax,'colormap',clr_topo(50:end,:))
    set(ax,'tickdir','in','fontsize',8)%,'ycolor',xyclr,'xcolor',xyclr)
    title(ax,{'m'},'fontweight','normal','fontsize',8);
    m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,[letters{ii},direc{ii}])

    % histogram # of cruises (month)
    subplot(2,3,2); hold off
    [~,m] = datevec(time_avg);
    h = histogram(m,0.5:1:12.5); hold on
    set(h,'EdgeColor',[0.1294 0.1294 0.1294])
    axis tight
    set(gca,'XLim',[0.5 12.5])
    set(gca,'Position',[0.38 0.56 0.23 0.36])
    y_lim = ylim;
    text(0.5,y_lim(2)+diff(y_lim)/15,'(b) Monthly distribution of cruises')
    ylabel('number of cruises')

    % mean temperature
    subplot(2,3,4); hold off
    mtemp = mean(temp_ref,3,'omitmissing');
    tlim = linspace(min(mtemp(:)),max(mtemp(:)),17);
    [CS,CH] = contourf(ax_ref,-depth_ref,mtemp,tlim,'LineStyle','none');
    clim([ceil(tlim(1)) floor(tlim(end))])
    colormap(gca,clr_T);
    set(gca,'color',[.7 .7 .7])
    hold on
    plot(ax_sample,0,'.k')
    ax=m_contfbar(0.82,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
    set(ax,'tickdir','in','fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
    set(ax,'ylim',[(tlim(1)) (tlim(end))])
    set(ax,'colormap',clr_T)
    title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')
    set(gca,'Position',[0.06 0.1 0.23 0.35])
    ylabel('Depth (m)')
    axt = ax_ref(1);
    if ii==1
        set(gca,'xdir','reverse')
    end
    if ii==1 || ii==4
       axt = ax_ref(end);
    end
       text(axt,50,'(d) Mean temperature')
  
    tick_lon = get(gca,'xtick');
    label_lon = strcat(num2str(abs(tick_lon)'),labE{ii});

    set(gca,'xtick',tick_lon,'xticklabel',label_lon,'XTickLabelRotation',0)


    % mean salinity
    msalt = mean(salt_ref,3,'omitmissing');
    slim = linspace(min(msalt(:)),max(msalt(:)),17);
    subplot(2,3,5); hold off
    [CS,CH] = contourf(ax_ref,-depth_ref,msalt,slim,'LineStyle','none');
    clim([slim(1) slim(end)])
    colormap(gca,clr_S);
    set(gca,'color',[.7 .7 .7])
    ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
    set(ax,'tickdir','in','ylim',[slim(1) slim(end)],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
    set(ax,'colormap',clr_S)
    title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')
    set(gca,'Position',[0.38 0.1 0.23 0.35])
    axt = ax_ref(1);
    if ii==1
        set(gca,'xdir','reverse')
    end
    if ii==1 || ii==4
       axt = ax_ref(end);
    end
    text(axt,50,'(e) Mean salinity')
    set(gca,'xtick',tick_lon,'xticklabel',label_lon,'XTickLabelRotation',0)

    % mean velocity
    mvel = mean(vel_ref,3,'omitmissing'); 
    subplot(2,3,6); hold off
    [CS,CH] = contourf(ax_ref,-depth_ref,mvel,linspace(-vlim(ii),vlim(ii),14),'LineStyle','none');
    clim([-vlim(ii) vlim(ii)]) 
    pos = get(gca,'Position');
    colormap(gca,clr_V2);
    set(gca,'color',[.7 .7 .7])
    xlim(Blim)

     ax=m_contfbar(1.2,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
     set(ax,'tickdir','in','ylim',[-vlim(ii) vlim(ii)],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
   %  set(ax,'ytick',[-vlim(ii):.2: vlim(ii)],'yticklabel',[-vlim(ii):.2: vlim(ii)])
    set(ax,'colormap',clr_V2)
    title(ax,'m/s','FontSize',10,'color',xyclr,'FontWeight','normal')
    set(gca,'Position',[0.7 0.1 0.23 0.35])
    axt = ax_ref(1);
    if ii==1
        set(gca,'xdir','reverse')
    end
    if ii==1 || ii==4
       axt = ax_ref(end);
    end

    text(axt,50,'(f) Mean velocity')

    xlim(blim(ii,:))
    tick_lon = get(gca,'xtick');
    tick_lon = tick_lon(tick_lon-round(tick_lon)==0); %just integers
    label_lon = strcat(num2str(abs(tick_lon)'),labE{ii});

    set(gca,'xtick',tick_lon,'xticklabel',label_lon,'XTickLabelRotation',0)

%end
%PIE

%CREATE X,Y AXES
timeYMD = datestr(floor(time_avg),'yyyymmdd');
data1 = str2num(timeYMD(:,5:6))'; %Month numbers
data2 = str2num(timeYMD(:,1:4))'; %year numbers

%GET DATA RANGE
axis1 = [min(data1):max(data1)]';
axis2 = [min(data2):max(data2)+1]';

%create month axis if wanted
leg1 = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};

%Add N number of spaces in the middle
N = 3;
C = double(npf'); %number of profiles/transect
leg2 = cell(1,length(axis2));
leg2{end} = '';
for i=1:length(axis2)-1
    aux = num2str(axis2(i));
    leg2{i} = aux(end-1:end);
end

%MAKE PIE PLOTS
%1: month, year
%clf
[g,col] = make_pie_subplot2(axis2,axis1,leg2,leg1,data2,data1,C,N,[.62 .51 .345 .415]);
title(col,'# prof')
text(-3.5,-.5,'(c)')

exportgraphics(gcf, ['figure',num2str(ii),direc{ii},'.png'], 'Resolution', 400);
end
