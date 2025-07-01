%%
%Packages needed: GSW, M_MAP, cmocean, TS_PACKAGE
addpath GSW_2015/
addpath m_map/
addpath TS_PACKAGE_THACKER_noyear_globe/

load temp_ax32_paper.mat
%Input: depth_xbt,temp_xbt,lat_xbt,lon_xbt,
%Output: time_xbt temp_ref lon_ref lat_ref depth_ref npf lon_xbt lat_xbt
disp('calc reference section')

% grid to reference section [lat_ref & lon_ref] and every 10 dbar [pres_ref]
kz = find(depth_xbt<=760);
depth_ref = depth_xbt(kz);
lat_ref = (32.5:0.1:40.5)'; % laittude grid for the reference section
nx_ref = length(lat_ref);
nz_ref = length(depth_ref);
nt = length(npf);
nz = length(depth_xbt);

% find the reference longitude based on averaged longitude from all transect after removing outliers
lon_ref = nan(nx_ref,nt);
for i=1:nt
    lon_ref(:,i) = interp1(lat_xbt{i},lon_xbt{i},lat_ref,'linear','extrap');
end
dlon=(lon_ref-mean(lon_ref,2,'omitmissing'));
kgood=find(mean(abs(dlon),1,'omitmissing')<=0.25);
mean_lon_ref=mean(lon_ref(:,kgood),2,'omitmissing');

% fill in gaps, then interpolate to reference latitude (lat_ref)
temp_ref = nan(length(depth_ref),nx_ref,nt);
for i = 1:nt
    tempi = temp_xbt{i};
    lati = lat_xbt{i};
    loni = lon_xbt{i};
    % if two profile are close together, which may induce anomous large velocity
    %   average to mid-latitude
    dist_y = deg2km(distance(lati(1:end-1),loni(1:end-1),lati(2:end),loni(2:end),'degrees'));
    ik = find(dist_y<10);
    if ~isempty(ik)
        for j = 1:length(ik)
            ik0 = ik(j);
            tempi(:,ik0) = mean(tempi(:,ik0:ik0+1),2,'omitmissing');
            lati(ik0) = mean(lati(ik0:ik0+1));
            loni(ik0) = mean(loni(ik0:ik0+1));
        end
        tempi(:,ik+1) = [];
        lati (ik+1) = [];
        loni (ik+1) = [];
    end
    for iz = 1:nz_ref
        % fill in gap
        ki = find(~isnan(tempi(iz,:)));
        if length(ki) >= length(lati)*0.8 & length(ki)<length(lati)
            tempi(iz,:) = interp1(lati(ki),tempi(iz,ki)',lati)';
        end
        % interpolate to reference grid
        ki = find(~isnan(tempi(iz,:)));
        temp_ref(iz,:,i) = interp1(lati(ki),tempi(iz,ki)',lat_ref)';
    end
end
%save('Output.mat','time_xbt','temp_ref','lon_ref','lat_ref','depth_ref')
save Output.mat time_xbt temp_ref lon_ref lat_ref depth_ref npf lon_xbt lat_xbt
disp('Done!')
clear 
%% =======================Marlos Salinity===========
%Input:time_xbt,temp_ref,lon_ref,lat_ref,depth_ref
%Output:salt_ref (Goes method),time_avg 
disp('calc sal')

load Output.mat
nt = length(npf);
nx_ref = length(lat_ref);
nz_ref = length(depth_ref);
time_avg = nan(nt,1);  %averaged observation time for the trensect
for it = 1:nt
     time_avg(it) = mean(time_xbt{it},'omitmissing');
end
timeYMD = datestr(time_avg','yyyymmdd');
time2 = str2num(timeYMD);
time2 = time2(:,ones(1,nx_ref))';

%Calc Salinity - All at once
auxT = reshape(temp_ref,[nz_ref nx_ref*nt]);
auxLon = reshape(lon_ref,[1 nx_ref*nt]);
auxLat = reshape(lat_ref(:,ones(1,nt)),[1 nx_ref*nt]);
auxTime = reshape(time2,[1 nx_ref*nt]);
salt_ref = Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(auxT,depth_ref,auxLat,auxLon,auxTime,0,depth_ref,'goes');

%Reshape back to initial size
salt_ref = reshape(salt_ref,[nz_ref nx_ref nt]);

save -append Output.mat salt_ref time_avg
disp('Done!')
clear

%% compute velocity (GSW toolbox)
% first compute pressure from depth at each grid
%load /phodnet/share/dong/gulfstream/oleander/temp_ax32_gridded.mat
%Input: temp_ref,depth_ref,salt_ref,lat_ref,lon_ref
%Output:vel_ref lat_vel lon_vel depth_ref dh0_ref
disp('compute velocity. Takes some time to finish')

load('Output.mat')

[nz_ref,nx_ref,nt] = size(temp_ref);
[~,ga_refInd] = min(abs(depth_ref-760));  %Reference for dynamic height output (760 m)
pres_ref = nan(nz_ref,nx_ref);
for ix = 1:nx_ref
    pres_ref(:,ix) = gsw_p_from_z(-depth_ref,lat_ref(ix));
end
vel_ref = nan(nz_ref,nx_ref-1,nt);
lon_vel = nan(nx_ref-1,nt);
dh0_ref = nan(nx_ref,nt);
for it = 1:nt
    [SA, in_ocean] = gsw_SA_from_SP(salt_ref(:,:,it),pres_ref,lon_ref(:,it),lat_ref);
    CT = gsw_CT_from_t(SA,temp_ref(:,:,it),pres_ref);
    ga = gsw_geo_strf_dyn_height(SA,CT,pres_ref,0);
    [vel0,mid_lat,mid_lon] = gsw_geostrophic_velocity(ga,lon_ref(:,it),lat_ref,pres_ref);
    lat_vel = mid_lat(1,:)';
    lon_vel(:,it) =  mid_lon(1,:)';
    vel_ref(:,:,it) = -vel0;
    dh0_ref(:,it) = -ga(ga_refInd,:)/9.81; %Convert from m2/s2 to meter at surface

end
    % vel_ref0 is velocity referenced to surface
    % convert to reference to 760 m or to bottom whichever is shallower
    disp('reference velocity')
for it = 1:nt
    for iy = 1:nx_ref-1
        kzi = find(~isnan(vel_ref(:,iy,it)),1,'last');
        if isempty(kzi)
            kzi = nz_ref;
        end
        vel_ref(:,iy,it) = vel_ref(:,iy,it) - vel_ref(kzi,iy,it);
    end
end

save -append Output.mat vel_ref lat_vel lon_vel depth_ref dh0_ref
disp('Done!')
clear

%%  Calculate Transport
%Input: depth_ref, temp_ref, lon_ref, lat_ref, vel_ref, lat_vel, lon_vel
%Output: GStransp, GSspd_z, GStransp_z; 
disp('Calc transport')

load('Output.mat')

%Define northern and southern boundary for transport calculations in degrees
BLim = [36.5 38.5];

[nz_ref,nx_ref,nt] = size(temp_ref);

% compute the area of each grid cell
 delt_z = diff(depth_ref(1:2));
 Area_GS = nan(nz_ref,nx_ref-1,nt);
 for i = 1:nt
     delt_y=deg2km(distance(lat_ref(1:nx_ref-1),lon_ref(1:nx_ref-1,i),lat_ref(2:nx_ref),lon_ref(2:nx_ref,i),'degrees'))*1000;
     Area_GS(:,:,i) = repmat(delt_y',[nz_ref 1])*delt_z;
 end

% Determine the maximum velocity at GS core (GSspd),
% and GS location (GSlat, GSlon) and volume transport
GSspd_z = nan(nz_ref,nt);
GStransp_z = nan(nz_ref,nt);
GSspd = nan(nt,1);
GSlat = GSspd;
GSlon = GSspd;
GStransp = GSspd;
GStransp_alt=GStransp;
indx_GS = find(lat_vel>=BLim(1) & lat_vel<=BLim(2));  %GS latitude range
for i = 1:nt
    transp = sum(vel_ref(:,:,i).*Area_GS(:,:,i),1,'omitmissing');
    vel0 = mean(vel_ref(1:25,:,i),1); % use top 50-m averaged velocity to determine GS core
    [~,j] = max(vel0(indx_GS)); % based on near surface max velocity
    indx0 = indx_GS(j);
    GSspd(i) = vel_ref(1,indx0,i);
    GSlat(i) = lat_vel(indx0);
    GSlon(i) = lon_vel(indx0,i);
    % find GS northern and southern boundary (where current changes sign)
    k1 = find(transp(1:indx0-1)<0,1,'last'); %southern boundary based on transport
    if isempty(k1)
        k1 = find(vel0(1:indx0-1)<0,1,'last'); %southern boundary based on velocity
    end
    k2 = find(transp(indx0:nx_ref-1)<0,1,'first');
    if isempty(k2)
        k2 = find(vel0(indx0:nx_ref-1)<0,1,'first');
    end
    k2 = k2 + indx0 - 1;
    GStransp(i) = sum(transp(k1+1:k2-1));
    GStransp_z(:,i) = sum(vel_ref(:,k1+1:k2-1,i).*Area_GS(:,k1+1:k2-1,i)/delt_z,2,'omitmissing');
    GSspd_z(:,i) = max(vel_ref(:,k1+1:k2-1,i),[],2);

end

save Output.mat -append GSlat GSlon GStransp 
disp('Done!')
clear

%% seasonal cycle
%Input: time_avg, GStransp_alt GSlat_alt GSspd_alt
%Output: mGSlat, mGSlat_alt, mGStransp, mGStransp_alt
disp('Calc seasonal means')
load Output.mat
[y_xbt,m_xbt,d_xbt] = datevec(time_avg);
Yday = datenum(0,m_xbt,d_xbt);
DayM = nan(12,1);
mGSlat = nan(12,2);
mGStransp = nan(12,2);
Yday2 = [Yday-366; Yday; Yday+366];
GSlat2 = [GSlat; GSlat; GSlat];
GStransp2 = [GStransp; GStransp; GStransp];
for i = 1:12
    DayM(i) = mean(datenum(0,i,1):datenum(0,i+1,1));
    ki = find(abs(Yday2-DayM(i))<=45);
    mGSlat(i,1) = mean(GSlat2(ki));
    mGSlat(i,2) = std(GSlat2(ki));
    mGStransp(i,1) = mean(GStransp2(ki));
    mGStransp(i,2) = std(GStransp2(ki));
 end
save Output.mat -append mGStransp mGSlat
clear

%%  set colormaps
disp('Plotting section')
%addpath /phodnet/share/mgoes/matlab/mapcolor/ 

clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
clr_S = cmocean('haline'); clr_S(1,:) = [.7 .7 .7];
clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
xyclr = [.4 .4 .4];
lon_lim = [-78 -61];                                   %Define lon lim for map
lat_lim = [31 44];                                     %Define lat lim for map

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
load Output.mat

label_lat = {'34\circN','','36\circN','','38\circN','','40\circN'};  %Define horizontal Label ticks

fig1 = myFigSize(1,10,6); clf % papersize (10,6)
subplot(2,3,1);hold off
m_proj('miller','lon',lon_lim,'lat',lat_lim);
[CS,CH]=m_contourf(lonb,latb,depth_btm,[-7000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
colormap(gca,clr_topo(50:end,:))
m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on
for i = 1:length(npf)
    m_plot(lon_xbt{i},lat_xbt{i},'.r','markersize',4)
end
m_plot(mean(lon_ref,2),lat_ref,'k','linewidth',1)
m_grid('tickdir','in','color',xyclr)
set(gca,'Position',[.01 .56 .34 .36])
ax=m_contfbar(0.27,[.57 .94],CS,CH,'axfrac',.025);
set(ax,'colormap',clr_topo(50:end,:))
set(ax,'tickdir','in','fontsize',8,'ycolor',xyclr,'xcolor',xyclr)
title(ax,{'m'},'fontweight','normal','fontsize',8);
m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,'(a) AX32')

% histogram # of cruises (month)
subplot(2,3,2); hold off
[~,m] = datevec(time_avg);
h = histogram(m,0.5:1:12.5); hold on
set(h,'EdgeColor',[0.1294 0.1294 0.1294])
set(gca,'XLim',[0.5 12.5],'YLim',[0 12])
set(gca,'Position',[0.38 0.56 0.23 0.36])
text(0.5,12.8,'(b) Monthly distribution of cruises')
ylabel('number of cruises')

% mean temperature
subplot(2,3,4); hold off
mtemp = mean(temp_ref,3,'omitmissing');
tmp = mtemp;
tmp(~isnan(tmp)) = 26;
tmp(isnan(tmp)) = -5;
contourf(lat_ref,-depth_ref,tmp,-6:32:26); hold on
colormap(gca,clr_T);
[CS,CH] = contourf(lat_ref,-depth_ref,mtemp,2:1:26,'LineStyle','none');
clim([0 26])
plot(lat_xbt{102},0,'.k')
ax=m_contfbar(0.82,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[2 26],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_T)
title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.06 0.1 0.23 0.35])
ylabel('Depth (m)')
set(gca,'xtick',34:1:40,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lat_ref(1),50,'(d) Mean temperature')

% mean salinity
msalt = mean(salt_ref,3,'omitmissing');
msalt(msalt<32.5) = 32.5;
tmp = msalt;
tmp(~isnan(tmp)) = 37;
tmp(isnan(tmp)) = 27;
subplot(2,3,5); hold off
contourf(lat_ref,-depth_ref,tmp,22:15:37); hold on
colormap(gca,clr_S);
[CS,CH] = contourf(lat_ref,-depth_ref,msalt,32:0.2:37,'LineStyle','none');
clim([30 37])
ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[32 37],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_S)
title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.38 0.1 0.23 0.35])
set(gca,'xtick',34:1:40,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lat_ref(1),50,'(e) Mean salinity')

% mean velocity
mvel = mean(vel_ref,3,'omitmissing');
tmp = mvel;
tmp(~isnan(tmp)) = 0.65;
tmp(isnan(tmp)) = -0.65;
subplot(2,3,6); hold off
contourf(lat_vel,-depth_ref,tmp,-0.7:1.4:0.7); hold on
colormap(gca,clr_V2);
[CS,CH] = contourf(lat_vel,-depth_ref,mvel,-0.6:0.05:0.6,'LineStyle','none');
contour(lat_vel,-depth_ref,mvel,[0 0],'color',[.7 .7 .7]);
clim([-0.65 0.65])
ax=m_contfbar(1.2,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[-0.6 0.6],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',-0.6:0.2:0.6,'yticklabel',-0.6:0.2:0.6)
set(ax,'colormap',clr_V2)
title(ax,'m/s','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.7 0.1 0.23 0.35])
set(gca,'xtick',34:1:40,'xticklabel',label_lat,'XTickLabelRotation',0)
text(lat_ref(1),50,'(f) Mean velocity')

%PIE
%CREATE X,Y AXES
timeYMD = datestr(time_avg','yyyymmdd');
data1 = str2num(timeYMD(:,5:6))'; %Month numbers
data2 = str2num(timeYMD(:,1:4))'; %year numbers

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
[g,col] = make_pie_subplot(axis2,axis1,[],leg1,data2,data1,C,N,[.62 .51 .345 .415]);
title(col,'# prof')
print -dpng -r400 figure1_ax32.png

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
