%%
clear
%Packages needed: GSW, M_MAP, cmocean, TS_PACKAGE
addpath /datasets/work/soop-xbt/work/UOT/programs/TS-lookup/datafiles/TS_PACKAGE_THACKER_noyear_globe_ftp/
addpath(genpath('/datasets/work/soop-xbt/work/cowley/matlab_localtools/toolbox/local/gsw_toolbox/'))
addpath /datasets/work/soop-xbt/work/UOT/programs/MATLAB/seawater_ver3_3/
%% load the data
% should be in this format:
  % Name           Size                Bytes  Class    Attributes
  % 
  % depth_xbt      1x450                3600  int64              
  % lat_xbt        1x139              103392  cell               
  % lon_xbt        1x139              103392  cell               
  % npf            1x139                1112  int64              
  % temp_xbt       1x139            40035656  cell               
  % time_xbt       1x139              103392  cell               
       
load PX30_data_cleaned.mat
% change depth_XBT to a double
depth_xbt = double(depth_xbt);
%% functions setup
function out = interp_fun(depths, temps, depth_ref)
    if numel(depths) > 1
        out = interp1(depths, temps, depth_ref, 'linear', NaN);
        % fill missing data near the surface
        inan = isnan(out(1:5));
        if sum(inan) > 0 && sum(inan) < 5
            % find the first instance of non-nan value and fill to surface
            igood = find(~inan);
            out(find(inan)) = out(igood(1));
        end
    else
        out = nan(size(depth_ref));
    end
end

function [groups,npf] = findTimestampGroups(timestamps, maxDiff)
    % findTimestampGroups - Identify groups of timestamps with differences < maxDiff
    %
    % Syntax: groups = findTimestampGroups(timestamps, maxDiff)
    %
    % Inputs:
    %   timestamps - Array of timestamps (datetime, datenum, or numeric)
    %   maxDiff    - Maximum time difference for grouping (default: 2 days)
    %
    % Outputs:
    %   groups - Cell array where each cell contains indices of grouped timestamps
    
    if nargin < 2
        maxDiff = 2; % Default: 2 days
    end
    
    % Sort timestamps and keep track of original indices
    [sortedTimestamps, sortIdx] = sort(timestamps);
    n = length(sortedTimestamps);
    
    % Initialize variables
    currentGroup = 1;
    groupAssignment = ones(n, 1);
    
    % Find groups based on time differences
    for i = 2:n
        timeDiff = days(sortedTimestamps(i) - sortedTimestamps(i-1));
        if timeDiff > maxDiff
            currentGroup = currentGroup + 1;
        end
        groupAssignment(i) = currentGroup;
    end
    
    % Create output groups with original indices
    numGroups = max(groupAssignment);
    groups = cell(numGroups, 1);
    npf = zeros(1,numGroups);
    for g = 1:numGroups
        groupMembers = find(groupAssignment == g);
        groups{g} = sort(sortIdx(groupMembers));
        npf(g) = length(groupMembers);
    end
end
%% create the grid and tidy
% grid onto depths by grouping by uniqueid
% [G, station_idx] = findgroups(uniqueid);
% 
% TEMP_cell = splitapply(@(d, t) {interp_fun(d, t, depth_xbt)}, ...
%     DEPTH, TEMP, G);
% 
% TEMP_grid = cell2mat(TEMP_cell');

%% clear some variables
% clear DEPTH LATITUDE LONGITUDE Station_ID TEMP TIME Callsign Cruise_ID ic
%% index each transect based on days between drops
% > 2days = new transect
% [grps, np] = findTimestampGroups(time, 2);
% [utrans,ib,ic] = unique(transect_number);
% npf = NaN*ones(1,length(utrans));
% for a = 1:length(npf)
%     npf(a) = sum(ic == a);
% end
%% focus on only data between 153 and 157.5 longitude and -26.75 to -25.75
% remove whole transects that have data outside the latitude bounds and any
% that have no data inside
lon_ref_30 = (153.13:0.1:166.83)';
lon_ref_31 = (166.93:0.1:178.33)';
lon_ref = [lon_ref_30;lon_ref_31];
xp = [153.1 153.1 165 178.4 178.4];
yp = [-22 -29 -26 -20 -14 ];
for a = 1:length(npf)
    lon = lon_xbt{a};
    lat = lat_xbt{a};
    % remove casts ti = time_xbt{a};outside specified region
    J = inpolygon(lon,lat,xp ,yp);
    if any(~J)
        ti = time_xbt{a};
        te = temp_xbt{a};
        lon(~J) = [];
        lat(~J)=[];
        ti(~J) = [];
        te(:,~J) = [];
        lon_xbt{a} = lon;
        lat_xbt{a} = lat;
        temp_xbt{a} = te;
        time_xbt{a} = ti;
        npf(a) = sum(J);
    end
end
% remove any transects with less than 10 profiles
irem = npf < 10;
lon_xbt(irem) = [];
lat_xbt(irem) = [];
temp_xbt(irem) = [];
time_xbt(irem) = [];
npf(irem) = [];

%%
%Input: depth,temp,xbt_lat,xbt_lon,
%Output: time_xbt temp_ref lon_ref lat_ref depth_ref npf LONGITUDE LATITUDE
disp('calc reference section')

% grid to reference section [lat_ref & lon_ref] and every 10 dbar [pres_ref]
% split the lon_ref for half-voyages
kz = find(depth_xbt<=760);
depth_ref = depth_xbt(kz);
nx_ref = length(lon_ref);
nz_ref = length(depth_ref);
nt = length(npf);
nz = length(depth_ref);

% find the reference latitude based on averaged latitude from all transect after removing outliers
lat_ref = nan(nx_ref,nt);
% fill in gaps, then interpolate to reference longitude (lon_ref)
temp_ref = nan(length(depth_ref),nx_ref,nt);
for i=1:nt
    % x values must be unique for interp1 therefore, remove any duplicates
    % by removing the dup profile

    tempi = temp_xbt{i};
    lati = lat_xbt{i};
    loni = lon_xbt{i};

    [loni, ie, id] = unique(loni,'stable');
    if length(id) ~= length(ie)
        disp(i)
        lati = lati(ie);
        tempi = tempi(:,ie);
    end
    if length(loni) < 2
        continue
    end
    lat_ref(:,i) = interp1(loni,lati,lon_ref,'linear');

    % if two profile are close together, which may induce anomous large velocity
    %   average to mid-longitude
    % dist_x = deg2km(distance(loni(1:end-1),lati(1:end-1),loni(2:end),lati(2:end),'degrees'));
    dist_x = sw_dist(lati,loni,'km');
    ik = find(dist_x<10);
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
        if length(ki) < 2
            continue
        end
        if length(ki) >= length(loni)*0.8 & length(ki)<length(loni)
            tempi(iz,:) = interp1(loni(ki),tempi(iz,ki),loni);
        end
        % interpolate to reference grid
        ki = find(~isnan(tempi(iz,:)));
        temp_ref(iz,:,i) = interp1(loni(ki),tempi(iz,ki),lon_ref)';
    end
    % backfill with NaNs where the range of longitude does not extend for
    % full lon_ref
    i31 = max(loni) >= min(lon_ref_31);
    i30 = min(loni) <= max(lon_ref_30);
    if ~(i31 && i30)
        if i30
            %just first part px30
            lat_ref(length(lon_ref_30)+1:end,i) = NaN;
            temp_ref(:, length(lon_ref_30)+1:end,i) = NaN;
        end
        if i31
                        %just px31
            lat_ref(1:length(lon_ref_30),i) = NaN;
            temp_ref(:,1:length(lon_ref_30),i) = NaN;
        end
    end
end
dlat=(lat_ref-mean(lat_ref,2,'omitmissing'));
kgood=find(mean(abs(dlat),1,'omitmissing')<=0.25);
mean_lat_ref=mean(lat_ref(:,kgood),2,'omitmissing');

%save('Output.mat','time_xbt','temp_ref','lon_ref','lat_ref','depth_ref')
save Output.mat time_xbt temp_ref lon_ref lat_ref depth_ref npf lon_xbt lat_xbt mean_lat_ref
disp('Done!')
clear 
%% =======================Marlos Salinity===========
%Input:time_xbt,temp_ref,lon_ref,lat_ref,depth_ref
%Output:salt_ref (Goes method),time_avg 
disp('calc sal')

load Output.mat
nt = length(npf);
ny_ref = length(lon_ref);
nz_ref = length(depth_ref);
time_avg = nan(nt,1);  %averaged observation time for the trensect
for it = 1:nt
     time_avg(it) = mean(time_xbt{it},'omitmissing');
end
timeYMD = datestr(time_avg','yyyymmdd');
time2 = str2num(timeYMD);
time2 = time2(:,ones(1,ny_ref))';

%Calc Salinity - All at once
auxT = reshape(temp_ref,[nz_ref ny_ref*nt]);
auxLat = reshape(lat_ref, [1 ny_ref*nt]);
auxLon = reshape(lon_ref(:,ones(1,nt)),[1 ny_ref*nt]);
auxTime = reshape(time2,[1 ny_ref*nt]);
salt_ref = Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(auxT,depth_ref,auxLat,auxLon,auxTime,0,depth_ref,'goes');

%Reshape back to initial size
salt_ref = reshape(salt_ref,[nz_ref ny_ref nt]);

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
% use mean_lat_ref for PX30
for ix = 1:nx_ref
    pres_ref(:,ix) = gsw_p_from_z(-depth_ref,mean_lat_ref(ix));
end
vel_ref = nan(nz_ref,nx_ref-1,nt);
lat_vel = nan(nx_ref-1,nt);
dh0_ref = nan(nx_ref,nt);
for it = 1:nt
    [SA, in_ocean] = gsw_SA_from_SP(salt_ref(:,:,it),pres_ref,lon_ref,lat_ref(:,it));
    CT = gsw_CT_from_t(SA,temp_ref(:,:,it),pres_ref);
    ga = gsw_geo_strf_dyn_height(SA,CT,pres_ref,0);
    [vel0,mid_lat,mid_lon] = gsw_geostrophic_velocity(ga,lon_ref,lat_ref(:,it),pres_ref);
    lat_vel(:,it) = mid_lat(1,:)';
    lon_vel =  mid_lon(1,:)';
    vel_ref(:,:,it) = vel0;
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

%Define eastern and western boundary for transport calculations in degrees
BLim = [153.1 158.5];

[nz_ref,ny_ref,nt] = size(temp_ref);

% compute the area of each grid cell
 delt_z = diff(depth_ref(1:2));
 Area_GS = nan(nz_ref,ny_ref-1,nt);
 for i = 1:nt
     % delt_y=deg2km(distance(lat_ref(1:ny_ref-1,i),lon_ref(1:ny_ref-1),lat_ref(2:ny_ref,i),lon_ref(2:ny_ref),'degrees'))*1000;
    delt_y = sw_dist(lat_ref(:,i),lon_ref,'km')*1000;  %In meters
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
indx_GS = find(lon_vel>=BLim(1) & lon_vel<=BLim(2));  %GS longitude range
for i = 1:nt
    transp = sum(vel_ref(:,:,i).*Area_GS(:,:,i),1,'omitmissing');
    vel0 = mean(vel_ref(1:25,:,i),1,'omitmissing'); % use top 50-m averaged velocity to determine GS core
    [~,j] = min(vel0(indx_GS)); % based on near surface max velocity
    indx0 = indx_GS(j);
    GSspd(i) = vel_ref(1,indx0,i);
    GSlat(i) = lat_vel(indx0,i);
    GSlon(i) = lon_vel(indx0);
    % find GS eastern and western boundary (where current changes sign)
    k1 = find(transp(1:indx0-1)>0,1,'last'); %southern boundary based on transport

    if isempty(k1)
        k1 = find(vel0(1:indx0-1)>0,1,'last'); %southern boundary based on velocity
    end
    if isempty(k1), k1=0;end
    k2 = find(transp(indx0:ny_ref-1)>0,1,'first');
    if isempty(k2)
        k2 = find(vel0(indx0:ny_ref-1)>0,1,'first');
    end
    k2 = k2 + indx0 - 1;
    if k2>k1+1
        GStransp(i) = sum(transp(k1+1:k2-1));
        GStransp_z(:,i) = sum(vel_ref(:,k1+1:k2-1,i).*Area_GS(:,k1+1:k2-1,i)/delt_z,2,'omitmissing');
        GSspd_z(:,i) = min(vel_ref(:,k1+1:k2-1,i),[],2);
    end
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
mGSlon = nan(12,2);
mGStransp = nan(12,2);
Yday2 = [Yday-366; Yday; Yday+366];
GSlon2 = [GSlon; GSlon; GSlon];
GStransp2 = [GStransp; GStransp; GStransp];
for i = 1:12
    DayM(i) = mean(datenum(0,i,1):datenum(0,i+1,1));
    ki = find(abs(Yday2-DayM(i))<=45);
    mGSlon(i,1) = mean(GSlon2(ki));
    mGSlon(i,2) = std(GSlon2(ki));
    mGStransp(i,1) = mean(GStransp2(ki));
    mGStransp(i,2) = std(GStransp2(ki));
 end
save Output.mat -append mGStransp mGSlon
clear

%%  set colormaps
disp('Plotting section')
addpath /datasets/work/soop-xbt/work/UOT/programs/xbt-products/figures_code/ 
addpath(genpath('/datasets/work/soop-xbt/work/UOT/programs/xbt-products/figures_code/m_map/'))
clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
clr_S = cmocean('haline'); clr_S(1,:) = [.7 .7 .7];
clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
xyclr = [.4 .4 .4];
lon_lim = [151 180];                                   %Define lon lim for map
lat_lim = [-28 -14];                                     %Define lat lim for map

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

tick_lon = 155:10:175;
label_lon = strcat(num2str(tick_lon'),'\circE');

fig1 = myFigSize(1,10,6); clf % papersize (10,6)
subplot(2,3,1);hold off
m_proj('miller','lon',lon_lim,'lat',lat_lim);
[CS,CH]=m_contourf(lonb,latb,depth_btm,[-7000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
colormap(gca,clr_topo(50:end,:))
m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on

for i = 1:length(npf)
    m_plot(lon_xbt{i},lat_xbt{i},'.r','markersize',4)
end

m_plot(lon_ref,mean(lat_ref,2,'omitmissing'),'k','linewidth',1)
m_grid('tickdir','in','color',xyclr)
set(gca,'Position',[.05 .56 .24 .36])
ax=m_contfbar(0.27,[.57 .94],CS,CH,'axfrac',.025);
set(ax,'colormap',clr_topo(50:end,:))
set(ax,'tickdir','out','fontsize',8,'ycolor',xyclr,'xcolor',xyclr)
title(ax,{'m'},'fontweight','normal','fontsize',8);
m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,'(a) PX30')
axis normal

% histogram # of cruises (month)
subplot(2,3,2); hold off
[~,m] = datevec(time_avg);
h = histogram(m,0.5:1:12.5); hold on
set(h,'EdgeColor',[0.1294 0.1294 0.1294])
set(gca,'XLim',[0.5 12.5],'YLim',[0 16])
set(gca,'Position',[0.38 0.56 0.23 0.36])
y_lim = ylim;
text(0.5,y_lim(2)+diff(y_lim)/15,'(b) Monthly distribution of cruises')
ylabel('number of cruises')

% mean temperature
subplot(2,3,4); hold off
mtemp = mean(temp_ref,3,'omitmissing');
tmp = mtemp;
tmp(~isnan(tmp)) = 26;
tmp(isnan(tmp)) = -5;

[CS,CH] = contourf(lon_ref,-depth_ref,mtemp,2:1:26,'LineStyle','none');
clim([0 26])
colormap(gca,clr_T);
set(gca,'color',[.7 .7 .7])
hold on
plot(lon_xbt{1},0,'.k')
ax=m_contfbar(0.82,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[2 26],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_T)
title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.06 0.1 0.23 0.35])
ylabel('Depth (m)')
set(gca,'xtick',tick_lon,'xticklabel',label_lon,'XTickLabelRotation',0)
text(lon_ref(1),50,'(d) Mean temperature')

% mean salinity
msalt = mean(salt_ref,3,'omitmissing');
msalt(msalt<32.5) = 32.5;
tmp = msalt;
tmp(~isnan(tmp)) = 37;
tmp(isnan(tmp)) = 27;
subplot(2,3,5); hold off
[CS,CH] = contourf(lon_ref,-depth_ref,msalt,32:0.2:37,'LineStyle','none');
clim([30 37])
colormap(gca,clr_S);
set(gca,'color',[.7 .7 .7])
ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[32 37],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',clr_S)
title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.38 0.1 0.23 0.35])
set(gca,'xtick',tick_lon,'xticklabel',label_lon,'XTickLabelRotation',0)
text(lon_ref(1),50,'(e) Mean salinity')

% mean velocity
mvel = mean(vel_ref,3,'omitmissing');
tmp = mvel;
tmp(~isnan(tmp)) = 0.65;
tmp(isnan(tmp)) = -0.65;
subplot(2,3,6); hold off
colormap(gca,clr_V2);
[CS,CH] = contourf(lon_vel,-depth_ref,mvel,-0.8:0.05:0.8,'LineStyle','none');
% contour(lon_vel,-depth_ref,mvel,[0 0],'color',[.7 .7 .7]);
clim([-0.8 0.8])
colormap(gca,clr_V2);
set(gca,'color',[.7 .7 .7])
ax=m_contfbar(1.2,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','in','ylim',[-0.85 0.85],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',-0.8:0.2:0.8,'yticklabel',-0.8:0.2:0.8)
set(ax,'colormap',clr_V2)
title(ax,'m/s','FontSize',10,'color',xyclr,'FontWeight','normal')
set(gca,'Position',[0.7 0.1 0.23 0.35])
xlim([153.1800 160])
tick_lon = 154:2:160;
label_lon = strcat(num2str(tick_lon'),'\circE');
set(gca,'xtick',tick_lon,'xticklabel',label_lon,'XTickLabelRotation',0)
text(lon_ref(1),50,'(f) Mean velocity')

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
C = double(npf'); %number of profiles/transect

%MAKE PIE PLOTS
%1: month, year
%clf
[g,col] = make_pie_subplot(axis2,axis1,[],leg1,data2,data1,C,N,[.62 .51 .345 .415]);
title(col,'# prof')
print -dpng -r400 figure1_px30.png

%%
%Redo auxiliary month dates
[y_xbt,m_xbt,d_xbt] = datevec(time_avg);
Yday = datenum(0,m_xbt,d_xbt);
for i = 1:12
   DayM(i) = mean(datenum(0,i,1):datenum(0,i+1,1));
end


tick_lon = 153:2:160;
label_GSpos = strcat(num2str(tick_lon'),'\circE');
% label_GSpos = {'153\circN','','154\circN','','155\circN','','156\circN'};
fig2 = myFigSize(2,10,6); clf

subplot(2,2,1); hold off
plot(Yday,GSlon,'.k'); hold on
plot(DayM,mGSlon(:,1),'Color',[0 .45 .75],'LineWidth',1.5)
set(gca,'xlim',[0 366],'ylim',[153 160],'ytick',tick_lon,'YTickLabel',label_GSpos)
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'position',[0.1 0.56 0.36 0.35])
text(0,160.5,'(a) EAC position (seasonal)')

subplot(2,2,2); hold off
plot(Yday,GStransp/1e6,'.k'); hold on
plot(DayM,mGStransp(:,1)/1e6,'Color',[0 .45 .75],'LineWidth',1.5)
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'xlim',[0 366])%,'ylim',[-40 0])
set(gca,'position',[0.56 0.56 0.36 0.35])
text(0,2,'(b) EAC transport (seasonal)')

subplot(2,2,3);hold off
plot(time_avg,GSlon,'color',[0 .45 .75],'LineWidth',1);hold on
%set(gca,'ylim',[36.5 38.7], ...
   set(gca,'xlim',[time_avg(1)-30 time_avg(end)+30])
set(gca,'xtick',datenum(1992:5:2022,1,1),'xticklabel',1992:5:2022,'XTickLabelRotation',0)
set(gca,'ytick',tick_lon,'YTickLabel',label_GSpos)
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'position',[0.1 0.1 0.36 0.35])
text(time_avg(1)-30,159.2,'(c) EAC position')

subplot(2,2,4);hold off
plot(time_avg,GStransp/1e6,'color',[0 .45 .75],'LineWidth',1);hold on
ylabel('Transport (Sv)')
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca,'ylim',[-40 0],'xlim',[time_avg(1)-30 time_avg(end)+30])
set(gca,'yaxislocation','left','xaxislocation','bottom')%,'color','none')
set(gca,'xtick',datenum(1992:5:2022,1,1),'xticklabel',1992:5:2022,'XTickLabelRotation',0)
set(gca,'position',[0.56 0.1 0.36 0.35])
text(time_avg(1)-30,1.774,'(d) EAC transport')

print -dpng -r400 figure2_px30.png


%% Hovmoller of dynamic height
fig3 = myFigSize(3,5,8); clf
tick_lon = 155:5:175;
label_lon = strcat(num2str(tick_lon'),'\circE');

nt = length(time_avg);
time_label=str2num(datestr(time_avg,'yyyy'));
%THIS USES THE REAL AXIS
%contourf(lat_ref,time_avg,dh0_ref',[0.6:.1:1.8]); hold on
%set(gca,'ytick',datenum(2010:2:2024,1,1),'yticklabel',2010:2:2024,'XTickLabelRotation',0)

%THIS USES A FAKE AXIS
contourf(lon_ref,1:nt,dh0_ref',[0.6:.1:1.8],'edgecolor','none'); hold on
set(gca,'ytick',10:10:nt,'yticklabel',time_label(10:10:nt))
%
colormap(gca,clr_V2);
xlim([154 177])
set(gca,'xtick',tick_lon,'xticklabel',label_lon,'XTickLabelRotation',0)
set(gca,'position',[0.13 0.2 0.7 0.7])
ax = colorbar('position',[.85 .2 .03 .7]);
set(ax,'tickdir','in','fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',1.3:0.1:1.8,'yticklabel',1.3:0.1:1.8)
set(ax,'colormap',clr_V2)
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
title(ax,'m','FontSize',10,'color',xyclr,'FontWeight','normal')
title('Surface Dynamic Height ')
print -dpng -r400 figure3_px30.png
