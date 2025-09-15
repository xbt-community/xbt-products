addpath('/data2/Movar_analises/nc_final/Figures_Shenfu/Figures_code/')

clear all, close all, clc

load /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/output_ax97.mat
temp_ref1=ncread('/data2/Movar_analises/arquivos_pro_sal/NETCDF/AX97_MOVAR.nc','CT');
temp_ref=permute(temp_ref1(1:end-1,1:end-1,:),[3 2 1]); clear temp_ref1
salt_ref1  =ncread('/data2/Movar_analises/arquivos_pro_sal/NETCDF/AX97_MOVAR.nc','SA');
salt_ref=permute(salt_ref1(1:end-1,1:end-1,:),[3 2 1]);  clear salt_ref1


%% compute velocity (BCW toolbox)
% first compute pressure from depth at each grid
%load /phodnet/share/dong/gulfstream/oleander/temp_ax32_gridded.mat
%Input: temp_ref,depth_ref,salt_ref,lat_ref,lon_ref
%Output:vel_ref lat_vel lon_vel depth_ref dh0_ref
% disp('compute velocity. Takes some time to finish')
% addpath ../PX30
% load('Output.mat')

[nz_ref,nx_ref,nt] = size(temp_ref); 
[~,ga_refInd] = min(abs(depth_ref-800));  %Reference for dynamic height output (760-800 m)
ref_dep = 500%depth_ref(ga_refInd);
pres_ref = nan(nz_ref,nx_ref);

% use mean_lon_ref for AX32
for ix = 1:nx_ref
    pres_ref(:,ix) = gsw_p_from_z(-depth_ref,lon_ref(ix));
end
vel_ref = nan(nz_ref,nx_ref-1,nt);
lon_vel = nan(nx_ref-1);
dh0_ref = nan(nx_ref,nt);
for it = 1:nt
    [SA, in_ocean] = gsw_SA_from_SP(salt_ref(:,:,it),pres_ref,lon_ref,lat_ref);
    CT = gsw_CT_from_t(SA,temp_ref(:,:,it),pres_ref);
    ga = gsw_geo_strf_dyn_height(SA,CT,pres_ref,0);

    inan = ~isnan(ga(1,:));
    DH_abs=ga*0;
    [DH_abs(:,inan),addep] = add_abs_gpan(ga(:,inan),lat_ref(inan),lon_ref(inan),depth_ref,ref_dep);
     ga = DH_abs;
    [vel0,mid_lat,mid_lon] = gsw_geostrophic_velocity(ga,lon_ref,lat_ref,pres_ref);
    lat_vel = mid_lat(1,:)';
    lon_vel =  mid_lon(1,:)';
    vel_ref(:,:,it) = vel0;
    dh0_ref(:,it) = ga(1,:)/9.81; %Convert from m2/s2 to meter at surface

end

save /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/output_ax97.mat -append vel_ref lat_vel lon_vel depth_ref dh0_ref salt_ref temp_ref
% disp('Done!')
% clear

%%  Calculate Transport
%Input: depth_ref, temp_ref, lon_ref, lat_ref, vel_ref, lat_vel, lon_vel
%Output: BCtransp, BCspd_z, BCtransp_z; 
% disp('Calc transport')

clear all, close all, clc

load /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/output_ax97.mat

% load('Output_ax97.mat')

%Define northern and southern boundary for transport calculations in degrees
BLim = [-40.9 -39];

[nz_ref,nx_ref,nt] = size(temp_ref);

% compute the area of each grid cell
 delt_z = diff(depth_ref(1:2));
 Area_BC = nan(nz_ref,nx_ref-1);
%  for i = 1:nt
     delt_y=deg2km(distance(lat_ref(1:nx_ref-1),lon_ref(1:nx_ref-1),lat_ref(2:nx_ref),lon_ref(2:nx_ref),'degrees'))*1000;
     Area_BC(:,:) = repmat(delt_y',[nz_ref 1])*delt_z;
%  end

% Determine the maximum velocity at BC core (BCspd),
% and BC location (BClat, BClon) and volume transport
BCspd_z = nan(nz_ref,nt);
BCtransp_z = nan(nz_ref,nt);
BCspd = nan(nt,1);
BClat = BCspd;
BClon = BCspd;
BCtransp = BCspd;
BCtransp_alt=BCtransp;
indx_BC = find(lon_vel>=BLim(1) & lon_vel<=BLim(2));  %BC latitude range
for i = 1:nt
    transp = sum(vel_ref(1:51,:,i).*Area_BC(1:51,:),1,'omitnan');
    vel0 = mean(vel_ref(1:6,:,i),1); % use top 50-m averaged velocity to determine BC core
    [~,j] = min(vel0(indx_BC)); % based on near surface max velocity
    indx0 = indx_BC(j);
    BCspd(i) = vel_ref(1,indx0,i);
    BClat(i) = lat_vel(indx0);
    BClon(i) = lon_vel(indx0);
    % find BC northern and southern boundary (where current changes sign)
    k1 = find(transp(1:indx0-1)>0,1,'last'); %southern boundary based on transport
    if isempty(k1)
        k1 = find(vel0(1:indx0-1)>0,1,'last'); %southern boundary based on velocity
    end
    if isempty(k1), k1 = 0; end
    
    k2 = find(transp(indx0:nx_ref-1)>0,1,'first');
    if isempty(k2)
        k2 = find(vel0(indx0:nx_ref-1)>0,1,'first');
    end
    k2 = k2 + indx0 - 1;
    BCtransp(i) = sum(transp(k1+1:k2-1));
    BCtransp_z(:,i) = sum(vel_ref(:,k1+1:k2-1,i).*Area_BC(:,k1+1:k2-1)/delt_z,2,'omitnan');
%     BCspd_z(:,i) = max(vel_ref(:,k1+1:k2-1,i),[],2);

end

save /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/output_ax97.mat -append -append BClat BClon BCtransp 
disp('Done!')
clear


%%  set colormaps
disp('Plotting section')
%addpath /phodnet/share/mgoes/matlab/mapcolor/ 
addpath('/data2/tayannepf/toolbox_movar/scripts_matlab/colorbar_yyw/')
addpath(genpath('/data2/tayannepf/toolbox_movar/scripts_matlab/'))
% addpath('/data2/tayannepf/toolbox_movar/scripts_matlab/m_map')

clr_T = jet;
clr_S = colormap_yyw('MPL_rainbow');
clr_V = colormap_yyw('GMT_polar');

clr_topo = m_colmap('blue',256);
clr_V2 = flipud(othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
%clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
%clr_S = cmocean('haline'); clr_S(1,:) = [.7 .7 .7];
%clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
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

% tick_lat = -42:2:-30;
% label_lat = strcat(num2str(tick_lat'),'\circW');
label_lat = {'42\circW','39\circW','36\circW','33\circW','30\circW'};  %Define horizontal Label ticks

%label_lat = {'34\circN','','36\circN','','38\circN','','40\circN'};  %Define horizontal Label ticks

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
%text(0.5,12.8,'(b) Monthly distribution of cruises')
y_lim = ylim;
text(0.5,y_lim(2)+diff(y_lim)/15,'(b) Monthly distribution of cruises')
ylabel('number of cruises')

% mean temperature
subplot(2,3,4); hold off
mtemp = mean(temp_ref,3,'omitnan');
tmp = squeeze(mtemp);
tmp(~isnan(tmp)) = 26;
tmp(isnan(tmp)) = -5;
contourf(lon_ref,-depth_ref,tmp(1:end,:),-6:32:26); hold on

[CS,CH] = contourf(lon_ref,-depth_ref,mtemp(:,:),2:1:26,'LineStyle','none'); hold on
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
set(gca,'color',[.7 .7 .7])

% mean salinity
msalt = mean(salt_ref,3,'omitnan');
msalt(msalt<32.5) = 32.5;
tmp = squeeze(msalt);
tmp(~isnan(tmp)) = 37;
tmp(isnan(tmp)) = 27;
subplot(2,3,5); hold off
set(gca,'color',[.7 .7 .7])

contourf(lon_ref,-depth_ref,tmp(1:end,:),22:15:37); hold on
[CS,CH] = contourf(lon_ref,-depth_ref,msalt(:,:),32:0.2:37,'LineStyle','none'); hold on
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
mvel = nanmean(vel_ref,3);
tmp = squeeze(mvel);
tmp(~isnan(tmp)) = 0.35;
tmp(isnan(tmp)) = -0.35;
subplot(2,3,6); hold off
contourf(lon_vel,-depth_ref,tmp,-0.7:1.4:0.7); hold on
colormap(gca,clr_V2);
[CS,CH] = contourf(lon_vel,-depth_ref,mvel,-0.20:0.01:0.20,'LineStyle','none');
contour(lon_vel,-depth_ref,mvel,[0 0],'color',[.7 .7 .7]);
caxis([-0.20 0.20])
ax=m_contfbar(1.2,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
% ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');

set(ax,'tickdir','in','ylim',[-0.20 0.20],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'ytick',-0.20:0.05:0.20,'yticklabel',-0.20:0.05:0.20)
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
% axis1 = [min(data1):max(data1)]';
% axis2 = [min(data2):max(data2)+1]';
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
title(col,'# prof')
%print -dpng -r400 figure1_ax32.png

print -dpng -r400 /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/figure1_ax97_new.png
%%
