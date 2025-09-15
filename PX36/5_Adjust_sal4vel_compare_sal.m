rmpath /phodnet/share/mgoes/TS_PACKAGE_THACKER_noyear_globe_ftp 
addpath /phodnet/share/mgoes/TS_PACKAGE_THACKER_noyear_globe 
close all
clear
cd('C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2\')
load Output.mat

z = depth_ref;
lat = lat_ref'; 
lon = lon_ref;   

 
time = time_avg; %202001*ones(3151,529); 
[nz,nx,nt] = size(temp_ref);

%Remove visible outliers in time 9
inan = true([1 145]);
inan([88 89 93]) = false;
t = temp_ref(:,:,9);
t(:,~inan) = nan;
s = salt_ref(:,:,9);
s(:,~inan) = nan;
t = interp1(lat_ref(inan),t(:,inan)',lat_ref);
s = interp1(lat_ref(inan),s(:,inan)',lat_ref);
temp_ref(:,:,9) = t';

%Slab layer approximation to remove nans from surfacce
for ii=1:nt    
    for jj=1:nx
       t = temp_ref(:,jj,ii);
      if isnan(t(1))
          oi = find(isnan(t));
          oi(oi>5)=[];
          t(oi) = t(oi(end)+1);
          temp_ref(:,jj,ii) = t;
      end
    end
end

%Load salinity and replace salt_ref
y1 = temp_ref*nan;y2=y1;
for ii=1:nt
    t = temp_ref(:,:,ii);
  [y1(:,:,ii)] = Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(t,z(:),lat,lon(:,ii)',time(ii*ones(nx,1)),0,[],'goes');
  [y2(:,:,ii)] = Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(t,z(:),lat,lon(:,ii)',time(ii*ones(nx,1)),0,[],'annual');
end

salt_ref = y1;  %replace salt_ref

%Remove salt where difference to annual is larger than 1
salt_ref((salt_ref-y2)>1) = y2((salt_ref-y2)>1);


%%
addpath GSW_2015/
addpath m_map/
addpath TS_PACKAGE_THACKER_noyear_globe/
addpath C:/Users/fanto/OneDrive - Uniparthenope/MATLAB/Seawater/
%addpath ../PX30

[nz_ref,nx_ref,nt] = size(temp_ref);
[~,ga_refInd] = min(abs(depth_ref-740));  %Reference for dynamic height output (740-800 m)
ref_dep = depth_ref(ga_refInd);
pres_ref = nan(nz_ref,nx_ref);

% use mean_lon_ref for AX32
for ix = 1:nx_ref
    pres_ref(:,ix) = gsw_p_from_z(-depth_ref,lat_ref(ix));
end
vel_ref = nan(nz_ref,nx_ref-1,nt);
lon_vel = nan(nx_ref-1,nt);
dh0_ref = nan(nx_ref,nt);
for it = 1:nt
    if it==25 
        continue 
    end
    % Comment GSW dynamic height. Use seawater instead
 %   [SA, in_ocean] = gsw_SA_from_SP(salt_ref(:,:,it),pres_ref,lon_ref(:,it),lat_ref);
 %   CT = gsw_CT_from_t(SA,temp_ref(:,:,it),pres_ref);
   % ga = gsw_geo_strf_dyn_height(SA,CT,pres_ref,0);
    [gpan]= sw_gpan(salt_ref(:,:,it),temp_ref(:,:,it),pres_ref);ga =-1*gpan;
    inan = ~isnan(ga(1,:));
    DH_abs=ga*0;
    [DH_abs(:,inan),addep] = add_abs_gpan(ga(:,inan),lat_ref(inan),lon_ref(inan,it),depth_ref,ref_dep);
     ga = DH_abs;
    [vel0,mid_lat,mid_lon] =  gsw_geostrophic_velocity(ga,lon_ref(:,it),lat_ref,pres_ref);%gsw
    vel0 = -1*sw_gvel(ga,lat_ref,lon_ref(:,it)); %replace seawater (optional)
    lat_vel = mid_lat(1,:)';
    lon_vel(:,it) =  mid_lon(1,:)';
    vel_ref(:,:,it) = vel0;
    dh0_ref(:,it) = ga(1,:)/9.81; %Convert from m2/s2 to meter at surface

end

%Copy Output to Output_new and append variables.

save('Output.mat','-append','dh0_ref','vel_ref','salt_ref','temp_ref')

%load('ssh4px36_all_cruises.mat')

clear
%% ======================= Plotting Section =======================
disp('Plotting section')
% Set up colormaps and limits for topography, temperature, salinity, and velocity
cd('C:/Users/fanto/OneDrive - Uniparthenope/Marlos_Goes_TS-lookup-master/V1_Scripts/V2/colorbar_yyw/');
clr_T = jet;
clr_S = colormap_yyw('MPL_rainbow');
clr_V = colormap_yyw('GMT_polar');

clr_topo = m_colmap('blue',256);
clr_V2 =flipud (othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
xyclr = [.4 .4 .4];
lon_lim = [165 185];  lat_lim = [-74 -44];
cd('C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2')
% Load and mask topography data
x = ncread('ETOPO1_Ice_g_gmt4.grd','x'); x = wrapTo360(x);
y = ncread('ETOPO1_Ice_g_gmt4.grd','y');
z = ncread('ETOPO1_Ice_g_gmt4.grd','z');
kx = find(x >= lon_lim(1)-0.25 & x <= lon_lim(2)+0.25);
ky = find(y >= lat_lim(1)-0.25 & y <= lat_lim(2)+0.25);
lonb = x(kx); latb = y(ky); depth_btm = z(kx,ky)';
depth_btm(depth_btm>0) = NaN;
clear kx ky x y z

%% Start Plot
load Output.mat

fig1 = myFigSize(1,10,6); clf % papersize (10,6)
subplot(2,3,1);hold off
m_proj('miller','lon',lon_lim,'lat',lat_lim);
[CS,CH]=m_contourf(lonb,latb,depth_btm,[-7000:1000:-1000 -500 -200 0],'edgecolor','none','LineStyle','none');
colormap(gca,clr_topo(50:end,:))
m_coast('patch',[.8 .8 .8],'edgecolor',[.8 .8 .8]);hold on
for i = 1:length(npf)
    lat_i = lat_xbt{i};
    lon_i = lon_xbt{i};

    % Applica il filtro
    idx = lat_i >= -72 & lat_i <= -46 & lon_i >= 170 & lon_i <= 180;

    % Plot solo i punti che soddisfano il filtro
    m_plot(lon_i(idx), lat_i(idx), '.r', 'markersize', 4);
end

m_plot(mean(lon_ref,2,'omitnan'),lat_ref,'k','linewidth',1)
m_grid('tickdir','out','color',xyclr,'xtick', 3, ...          % ogni 5° di longitudine
    'ytick', 10, ...          % ogni 5° di latitudine (opzionale)
    'fontsize', 8);
set(gca,'Position',[.01 .56 .34 .36])
ax=m_contfbar(0.27,[.57 .94],CS,CH,'axfrac',.025);
set(ax,'colormap',clr_topo(50:end,:))
set(ax,'tickdir','out','fontsize',8,'ycolor',xyclr,'xcolor',xyclr)
title(ax,{'m'},'fontweight','normal','fontsize',8);
m_text(lon_lim(1),lat_lim(2)+diff(lat_lim)/15,'(a) PX36')

% histogram # of cruises (month)
subplot(2,3,2); hold off
[~,m] = datevec(time_avg);
h = histogram(m,0.5:1:21.5); hold on
set(h,'EdgeColor',[0.1294 0.1294 0.1294])
set(gca,'XLim',[0.5 12.5],'YLim',[0 21])
set(gca,'Position',[0.38 0.56 0.23 0.36])
text(0.5, 1.05, '(b) Monthly distribution of cruises', ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'normal');
ylabel('number of cruises')

% mean temperature
subplot(2,3,4); hold off
mtemp = mean(temp_ref,3,'omitnan');
tmp = mtemp;
tmp(~isnan(tmp)) = 12;
tmp(isnan(tmp)) = -2;
contourf(lat_ref,-depth_ref,tmp,-2:0.5:12, 'LineStyle', 'none');
hold on
colormap(gca,clr_T);
[CS,CH] = contourf(lat_ref,-depth_ref,mtemp,-2:0.5:12,'LineStyle','none');
caxis([-3 12]) % Imposta scala colori tra -2 e 12
lat_all = vertcat(lat_xbt{4});
lat_all = lat_all(lat_all >= -72 & lat_all <= -46);
hold on
plot(lat_all,0,'.k')% Plotti tutte lat di xbt della prima campagna
ax = m_contfbar(0.82,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','out','ylim',[-3 12],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(gca,'Position',[0.06 0.1 0.23 0.35])
ylabel('Depth (m)')
set(ax,'colormap',clr_T)
title(ax,'\circC','FontSize',10,'color',xyclr,'FontWeight','normal')

% Imposta ticks con valori NEGATIVI (coerenti con lat_ref)
xticks(-72:2:-46)
% Labels con passo 2 (solo valori pari, senza segno meno)
vals = -72:2:-46;
labels = cell(1,length(vals));
for k = 1:length(vals)
    if mod(abs(vals(k)),4) == 0
        labels{k} = sprintf('%d\\circS', abs(vals(k)));
    else
        labels{k} = '';
    end
end
xticklabels(labels)

% Inverti l'asse X per mostrare latitudine in ordine crescente (da sud a
% nord)
set(gca,'XDir','reverse')

% Titolo sopra il grafico al centro
title('(d) Mean temperature', 'FontSize', 12, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.05, 0])



% mean salinity
msalt = median(salt_ref,3,'omitnan');
%msalt(msalt>32.5) = 32.5;   %% Perchè?
tmp = msalt;
%tmp(~isnan(tmp)) = 35;
%tmp(isnan(tmp)) = 32;
%tmp(1:,138:145)=33.66;

subplot(2,3,5); hold off
contourf(lat_ref,-depth_ref,tmp,32:0.1:35); hold on
colormap(gca,clr_S);
[CS,CH] = contourf(lat_ref,-depth_ref,msalt,32:0.1:35,'LineStyle','none');
clim([32 35])
ax=m_contfbar(1,[.05 0.9],CS,CH,'axfrac',.035,'endpiece','no');
set(ax,'tickdir','out','ylim',[32 35],'fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(gca,'Position',[0.38 0.1 0.23 0.35])
ylabel('Depth (m)')
set(ax,'colormap',clr_S)
title(ax,'psu','FontSize',10,'color',xyclr,'FontWeight','normal')

% Imposta ticks con valori NEGATIVI (coerenti con lat_ref)
xticks(-72:2:-46)
% Labels con passo 2 (solo valori pari, senza segno meno)
vals = -72:2:-46;
labels = cell(1,length(vals));
for k = 1:length(vals)
    if mod(abs(vals(k)),4) == 0
        labels{k} = sprintf('%d\\circS', abs(vals(k)));
    else
        labels{k} = '';
    end
end
xticklabels(labels)

% Inverti l'asse X per mostrare latitudine in ordine crescente (da sud a
% nord)
set(gca,'XDir','reverse')
% Titolo sopra il grafico al centro
title('(e) Mean salinity', 'FontSize', 12, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.05, 0])



% mean velocity
%vel_ref(vel_ref > 4) = NaN;

%vel_ref(vel_ref < -4) = NaN;


mvel = mean(vel_ref, 3, 'omitnan');
tmp = mvel;
tmp(~isnan(tmp)) = 0.1;
tmp(isnan(tmp)) = -0.1;

subplot(2,3,6); hold off

% Prima mappa sfondo con tutti i valori colorati (riempitivo)
contourf(lat_vel, -depth_ref, tmp, -0.1:0.01:0.1, 'LineStyle', 'none'); hold on

% Colormap dettagliata
colormap(gca, (clr_V2));

% Sovrapponi mappa reale
[CS,CH] = contourf(lat_vel, -depth_ref, mvel, -0.15:0.01:0.15, 'LineStyle', 'none');
contour(lat_vel, -depth_ref, mvel, [0 0], 'color', [0.7 0.7 0.7]);

clim([-0.15 0.15])

% Colorbar dettagliata
ax = m_contfbar(1.2, [.05 0.9], CS, CH, 'axfrac', .035, 'endpiece', 'no');
ticks = (-0.15:0.05:0.15);
set(ax, 'tickdir', 'in', 'ylim', [-0.15 0.15], ...
    'fontsize', 10, 'ycolor', xyclr, 'xcolor', xyclr, ...
    'ytick', ticks, ...
    'yticklabel', arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false))
set(ax, 'colormap', clr_V2)
title(ax, 'm/s', 'FontSize', 10, 'color', xyclr, 'FontWeight', 'normal')

set(gca, 'Position', [0.7 0.1 0.23 0.35])

% Asse X con latitudine
xticks(-72:2:-46)
vals = -72:2:-46;
labels = cell(1,length(vals));
for k = 1:length(vals)
    if mod(abs(vals(k)),4) == 0
        labels{k} = sprintf('%d\\circS', abs(vals(k)));
    else
        labels{k} = '';
    end
end
xticklabels(labels)
set(gca, 'XDir', 'reverse')

% Titolo
title('(f) Mean velocity', 'FontSize', 12, 'FontWeight', 'normal', ...
    'Units', 'normalized', 'Position', [0.5, 1.05, 0])
axes_handles = findall(gcf, 'type', 'axes');
set(axes_handles, 'TickDir', 'out');

%PIE CREATE X,Y AXES
timeYMD = datestr(time_avg','yyyymmdd');
data1 = str2num(timeYMD(:,5:6))'; %Month numbers
data1(end+1)=NaN;
data2 = str2num(timeYMD(:,1:4))'; %year numbers
data2(end+1)=NaN;

%GET DATA RANGE
axis1 = [min(data1):max(data1)]';
axis2 = [min(data2):max(data2)]';
axis2(end+1)=NaN;

%create month axis if wanted
leg1 = {'jan','feb','mar','apr','may','jun','jul','ago','sep','oct','nov','dec'};

%Add N number of spaces in the middle
N = 3;
C = npf; %number of profiles/transect
C(end+1)=0;
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
print -dpng -r400 Figure1_px36.png
