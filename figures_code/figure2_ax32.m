addpath /phodnet/share/mgoes/matlab/GSW_2015/
%Output: GStransp, GSspd_z, GStransp_z; 
% Extra Output: GStransp_alt, GSlat_alt, GSspd_alt
close all
clear

load Output.mat

% load ssh for the reference section (matching the XBT grid)
%load /phodnet/share/dong/gulfstream/oleander/
load ssh4ax32_2025.mat ssh_ref
[nz_ref,nx_ref,nt] = size(temp_ref);

% compute the area of each grid cell
delt_z = diff(depth_ref(1:2));
Area_GS = nan(nz_ref,nx_ref-1,nt);
vel_alt = nan(length(lat_vel),nt);
for i = 1:nt
    % delt_y=deg2km(distance(lat_ref(1:nx_ref-1),lon_ref(1:nx_ref-1,i),lat_ref(2:nx_ref),lon_ref(2:nx_ref,i),'degrees'))*1000;
    delt_y = sw_dist(lat_ref(:,i),lon_ref,'km')*1000;
    Area_GS(:,:,i) = repmat(delt_y',[nz_ref 1])*delt_z;
    vel_alt(:,i)=-gsw_grav(lon_vel)./gsw_f(lon_vel).*(diff(ssh_ref(:,i))./delt_y);
end

% Determine the maximum velocity at GS core (GSspd),
% and GS location (GSlat, GSlon) and volume transport
GSspd_z = nan(nz_ref,nt);
GStransp_z = nan(nz_ref,nt);
GSspd = nan(nt,1);
GSlat = GSspd;
GSlon = GSspd;
GStransp = GSspd;
GSspd_alt = GSspd;
GSlat_alt = GSspd;
GStransp_alt=GStransp;
indx_GS = find(lat_vel>=36.5 & lat_vel<=38.5);  %GS latitude range

for i = 1:nt
    transp = sum(vel_ref(:,:,i).*Area_GS(:,:,i),1,'omitmissing');
    vel0 = mean(vel_ref(1:25,:,i),1); % use top 50-m averaged velocity to determine GS core
    [~,j] = max(vel0(indx_GS)); % based on near surface max velocity
    indx0 = indx_GS(j);
    GSspd(i) = vel_ref(1,indx0,i);
    GSlat(i) = lat_vel(indx0);
    GSlon(i) = lon_vel(indx0,i);
    % find GS northern and southern boundary (where current changes sign)
    k1 = find(vel0(1:indx0-1)<0,1,'last'); %southern boundary based on transport
    if isempty(k1)
        k1 = find(transp(1:indx0-1)<0,1,'last'); %southern boundary based on velocity
    end
    if isempty(k1), k1 = 0; end
    k2 = find(vel0(indx0:nx_ref-1)<0,1,'first');
    if isempty(k2)
        k2 = find(transp(indx0:nx_ref-1)<0,1,'first');
    end
    k2 = k2 + indx0 - 1;
    GStransp(i) = sum(transp(k1+1:k2-1));
    GStransp_z(:,i) = sum(vel_ref(:,k1+1:k2-1,i).*Area_GS(:,k1+1:k2-1,i)/delt_z,2,'omitmissing');
    GSspd_z(:,i) = max(vel_ref(:,k1+1:k2-1,i),[],2);
    GSlat_SB(i) = lat_ref(k1+1);
    GSlat_NB(i) = lat_ref(k2);
    % from altimeter SSH
    if ~isnan(mean(vel_alt(indx_GS,i)))
        [GSspd_alt(i),j] = max(vel_alt(indx_GS,i));
        GSlat_alt(i)=lat_vel(indx_GS(j));
        GStransp_alt(i) = ssh_ref(k1+1,i)-ssh_ref(k2-1,i);
    end
end

%save Output.mat -append GStransp_alt GSlat_alt GSspd_alt GSlat GSlon GStransp

%Calculate Correlation between altimetry and XBT
inan = ~isnan(GStransp_alt+GStransp);
C = corr(GStransp_alt(inan),GStransp(inan))

%%
fig3 = myFigSize(3,7,3); clf % papersize (10,10)

plot(time_avg,GStransp/1e6,'color',[0 .45 .75],'LineWidth',1);
hold on
plot(time_avg,GStransp/1e6,'.','color',[0 .45 .75],'markersize',6);
ylabel('Transport (Sv)')
%set(gca,'position',[0.53 0.08 0.4 0.23])
set(gca,'position',[0.1 0.1 0.8 0.8])
h1=gca;h2=axes('position',get(h1,'position'));
plot(time_avg,GStransp_alt,'color','r','LineWidth',1);hold on
ylabel('cross-front \DeltaSSH (m)')
set(h2,'yaxislocation','right','xaxislocation','top','color','none')
set(h2,'xcolor','k','ycolor','r','box','off','xticklabel',[])
set(h2,'ylim',[0.8 1.7],'xlim',[time_avg(1)-30 time_avg(end)+30])
set(h2,'xtick',datenum(2010:2:2024,1,1),'xticklabel',[])
set(h1,'yaxislocation','left','xaxislocation','bottom','color','none')
set(h1,'xcolor','k','ycolor',[0 0.45 0.75],'box','off')
set(h1,'xtick',datenum(2010:2:2024,1,1),'xticklabel',2010:2:2024)
set(h1,'ylim',[26 44],'xlim',[time_avg(1)-30 time_avg(end)+30])
text(time_avg(70),1,'XBT','Color',[0 .45 .75])
text(time_avg(70),0.88,'Altimeter','Color','r')
text(time_avg(1)-30,1.774,'(h) GS transport')
text(time_avg(1)+90,1.65,sprintf('R = %4.2f',C))

print -dpng -r600 figure2.AX32_GS.png


