%%  Calculate Transport
%Input: depth_ref, temp_ref, lon_ref, lat_ref, vel_ref, lat_vel, lon_vel
%Output: GStransp, GSspd_z, GStransp_z; 
cd('C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2\')

disp('Calc transport')

load('Output.mat')
% ssh_ref=ssh_xbt;
load('C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2\SSH_inerp_all_curises_Mod.mat')


%Define northern and southern boundary for transport calculations in
%degrees for whole ACC
BLim = [-64.50  -52.25];

[nz_ref,nx_ref,nt] = size(temp_ref);

% compute the area of each grid cell
 delt_z = diff(depth_ref(1:2));
 Area_NSAF = nan(nz_ref,nx_ref-1,nt);
 vel_alt = nan(length(lat_vel),nt);

 for i = 1:nt
     delt_y=deg2km(distance(lat_ref(1:nx_ref-1),lon_ref(1:nx_ref-1,i),lat_ref(2:nx_ref),lon_ref(2:nx_ref,i),'degrees'))*1000;
     Area_NSAF(:,:,i) = repmat(delt_y',[nz_ref 1])*delt_z;
         if i==nt
    continue
    end
    vel_alt(:,i)=-gsw_grav(lat_vel)./gsw_f(lat_vel).*(diff(ssh_ref(:,i))./delt_y);

 end

% Determine the maximum velocity at NSAF core (NSAFspd),
% and NSAF location (NSAFlat, NSAFlon) and volume transport
NSAFspd_z = nan(nz_ref,nt);
NSAFtransp_z = nan(nz_ref,nt);
NSAFspd = nan(nt,1);
NSAFlat = NSAFspd;
NSAFlon = NSAFspd;
NSAFtransp = NSAFspd;
NSAFtransp_alt=NSAFtransp;
indx_NSAF = find(lat_vel>=BLim(1) & lat_vel<=BLim(2));  %NSAF latitude range
% vel_ref(vel_ref>2)=NaN;
% vel_ref(vel_ref<-2)=NaN;
% 
for i = 1:nt
    if i==25 || i==16 || i==18 || i==28|| i==31 || i==36 || i==20 || i==26  %% Discarded cruises due to low number of XBT cast related
%  to the applyed Geospatial filter
        continue
    end

    transp = sum(vel_ref(:,:,i).*Area_NSAF(:,:,i),1,'omitnan');
    vel0 = mean(vel_ref(1:25,:,i),1,'omitnan'); % use top 50-m averaged velocity to determine NSAF core
    [~,j] = max(vel0(indx_NSAF)); % based on near surface max velocity
    indx0 = indx_NSAF(j);
    NSAFspd(i) = vel_ref(1,indx0,i);
    NSAFlat(i) = lat_vel(indx0);
    NSAFlon(i) = lon_vel(indx0,i);
    % find NSAF northern boundary where we find first positive value at northern latitude 
    % respect to the max vel index and  at southern latitudes respect to the Blim(1)
      k1 = find(transp(1:indx0-1)>0,1,'first'); 
    if isempty(k1)
        k1 = find(vel0(1:indx0-1)>0,1,'first'); 
    end
     if k1<indx_NSAF(1)
        k1=indx_NSAF(1);
     end  
     
    % find southern boundary where we find first negative value at southern
    % latitudes respecto the max vel index and at northern latitudes respect to the Blim(2)

   
    k2 = find(transp(indx0:nx_ref-1)<0,1,'first');
    if isempty(k2)
        k2 = find(vel0(indx0:nx_ref-1)<0,1,'first');
    end
    k2 = k2 + indx0 - 1;
     if k2 >indx_NSAF(end)
        k2=indx_NSAF(end);
    end
    NSAFtransp(i) = sum(transp(k1+1:k2-1),'omitnan');
    NSAFtransp_z(:,i) = sum(vel_ref(:,k1+1:k2-1,i).*Area_NSAF(:,k1+1:k2-1,i)/delt_z,2,'omitnan');
    NSAFspd_z(:,i) = max(vel_ref(:,k1+1:k2-1,i),[],2);
    if ~isnan(mean(vel_alt(indx_NSAF,i)))
        [NSAFspd_alt(i),j] = max(vel_alt(indx_NSAF,i));
        NSAFlat_alt(i)=lat_vel(indx_NSAF(j));
        NSAFtransp_alt(i) = ssh_ref(k1+1,i)-ssh_ref(k2-1,i);
    end
clear vel0 transp
end

%save Output.mat -append NSAFlat NSAFlon NSAFtransp 
disp('Done!')
%clear

%Calculate Correlation between altimetry and XBT
inan = ~isnan(NSAFtransp_alt+NSAFtransp);
C = corr(NSAFtransp_alt(inan),NSAFtransp(inan))

%%
fig3 = myFigSize(3,7,3); clf % papersize (10,10)
% Plot only row indexes different from nan
valid_idx = ~isnan(NSAFtransp)
%
plot(time_avg(valid_idx),NSAFtransp(valid_idx)/1e6,'color',[0 .45 .75],'LineWidth',1);
hold on
plot(time_avg(valid_idx),NSAFtransp(valid_idx)/1e6,'.','color',[0 .45 .75],'markersize',6);
ylim([-2 10])
ylabel('Transport (Sv)')
%set(gca,'position',[0.53 0.08 0.4 0.23])
set(gca,'position',[0.1 0.1 0.8 0.8])
h1=gca;h2=axes('position',get(h1,'position'));
valid_idx = ~isnan(NSAFtransp_alt)
plot(time_avg(valid_idx),NSAFtransp_alt(valid_idx),'color','r','LineWidth',1);
hold on
plot(time_avg(valid_idx),NSAFtransp_alt(valid_idx),'.','color','r','markersize',6);
ylabel('cross-front \DeltaSSH (m)')
set(h2,'yaxislocation','right','xaxislocation','top','color','none')
set(h2,'xcolor','k','ycolor','r','box','off','xticklabel',[])
set(h2,'ylim',[0.5 1.8],'xlim',[time_avg(1)-30 time_avg(end)+30])
set(h2,'xtick',datenum(1994:4:2024,1,1),'xticklabel',[])
set(h1,'yaxislocation','left','xaxislocation','bottom','color','none')
set(h1,'xcolor','k','ycolor',[0 0.45 0.75],'box','off')
set(h1,'xtick',datenum(1994:4:2024,1,1),'xticklabel',1994:4:2024)
set(h1,'ylim',[0 100],'xlim',[time_avg(1)-30 time_avg(end)+30])
text(time_avg(end-10),-0.07,'XBT','Color',[0 .45 .75])
text(time_avg(end-9)-140,0.27,'Altimeter','Color','r')
text(time_avg(1)-30,1.9,'(h)ACC transport')
text(time_avg(1)+90,1.84,sprintf('R = %4.2f',C))


print -dpng -r600 Figure2_px36_ACC.png
