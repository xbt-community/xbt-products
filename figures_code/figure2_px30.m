addpath(genpath('/datasets/work/soop-xbt/work/cowley/matlab_localtools/toolbox/local/gsw_toolbox/'))
addpath /datasets/work/soop-xbt/work/UOT/programs/MATLAB/seawater_ver3_3/
%Output: GStransp, GSspd_z, GStransp_z; 
% Extra Output: GStransp_alt, GSlat_alt, GSspd_alt
close all
clear
%%
load Output.mat

% load ssh for the reference section (matching the XBT grid)
%load /phodnet/share/dong/gulfstream/oleander/
load ssh4px30.mat ssh_ref_px30
[nz_ref,ny_ref,nt] = size(temp_ref);

% compute the area of each grid cell
delt_z = diff(depth_ref(1:2));
Area_GS = nan(nz_ref,ny_ref-1,nt);
vel_alt = nan(length(lon_vel),nt);
for i = 1:nt
    % delt_y=deg2km(distance(lat_ref(1:nx_ref-1),lon_ref(1:nx_ref-1,i),lat_ref(2:nx_ref),lon_ref(2:nx_ref,i),'degrees'))*1000;
    delt_y = sw_dist(lat_ref(:,i),lon_ref,'km')*1000;
    Area_GS(:,:,i) = repmat(delt_y',[nz_ref 1])*delt_z;
    vel_alt(:,i)=gsw_grav(lat_vel(:,i))./gsw_f(lat_vel(:,i)).*(diff(ssh_ref_px30(:,i))./delt_y);
end
%%
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
indx_GS = find(lon_vel>=153.1 & lon_vel<=160);  %GS longitude range to get the edge of EAC

for i = 1:nt
    transp = sum(vel_ref(:,:,i).*Area_GS(:,:,i),1,'omitmissing');
    vel0 = mean(vel_ref(1:25,:,i),1,'omitmissing'); % use top 50-m averaged velocity to determine GS core
    % if all(isnan(vel0))
    %     continue
    % end
    [~,j] = min(vel0(indx_GS)); % based on near surface max velocity
    indx0 = indx_GS(j);
    GSspd(i) = vel_ref(1,indx0,i);
    GSlat(i) = lat_vel(indx0,i);
    GSlon(i) = lon_vel(indx0);
    % find GS eastern and western boundary (where current changes sign)
    k1 = find(vel0(1:indx0-1)>0,1,'last'); %eastern boundary based on velocity
    if isempty(k1)
        k1 = find(transp(1:indx0-1)>0,1,'last'); %eastern boundary based on transport
    end
    if isempty(k1), k1 = 1; end
    k2 = find(vel0(indx0:indx_GS(end)-1)>0,1,'first');
    if isempty(k2)
        k2 = find(transp(indx0:indx_GS(end)-1)>0,1,'first');
    end
    k2 = k2 + indx0 - 1;
    % k1 = 0; k2 = indx_GS(end);
    if k2>k1+1
        GStransp(i) = sum(transp(k1+1:k2-1),"all","omitmissing");
        GStransp_z(:,i) = sum(vel_ref(:,k1+1:k2-1,i).*Area_GS(:,k1+1:k2-1,i)/delt_z,2,'omitmissing');
        GSspd_z(:,i) = min(vel_ref(:,k1+1:k2-1,i),[],2);
        GSlon_EB(i) = lon_ref(k1+1);
        GSlon_WB(i) = lon_ref(k2);

            % from altimeter SSH
        if ~isnan(mean(vel_alt(indx_GS,i),'omitmissing'))
            [GSspd_alt(i),j] = min(vel_alt(indx_GS,i));
            GSlat_alt(i)=lat_vel(indx_GS(j));
            if isnan(ssh_ref_px30(k1+1,i))
                inan = find(~isnan(ssh_ref_px30(:,i)));
                % if isempty(inan)
                %     continue
                % end
                k1 = inan(1);
            end
            GStransp_alt(i) = ssh_ref_px30(k1+1,i)-ssh_ref_px30(k2-1,i);
        end
    end
end
%%
save Output.mat -append GStransp_alt GSlat_alt GSspd_alt GSlat GSlon GStransp

%Calculate Correlation between altimetry and XBT
inan = ~isnan(GStransp_alt+GStransp);
C = corr(GStransp_alt(inan),GStransp(inan));

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
set(h2,'ylim',[-0.9 0.2],'xlim',[time_avg(1)-30 time_avg(end)+30])
set(h2,'xtick',datenum(1992:3:2024,1,1),'xticklabel',[])
set(h1,'yaxislocation','left','xaxislocation','bottom','color','none')
set(h1,'xcolor','k','ycolor',[0 0.45 0.75],'box','off')
set(h1,'xtick',datenum(1992:3:2024,1,1),'xticklabel',1992:3:2024)
set(h1,'ylim',[-35 2],'xlim',[time_avg(1)-30 time_avg(end)+30])
text(time_avg(50),-0.05,'XBT','Color',[0 .45 .75])
text(time_avg(60),-0.65,'Altimeter','Color','r')
text(time_avg(1)+30,0.15,'(h) GS transport')
text(time_avg(1)+90,0.05,sprintf('R = %4.2f',C))

print -dpng -r600 figure2.PX30_GS.png


