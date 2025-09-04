clear all, close all, clc
addpath('/home2/tayannepf/matlab/gsw/')
addpath('/home2/tayannepf/matlab/gsw/library/')
addpath('/data2/Movar_analises/nc_final/Figures_Shenfu/Figures_code/')


load /data/MOVAR/sumario_cruzeiros/02_proc_CH14/files_netcdf/datas_cruz_CB_2004_2024_datetime.mat
clear dini dend

%%%% Organizando o SSH
load /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/ssh4ax97.mat ssh_ax97 day_ssh month_ssh year_ssh lon_ax97
data_ssh=datetime(year_ssh,month_ssh,day_ssh); clear year_ssh month_ssh day_ssh

ddd=ismember(data_ssh,diaCB);
ind=find(ddd==1);

% new_data=data_ssh(ddd);
for iii=1:length(diaCB)-1;
   
ssh_ref(:,iii)=nanmean(ssh_ax97(5:end,ind(iii)-1:ind(iii)+1),2);
% ssh_ref1(:,iii)=nanmean(ssh_ax97(5:end,ind(iii)-1:ind(iii)+1),2);
% ssh_ref3(:,iii)=nanmean(ssh_ax97(5:end,ind(iii)-3:ind(iii)+3),2);

end
%%%%

load /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/output_ax97.mat

[nz_ref,nx_ref,nt] = size(temp_ref);



%%  Calculate Transport
%Input: depth_ref, temp_ref, lon_ref, lat_ref, vel_ref, lat_vel, lon_vel
%Output: BCtransp, BCspd_z, BCtransp_z; 
disp('Calc transport')

% load('Output.mat')

%Define northern and southern boundary for transport calculations in degrees
BLim = [-41 -39];

% compute the area of each grid cell
 delt_z = diff(depth_ref(1:2));
 Area_BC = nan(nz_ref,nx_ref-1);
 vel_alt = nan(length(lon_ref(1:end-1)),nt);

     delt_y=deg2km(distance(lat_ref(1:nx_ref-1),lon_ref(1:nx_ref-1),lat_ref(2:nx_ref),lon_ref(2:nx_ref),'degrees'))*1000;
     Area_BC(:,:) = repmat(delt_y',[nz_ref 1])*delt_z;
     vel_alt(:,:)=-gsw_grav(lon_ref(1:end-1))./gsw_f(lon_ref(1:end-1)).*(diff(ssh_ref(:,:))./delt_y);

% Determine the maximum velocity at BC core (BCspd),
% and BC location (BClat, BClon) and volume transport
BCspd_z = nan(nz_ref,nt);
BCtransp_z = nan(nz_ref,nt);
BCspd = nan(nt,1);
BClat = BCspd;
BClon = BCspd;
BCtransp = BCspd;
BCtransp_alt=BCtransp;
indx_BC = find(lon_ref>=BLim(1) & lon_ref<=BLim(2));  %BC latitude range
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

% from altimeter SSH
    if ~isnan(mean(vel_alt(indx_BC,i)))
        [BCspd_alt(i),j] = max(vel_alt(indx_BC,i));
        BClat_alt(i)=lon_ref(indx_BC(j));
        BCtransp_alt(i) = ssh_ref(k1+1,i)-ssh_ref(k2-1,i);
    end

end

disp('Done!')

%Calculate Correlation between altimetry and XBT
inan = ~isnan(BCtransp_alt+BCtransp);
C = corr(BCtransp_alt(inan),BCtransp(inan))


fig3 = myFigSize(3,7,3); clf % papersize (10,10)

plot(diaCB(1:end-1),BCtransp/1e6,'color',[0 .45 .75],'LineWidth',1);
hold on
plot(diaCB(1:end-1),BCtransp/1e6,'.','color',[0 .45 .75],'markersize',6);
ylabel('Transport (Sv)')
%set(gca,'position',[0.53 0.08 0.4 0.23])
set(gca,'position',[0.1 0.1 0.8 0.8])
h1=gca;h2=axes('position',get(h1,'position'));
plot(diaCB(1:end-1),BCtransp_alt,'color','r','LineWidth',1);hold on
ylabel('cross-front \DeltaSSH (m)')
set(h2,'yaxislocation','right','xaxislocation','top','color','none')
set(h2,'xcolor','k','ycolor','r','box','off','xticklabel',[])
set(h2,'ylim',[-.25 .1],'xlim',[diaCB(1)-30 diaCB(end-1)+30])
% set(h2,'xtick',datenum(2004:2:2024,1,1),'xticklabel',[])
set(h1,'yaxislocation','left','xaxislocation','bottom','color','none')
set(h1,'xcolor','k','ycolor',[0 0.45 0.75],'box','off')
% set(h1,'xtick',datenum(2004:2:2024,1,1),'xticklabel',2004:2:2024)
set(h1,'ylim',[-15 5],'xlim',[diaCB(1)-30 diaCB(end-1)+30])
text(diaCB(70),0,'XBT','Color',[0 .45 .75])
text(diaCB(70),-.17,'Altimeter','Color','r')
text(diaCB(1),0.12,'(h) BC transport')
text(diaCB(1),0.08,sprintf('R = %4.2f',C))

% print -dpng -r600 /data2/Movar_analises/nc_final/Figures_Shenfu/ax97/figure2.AX97_BC.png




