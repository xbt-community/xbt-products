%Figure 7 - Comparison GOXBT with ARGO RG
%Marlos Goes Sep 29, 2025

clear
close all

argo = '../Argo/';
%argo = 'Argo/';  %Download the RG Argo folder first


    %% PLOTTING section
    close all

   direc = {'AX32','AX97','PX30','PX36','PX40'};
    letters={'a) ','b) ','c) ','d) ','e) '};
     vlim = [.8 .2 .7 .2 .8];
     labE = {'N','W','E','S','E'};
     stp = [0.05 0.01 0.05 0.05 .05];
     alim = [ 32.5500  40.4500;...
             -40.8167 -38;...
             153.1800  160.00;...
             -71.8300 -46.0900;...  
              139.00   150.00];

    fig3 = myFigSize(1,7,10); clf
    count = 0;
         for ii = 1:length(direc)
             file = [direc{ii},'/',direc{ii},'_gridded_new.nc'];
           if ~exist(file,'file')
              file = [direc{ii},'/',direc{ii},'_gridded.nc'];
           end

    vel_ref = ncread(file,'velocity');
    GStransp = ncread(file,'geostrophic_transport');
    time_avg = ncread(file,'time')+datetime(1950,1,1);
    lon_vel = ncread(file,'longitude');
    lat_vel = ncread(file,'latitude');
    depth_ref = ncread(file,'depth');
    mid_lat = lat_vel(2:end,:);


        mvel = median(vel_ref,3,'omitmissing');
        mask = double(~isnan(mvel)); mask(mask==0)=nan;

        fileargo = [argo,'argo_',direc{ii},'_1x1.mat'];
        load(fileargo)
        if any(ii==[2 3 5])
            lat_vel = mean(lon_vel,2,'omitmissing');
            mid_lat = mid_lon;
        end
        nx = size(mid_lat,2);
        depth_ref = double(depth_ref);
        vel0 = interp2(mid_lat,Z(:,ones(1,nx)),vel0,lat_vel,depth_ref(:)');
        Z = depth_ref;mid_lat = lat_vel;
        vel0 = vel0.*mask;


        subplot(5,2,2*ii-1)
        contourf(lat_vel,-depth_ref,mvel,-0.8:stp(ii):0.8,'LineStyle','none');
        hold on
        clim([-vlim(ii) vlim(ii)]) % = clim;
        contour(lat_vel,-depth_ref,mvel,[0 0],'color',[.7 .7 .7]);
        set(gca,'color',[.7 .7 .7])
        set(gca,'fontsize',14)
        if ii==1
            set(gca,'xdir','reverse')
        end
        title([letters{ii},direc{ii}])
        xlim(alim(ii,:))
        ylabel('Depth (m)')

        tick_lat = get(gca,'xtick');
        label_lat = strcat(num2str(abs(tick_lat)'),labE{ii});
        set(gca,'fontsize',10)
        set(gca,'ytick',-800:200:0)
        set(gca,'xtick',tick_lat,'xticklabel',label_lat,'XTickLabelRotation',0)


        subplot(5,2,2*ii)
        [CS,CH] = contourf(mid_lat,-Z,vel0,-0.8:stp(ii):0.8,'LineStyle','none');
        shading flat
        clim([-vlim(ii) vlim(ii)]) 
        hold on
        contour(mid_lat,-Z,vel0,[0 0],'color',[.7 .7 .7]);
        set(gca,'color',[.7 .7 .7])
        set(gca,'fontsize',12)
        if ii==1
            set(gca,'xdir','reverse')
        end

        xlim(alim(ii,:))
        title('RG Argo')
        divergingmap(17);
        pos = get(gca,'Position');
        ax = colorbar('position',[.92 pos(2)+.01 .01 .1]);
        title(ax,'m/s','FontSize',10,'FontWeight','normal')
        set(gca,'fontsize',10)
        set(gca,'ytick',-800:200:0)
        set(gca,'xtick',tick_lat,'xticklabel',label_lat,'XTickLabelRotation',0)
        end
        exportgraphics(gcf, 'Figure7.png', 'Resolution', 400);

