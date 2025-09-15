[ny_ref,nz_ref,nt] = size(temp_ref);

% compute the area of each grid cell
 delt_z = diff(depth_ref(1:2));
 Area_GS = nan(nz_ref,ny_ref-1,nt);
 vel_alt = nan(length(lon_vel),nt);
 
 for i = 1:nt
     delt_y=deg2km(distance(lat_ref(1:ny_ref-1,i),lon_ref(1:ny_ref-1,i),lat_ref(2:ny_ref,i),lon_ref(2:ny_ref,i),'degrees'))*1000;
     Area_GS(:,:,i) = repmat(delt_y,[1 nz_ref])'*delt_z;
     vel_alt(:,i)=gsw_grav(lat_vel(:,i))./gsw_f(lat_vel(:,i)).*(diff(ssh_ref_px40(:,i))./delt_y);
 end


%%
% Determine the maximum velocity at GS core (GSspd),
% and GS location (GSlat, GSlon) and volume transport
Blim = [139.75 155]%144]%150]%165];
GSspd_z = nan(nz_ref,nt);
GStransp_z = nan(nz_ref,nt);
GSspd = nan(nt,1);
GSlat = GSspd;
GSlon = GSspd;
GStransp = GSspd;
GSspd_alt = GSspd;
GSlon_alt = GSspd;
GStransp_alt=GStransp;
indx_GS = find(lon_vel(:,1)>=Blim(1) & lon_vel(:,1)<=Blim(2));  %GS latitude range
%indx_max = find(lon_vel(:,1)<=160,1,'last');
%%
for i = 1:nt
    mask = double(vel_ref(:,:,i)>0);
    transp = sum(vel_ref(:,:,i).*Area_GS(:,:,i).*mask,1,'omitmissing');
    vel0 = mean(vel_ref(1:25,:,i),1); % use top 50-m averaged velocity to determine GS core
    if all(isnan(vel0))
        continue
    end
    [~,j] = max(vel0(indx_GS)); % based on near surface max velocity
    indx0 = indx_GS(j);
 %   indx0 = min(indx0,indx_max);

    GSspd(i) = vel_ref(1,indx0,i);
    GSlat(i) = lat_vel(indx0,i);
    GSlon(i) = lon_vel(indx0,i);
    % find GS northern and southern boundary (where current changes sign)
    k1 = find(vel0(1:indx0-1)<0,1,'last'); %eastern boundary based on velocity
    if isempty(k1)
        k1 = find(transp(1:indx0-1)<0,1,'last'); %eastern boundary based on transport
    end
    if isempty(k1), k1 = 1; end
    k2 = find(vel0(indx0:ny_ref-1)<0,1,'first');   %First location of vel< 0
    if isempty(k2)
        k2 = find(transp(indx0:ny_ref-1)<0,1,'first');
    end
  %  k2 = min(k2 + indx0 - 1,indx_max);
    k2 = k2 + indx0 - 1;  %original
    [k1 k2]
    if k2>k1+1
        GStransp(i) = sum(transp(k1+1:k2-1));
        GStransp_z(:,i) = sum(vel_ref(:,k1+1:k2-1,i).*Area_GS(:,k1+1:k2-1,i)/delt_z,2,'omitmissing');
        GSspd_z(:,i) = max(vel_ref(:,k1+1:k2-1,i),[],2);
        GSlon_EB(i) = lon_ref(k1+1);
        GSlon_WB(i) = lon_ref(k2);
        % from altimeter SSH
        if ~isnan(mean(vel_alt(indx_GS,i),'omitmissing'))

            [GSspd_alt(i)] = max(vel_alt(indx_GS,i));
            [GSspd_alt(i),j] = max(vel_alt(indx_GS,i));
         %   GSlon_alt(i)=lon_vel(indx_GS(j));
            inan = find(~isnan(ssh_ref_px40(k1+1:k1+10,i)),1,'first'); %if there are nans in SSH, choose next
            GStransp_alt(i) = ssh_ref_px40(k2-1,i)-ssh_ref_px40(k1+inan,i);
        end
    end
end