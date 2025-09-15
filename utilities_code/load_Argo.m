%Load Argo data and interpolates tp the lat,lon specified.
%ADDEP output is in cm/dyn
%by Marlos Goes  Aug 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [addep,z,temp,sal]= load_Argo(lon,lat,month)
direc = './';%[preffix,'IPRC/'];

%filename = [direc,'argo_CLIM_2005-2020_grd.nc'];
%filename = [direc,'argo_CLIM_1x1.nc'];
filename = ['ArgoData.nc'];    %ANNUAL CLIMATOLOGY

if nargin<3
    %For annual data.
    sal  = ncread(filename,'SALT');
    sal = permute(sal,[2 1 3]);
    temp  = ncread(filename,'TEMP');
    temp = permute(temp,[2 1 3]);
    addep = ncread(filename,'ADDEP');
    addep = permute(addep,[2 1 3]);
else
    %For monthly data.
    %filename = [direc,'argo_CLIM_2005-2016_grd.nc'];

    sal  = ncread(filename,'SALT',[1 1 1 month],[Inf Inf Inf 1]);
    temp = ncread(filename,'TEMP',[1 1 1 month],[Inf Inf Inf 1]);
    addep = ncread(filename,'ADDEP',[1 1 1 month],[Inf Inf Inf 1]);
    % addep = permute(addep,[2 1 3]);
    % sal = permute(sal,[2 1 3]);
    % temp = permute(temp,[2 1 3]);

end


MISS = min(addep(:));
if MISS<-100
    sal(sal == MISS) = nan;
    temp(temp == MISS) = nan;
    %ptemp(ptemp == MISS) = nan;
    addep(addep == MISS) = nan;
end

x = double(ncread(filename,'LONGITUDE'));
y = double(ncread(filename,'LATITUDE'));
z = ncread(filename,'LEVEL');

x(x>180) = x(x>180)-360; %FROM 0-360 to -180-180
[~,ilon] = sort(x);

if size(addep,1) == length(lon)
    addep = addep(ilon,:,:,:);  %SOMHOW WAS INVERTED
    temp  = temp(ilon,:,:,:);   %SOMHOW WAS INVERTED
    sal   = sal(ilon,:,:,:);    %SOMHOW WAS INVERTED
else
    addep = addep(:,ilon,:,:);  %SOMHOW WAS INVERTED
    temp  = temp(:,ilon,:,:);   %SOMHOW WAS INVERTED
    sal   = sal(:,ilon,:,:);
end
x = x(ilon);



%Use intepolation instead
[X,Y] = meshgrid(x,y);
data_ini = nan*lat;
st = size(addep,3);

for qt = 1:st
       sat = addep(:,:,qt);
     % disp(['interpolating cast #',num2str(i)])
      acha  = find(~isnan(sat));
      data_ini(:,qt) = griddata(X(acha),Y(acha),sat(acha),lon,lat,'cubic');

end
 addep = data_ini'*100/10;  %Transform to dyn cm and geopotential anomaly

%addep = interp1(x(vec),ad,lon)*100/10;

return

