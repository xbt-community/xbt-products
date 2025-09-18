%% Packages needed: GSW, M_MAP, cmocean, TS_PACKAGE
% Add paths for required toolboxes and libraries
cd('C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2\')

addpath GSW_2015/
addpath m_map/
addpath TS_PACKAGE_THACKER_noyear_globe/

%% Load Data
% Load dataset
load('C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2\temp_px36_paper.mat');
% Inputs:
%   depth_xbt, temp_xbt, lat_xbt, lon_xbt
% Outputs (initial):
%   time_xbt, temp_ref
%   lon_ref, lat_ref, depth_ref, npf, lon_xbt, lat_xbt

%% == Reference Section Retrieval
disp('calc reference section')

% Define depth reference grid (<= 760 m)
kz = find(depth_xbt <= 760);
depth_ref = depth_xbt(kz);

% Define latitude grid for the reference section (~20 km steps ~0.18°)
lat_ref = (-46:-0.18:-72)';  % latitude grid from -46° to -72° (southward)
nx_ref = length(lat_ref);
nz_ref = length(depth_ref);
nt = length(npf);
nz = length(depth_xbt);

% Initialize longitude reference array
lon_ref = nan(nx_ref, nt);

% Loop through each transect to build reference longitude
for i = 1:nt
    lat_i = lat_xbt{i};
    lon_i = lon_xbt{i};
    temp_i = temp_xbt{i};  % associated temperature

    % Remove duplicate points and handle epsilon injection
    new_lat = [];
    new_lon = [];
    new_temp = [];
    j = 1;
    while j <= length(lat_i)
        if j < length(lat_i)
            lat_equal = (lat_i(j) == lat_i(j+1));
            lon_equal = (lon_i(j) == lon_i(j+1));
            if lat_equal && lon_equal
                % Duplicate point pair: average temperature, keep one coord
                avg_lat = lat_i(j);
                avg_lon = lon_i(j);
                avg_temp = mean([temp_i(j), temp_i(j+1)], 'omitnan');
                % Keep only points within the study bounds
                if avg_lat >= -72 && avg_lat <= -46 && avg_lon >= 170 && avg_lon <= 180
                    new_lat(end+1) = avg_lat;
                    new_lon(end+1) = avg_lon;
                    new_temp(end+1) = avg_temp;
                end
                j = j + 2;
            elseif lat_equal && ~lon_equal
                % Same latitude, inject tiny epsilon to second lat point
                new_lat(end+1) = lat_i(j);
                new_lon(end+1) = lon_i(j);
                new_temp(end+1) = temp_i(j);
                new_lat(end+1) = lat_i(j+1) + 1e-11;
                new_lon(end+1) = lon_i(j+1);
                new_temp(end+1) = temp_i(j+1);
                j = j + 2;
            elseif ~lat_equal && lon_equal
                % Same longitude, inject tiny epsilon to second lon point
                new_lat(end+1) = lat_i(j);
                new_lon(end+1) = lon_i(j);
                new_temp(end+1) = temp_i(j);
                new_lat(end+1) = lat_i(j+1);
                new_lon(end+1) = lon_i(j+1) + 1e-11;
                new_temp(end+1) = temp_i(j+1);
                j = j + 2;
            else
                % No duplicates, keep current point if in bounds
                if lat_i(j) >= -72 && lat_i(j) <= -46 && lon_i(j) >= 170 && lon_i(j) <= 180
                    new_lat(end+1) = lat_i(j);
                    new_lon(end+1) = lon_i(j);
                    new_temp(end+1) = temp_i(j);
                end
                j = j + 1;
            end
        else
            % Last single point
            if lat_i(j) >= -72 && lat_i(j) <= -46 && lon_i(j) >= 170 && lon_i(j) <= 180
                new_lat(end+1) = lat_i(j);
                new_lon(end+1) = lon_i(j);
                new_temp(end+1) = temp_i(j);
            end
            j = j + 1;
        end
    end

    % Consolidate duplicate latitudes by averaging longitudes
    [unique_lat, ~, ic] = unique(new_lat);
    unique_lon = zeros(size(unique_lat));
    for k = 1:length(unique_lat)
        idxs = find(ic == k);
        unique_lon(k) = mean(new_lon(idxs),'omitnan');
    end

    % Interpolate lon_ref onto lat_ref with linear interpolation
    lon_ref(:, i) = interp1(unique_lat, unique_lon, lat_ref, 'linear');
    % Fill missing values at boundaries
    lon_ref(:, i) = fillmissing(lon_ref(:, i), 'linear', 'EndValues', 'nearest');
end

% Compute longitude deviation and filter outliers
dlon = (lon_ref - mean(lon_ref, 2, 'omitnan'));
%kgood = find(mean(abs(dlon), 1, 'omitnan') <= 0.25);
mean_lon_ref = mean(lon_ref, 2, 'omitnan');

%% == Temperature Reference Grid Interpolation
% Initialize temperature reference array: [depth x lat x time]
temp_ref = nan(nz_ref, nx_ref, nt);
for i = 1:nt
    tempi = temp_xbt{i};
    lati = lat_xbt{i};
    loni = lon_xbt{i};

    % Identify close-by profiles (<10 km apart)
    dist_y = deg2km(distance(lati(1:end-1), loni(1:end-1), lati(2:end), loni(2:end), 'degrees'));
    close_idx = find(dist_y < 10);

    if isempty(close_idx)
        % No clustering needed
        new_lat = lati;
        new_lon = loni;
        new_temp = tempi;
    else
        % Cluster adjacent profiles
        diff_idx = diff(close_idx);
        cluster_breaks = [0, find(diff_idx > 1)', length(close_idx)];
        new_lat = [];
        new_lon = [];
        new_temp = [];
        start_idx = 1;
        for c = 1:length(cluster_breaks)-1
            cluster_range = close_idx(start_idx):(close_idx(cluster_breaks(c+1))+1);
            start_idx = cluster_breaks(c+1) + 1;
            avg_lat = mean(lati(cluster_range));
            avg_lon = mean(loni(cluster_range));
            avg_temp = mean(tempi(:, cluster_range), 2, 'omitnan');
            new_lat(end+1) = avg_lat;
            new_lon(end+1) = avg_lon;
            new_temp(:, end+1) = avg_temp;
        end
        % Append single points not in clusters
        all_points = 1:length(lati);
        clustered_pts = unique([cluster_range]);
        single_pts = setdiff(all_points, clustered_pts);
        for sp = single_pts
            new_lat(end+1) = lati(sp);
            new_lon(end+1) = loni(sp);
            new_temp(:, end+1) = tempi(:, sp);
        end
        % Sort by latitude
        [new_lat, sort_idx] = sort(new_lat);
        new_lon = new_lon(sort_idx);
        new_temp = new_temp(:, sort_idx);
    end


    % Prevent any exact latitude duplicates
    epsilon = 1e-6;  %
    % ic(k) = index of the unique value corresponding to new_lat(k)
    [~, ~, ic] = unique(new_lat, 'stable');
    % Count how many times each unique value occurs
    counts = histcounts(ic, 1:(max(ic)+1));
    % For each unique value that occurs more than once...
    for valIdx = find(counts > 1)
        % find all positions of that value
        dupPositions = find(ic == valIdx);
        % keep the first unchanged; add epsilon*(n-1) to the others
        for n = 2:length(dupPositions)
            k = dupPositions(n);
            new_lat(k) = new_lat(k) + epsilon * (n-1);
        end
    end


    % Interpolate temperature onto lat_ref for each depth
    for iz = 1:nz_ref
        valid_idx = find(~isnan(new_temp(iz, :)));
        if isempty(valid_idx)
            temp_ref(iz, :, i) = NaN;
        else
            temp_ref(iz, :, i) = interp1(new_lat(valid_idx), new_temp(iz, valid_idx), lat_ref, 'linear', NaN);
        end
    end
end

% Save intermediate results\...Output.mat time_xbt temp_ref lon_ref lat_ref depth_ref npf lon_xbt lat_xbt
save Output.mat time_xbt temp_ref lon_ref lat_ref depth_ref npf lon_xbt lat_xbt

disp('Reference section smoothing complete')

%% ======================= Salinity Calculation =======================
disp('calc salinity')

load Output.mat
nt = length(npf);
nx_ref = length(lat_ref);
nz_ref = length(depth_ref);

% Compute average observation time for each transect
time_avg = nan(nt,1);
for it = 1:nt
    time_avg(it) = mean(time_xbt{it}, 'omitnan');
end
timeYMD = datestr(time_avg', 'yyyymmdd');
time2 = str2num(timeYMD);
time2 = time2(:, ones(1, nx_ref))';

% Calculate salinity all at once using GOES method\auxT = reshape(temp_ref, [nz_ref, nx_ref*nt]);
auxLon = reshape(lon_ref, [1, nx_ref*nt]);
auxLat = reshape(lat_ref(:, ones(1, nt)), [1, nx_ref*nt]);
auxTime = reshape(time2, [1, nx_ref*nt]);
auxT = reshape(temp_ref,[nz_ref nx_ref*nt]);
salt_ref = Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(auxT, depth_ref, auxLat, auxLon, auxTime, 0, depth_ref, 'goes');

% Reshape salinity back to original grid
salt_ref = reshape(salt_ref, [nz_ref, nx_ref, nt]);

save -append Output.mat salt_ref time_avg
disp('Salinity calculation complete')


%%  ======================= Compute velocity (GSW toolbox)  =======================
% first compute pressure from depth at each grid
%load /phodnet/share/dong/gulfstream/oleander/temp_ax32_gridded.mat
%Input: temp_ref,depth_ref,salt_ref,lat_ref,lon_ref
%Output:vel_ref lat_vel lon_vel depth_ref dh0_ref
disp('compute velocity. Takes some time to finish')
addpath C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2
load('Output.mat')

[nz_ref,nx_ref,nt] = size(temp_ref);
[~,ga_refInd] = min(abs(depth_ref-740));  %Reference for dynamic height output (760-800 m)
ref_dep = depth_ref(ga_refInd);
pres_ref = nan(nz_ref,nx_ref);

% use mean_lon_ref
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
    [SA, in_ocean] = gsw_SA_from_SP(salt_ref(:,:,it),pres_ref,lon_ref(:,it),lat_ref);
    CT = gsw_CT_from_t(SA,temp_ref(:,:,it),pres_ref);
    ga = gsw_geo_strf_dyn_height(SA,CT,pres_ref,0);
    inan = ~isnan(ga(1,:));
    DH_abs=ga*0;
    [DH_abs(:,inan),addep] = add_abs_gpan(ga(:,inan),lat_ref(inan),lon_ref(inan,it),depth_ref,ref_dep);
    ga = DH_abs;
    [vel0,mid_lat,mid_lon] = gsw_geostrophic_velocity(ga,lon_ref(:,it),lat_ref,pres_ref);
    lat_vel = mid_lat(1,:)';
    lon_vel(:,it) =  mid_lon(1,:)';
    vel_ref(:,:,it) = -vel0;
    dh0_ref(:,it) = ga(1,:)/9.81; %Convert from m2/s2 to meter at surface

end

S10=salt_ref(:,:,10);
S11=salt_ref(:,:,11);
S12=salt_ref(:,:,12);
V10=vel_ref(:,:,10);
V11=vel_ref(:,:,11);
V12=vel_ref(:,:,12);



save -append Output.mat vel_ref lat_vel lon_vel depth_ref dh0_ref
disp('Done!')
clear
