%% Script for processing and managing AX97 cruise data
% In the AX97 section, we perform two preprocessing steps before this script,
% which corresponds to Step 3. In Step 1, we remove poor-quality profiles based 
% on the RMSE relative to the climatology. In Step 2, we apply the fall-rate equation
% proposed by Cheng (2014). 

% Here, in Step 3, we interpolate the data onto a reference transect defined by a straight 
% line between Rio de Janeiro and Trindade Island, which follows the standard path of the cruises. 
% A spatial filter is then applied (using the smoothing function) to reduce noise along the section.
% Salinity is subsequently estimated using the method proposed by Goes, and the geostrophic 
% velocity is calculated based on the resulting temperature and salinity fields.

%% Start the script

% Close all open figure windows, clear all variables from the workspace, and clear the command window
close all, clear all, clc

% Add necessary toolboxes to the MATLAB path
addpath(genpath('/data2/tayannepf/toolbox_movar/'));
addpath(genpath('/home2/tayannepf/matlab/seawater/'));

% Set directory containing the netCDF files of the processed cruises
dire = '/data2/tayannepf/sumario_cruzeiros/02_proc_CH14/files_netcdf/';

% Load bathymetric mask (used later for masking land or invalid regions)
load /data2/tayannepf/sumario_cruzeiro/03_smooth2d/mascara.mat

%% Load netCDF data from AX97 line until 2023
names = dir([dire, '*.nc']);
le = length(names);

% Loop through each netCDF file and read basic variables
for ii = 1:le
    dados_FRE(ii).long = ncread([dire, names(ii).name],'longitude');
    dados_FRE(ii).prof = ncread([dire, names(ii).name],'depth');
    dados_FRE(ii).temp = ncread([dire, names(ii).name],'temperature');
end

%% Create reference grid for interpolation

% Create a fixed longitude grid (radial definition)
longitude = [-40.9:0.16667:-40.1 -40:0.25:-29.6 -29.65]';

% Estimate the corresponding latitudes along the transect using linear regression
[p,~] = polyfit([-40.9 -29.65], [-23 -20.55], 1);
latitude = polyval(p, longitude); clear p

% Define vertical depth levels (0 to 800 m with 10 m spacing)
depth = 0:10:800;

% Create mesh grid for interpolation
[A,B] = meshgrid(longitude, depth);

% Initialize storage and timing
k = 1;
date(k,2) = datetime;
tic

% Process each cruise profile
for i = 1:length(dados_FRE)
    
    % Extract and convert data to double precision
    long = double(dados_FRE(i).long);
    prof = double(dados_FRE(i).prof);
    temp = double(dados_FRE(i).temp);

    % Replace zero values with NaN
    temp(temp == 0) = NaN;
    prof(prof == 0) = NaN;

    % Adjust longitude array dimensions if necessary
    if ~isequal(size(temp), size(long))
        disp('TEMP and LONG dimensions do not match, correcting...')
        long = repmat(long,1,size(temp,2));
    end

    % Adjust depth array dimensions if necessary
    if ~isequal(size(temp), size(prof))
        disp('TEMP and PROF dimensions do not match, correcting...')
        prof = repmat(prof,2,size(temp,1))';
    end

    % Find valid (non-NaN) temperature data
    j = find(~isnan(temp));
    y = find(~isnan(prof));

    % Fix edge case with NaNs at the first row/column
    if isnan(prof(1,1)) || isnan(prof(1,4))
        prof(:,1) = 0; 
        prof(1,:) = 0;
    end

    %% Interpolate temperature to the reference grid
   
    % Use nearest-neighbor interpolation
    F = scatteredInterpolant(long(j), prof(j), temp(j), 'nearest');

    % Interpolate onto the regular grid
    T1 = F(A,B);

    % Apply LOESS smoothing to interpolated temperature
    [T,~] = smooth2d_loess(T1, longitude, depth, 1, 50, longitude, depth);

    % Extract cruise date from filename
    time = repmat(str2double(strcat(names(i).name(20:end-3))), length(latitude), 1);

    % Parameters for salinity estimation
    Pout = [depth]; 
    Pad = 0; % Do not add climatology
    method = 'Goes';

    %% Estimate salinity using empirical method
    [S2, S3, TT, PP] = Calc_sal_Thacker_Goes_EmDr_Stom_svd_globe(T, depth, latitude, longitude, time(1), Pad, Pout, method);

    %% Compute dynamic height and geostrophic velocity
    ref_depth = 500;
    [dynH, vg] = calc_dynH_fixlev_clim(longitude', latitude', PP, TT, S2, ref_depth, time(1));

    % Store processed variables in arrays
    SAL2(:,:,k) = S2;
    SAL3(:,:,k) = S3;
    TEMP(:,:,k) = TT;
    DH(:,:,k) = dynH;
    
    % Compute SSH from dynamic height
    [ll, dd] = meshgrid(latitude, depth);
    SSH_movar(:,:,k) = dynH ./ -sw_g(ll, dd); clear ll dd

    VEL(:,:,k) = vg';
    DEPTH(:,k) = PP;
    LAT(:,k) = latitude;
    LON(:,k) = longitude;
    TIME(:,k) = time;

    k = k + 1;
end

toc

% Clean up temporary variables
clear i j m n k F dados p Pout PP prof temp long depth dynH T vg S2 S3 TT time A B

% Compute pressure field
PRESS = sw_pres(repmat(DEPTH(:,1),1,48), LAT(:,1)');

% Compute density from salinity, temperature, and pressure
for i = 1:length(date)
    DENS(:,:,i) = sw_dens(SAL3(:,:,i), TEMP(:,:,i), PRESS); 
end

%% Remove specific cruises with poor data quality
bad_cruises = [21 51 53 24 27];

% Remove data from all variables corresponding to bad cruises
VEL(:,:,bad_cruises) = [];
TIME(:,bad_cruises) = [];
LON(:,bad_cruises) = [];
LAT(:,bad_cruises) = [];
TEMP(:,:,bad_cruises) = [];
SAL2(:,:,bad_cruises) = [];
SAL3(:,:,bad_cruises) = [];
DH(:,:,bad_cruises) = [];
DEPTH(:,bad_cruises) = [];
SSH_movar(:,:,bad_cruises) = [];

% Save all processed results to a .mat file
save('/data2/tayannepf/sumario_cruzeiro/03_smooth2d/AX97_pos_proc_temp_sal_vel_smo_suavizado.mat', ...
    'DENS','VEL','date','LON','LAT','DEPTH','TIME','DH','SAL2','SAL3','SSH_movar','TEMP')


