% Load Output matrix

load('Output.mat')
clearvars -except time_avg
Time_cruise=datevec(time_avg)

% Load SSH matrix
load('ssh4px36_all_cruises.mat')

%% Remove XBT casts with time gaps caused by upstream latitude filtering applied to the dataset
lat_xbt{2}(97:98,:)=[];
lon_xbt{2}(97:98,:)=[];
time_xbt{2}(97:98,:)=[];
ssh_xbt{2}(97:98,:)=[];

lat_xbt{4}(111:116,:)=[];
lon_xbt{4}(111:116,:)=[];
time_xbt{4}(111:116,:)=[];
ssh_xbt{4}(111:116,:)=[];

lat_xbt{13}([34,37,40,42,44,45,47,49,50],:)=[];
lon_xbt{13}([34,37,40,42,44,45,47,49,50],:)=[];
time_xbt{13}([34,37,40,42,44,45,47,49,50],:)=[];
ssh_xbt{13}([34,37,40,42,44,45,47,49,50],:)=[];


lat_xbt{15}([80,82,84,86,87,88,90,91,94,96,97,99,101,104,105,109,110],:)=[];
lon_xbt{15}([80,82,84,86,87,88,90,91,94,96,97,99,101,104,105,109,110],:)=[];
time_xbt{15}([80,82,84,86,87,88,90,91,94,96,97,99,101,104,105,109,110],:)=[];
ssh_xbt{15}([80,82,84,86,87,88,90,91,94,96,97,99,101,104,105,109,110],:)=[];

save("ssh4px36_all_cruises_Modified.mat", "ssh_xbt","lat_xbt","lon_xbt","time_xbt")


%% Interpolate SSH onto the same latitude grid as the XBT data
load('ssh4px36_all_cruises_Modified.mat')
% Define latitude grid for the reference section (~20 km steps ~0.18°)
lat_ref = (-46:-0.18:-72)';  % latitude grid from -46° to -72° (southward)
ssh_interp = NaN(length(lat_ref), size(ssh_xbt,2));  % Preallocate array

for i = 1:size(ssh_xbt,2)% Small value to perturb duplicate entries

    % Modifiable copy of latitude values
    new_lat = lat_xbt{i};

    % Find unique values and their indices
    [~, ~, ic] = unique(new_lat, 'stable');
    counts = histcounts(ic, 1:(max(ic)+1));

    % Handle duplicates
    for valIdx = find(counts > 1)
        dupPositions = find(ic == valIdx);
        for n = 2:length(dupPositions)
            k = dupPositions(n);
            new_lat(k) = new_lat(k) + epsilon * (n-1);  % Aggiunge piccola perturbazione
        end
    end

    % Interpolation using the corrected latitude values
    ssh_interp(:,i) = interp1(new_lat, ssh_xbt{i}, lat_ref, 'linear', NaN);
end
clear ssh_xbt
clear ssh_ref
ssh_ref=ssh_interp;
save("SSH_interp_all_curises_Mod.mat", 'ssh_ref')

