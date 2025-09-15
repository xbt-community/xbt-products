% Carichiamo Output matrice ricavata da MArlos_adapted_routin

load('Output.mat')
clearvars -except time_avg
Time_cruise=datevec(time_avg)

% Carichiamo matrice di SSH fornita da Shenfu
load('ssh4px36_all_cruises.mat')


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


%% Interpoliamo la ssh sulla stessa griglia di lat dei dati xbt
load('ssh4px36_all_cruises_Modified.mat')
% Define latitude grid for the reference section (~20 km steps ~0.18°)
lat_ref = (-46:-0.18:-72)';  % latitude grid from -46° to -72° (southward)
ssh_interp = NaN(length(lat_ref), size(ssh_xbt,2));  % Prealloca

for i = 1:size(ssh_xbt,2)
    epsilon = 1e-6;  % Piccolo valore per perturbare i duplicati
    
    % Copia modificabile della latitudine
    new_lat = lat_xbt{i};
    
    % Trova valori unici e i loro indici
    [~, ~, ic] = unique(new_lat, 'stable');
    counts = histcounts(ic, 1:(max(ic)+1));
    
    % Gestione duplicati
    for valIdx = find(counts > 1)
        dupPositions = find(ic == valIdx);
        for n = 2:length(dupPositions)
            k = dupPositions(n);
            new_lat(k) = new_lat(k) + epsilon * (n-1);  % Aggiunge piccola perturbazione
        end
    end

    % Interpolazione usando la latitudine corretta
    ssh_interp(:,i) = interp1(new_lat, ssh_xbt{i}, lat_ref, 'linear', NaN);
end
clear ssh_xbt
clear ssh_ref
ssh_ref=ssh_interp;
save("SSH_inerp_all_curises_Mod.mat", 'ssh_ref')

