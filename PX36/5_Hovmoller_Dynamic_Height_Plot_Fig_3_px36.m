
cd('C:\Users\fanto\OneDrive - Uniparthenope\Marlos_Goes_TS-lookup-master\V1_Scripts\V2')
load('Output.mat')
%% Hovmoller of dynamic height

disp('Plotting section')
%addpath /phodnet/share/mgoes/matlab/mapcolor/
%addpath ../Colormap/colorbar_yyw/

clr_topo = m_colmap('blue',256);
clr_V2 = (othercolor('RdBu11',256)); clr_V2(1,:) = [.7 .7 .7];
%clr_T = cmocean('thermal'); clr_T(1,:) = [.7 .7 .7];
%clr_S = cmocean('haline'); clr_S(1,:) = [.7 .7 .7];
%clr_V = cmocean('balance'); clr_V(1,:) = [.7 .7 .7];
xyclr = [.4 .4 .4];

fig3 = myFigSize(3,5,8); clf
vals = -72:2:-46;
label_lat = cell(1,length(vals));
for k = 1:length(vals)
    if mod(abs(vals(k)),4) == 0
        label_lat{k} = sprintf('%d\\circS', abs(vals(k)));
    else
        label_lat{k} = '';
    end
end
nt = length(time_avg);
time_label=str2num(datestr(time_avg,'yyyy'));
%THIS USES THE REAL AXIS contourf(lat_ref,time_avg,dh0_ref',[0.6:.1:1.8]);
%hold on
%set(gca,'ytick',datenum(2010:2:2024,1,1),'yticklabel',2010:2:2024,'XTickLabelRotation',0)

%THIS USES A FAKE AXIS
%PEr evitare che contourf non ci plotti l'ultima Colonna dobbiamo
%aggiungere l'ultima colonna di nan alla fine di dh0_ref
dupl_col=dh0_ref(:,end);
dh0_ref = [dh0_ref, dupl_col];  % duplica ultima riga
contourf(lat_ref,1:nt+1,dh0_ref',[-3:.05:0.8],'edgecolor','none'); hold on
set(gca,'ytick',1:2:nt,'yticklabel',time_label(1:2:nt))

%% Create a custom symmetric colorbar with different limits so that 0 = white

nColors = 256;
clr = flipud(othercolor('RdBu11', nColors));  % o cmocean('balance')

% Define Range
vmin = -1.2;
vmax = 0.8;

% Normalize values on a symmetric scale around zero
max_abs = max(abs([vmin, vmax]));

% Resample the colormap on the symmetric scale
x = linspace(-max_abs, max_abs, nColors);  % scala simmetrica
x_actual = linspace(vmin, vmax, nColors);  % scala reale

% Interpolate the colormap onto the actual scale
clr_interp = interp1(x, clr, x_actual, 'linear');


colormap(gca,(clr_interp));

xlim([-72 -46])
set(gca,'xtick',-72:2:-46,'xticklabel',label_lat,'XTickLabelRotation',0)
set(gca,'position',[0.13 0.2 0.7 0.7])
ax = colorbar('position',[.85 .2 .03 .7]);
set(ax,'tickdir','out','fontsize',10,'ycolor',xyclr,'xcolor',xyclr)
set(ax,'colormap',(clr_interp))
set(gca,'xcolor',xyclr,'ycolor',xyclr,'fontsize',10)
set(gca, 'XDir', 'reverse')
caxis([vmin vmax])
title(ax,'m','FontSize',10,'color',xyclr,'FontWeight','normal')
title('Surface Dynamic Height')
print -dpng -r400 Figure3_px36.png

