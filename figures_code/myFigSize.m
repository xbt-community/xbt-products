function myFig = myFigSize(Fig_number,width,height)
% 
% Resize figure on screen according to papersize.
%     
figure(Fig_number);
set(gcf, 'PaperUnits', 'inches'); 
set(gcf, 'PaperSize', [width height]); 
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperPosition', [0 0 width height]);

% wysiwyg;
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(3:4) = ppos(3:4);
% pos(1:2) = [1 1];
set(gcf,'position',pos);
set(gcf,'units',unis);

myFig = figure(Fig_number);