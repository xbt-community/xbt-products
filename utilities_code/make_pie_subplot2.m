
function [g,col] = make_pie_subplot(varargin)
%vector1,vector2,names1,names2,mesa,anos,C,N)
%This code makes a pie chart using two axes
% and different colors dependent on C values
%Usage: make_pie(vector1,vector2,names1,names2,mesa,anos,C,N)
% INPUT: vector1 - Inner axis values - ex: 1:12 for months
%        vector2 - outer axis values - ex: 2010:2020 for years
%        names1  - [if needed] labels for axis 1 - ex: {'jan','feb',...}
%        names2  - [if needed] labels for axis 2 - ex: {'2010','2011',...}
%        mesa[1:N]  - vector values for axis 1 - ex: [1 2 12 4 ...]
%        anos[1:N]  - vector values for axis 2 - ex: [2009 2020 2018...]
%        C[1:N] - values to define colors - ex: [100 50 20 30...]
%        N - number of boxes to skip in the middle - ex: 3
% OUTPUT: g - handle for vertical axis
%         col - handle for colorbar
%By: Marlos Goes - May, 2025

InpVars = cell(1,9);
InpVars(1:nargin) = varargin;
[vector1,vector2,names1,names2,mesa,anos,C,N,eixo] = InpVars{:};

if isempty(C), C = ones(size(mesa));end
if isempty(N), N=3;end
if isempty(eixo), eixo=[0 0 1 1];end

    C0 = C;
%ADD N CENTRAL SPACES TO RADIUS VECTOR
vector2 = [vector2; vector2(end)+(vector2(end)-vector2(end-1))*[1:N]'];
valLIM1 = [min(mesa) max(mesa)];
valLIM2 = [min(anos) max(anos)];
valLIMC = [min(C) max(C)];

for ii=length(names2)+1:length(vector2)
    aux = vector2(ii);
    if aux < valLIM2(1)||aux > valLIM2(2)
        aux = '';
    end
    aux2 = num2str(aux);
    names2{ii} = aux2;%(max(1,end-1):end);


end

if isempty(names1)
    for ii=1:length(vector1)
        aux = vector1(ii);
        if aux<valLIM1(1)||aux>valLIM1(2)
            aux = '';
        end
              aux2 = num2str(aux);
        names1{ii} = aux2;%(end-1:end);
    end
end
%NORMALIZE VECTORS
ano = vector2-min(anos)+1;
mes = vector1-min(mesa)+1;
anos = anos-min(anos)+1;
mesa = mesa-min(mesa)+1;
%C = C-min(C);C=C/max(C);

stpa = 1/diff([min(ano) max(ano)])/2;
stpa = (eixo(3))/diff([min(ano) max(ano)])/2;

vec1 = ones(1,length(mes));
vec2 = zeros(1,length(mes));

label1 = names1; %{'jan','feb','mar','apr','may','jun','jul','ago','sep','oct','nov','dec'};
label2 = names2;

for nn=1:length(label1)
    null{nn} = '';
end
valid1 = logical(ones(1,length(mes))*0);

%DEFINE COLORSCALE FROM C
Ccolors = spring(24);%autumn(24);%length(C));
Ccolors = hsv(length(vec1));

% colors = interp1(linspace(0,1,length(Ccolors)),Ccolors,C);

colors = interp1(linspace(valLIMC(1),valLIMC(2),length(Ccolors)),Ccolors,C);
Cax1 = linspace(valLIMC(1),valLIMC(2),length(vec1)); 
C = round(interp1(Cax1,1:length(Cax1),C));


%MAKE FIGURE
%figure
%clf
for ii=0:length(ano)-N
    ii
    valid = valid1;
    indices = find(anos==ano(ii+1));
    [~,ia]=unique(mesa(indices));  %REMOVE REPETITIONS
    indices = indices(ia);
    valid(mesa(indices))=1;  %-valLIM1(1)+1
    c = (C(indices)); 
    C2 = double(valid);
    C2(C2==1)=c;
    cor = colors(indices,:);
    loc = eixo+[ii*stpa ii*stpa -ii*2*stpa -ii*2*stpa]; %location os axis
    a{ii+1} = axes('position',loc); %[0 0 1 1])
    if ii==0
        h = pie(vec1,vec2,label1);
    else
        h = pie(vec1,vec2,null);
    end
    colormap(gca,Ccolors)
    count = 0;
    countc = 0;
    for jj=1:2:length(mes)*2
        count = count+1;
        vali = valid(count);
        h(jj).EdgeColor=[.7 .7 .7];
        if vali==0
            h(jj).FaceColor=[1 1 1];
 %            h(jj).EdgeAlpha=.5;
        else
            countc = countc+1;
            h(jj).FaceColor=cor(countc,:)  ;
            ax = h(jj).Parent ;
            h(jj).CData = h(jj).CData*0+C2(count);
        end
       %Adjust fontsize if wanted:
       set(h(jj+1),'fontsize',8)

    end
    clim(a{ii+1},[1 length(Cax1)]);%length(mes)])
    if ii==length(ano)-N
        h = pie(1,1,'');
        h(1).FaceColor=[1 1 1];
    end
    if ii==0
        %Add Colorbar
        ax = get(h(1),'Parent');
        col = colorbar(ax,'Position',[(eixo(1)+eixo(3))-.05 (eixo(2)+eixo(4))*.85 .01 (eixo(2)+eixo(4))*.2]);
        set(col,'TickLabels',round(Cax1(col.Ticks)))
    end
end


%ADD INNER AXIS LABELS
a0 = a{1}.Position;
%g=axes('position',a0+[a0(3)/2 a0(4)/2+(N+.5)*stpa -a0(3)+0.02 -a0(4)/1.68-(N+.5)*stpa]);%-N*stpa])%+.345/2
g=axes('position',a0+[a0(3)/2-0.02 a0(4)/2+(N+.5)*stpa -a0(3)+0.02 -a0(4)/1.68-(N+.5)*stpa]);%-N*stpa])%+.345/2


plot(zeros(1,length(ano)-N),linspace(0,1,length(ano)-N),'color','none')
xlim([0 1])
yax = linspace(0,1,length(ano)-N+1);
    yaxL = yax(2:2:end);
if mod(length(yax),2)~=0
    yaxL(end) = yax(end);
end
set(g,'YAxisLocation','right')    
set(g,'ydir','reverse','ytick',yaxL,'fontsize',8,'fontweight','bold')
set(g,'color','none','xcolor','none');box off
set(g,'yticklabel',label2(2:2:length(ano)-N+1),'xticklabel','')

% g2=axes('position',a0+[a0(3)/2 a0(4)/2+(N+.5)*stpa -a0(3)+0.02 -a0(4)/1.68-(N+.5)*stpa]);%-N*stpa])%+.345/2
% plot(zeros(1,length(ano)-N),linspace(0,1,length(ano)-N),'color','none')
% xlim([0 1])
% set(g2,'ydir','reverse','ytick',yaxL,'fontsize',8)
% set(g2,'color','none','xcolor','none');box off
% set(g2,'yticklabel',label2(2:2:length(ano)-N+1),'xticklabel','')


disp('Bring me some pies!')
end %function
