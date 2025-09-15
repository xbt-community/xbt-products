%This funtion interpolates the altimetry data to the 
%positions of the xbt data
%Now it accepts 3-D data, as long as you have in the format SAT(lat,lon,z)
%lon and lat have to be matrices (meshgrid)
%enter method(nargin=6) the method: linear,oi,fast
%usage: [lon,lat,data_ini,e_oi]=interp_sat2xbt(vecx,vecy,SAT,lon,lat,method)
function [lon,lat,data_ini,e_oi]=interp_sat2xbt(vecx,vecy,SAT,lon,lat,varargin)
path_i = pwd;
ai = strfind(path_i,'/');
%addpath([path_i(1:ai(2)),'mgoes/matlab/'])

if ischar(vecx),vecx = decyear(str2num(vecx));end
if ischar(vecy),vecy = decyear(str2num(vecy));end
if ischar(lon), vecx = decyear(str2num(lon)); end
if ischar(lat), vecx = decyear(str2num(lat)); end
%Select between cubic spline and optimal interpolation
%Due to the size constraint, the 'oi' divides the depth into 2 parts

paramNames = {'method','modo'};
defaults = {'fast','linear'};
nv = length(varargin);
if nv<=1
  varargin{2} = defaults{2};
end 
if isempty(varargin)
    varargin = defaults;
end
if isempty(varargin{1})
    varargin{1} = defaults{1};
end
%if isempty(varargin{2})
%    varargin{2} = defaults{2};
%end
%keyboard
[method,modo] = internal.stats.parseArgs(paramNames, varargin);

e_oi = 0;  %error of the optimal interpolation
%dbstop('86')
% if (nargin < 6)
% method = 'fast'; %'oi','linear','fast'
% %method = 'linear'; %'oi','linear','fast'
% end
% if (nargin<7)
%     modo = 'linear'
% end
    
sx = length(vecx); sy = length(vecy);
tam_int = size(lon);
lon = lon(:); lat = lat(:);
vecx = vecx(:); vecy = vecy(:);

%Allow SAT with 3 and 4 dimensions(time and depth)%%%%%%%%
%Working for 3-d so far
tam = ndims(SAT);
ax = (size(SAT) == sx);
ay = (size(SAT) == sy);
at = 1; st = 1;
if tam > 2
at = find(~(ax+ay))
st = size(SAT,at);
% if (sum(at) > 3), az = find(at); at = az(1); az = az(2);
%  sz=size(SAT,az);st=size(SAT,at);
% end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(vecx,vecy);
if find(ax) ==1;X=X';Y=Y';end   %ADDED TO AVOID PROBLEMS
data_ini = nan*lat;
switch method

case 'linear',
disp(' Interpolating using linear')
for qt = 1:st
       sat = SAT(:,:,qt);
 %for i=1:length(lat)
     disp(['interpolating cast #',num2str(i)])
%     acha  = find(~isnan(sat) & abs(lat(i)-Y) < 10 & abs(lon(i)-X) <10);
%     dist = sqrt((lat(i)-Y(acha)).^2 + (lon(i)-X(acha)).^2);
 %     data_ini(i,qt) = griddata(X(acha),Y(acha),sat(acha),lon(i),lat(i),'cubic');
      acha  = find(~isnan(sat));% & abs(lat(:)-Y) < 10 & abs(lon(:)-X) <10);
      data_ini(:,qt) = griddata(X(acha),Y(acha),sat(acha),lon,lat,'cubic');

      %end
end
 
case 'fast'
%     disp(' Interpolating using fast interpolation')
      for qt = 1:st
       sat = SAT(:,:,qt);
    if find(ax)==1;
     acha  = find(~isnan(sat(1:sx,1:sy)));
    else
     acha  = find(~isnan(sat(1:sy,1:sx)));
    end
 %     F = TriScatteredInterp(X(acha),Y(acha),sat(acha));
      F = scatteredInterpolant(X(acha),Y(acha),sat(acha),modo);%,'none')%'nearest');

      data_ini(:,qt) = F(lon,lat);
       end

case 'oi',
    disp(' Interpolating using oi')

%%%%%%%%%%%%%%OI parameters%%%%%%%%%%%%%%%%%%%%
%scales = [4*mean(diff(vecx)),4*mean(diff(vecy))];
scales = [2 2];
   err = 0.01; %was 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for qt = 1:st
        sat = SAT(:,:,qt);  
        data_aux = nan* ones(1,length(lat));
     disp(['interpolating depth #',num2str(qt)])
% for i=1:length(lat)
%    acha  = find(~isnan(sat) & abs(lat(i)-Y) < 10 & abs(lon(i)-X) <10);
%[data_aux(i),e] = objmap(X(acha),Y(acha),sat(acha),lon(i),lat(i),...
%[x,y]=meshgrid(
    acha  = find(~isnan(sat)&X<max(lon)+2&...
      X>min(lon)-2&Y<max(lat)+2&Y>min(lat)-2);
[data_aux,e] = objmap(X(acha),Y(acha),sat(acha),lon,lat,...
    scales,err,[32 8 500 1/3 0]);
% end
data_ini(:,qt) = data_aux; 
e_oi(1:length(e),qt) = e;
end

end
lon = lon';
lat = lat';
data_ini = reshape(data_ini,[tam_int qt]);data_ini = squeeze(data_ini);
data_ini = data_ini';   %MAY GIVE PROBLEMS
