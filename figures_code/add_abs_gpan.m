function [gpan,addep2] = add_abs_gpan(gpan,lat2,lon2,Pn,pref,month_in)
%NOW ACCEPTS 3-D filed
%DEPTH HAS TO BE FIRST DIMENSION
%by Marlos Goes Aug 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizeg = size(gpan);
[~,dep] = min(abs(Pn-pref));

ndep = length(Pn);

if exist('month_in','var')

    month = month_in;
    
    addep = load_abs_gpan_section(lat2,lon2,Pn,month);
    
    
    
    
    
    
else
    
    
    addep = load_abs_gpan_section(lat2,lon2,Pn);
    
    
    
    
    
    
end

disp('add') %,keyboard
gpan = ones(ndep,1)*gpan(1,:)-gpan(:,:);              
gpanDep = gpan(dep,:);
gpanDep = fillmissing(gpanDep,'linear');       % disp('Updated')
gpan = ones(ndep,1)*gpanDep-gpan(:,:); 
gpan = gpan + ones(ndep,1)*addep(dep,:);

gpan = reshape(gpan,sizeg);
addep2 = reshape(addep,sizeg);

return
