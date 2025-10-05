function make_mean
transect = {'AX32','AX97','PX30','PX36','PX40'}

%files =  {'Output_final.mat','Output002t.mat','Output.mat','Output_new.mat','Output001.mat'}
files =  {'Output_final.mat','Output002t.mat','Output.mat','Output.mat','Output001.mat'}
%  load([inpath '/Output.mat'])  %PX30 good
%    load([inpath '/Output_new.mat'])  %PX36 good
  %  load([inpath '/Output_final.mat']) %AX32
   % load([inpath '/Output002t.mat'])%AX97  %t for tayanne
  %  load([inpath '/Output001.mat'])%PX40
for ii=3 %5
    ii
      load(['../',transect{ii},'/',files{ii}])
      maxS0 = 10;
      if ii==5
         lat_xbt = OG_lat([1:57 60:3:253],:);
         lon_xbt = OG_lon([1:57 60:3:253],ones(1,size(OG_lat,2)));
      end
      if iscell(lon_xbt);
          nt = length(lon_xbt);
      else 
          nt = size(lon_xbt,2);
      end
      for jj=1:nt
          try
             maxS = length(lon_xbt{jj});
             str1 = '= lon_xbt{select};';
             str2 = '= lat_xbt{select};';
          catch
             maxS = length(lon_xbt(:,jj)); 
             str1 = '= lon_xbt(:,select);';
             str2 = '= lat_xbt(:,select);';
          end
          if maxS>maxS0
              select = jj;
              maxS0=maxS;
          end
      end
           eval(['lon_',transect{ii},str1])
           eval(['lat_',transect{ii},str2])
      %     lon_PX40 = lon_PX40([1:57 60:3:253]);
      %     lat_PX40 = lat_PX40([1:57 60:3:253]);
      
      
      if ii==1
         save("latlon_goxbt.mat",['lon_',transect{ii}],['lat_',transect{ii}])
      else
         save('-append',"latlon_goxbt.mat",['lon_',transect{ii}],['lat_',transect{ii}])
      end         
end
