function make_npf
transect = {'AX32','AX97','PX30','PX36','PX40'}

%files =  {'Output_final.mat','Output002t.mat','Output.mat','Output_new.mat','Output001.mat'}
files =  {'Output_final.mat','Output002t.mat','Output.mat','Output.mat','Output001.mat'}
%  load([inpath '/Output.mat'])  %PX30 good
%    load([inpath '/Output_new.mat'])  %PX36 good
  %  load([inpath '/Output_final.mat']) %AX32
   % load([inpath '/Output002t.mat'])%AX97  %t for tayanne
  %  load([inpath '/Output001.mat'])%PX40
for ii=1:5
      load(['../',transect{ii},'/',files{ii}])
      npf = npf(:)';
      eval(['npf_',transect{ii},'= npf'])
      if ii==1
         save("npf_goxbt.mat",['npf_',transect{ii}])
      else
         save('-append',"npf_goxbt.mat",['npf_',transect{ii}])
      end         
end
