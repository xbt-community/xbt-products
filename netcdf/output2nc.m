function output2nc(line_id,inpath,outfname)
% function output2nc(line_id,infname,repopath,outfname)
% Inputs:   line_id (eg, PX30) as character array
%           inpath: full path to location of 'Output.mat' file from XBT
%           gridding
%           outfname (optional): full path and output the netcdf gridded
%           data to. Default is the same as inpath.
%
% Loads the Output.mat file containing D_pr structure with gridded data 
%ll_grid (lat/long grid) and pr_grid (pressure grid)
% created from figure1_<line_id>.m
% Outputs gridded data to CF compliant NETCDF format
%
% Using text files for easy editing of metadata attributes.
% Bec Cowley, August, 2025.

% the file is created in the following order
%
% 1. global attributes
% 2. dimensions / coordinate variables
% 3. variable definitions
% 4. data

% directory in upper case, filname in lower case (as of ver.1)
LINE = upper(line_id);

narginchk(2,3);
if nargin < 3
    outfname = inpath;
end
outfname = [outfname '/' LINE '_gridded.nc'];

%load the input file:
try
    load([inpath '/Output.mat'])
catch
    error(['File ' inpath '/Output.mat does not exist'])
end

%Netcdf file creation
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLOBBER'));
fidnc = netcdf.create(outfname, cmode);
if fidnc == -1, error(['Could not create ' filename]); end

% we don't want the API to automatically pre-fill with FillValue, we're
% taking care of it ourselves and avoid 2 times writting on disk
netcdf.setFill(fidnc, 'NOFILL');

%% get the global attributes template:
globalatts = parseNCTemplate('global_attributes_gridded.txt');

%populate empty fields where we can:
%list of fields
flds = fieldnames(globalatts);

%add the line name to the end of the text
%handle the 'goship_woce_line_id' here also
globalatts.title = [globalatts.title ' ' LINE];

%'date_issued'
globalatts.date_issued = datestr(now,'yyyymmdd');
%'soop_line_id'
%input into function
globalatts.soop_line_id = LINE;
% 'geospatial_bounds'
%also handle all the geospatial* etc here

%Final output needs to be -180 to 180 longitude, do conversion here
ineg = lon_ref > 180;
lon_ref(ineg) = lon_ref(ineg)-360;
% decide if this is a transect that runs north-south or east-west
[latm,latn] = size(lat_ref);
[lonm, lonn] = size(lon_ref);
if lonn < latn
    % east-west line
    ilatlon = 1;
else
    % north-south line
    ilatlon = 2;
end

%assign the global attributes:
globalatts.geospatial_lat_min = min(min(lat_ref));
globalatts.geospatial_lat_max = max(max(lat_ref));
globalatts.geospatial_lon_min = min(min(lon_ref));
globalatts.geospatial_lon_max = max(max(lon_ref));
globalatts.geospatial_vertical_min = min(depth_ref);
globalatts.geospatial_vertical_max = max(depth_ref);
%also handle the other time* fields here
mintime = NaN*time_avg;
maxtime = mintime;
for a = 1:length(time_xbt)
    mintime(a) = min(time_xbt{a});
    maxtime(a) = max(time_xbt{a});
end
globalatts.time_coverage_start = datestr(min(mintime),'dd-mm-yyyy HH:MM:SS');
globalatts.time_coverage_end = datestr(max(maxtime),'dd-mm-yyyy HH:MM:SS');

%write out global attributes:
vid = netcdf.getConstant('NC_GLOBAL'); %get the global attributes reference
flds = fieldnames(globalatts);
for a = 1:length(flds)
    name = flds{a};
    netcdf.putAtt(fidnc, vid, name, globalatts.(flds{a}));
end
%% get the dimensions set up
%dimensions
%set up the section data:
sect = 1:length(time_avg);

%if east-west line: 
if ilatlon == 1
    dimnames = {'time','longitude','depth'};
    dimdata = {'time_avg','lon_ref','depth_ref'};
else
    %north-south
    dimnames = {'time','latitude','depth'};
    dimdata = {'time_avg','lat_ref','depth_ref'};
end
time_avgatts = parseNCTemplate('time_attributes.txt');
depth_refatts = parseNCTemplate('depth_attributes_gridded.txt');
lon_refatts = parseNCTemplate('longitude_attributes.txt');
lat_refatts = parseNCTemplate('latitude_attributes.txt');

for m=1:length(dimnames)
    eval(['data = ' dimdata{m} ';']);
    % create dimension
    did(m) = netcdf.defDim(fidnc, dimnames{m}, length(data));
    % create coordinate variable and attributes
    eval(['atts = ' dimdata{m} 'atts;']);
    vid(m) = netcdf.defVar(fidnc, dimnames{m}, 'NC_DOUBLE', did(m));
    fldn = fieldnames(atts);
    for b = 1:length(fldn)
        netcdf.putAtt(fidnc,vid(m),fldn{b},atts.(fldn{b}))
    end
end
%% now for each variables attributes
if ilatlon ==1
    % east-west
    loc_var = 'latitude';
    loc_name = 'lat_ref';
    loc_units = 'degrees_east';
    loc_range = [-90, 90];
else
    %north-south
    loc_var = 'longitude';
    loc_name = 'lon_ref';
    loc_units = 'degrees_north';
    loc_range = [-180, 180];
end
%populate for each one:
varname = {loc_var,'temperature','reference_salinity','velocity','surface_dynamic_height',...
    'sea_surface_height','geostrophic_transport','altimetric_ssh_gradient'...
    };
dataname = {loc_name,'temp_ref','salt_ref','vel_ref','dh0_ref','ssh_ref_px30',...
    'GStransp', 'GStransp_alt'};
stdn = {loc_var,'sea_water_temperature','sea_water_practical_salinity',...
    'geostrophic_northward_sea_water_velocity','','sea_surface_height','',...
    ''};
long_name = {loc_var,'sea_water_temperature','sea_water_practical_salinity',...
    'geostrophic_northward_sea_water_velocity',...
    'absolute_dynamic_height_at_surface','sea_surface_height_altimetry','mass_transport_from_XBT',...
    'cross_front_delta_ssh'};
units = {loc_units,'degC','1','m s-1','m','m','Sv','m'}; 
refscale = {'WGS84 geographic coordinate system','ITS-90','PSS-78','','','','',''};
vmin = [loc_range(1),-2.5,2.0,-10,-10,-10,-50,-10];
vmax = [loc_range(2),40.0,41.0,10,10,10,50,10];
for a = 1:length(varname)
    clear varatts
    if ~isempty(stdn{a})
        varatts.standard_name = stdn{a};
    end
    varatts.long_name = long_name{a};
    if a <=4
        varatts.coordinates = [dimnames{1},' ', dimnames{2},' ' dimnames{3}] ;
    else
        varatts.coordinates = dimnames{1};
    end
    varatts.units = units{a};
    varatts.valid_min = vmin(a);
    varatts.valid_max = vmax(a);    
    varatts.FillValue = NaN;
    if ~isempty(refscale{a}) %do for all but time and oxygen
        varatts.reference_scale = refscale{a};
    end
    if a == 7 | a == 8% time only dimensional var
        varatts.coordinates = dimnames{1};
        didv = did(1);
    elseif a == 1 |a == 5 | a ==6 %variables,2 dims
        varatts.coordinates = [dimnames{1},' ', dimnames{2}] ;
        didv = did([2,1]);
    else % all three dimensions
        varatts.coordinates = [dimnames{1},' ', dimnames{2}, ' ', dimnames{3}] ;
        didv = fliplr(did);
    end
    
    %create the variable:
    vidv(a) = netcdf.defVar(fidnc,varname{a},'NC_DOUBLE',didv);
    fldn = fieldnames(varatts);
    %now write out the variable atts:
    for b = 1:length(fldn)
        name = fldn{b};
        if strcmpi(name, 'FillValue')
            netcdf.defVarFill(fidnc, vidv(a), false, NaN); % false means noFillMode == false
        else
            netcdf.putAtt(fidnc, vidv(a), name, varatts.(fldn{b}));
        end
    end
end
% we're finished defining dimensions/attributes/variables
netcdf.endDef(fidnc);

%% Now the data
% dimension data
for a = 1:length(dimnames)
    eval(['data = ' dimdata{a} ';']);
    netcdf.putVar(fidnc, vid(a), data);
end

%variable data
for a = 1:length(varname)
    eval(['data = ' dataname{a} ';']);
    disp(dataname{a})
       
    % reshape velocity, add one more time stamp of empty values at the
    % start to match the lon_ref or lat_ref
    if a == 4
       [row,col,pages]=size(data);
       nan_row = NaN(row,1,pages);
       data = [nan_row,data];
    end

    %data needs to be kept to 3 decimal places
    data = round(data,3);

    netcdf.putVar(fidnc, vidv(a), data);
end
netcdf.close(fidnc)
end
