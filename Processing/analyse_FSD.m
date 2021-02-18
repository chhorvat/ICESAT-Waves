% function analyse_FSD(file_path,outdir,gridname)
%%
% A pointer for the .mat file with IS2 data.
IS2data = matfile(file_path);
% [ntracks,nfields] = size(IS2data.fields);

% Output -
% frac_under_geo nlat x nlon - fraction (by length) of measurements made in
% each lat-lon bin that are wave affected. Equal to twice the length of the
% measurements that are negative divided by the length of all measurements

% total_under_geo - nlat x nlon - fraction (of number) " . Equal to twice
% the number of all segments in each bin that are negative.

% field_names = ...
%     { ... % Get the ice segment information
%     'freeboard_beam_segment/height_segments/height_segment_height' ...
%     'freeboard_beam_segment/height_segments/height_segment_length_seg', ...
%     'freeboard_beam_segment/height_segments/latitude', ...
%     'freeboard_beam_segment/height_segments/longitude', ...
%     'freeboard_beam_segment/height_segments/photon_rate', ...
%     'freeboard_beam_segment/beam_freeboard/beam_fb_height', ...
%     'freeboard_beam_segment/height_segments/height_segment_type', ...
%     'freeboard_beam_segment/height_segments/height_segment_ssh_flag', ...
%     'freeboard_beam_segment/height_segments/ice_conc', ...
%     };


% Options

% Which beams will we choose? Index is [gt1r,gt1l,gt2r,gt2l,gt3r,gt3l]

% Use CESM grid
% lat = ncread('cesm_grid_area.nc','TLAT'); lon = ncread('cesm_grid_area.nc','TLON'); lat(isnan(lat)) = 1e6; lon(isnan(lon)) = 1e6;

%% Because we don't care about the tracks themselves, we want to create a
% long vector which has lat/lon/elevation. We do this by concatenating all
% of the beams together.
disp(['Loading...']);

conc = IS2data.fields(:,end);
% Looking now but not using ssh fields
fieldmat = IS2data.fields(:,1:end-2);
load(file_path,'timer');

disp('Loaded');

% Fieldmat is nseg x nfields - for every segment we have
% latitude/longitude/elevation/segment length/photon rate/freeboard

%% Remove places with local low surface heights

% Identify local moving average
window_1k = 1000;
window_10k = 10000; % meters - size of moving window
window_50k = 50000;
max_seg_size = 200;

earthellipsoid = referenceSphere('earth','m');
sortvec = [];
binary_floelength = []; 
binary_floeid = []; 
moving_floenum = []; 
moving_floelength = []; 
moving_RCL = []; 
moving_MCL = []; 

numtracks = length(timer);

%%
%
% Compute distance along the track, in units of m.
for i = 1:numtracks
    %%
    lat = fieldmat{i,3};
    lon = fieldmat{i,4};
    
    if length(lat) > 1
        % along-track distance
        dist = distance([lat(1:end-1) lon(1:end-1)],[lat(2:end) lon(2:end)],earthellipsoid);
        
    else
        
        dist = [];
        
    end
    
    if ~isempty(dist)
        
        dist = [0; cumsum(dist)];
        
    end
    
    % Preprocess
    dupes = find(diff(dist)==0)+1;
    dist(dupes) = dist(dupes) + .01*abs(rand(size(dupes)));
    [dist,b] = sort(dist);
    
    %     ssh_flag = fieldmat{i,8};
    %     ssh_flag = ssh_flag(b);
    
    isice = fieldmat{i,7};
    isice = isice(b);
    
    isocean = isice > 1;
    
    %% Objects that are evaluated locally
    
    % Sorting to put in order
    sortvec = cat(1,sortvec,b+length(sortvec));
    
    
    %% FSD-related. 
    up = strfind([0,isocean'],[0 1]);
    down = strfind([isocean',0],[1 0]);
    
    if ~isempty(up)
        
        toosmall = intersect(up,down); 
        up = setxor(up,toosmall)'; 
        down = setxor(down,toosmall)'; 
        
        Uloc = up - 1; 
        Uloc(Uloc==0) = 1; 
        
        floelen = dist(down) - dist(up); 
        floeloc = .5*(dist(down) + dist(up));        
        floeind = round(.5*(down + up));
        floe_seglength = floelen./(down - up); 
                
    end
    
    usable_floe = logical((floelen > 5).*(floe_seglength < 100)); 
    usable_floe(1) = 0; % exclude endpoints
    usable_floe(end) = 0; % exclude endpoints
    
    floelen = floelen(usable_floe); 
    floe_seglength = floe_seglength(usable_floe); 
    floeloc = floeloc(usable_floe); 
    floeind = floeind(usable_floe); 
    
    % One at the location of a floe
    % Zero otherwise
    hasfloe = 0*dist;
    hasfloe(floeind) = 1; 
    
    % Length of floe, at the center location of each floe. Should find a
    % way to add these up to avoid undercounting... 
    floe_lengths = 0*dist; 
    floe_lengths(floeind) = floelen; 
       
    % Map of floe centers and floe lengths. 
    binary_floeid = cat(1,binary_floeid,hasfloe); 
    binary_floelength = cat(1,binary_floelength,floe_lengths);
    moving_floenum = cat(1,moving_floenum,movsum(hasfloe,window_50k,'samplepoints',dist));
    moving_floelength = cat(1,moving_floelength,movsum(floe_lengths,window_50k,'samplepoints',dist));
    
    % Moving mean chord length (int r * N)/ int(N)
    moving_MCL = cat(1,moving_MCL, ...
        movsum(floe_lengths,window_50k,'samplepoints',dist) ...
        ./ movsum(hasfloe,window_50k,'samplepoints',dist)); 
    
    % Moving representative chord length (int r * f / int f)
    moving_RCL = cat(1,moving_RCL, ...
        movsum(floe_lengths.*floe_lengths.^2,window_50k,'samplepoints',dist) ...
        ./ movsum(floe_lengths.^2,window_50k,'samplepoints',dist));
    
    
end

%%

% Get into a matrix and then sort it!
fieldmat = cell2mat(fieldmat);
fieldmat = fieldmat(sortvec,:);

%% FSD Code First

% Criteria for inclusion of a location. Kind of already did it. 
fprintf('\nMeasured %2.0f thousand chord lengths \n',sum(binary_floeid)/1000);
fprintf('Total chord length is %2.0f km \n',sum(binary_floelength)/1000);
fprintf('Average chord length (number) is %2.0f m \n',sum(binary_floelength)/sum(binary_floeid));

smallfloe_MCL = moving_MCL < 100; 
midfloe_MCL = (moving_MCL >= 100) & (moving_MCL <= 1000); 
bigfloe_MCL = (moving_MCL > 1000); 

smallfloe_RCL = moving_RCL < 100; 
midfloe_RCL = (moving_RCL >= 100) & (moving_RCL <= 1000); 
bigfloe_RCL = (moving_RCL > 1000); 

%% Now we want to figure out where the waves are

numsegs = size(fieldmat,1);

fprintf('Have %d total segments to analyse \n',numsegs)

% What grid will we bin the lat-lon values into

kdloc = ['../Processing/KDTrees/KDTree_' gridname];

load(kdloc,'lat_X','lon_X','KDTree');

disp('Finding Locations');

% K nerest neighbor search into the loaded grid to find grid locations
posloc = knnsearch(KDTree,[fieldmat(:,3) fieldmat(:,4)],'K',1);

%% Accumulate all values that are below zero into the chosen lat-lon array

floenum_geo = accumarray(posloc,binary_floeid,[numel(lat_X) 1],@sum);
floelength_geo = accumarray(posloc,binary_floelength,[numel(lat_X) 1],@sum);
MCL_geo = floelength_geo ./ floenum_geo; 

smallfloe_MCL_geo = accumarray(posloc,smallfloe_MCL,[numel(lat_X) 1],@mean);
midfloe_MCL_geo = accumarray(posloc,midfloe_MCL,[numel(lat_X) 1],@mean);
bigfloe_MCL_geo = accumarray(posloc,bigfloe_MCL,[numel(lat_X) 1],@mean);

smallfloe_RCL_geo = accumarray(posloc,smallfloe_RCL,[numel(lat_X) 1],@mean);
midfloe_RCL_geo = accumarray(posloc,midfloe_RCL,[numel(lat_X) 1],@mean);
bigfloe_RCL_geo = accumarray(posloc,bigfloe_RCL,[numel(lat_X) 1],@mean);

%%

disp(['Saving ' outdir]);

% save(outdir,'lat_X','lon_X','*_geo','numsegs','numtracks');