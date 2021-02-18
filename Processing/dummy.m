% function analyse_wave_frac_affected(file_path,outdir,gridname)
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

%% Wave-related things
moving_wavy = [];

sortvec = [];
moving_ssh_no = [];


%% FSD-related things
binary_floelength = [];
binary_floeid = [];
moving_floenum = [];
moving_floelength = [];
moving_CLD_mom_0 = [];
moving_CLD_mom_1 = [];
moving_CLD_mom_2 = [];
moving_CLD_mom_3 = [];

%%
numtracks = length(timer);

%
% Compute distance along the track, in units of m.
for i = 1:50;%numtracks
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
    
    is_ice = fieldmat{i,7};
    is_ice = is_ice(b);
    
    is_ocean = is_ice > 1;
    
    height = fieldmat{i,1};
    height = height(b);
    
    seg_len = fieldmat{i,2};
    seg_len = seg_len(b);
    
    % check to make sure all ice points are close to open water points
    % close_to_ssh = cellfun(@length,rangesearch(dist(isocean),dist,window_ssh));
    
    % This is the index of the nearest neighbor to each point that has open
    % water
    
    if isempty(dist(is_ocean))
        
        dist_to_ssh = nan*dist;
        
    else
        
        [~,dist_to_ssh] = knnsearch(dist(is_ocean),dist);
    end
    
    if sum(is_ocean) > 1
        
        % This is the interpolated SSH field, on all of the dist places
        ssh_interp = interp1(dist(is_ocean),height(is_ocean),dist);
        
        
    else
        
        ssh_interp = nan*dist;
        
    end
    
    %% Wave Diagnostics - many repeated in analyse_wave_frac_affected
    
    % First eliminate all segments larger than 1 km - may be artificially big
    not_too_long = seg_len < max_seg_size;
    
    ssh_window_no = movsum(is_ocean,window_10k,'samplepoints',dist);
    ssh_window_std = movstd(ssh_interp,window_10k,'samplepoints',dist);
    height_window_std = movstd(height,window_10k,'samplepoints',dist);
    isneg = height - ssh_interp < 0;
    
    window_pos = movsum(~isneg,window_1k,'samplepoints',dist); 
    window_neg = movsum(isneg,window_1k,'samplepoints',dist);

    
    % Need to have two ssh points within 10 km
    close_to_ssh = (ssh_window_no >= 2);%  .* (ssh_neighbors <= window_10k);
    
    % Minimum threshold on moving std
    wave_cutoff_ssh = max(ssh_window_std,.1);
    wave_cutoff_height = max(height_window_std,.1);
    both_cutoff_height = max(wave_cutoff_ssh,wave_cutoff_height);
    
    % adjust for deviation from local ssh
    height_adjusted = height - ssh_interp;
    
    % Has positive points nearby
    close_to_positive = window_pos >= 2;
    % Has negative points nearby - only for those values that are negative
    close_to_negative = window_neg >= 2;
    
    % Included points have positive points nearby, aren't too long, and are
    % identified as ice. Criteria I1-I2.
    is_included = logical(not_too_long .* is_ice .* close_to_positive);
    
    % Wave candidates are those i|ncluded segs close to ssh points, and
    % close to other negative values.
    is_wave_candidate = logical(is_included .*close_to_ssh.* close_to_negative); % - XX last part kills most tracks!?
    
    % Height is negative beyond ssh or height variance
    is_under_both_var = logical((height_adjusted < -both_cutoff_height).*is_wave_candidate);
    
    moving_wavy = cat(1,moving_wavy,is_under_both_var);    
    
    %% FSD-related things
    up = strfind([0,is_ocean'],[0 1]);
    down = strfind([is_ocean',0],[1 0]);
    
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
    
    % Sum up moments of the floe length distsribution in each 50k window
    mom0 =  movsum(hasfloe,window_50k,'samplepoints',dist);
    mom1 = movsum(floe_lengths,window_50k,'samplepoints',dist);
    mom2 = movsum(floe_lengths.^2,window_50k,'samplepoints',dist);
    mom3 = movsum(floe_lengths.^3,window_50k,'samplepoints',dist);
    
    
    moving_CLD_mom_0 = cat(1,moving_CLD_mom_0,mom0);
    moving_CLD_mom_1 = cat(1,moving_CLD_mom_1,mom1);
    moving_CLD_mom_2 = cat(1,moving_CLD_mom_2,mom2);
    moving_CLD_mom_3 = cat(1,moving_CLD_mom_3,mom3);
    
end

%% FSD Code First

% Criteria for inclusion of a location. Kind of already did it.
fprintf('\nMeasured %2.0f thousand chord lengths \n',sum(binary_floeid)/1e3);
fprintf('Total chord length is %2.0f km \n',sum(binary_floelength)/1000);
fprintf('Average chord length (number) is %2.0f m \n',sum(binary_floelength)/sum(binary_floeid));

moving_MCL = moving_CLD_mom_1 ./ moving_CLD_mom_0;
moving_RCL = moving_CLD_mom_3 ./ moving_CLD_mom_0;

smallfloe_MCL = moving_MCL < 100;
midfloe_MCL = (moving_MCL >= 100) & (moving_MCL <= 1000);
bigfloe_MCL = (moving_MCL > 1000);

smallfloe_RCL = moving_RCL < 100;
midfloe_RCL = (moving_RCL >= 100) & (moving_RCL <= 1000);
bigfloe_RCL = (moving_RCL > 1000);

%% Now we want to figure out where the places go

numsegs = size(fieldmat,1);

fprintf('Have %d total segments to analyse \n',numsegs)

% What grid will we bin the lat-lon values into

kdloc = ['../Processing/KDTrees/KDTree_' gridname];

load(kdloc,'lat_X','lon_X','KDTree');

disp('Finding Locations');

% K nerest neighbor search into the loaded grid to find grid locations
posloc = knnsearch(KDTree,[fieldmat(:,3) fieldmat(:,4)],'K',1);

%% Wave-FSD Geographic Fields
 
% Add up number of places with wave index = 1
FSD_wave_geo = accumarray(posloc,moving_wavy,[numel(lat_X) 1],@sum);
% Number of segments in each
FSD_segnum_geo = accumarray(posloc,1 + 0*moving_wavy,[numel(lat_X) 1],@sum);
% Fraction of segments with WAVY
FSD_wavefrac_geo = FSD_wave_geo ./ FSD_segnum_geo; 
%% FSD Geographic Fields

% Number of floes at each point
floenum_geo = accumarray(posloc,binary_floeid,[numel(lat_X) 1],@sum);
% Total floe length at each point
floelength_geo = accumarray(posloc,binary_floelength,[numel(lat_X) 1],@sum);
% Mean chord length at each point

CLD_mom_0_geo = accumarray(posloc,moving_CLD_mom_0,[numel(lat_X) 1],@sum);
CLD_mom_1_geo = accumarray(posloc,moving_CLD_mom_1,[numel(lat_X) 1],@sum);
CLD_mom_2_geo = accumarray(posloc,moving_CLD_mom_2,[numel(lat_X) 1],@sum);
CLD_mom_3_geo = accumarray(posloc,moving_CLD_mom_3,[numel(lat_X) 1],@sum);

MCL_geo = CLD_mom_1_geo ./ CLD_mom_0_geo;
RCL_geo = CLD_mom_3_geo ./ CLD_mom_2_geo;

% Breakdown of % where made up of small, medium, or large floes. Useful
% when comparing areas with lots of tracks on top of them.
smallfloe_MCL_geo = accumarray(posloc,smallfloe_MCL,[numel(lat_X) 1],@mean);
midfloe_MCL_geo = accumarray(posloc,midfloe_MCL,[numel(lat_X) 1],@mean);
bigfloe_MCL_geo = accumarray(posloc,bigfloe_MCL,[numel(lat_X) 1],@mean);

smallfloe_RCL_geo = accumarray(posloc,smallfloe_RCL,[numel(lat_X) 1],@mean);
midfloe_RCL_geo = accumarray(posloc,midfloe_RCL,[numel(lat_X) 1],@mean);
bigfloe_RCL_geo = accumarray(posloc,bigfloe_RCL,[numel(lat_X) 1],@mean);

%%

disp(['Saving ' outdir]);

save(outdir,'lat_X','lon_X','*_geo','numsegs','numtracks');
