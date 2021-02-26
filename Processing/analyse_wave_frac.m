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
height_moving_avg = [];
height_moving_std = [];
moving_en = [];
ssh = [];
ssh_moving_std = [];
ssh_neighbors = [];
moving_neg = []; 
moving_pos = []; 

sortvec = [];
moving_ssh_no = []; 

lenct = 0; 

%%
numtracks = length(timer); 

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
    % Remove duplicate values
    dupes = find(diff(dist)<0.5)+1;
    dist(dupes) = []; 
    
    % dist(dupes) = dist(dupes) + .01*abs(rand(size(dupes)));
    [dist,b] = sort(dist);
    
    % 
    % Dedupe and sort ice vector
    is_ice = fieldmat{i,7};
    is_ice(dupes) = []; 
    is_ice = is_ice(b);
    
    % Ocean is the stuff that isn't ice. 
    is_ocean = is_ice > 1;
    
    % Sorting to put in order. Have indices of duplicate values first
    % Keep in mind where in the initial vector we had a duplicate
    dupevec = cat(1,dupevec,dupes + lenct); 
   
    % This adds the dupe vector by the size of the field
    lenct = lenct + size(fieldmat{i,1},1); 

    is_ocean = is_ice > 1;
    
    height = fieldmat{i,1};
    height(dupes) = []; 
    height = height(b);
    
    seg_len = fieldmat{i,2};
    seg_len(dupes) = []; 
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
    
    %% Objects that are evaluated locally
    
    % How far away is the nearest SSH point. 
    ssh_neighbors = cat(1,ssh_neighbors,dist_to_ssh);
    
    % Sorting to put in order
    sortvec = cat(1,sortvec,b+length(sortvec));

    % Index of points in this window
    pts = length(ssh_neighbors)-length(seg_len)+1:length(ssh_neighbors); 
    
    %% Objects interpolated on the 50km moving window 
    % Total wave energy and SSH
    
    % SSH 
    ssh = cat(1,ssh,ssh_interp); % movmean(ssh_interp,window_50k,'samplepoints',dist));

    % Moving average wave energy
    % Number of SSH points within that window 
    
    moving_en = cat(1,moving_en,movmean((height-ssh(pts)).^2 .* seg_len,window_50k,'samplepoints',dist));

    
    %% Objects interpolated on the 10km moving window 
    % Standard deviation of Height
    
    % height_moving_avg = cat(1,height_moving_avg,movmean(height,window_variance,'samplepoints',dist));
    height_moving_std = cat(1,height_moving_std,movstd(height,window_10k,'samplepoints',dist));

    moving_ssh_no = cat(1,moving_ssh_no,movsum(is_ocean,window_10k,'samplepoints',dist));
    
    ssh_moving_std = cat(1,ssh_moving_std,movstd(ssh_interp,window_10k,'samplepoints',dist));

    
    %% Objects interpolated on the 1km moving window

    % Negative values of adjusted SSH
    isneg = height - ssh(pts) < 0; 
    
    moving_neg = cat(1,moving_neg,movsum(isneg,window_1k,'samplepoints',dist));
    moving_pos = cat(1,moving_pos,movsum(~isneg,window_1k,'samplepoints',dist)); 
        
     
end

%%

% Get into a matrix and then sort it!
fieldmat = cell2mat(fieldmat);
fieldmat(dupevec,:) = []; 
fieldmat = fieldmat(sortvec,:);

%% Now Wave Code

% Multiple for WAL
Mval = (.5 - (1/pi)*asin(ssh./moving_en + 1/sqrt(2))).^(-1); 
Mval(abs(imag(Mval)) > 0) = 0; 
Mval(isnan(Mval)) = 0; 
Mval = min(Mval,10); 

% Compute the fraction of all measurements (by length or number) that are
% "wave-affected"

% First eliminate all segments larger than 1 km - may be artificially big
not_too_long = fieldmat(:,2)<max_seg_size;

% Need to have two ssh points within 10 km
close_to_ssh = (moving_ssh_no >= 2);%  .* (ssh_neighbors <= window_10k);

% Minimum threshold on moving std
wave_cutoff_ssh = max(ssh_moving_std,.1); 
wave_cutoff_height = max(height_moving_std,.1); 
both_cutoff_height = max(wave_cutoff_ssh,wave_cutoff_height); 

% ice if tagged as sea ice by ATL07
is_ice = fieldmat(:,7) == 1;

% adjust for deviation from local ssh
height_adjusted = fieldmat(:,1) - ssh;

% Not too long, is ice, and has positive moving average
naive_reasonable = logical(not_too_long .* is_ice);

% Has positive points nearby
close_to_positive = moving_pos >= 2; 
% Has negative points nearby - only for those values that are negative
close_to_negative = moving_neg >= 2; 

% Included points have positive points nearby, aren't too long, and are
% identified as ice. Criteria I1-I2.
is_included = logical(not_too_long .* is_ice .* close_to_positive);

% Wave candidates are those i|ncluded segs close to ssh points, and
% close to other negative values. 
is_wave_candidate = logical(is_included .*close_to_ssh.* close_to_negative); % - XX last part kills most tracks!?

% Different metrics for being wave-affected

% Just negative, ice, and not too long. 
naive_under = logical((fieldmat(:,1) < 0).*is_wave_candidate);

% Adjusted height is negative
is_under = logical((height_adjusted < 0).*is_wave_candidate);

% Height is negative beyond ssh or height variance
is_under_ssh_var = logical((height_adjusted < -wave_cutoff_ssh).*is_wave_candidate);
is_under_height_var = logical((height_adjusted < -wave_cutoff_height).*is_wave_candidate);
is_under_both_var = logical((height_adjusted < -both_cutoff_height).*is_wave_candidate);

% Take total length of negative elevations divided by total length of all
% segments. Multiply by two because all troughs have crests.
wave_area_frac_naive = 2 * sum(fieldmat(naive_under,2)) / sum(fieldmat(is_ice,2));
wave_area_frac = 2 * sum(fieldmat(is_under,2)) / sum(fieldmat(is_ice,2));
wave_area_frac_ssh = 2 * sum(fieldmat(is_under_ssh_var,2)) / sum(fieldmat(is_ice,2));
wave_area_frac_height = 2 * sum(fieldmat(is_under_height_var,2)) / sum(fieldmat(is_ice,2));
wave_area_frac_both = 2 * sum(fieldmat(is_under_both_var,2)) / sum(fieldmat(is_ice,2));
wave_area_frac_M_height = 2 * sum(Mval(is_under_height_var).*fieldmat(is_under_height_var,2)) / sum(fieldmat(is_ice,2));

% Take total number of negative elevations divided by total number of all
% segments
wave_num_frac_naive = 2 * sum(naive_under)/sum(is_ice);
wave_num_frac = 2 * sum(is_under)/sum(is_included);
wave_num_frac_ssh = 2 * sum(is_under_ssh_var)/sum(is_included);
wave_num_frac_height = 2 * sum(is_under_height_var)/sum(is_included);
wave_num_frac_both = 2 * sum(is_under_both_var)/sum(is_included);
wave_num_frac_M_height = 2 * sum(Mval.*is_under_height_var)/sum(is_included);

fprintf('\nTotal number of ice segs: %2.1f million \n',sum(is_ice)/1e6)
fprintf('Total number of not-too-long ice segs: %2.1f million \n',sum(naive_reasonable)/1e6)
fprintf('Total number of not-too-long ice segs with positive MA: %2.1f million \n',sum(is_included)/1e6)
fprintf('Total number of wave_candidate ice segs: %2.1f million \n',sum(is_wave_candidate)/1e6)
fprintf('%2.1f percent of all ice segs (by number) have negative heights \n',wave_num_frac_naive*100);
fprintf('%2.1f percent of suitable segs (by number) have negative heights \n',wave_num_frac*100);
fprintf('%2.1f percent of suitable segs (by number) have statistically outlier negative heights vs SSH variance \n',wave_num_frac_ssh*100);
fprintf('%2.1f percent of suitable segs (by number) have statistically outlier negative heights vs height variance \n',wave_num_frac_height*100);
fprintf('%2.1f percent of suitable segs (by number) have statistically outlier negative heights vs both variance \n',wave_num_frac_both*100);
fprintf('%2.1f percent of M-adjusted suitable segs (by number) have statistically outlier negative heights \n',wave_num_frac_M_height*100);

fprintf('\nTotal length of ice segs: %2.1f km \n',sum(fieldmat(is_ice,2))/1000)
fprintf('Total length of not-too-long ice segs: %2.1f km \n',sum(fieldmat(logical(is_ice.*not_too_long),2))/1000)
fprintf('Total length of not-too-long ice segs with positive MA: %2.1f km \n',sum(fieldmat(is_included,2))/1000)
fprintf('Total length of wave candidate ice segs: %2.1f km \n',sum(fieldmat(is_wave_candidate,2))/1000)
fprintf('%2.1f percent of all ice segs (by length) have negative heights \n',wave_area_frac_naive*100);
fprintf('%2.1f percent of suitable segs (by length) have negative heights \n',wave_area_frac*100);
fprintf('%2.1f percent of suitable segs (by length) have statistically outlier negative heights vs SSH variance \n',wave_area_frac_ssh*100);
fprintf('%2.1f percent of suitable segs (by length) have statistically outlier negative heights vs Height variance \n',wave_area_frac_height*100);
fprintf('%2.1f percent of suitable segs (by length) have statistically outlier negative heights vs Height variance \n',wave_area_frac_both*100);
fprintf('%2.1f percent of M-adjusted suitable segs (by length) have statistically outlier negative heights \n',wave_area_frac_M_height*100);

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

if ~isempty(conc)
    
    conc_geo = accumarray(posloc,is_included,[numel(lat_X) 1],@mean);
    
end

% This adds up all measurements that are reasonable
num_meas_geo = accumarray(posloc,is_included,[numel(lat_X) 1],@sum);
naive_num_meas_geo = accumarray(posloc,naive_reasonable,[numel(lat_X) 1],@sum);

% This adds up the length of all measurements
len_meas_geo = accumarray(posloc,is_included.*fieldmat(:,2),[numel(lat_X) 1],@sum);
naive_len_meas_geo = accumarray(posloc,naive_reasonable.*fieldmat(:,2),[numel(lat_X) 1],@sum);

% This is all measurements below zero (times 2)
num_under_geo = 2*accumarray(posloc,is_under,[numel(lat_X) 1],@sum);
naive_num_under_geo = 2*accumarray(posloc,naive_under,[numel(lat_X) 1],@sum);
num_ssh_under_geo = 2*accumarray(posloc,is_under_ssh_var,[numel(lat_X) 1],@sum);
num_height_under_geo = 2*accumarray(posloc,is_under_height_var,[numel(lat_X) 1],@sum);
num_both_under_geo = 2*accumarray(posloc,is_under_both_var,[numel(lat_X) 1],@sum);

% Adds up length of all that are negative (times 2) in each cell
len_under_geo = 2*accumarray(posloc,is_under.*fieldmat(:,2),[numel(lat_X) 1],@sum);
naive_len_under_geo = 2*accumarray(posloc,naive_under.*fieldmat(:,2),[numel(lat_X) 1],@sum);
len_ssh_under_geo = 2*accumarray(posloc,is_under_ssh_var.*fieldmat(:,2),[numel(lat_X) 1],@sum);
len_height_under_geo = 2*accumarray(posloc,is_under_height_var.*fieldmat(:,2),[numel(lat_X) 1],@sum);
len_both_under_geo = 2*accumarray(posloc,is_under_both_var.*fieldmat(:,2),[numel(lat_X) 1],@sum);


% Don't need to multiply by 2 with M
len_M_ssh_under_geo = accumarray(posloc,is_under_ssh_var.*fieldmat(:,2).*Mval,[numel(lat_X) 1],@sum);
len_M_height_under_geo = accumarray(posloc,is_under_height_var.*fieldmat(:,2).*Mval,[numel(lat_X) 1],@sum);
len_M_both_under_geo = accumarray(posloc,is_under_both_var.*fieldmat(:,2).*Mval,[numel(lat_X) 1],@sum);

% 2-D map of fraction by length in each cell
perc_under_geo = len_under_geo ./ len_meas_geo;
perc_under_geo(num_meas_geo==0) = nan;

perc_ssh_under_geo = len_ssh_under_geo ./ len_meas_geo;
perc_ssh_under_geo(num_meas_geo==0) = nan;

perc_height_under_geo = len_height_under_geo ./ len_meas_geo;
perc_height_under_geo(num_meas_geo==0) = nan;

perc_M_ssh_under_geo = len_M_ssh_under_geo ./ len_meas_geo;
perc_M_ssh_under_geo(num_meas_geo==0) = nan;

perc_M_height_under_geo = len_M_height_under_geo ./ len_meas_geo;
perc_M_height_under_geo(num_meas_geo==0) = nan;

naive_perc_under_geo = naive_len_under_geo ./ naive_len_meas_geo;
naive_perc_under_geo(naive_num_meas_geo==0) = nan;

% ssh anomaly
ssh_anom_geo = accumarray(posloc,ssh,[numel(lat_X) 1],@mean);
ssh_std_geo = accumarray(posloc,ssh_moving_std,[numel(lat_X) 1],@mean);

% Height variance 
height_std_geo = accumarray(posloc,height_moving_std,[numel(lat_X) 1],@mean);

% Accumulate wave energy into a single matrix
wave_energy_geo = accumarray(posloc,moving_en,[numel(lat_X) 1],@sum);

%%

disp(['Saving ' outdir]);

save(outdir,'lat_X','lon_X','*_geo','numsegs','numtracks');
