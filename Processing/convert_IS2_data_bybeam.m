function convert_IS2_data_bybeam(year,month,beamind,filedir,savedir)
%% Code that loads in data from the ATL10 hdf5 files - bulky, yes, but necessary
% to susequently use all Matlab processing. The saved data is contained in
% two arrays which have consistent dimensions that can be segmented by
% common matrix commands.

% Optional input - save_path - a filename where the .mat file will be saved
% will be defaulted to Data.mat if not applied

% Output - a file $savename$.mat, which contains
% fields - a cell array Ntracks x Nbeams x Nfields for segment-based data
% lead_fields - a cell array Ntracks x Nbeams x Nfields for lead-based data
% timer - a character array Ntracks x 28 in size - each row describes the
% data on which the track began.

% field_names - a cell array with the name of each field in fields
% lead_field_names - " for lead_fields

%% Options are listed in this code block

% Directory of the files to consider - will load in all .h5 files in that
% directory

% Which beams to load in
beam_names = {'gt1r','gt1l','gt2r','gt2l','gt3r','gt3l'};

% Data that is read in on a per-segment basis
field_names = ...
    { ... % Get the ice segment information
    'freeboard_beam_segment/height_segments/height_segment_height' ...   
    'freeboard_beam_segment/height_segments/height_segment_length_seg', ...
    'freeboard_beam_segment/height_segments/latitude', ...
    'freeboard_beam_segment/height_segments/longitude', ...
    'freeboard_beam_segment/height_segments/photon_rate', ...
    'freeboard_beam_segment/beam_freeboard/beam_fb_height', ...
    'freeboard_beam_segment/height_segments/height_segment_type', ...
    'freeboard_beam_segment/height_segments/height_segment_ssh_flag', ...
    'freeboard_beam_segment/height_segments/ice_conc',
    };

yrstr = num2str(year); 

mostr=num2str(month);
if month<10, mostr=['0',mostr]; end

% Obtains all .h5 files in that directory
[filedir '*ATL10*' yrstr mostr '*.nc']
files = dir([filedir '*ATL10*' yrstr mostr '*.nc']);

%% Initialize the output fields

fields = cell(length(files),length(field_names));
timer = cell(length(files),1);

%%
for fileind = 1:length(files)
    
    if mod(fileind,100) == 1

    fprintf('File %d of %d \n',fileind,length(files));

    end    
    
    filename = [filedir files(fileind).name];
    
    % Load in start of track
    try 
	timer{fileind} = h5readatt(filename,'/','time_coverage_start');
	catch timerrread
	disp('WE HAVE A GODDAMN ERROR RIGHT HERE');
	filename
	disp('END OF ERROR'); 
	throw(timerrread)
end        

        % For all segment-indexed fields
        for fieldind = 1:length(field_names)
            
            try
                
            fields{fileind,fieldind} = double(h5read(filename,['/' beam_names{beamind} '/' field_names{fieldind}]));
           
            catch errread
                
                fields{fileind,fieldind} = double(zeros(0,1));
                
            end
            
        end
    
end

%%

% Convert the timing cell array to a character array as all inputs have the
% same dimensions
%timer = cell2mat(timer);

if length(files) > 0

save_path = [savedir yrstr mostr '-beam-' beam_names{beamind} '.mat']
save(save_path,'field_names','fields','timer','beam_names','-v7.3');

end

end
