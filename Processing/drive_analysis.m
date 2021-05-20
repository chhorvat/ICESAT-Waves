%% Analyse_frac_wave_affected
% This code takes the collected track data, the subsetted data, bins it
% into lat-lon bins, and then produces maps of the fraction of those bins
% that are affected by waves.

% Input require save_path - which is the directory of the .mat file created
% by previous code.

clear

file_heading = './';

hemi_dir = {'NH','SH',};

DO_WAVES = 1;
DO_FSD = 0;
DO_REPLACE = 0;

% try
%  parpool()
% catch err
%   err
% end
%%
gridname = '100km';


addpath('../Processing/');

for i = 1:2 % length(hemi_dir)
    
    file_dir = [file_heading hemi_dir{i}];
    
    files = dir([file_dir '/*.mat']);
    
    % Overall place where we will save each file
    save_dir = [file_heading 'Processed/' hemi_dir{i} '-' gridname '/'];
    
    mkdir(save_dir);
    mkdir([save_dir 'WAVES/']);
    mkdir([save_dir 'FSD/']);
    
    % par
    for filename = 1:length(files)% :-1:1
        
        disp('-----------------------------------------');
        disp(['File at ' file_dir files(filename).name(1:end-4)]);
        
        % Specific file save name
        save_loc_waves = [save_dir 'WAVES/' files(filename).name];
        save_loc_fsd = [save_dir 'FSD/' files(filename).name];
        
        
        if (exist([save_loc_waves]) == 2) & (~DO_REPLACE)
            
            make_waves = 0;
            disp('Already There');
            
        else
            
            make_waves = 1;
            disp(['Saving wave files to ' save_loc_waves]);
            
        end
       
    
    if DO_FSD
        
        if (exist([save_loc_fsd]) == 2) & (~DO_REPLACE)
            
            make_fsd = 0;
            disp('Already There');
            
        else
            
            make_fsd = 1;
            
            disp(['Saving FSD files to ' save_loc_fsd]);
            
        end
        
    end

if make_fsd + make_waves > 0
    
    analyse_waves_and_FSD([file_dir '/' files(filename).name], ...
        save_loc_waves,save_loc_fsd,gridname,make_waves,make_fsd);

end    

    end
    
end
