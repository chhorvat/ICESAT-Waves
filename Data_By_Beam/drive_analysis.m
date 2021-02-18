%% Analyse_frac_wave_affected
% This code takes the collected track data, the subsetted data, bins it
% into lat-lon bins, and then produces maps of the fraction of those bins
% that are affected by waves.

% Input require save_path - which is the directory of the .mat file created
% by previous code.

clear

file_heading = './';

hemi_dir = {'NH','SH',};

DO_WAVES = 0;
DO_FSD = 1;
DO_REPLACE = 0;

% try
%  parpool()
% catch err
%   err
% end
%%
gridname = '200km';


addpath('../Processing/');

for i = 1:1 % length(hemi_dir)
    
    file_dir = [file_heading hemi_dir{i}];
    
    files = dir([file_dir '/*.mat']);
    
    % Overall place where we will save each file
    save_dir = [file_heading 'Processed/' hemi_dir{i} '-' gridname '/'];
    
    mkdir(save_dir);
    mkdir([save_dir 'WAVES/']); 
    mkdir([save_dir 'FSD/']); 
    
    % par
    for filename = 1:length(files)
        
        disp('-----------------------------------------');
        disp(['File at ' file_dir files(filename).name]);
        
        % Specific file save name
        save_loc_waves = [save_dir 'WAVES/' files(filename).name];
        save_loc_fsd = [save_dir 'FSD/' files(filename).name];
        
        if DO_WAVES
            
            if (exist([save_loc_waves]) == 2) & (~DO_REPLACE)
                
                disp('Already There');
                
            else
                
                disp(['Saving to ' save_loc_waves]);
                analyse_wave_frac_affected([file_dir '/' files(filename).name],save_loc_waves,gridname);
                
            end
            
        end
        
        if DO_FSD
            
            if (exist([save_loc_fsd]) == 2) & (~DO_REPLACE)
                
                disp('Already There');
                
            else
                
                disp(['Saving to ' save_loc_fsd]);
                analyse_fsd([file_dir '/' files(filename).name],save_loc_fsd,gridname);
                
            end
            
        end
        
    end
    
end
