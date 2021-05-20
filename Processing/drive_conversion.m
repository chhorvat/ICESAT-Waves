data_loc = '/gpfs/data/epscor/chorvat/IS2/All_Tracks_NetCDF/'

filedirs = {[data_loc '/NH/'],[data_loc '/SH/']};

savedirs = {'/gpfs/data/epscor/chorvat/IS2/Data_By_Beam/NH/','/gpfs/data/epscor/chorvat/IS2/Data_By_Beam/SH/'};

% parpool()

beam_names = {'gt1r','gt1l','gt2r','gt2l','gt3r','gt3l'};

DO_REPLACE = 0; 

for i = 1:2
    for yr = 2018:2021
        for mo = 1:12
            for beamind = 1:6


		yrstr = num2str(yr);

		mostr=num2str(mo);
		if mo<10, mostr=['0',mostr]; end

		save_str = [savedirs{i} yrstr mostr '-beam-' beam_names{beamind} '.mat']

		try

			MF = matfile(save_str); 
 			p = size(MF,'fields')
%		if (exist(save_str) ==2) & ~DO_REPLACE
			disp(['Already exists at ' save_str])
%		else
		catch err

            		convert_IS2_data_bybeam(yr,mo,beamind,filedirs{i},savedirs{i});
		end	
		

            end

        end

    end

end




