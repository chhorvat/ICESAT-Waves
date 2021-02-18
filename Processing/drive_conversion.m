data_loc = '/gpfs/data/epscor/chorvat/IS2/All_Tracks_NetCDF/'

filedirs = {[data_loc '/NH/'],[data_loc '/SH/']};

savedirs = {'/gpfs/data/epscor/chorvat/IS2/Data_By_Beam/NH/','/gpfs/data/epscor/chorvat/IS2/Data_By_Beam/SH/'};

parpool()

for i = 2
    for yr = 2018:2020
        parfor mo = 1:12
            for beamind = 1:6

            convert_IS2_data_bybeam(yr,mo,beamind,filedirs{i},savedirs{i});

            end

        end

    end

end




