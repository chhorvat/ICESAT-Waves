clear
% in NH
lat = -90:90; 
lon = -180:180; 

[LAT,LON] = meshgrid(lat,lon); 

load('bybeam/201905-beam-gt1l.mat'); 

%% Do with KNN
lat_X = [LAT(:)];  
lon_X = [LON(:)]; 

M = createns([lat_X,lon_X]); 

lat_Y = cell2mat(fields(:,3)); 
lon_Y = cell2mat(fields(:,4)); 

ID = knnsearch(M,[lat_Y lon_Y],'K',1); 

%% Do in old way 

F = scatteredInterpolant(lat_X, lon_X(:),(1:numel(lat_X)).', 'nearest');

% Linear index for every one.
% posloc = zeros(size(fieldmat,1),1);

% This can take quite a long time to identify the location of every segment
% within the CESM grid.

numsegs = length(lat_Y)

fprintf('Have %d total segments to analyse \n',numsegs)

% The function (F) evaluation can die if too many inputs are sent.
% So we send one million at a time. 
% This is way, way faster than doing each function call individually! 

num_mil = ceil(numsegs/1e6);

posloc = zeros(numsegs,1); 

for i = 1:num_mil
    
    index = (i-1)*1e6+1:i*1e6;
    index(index > numsegs) = [];
    
    posloc(index) = F(lat_Y(index),lon_Y(index));
    
    fprintf('%d million out of %2.1f million \n',i,numsegs/1e6)
    
    
end

%%
earthellipsoid = referenceSphere('earth','km');


latbounds = LAT(1,:)'; 
lonbounds = LON(:,1); 
latbounds = [latbounds; latbounds(1);]; 
lonbounds = [lonbounds; lonbounds(1);]; 

%%

AREA = 0*LAT;

for i = 1:size(LON,1)
    for j = 1:size(LAT,2)
        
        AREA(i,j) = areaquad(latbounds(j),lonbounds(i),latbounds(j+1),lonbounds(i+1),earthellipsoid);
        
    end
end


grid_area = double(grid_area);