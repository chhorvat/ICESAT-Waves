clear
% in NH
fid = fopen('Grid_files/psn25lats_v3.dat','r'); 
lat_NH = fread(fid,'integer*4')/1e5; 
lat_NH = reshape(lat_NH, [304 448]); 
lat_NH = lat_NH(3:4:end,3:4:end); 

fid = fopen('Grid_files/psn25lons_v3.dat','r'); 
lon_NH = fread(fid,'integer*4')/1e5;  
lon_NH = reshape(lon_NH,[304 448]); 
lon_NH = lon_NH(3:4:end,3:4:end); 

fid = fopen('Grid_files/pss25lats_v3.dat','r'); 
lat_SH = fread(fid,'integer*4')/1e5; 
lat_SH = reshape(lat_SH,[316 332]); 
lat_SH = lat_SH(3:4:end,3:4:end); 

fid = fopen('Grid_files/pss25lons_v3.dat','r'); 
lon_SH = fread(fid,'integer*4')/1e5; 
lon_SH = reshape(lon_SH,[316 332]); 
lon_SH = lon_SH(3:4:end,3:4:end); 

%%
% Area a little more delicate
fid = fopen('Grid_files/psn25area_v3.dat','r'); 
area_NH = fread(fid,'integer*4')/1e3; 
area_NH = reshape(area_NH,[304 448]); 
area_NH = reshape(area_NH,[4 304/4 4 448/4]); 
area_NH = squeeze(sum(sum(area_NH,1),3)); 

fid = fopen('Grid_files/pss25area_v3.dat','r'); 
area_SH = fread(fid,'integer*4')/1e3; 
area_SH = reshape(area_SH,[316 332]); 
area_SH = reshape(area_SH,[4 316/4 4 332/4]); 
area_SH = squeeze(sum(sum(area_SH,1),3)); 

%% Do with KNN
lat_X = [lat_NH(:); lat_SH(:)];  
lon_X = [lon_NH(:); lon_SH(:)]; 

KDTree = createns([lat_X,lon_X]); 

save('KDTrees/KDTree_100km.mat','KDTree','lat_X','lon_X','*_SH','*_NH'); 
%%


load('bybeam/201905-beam-gt1l.mat'); 

%% Do with KNN

lat_Y = cell2mat(fields(:,3)); 
lon_Y = cell2mat(fields(:,4)); 

ID = knnsearch(KDTree,[lat_Y lon_Y],'K',1); 

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