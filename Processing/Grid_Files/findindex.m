clear
% in NH
fid = fopen('Grid_files/psn25lats_v3.dat','r'); 
lat_NH = fread(fid,'integer*4')/1e5; 
lat_NH = reshape(lat_NH,[304 448]); 

fid = fopen('Grid_files/psn25lons_v3.dat','r'); 
lon_NH = fread(fid,'integer*4')/1e5; 
lon_NH = reshape(lon_NH,[304 448]); 

fid = fopen('Grid_files/psn25area_v3.dat','r'); 
area_NH = fread(fid,'integer*4')/1e5; 
area_NH = reshape(area_NH,[304 448]); 

fid = fopen('Grid_files/pss25lats_v3.dat','r'); 
lat_SH = fread(fid,'integer*4')/1e5; 
lat_SH = reshape(lat_SH,[316 332]); 

fid = fopen('Grid_files/pss25lons_v3.dat','r'); 
lon_SH = fread(fid,'integer*4')/1e5; 
lon_SH = reshape(lon_SH,[316 332]); 

fid = fopen('Grid_files/pss25area_v3.dat','r'); 
area_SH = fread(fid,'integer*4')/1e5; 
area_SH = reshape(area_SH,[316 332]); 

lat_X = [lat_NH(:); lat_SH(:)];  
lon_X = [lon_NH(:); lon_SH(:)]; 

save('grid-25km.mat','lat_X','lon_X','*_SH','*_NH')