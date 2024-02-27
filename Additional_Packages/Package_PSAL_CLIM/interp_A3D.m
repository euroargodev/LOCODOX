function [tmp]=interp_A3D(lon_float1,lat_float1,depth_float,juld_float1,rep,filename)
% Author Thierry Reynaud
% 2024/02/07

% Generate Armor-3D interpolated files
% Interpolates ARMOD-3D fields on the float time and depth.

% Read times
fname=fullfile(rep,filename);
time_file=double(ncread(fname,'time'));
scale_units = ncreadatt(fname,'time','units');
dateref=sscanf(scale_units,'hours since %f-%f-%f %f:%f:%f UTC');
time_A3D=time_file/24+datenum(dateref');

% Read longitude -180 to 180
lon_A3D=double(ncread(fname,'longitude'));
%Read latitude
lat_A3D=double(ncread(fname,'latitude'));
%Read Depth
z_A3D=double(ncread(fname,'depth'));

if strfind(filename,'PSAL')
    var_A3D=double(ncread(fname,'so'));% lon x lat * depth x time
else
    var_A3D=double(ncread(fname,'to'));
end

tmp = interpn(lon_A3D, lat_A3D,z_A3D,time_A3D,var_A3D,lon_float1,lat_float1,depth_float,juld_float1);
tmp=squeeze(tmp);

return;
end