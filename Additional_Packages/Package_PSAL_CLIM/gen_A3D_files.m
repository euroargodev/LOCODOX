function [out]=gen_A3D_files(juld_float,lat_float,lon_float)
% Author Thierry Reynaud
% 2024/02/07

% Generate Armor-3D interpolated files
% Interpolates ARMOD-3D fields on float a specific float position.


%lon_float=-24.1947183333333;
%lat_float=52.4994333333333;
%juld_float=7.372402719097222e+05;


% Set directories and output filesnames
tic;
file_PSAL='dataset_subset_PSAL.nc';
file_TEMP='dataset_subset_TEMP.nc';

out.dir=pwd;
out.file_PSAL=file_PSAL;
out.file_TEMP=file_TEMP;

if exist(file_PSAL,'file')
    delete(file_PSAL);
end

if exist(file_TEMP,'file')
    delete(file_TEMP);
end

%Access to copernicusmarine services.
% Read login and password
filename_access='~/.coperniscus-access';
fid = fopen(filename_access,'r');
access.login=fscanf(fid,'%s/n');
access.pwd=fscanf(fid,'%s/n');
fclose(fid);

wline1='#!/bin/csh -f';
wline2='source ~/.cshrc';
wline3='conda activate cmt_1.0';
wline4='copernicusmarine --version';


filename='copernicus.scr';
if exist(filename,'file')
    delete(filename);
end

% Build the CLI subset lines
% CLI subset do not interpolate times
% get 2 profiles 1 before and 1 after
delta=21;% Days
juld_min=juld_float-delta;
cjuld_min=datestr(juld_min,'yyyy-mm-ddTHH:MM:SS');
juld_max=juld_float+delta;
cjuld_max=datestr(juld_max,'yyyy-mm-ddTHH:MM:SS');

% La re-analyse:
% dataset_id="dataset-armor-3d-rep-weekly", # 1993-01-06 / 2022-12-28
% Le real-time:
% dataset_id="dataset-armor-3d-nrt-weekly", # 2019-01-02 / 2024-01-24
juld_copernicus_rep=datenum(2022,12,28,0,0,0);

if juld_max<=juld_copernicus_rep
    dataset_id="dataset-armor-3d-rep-weekly";
else 
    dataset_id="dataset-armor-3d-nrt-weekly"
end
%wline_part1_PSAL='copernicusmarine subset -i dataset-armor-3d-rep-weekly -v so';
%wline_part1_TEMP='copernicusmarine subset -i dataset-armor-3d-rep-weekly -v to';

wline_part1_PSAL='copernicusmarine subset -i dataset_id -v so';
wline_part1_PSAL=strrep(wline_part1_PSAL,'dataset_id',dataset_id);

wline_part1_TEMP='copernicusmarine subset -i dataset_id -v to';
wline_part1_TEMP=strrep(wline_part1_TEMP,'dataset_id',dataset_id);


wline_part2='-t 2022-01-01T00:00:00 -T 2022-01-01T00:00:00';% Original
wline_part2='-t cjuld_min -T cjuld_max';% Modified
wline_part2=strrep(wline_part2,'cjuld_min',cjuld_min);
wline_part2=strrep(wline_part2,'cjuld_max',cjuld_max);


%wline_part3='-x -6.17 -X -6.17 -y 35.75 -Y 35.75 -z 0.0 -Z 4000.0';Original
wline_part3='-x lonmin -X lonmax -y latmin -Y latmax -z 0.0 -Z 4000.0';
delta_X=0.25;
tmp=sprintf('%13.8f',lon_float-delta_X);
tmp=deblank(tmp);
tmp=deblank(fliplr(tmp));
tmp=fliplr(tmp);
wline_part3=strrep(wline_part3,'lonmin',tmp);

tmp=sprintf('%13.8f',lon_float+delta_X);
tmp=deblank(tmp);
tmp=deblank(fliplr(tmp));
tmp=fliplr(tmp);
wline_part3=strrep(wline_part3,'lonmax',tmp);

%delta_Y=0.5*cosd(lat_float);
delta_Y=0.5;
tmp=sprintf('%13.8f',lat_float-delta_Y);
tmp=deblank(tmp);
tmp=deblank(fliplr(tmp));
tmp=fliplr(tmp);
wline_part3=strrep(wline_part3,'latmin',tmp);

tmp=sprintf('%13.8f',lat_float+delta_Y);
tmp=deblank(tmp);
tmp=deblank(fliplr(tmp));
tmp=fliplr(tmp);
wline_part3=strrep(wline_part3,'latmax',tmp);

wline_part4_PSAL='-f dataset_subset_PSAL.nc --disable-progress-bar --overwrite --force-download';
wline_part4_TEMP='-f dataset_subset_TEMP.nc --disable-progress-bar --overwrite --force-download';
wline_part5='--username LOGIN --password PWD';
wline_part5=strrep(wline_part5,'LOGIN',access.login);
wline_part5=strrep(wline_part5,'PWD',access.pwd);

fid = fopen(filename,'w');

fprintf(fid,'%s\n',wline1);
fprintf(fid,'%s\n',wline2);
fprintf(fid,'%s\n',wline3);
fprintf(fid,'%s\n',wline4);
fprintf(fid,'%s %s %s %s %s\n',wline_part1_PSAL,wline_part2,wline_part3,wline_part4_PSAL,wline_part5);
fprintf(fid,'%s %s %s %s %s\n',wline_part1_TEMP,wline_part2,wline_part3,wline_part4_TEMP,wline_part5);

fclose(fid);

line_chmod=['chmod u+x ',filename];
system(line_chmod);
system([filename,'>&output']);

out.error_file='output';
toc;
return
end