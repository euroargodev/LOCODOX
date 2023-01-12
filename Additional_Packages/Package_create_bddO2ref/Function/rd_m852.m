function [lat,lon,juld,pres,deph,temp,psal,doxy,sig0] = rd_m852(filename)
%
%function [lat,lon,date,pres,deph,temp,psal,doxy,sig0] = rd_m852(inpath,fname)
% Lecture des fichiers de la campagne m852 de Johannes Karstensen
%

fid=fopen([filename],'r');
for i=1:4
    vline=fgetl(fid);
end
vline=fgetl(fid);
ieg=findstr(vline,'=');
lat=str2num(vline(ieg+1:end));
vline=fgetl(fid);
ieg=findstr(vline,'=');
lon=str2num(vline(ieg+1:end));
vline=fgetl(fid);
ieg=findstr(vline,'=');
vdate=(vline(ieg+1:end));
vline=fgetl(fid);
ieg=findstr(vline,'=');
vtime=(vline(ieg+1:end));

juld=datenum(str2double(vdate(1:5)),str2double(vdate(7:8)),str2double(vdate(10:11)));

for i=9:48
    vline=fgetl(fid);
end
ind=0;
while ischar(vline) == 1
    ind=ind+1;
    vline=str2num(fgetl(fid));
    pres(ind)=vline(1);
    deph(ind)=vline(2);
    temp(ind)=vline(3);
    psal(ind)=vline(4);
    doxy(ind)=vline(5);
    sig0(ind)=vline(7);
    vline=fgetl(fid);
end

fclose(fid)



