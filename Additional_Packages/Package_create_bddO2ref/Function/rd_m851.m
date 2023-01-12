function [lat,lon,juld,pres,temp,tpot,psal,doxy,sig0] = rd_m851(filename)
%
%function [lat,lon,date,pres,deph,temp,psal,doxy,sig0] = rd_m852(inpath,fname)
% Lecture des fichiers de la campagne m852 de Johannes Karstensen
%

fid=fopen([filename],'r');
for i=1:6
    vline=fgetl(fid);
end
vline=fgetl(fid);
ieg=findstr(vline,'=');
vdate=(vline(ieg+1:end));

juld=datenum(str2double(vdate(1:5)),str2double(vdate(7:8)),str2double(vdate(10:11)));


vline=fgetl(fid);
ieg=findstr(vline,'=');
vtime=(vline(ieg+1:end));

vline=fgetl(fid);
ieg=findstr(vline,'=');
lat=str2num(vline(ieg+1:ieg+3))+str2num(vline(ieg+5:ieg+9))/60;

vline=fgetl(fid);
ieg=findstr(vline,'=');
lon=-(str2num(vline(ieg+1:ieg+5))+str2num(vline(ieg+6:ieg+10))/60);



for i=9:48
    vline=fgetl(fid);
end
vlin=fgetl(fid);
ind=0;
while ischar(vlin) == 1
    ind=ind+1;
    vline=str2num(vlin);
    pres(ind)=vline(1);
    temp(ind)=vline(2);
    tpot(ind)=vline(3);
    psal(ind)=vline(4);
    sig0(ind)=vline(5);
    oxyl(ind)=vline(6);
    doxy(ind)=convert_oxygen(oxyl(ind),'mL/L','mumol/kg',sig0(ind));
    vlin=fgetl(fid);
end

fclose(fid)



