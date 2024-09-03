function data=get_pangaea_ctd(filename)

fid=fopen(filename,'rt');
iread=1;
while iread
    line=fgetl(fid);
    if ~isempty(strfind(line,'Columns'))
        iread=0;%Start Reading data
    end
    if ~isempty(strfind(line,'Station'))
        idx=strfind(line,'=');
        cval=line(idx+1:end);
        data.nsta=str2num(cval);
    end
    if ~isempty(strfind(line,'Profile'))
        idx=strfind(line,'=');
        cval=line(idx+1:end);
        data.nprof=str2num(cval);
    end
    if ~isempty(strfind(line,'Date'))
        if isempty(strfind(line,'FileDat'))
            idx=strfind(line,'=');
            cdate=line(idx+1:end);
        end
    end
    if ~isempty(strfind(line,'Time'))
        idx=strfind(line,'=');
        ctime=line(idx+1:end);
    end
    if ~isempty(strfind(line,'Latitude'))
        idx=strfind(line,'=');
        cval=line(idx+1:end);
        val=sscanf(cval(1:end-1),'%f %f');
        lat=val(1)+val(2)/60;
        if ~isempty(strfind(cval,'S'))
            lat=-lat;
        end
        data.lat=lat; clear lat;

    end
    if ~isempty(strfind(line,'Longitude'))
        idx=strfind(line,'=');
        cval=line(idx+1:end);
        val=sscanf(cval(1:end-1),'%f %f');
        lon=val(1)+val(2)/60;
        if ~isempty(strfind(cval,'W'))
            lon=-lon;
        end
        data.lon=lon;clear lon;
    end
end
data.juld=datenum([cdate ctime],'yyyy/mm/dd HH:MM');
tmp = textscan(fid,'%f %f %f %f %f %f','Delimiter','\n');
fclose(fid);


data.pres=tmp{:,1};%Pressure
data.temp=tmp{:,2};%Temperature
data.psal=tmp{:,4};%Salinity
data.doxy=tmp{:,6};%O2



return
end
