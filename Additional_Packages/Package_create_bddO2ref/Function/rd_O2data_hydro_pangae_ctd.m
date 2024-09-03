function [refID,lon,lat,sta_ctd,juld,presi,tempi,theta0i,psali,sig0i,doxyi,o2ok] = rd_O2data_hydro_pangae_ctd(ficcamp,datatype,presi,lid,scar,addzero )
%
% Fonction qui lit les donnees oxygene calibrees (campagne et bouteille)
% pour calibrer les donnees flotteurs
% Input:
%        ficcamp : nom du fichier HYDRO (avec chemin complet)
%        datatype = type de donnees lues
%                 = 'ctd' ou 'btl'
%
%
% addpath /Users/thierry_reynaud/IFREMER/MATLAB/seawater_330_its90_lpo
%

list_files=dir(strcat(ficcamp,'*.ctd'));

%[num,txt,raw] =xlsread(ficcamp);

tempi=NaN*ones(length(list_files),length(presi));
theta0i=NaN*ones(length(list_files),length(presi));
psali=NaN*ones(length(list_files),length(presi));
sig0i=NaN*ones(length(list_files),length(presi));
doxyi=NaN*ones(length(list_files),length(presi));

%refID=[];

if strcmp(datatype,'ctd') == 1
    % Lecture des donnees campagne CTD
    o2ok = 1;

    %     [refID, lat,lon,juld,platform,temp,theta0,psal,doxy,pres,sig0,sta_ctd] = ...
    %         read_nc_base_hydro(ficcamp,'TEMP','TPOT','PSAL','OXYK','PRES','SIG0','STATION_NUMBER');
    %     temp(temp == -9999) = NaN;
    %     theta0(theta0 == -9999) = NaN;
    %     psal(psal == -9999) = NaN;
    %     doxy(doxy == -9999) = NaN;
    %     pres(pres == -9999) = NaN;
    %     sig0(pres == -9999) = NaN;
    for i=1:length(list_files)
        filename=fullfile(list_files(i).folder,list_files(i).name);
        data=get_pangaea_ctd(filename);
        %refID{i}=sprintf('%s%i_%i',scar,data.nsta,data.nprof);
        refID{i}=sprintf('%s_%i',scar,data.nprof);
        sta_ctd(i)=i;
        juld(i)=data.juld;
        lat(i)=data.lat;
        lon(i)=data.lon;
        pres=data.pres;
        temp=data.temp;
        psal=data.psal;
        doxy=data.doxy;

        pr=0;
        theta0=sw_ptmp(psal,temp,pres,pr);
        sig0=sw_pden(psal,temp,pres,pr)-1000;

        iok = find(isfinite(pres));
        tempi(i,:) = interp1(pres(iok),temp(iok),presi);
        theta0i(i,:) = interp1(pres(iok),theta0(iok),presi);
        psali(i,:) = interp1(pres(iok),psal(iok),presi);
        sig0i(i,:) = interp1(pres(iok),sig0(iok),presi);
        doxyi(i,:) = interp1(pres(iok),doxy(iok),presi);
    end

elseif strcmp(datatype,'btl') == 1
    % Lecture des donnees campagne bouteilles
    o2ok = 0;
    pres = NaN;
    temp = NaN;
    theta0 = NaN;
    psal = NaN;
    sig0 = NaN;
    oxyl = NaN;
    doxy = NaN;
end

return
end


