function [refID,lon,lat,sta_ctd,juld,presi,tempi,theta0i,psali,sig0i,doxyi,o2ok] = rd_O2data_hydro_pangae(ficcamp,ctd_sta,datatype,presi,lid,scar,addzero )
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
% % 1900943
% % Stations: CTD128 and CTD129
% ctd_sta=128:129;
% 
% datatype='ctd';

% lid=1 ==> built refid
% 
% ficcamp='/Users/thierry_reynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_DATA/Convert_CTD/data/PANGAE_POS348_ctd_2007.tab.xlsx';
[num,txt,raw] =xlsread(ficcamp);

tempi=NaN*ones(length(ctd_sta),length(presi));
theta0i=NaN*ones(length(ctd_sta),length(presi));
psali=NaN*ones(length(ctd_sta),length(presi));
sig0i=NaN*ones(length(ctd_sta),length(presi));
doxyi=NaN*ones(length(ctd_sta),length(presi));

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
    for i=1:length(ctd_sta)
        sta_ctd(i)=ctd_sta(i);
        scar0=scar;
        
        if sta_ctd(i)<10&addzero
           scar0=[scar0,'0'];
        end
        CTD_NUM=[scar0,num2str(ctd_sta(i))];
        idx=find(strcmpi(raw(:,2), CTD_NUM));
        
        dateC=cell2mat(raw(idx(1),3));
        juld(i)=datenum(dateC,'yyyy-mm-ddTHH:MM');
        
        latC=cell2mat(raw(idx(1),4));
        lat(i)=double(latC);
        
        lonC=cell2mat(raw(idx(1),5));
        lon(i)=double(lonC);
        
        tmp=cell2mat(raw(idx(1),1));
        jdx=strfind(tmp,'_');
        if length(jdx>0)
        jdx=jdx(1);
        else
            jdx=length(tmp)+1;
        end
        
        refID{i}=[strcat(tmp(1:jdx-1),'_',num2str(ctd_sta(i),'%i'))];
        
        pres=cell2mat(raw(idx,6));
        temp=cell2mat(raw(idx,8));
        psal=cell2mat(raw(idx,9));
        doxy=cell2mat(raw(idx,10));
        
        pr=0;
        theta0=sw_ptmp(psal,temp,pres,pr);
        sig0=sw_pden(psal,temp,pres,pr)-1000;
        
        %sta_ctd = sta_ctd';
        
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




