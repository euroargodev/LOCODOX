% =========================================================================
% Create reference database (structure REF) for O2 correction with LOCODOX
% Create the structure .mat used as input in the software
%
% adaptation by M.GALLIAN (ALTRAN) in february 2019 of the function gen_bdd_O2ref_all.m
% (Virginie Thierry and Emilie Brion (ALTRAN))
%
% Explaination :
% The user should change the "path_file" and set the corresponding function to decode the netcdf file ("name_func").
% Run this code as many netcdf files you have to add to the database, and select "extend" after the first run.
% Be carefull to not write 2 times the same file in the reference matrix, it will induce an error in LOCODOX

% =========================================================================

clear all;
close all;

% =========================================================================
%% PARAMETERS
% =========================================================================
% data type :
% - Hydro LOPS : name_func=rd_O2data_hydro_LOPS
% - Strass : name_func=rd_O2data_hydro_strass
% - Ovide2011me : name_func= rd_O2data_hydro_ovid11
% - Ovide2011di : name_func= rd_O2data_hydro_ovid11
    
%Define the function to use to decode the netcdf file and its directory
%The function has to return : refId,lon,lat,juld,pres,temp,psal,sig0,doxy
%Path of the netcdf file to decode or name (=input 1 of the function)

path_func='Function/';

rep_data_lpo='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/Convert_CTD/data_lpo/';

for ifile=15:15
%for ifile=10:10    
    
%     if i==1; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/ovid10_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
%     if i==2; path_file='//Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/cat12_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
%     if i==3; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/geov_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
%     if i==4; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/OVIDE/bo16_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
%     if i==5; path_file='/Volumes/qlpo5/HYDROCEAN/OVIDE18_TEMP/MLT/ov18_d_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
%     if i==6; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/RREX/RREX15/rr15_PRES.nc';name_func='rd_O2data_hydro_LOPS';end
%     if i==7; path_file='/Volumes/qlpo5/HYDROCEAN/MLT_NC/LPO/RREX/RREX17/rr17_PRES.nc';name_func='rd_O2data_hydro_LOPS';end

    if ifile==1; path_file=strcat(rep_data_lpo,'ovid10_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    if ifile==2; path_file=strcat(rep_data_lpo,'cat12_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    if ifile==3; path_file=strcat(rep_data_lpo,'geov_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    if ifile==4; path_file=strcat(rep_data_lpo,'bo16_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    if ifile==5; path_file=strcat(rep_data_lpo,'ov18_d_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    if ifile==6; path_file=strcat(rep_data_lpo,'rr15_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    if ifile==7; path_file=strcat(rep_data_lpo,'rr17_PRES.nc');name_func='rd_O2data_hydro_LOPS';end


    if ifile==8; path_file='strass';name_func='rd_O2data_hydro_strass';end
    if ifile==9; path_file='ov11me';name_func= 'rd_O2data_hydro_ovid11';end
    if ifile==10; path_file='ov11di';name_func= 'rd_O2data_hydro_ovid11';end
    if ifile==11; path_file='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/Convert_CTD/data/PANGAE_POS348_ctd_2007.tab.xlsx';name_func= 'rd_O2data_hydro_pangae';end
    if ifile==12; path_file='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/Convert_CTD/data/MSM_08_1_phys_oce.tab_tr.xlsx';name_func= 'rd_O2data_hydro_pangae';end
    if ifile==13; path_file='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/Convert_CTD/data/ATA03_phys_oce.tr.xlsx';name_func= 'rd_O2data_hydro_pangae';end
    if ifile==14; path_file=strcat(rep_data_lpo,'ovid21_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    if ifile==15; path_file=strcat(rep_data_lpo,'bo23_PRES.nc');name_func='rd_O2data_hydro_LOPS';end
    display(path_file);
    
    %Path of the .mat structure to update, or create
    %BDD_O2_REF= 'bddo2ref_avecov18_temp.mat';
    %BDD_O2_REF= 'bddo2ref_ov18_temp.mat';
    %BDD_O2_REF= 'bddo2ref_vracape.mat';
    BDD_O2_REF= 'bddo2ref_all_TR.mat';

    % =========================================================================
    %% Initialisation
    % =========================================================================
    
    %Add path of the function
    addpath(path_func)
    
    %Initialise the REF structure
    REF.id   =   [];
    presi = 0:10:4000;
    REF.lon  =   [];
    REF.lat  =   [];
    REF.juld =   [];
    REF.temp =   [];
    REF.psal =   [];
    REF.sig0 =   [];
    REF.doxy =   [];
    REF.pres =   [];
    REF.psat =   [];
    
    % =========================================================================
    %% Load reference matrix or overwrite
    % =========================================================================
    
    if exist(BDD_O2_REF,'file')
        answer = questdlg('Reference database exists ! Do you want to overwrite or extend it ?','Refrence database','Overwrite','Extend','Extend');
        if strcmp(answer,'Overwrite')
            if exist(BDD_O2_REF,'file')
                eval(['!\rm ' BDD_O2_REF]);
            end
            REF.id   =   [];
            presi = 0:10:4000;
            REF.lon  =   [];
            REF.lat  =   [];
            REF.juld =   [];
            REF.temp =   [];
            REF.psal =   [];
            REF.sig0 =   [];
            REF.doxy =   [];
            REF.pres =   [];
            REF.psat =   [];
        else
            load(BDD_O2_REF) ;
        end
    end
    
    
    % =========================================================================
    %% Create reference database
    % =========================================================================
    
    %Extract data
    if ifile==11
        addpath /Users/treynaud/IFREMER/MATLAB/seawater_330_its90_lpo
        ctd_sta=128:129; % Float 1900943
        scar='CTD';%Search Character
        lid=1;
        addzero=1;
        [refId,lon,lat,~,juld,~,tempi,~,psali,sig0i,doxyi,~] = feval(name_func,path_file,ctd_sta,'ctd',presi,lid,scar,addzero);
    elseif ifile==13
        addpath /Users/treynaud/IFREMER/MATLAB/seawater_330_its90_lpo
        ctd_sta=[7,35,40]; % Float 6900627 6900628 6900630
        lid=0;
        scar='ATA3_';
        addzero=1;
        [refId,lon,lat,~,juld,~,tempi,~,psali,sig0i,doxyi,~] = feval(name_func,path_file,ctd_sta,'ctd',presi,lid,scar,addzero);
    elseif ifile==12
        addpath /Users/treynaud/IFREMER/MATLAB/seawater_330_its90_lpo
        ctd_sta=[7,9,20,24,25,26,39,40]; % Float 6900524 6900525 6900629
        scar='1_CTD-RO';
        lid=1;
        addzero=0;
        [refId,lon,lat,~,juld,~,tempi,~,psali,sig0i,doxyi,~] = feval(name_func,path_file,ctd_sta,'ctd',presi,lid,scar,addzero);
    else
        [refId,lon,lat,~,juld,pres,temp,~,psal,sig0,doxy,~] = feval(name_func,path_file,'ctd');
    end
    
    %Interpolate data on pressure
    if ~isempty(lon) && ~strcmp(name_func,'rd_O2data_hydro_LOPS')&& ~strcmp(name_func,'rd_O2data_hydro_pangae')
        for i = 1:length(lon)
            iok = find(isfinite(pres(i,:)));
            pres=double(pres);
            temp=double(temp);
            psal=double(psal);
            sig0=double(sig0);
            doxy=double(doxy);
            tempi(i,:) = interp1(pres(i,iok),temp(i,iok),presi);
            psali(i,:) = interp1(pres(i,iok),psal(i,iok),presi);
            sig0i(i,:) = interp1(pres(i,iok),sig0(i,iok),presi);
            doxyi(i,:) = interp1(pres(i,iok),doxy(i,iok),presi);
        end
    elseif ~isempty(lon) && strcmp(name_func,'rd_O2data_hydro_LOPS')
        tempi=temp;
        psali=psal;
        sig0i=sig0;
        doxyi=doxy;
    end
    
    jdate=find(juld>=datenum(1950,1,1));
    if ~isempty(jdate)
        juld(jdate)=single(juld(jdate)-datenum(1950,1,1));
    end
    juld=single(juld);
    tempi=single(tempi);
    psali=single(psali);
    sig0i=single(sig0i);
    doxyi=single(doxyi);
    %Add reference data to the REF structure
    if ~isempty(lon)
        tmp=[REF.lon;reshape(lon,length(lon),1)];
        REF.lon = tmp;clear tmp;

        tmp=[REF.lat;reshape(lat,length(lat),1)];
        REF.lat = tmp;clear tmp;

        tmp=[REF.juld;reshape(juld,length(juld),1)];
        REF.juld = tmp;clear tmp;

        REF.temp = [REF.temp;tempi];
        REF.psal = [REF.psal;psali];
        REF.sig0 = [REF.sig0;sig0i];
        REF.doxy = [REF.doxy;doxyi];

        tmp=[REF.id;cellstr(reshape(refId,length(refId),1))];
        REF.id = tmp;clear tmp;

        tmp=ones(1,length(REF.lon)).*presi';       
        REF.pres = tmp';clear tmp;

        tmpsat = sw_satO2(double(REF.psal),double(REF.temp));
        sat = DOXY_convert(tmpsat,'mL/L','mumol/kg',double(REF.sig0));
        REF.psat = single(100*double(REF.doxy)./sat);
        
        clear lon lat juld temp psal sig0 doxy doxyi tempi psali sig0i pres psat 
        % Save the data
        save(BDD_O2_REF,'REF');
%         display('<===============>')
%         ifile
%         REF

        %refId_old=refId;
    end
    
    
end
