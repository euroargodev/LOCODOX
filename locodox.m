
% locodox apply correction on the Argo Oxygen data
%
% SYNTAX
% [] = locodox
%
% DESCRIPTION
% locodox apply correction on the Argo Oxygen data, using the ATLN
% method.
%
% INPUT
%   config_file     (Optional) the full path of the configuration file.
%
% OUTPUT
%   Adjusted file written
%
% CALL
%   DOXY_NCR_WOA, DOXY_argo_read, DOXY_argo_prep, DOXY_MAP, DOXY_PLOT_QC,
%   DOXY_interp_WOA_main, grep -ni doxy_drift *.m, DOXY_corr_compute,
%   DOXY_corr_apply_main,DOXY_interp_REF_main,DOXY_update_fields,
%
% SEE ALSO
%

% HISTORY
%   $created: //2009 $author: Mathieu Le Steun, Thomas Bouinot, LPO
%   $Revision: version $Date: $author:
%       v1.2  /
%       v2 18/11/2015   Emilie Brion, ALTRAN OUEST
%                       ergonomics
%       v3 25/05/2016   Emilie Brion, ALTRAN OUEST
%                       argo 3.1 adjustement, Vertical Sampling Scheme
%                       management
%       v4 26/01/2017   Anne Piron, Emilie Brion, Altran Ouest
%                       argo - ameliorations 2017 phase 1
%       v4 21/02/2019   Marine GALLIAN, Altran Ouest
%                       Take into account when there is no PPOX data, the 
%                       IN AIR correction is not possible
%                       Add pressure effect correction and display
%                       corrected data 
%                       Write configuration and choices of the user in the 
%                       summary file
%      v5 07/02/2020    Access to M_MAP added by T. Reynaud
%      V5 12/03/2020    Bug corrected by Thierry Reynaud for NC/MAT directories for REF case
%                       (dirout)
%      v5 05.04.2020    Modifications by T. Reynaud for plotting and
%                       generating corrected oxy profiles figures without   
%                       a REF profile
%      v6 24.04.2020    DOXY_corr_main.m renamed as locodox.m    
%      v7 19.01.2021    Bug corrected by VT+TR    
%      v8 12.04.2021    TD depth modified by TR    
%      v9 07.07.2021    DRIFT/SLOPE ==> test added for TD Degree 2 (T.Reynaud)
%     v10 08.07.2021    DRIFT2 added for TD Degree 2 (T.Reynaud)
%     v11 11.02.2022    Mutliple REF profiles allowed (T.Reynaud)
%     v12 11.04.2022    Bug corrected for Work.launchdate (T.Reynaud)
%     v13 31.10.2023    Online or Dialog Box questions (T. Reynaud and C.Kermabon)
%     v14 10.11.2023    Hybrid DM PSAL and RT PSAL used
%     v15 22.03.2024    Bug correction line #279 ==> biofiles T. Reynaud
%     v16 26.04.2024    PWLF option added for Time Drift Gain Correction (T. Reynaud)
%     v16 17.08.2024    Save matlab Figures as .fig as well

function [] = locodox(config_prog)

close all
% =========================================================================
%% *Initialisation*
% =========================================================================
goProg = 1;

% Read the configuration file, addpath and load files
if nargin == 0
    config_prog = '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX/locodox_config';
end
% if nargin == 0
%     config_prog = '/Users/vthierry/matlab/GITHUB_LOCODOX/LOCODOX/locodox_config';
% end
[config_dir,config_prog] = fileparts(config_prog);
cd(config_dir);
CONFIG = feval(config_prog);

% Prepare results directories 
CONFIG.saveDataDir=fullfile(CONFIG.resultsDir,'data',filesep);
CONFIG.savePlotDir = fullfile(CONFIG.resultsDir,'plots',filesep);
CONFIG.saveLogDir = fullfile(CONFIG.resultsDir,'logs',filesep);
% CONFIG.NCEPDataDir ==> move in DOXY_config.m by TR 09.04.2020
% CONFIG.NCEPDataDir = fullfile(CONFIG.LocodoxMainDir,'data_input',filesep);

% Initialise Work Structure
Work = load('readme');
Work.R2threshold = CONFIG.R2threshold;
Work.adjusted_error_abs = CONFIG.adjusted_error_abs;
Work.adjusted_error_rel = CONFIG.adjusted_error_rel;

% Units
% -------------------------------------------------
% Units shortname = {'mumol/kg','mumol/L','mL/L'};
% Units longname = {'micromole per kilogram','micromole per liter','milliliter per liter'};
Work.unit = 'mumol/kg';
Work.refUnit = CONFIG.refUnit;

% General
Work.Wdata.DM_psal = CONFIG.DM_psal;
Work.Wdata.DM_temp = CONFIG.DM_temp;
Work.Wdata.DM_pres = CONFIG.DM_pres;

% Added by T. Reynaud 10.11.2023 Hybrid DM+RT SAL data use
if isfield(CONFIG,'RT_psal_cycle')
    Work.Wdata.RT_psal_cycle=CONFIG.RT_psal_cycle;
end    

Work.presEff = CONFIG.presEff;
Work.QC_O = CONFIG.QC_O;
Work.QC_P = CONFIG.QC_P;
Work.QC_T = CONFIG.QC_T;
Work.QC_S = CONFIG.QC_S;
Work.makePlot = 1;
Work.savePlot = CONFIG.savePlot;
Work.savePlotFig = CONFIG.savePlotFig;%Added by T. Reynaud 17.08.2024
Work.formattype = CONFIG.formattype;
Work.drift_fitPolynomialDegree = CONFIG.drift_fitPolynomialDegree;
Work.drift_spec = CONFIG.drift_spec;
if CONFIG.drift_spec == 1
    Work.drift_fittype = CONFIG.drift_fittype;
end
Work.fontsize = CONFIG.fontsize;
Work.resol = sprintf('-r%d',CONFIG.resolution);
Work.history.history_software = CONFIG.history_software;
Work.history.history_reference = CONFIG.history_reference;
Work.history.history_software_release=CONFIG.history_software_release;
Work.isokC=CONFIG.isokC;

% Modified by T. Reynaud 26/04/2024- Piece Wise Linear Fitting for Time
% Drift
Work.drift_PWLF_N=CONFIG.drift_PWLF_N;

initialWork = Work;

% =====================================================================
%% *Read the WOA and convert in the usefull unit*
% =====================================================================
WOA.init = DOXY_WOA_read(CONFIG.WOAfile, Work.unit);
if CONFIG.presEff == 1
    [WOA.init.psatwoa] = WOA.init.psatwoa_preswoa;
    WOA.init = rmfield(WOA.init,'psatwoa_preswoa');
else
    WOA.init = rmfield(WOA.init,'psatwoa_preswoa');
end


% =====================================================================
%% *Doxy correction*
% =====================================================================
while goProg

    % Modofied by TR 12/04/2021
    Work.min_drift_depth_deep=CONFIG.min_drift_depth_deep;
    Work.max_drift_depth_deep=CONFIG.max_drift_depth_deep;
    Work.step_drift_depth_deep=CONFIG.step_drift_depth_deep;

    Work.min_drift_depth_surf=CONFIG.min_drift_depth_surf;
    Work.max_drift_depth_surf=CONFIG.max_drift_depth_surf;
    Work.step_drift_depth_surf=CONFIG.step_drift_depth_surf;


    %Modified by T. Reynaud and C. Kermabon 26/04/2024
    Work.onlineq=CONFIG.onlineq;

    options.WindowStyle = 'normal';
    wmoIn = inputdlg('WMO of the float: ','ARGO Float DOXY correction',1,cellstr(''),options);
    if isempty(wmoIn)
        msgbox({'  Good Bye !  ';''},sprintf('%s',CONFIG.history_software),'custom',imread(CONFIG.logo));
        return
    end
    Work.wmo=str2double(wmoIn);
    Work.replist=fullfile(cell2mat(wmoIn),'/profiles/');
    Work.metafile=fullfile(cell2mat(wmoIn),'/');
    close all

%     a=questdlg('Do you want to use PWLF for Time Drift gain calculation ?',sprintf('%d',Work.wmo),'Yes','No','No');
%     if strcmp(a,'Yes')
%         Work.drift_PWLF=1;
%     else
%         Work.drift_PWLF=0;
%     end

    % =================================================================
    %% Chose the correction option : WOA, REF or INAIR
    % =================================================================
    CONFIG.corrTypeDesc = {'WOA climatology','Reference In-Situ database','In-Air correction'};
    CONFIG.corrTypeShort = {'WOA','REF','INAIR'};
    nbr_tmp = listdlg('Name','METHOD FOR THE DOXY GAIN CALCULATION',...
    'PromptString','Which method do you want to use for your DOXY Gain calculation (slope/offset) ?',...
    'SelectionMode','single',...
    'Listsize',[500,150],...
    'ListString',CONFIG.corrTypeDesc);
    
    if isempty(nbr_tmp) 
        Work.whichCorr = questdlg({'No correction method has been choosen.';...
            'Do you want to restart ?'},...
            sprintf('%s: Cancel method correction choice',CONFIG.history_software'),...
            'YES','NO','NO');
        if strcmp(Work.whichCorr,'YES')
            continue
        else
            msgbox({'  Good Bye !  ';''},sprintf('%s',CONFIG.history_software),'custom',imread(CONFIG.logo));
            return
        end        
    end
    Work.whichCorr = CONFIG.corrTypeShort{nbr_tmp};
    clear nbr_tmp
    
    
    fprintf('INFO >>>>>    Apply correction from %s\n',Work.whichCorr);
    
    % Create directory for plots and logs
    if ~exist(fullfile(CONFIG.savePlotDir,Work.whichCorr,num2str(Work.wmo)),'dir')
        mkdir(fullfile(CONFIG.savePlotDir,Work.whichCorr,num2str(Work.wmo)));
    end
    Work.dirPlot = fullfile(CONFIG.savePlotDir,Work.whichCorr,num2str(Work.wmo));
    
    if ~exist(fullfile(CONFIG.saveLogDir,Work.whichCorr,num2str(Work.wmo)),'dir')
        mkdir(fullfile(CONFIG.saveLogDir,Work.whichCorr,num2str(Work.wmo)));
    end
    Work.dirLog = fullfile(CONFIG.saveLogDir,Work.whichCorr,num2str(Work.wmo));
    
    %% --------   SUMMARY FILE   --------------
    % CREATE THE DOC SUMMARIZING ALL THE OPTION CONSIDERED BY THE USER
    summary_file=fullfile(Work.dirLog,sprintf('Summary_correction_options_%d_%s.txt',Work.wmo,datestr(now,'yyyymmdd_HHMM')));
    log_summ= fopen(summary_file,'w');
    DOXY_print_config_in_logfile(log_summ,CONFIG,Work.wmo);
    Work.log_summ= log_summ;

    fprintf(log_summ,'================================================================================\n');
    fprintf(log_summ,'---------------------    PART 2 : CORRECTIONS APPLIED     -----------------------\n');
    fprintf(log_summ,'================================================================================\n');
    fprintf(log_summ,'\n');
    fprintf(log_summ,sprintf('The float correction is based on the %s method\n',Work.whichCorr));
    fprintf(log_summ,'\n');
    fprintf(log_summ,'----------------------------------------------------\n');
    fprintf(log_summ,'\n');
    
    % =====================================================================
    %% *Read the Reference In-Situ data if the correction used is REF or
    %% if the Reference In-situ data exist
    % =====================================================================
    if strcmp(Work.whichCorr,'REF')
        load(CONFIG.bddFile);
        ff=fopen(CONFIG.RefArgoFile);
        linewmo=textscan(ff,'%s %s %s %s');  
        ii = find(strcmp(linewmo{1},wmoIn{1}));
        if isempty(ii)
            warndlg({sprintf('The wmo %d does not exist in the reference list %s.',...
                str2double(wmoIn{1}),CONFIG.RefArgoFile); '';...
                '=> Choose another one or use WOA or IN AIR correction'},...
                sprintf('%s: Uncorrect WMO',CONFIG.history_software'));
            continue
        else
            %Creation de la structure REF_ARGO
            %Modified TR 01/02/2022 for multiple reference profiles
            for k=1:size(ii,1)
                REF_ARGO.wmo(k)   = str2double(cell2mat((linewmo{1}(ii(k)))));
                REF_ARGO.cycle(k) = str2double(cell2mat(linewmo{2}(ii(k))));
                REF_ARGO.refId{k} = [cell2mat(linewmo{3}(ii(k))) '_' cell2mat(linewmo{4}(ii(k)))];
            end
        end
    elseif exist(CONFIG.RefArgoFile,'file')==2
        load(CONFIG.bddFile);
        ff=fopen(CONFIG.RefArgoFile);
        linewmo=textscan(ff,'%s %s %s %s');
        ii = find(strcmp(linewmo{1},wmoIn{1}));
        if ~isempty(ii)
            %Modified TR 25/02/2022 for multiple reference profiles
            for k=1:size(ii,1)
                REF_ARGO.wmo(k)   = str2double(cell2mat((linewmo{1}(ii(k)))));
                REF_ARGO.cycle(k) = str2double(cell2mat(linewmo{2}(ii(k))));
                REF_ARGO.refId{k} = [cell2mat(linewmo{3}(ii(k))) '_' cell2mat(linewmo{4}(ii(k)))];
            end
            CONFIG.REF_ARGO=REF_ARGO;
            %--MG ajouter recuperation profil de reference
        end
    end
    
    
    % =================================================================
    %% Open and prepare the argo data
    % =================================================================
    % ******************************************************
    % check if files exist. If not, warn, and reload LOCODOX to the wmo
    % asking for.
    % ******************************************************
    biofiles = dir(fullfile(CONFIG.DataDir,Work.replist,'B*nc'));
    core_files = dir(fullfile(CONFIG.DataDir,Work.replist,'/*nc'));
    core_files = {core_files(:).name}';
    iscore = ~cellfun('isempty',regexp(core_files,'^R[0-9]*')) | ~cellfun('isempty',regexp(core_files,'^D[0-9]*'));
    core_files = core_files(iscore);
    
    if isempty(core_files) && isempty(biofiles) % Corrected by TR: 22.03.2024
    %if isempty(core_files)
        warndlg(sprintf('No Bio nor Core file in the directory : %s\n \n=> Check CONFIG.DataDir \n',fullfile(CONFIG.DataDir,Work.replist)),'Warning');
        continue
    elseif isempty(biofiles)
        warndlg(sprintf('No Bio file in the directory : %s\n \n=> Check CONFIG.DataDir \n',fullfile(CONFIG.DataDir,Work.replist)),'Warning');
        continue
    elseif isempty(core_files)
        warndlg(sprintf('No Core file in the directory : %s\n \n=> Check CONFIG.DataDir \n',fullfile(CONFIG.DataDir,Work.replist)),'Warning');
        continue
    end
    
    % *****************************************************************
    % Read profile data : if no profil data, warn and stop
    % *****************************************************************
    [DATA.primary.argo, DATA.nearsurf.argo, DATA.secondary.argo,DATA.other.argo, ...
        DATA.primary.Dim, DATA.nearsurf.Dim, DATA.secondary.Dim,DATA.other.Dim,...
        argo, ~] = DOXY_argo_read(CONFIG.DataDir,Work.wmo,Work.replist,Work.whichCorr,Work.dirLog,Work.metafile);
    if isfield(DATA.primary.argo,'no_inair_data') && strcmp(Work.whichCorr,'INAIR')
        warndlg({sprintf('The wmo %d can''t be corrected with IN AIR method, ',...
            str2double(wmoIn{1})); '';...
            '=> Choose another one or try WOA or REF correction'},...
            'IN AIR correction not possible');
        continue
    end

    
    %Prepare data to compute drift on NCEP %marine 20/06/19
    if CONFIG.ok_inair_drift==1 && ~strcmp(Work.whichCorr,'INAIR')
        Work_NCEP=Work;
        Work_NCEP.whichDrift='NCEP';          
        Work_NCEP.is_not_inair_corr=1;
        [DATA_NCEP.primary.argo, DATA_NCEP.nearsurf.argo, DATA_NCEP.secondary.argo,DATA_NCEP.other.argo, ...
        DATA_NCEP.primary.Dim, DATA_NCEP.nearsurf.Dim, DATA_NCEP.secondary.Dim,DATA_NCEP.other.Dim,...
        argo_NCEP, ~] = DOXY_argo_read(CONFIG.DataDir,Work_NCEP.wmo,Work_NCEP.replist,Work_NCEP.whichDrift,Work_NCEP.dirLog,Work_NCEP.metafile);  
 
        if isfield(DATA_NCEP.primary.argo,'no_inair_data')
            no_inair_drift=1;
            h = warndlg({sprintf('No PPOX data or in-air measurement for %d !!!',Work.wmo),...
                ' ',...
                'Impossible to compute drift with NCEP',...
                '=> Drift will be computed on WOA'});
            uiwait(h);        
            clear DATA_NCEP
        end
    end   
    
    Work.launchdate=argo.launchdate;
    % Added by T. Reynaud 11.04.2022:
    save('tmp.mat','Work');
    
    % *****************************************************************
    % Apply a pressure effect correction 
    % *****************************************************************    

    a=questdlg('Do you want to apply a pressure effect correction on this float ?',sprintf('PRESSURE EFFECT CORRECTION - %d',Work.wmo),'Yes','No','No');   
    if strcmp(a,'Yes')
        %figure
        m_outerpos = [0.52 0.15 0.24 0.50];
        fig=figure('unit','normalized','Outerposition',m_outerpos,...
            'Name', sprintf('PRESSURE EFFECT CORRECTION - %d',Work.wmo),'NumberTitle','off','visible','off');
        hold on; grid on; hg1=hggroup;
        if strcmp(argo.VSS,'Near-surface sampling')
            if isempty(DATA.primary.argo.doxy.data)
                plot(DATA.secondary.argo.doxy.data',DATA.secondary.argo.pres.data','b','Parent',hg1)                
            else
                plot(DATA.primary.argo.doxy.data',DATA.primary.argo.pres.data','b','Parent',hg1)
            end
        else
            plot(argo.doxy.data',argo.pres.data','b','Parent',hg1)
        end

        %Apply the pressure effect correction to the data
        Work.anteriorcorr=1;
        Work.coeff_corr = cell2mat(inputdlg(sprintf('Enter coefficient value for pressure effect correction with the formula : \n    "DOXY_corr=DOXY*(1+coeff*PRES/1000)" \n '),sprintf('PRESSURE EFFECT CORRECTION - %d',Work.wmo)));
        if isempty(Work.coeff_corr)
            Work.anteriorcorr=0;
            fprintf('INFO >>>>>   No pressure effect applied \n')
            % ----------    SUMMARY FILE     -------------
            fprintf(log_summ,'No presure effect applied.\n');
            fprintf(log_summ,'\n');
            fprintf(log_summ,'----------------------------------------------------\n');
            fprintf(log_summ,'\n');
        else
            fprintf('INFO >>>>>   A pressure effect is applied \n')
            %apply presure effect
            DATA.primary.argo.doxy.data=DATA.primary.argo.doxy.data.*(1+DATA.primary.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            DATA.secondary.argo.doxy.data=DATA.secondary.argo.doxy.data.*(1+DATA.secondary.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            DATA.nearsurf.argo.doxy.data=DATA.nearsurf.argo.doxy.data.*(1+DATA.nearsurf.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            DATA.other.argo.doxy.data=DATA.other.argo.doxy.data.*(1+DATA.other.argo.pres.data.*str2double(Work.coeff_corr)/1000);
            argo.doxy.data=argo.doxy.data.*(1+argo.pres.data.*str2double(Work.coeff_corr)/1000);
            
            %Add corrected data to the figure
            hg2=hggroup;
            if strcmp(argo.VSS,'Near-surface sampling')
                if isempty(DATA.primary.argo.doxy.data)
                    plot(DATA.secondary.argo.doxy.data',DATA.secondary.argo.pres.data','r','Parent',hg2)
                else
                    plot(DATA.primary.argo.doxy.data',DATA.primary.argo.pres.data','r','Parent',hg2)
                end
            else
                plot(argo.doxy.data',argo.pres.data','r','Parent',hg2)
            end
           
          %  plot(DATA.primary.argo.doxy.data',DATA.primary.argo.pres.data','r','Parent',hg2)
            set(get(get(hg1,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            set(get(get(hg2,'Annotation'),'LegendInformation'),...hD
                'IconDisplayStyle','on');
            legend('Raw data','Raw data corrected from pressure effect','Location','Southwest')
            title(['Pressure effect correction with coeff=' num2str(Work.coeff_corr)])
            set(gca,'YDir','reverse');
            set(gca,'fontweight','bold');
            xlabel('[O2] (umol/kg)','fontweight','bold','fontsize',Work.fontsize+2);
            ylabel('Pres (dbar)','fontweight','bold','fontsize',Work.fontsize+2);
            set(fig,'visible','on')
            
            warning('From now on, every "raw" data presented will be corrected of the pressure effect');
            
            % ----------    SUMMARY FILE     -------------
            fprintf(log_summ,sprintf('A presure effect was applied : \n  DOXY_corr=DOXY*(1+coeff*PRES/1000)  with coeff=%s\n',num2str(Work.coeff_corr)));
            fprintf(log_summ,'\n');            
            fprintf(log_summ,'----------------------------------------------------\n');
            fprintf(log_summ,'\n');
            
            if Work.savePlot == 1
                hFig=gcf;
                saveFile = fullfile(Work.dirPlot,sprintf('DOXY_pressure_effect_coeff%s_%d',num2str(Work.coeff_corr),Work.wmo));
                DOXY_PLOT_settingsToPrint(hFig,Work,saveFile);       
            end
        end
    else
       Work.anteriorcorr=0;  
       Work.coeff_corr='';% Modified by VT+TR 19.01.2021
       fprintf('INFO >>>>>    No pressure effect applied \n')   
       
       % ----------    SUMMARY FILE     -------------  
       fprintf(log_summ,'No presure effect applied.\n');   
       fprintf(log_summ,'\n');         
       fprintf(log_summ,'----------------------------------------------------\n');     
       fprintf(log_summ,'\n');  
    end
    
      DATA.primary.Work = Work;
      DATA.nearsurf.Work = Work;
      DATA.secondary.Work = Work;
      DATA.other.Work = Work;
    
    if argo.n_prof == 0
        warnstr = {sprintf('%d : the main profile is empty (here : %s). ',REF_ARGO.wmo,argo.VSS);...
            ' ';...
            '=> No correction could be applied.';' '};
        h = warndlg(warnstr,'DOXY_argo_read WARNING');
        uiwait(h);
        continue
    end
    
    
    [argo1Struct, argo2Struct, argo3Struct, argo4Struct, ...
        argoWork, Work] = DOXY_argo_prepare_main(DATA, Work, argo);

    % *****************************************************************
    % Read trajectory data : if no in-air measurement, warn and stop
    % *****************************************************************
    noIAmeas = 0;
    if strcmp(Work.whichCorr,'INAIR')
        [~, ~, argoTrajWork] = DOXY_argoTraj_read(CONFIG,Work.wmo);
        Work.trajVar = argoTrajWork.trajVar;
        argoTrajWork = rmfield(argoTrajWork,'trajVar');
        
        if ~isfield(argoTrajWork,'ppox_doxy_adjusted')
            noIAmeas = 1;        
            h = warndlg({sprintf('No PPOX data for %d !!!',Work.wmo),...
                ' ',...
                'Impossible to use the INAIR correction',...
                '=> Try with WOA or REF correction'},'DOXY_argoTraj_read : no PPOX data');
            uiwait(h);        
            continue
        
        elseif all(cellfun('isempty',argoTrajWork.juld.data))
            noIAmeas = 1;
            h = warndlg({sprintf('No Inair Measurement for %d !!!',Work.wmo),...
                ' ',...
                sprintf('INFO: Measurement code for in-air measurement get in [%s]',num2str(CONFIG.inAirMC)),...
                '',...
                'Impossible to use the INAIR correction',...
                '=> Try with WOA or REF correction'},'DOXY_argoTraj_read : no Inair Measurement');
            uiwait(h);
            continue
        end
    elseif CONFIG.ok_inair_drift==1 && ~exist('no_inair_drift','var') %Prepare data to compute drift on NCEP %marine 20/06/19
        [~, ~, argoTrajWork_NCEP] = DOXY_argoTraj_read(CONFIG,Work_NCEP.wmo);
        Work_NCEP.trajVar = argoTrajWork_NCEP.trajVar;
        argoTrajWork_NCEP = rmfield(argoTrajWork_NCEP,'trajVar');
        
        if ~isfield(argoTrajWork_NCEP,'ppox_doxy_adjusted') || all(cellfun('isempty',argoTrajWork_NCEP.juld.data))
            noIAmeas_NCEP = 1;        
            h = warndlg({sprintf('No PPOX data or in-air measurement for %d !!!',Work.wmo),...
                ' ',...
                'Impossible to compute drift with NCEP',...
                '=> Drift will be compute on WOA'});
            uiwait(h);        
            clear argoTrajWork_NCEP Work_NCEP noIAmeas_NCEP
            no_inair_drift=1;
        else
            Work_NCEP.makePlot=0;  
            Work_NCEP.is_not_inair_corr=1;
            DATA_NCEP.primary.Work = Work_NCEP;
            DATA_NCEP.nearsurf.Work = Work_NCEP;
            DATA_NCEP.secondary.Work = Work_NCEP;
            DATA_NCEP.other.Work = Work_NCEP;
            [argo1Struct_NCEP, argo2Struct_NCEP, argo3Struct_NCEP, argo4Struct_NCEP, ...
        argoWork_NCEP, Work_NCEP] = DOXY_argo_prepare_main(DATA_NCEP, Work_NCEP, argo_NCEP);
        end
    end
    % =====================================================================
    %% Plot argo trajectory
    % =====================================================================
    if CONFIG.M_MAP_ACTIVE==0
        if strcmp(Work.whichCorr,'WOA') || strcmp(Work.whichCorr,'INAIR')
            DOXY_MAP(argo,Work,1)
        else
            iok = strcmp(REF.id,REF_ARGO.refId);
            DOXY_MAP(argo,Work,3,REF,iok)
        end
    end
    
    if CONFIG.M_MAP_ACTIVE>0
        if strcmp(Work.whichCorr,'WOA') || strcmp(Work.whichCorr,'INAIR')
            DOXY_MAP2(CONFIG,argo,Work,1)
        else
            iok=logical(zeros(size(REF.id)));%Added by T. Reynaud 01/02/2022
            for k=1:size(REF_ARGO.refId,2)
                iok(strcmp(REF.id,REF_ARGO.refId{k})==1)=1;
            end
            % End modifications TR.
            DOXY_MAP2(CONFIG,argo,Work,3,REF,iok)
        end
    end
    
    % =====================================================================
    %% Choose the vertical scale : pressure or density ?
    % =====================================================================
    Work.PresOrDens = 'Pressure';

    % =====================================================================
    %% Correction
    % =====================================================================
    % Probably : the DOXY_corr_compute and DOXY_corr_apply_main is the
    % same for woa and ref : keep it outside, and replace DOXY_woa_corr
    % and DOXY_ref_corr by DOXY_woa_corr_prepare and
    % DOXY_ref_corr_prepare
    
    argo1Structini=argo1Struct;
    argo2Structini=argo2Struct;
    argo3Structini=argo3Struct;
    argo4Structini=argo4Struct;
    
    switch Work.whichCorr
        case {'WOA','REF'}
            if CONFIG.ok_inair_drift==1 && ~exist('no_inair_drift','var') 
                NCEP.init = DOXY_NCEP_read(CONFIG,CONFIG.ncepYears);
                [Work_NCEP] = DOXY_inair_prepare_for_drift(CONFIG, WOA,...
                    NCEP, argo1Struct_NCEP, argo2Struct_NCEP, argo3Struct_NCEP, argo4Struct_NCEP, argo_NCEP, argoWork_NCEP, Work_NCEP, argoTrajWork_NCEP);
                clear argo1Struct_NCEP argo2Struct_NCEP argo3Struct_NCEP argo4Struct_NCEP 
                Work.NCEP_drift.Work=Work_NCEP;
                Work.NCEP_drift.argo=argo_NCEP;
                Work.NCEP_drift.argoWork=argoWork_NCEP;
                Work.NCEP_drift.argoTrajWork=argoTrajWork_NCEP;
                clear Work_NCEP argo_NCEP argoWork_NCEP argoTrajWork_NCEP
            end
            
            if strcmp(Work.whichCorr,'WOA')
                [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA,Work,goProg] = DOXY_woa_corr(CONFIG, ...
                WOA, argo1Struct, argo2Struct, argo3Struct, argo4Struct, argo, argoWork, Work);
            else 
                [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA,Work,goProg] = DOXY_ref_corr(CONFIG, ...
                WOA, REF_ARGO, argo1Struct, argo2Struct, argo3Struct, argo4Struct,argo, argoWork, Work);
            end
            
        case 'INAIR'
            % Read the NCEP data
            NCEP.init = DOXY_NCEP_read(CONFIG,CONFIG.ncepYears);
            [argo1Struct,argo2Struct,argo3Struct,argo4Struct,WOA,NCEP,Work,goProg] = DOXY_inair_corr(CONFIG,...
                WOA, NCEP, argo1Struct, argo2Struct, argo3Struct, argo4Struct, argo, argoWork, Work, argoTrajWork);
    end
    if goProg == 0
        continue;
    end
    
    % =====================================================================
    %% Apply correction
    % =====================================================================
    
    if Work.DODRIFT==true
        %Work.drift_val(num_fic)
        if Work.drift_fitPolynomialDegree == 1
            if ~Work.drift_PWLF % Added T.Reynaud 18.04.2024
                a=Work.PPOX_DRIFTCORR_COEF(1);%at+b
                b=Work.PPOX_DRIFTCORR_COEF(2);
                c=0;%08.07.2021 added by T. Reynaud
            else
                a=Work.PPOX_DRIFTCORR_COEF(1,:);%at+b
                b=Work.PPOX_DRIFTCORR_COEF(2,:);
                c=zeros(size(b));%18.04.2022 added by T. Reynaud
            end
        elseif Work.drift_fitPolynomialDegree == 2
            c=Work.PPOX_DRIFTCORR_COEF(1);%ct^2+at+b
            a=Work.PPOX_DRIFTCORR_COEF(2);
            b=Work.PPOX_DRIFTCORR_COEF(3);%08.07.2021 added by T. Reynaud            
        else
            disp('Problem, drift equation not in the type at+b!');
        end
    else
        a=0;
        b=1;
        c=0;%added by T. Reynaud 08.07.2021
    end
    if ~isfield(Work,['OFFSET_' Work.whichO2quantity])
        G=Work.(['FIT_slope_' Work.whichO2quantity]);
        OFFSET=Work.(['FIT_intercept_' Work.whichO2quantity]);
    else
        G=1;
        OFFSET=Work.(['OFFSET_' Work.whichO2quantity]);
    end
    
    % Added by T. Reynaud 06/07/2021
    if Work.drift_fitPolynomialDegree == 1
        SLOPE=b*G;
        DRIFT=100*365*G.*a./SLOPE;
        DRIFT2=0;%added by T. Reynaud 08.07.2021
    elseif Work.drift_fitPolynomialDegree == 2
        SLOPE=b*G;
        DRIFT=100*365*G*a/SLOPE;
        DRIFT2=100*(365)^2*G*c/SLOPE;%added by T. Reynaud 08.07.2021
    else
        SLOPE=0;
        DRIFT=0;
        DRIFT2=0;
    end

%  For H. Bittig/Schmichtig calculations
%     OFFSET=0;
%     SLOPE=1.0740;
%     DRIFT=1.110;
%     INCLINE_T=0;
%     DRIFT2=0;


    Work.SLOPE=SLOPE;
    Work.DRIFT=DRIFT;
    Work.DRIFT2=DRIFT2;
    if ~isfield(Work,'drift_PWLF')%Added by T. Reynaud 10.01.2025
        Work.drift_PWLF=0;
    end
    if Work.drift_PWLF
        Work.TIME_PWLF=Work.PPOX_DRIFTCORR_SEG(2);
    else
        Work.TIME_PWLF=0;
    end
    TIME_PWLF=Work.TIME_PWLF;
    % end TR
    INCLINE_T=0;
    Work.INCLINE_T=INCLINE_T;
    Work.OFFSET=OFFSET;
    
    disp('La correction est ')
    disp(['OFFSET = ', num2str(OFFSET)]);
    disp(['SLOPE = ', num2str(SLOPE)]);
    disp(['DRIFT = ', num2str(DRIFT)]);
    disp(['DRIFT2 = ', num2str(DRIFT2)]);
    disp(['INCLINE_T = ', num2str(INCLINE_T)]);
    disp(['TIME_PWLF = ', num2str(TIME_PWLF)]);
    
    for n = 1:4
    argoStruct = eval(['argo' num2str(n) 'Struct;']);  
    argoStructini = eval(['argo' num2str(n) 'Structini;']); 
    if argoStruct.argo.n_prof > 0 && isempty(argoStructini.argo.juld.data) == 0
        diffday=argoStructini.argo.juld.data-argo.launchdate*ones(1,size(argoStruct.argoWork.ppox_adjusted.data,2));
        
        GNTD=1;
        GNTD=GNTD+DRIFT(1)/100*diffday/365;
        GNTD=GNTD+DRIFT2/100.*diffday.*diffday/365/365;
        GNTD=GNTD+INCLINE_T.*double(argoStruct.argoWork.temp_adjusted.data);% double added by TR

        SLOPE_TIME=SLOPE(1)*ones(size(GNTD));

        if Work.drift_PWLF % Added T. Reynaud 19.04.2024

            GNTD2=1;
            GNTD2=GNTD2+DRIFT(2)/100*diffday/365;
            GNTD2=GNTD2+DRIFT2/100.*diffday.*diffday/365/365;
            GNTD2=GNTD2+INCLINE_T.*double(argoStruct.argoWork.temp_adjusted.data);% double added by TR

            idx=find(diffday(:,1)>=Work.PPOX_DRIFTCORR_SEG(2));
            SLOPE_TIME(idx,:)=SLOPE(2);

%             x1=diffday(:,1);
%             y1=double(GNTD(:,1));
%             x1=x1(~isnan(y1));
%             y1=y1(~isnan(y1));
%             [xs,ii]=sort(x1);
%             x1=xs;
%             y1=y1(ii);
%             p1=polyfit(x1,y1,1);          
            %test1=SLOPE(1)*(1+(DRIFT(1)/100*Work.PPOX_DRIFTCORR_SEG(2)/365))
            %test2=SLOPE(2)*(1+(DRIFT(2)/100*Work.PPOX_DRIFTCORR_SEG(2)/365))
            %time_check=Work.PPOX_daydiff;
            %gain_check=G*Work.PPOX_fitRegression;
            % test1=test2=gain_check(10) ==> C bon

            GNTD3=GNTD;
            GNTD3(idx,:)=GNTD2(idx,:);

%             h0=plot(Work.PPOX_daydiff,G*Work.PPOX_fitRegression,'-+g');hold on
%             h3=plot(diffday(:,1),SLOPE_TIME(:,1).*GNTD3(:,1),'*m');
%             ylabel('Gain');
%             xlabel('Days');
%             title(num2str(Work.wmo));
%             hleg=legend([h0,h3],'From Doxy\_drift','Gain Final');
%             hleg.FontSize=15;

            GNTD=GNTD3;
            clear GNTD2 GNTD3
        end
        
        if isfield(Work,'ind_drift_stop')
           idx=find(diffday>=Work.ind_drift_stop(2), 1, 'first');
           for ic=idx:size(GNTD,1)
             GNTD(ic,:)=GNTD(idx,:);
           end
        end
        
        argoStruct.argoWork.ppox_adjusted.data = OFFSET+SLOPE_TIME.*GNTD.* ...
                                                 argoStructini.argoWork.ppox_adjusted.data;
        
        if Work.presEff == 1, P = argoStructini.pres_adjusted.data; else, P = 0; end
        
        tmpDoxy=O2ptoO2c(argoStruct.argoWork.ppox_adjusted.data,argoStruct.argoWork.temp_adjusted.data,argoStruct.argoWork.psal_adjusted.data,P);                                    
        argoStruct.argoWork.doxy_adjusted.data = DOXY_convert(tmpDoxy,'mumol/L',argoStruct.argoWork.doxy_adjusted.units,argoStruct.argoWork.an_dens.data);
        argoStruct.argoWork.psat_adjusted.data = O2ptoO2s(argoStruct.argoWork.ppox_adjusted.data,argoStruct.argoWork.temp_adjusted.data,argoStruct.argoWork.psal_adjusted.data,P);
       
    end
    eval(['argo' num2str(n) 'Struct = argoStruct;']);
    end
    


    
    % =====================================================================
    %% Save structures for wmo
    % =====================================================================
    %marine
    if Work.presEff
        presEffStr = 'okpreseff'; 
    else
        presEffStr = 'nopreseff'; 
    end
    % Modified T. Reynaud 19.01.2021
    %if isempty(Work.coeff_corr) || Work.coeff_corr==0
    if isempty(Work.coeff_corr) || str2double(Work.coeff_corr)==0
        presCorrStr='nopresCorr';
    else
        % Modified T. Reynaud 19.01.2021
        %presCorrStr=['presCorr' num2str(Work.coeff_corr)];
        presCorrStr=['presCorr' (Work.coeff_corr)];
    end
        
    if Work.DODRIFT
        driftStr = ['drifton' Work.whichDrift];
    else
        driftStr = 'nodrift';
    end
    if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0, offsetStr = 'nooffset';
    else, offsetStr = 'okoffset';
    end
    if isfield(argo1Struct.Work,['OFFSET_' argo1Struct.Work.whichO2quantity])
        corr_type=[Work.whichCorr 'constant'];
    else
        corr_type=Work.whichCorr;
    end
    Work.dirout=sprintf('%s_%s_%s_%s_%s_%s',corr_type,presEffStr,presCorrStr,driftStr,offsetStr,Work.whichO2quantity);
    
    if ~exist(fullfile(CONFIG.saveDataDir,'MAT',Work.dirout),'dir')
        mkdir(fullfile(CONFIG.saveDataDir,'MAT',Work.dirout));
    end 
    save(fullfile(CONFIG.saveDataDir,'MAT',Work.dirout,...
        sprintf('%d_corr',Work.wmo)),'argo','argoWork','Work',...
        'argo1Struct','argo2Struct','argo3Struct','argo4Struct');
    
    % =====================================================================
    %% Write NetCDF corrected
    % =====================================================================
    % Prepare the data to be writen
    [argo1Struct, argo2Struct, argo3Struct, argo4Struct, CONFIG] = DOXY_write_prepare(argo1Struct, ...
        argo2Struct, argo3Struct, argo4Struct, CONFIG, Work);
      
    % Write the new monoprofiles
    if strcmp(Work.whichCorr,'REF')
        DOXY_argo_write(CONFIG.DataDir,Work.replist,CONFIG.saveDataDir,...
        CONFIG.prefix,argo1Struct, argo2Struct, argo3Struct, argo4Struct, Work,REF_ARGO);  
    else
        DOXY_argo_write(CONFIG.DataDir,Work.replist,CONFIG.saveDataDir,...
        CONFIG.prefix,argo1Struct, argo2Struct, argo3Struct, argo4Struct, Work);  
    end
    
    % --------   SUMMARY FILE   --------------
    fprintf(log_summ,'\n');
    fprintf(log_summ,'             CORRECTION COMPLETED\n');
    
    
    % close the summary file
    fclose('all');

    % Added by T. Reynaud 11.04.2022
    if exist('tmp.mat','file')
        delete('tmp.mat');
    end

    newWMO = questdlg('Do you want to manage a new wmo ?',...
        sprintf('%s: Continue the correction', CONFIG.history_software),'YES',...
        sprintf('NO, STOP %s',CONFIG.history_software),sprintf('NO, STOP %s',CONFIG.history_software));
    if goProg == 0, 
        goProg =1; 
    end
    if strcmp(newWMO,sprintf('NO, STOP %s',CONFIG.history_software))
        goProg = 0;
        %msgbox({'  Good Bye !  ';''},sprintf('%s',CONFIG.history_software),'custom',imread(CONFIG.logo));
    else
        if strcmp(Work.whichCorr,'INAIR') && ~noIAmeas
            NCEP = rmfield(NCEP,'coloc');
        end
        Work = initialWork;
        if length(fieldnames(WOA)) == 2
            WOA = rmfield(WOA,'interp');
            clear DATA argo1Struct argo2Struct argo3Struct argo argoWork;
        end
    end
end


