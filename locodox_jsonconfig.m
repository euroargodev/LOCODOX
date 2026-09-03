% LOCODOX_json_config check the user configuration parameters usefull for the oxygen
% correction LOCODOX processing chain and merged with the tool configuration - 
% The user configuration come from the user json file
%
% SYNTAX
% [] = LOCODOX_json_config()
%
% DESCRIPTION
% LOCODOX_json_config return the configuration parameters usefull for the oxygen
% correction LOCODOX processing chain - need .json 
%
% INPUT
%   full path of the locodox_config.json
%
% OUTPUT
%   CONFIG (struct)       Structure bearing the configuration information.
%                           CONFIG = 
%                            LocodoxMainDir : 'C:\Users\...\LOCODOX\...'
%                                    DataDir: 'C:\Users\...\LOCODOX\...'
%                                NCEPDataDir: 'C:\Users\...\LOCODOX\...'
%                                resultsDir :  'C:\Users\...\LOCODOX\...'
%                                    WOAfile: 'WOA_monthly_5500_1deg.nc'
%                                    bddFile: 'bddo2ref.mat'
%                                RefArgoFile: 'bdd_ref_argo.txt'
%                                  bathyFile: 'etopo5.cdf'
%                                   varBathy: {'etopo05_x'  'etopo05_y'  'rose'}
%                                   maskFile: 'landsea_masks.cdf'
%                                    varMask: 'landsea'
%                                       logo: 'C:\Users\...\SOFT_2017...'
%                                    ncepFtp: 'ftp.cdc.noaa.gov'
%                                 ncepFtpDir: 'Datasets/ncep.reanalysis/'
%                              ncepFtpSubDir: {'surface'  'surface'  'surface'}
%                                  ncepFiles: {'slp'  'air.sig995'  'rhum.sig995'}
%                               ncepDoUpdate: 0
%                                  ncepYears: [1xn double]
%                               corrTypeDesc: {1x3 cell}
%                              corrTypeShort: {'WOA'  'REF'  'INAIR'}
%                                   refUnit : 'mumol/kg'
%                                   presEff : 0
%                                     isokC : 0
%                                    DM_pres: 1
%                                    DM_temp: 1
%                                    DM_psal: 0
%                                    QC_[..]: [1 2]
%                 drift_fitPolynomialDegree : 1
%                                drift_spec : 0
%                                 inWaterMC : [690 710]
%                                   inAirMC : [699 711 799]
%                         inAirMaxPresForTS : 20
%                                 R2treshold: 0.8000
%                         adjusted_error_raw: 0
%                         adjusted_error_rel: 1
%                               trajSpeedLim: 0.6000
%                           history_software: 'LOCODOX'
%                  history_software_release : '3.0'
%                          history_reference: 'LOCODOX2016'
%                                     prefix: 'BD'
%                                   makePlot: 1
%                                   savePlot: 0
%                                   fontsize: 8
%                                resolution : 100
%                                formattype : {'-dpng'}
%
% CALL
%   'configuration/locodox_default_user_config.json'
%
% SEE ALSO
%
%
% HISTORY
%   $created: //2026 $author: Virginie Racape, pokapok
%                   distinction between the :
%                       user config coming from locodox_user_config.json and from locodox_default_user_config.json
%                       tool config from locodox_tool_config.json
%   $Revision: 

function [CONFIG] = locodox_jsonconfig()

% -------------------------------------------------------------------------
% USER CONFIG FILE
% -------------------------------------------------------------------------
user_config_json='configuration/locodox_user_config.json';
if ~exist(user_config_json,'file')
    fprintf('STOP: %s does not exist, please copy/paste configuration/locodox_default_user_config.json, change its name by locodox_user_config.json and fill with your configuration\n', user_config_json)
    CONFIG = [];
    return
end

% -------------------------------------------------------------------------
% DEFINED FUNCTION FOR CONFIG
% -------------------------------------------------------------------------
% function to manage srtuc from jsons
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
% function to check config structure
checkEmpty = @(s)isempty(s);
checkNum = @(s)isnumeric(s);
% function to update fields
function S = updatestruct(S, U)
    % Met à jour les champs de S avec ceux de U.
    % Si un champ existe déjà dans S, sa valeur est écrasée par celle de U.
    % Si un champ n'existe pas dans S, il est ajouté.
    noms = fieldnames(U);
    for i = 1:numel(noms)
        S.(noms{i}) = U.(noms{i});
    end
end

% -------------------------------------------------------------------------
% LOAD JSON - user configuration
% -------------------------------------------------------------------------
% To be compatible with R2016
jsonText = fileread(user_config_json);
user_jsondata = jsondecode(jsonText);

% -------------------------------------------------------------------------
% LOAD JSON - default user configuration (template)
% -------------------------------------------------------------------------
% To be compatible with R2016
jsonText = fileread('configuration/locodox_default_user_config.json');
default_jsondata = jsondecode(jsonText);

% -------------------------------------------------------------------------
% LOAD JSON - tool configuration
% -------------------------------------------------------------------------
% To be compatible with R2016
jsonText = fileread('configuration/locodox_tool_config.json');
tool_jsondata = jsondecode(jsonText);

% -------------------------------------------------------------------------
% DIRECTORIES
% -------------------------------------------------------------------------
% check if the fields exist
if ~isfield(user_jsondata, 'directories')
    fprintf('STOP: directory section is mandatory in json\n')
    CONFIG = [];
    return
end

% check user content key values
default_jsondata.directories=rmfield(default_jsondata.directories,"comment");
dirnames = fieldnames(default_jsondata.directories);
user_dirnames = fieldnames(user_jsondata.directories);
isOk = ismember(dirnames,user_dirnames);
if any(isOk==0)
    fprintf('STOP: %s directory is mandatory in json\n',dirnames{isOk == 0})
    CONFIG = [];
    return
end

% check user content empty values
isOk = ismember(user_dirnames,dirnames);
isEmpty = structfun(checkEmpty, user_jsondata.directories);
if any(isEmpty == 1 & isOk == 1)
    fprintf('STOP: %s is mandatory and the path should not be empty\n',user_dirnames{[isEmpty == 1 & isOk == 1] == 1})
    CONFIG = [];
    return
end
if any(isOk==0)
    fprintf('INFO: %s should be removed from the user configuration \n',user_dirnames{isOk == 0})
    user_jsondata.directories = rmfield(user_jsondata.directories, user_dirnames(isOk==0));
    user_dirnames(isOk==0)=[];
end

% update CONFIG
CONFIG = user_jsondata.directories;

mytools = fullfile(CONFIG.LocodoxMainDir,'share','MyTools',filesep);
addpath(mytools);
hydcalDir = fullfile(CONFIG.LocodoxMainDir,'share','hydcal',filesep);
seawaterDir = fullfile(CONFIG.LocodoxMainDir,'share','seawater','seawater_330_its90_lpo',filesep);
addpath(hydcalDir);
addpath(seawaterDir);
addpath(fullfile(CONFIG.LocodoxMainDir,'share',filesep));
addpath(fullfile(CONFIG.LocodoxMainDir,'doxy_corr',filesep));
addpath(fullfile(CONFIG.RefBddDir,'data_input',filesep));

% -------------------------------------------------------------------------
% FILES
% -------------------------------------------------------------------------
% Mask file and the variable to be read in the mask file
CONFIG = mergestructs(CONFIG, tool_jsondata.mask_file);

% Logo file
CONFIG.logo = fullfile(CONFIG.LocodoxMainDir,tool_jsondata.logo_file);

% -------------------------------------------------------------------------
% DATA SELECTION
% -------------------------------------------------------------------------
key_val = {'data_mode_selection','QC_selection'};
for ikey = 1:length(key_val)

    % check user content key values
    default_jsondata.(key_val{ikey})=rmfield(default_jsondata.(key_val{ikey}),"comment");

    % check if the fields exist
    if ~isfield(user_jsondata, key_val{ikey})
        fprintf('INFO: %s is mandatory in json, it will be updated with default value\n',key_val{ikey})
        user_jsondata.(key_val{ikey}) = default_jsondata.(key_val{ikey});       
    end

    dirnames = fieldnames(default_jsondata.(key_val{ikey}));
    user_dirnames = fieldnames(user_jsondata.(key_val{ikey}));

    % is avail or nor
    isOk = ismember(dirnames,user_dirnames);
    if any(isOk==0)
        fprintf('WARN: %s is missing in json, it will be replace by the default value\n',dirnames{isOk == 0})
        U = rmfield(default_jsondata.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.(key_val{ikey})), dirnames(isOk == 0)));
        user_jsondata.(key_val{ikey}) = updatestruct(user_jsondata.(key_val{ikey}), U);
        user_dirnames = fieldnames(user_jsondata.(key_val{ikey}));
    end

    % check user content empty values or different from numeric
    isOk = ismember(user_dirnames,dirnames);
    isEmpty = structfun(checkEmpty, user_jsondata.(key_val{ikey}));
    if any(isEmpty == 1 & isOk == 1)
        fprintf('WARN: %s is empty, it will be replace by the default value\n',user_dirnames{[isEmpty == 1 & isOk == 1] == 1})
        U = rmfield(default_jsondata.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.(key_val{ikey})), user_dirnames{[isEmpty == 1 & isOk == 1] == 1}));
        user_jsondata.(key_val{ikey}) = updatestruct(user_jsondata.(key_val{ikey}), U);
    end
    
    if any(isOk==0)
        fprintf('INFO: %s should be removed from the user configuration \n',user_dirnames{isOk == 0})
        user_jsondata.(key_val{ikey}) = rmfield(user_jsondata.(key_val{ikey}), user_dirnames(isOk==0));
        user_dirnames(isOk==0)=[];
    end

    % numeric or not
    isNum = structfun(checkNum, user_jsondata.(key_val{ikey}));
    if any(isNum==0)
        fprintf('WARN: %s should be a numeric, it will be replace by the default value\n',user_dirnames{isNum == 0})
        U = rmfield(default_jsondata.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.(key_val{ikey})), user_dirnames{isNum == 0}));
        user_jsondata.(key_val{ikey}) = updatestruct(user_jsondata.(key_val{ikey}), U);
    end

    CONFIG = mergestructs(CONFIG, user_jsondata.(key_val{ikey}));
end

disp(['INFO: if the profile is Near-Surface, LOCODOX chooses PSAL', ...
       'instead of PSAL_ADJUSTED (whatever the user choice), because', ...
       'unpumped data get artificially a QC = 4 (so data = NaN), not', ... 
       'justified for DOXY calculation'])


CONFIG.QC_O = CONFIG.QC_O'; % DOXY
CONFIG.QC_P = CONFIG.QC_P'; % PRES
CONFIG.QC_T = CONFIG.QC_T'; % TEMP
CONFIG.QC_S = CONFIG.QC_S'; % PSAL

% -------------------------------------------------------------------------
% FOR CORRECTION
% -------------------------------------------------------------------------
key_val_num = {'general','drift','from_inair'};
key_val_str = {'from_reference','from_clim'};
key_val = [key_val_num, key_val_str];
for ikey = 1:length(key_val)

    % kind of key_val
    if ismember(key_val{ikey}, key_val_num)
        val = 1; % numeric
    elseif ismember(key_val{ikey}, key_val_str)
        val = 2 ;% string
    end

    % check user content key values
    default_jsondata.doxy_correction.(key_val{ikey})=rmfield(default_jsondata.doxy_correction.(key_val{ikey}),...
        "comment");
    
    % check if the fields exist
    if ~isfield(user_jsondata.doxy_correction, key_val{ikey})
        fprintf('INFO: %s is missing in json, it will be updated with default value\n',key_val{ikey})
        user_jsondata.doxy_correction.(key_val{ikey}) = default_jsondata.doxy_correction.(key_val{ikey});       
    end

    dirnames = fieldnames(default_jsondata.doxy_correction.(key_val{ikey}));
    user_dirnames = fieldnames(user_jsondata.doxy_correction.(key_val{ikey}));

    % is avail or nor
    isOk = ismember(dirnames,user_dirnames);
    if any(isOk==0)
        fprintf('WARN: %s is missing in json, it will be replace by the default value\n',dirnames{isOk == 0})
        U = rmfield(default_jsondata.doxy_correction.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.doxy_correction.(key_val{ikey})), dirnames(isOk == 0)));
        user_jsondata.doxy_correction.(key_val{ikey}) = updatestruct(user_jsondata.doxy_correction.(key_val{ikey}), U);
        user_dirnames = fieldnames(user_jsondata.doxy_correction.(key_val{ikey}));
    end

    % check user content empty values for expected numeric fields
    isOk = ismember(user_dirnames,dirnames);
    isEmpty = structfun(checkEmpty, user_jsondata.doxy_correction.(key_val{ikey}));
    if any(isEmpty == 1 & isOk==1) & val == 1
        fprintf('WARN: %s is empty, it will be replace by the default value\n',user_dirnames{[isEmpty == 1 & isOk==1]==1})
        U = rmfield(default_jsondata.doxy_correction.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.doxy_correction.(key_val{ikey})), user_dirnames{[isEmpty == 1 & isOk==1]==1}));
        user_jsondata.doxy_correction.(key_val{ikey}) = updatestruct(user_jsondata.doxy_correction.(key_val{ikey}), U);
    elseif any(isEmpty == 1 & isOk==1) & val == 2
        fprintf('WARN: %s is empty\n',user_dirnames{[isEmpty == 1 & isOk==1]==1})
        fprintf('WARN, in this way, the tool cannot use the %s option for doxy correction\n',key_val{ikey})
    end
    
    if any(isOk==0)
        fprintf('INFO: %s should be removed from the user configuration \n',user_dirnames{isOk == 0})
        user_jsondata.doxy_correction.(key_val{ikey}) = rmfield(user_jsondata.doxy_correction.(key_val{ikey}), user_dirnames(isOk==0));
        user_dirnames(isOk==0)=[];
    end

    % numeric or not
    isNum = structfun(checkNum, user_jsondata.doxy_correction.(key_val{ikey}));
    % exception managment
    if strcmp(key_val{ikey},'drift')
        isNum(strcmp(fieldnames(user_jsondata.doxy_correction.(key_val{ikey})),'drift_fittype')) = 1;
    end
    if any(isNum==0) & val == 1
        fprintf('WARN: %s should be a numeric, it will be replace by the default value\n',user_dirnames{isNum == 0})
        U = rmfield(default_jsondata.doxy_correction.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.doxy_correction.(key_val{ikey})), ...
                user_dirnames{isNum == 0}));
        user_jsondata.doxy_correction.(key_val{ikey}) = updatestruct(user_jsondata.doxy_correction.(key_val{ikey}), U);
    end

    CONFIG = mergestructs(CONFIG, user_jsondata.doxy_correction.(key_val{ikey}));
end

%update information
% for clim
CONFIG.WOAfile = [CONFIG.WOADataDir, CONFIG.WOA_file];
CONFIG = rmfield(CONFIG, 'WOA_file');
if isnumeric(CONFIG.WOA_YEAR)
    CONFIG.WOA_YEAR = num2str(CONFIG.WOA_YEAR)
end

% for in air
CONFIG = mergestructs(CONFIG, tool_jsondata.doxy_correction.from_inair);
CONFIG.ncepFtpSubDir = CONFIG.ncepFtpSubDir';
CONFIG.ncepFiles = CONFIG.ncepFiles';
CONFIG.ncepYears = [CONFIG.ncepYearStart:CONFIG.ncepYearEnd];
CONFIG.ncepGetYears = CONFIG.ncepYears;

CONFIG.inWaterMC = CONFIG.inWaterMC';
CONFIG.inAirMC = CONFIG.inAirMC';

% -------------------------------------------------------------------------
% FOR NetCDF writing and FOR HISTORY
% -------------------------------------------------------------------------
CONFIG = mergestructs(CONFIG, tool_jsondata.tool_version);
CONFIG.history_reference = CONFIG.WOA_YEAR;

% -------------------------------------------------------------------------
% GRAPHICAL INFORMATION / M_MAP INFORMATION / Ask online questions or not
% -------------------------------------------------------------------------
key_val = {'graphical_information','mmap_information','online_requeriment'};
for ikey = 1:length(key_val)

    % check user content key values
    default_jsondata.(key_val{ikey})=rmfield(default_jsondata.(key_val{ikey}),"comment");

    % check if the fields exist
    if ~isfield(user_jsondata, key_val{ikey})
        fprintf('INFO: %s is mandatory in json, it will be updated with default value\n',key_val{ikey})
        user_jsondata.(key_val{ikey}) = default_jsondata.(key_val{ikey});       
    end

    dirnames = fieldnames(default_jsondata.(key_val{ikey}));
    user_dirnames = fieldnames(user_jsondata.(key_val{ikey}));

    % is avail or nor
    isOk = ismember(dirnames,user_dirnames);
    if any(isOk==0)
        fprintf('WARN: %s is missing in json, it will be replace by the default value\n',dirnames{isOk == 0})
        U = rmfield(default_jsondata.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.(key_val{ikey})), dirnames(isOk == 0)));
        user_jsondata.(key_val{ikey}) = updatestruct(user_jsondata.(key_val{ikey}), U);
        user_dirnames = fieldnames(user_jsondata.(key_val{ikey}));
    end

    % check user content empty values or different from numeric
    isOk = ismember(user_dirnames,dirnames);
    isEmpty = structfun(checkEmpty, user_jsondata.(key_val{ikey}));
    if any(isEmpty == 1 & isOk == 1)
        fprintf('WARN: %s is empty, it will be replace by the default value\n',user_dirnames{[isEmpty == 1 & isOk == 1] == 1})
        U = rmfield(default_jsondata.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.(key_val{ikey})), user_dirnames{[isEmpty == 1 & isOk == 1] == 1}));
        user_jsondata.(key_val{ikey}) = updatestruct(user_jsondata.(key_val{ikey}), U);
    end
    
    if any(isOk==0)
        fprintf('INFO: %s should be removed from the user configuration \n',user_dirnames{isOk == 0})
        user_jsondata.(key_val{ikey}) = rmfield(user_jsondata.(key_val{ikey}), user_dirnames(isOk==0));
        user_dirnames(isOk==0)=[];
    end

    % numeric or not
    isNum = structfun(checkNum, user_jsondata.(key_val{ikey}));
    if strcmp(key_val{ikey},'graphical_information')
        isNum(strcmp(fieldnames(user_jsondata.(key_val{ikey})),'formattype')) = 1;
    end
    if any(isNum==0)
        fprintf('WARN: %s should be a numeric, it will be replace by the default value\n',user_dirnames{isNum == 0})
        U = rmfield(default_jsondata.(key_val{ikey}), ...
            setdiff(fieldnames(default_jsondata.(key_val{ikey})), user_dirnames{isNum == 0}));
        user_jsondata.(key_val{ikey}) = updatestruct(user_jsondata.(key_val{ikey}), U);
    end

    CONFIG = mergestructs(CONFIG, user_jsondata.(key_val{ikey}));
end

CONFIG.formattype=CONFIG.formattype';
mmp = fullfile(CONFIG.LocodoxMainDir,'share','m_map1.4m',filesep);
CONFIG.M_MAPDir=mmp;

disp('')
disp(' ---- USER CONFIG ------')
disp(CONFIG)
end