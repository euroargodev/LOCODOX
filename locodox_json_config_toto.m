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
%                       user config coming from locodox_config.json
%                       tool confing available here
%                   used from version 5.0
%   $Revision: 

function [CONFIG] = locodox_json_config(config_json)

% -------------------------------------------------------------------------
% DEFINED FUNCTION FOR CONFIG
% -------------------------------------------------------------------------
% function to merge srtuc from jsons
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
% function to check config structure
checkEmpty = @(s)isempty(s);
checkNum = @(s)isnumeric(s);

% -------------------------------------------------------------------------
% LOAD JSON - user configuration
% -------------------------------------------------------------------------
% To be compatible with R2016
jsonText = fileread(config_json);
user_jsondata = jsondecode(jsonText);

% -------------------------------------------------------------------------
% LOAD JSON - default configuration
% -------------------------------------------------------------------------
% To be compatible with R2016
jsonText = fileread('configuration/locodox_default_user_config.json');
default_jsondata = jsondecode(jsonText);

% -------------------------------------------------------------------------
% LOAD JSON - tool configuration
% -------------------------------------------------------------------------
% To be compatible with R2016
jsonText = fileread('configuration/locodox_config.json');
jsondata = jsondecode(jsonText);

% -------------------------------------------------------------------------
% DIRECTORIES
% -------------------------------------------------------------------------
% check user content key values
dirnames = fieldnames(default_jsondata.directories);
user_dirnames = fieldnames(user_jsondata.directories);
isOk = ismember(dirnames,user_dirnames);
if any(isOk==0)
    fprintf('%s directory/ies is/are mandatory in json\n',dirnames{isOk == 0})
    return
end

% check user content empty values
isEmpty = structfun(checkEmpty, user_jsondata.directories);
if any(isEmpty)
    fprintf('%s is/are mandatory and not empty\n',user_dirnames{isEmpty == 1})
    return
end

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
CONFIG = mergestructs(CONFIG, jsondata.mask_file);

% Logo file
CONFIG.logo = fullfile(CONFIG.LocodoxMainDir,jsondata.logo_file);

% -------------------------------------------------------------------------
% DATA SELECTION
% -------------------------------------------------------------------------
% check user content key values
key_val = {'data_mode_selection','QC_selection'}
for ikey = 1:length(key_val)
    dirnames = fieldnames(default_jsondata.(key_vale{ikey}));
    user_dirnames = fieldnames(user_jsondata.(key_vale{ikey}));
    % is avail or nor
    isOk = ismember(dirnames,user_dirnames);
    if any(isOk==0)
        fprintf('WARN: %s is/are mandatory in json\n',dirnames{isOk == 0})
        fprintf('The tool replaces by defaults values %s\n',default_jsondata.(key_vale{ikey}).dirnames{isOk == 0})
    end
    % numeric or not
    % isOk = structfun(dirnames,user_dirnames);
    % if any(isOk==0)
    %     fprintf('WARN: %s is/are mandatory in json\n',dirnames{isOk == 0})
    %     fprintf('The tool replaces by defaults values %s\n',default_jsondata.(key_vale{ikey}).dirnames{isOk == 0})
    % end

% Data mode selection
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Do you want to force the use of the 'Real Time' fields even if the 
% 'Delayed Mode' fields exist?
% If yes, set the DM_fields to 0 (By default: DM_fields = 1 and the 
% 'Delayed Mode' fields are privileged) 
% (CONFIG.DM : 1=ADJUSTED ; 0=RAW)
% Note: WARNING, if the profile is Near-Surface, LOCODOX chooses PSAL 
%       instead of PSAL_ADJUSTED (whatever the user choice), because 
%       unpumped data get artificially a QC = 4 (so data = NaN), not 
%       justified for DOXY calculation 
CONFIG = mergestructs(CONFIG, jsondata.data_mode_selection);
 
% QC selection: 
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Only the data whose QC equal to these values are 
% retained for the calculation of the correction (for DOXY, PRES, TEMP and PSAL).
% Format: Example: CONFIG.QC_x = [1 2 3 4] (possible values: from 0 to 9)
CONFIG = mergestructs(CONFIG, jsondata.QC_selection);
CONFIG.QC_O = CONFIG.QC_O'; % DOXY
CONFIG.QC_P = CONFIG.QC_P'; % PRES
CONFIG.QC_T = CONFIG.QC_T'; % TEMP
CONFIG.QC_S = CONFIG.QC_S'; % PSAL

% -------------------------------------------------------------------------
% FOR CORRECTION
% -------------------------------------------------------------------------

% General information
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% include
% - Conversion DOXY/PSAT/PPOX : take into account pressure effect or not : 0/1 option unactivated/activated
% - Linear correction R2 coefficient : treshold under which LOCODOX suggests
% to apply Constant correction
CONFIG = mergestructs(CONFIG, jsondata.doxy_correction.general);

% % For FSD floats: Hybrid DM a,d RT
% % Comment otherwise
% % CONFIG.RT_psal_cycle=158;% Pour 6902806 FSD use RT PSAL data from specified cycle ==> DM_PSAL=1
% % CONFIG.QC_S = [1 2 3 4]; % PSAL FSD 6902806 ==> DP_psal=1


% Time Drift correction
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% - Compute inair drift if possible for WOA (default) and REF correction (in addition)
% - near surface (Introduced by TR 12/04/2021) & deep waters (modified by TR 12/04/2021) min, max, step info 
% - Added by T. Reynaud for Piece wise linear fitting for Time Drift 
% 26/04/2024
CONFIG = mergestructs(CONFIG, jsondata.doxy_correction.drift);
% By default, the data drift is computed using the polynomial fitting
% (order 1, y ~ ax+b). The polynomial degree can be increased up to 3. 
% A new equation could be used if the parameter CONFIG.drift_spec is 
% activated and if the parameter CONFIG.drift_fittype is filled.
CONFIG.drift_fitPolynomialDegree = 1;
CONFIG.drift_spec = 0;
if CONFIG.drift_spec == 1
    % If the fitting equation is a classical one (ex : a*exp(b*x)), enter
    % the name. If it is not a classical one, create the new equation using
    % the matlab "fittype" fonction defining way (see the matlab help).
    % Example : g = fittype('a*u+b*exp(n*u)','problem','n','independent','u')
    CONFIG.drift_fittype = fittype('exp1');
    %CONFIG.drift_fittype = fittype('exp1');
end

% from reference data set
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% - BDD files including list assiciated wmo with reference ctd
% - Reference data unit : precise the unit of your reference data (most of the
% time in-situ CTD)
% Possible options:
%         mL/L, milliliter per liter, milliliter/L
%         mumol/m3, micromole per m3, mmol/m3, micromole/m3
%         mumol/L, micromole per liter, mmol/L, micromole/L
%         mg/L, milligram per liter, milligram/L
%         mumol/kg, micromole per kilo, mmol/kg, micromole per kilogram, micromole/kg 
CONFIG = mergestructs(CONFIG, jsondata.doxy_correction.from_reference);

% from in air
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% - Carry over parameter. If isokC is set to 0, we suppose that in air oxygen measurements are not biased by splash of water; C is set to 1 otherwise
% See "Oxygen Optode Sensors : Principle Characterization, Calibration, and
% Application in the Ocean", Henry Bittig and al (2018)
% - ncep update info
CONFIG = mergestructs(CONFIG, jsondata.doxy_correction.from_inair);

% additional needs
% NCEP files and FTP address
% ncepFtp : the ftp website of NCEP
% ncepFtpDir : the path where to find the NCEP data in the ftp website
% ncepFtpSubDir : the sub directory where to find the NCEP data in the ftp website
% ncepFiles : the NCEP files to be read
% ncepYears : read the NCEP data for the years specified
CONFIG.ncepFtp = 'ftp.cdc.noaa.gov';
CONFIG.ncepFtpDir = 'Datasets/ncep.reanalysis/';
CONFIG.ncepFtpSubDir = {'surface','surface','surface'};
CONFIG.ncepFiles = {'slp','air.sig995','rhum.sig995'};
CONFIG.ncepYears = [CONFIG.ncepYearStart:CONFIG.ncepYearEnd];% To be used
CONFIG.ncepGetYears = str2double(datestr(now,'YYYY'));% To be downloaded
%CONFIG.ncepGetYears = [2021:2022];% To be downloaded

% trajectory data : measurement codes for In-Air and Near-surface samples
% Measurement Code which carries the InAir measurement :
%  - CONFIG.inAirMC : MC for the in-air samples, part of the surface
%  sequence (double, array)
%  - CONFIG.inWaterMC : MC for the in-water samples, part of the surface
%  sequence (or the profile sequence).(double, array)
CONFIG.inWaterMC = [690 710];
CONFIG.inAirMC = [699 711 799];

% from climatology
% WOA climatology file
WOA_YEAR='2018';% Release Year
CONFIG.WOAfile = [CONFIG.WOADataDir,'WOA2018_DECAV_monthly_5500_1deg.nc'];

% % Error on adjusted field
% % '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% % abs (absolute): adjusted_error is a fixed value
% % rel (relative): adjusted error is a % value:  doxy_adjusted + adjusted_error_rel*doxy_adjusted/100
CONFIG = mergestructs(CONFIG, jsondata.doxy_correction.general.error);

% -------------------------------------------------------------------------
% FOR NetCDF writing and FOR HISTORY
% -------------------------------------------------------------------------
CONFIG.history_software = 'LOCODOX';
CONFIG.history_reference = ['LOPS2020_WOA',WOA_YEAR];
CONFIG.history_software_release = '5.0';
CONFIG.prefix = 'BD';

% -------------------------------------------------------------------------
% GRAPHICAL INFORMATION
% -------------------------------------------------------------------------
% Save Plots        : 1=yes, 0=no
% Font size         : tunes the font size on the plot, for both screen display and
%                     saving plots
% Resolution in dpi : 0 : screen resolution, 100 : mid resolution, 300 : high resolution
%                     warning: the higher the resolution setting, the
%                     longer it takes to render your figure
% Format for save   : '-dpng', '-jpeg', '-ps',... all the fromattype
%                      allowed by matlab for the print function (doc print
%                      for more information). Cell array of string.
% ------
CONFIG = mergestructs(CONFIG, jsondata.graphical_information);
% CONFIG.formattype = {'-dpng'}; 

% -------------------------------------------------------------------------
% M_MAP INFORMATION
% -------------------------------------------------------------------------
% Added by Thierry Reynaud 06/02/2020
% The M_MAP Library is introduced to replace Google Plot Maps which requered internet.
% The M_MAP path is defined here:
CONFIG.M_MAP_ACTIVE=1; % set to 0 ==> google_map plots
mmp = fullfile(CONFIG.LocodoxMainDir,'share','m_map1.4m',filesep);
CONFIG.M_MAPDir=mmp;
CONFIG.M_MAP_PLOT_BATHY=1;% PLOTTING ==> Reading ETOPO2 File

% -------------------------------------------------------------------------
% Ask online questions or not
% -------------------------------------------------------------------------
% By T. Reynaud and C. Kermabon 31/10/2023
CONFIG.onlineq=0;%Open dialog box for questions
%CONFIG.onlineq=1;%matlab prompt questions

% -------------------------------------------------------------------------
% climatological PSAL subtitution
% Added by T. Reynaud 14/12/2023
% -------------------------------------------------------------------------
% CONFIG.PSAL_REPLACE=false;
% CONFIG.PSAL_REPLACE_CLIM='ISAS';
% CONFIG.PSAL_REPLACE_DIR='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/ISAS/';
% CONFIG.PSAL_REPLACE_CLIM_file='isas17_PSAL.mat';
% %Config.PSAL_REPLACE_CLIM_file='isas15_PSAL.mat';
% CONFIG.PSAL_REPLACE_plot=false;
% CONFIG.PSAL_REPLACE_cycle_beg=6;
% CONFIG.PSAL_REPLACE_cycle_end=90;


disp(CONFIG)