% DOXY_config define the configuration parameters usefull for the oxygen
% correction LOCODOX processing chain
%
% SYNTAX
% [] = DOXY_config
%
% DESCRIPTION
% DOXY_config define the configuration parameters usefull for the oxygen
% correction LOCODOX processing chain
%
% INPUT
%   no inputs
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
%
% SEE ALSO
%
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
%       v4 26/01/2017   Anne Piron, Emilie Brion, ALTRAN OUEST
%                       argo - ameliorations 2017 phase 1
%       v4.1 19/09/2017 Emilie Brion, ALTRAN OUEST
%                       add the "unselect" section
%       v4.2 31/07/2018 Emilie Brion, ALTRAN OUEST
%                       - add format type for figure saving
%                       - add the pressure effect section : CONFIG.presEff
%                       - tune the measurement_code for in-air and in-water
%                       samples in the trajectory
%            20/11/2018 Marine GALLIAN, ALTRAN OUEST
%                       New format for reference data : we use now a unique
%                       txt file while we were using two .mat files before.
%
%       v3.3 14/08/2018 Marine GALLIAN, ALTRAN OUEST
%       v3.4 09.04.2020 Thierry Reynaud -> WOA+M_MAP+Share rearrangment
%       v3.4 10.04.2020 Thierry Reynaud -> add CONFIG.ncepGetYears = str2double(datestr(now,'YYYY'));% To be downloaded
%       v4.0 18.01.2021 Thierry Reynaud -> add path selection upon username
%       v5.0 12.04.2021 Thierry Reynaud --> Time Drift Depth parameters
%       modified

function [CONFIG] = locodox_config

% -------------------------------------------------------------------------
% DIRECTORIES
% -------------------------------------------------------------------------

if ismac
    % Code to run on Mac platform
    username=deblank(getenv('USER'));% For Mac and Linux
elseif isunix
    % Code to run on Linux platform
    username=deblank(getenv('USER'));% For Mac and Linux
elseif ispc
    % Code to run on Windows platform
    username=deblank(getenv('USERNAME'));% For windows   
else
    disp('Platform not supported')
    exit;
end


if strfind(username,'treynaud')
    % Main directory of LOCODOX
    CONFIG.LocodoxMainDir = '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX/';
    % Directory of the Argo NetCDF data
    CONFIG.DataDir = '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_FLOAT_DATA/DMQC_PSAL/coriolis/';
    %CONFIG.DataDir = '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_FLOAT_DATA/CORIOLIS/coriolis/';% Pour CK 2023.04.06
    % External data directories : NCEP, WOA and TOPOGRAPHY
    CONFIG.ExtDataDir = '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_DATA/';
    % added by Thierry Reynaud 23.04.2020
    CONFIG.LOPSDataDir = '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/';
    % added by Thierry Reynaud 06.02.2023
    CONFIG.NCEPDataDir= '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_DATA/NCEP/';
    % added by Thierry Reynaud 23.04.2020
    CONFIG.WOADataDir= '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_DATA/WOA/';
    % added by Thierry Reynaud 23.04.2020
    CONFIG.BathyDataDir= '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_DATA/TOPOGRAPHY/';
    % Directory for saving
    CONFIG.resultsDir= '/Users/treynaud/IFREMER/MATLAB/LOCODOX/results/';
    
elseif strfind(username,'vthierry')
    % Main directory of LOCODOX
    CONFIG.LocodoxMainDir = '/Users/vthierry/matlab/GITHUB_LOCODOX/LOCODOX/';
    % Directory of the Argo NetCDF data
    %CONFIG.DataDir = '/Volumes/Virginie_data/ARGOGDAC/dac/coriolis/';
    CONFIG.DataDir = '/Users/vthierry/Desktop/';
    % External data directories : NCEP, WOA and TOPOGRAPHY
    CONFIG.ExtDataDir = '/Users/vthierry/matlab/GITHUB_LOCODOX/LOCODOX/LOCODOX_EXTERNAL_DATA/';
    % added by Thierry Reynaud 23.04.2020
    CONFIG.LOPSDataDir = '/Users/vthierry/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/';
    % added by Thierry Reynaud 06.02.2023
    CONFIG.NCEPDataDir= '/Users/vthierry/matlab/GITHUB_LOCODOX/LOCODOX/LOCODOX_EXTERNAL_DATA/NCEP/';
    % added by Thierry Reynaud 23.04.2020
    CONFIG.WOADataDir= '/Users/vthierry/matlab/GITHUB_LOCODOX/LOCODOX/LOCODOX_EXTERNAL_DATA/WOA/';
    % added by Thierry Reynaud 23.04.2020
    CONFIG.BathyDataDir= '/Users/vthierry/matlab/GITHUB_LOCODOX/LOCODOX/LOCODOX_EXTERNAL_DATA/TOPOGRAPHY/';
    % Directory for saving
    CONFIG.resultsDir= '/Users/vthierry/matlab/GITHUB_LOCODOX/LOCODOX/results/';
else
    disp('User unknown ==> Define path in locodox_config.m')
    exit;
end

% Add the useful paths
mytools = fullfile(CONFIG.LocodoxMainDir,'share','MyTools',filesep);
addpath(mytools);
hydcalDir = fullfile(CONFIG.LocodoxMainDir,'share','hydcal',filesep);
seawaterDir = fullfile(CONFIG.LocodoxMainDir,'share','seawater','seawater_330_its90_lpo',filesep);
addpath(hydcalDir);
addpath(seawaterDir);
addpath(fullfile(CONFIG.LocodoxMainDir,'share',filesep));
addpath(fullfile(CONFIG.LocodoxMainDir,'doxy_corr',filesep));
addpath(fullfile(CONFIG.LOPSDataDir,'data_input',filesep));% TR 06.02.2023
% woa_pth added by Thierry Reynaud 09.04.2020
%woa_pth=fullfile(CONFIG.ExtDataDir,'WOA',filesep);
%addpath(woa_pth);

% -------------------------------------------------------------------------
% FILES
% -------------------------------------------------------------------------
% WOA climatology file
WOA_YEAR='2018';% Release Year
CONFIG.WOAfile = [CONFIG.WOADataDir,'WOA2018_DECAV_monthly_5500_1deg.nc'];

%WOA_YEAR='2009';% Release Year
%CONFIG.WOAfile = [CONFIG.WOADataDir,'WOA2009_monthly_5500_1deg.nc'];


% Reference In-Situ DataBase
%CONFIG.bddFile = 'bddo2ref_avecov18_temp.mat';
%CONFIG.bddFile='bddo2ref_ov18_temp.mat'; % Only for floats no 6901763 6902818 6902881 6902882 6901601
%CONFIG.bddFile='bddo2ref_ov18_temp.mat'; % Only for floats no 6901763 6902818 6902881 6902882 6901601 6902800
%CONFIG.bddFile='bddo2ref.mat'; % rr15+rr17
%CONFIG.bddFile='bddo2ref_vracape.mat'; % 1900943 6900629
CONFIG.bddFile='bddo2ref_all_TR_2024.mat';


% Reference data associated to a wmo
CONFIG.RefArgoFile='bdd_REF_ARGO.txt';%

% Mask file and the variable to be read in the mask file
CONFIG.maskFile = 'landsea_masks.cdf';
CONFIG.varMask = 'landsea';

% Logo file
CONFIG.logo = fullfile(CONFIG.LocodoxMainDir,'locodox_logo.jpg');

% NCEP files and FTP address
% ncepDoUpdate : read the ftp website to update your NCEP 
% ncepFtp : the ftp website of NCEP
% ncepFtpDir : the path where to find the NCEP data in the ftp website
% ncepFtpSubDir : the sub directory where to find the NCEP data in the ftp website
% ncepFiles : the NCEP files to be read
% ncepYears : read the NCEP data for the years specified
CONFIG.ncepDoUpdate = 1;  
CONFIG.ncepFtp = 'ftp.cdc.noaa.gov';
CONFIG.ncepFtpDir = 'Datasets/ncep.reanalysis/';
CONFIG.ncepFtpSubDir = {'surface','surface','surface'};
CONFIG.ncepFiles = {'slp','air.sig995','rhum.sig995'};
CONFIG.ncepYears = [2014:2023];% To be used
CONFIG.ncepGetYears = str2double(datestr(now,'YYYY'));% To be downloaded
%CONFIG.ncepGetYears = [2021:2022];% To be downloaded


% -------------------------------------------------------------------------
% FOR CORRECTION
% -------------------------------------------------------------------------

% Reference data unit : precise the unit of your reference data (most of the
% time in-situ CTD)
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Possible options:
%         mL/L, milliliter per liter, milliliter/L
%         mumol/m3, micromole per m3, mmol/m3, micromole/m3
%         mumol/L, micromole per liter, mmol/L, micromole/L
%         mg/L, milligram per liter, milligram/L
%         mumol/kg, micromole per kilo, mmol/kg, micromole per kilogram, micromole/kg 
CONFIG.refUnit = 'mumol/kg';

% Conversion DOXY/PSAT/PPOX : take into account pressure effect or not
CONFIG.presEff = 0;   % 0/1 option unactivated/activated

% Carry over parameter. If isokC is set to 0, we suppose that in air oxygen measurements are not biased by splash of water; C is set to 1 otherwise
% See "Oxygen Optode Sensors : Principle Characterization, Calibration, and
% Application in the Ocean", Henry Bittig and al (2018)
CONFIG.isokC=1;

%CONFIG.isokC=0;% Pour CK 2023.04.06
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
CONFIG.DM_pres = 1;
CONFIG.DM_temp = 1;
CONFIG.DM_psal = 1; 
 
% CONFIG.DM_pres = 0;% Pour CK
% CONFIG.DM_temp = 0;% Pour CK
% CONFIG.DM_psal = 0;% Pour CK

% QC selection: 
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Only the data whose QC equal to these values are 
% retained for the calculation of the correction (for DOXY, PRES, TEMP and PSAL).
% Format: Example: CONFIG.QC_x = [1 2 3 4] (possible values: from 0 to 9)
CONFIG.QC_O = [1 2 3]; % DOXY
CONFIG.QC_P = [1 2]; % PRES
CONFIG.QC_T = [1 2]; % TEMP
CONFIG.QC_S = [1 2 3]; % PSAL

% For FSD floats: Hybrid DM a,d RT
% Comment otherwise
% CONFIG.RT_psal_cycle=158;% Pour 6902806 FSD use RT PSAL data from specified cycle ==> DM_PSAL=1
% CONFIG.QC_S = [1 2 3 4]; % PSAL FSD 6902806 ==> DP_psal=1

% Time Drift correction
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
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
%Compute inair drift if possible for WOA and REF correction 
CONFIG.ok_inair_drift=1;

% Added by T. Reynaud for Piece wise linear fitting for Time Drift 
% 26/04/2024
CONFIG.drift_PWLF_N = 2;% Number of linear segments

%Minimum depth for calculating drift (default value = 1500m)
%CONFIG.min_drift_depth=800;
%Name modified by TR 12/04/2021
% In deeper water
%CONFIG.min_drift_depth_deep=2100;
%CONFIG.max_drift_depth_deep=2500;
%CONFIG.step_drift_depth_deep=100;

CONFIG.min_drift_depth_deep=1700;
%CONFIG.max_drift_depth_deep=1800;% Normal
CONFIG.max_drift_depth_deep=4100;
CONFIG.step_drift_depth_deep=100;

%CONFIG.min_drift_depth_deep=2500;
%CONFIG.max_drift_depth_deep=4000;
%CONFIG.step_drift_depth_deep=100;

%CONFIG.min_drift_depth_deep=450;
%CONFIG.max_drift_depth_deep=550;
%CONFIG.step_drift_depth_deep=50;

%CONFIG.min_drift_depth_deep=800;
%CONFIG.max_drift_depth_deep=1000;
%CONFIG.step_drift_depth_deep=100;

% Near surface
% Introduced by TR 12/04/2021
CONFIG.min_drift_depth_surf= 0;
%CONFIG.max_drift_depth_surf=25;%Normal
CONFIG.max_drift_depth_surf=30;
CONFIG.step_drift_depth_surf=5;

% WOA or REF correction
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Linear correction R2 coefficient : treshold under which LOCODOX suggests
% to apply Constant correction
CONFIG.R2threshold = 0.80;

% In-air correction
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% trajectory data : measurement codes for In-Air and Near-surface samples
% Measurement Code which carries the InAir measurement :
%  - CONFIG.inAirMC : MC for the in-air samples, part of the surface
%  sequence (double, array)
%  - CONFIG.inWaterMC : MC for the in-water samples, part of the surface
%  sequence (or the profile sequence).(double, array)
CONFIG.inWaterMC = [690 710];
CONFIG.inAirMC = [699 711 799];

% maximum depth for TEMP and PSAL selection as "surface temperature" and
% "surface salinity"
CONFIG.inAirMaxPresForTS = 20;

% Error on adjusted field
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% abs (absolute): adjusted_error is a fixed value
% rel (relative): adjusted error is a % value:  doxy_adjusted + adjusted_error_rel*doxy_adjusted/100
CONFIG.adjusted_error_abs = 0;
CONFIG.adjusted_error_rel = 3;  

% -------------------------------------------------------------------------
% FOR NetCDF writing and FOR HISTORY
% -------------------------------------------------------------------------
CONFIG.history_software = 'LOCODOX';
CONFIG.history_reference = ['LOPS2020_WOA',WOA_YEAR];
CONFIG.history_software_release = '4.4';
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
CONFIG.savePlot = 1;
CONFIG.savePlotFig = 1;%Save figures as matlab .fig files
CONFIG.fontsize = 14;
CONFIG.resolution = 100;
CONFIG.formattype = {'-dpng'}; 

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


