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
%   $created: //2025.01.16 $author: Thierry Reynaud, LOPS
%   $Revision: version $Date: $author:


function [CONFIG_CLIM] = DOXY_PSAL_clim_config(cwmo)

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
    % Main directory of DOXY computation
    CONFIG_CLIM.FloatMainDir = '/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_FLOAT_DATA/DMQC_STEP1/coriolis/';
    % Directory for saving
    CONFIG_CLIM.resultsDir= '/Users/treynaud/IFREMER/MATLAB/LOCODOX/results/';

    % CONFIG_CLIM.PSAL_REPLACE_CLIM='ISAS';
    % CONFIG_CLIM.PSAL_REPLACE_DIR='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/ISAS/';
    % CONFIG_CLIM.PSAL_REPLACE_CLIM_file='isas17.mat';

    % CONFIG_CLIM.PSAL_REPLACE_CLIM='WOA';
    % CONFIG_CLIM.PSAL_REPLACE_DIR='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/WOA/';
    % CONFIG_CLIM.PSAL_REPLACE_CLIM_file='woa18.mat';

    CONFIG_CLIM.PSAL_REPLACE_CLIM='ARMOR-3D';
    CONFIG_CLIM.PSAL_REPLACE_DIR='/Users/treynaud/IFREMER/MATLAB/LOCODOX/Tools_treynaud/';
    CONFIG_CLIM.PSAL_REPLACE_CLIM_file='copernicus';

    CONFIG_CLIM.PSAL_REPLACE_plot=true;% plot figures for each individual cycle
    CONFIG_CLIM.PSAL_REPLACE_plot_close=false;
    CONFIG_CLIM.DOXY_PLOT_all=false; % plot figures all cumulatives cycles    
    CONFIG_CLIM.DOXY_PLOT_manual=true; % generate figures 1 by 1   
    CONFIG_CLIM.DOXY_QUEST_TOP=false; % flag top DOXY CLIM values    

else
    disp('User unknown ==> Define path in DOXY_PSAL_CLIM_config.m')
    exit;
end