% DOXY_NCEP_colocalize colocalizes the NCEP data over the argoTrajWork time/lat/lon
% regridded in a 6h daily grid.
%
% SYNTAX
% [NCEP] = DOXY_NCEP_colocalize(NCEP,argoTrajWork,CONFIG)
%
% DESCRIPTION
% DOXY_NCEP_colocalize colocalizes the NCEP data over the argoTrajWork time/lat/lon
% regridded in a 6h daily grid.
%
% INPUT
%   NCEP (structure)     filled with the variable and its attributes found
%                        in the NCEP climatology files
%                        Example: NCEP
%                              slp: [1x1 struct]
%                              air: [1x1 struct]
%                             rhum: [1x1 struct]
%
%                        with NCEP.slp = 
%                         dimorder: 'C'
%                              lat: [1x1 struct]
%                              lon: [1x1 struct]
%                              slp: [1x1 struct]
%                           recdim: 'time'
%                              obj: 'ObsInSitu'
%                        fillisnan: 1
%                             juld: [1x1 struct]
%                             time: 23884
%                     firstdimname: 'time'
%
%                       and NCEP.slp.slp =
%                            name: 'slp'
%                             dim: {'time'  'lat'  'lon'}
%                            data: [23884x73x144 single]
%                       long_name: '4xDaily Sea Level Pressure'
%                           units: 'mBar'
%                       precision: 0
%                                   ..........
%
%  argoTrajWork (structure)     float data structure (argoTrajWork) choose for the Doxy
%                               correction.
%                       Example: argoTrajWork
%                          float_serial_no: [1x1 struct]
%                            wmo_inst_type: [1x1 struct]
%                                     juld: [1x1 struct]
%                                 latitude: [1x1 struct]
%                                longitude: [1x1 struct]
%                                   ..........
%
% CONFIG (struct)       Configuration structure with data path,
%                       operator choices, ...
%                             CONFIG = 
%                        DataDir: '/home/oo26/coriolis/co05/co0508/dac/coriolis/'
%                        NCEPDataDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/data_input/'
%                        saveDataDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/results/data/'
%                        savePlotDir: '/home4/homedir4/perso/mgallian/LOCODOX2/LOCODOX/results/plots/'
%                                   ..........
%                             
% OUTPUT
%   NCEP (structure)     the field coloc is filled with colocalized
%                        variables over the regridded argoTrajWork data (6h daily)
%                        Example: NCEP.coloc
%                                     juld: [1522x1 double]
%                                      slp: [1522x1 double]
%                                      air: [1522x1 double]
%                                     rhum: [1522x1 double]
%
% CALL :
%
% SEE ALSO
%   DOXY_corr_main, DOXY_NCEP_read

% HISTORY
%   $created: 09/05/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       22/03/2019      Marine GALLIAN, Altran Ouest
%                       Update function in order to manage new MC
%                       measurement code for in air data
%
% 27/04:2023 : C. Kermabon.
%                       Utilisation interpn + prise en compte ERA5.
%

function [NCEP] = DOXY_NCEP_colocalize(NCEP,argoTrajWork,CONFIG)

% Take inAir argoTrajWork time, lat and lon --MG
 for i = 1:length(argoTrajWork.ppox_doxy_adjusted.data)
     isInAir = ismember(argoTrajWork.measurement_code.data{i},CONFIG.inAirMC);
     if ~isempty(isInAir)
         argoTrajWork.juld.data{i} = argoTrajWork.juld.data{i}(isInAir);
         argoTrajWork.latitude.data{i} = argoTrajWork.latitude.data{i}(isInAir)  ;      
         argoTrajWork.longitude.data{i} = argoTrajWork.longitude.data{i}(isInAir);
     end
 end

% =========================================================================
% Interpolate argoTrajWork time, lat and lon over a 6h daily grid
% =========================================================================
if iscell(argoTrajWork.juld.data)
    juld_argo = cellfun(@(x) x+datenum(strrep(argoTrajWork.juld.units,...
            'days since',''),'yyyy-mm-dd HH:MM:SS'),...
            argoTrajWork.juld.data,'UniformOutput',false);
    juld_argo = cell2mat(juld_argo);
else
    juld_argo = argoTrajWork.juld.data + ...
        datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');
end
 

% =========================================================================
% Interpolate data to argoTrajWork profile date
% Finalize data
% =========================================================================
varNames = fieldnames(NCEP.init);
varNames = varNames(~ismember(varNames,'nodata'));
NCEP.coloc.juld = cell2mat(argoTrajWork.juld.data);
juld_ncep =  NCEP.init.(varNames{1}).juld.data + ...
        datenum(strrep(NCEP.init.(varNames{1}).juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');

if max(juld_ncep)<max(juld_argo)
        fprintf('\n\nAttention : INAIR DATA are loaded until %s.\nARGO are loaded until %s.\nUpdate the INAIR Data.',datestr(max(juld_ncep)),datestr(max(juld_argo)));
end
if min(juld_ncep)>min(juld_argo)
        fprintf('\n\nAttention : INAIR DATA are loaded from %s.\nARGO are loaded from %s.\nUpdate INAIR DATA.\n',datestr(min(juld_ncep)),datestr(min(juld_argo)));
end
%
% Verification que certaines annees, exceptee la derniere, ne sont pas
% tronquees.
%
diff_juld_ncep = diff(juld_ncep);
unic_diff = uniquetol(diff_juld_ncep,1/48);
if length(unic_diff)~=1
    fprintf('Attention : Error on INAIR DATA between\nThe difference between data are not unic.\n');
end
for i = 1:length(varNames)    
    easyVar = varNames{i};
    % Longitude between 0 and 360.
    %
    lon_ncep = double(NCEP.init.(easyVar).lon.data);
    lon_ncep(lon_ncep<0) = lon_ncep(lon_ncep<0)+360;
    %[lon_sort,IndLon]=sort(lon_ncep);
    lon_argo = cell2mat(argoTrajWork.longitude.data);
    lon_argo(lon_argo<0) = lon_argo(lon_argo<0) + 360;
    juld_ncep =  NCEP.init.(easyVar).juld.data + ...
        datenum(strrep(NCEP.init.(easyVar).juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS');
    
%     NCEP.coloc.(easyVar) = interp1(juld_ncep,dataNcepOnArgo.(easyVar).data,...
%         NCEP.coloc.juld + datenum(strrep(argoTrajWork.juld.units,'days since',''),'yyyy-mm-dd HH:MM:SS'));
%
% NCEP ne va que jusqu'à la longitude 357.5. 
% On boucle le tableau sur 360 degre.
%
% EX : 
% si la longitude va de 5° à 355°, ce qui suit permet d'avoir des longitudes de 5° à
% 365°. Toutes les longitudes sont ainsi couvertes.
% Les longitudes ARGO superieures à 355° seront correctement
% interpolées. Les longitudes ARGO inferieures à 5°, il faut alors leur
% rajouter 360° pour qu'elles puissent etre interpolees correctement.
%
bid = NaN(length(lon_ncep)+1,1);
bid(1:end-1) = lon_ncep;
bid(end) = min(lon_ncep) + 360; % min(lon_ncep) proche de 0 car lon entre 0 et 360.
lon_ncep = bid;
[lon_sort,IndLon]=sort(lon_ncep);
bid = NaN(size(NCEP.init.(easyVar).(easyVar).data,1),size(NCEP.init.(easyVar).(easyVar).data,2),size(NCEP.init.(easyVar).(easyVar).data,3)+1);
bid(:,:,1:end-1) = NCEP.init.(easyVar).(easyVar).data(:,:,1:end);
ind_min_lon = find(lon_ncep==min(lon_ncep));
bid(:,:,end) = NCEP.init.(easyVar).(easyVar).data(:,:,ind_min_lon);

%
% Prise en compte du cas où la longitude ARGO est inferieure à la longitude
% NCEP/ERA5.
%
isbad = find(lon_argo<min(lon_ncep));
lon_argo(isbad) = lon_argo(isbad) + 360;

     NCEP.coloc.(easyVar) =interpn(juld_ncep,double(NCEP.init.(easyVar).lat.data),lon_ncep(IndLon),bid(:,:,IndLon),juld_argo,cell2mat(argoTrajWork.latitude.data),lon_argo);
end