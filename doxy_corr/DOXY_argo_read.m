% DOXY_argo_read  read data
%
% SYNTAX
% [argo_primary, argo_nearsurf, argo_secondary,DimP, DimNS, DimS, argo, Dim] = ...
%    DOXY_argo_read(DataDir,wmo,replist,whichCorr,LogDir,metafiledir)
%
% DESCRIPTION
% DOXY_argo_read read data : 
% * Read data, taking into account the vertical sampling scheme
% * Extract the Primary, Secondary and NEar-Surface profile
% * Choose the pumped or unpumped profile depending on the correction type
% * Get the ascendants profile
%
% INPUT
%   DataDir (string)        full path of the data.
%   
%   wmo                     wmo of the float 
%   
%   replist (string)        wmo/profile/
%
%   whichCorr (string)      Indicate if the correction is made respect to
%                           climatology (WOA), InSitu reference database
%                           (REF) or In-air measurement (INAIR)
%                           Or if the drift correction is based on NCEP
%                           data
% 
%   LogDir (string)         the full path of the Log directory (defined in the
%                           configuration file)
%
%   ind_inair               ind_inair=1 if correction is based on in air method,
%                           ind_inair=0 if correction not based on in air
%                           method, but we need inair data to compute the drift
%
% OUTPUT
%     argo_primary (structure)    argo float "multiprofile" structure for
%     argo_nearsurf (structure)   the pumped data (Primary sampling) and
%     argo_secondary (structure)  the other vertical sampling scheme
%                                 (Near-surface and Secondary sampling).
%                                 Example :
%                                     argo_primary.temp.name = 'TEMPERATURE' name of the variable
%                                     argo_primary.temp.dim  = {'N_PROF','N_LEVEL'}
%                                     argo_primary.temp.data =  n_prof x n_level value of temperature
%                                     ...
%                                 If the structure.n_prof = 0, then the
%                                 argo float has no corresponding vertical
%                                 sampling scheme.
%                                 Example: argo_secondary.n_prof == 0 : No
%                                 secondary vertical sampling scheme in the
%                                 float.
%
%     DimP (structure)            argo float "multiprofile" structure for
%     DimNS(structure)            the pumped dimension (Primary sampling) and
%     DimS(structure)             the  other vertical sampling scheme
%                                 (Near-surface and Secondary sampling).
%                                 Example :
%                                 DimP = 
%                                         date_time: [1x1 struct]
%                                         string256: [1x1 struct]
%                                          string64: [1x1 struct]
%                                          string32: [1x1 struct]
%                                          ...
%
%     argo, Dim (structure)       float data structure (argo) and float Dim
%                                 structure (Dim) chosen for the Doxy
%                                 correction.
%
% CALL : DOXY_mono_to_multi, fillValue_2_nan, extract_profile_dim
%
% SEE ALSO
%   DOXY_corr_main

% HISTORYD
%   $created: 25/11/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%   23.03.2023 Kermabon Catherine and T. Reynaud ==> read properly the optode height 

function [argo_primary, argo_nearsurf, argo_secondary, argo_other, DimP, DimNS, DimS, DimO, argo, Dim] = ...
    DOXY_argo_read(DataDir,wmo,replist,whichCorr,LogDir,metafiledir)

% =========================================================================
%% Read data, taking into account the vertical sampling scheme
%   Primary sampling
%   Near-surface sampling
% =========================================================================

[argo_multiProfil, DimM] = DOXY_mono_to_multi(DataDir,wmo,replist,[],LogDir);
argo_multiProfil = fillValue_2_nan(argo_multiProfil,'FillValue_');
if argo_multiProfil.n_prof == 0
    fprintf('No data for the float %d\n',wmo);
end


% =========================================================================
%% Create vector to identifie which files contain DOXY %28/06/19 marine
% =========================================================================
cycle_nbr=unique(argo_multiProfil.cycle_number.data);
file_hereDoxy=zeros(size(argo_multiProfil.cycle_number.data)); % one if DOXY in the file containing the profile
for i=1:1:length(cycle_nbr)
   %Cycle number i with direction D
   prof_down=(argo_multiProfil.cycle_number.data==cycle_nbr(i)) & (argo_multiProfil.direction.data=='D');
   if any(argo_multiProfil.hereDoxy.data & prof_down)
       file_hereDoxy=file_hereDoxy | prof_down;
   end   
   
   %Cycle number i with direction A   
   prof_up=(argo_multiProfil.cycle_number.data==cycle_nbr(i)) & (argo_multiProfil.direction.data=='A');
   if any(argo_multiProfil.hereDoxy.data & prof_up)
       file_hereDoxy=file_hereDoxy | prof_up;
   end       

end

% =========================================================================
%% Extract the Primary, Secondary and NEar-Surface profile
%  The Primary is generally unique.
%  The Secondary and the Near-Surface profile could be multiple. The field
%  "hereDoxy" indicates which profiles that carries the doxy parameter has
%  been selected. 
%  All the other profiles are then stocked together in the structure named
%  "other".
% =========================================================================
VSS = {'Primary','Secondary','Near-surface'};
VSSstr = {'primary','secondary','nearsurf'};
VSSdim = {'P','S','NS'};
isVSSfinal = false(size(argo_multiProfil.vertical_sampling_scheme.data,1),1);
for iVSS = 1:length(VSS)
    isVSS = (findstr_tab(argo_multiProfil.vertical_sampling_scheme.data,VSS{iVSS})) & file_hereDoxy; %28/06/19 marine 
    if any(argo_multiProfil.hereDoxy.data) && strcmp(VSS{iVSS},'Primary') ~= 1
        % Primary is always unique in the vertical sampling scheme
        isVSS = isVSS & argo_multiProfil.hereDoxy.data ;        
        isVSSfinal = (isVSS | isVSSfinal);  
    else
        isVSSfinal = (isVSS | isVSSfinal);
    end

    eval(['[argo_' VSSstr{iVSS} ',Dim' VSSdim{iVSS} '] = extract_profile_dim(argo_multiProfil,DimM,''N_PROF'',isVSS);']);
    eval(['argo_' VSSstr{iVSS} '.VSS = [VSS{iVSS} '' sampling''];']);
end

% other case : nor Primary, nor secondary first choice carrying DOXY, nor
% near-surface first choice carrying DOXY
isother = ~(isVSSfinal) & argo_multiProfil.hereDoxy.data;
[argo_other, DimO] = extract_profile_dim(argo_multiProfil,DimM,'N_PROF',isother);
argo_other.VSS = 'Other remaining sampling';

% =========================================================================
%% Choose the pumped or unpumped profile depending on the correction type
% The profile pumped or unpumped used to compute the DOXY correction is
% choosen depending on the correction type :
%   - WOA and ATLN : Pumped
%   - In-Air : un-pumped 
% =========================================================================

% Depending on the correction type, choose the main profile from vertical
% sampling scheme
if ismember(whichCorr,{'WOA','REF'})
    if any(argo_primary.hereDoxy.data)
        argo = argo_primary;
        argo.VSS = 'Primary sampling';
        Dim = DimP;
    elseif any(argo_secondary.hereDoxy.data)
        argo = argo_secondary;
        argo.VSS = 'Secondary sampling';
        Dim = DimS;
    else
        warnstr = {sprintf('%d : no doxy data in both Primary and Secondary profile. ',wmo);...
            ' ';...
            sprintf('=> No correction could be applied for the %s correction.',whichCorr)};
        h = warndlg(warnstr,'DOXY_argo_read WARNING');
        uiwait(h);
        argo.n_prof = 0;
        Dim = [];
    end
elseif strcmp(whichCorr,'INAIR') || strcmp(whichCorr,'NCEP')
    if any(argo_nearsurf.hereDoxy.data)
            argo = argo_nearsurf;
            argo.VSS = 'Near-surface sampling';
            Dim = DimNS;
    else
        if strcmp(whichCorr,'INAIR')
            warnstr = {sprintf('%d : no doxy data in the Near-surface profile. ',wmo);...
                ' ';...
                sprintf('=> No correction could be applied for the %s correction.',whichCorr)};
            h = warndlg(warnstr,'DOXY_argo_read WARNING');
            uiwait(h);
            argo.n_prof = 0;
            Dim = [];
            argo_primary.no_inair_data=1;
            return
        elseif strcmp(whichCorr,'NCEP')  
            warnstr = {sprintf('%d : no doxy data in the Near-surface profile. ',wmo);...
                ' ';...
                sprintf('=> No drift correction based on NCEP could be applied for the %s correction.',whichCorr)};
            h = warndlg(warnstr,'DOXY_argo_read WARNING');
            uiwait(h);
            argo.n_prof = 0;
            Dim = [];
            argo_primary.no_inair_data=1;
            return        
        end 
    end
    
end

% =========================================================================
%% Read metadata file
% =========================================================================
metafile=[DataDir metafiledir num2str(wmo) '_meta.nc']
metaval=NCR_file(metafile);
yyyy=str2num(metaval.launch_date.data(1:4));
mm=str2num(metaval.launch_date.data(5:6));
dd=str2num(metaval.launch_date.data(7:8));
hh=str2num(metaval.launch_date.data(9:10));
mn=str2num(metaval.launch_date.data(11:12));
ss=str2num(metaval.launch_date.data(13:14));
argo.launchdate=datenum(yyyy,mm,dd,hh,mn,ss)-datenum(1950,1,1,0,0,0);

% Added by T. Reynaud and K. Kermabon : 23.03.2023
%config_launch = ncread(metafile,'LAUNCH_CONFIG_PARAMETER_NAME');
%config_launch = config_launch';
%config_launch = cellstr(config_launch);
%isinfo_optode = find(strcmp(config_launch,'CONFIG_OptodeVerticalPressureOffset_dbar'));

argo.heightoptode = 0.2; % 20 cm assumed as default height of optode above sea surface, if not given in 'CONFIG_OptodeVerticalPressureOffset_dbar' in meta file
if isfield(metaval,'OptodeVerticalPressureOffset_dbar')
    argo.heightoptode=abs(-metaval.OptodeVerticalPressureOffset_dbar);
%
% 27/03/2023 : C. Kermabon - T. Reynaud : Reading OptodeVerticalPressureOffset.
%
elseif isfield(metaval,'launch_config_parameter_name') 
    if ~isempty(find(strcmp(cellstr(metaval.launch_config_parameter_name.data),'CONFIG_OptodeVerticalPressureOffset_dbar')))
     isinfo_optode = find(strcmp(cellstr(metaval.launch_config_parameter_name.data),'CONFIG_OptodeVerticalPressureOffset_dbar'));
     argo.heightoptode = metaval.launch_config_parameter_value.data(isinfo_optode); 
    end
end
% End 27.03.2023

if argo.heightoptode<0, argo.heightoptode=-argo.heightoptode; disp(['  Assume optode to be ' num2str(argo.heightoptode) ' m above the pressure sensor, not below'])
elseif argo.heightoptode==0, argo.heightoptode=1e-4; % to be able to use same equation below; equal to sea-surface water vapor pressure
end 

% =========================================================================
%% Get the ascendants profile
% Only ascendant ones are used to compute the Doxy correction. The
% correction is then applied to each sort of profile.
% =========================================================================
isAprofile = strcmp(cellstr(argo.direction.data),'A');
[argo,Dim] = extract_profile_dim(argo,Dim,'N_PROF',isAprofile);

