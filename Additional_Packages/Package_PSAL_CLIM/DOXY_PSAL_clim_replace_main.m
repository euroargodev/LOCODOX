% DOXY_PSAL_clim_replace_main replaces float PSAL values by climatological ones.
%
% SYNTAX
% [DATA]=DOXY_PSAL_clim_replace(DATA,WOA,WORK)
%
% DESCRIPTION
% DOXY_PSAL_clim_replace replaces float PSAL values by climatological ones.
%
% INPUT
%     DATA (structure)    argo float "multiprofile" structure gathering the
%                         three kind of vertical sampling scheme. Subfields
%                         are argo (the data) and Dim (the dimension).
%                         Example :
%                             DATA =
%                                   primary: [1x1 struct]
%                                  nearsurf: [1x1 struct]
%                                 secondary: [1x1 struct]
%                                     other: [1x1 struct]
%
%                             DATA.primary =
%                                 argo: [1x1 struct]
%                                  Dim: [1x1 struct]
%
%                             DATA.primary.argo =
%                                   dimorder: 'C'
%                                  data_type: [1x1 struct]
%                                          ...
%                                       pres: [1x1 struct]
%                                       doxy: [1x1 struct]
%                                    doxy_qc: [1x1 struct]
%                              doxy_adjusted: [1x1 struct]
%                           doxy_adjusted_qc: [1x1 struct]
%                        doxy_adjusted_error: [1x1 struct]
%   WOA (struct)          Climatology WOA structure.
%                         Example:
%                          WOA =
%                                init: [1x1 struct]
%                              interp: [1x1 struct]
%                          WOA.init
%                                dimorder: 'C'
%                                latitude: [1x1 struct]
%                               longitude: [1x1 struct]
%                                   depth: [1x1 struct]
%                                    time: [1x1 struct]
%                                 doxywoa: [1x1 struct]
%                                           ...
%    Work (structure)     Doxy correction working structure, issued and
%                         computed from argo float data.
%                         Example:
%                             Work =
%                                   readme: [1x1 struct]
%                                     unit: 'mumol/kg'
%                                      wmo: 1901205
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
% HISTORY
%   $created: 12/12/2023 $author: Thierry Reynaud
%   $Revision: version $Date: $author:
%       v1.0 14/12/2023   Thierry Reynaud, LOPS
%                         initial
%

function [DATA, out_std,CLIM]=DOXY_PSAL_clim_replace_main(ptype,DATA,WOA,Work,CLIM)
%

% Read climatological fields if needed
%

if ~isfield(CLIM,'psalclim')

    if strcmp(Work.PSAL_REPLACE_CLIM,'ISAS')
        load(strcat(Work.PSAL_REPLACE_DIR,Work.PSAL_REPLACE_CLIM_file));
        %
        % on raboute le mois de decembre avec le mois de janvier et le mois de
        % janvier avec le mois de decembre.
        %    bid = NaN(length(timeclim)+2,1);


        bid = NaN(length(timeclim)+2,1);
        bid(2:end-1)=timeclim;
        bid(1) = timeclim(end)-365.25;
        bid(end) = timeclim(1)+365.25;
        timeclim = bid; % [14 x 1]
    else
        %
        % on raboute le mois de decembre avec le mois de janvier et le mois de
        % janvier avec le mois de decembre.
        %
        load(strcat(Work.PSAL_REPLACE_DIR,Work.PSAL_REPLACE_CLIM_file));

        timeclim=WOA.init.time.data'; % [12 x 1]
        bid = NaN(length(timeclim)+2,1);
        bid(2:end-1)=timeclim;
        bid(1) = timeclim(end)-365.25;
        bid(end) = timeclim(1)+365.25;
        timeclim = bid; % [14 x 1]

        lonclim=WOA.init.longitude.data'; % [360 x 1]
        latclim=WOA.init.latitude.data'; %  [180 x 1]
        %presclim=WOA.init.preswoa.data; %   [12x102x180x360]
        depthclim=WOA.init.depth.data';%    [102x1]
        psalclim=WOA.init.psal_woa.data; %   [12x102x180x360]
        tempclim=WOA.init.temp_woa.data; %   [12x102x180x360]
    end

    bid = NaN(length(timeclim),length(depthclim),length(latclim),length(lonclim));
    bid(1,:,:,:) = psalclim(end,:,:,:);
    bid(end,:,:,:) = psalclim(1,:,:,:);
    bid(2:end-1,:,:,:)=psalclim;
    psalclim = bid;

    bid = NaN(length(timeclim),length(depthclim),length(latclim),length(lonclim));
    bid(1,:,:,:) = tempclim(end,:,:,:);
    bid(end,:,:,:) = tempclim(1,:,:,:);
    bid(2:end-1,:,:,:)=tempclim;
    tempclim = bid;

    %psalclim = fillmissing(psalclim,'nearest'); % A commenter ou non selon les
    %cas --> 900 secondes pour isas
    %%%%%%%%%%%%%%%%%%%%%
    CLIM.timeclim=timeclim;
    CLIM.depthclim=depthclim;
    CLIM.latclim=latclim;
    CLIM.lonclim=lonclim;
    CLIM.psalclim=psalclim;
    CLIM.tempclim=tempclim;

    if exist('std_psalclim','var')
        CLIM.std_psalclim=std_psalclim;
    end

    if exist('std_tempclim','var')
        CLIM.std_tempclim=std_tempclim;
    end

end

% Replacing Float PSAL values by climatological values

%ptype='primary';% Select profile type
[DATA, out_std]=DOXY_PSAL_clim_replace_apply(ptype,DATA,Work,CLIM,Work.PSAL_REPLACE_cycle_beg,Work.PSAL_REPLACE_cycle_end);

% ptype='nearsurf';% Select profile type
% [DATA, outstd]=DOXY_PSAL_clim_replace_apply(ptype,DATA,Work,timeclim, depthclim, latclim, lonclim,psalclim,Work.PSAL_REPLACE_cycle_beg,Work.PSAL_REPLACE_cycle_end);


return
end
