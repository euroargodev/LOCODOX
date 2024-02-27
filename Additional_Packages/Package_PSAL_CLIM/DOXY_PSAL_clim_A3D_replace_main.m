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
%       v1.0 7/2/2024   Thierry Reynaud, LOPS
%                         initial
%

function [DATA, out_std,CLIM]=DOXY_PSAL_clim_A3D_replace_main(ptype,DATA,WOA,Work,CLIM)
%

% Replacing Float PSAL values by ARMOR-3D climatological values

%ptype='primary';% Select profile type
[DATA, out_std]=DOXY_PSAL_clim_A3D_replace_apply(ptype,DATA,Work,CLIM,Work.PSAL_REPLACE_cycle_beg,Work.PSAL_REPLACE_cycle_end);

% ptype='nearsurf';% Select profile type
% [DATA, outstd]=DOXY_PSAL_clim_replace_apply(ptype,DATA,Work,timeclim, depthclim, latclim, lonclim,psalclim,Work.PSAL_REPLACE_cycle_beg,Work.PSAL_REPLACE_cycle_end);


return
end
