% DOXY_interp_REF_2_argo interpolate REF to argo.
%
% SYNTAX
% [REF] = DOXY_interp_REF_2_argo(REF,iok,argoWork,icyc, PorD)
%
% DESCRIPTION
% DOXY_interp_REF_2_argo  interpolate REF data over argo horizontal and
% vertical grid. The vertical grid could be pressure (interpP) and over
% density (interpD).
%
% INPUT
%
%   REF (structure)      filled with the variable and its attributes found
%                        in the '.mat' insitu reference data file, and
%                        variables computed for the DOXY correction.
%                        Example:
%                        REF = 
%                              id: {438x1 cell}
%                           presi: [1x401 double]
%                             lon: [438x1 double]
%                             lat: [438x1 double]
%                             sta: [438x1 single]
%                            juld: [438x1 single]
%                            temp: [438x401 single]
%                          theta0: [438x401 single]
%                            psal: [438x401 single]
%                            sig0: [438x401 single]
%                            doxy: [438x401 single]
%                            pres: [438x401 double]
%                             sat: [438x401 single]
%                            psat: [438x401 single]
%                         doxy_CV: [361x1 single]
%                               ...
%
%   iok (boolean)         index of corresponding reference data for the
%                         argo float under correction
%
%   argowork (structure)  float intermediate working structure
%                         Example:
%                         argoWork = 
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%                                      psat: [1x1 struct]
%                                      ppox: [1x1 struct]
%
%   icyc (double)         index of the float under correction, in argoWork
%
%   PorD (cell string)    cell array of string, indicating if the level
%                         scale is the pressure or the density
%                         Example: {'pres','dens'}
%
% OUTPUT
%   REF (structure)      filled with the variable and its attributes found
%                        in the '.mat' insitu reference data file, and
%                        variables computed for the DOXY correction.
%                        Example:
%                        REF = 
%                              id: {438x1 cell}
%                           presi: [1x401 double]
%                             lon: [438x1 double]
%                             lat: [438x1 double]
%                             sta: [438x1 single]
%                            juld: [438x1 single]
%                            temp: [438x401 single]
%                          theta0: [438x401 single]
%                            psal: [438x401 single]
%                            sig0: [438x401 single]
%                            doxy: [438x401 single]
%                            pres: [438x401 double]
%                             sat: [438x401 single]
%                            psat: [438x401 single]
%                         doxy_CV: [361x1 single]
%                               ...
%
% CALL
%
% SEE ALSO
% DOXY_corr_main

% HISTORY
%   $created: 03/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%     v02 11.02.2022    Modified for Mutliple REF profiles use (T.Reynaud)
%                       see idxref

function [REF] = DOXY_interp_REF_2_argo(REF,iok,argoWork,icyc, PorD,presEff)


% =========================================================================
%% Initialisation
% =========================================================================

if any(~cellfun('isempty',strfind(PorD, 'pres')))
    letterPorD = 'P';
elseif any(~cellfun('isempty',strfind(PorD, 'dens')))
    letterPorD = 'D';
else
    letterPorD = 'P';
end

VarInRef = {'ppox','psat','doxy_CV','temp','psal',PorD{1}};
VarInArgo = {'ppox_adjusted','psat','doxy_adjusted','temp_adjusted','psal_adjusted',PorD{2}};

% ---------------------------------------------------------------------
% Compute ppox from REF data
% ---------------------------------------------------------------------

tmpDOXY=DOXY_convert(REF.doxy,'mumol/kg','mumol/L',REF.sig0);
if presEff == 1
    P = REF.pres; 
else
    P = 0; 
end
REF.ppox = O2ctoO2p(tmpDOXY,REF.temp,REF.psal,P);

% ---------------------------------------------------------------------
% Initialize the REF interpolated variables
% ---------------------------------------------------------------------
cmpl = sprintf('_interp%s',letterPorD);
for v = 1:length(VarInRef)
    REF.([VarInRef{v} cmpl]) = NaN(size(icyc,2),size(argoWork.(VarInArgo{v}).data,2));% Added T. Reynaud 07.02/2022
end

% =====================================================================
%% Vertical interpolation over pressure
% =====================================================================
idxref=REF.idxRef;%index in bddref - T. Reynaud 11.02.2022 - SORTED as bdd_REF_ARGO.txt file
for k=1:size(icyc,2)
    Z = REF.(PorD{1})(idxref(k),:);% from REF data base. TR 07.02.2022
    isok = ~isnan(Z);
    Z = Z(isok);
    if all(~isok)
        REF.(['ppox' cmpl])(iok,:) = NaN;
        REF.(['psat' cmpl])(iok,:) = NaN;
        REF.(['doxy_CV' cmpl])(iok,:) = NaN;
    else
        Zq = argoWork.(PorD{2}).data(icyc(:,k),:);
        for v=1:length(VarInRef)
            REF.([VarInRef{v} cmpl])(k,:) = interp1(Z',REF.(VarInRef{v})(idxref(k),isok),Zq,'linear');
        end
    end
end
