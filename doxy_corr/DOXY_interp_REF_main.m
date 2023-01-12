% DOXY_interp_REF_main interpolate REF data over argo horizontal and
% vertical grid.
%
% SYNTAX
% [REF, Work] = DOXY_interp_REF_main(bddFile, REF_ARGO, argo,argoWork, Work)
%
% DESCRIPTION
% DOXY_interp_REF_main interpolate Reference In Situ data over argo
% horizontal and vertical grid. The vertical grid could be pressure
% (interpP) and over density (interpD). The working structure is completed.
%
% INPUT
%
%   bddFile (structure)     file of data file of insitu Reference data (.mat).
%                           Once loaded, it's a structure.Example:
%                           REF =
%                                      id: {438x1 cell}
%                                   presi: [1x401 double]
%                                     lon: [438x1 double]
%                                     lat: [438x1 double]
%                                     sta: [438x1 single]
%                                    juld: [438x1 single]
%                                    temp: [438x401 single]
%                                  theta0: [438x401 single]
%                                    psal: [438x401 single]
%                                    sig0: [438x401 single]
%
%   REF_ARGO (structure)    List the WMO of argo float with DOXY sensor,
%                           and useful information : sensor, DAC,
%                           directory, reference profile. Example:
%                            REF_ARGO = 
%                                 cycle : [1x1 double]
%                                 refId : [1x1 string]
%                                 wmo : [1x1 double]
%
%   argo (structure)      Argo float structure (get directly from data).
%                         Example :
%                         argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%                               argo.pres =
%                                       name: 'DOXY'
%                                        dim: {'N_PROF'  'N_LEVELS'}
%                                       data: [85x120 single]
%                                  long_name: 'DISSOLVED OXYGEN'
%                              standard_name: 'moles_of_oxygen_per_unit_mass_in_sea_water'
%                                 FillValue_: 99999
%                                      units: 'micromole/kg'
%                                  valid_min: 0
%                                  valid_max: 650
%                                   C_format: '%9.3f'
%                             FORTRAN_format: 'F9.3'
%                                 resolution: 0.0010
%                                       type: 5
%
%   argowork (structure)  float working structure issued from the argo data
%                         Example:
%                         argoWork =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                   an_dens: [1x1 struct]
%                                   density: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%                                       sat: [1x1 struct]
%                                      psat: [1x1 struct]
%                             doxy_adjusted: [1x1 struct]
%
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA)
%                           Example:
%                        Work =
%                                     readme: [1x1 struct]
%                                 R2treshold: '0.80'
%                         adjusted_error_abs: '0'
%                         adjusted_error_rel: '5'
%                                       unit: Work.refUnit
%                                    dirPlot: '/home1/homedir4/perso/ebrion/NAOS/NAOS_2016(1)/plots2/REF'
%                                   makePlot: 1
%                                   savePlot: 1
%                                   fontsize: 8
%                                    history: [1x1 struct]
%                                        wmo: 5902302
%                                     sensor: 'Optode'
%                                  whichCorr: 'REF'
%                                   DOXY_RAW: [175x110 single]
%                                      timar: [175x1 datetime]
%                                      datat: [175x1 double]
%                                    DENSITY: [175x110 single]
%                                      DEPTH: [175x110 double]
%                                    DOXY_QC: [175x110 char]
%                                       DENS: [175x110 single]
%                                       PSAT: [175x110 single]
%                                       DOXY: [175x110 single]
%                                       PSAL: [175x110 single]
%                                       TEMP: [175x110 single]
%                                       PRES: [175x110 single]
%                                 PresOrDens: 'Pressure'
%                                      DRIFT: -3.1841
%                                    DODRIFT: 0
%
%
% OUTPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA),
%                        completed with the annual drift.
%                        Example:
%                        Work =
%                                     readme: [1x1 struct]
%                                 R2treshold: '0.80'
%                         adjusted_error_abs: '0'
%                         adjusted_error_rel: '5'
%                                       unit: Work.refUnit
%                                    dirPlot: '/home1/homedir4/perso/ebrion/NAOS/NAOS_2016(1)/plots2/REF'
%                                   makePlot: 1
%                                   savePlot: 1
%                                   fontsize: 8
%                                    history: [1x1 struct]
%                                        wmo: 5902302
%                                  whichCorr: 'REF'
%                                   DOXY_RAW: [175x110 single]
%                                      timar: [175x1 datetime]
%                                      datat: [175x1 double]
%                                    DENSITY: [175x110 single]
%                                      DEPTH: [175x110 double]
%                                    DOXY_QC: [175x110 char]
%                                       DENS: [175x110 single]
%                                       PSAT: [175x110 single]
%                                       DOXY: [175x110 single]
%                                       PSAL: [175x110 single]
%                                       TEMP: [175x110 single]
%                                       PRES: [175x110 single]
%                                 PresOrDens: 'Pressure'
%                                      DRIFT: -3.1841
%                                    DODRIFT: 0
%                           DOXY_REF_interpP: [1x110 double]
%                           PSAT_REF_interpP: [1x110 double]
%                           PSAL_REF_interpP: [1x110 double]
%                           TEMP_REF_interpP: [1x110 double]
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
% CALL : 
%   DOXY_drift, DOXY_interp_REF_2_argo, DOXY_PLOT_interpolation
%
% SEE ALSO
%   DOXY_corr_main, DOXY_ref_corr

% HISTORY
%   $created: 07/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v1.1 25/05/2016   Emilie Brion, ALTRAN OUEST
%                         improved
%       v1.2 06/09/2017   Emilie Brion, ALTRAN OUEST
%                         - reference cycle index properly managed
%                         - update description
%       v1.3 20/09/2018   Emilie Brion, ALTRAN OUEST
%                         - REF PSAT conversion properly manage, from
%                         Gordon and Garcia formula.
%            04/03/2019   Marine GALLIAN, ALTRAN OUEST
%                         - Reference cycle saved if IN AIR correction
%       v1.4 11.02.2022   - Modified for Mutliple REF profile use (T.Reynaud)


function [REF, Work] = DOXY_interp_REF_main(bddFile, REF_ARGO, argo, ...
    argoWork, Work)

% =========================================================================
%% Initialize
% =========================================================================
% Reference data
load(bddFile); 

% Reference PSAT computation
tmpDoxy = DOXY_convert(REF.doxy,Work.refUnit,'mumol/L',REF.sig0);
if Work.presEff == 1, P = REF.pres; else, P = 0; end
REF.psat = O2ctoO2s(tmpDoxy,REF.temp,REF.psal,P);
% tmpDoxy = DOXY_convert(REF.doxy,Work.refUnit,'mL/L',REF.sig0);
% REF.sat = (tmpDoxy./REF.psat)*100;
REF.doxy_CV = DOXY_convert(REF.doxy,Work.refUnit,Work.unit,REF.sig0);

% density
REF.density = REF.sig0 + 1000;


% =========================================================================
%% display the date of the reference data
% =========================================================================

% useful index

%Added by T. Reynaud 07/02/2022
iokRef=logical(zeros(size(REF.id)));
idxRef=zeros(size(REF_ARGO.refId,2),1);% Ref index order Reynaud 10.02.2022
for k=1:size(REF_ARGO.refId,2)
    iokRef(strcmp(REF.id,REF_ARGO.refId{k})==1)=1;
    idxRef(k)=find(strcmp(REF.id,REF_ARGO.refId{k})==1);
end
REF_ARGO.idxRef=idxRef;% Added by T. Reynaud 10.02.2022
REF.idxRef=idxRef;% Added by T. Reynaud 11.02.2022

iokCyc=logical(zeros(size(argo.cycle_number.data,1),size(REF_ARGO.refId,2)));%Added by T. Reynaud 07/02/2022
for k=1:size(REF_ARGO.refId,2)
    iokCyc(:,k) = argo.cycle_number.data == REF_ARGO.cycle(k) & argo.direction.data == 'A';% Comumn following the order in the bdd_REF_ARGO.txt file
end
% End of modifications T. Reynaud 07.02.2022

if isfield(argo,'is_main_profile')
    % Reference date
    refDate = datestr(double(REF.juld(idxRef)) + datenum(1950,1,1,0,0,0));% T.Reynaud 10.02.2022

    for k=1:size(REF_ARGO.refId,2)
        tmp = datestr(argo.juld.data(iokCyc(:,k)) + datenum(1950,1,1,0,0,0));
        if k==1
            argoDate=tmp;
        else
            argoDate=[argoDate;tmp];
        end
        fprintf('\t argo date %s, reference date %s\n', argoDate(k,:), refDate(k,:));

        % Reference DOXY range for the profile
        range = [nanmin(REF.doxy(idxRef(k),:)) nanmax(REF.doxy(idxRef(k),:))];% T. Reynaud 10.02.2022
        fprintf('\t the doxy reference data is : [%2.1f - %2.1f] \n', range(1), range(2))

    end

end
% =========================================================================
%% Interpolation : interpolate REF over argo
% =========================================================================
% The interpolation is made over pressure (interpP) and over density
% (interpD)
PorD = {{'pres','pres_adjusted'},{'density','density'}};
cmpl = {'_interpP','_interpD'};

[REF] = DOXY_interp_REF_2_argo(REF,iokRef,argoWork,iokCyc,PorD{1},Work.presEff);

Work.(['DOXY_REF' cmpl{1}]) = REF.(['doxy_CV' cmpl{1}]);
Work.(['PSAT_REF' cmpl{1}]) = REF.(['psat' cmpl{1}]);
Work.(['PSAL_REF' cmpl{1}]) = REF.(['psal' cmpl{1}]);
Work.(['TEMP_REF' cmpl{1}]) = REF.(['temp' cmpl{1}]);
Work.(['PPOX_REF' cmpl{1}]) = REF.(['ppox' cmpl{1}]);

if isfield(argo,'is_main_profile')
    vrf=0;% Cycle number without data
    for k=1:size(REF_ARGO.refId,2)
        if all(isnan(argoWork.doxy_adjusted.data(iokCyc(:,1),:)))
            vrf=k;
        end
    end    
    if vrf>0
        h = warndlg({sprintf('No DOXY data in one of the reference cycles (cycle number %s) \n \n --> Change the reference cycle in the ''bdd_REF_ARGO'' file, or choose another correction !!!',num2str(REF_ARGO.cycle(k)))});
        uiwait(h);     
        Work.goprog=0;
        return
    end
end


% =========================================================================
%% Corrections applied on QC <=2
% =========================================================================

for k=1:size(iokCyc,2)
    cqc=argoWork.doxy_adjusted_qc.data(iokCyc(:,k),:);
    cqc(strfind(cqc,' '))='0';
    cqc=cqc';% Should be a Column Vector

    isnok = find(str2num(cqc) > 1); %#ok<ST2NM>
    if ~isempty(isnok)
        REF.(['doxy_CV' cmpl{1}])(k,isnok) = NaN;
        REF.(['psat' cmpl{1}])(k,isnok) = NaN;
        REF.(['ppox' cmpl{1}])(k,isnok) = NaN;
    end
end
% TR 08.02.2022


Work.(['DOXY_REF' cmpl{1}]) = REF.(['doxy_CV' cmpl{1}]);
Work.(['PSAT_REF' cmpl{1}]) = REF.(['psat' cmpl{1}]);
Work.(['PPOX_REF' cmpl{1}]) = REF.(['ppox' cmpl{1}]);

% Take the REF cycle if INAIR corr --MG
if strcmp(Work.whichCorr,'INAIR') || strcmp(Work.whichCorr,'WOA')
    Work.refCycle=REF_ARGO.cycle;
end

% =========================================================================
%% Control plot
% =========================================================================
if Work.makePlot
    DOXY_PLOT_interpolation(Work, argoWork, 2, iokCyc, REF_ARGO.cycle);
end
