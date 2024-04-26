% DOXY_drift corrects data from sensor drift
%
% SYNTAX
% [Work,DRIFT] = DOXY_drift(Work, argoWork, argo, refCyc, whichCorr, argoTrajWork)
%
% DESCRIPTION
% DOXY_drift compare the doxy argo data to the WOA data or NCEP data and estimates a
% sensor drift. If needed, a correction is applied. The drift estimation is
% made at depth below the user choice for WOA data. The WOA data and the argo data are both
% interpolate on a pressure vertical grid every 100m from the depth chosen by the user
% (see the configuration file) to 6000m.
% For drift computes on NCEP data, in-air measurement are used.
%
% INPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA)
%                           Example:
%                           Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                           where
%                           Work.DOXY_WOA =
%                           name: 'DOXY_ADJUSTED'
%                            dim: {'N_PROF'  'N_LEVELS'}
%                           data: [85x120 double]
%                      long_name: 'DISSOLVED OXYGEN adjusted and converted'
%                  standard_name: 'moles_of_oxygen_per_unit_mass_in_sea_water'
%                     FillValue_: 99999
%                          units: 'mumol/kg'
%                      valid_min: 0
%                      valid_max: 650
%                       C_format: '%9.3f'
%                 FORTRAN_format: 'F9.3'
%                     resolution: 0.0010
%                           type: 5
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
%   argo (structure)      Argo float structure (get directly from data).
%                         Example :
%                         argo =
%                                      pres: [1x1 struct]
%                             pres_adjusted: [1x1 struct]
%                                      doxy: [1x1 struct]
%                                   doxy_qc: [1x1 struct]
%
%                               argo.pres
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
%  whichCorr (string)      Indicates which kind of correction is used
%                          Example: whichCorr='REF', whichCorr='WOA',
%                          whichCorr = 'INAIR';
%
%  argoTrajWork (struct)   Float trajectory data structure
%
% OUTPUT
%   Work (structure)     float working structure, issued and computed from
%                        the argo data and the climatology data (WOA),
%                        completed with the annual drift.
%                        Example:
%                        Work =
%                             pres_adjusted: [1x1 struct]
%                             temp_adjusted: [1x1 struct]
%                             psal_adjusted: [1x1 struct]
%                                  DOXY_WOA: [1x1 struct]
%                                     DRIFT: 4.8400
%
%   DRIFT (structure)    Drift estimating structure, from argo data and
%                        climatology data.
%                        Example:
%                        DRIFT =
%                                      doxy: [46x85 double]
%                                   doxyref: [46x85 double]
%                                   deltaO2: [46x85 double]
%                              deltaO2_Mean: [1x85 double]
%                            doxy_corrected: [85x120 double]
%                          fitRegression: [1x85 double]
%
%
%
% CALL :
%   DOXY_PLOT_drift, O2ptoO2c, O2ptoO2s, O2stoO2c
%
% SEE ALSO
%   DOXY_corr_main
%
% HISTORY
%   $created: 07/12/2015 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%             25/03/2019        Marine GALLIAN, Altran Ouest
%                               Fill the summary file
%                               Manage onvertion for drift computed on NCEP
%                               and WOA
%                               Improve message box
%             30/03/2019        Marine GALLIAN, Altran Ouest
%                               Now drift is computed using PRES and no
%                               more using the DEPTH variable
%       v3.4 10.04.2020         T.Reynaud temporal ==> time
%            23.11.2020         V. Thierry added Time Drift Correction on PSAT
%            23.11.2020         T. Reynaud driftonDOXY and driftonPSAT read
%                               from Work
%      v3.5  09.07.2021         For NCEP time drift the constant time drift
%                               index recalculated
%      v3.6  26.04.2024         PWLF option added for Time Drift Gain Correction (T. Reynaud)


function [Work, DRIFT] = DOXY_drift_apply(Work, argoWork, argo, argoTrajWork, DRIFT)

% =========================================================================
%% apply drift to data
% =========================================================================
switch DRIFT.applyDriftC
    case 'YES'
        fprintf('\t The drift is applied \n')
        % =========================================================================
        %% Compute drift correction over profile and trajectory
        % =========================================================================
        Work.DODRIFT = true;
        % -------------------------------------------------------------------------
        % Compute the profile correction
        % -------------------------------------------------------------------------
        
        %Variable Initialisation
        okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');
        dim = size(argoWork.pres_adjusted.data);
        DRIFT.ppox_corrected = NaN(dim);
        isub=0;% Added by T. Reynaud 16/07/2021

        
        %     if size(DRIFT.fitRegression,2) == size(Work.DOXY,1)
        %         DRIFT.fitRegression = DRIFT.fitRegression';
        %     end
        
        
        %Compute constant drift
        if isfield(Work,'ind_drift_stop')
            % Modified by T. Reynaud 09.07.2021
            ibeg=Work.ind_drift_stop(1);
            if strcmpi(Work.whichDrift,'NCEP')
                tmp_deltat=argo.juld.data-argo.launchdate;
                ibeg=find(tmp_deltat>=Work.ind_drift_stop(2),1,'first');
            end
            DRIFT.fitRegression(ibeg:end)=DRIFT.fitRegression(ibeg);
            tmp_deltat2=DRIFT.daydiff; % time base interpolated coefficient
            %DRIFT.fitRegression(Work.ind_drift_stop(1):end)=DRIFT.fitRegression(Work.ind_drift_stop(1));
            if size(DRIFT.fitRegression,1)<size(Work.PPOX,1) % Added by T.Reynaud 16.07.2021
                tmp_deltat2=DRIFT.daydiff; % time base interpolated coefficient
                DRIFT.fitRegression2=interp1(tmp_deltat2,DRIFT.fitRegression,tmp_deltat);
                isub=1;
            end
        end
        
        %Apply drift ==> Modified by T. Reynaud 16/07/2021
        for i = 1:dim(okdim)
            if isub
                DRIFT.ppox_corrected(i,:) = Work.PPOX(i,:) .* DRIFT.fitRegression2(i);
            else
                if i<=length(DRIFT.fitRegression)
                    DRIFT.ppox_corrected(i,:) = Work.PPOX(i,:) .* DRIFT.fitRegression(i);
                else
                    DRIFT.ppox_corrected(i,:) = Work.PPOX(i,:) .* DRIFT.fitRegression(end);% Debug TR 26/08/2022
                end
            end
        end
        
        %reconvert PPOX to DOXY and PSAT
        if Work.presEff == 1, P = argoWork.pres_adjusted.data; else, P = 0; end
        tmpDoxy = O2ptoO2c(DRIFT.ppox_corrected,argoWork.temp_adjusted.data, argoWork.psal_adjusted.data,P);
        DRIFT.doxy_corrected= DOXY_convert(tmpDoxy,'mumol/L',argoWork.doxy_adjusted.units,argoWork.an_dens.data);
        DRIFT.psat_corrected=O2ptoO2s(DRIFT.ppox_corrected,argoWork.temp_adjusted.data, argoWork.psal_adjusted.data,P);
        drift_on=Work.whichDrift;
        
        Work.PPOX_DRIFTCORR_COEF = DRIFT.coeffFit_mean;
        Work.PPOX_DRIFTCORR = DRIFT.ppox_corrected;
        Work.PPOX_fitRegression=DRIFT.fitRegression;
        Work.PPOX_daydiff=DRIFT.daydiff;

        if Work.drift_PWLF % Added T.Reynaud 16.04.2024
            Work.PPOX_DRIFTCORR_SEG=DRIFT.daydiff_Seg;
            num_seg=Work.drift_PWLF_N;
        else
            num_seg=1;
        end
        
        if strcmp(Work.whichCorr,'INAIR')
            
            if strcmp(Work.whichDrift,'WOA')
                
                %VT
                
                if Work.NSIAfloat
                    %daymin_airwater = min(Work.surface.argo.juld_air);
                    daydiff_airwater = Work.surface.argo.juld_air - argo.launchdate;
                else
                    %daymin_airwater = min(min(argoTrajWork.juld.data));
                    daydiff_airwater = argoTrajWork.juld.data - argo.launchdate;
                end
                %Recompute drift
                if Work.drift_spec == 0
                    if ~Work.drift_PWLF % Added T.Reynaud 16.04.2024
                        DRIFT.fitRegression_airwater = polyval(DRIFT.coeffFit_mean,double(daydiff_airwater));
                    else
                        DRIFT.fitRegression_airwater = polyval_PWLF(DRIFT.coeffFit_mean,DRIFT.daydiff_Seg,double(daydiff_airwater));
                    end
                else
                    listCoef_airwater = [];
                    for i=1:length(DRIFT.coeffFit_mean)
                        listCoef_airwater = [listCoef_airwater,sprintf('DRIFT.coeffFit_mean(%d),',i)];
                    end
                    %VT on uitilise l'offset
                    %listCoef_airwater(end) = '';
                    eval(['c = cfit(Work.drift_fittype,' listCoef_airwater ');']);
                    DRIFT.fitRegression_airwater = feval(c,double(daydiff_airwater));
                end
                if size(DRIFT.fitRegression,2) == size(Work.DOXY,1)
                    DRIFT.fitRegression_airwater = DRIFT.fitRegression_airwater';
                end
                
                %Apply constant drift
                if isfield(Work,'ind_drift_stop')
                    for i=1:length(daydiff_airwater)
                        if daydiff_airwater(i)>Work.ind_drift_stop(2) && Work.ind_drift_stop(2)>-1
                            Work.ind_drift_stop_airwater=[i,Work.ind_drift_stop(2)];
                            break
                        end
                    end
                    DRIFT.fitRegression_airwater(Work.ind_drift_stop_airwater(1):end)=DRIFT.fitRegression_airwater(Work.ind_drift_stop_airwater(1));
                end
            end
            
            DRIFT.pO2_corrected.inair = Work.surface.argo.pO2_air .* DRIFT.fitRegression_airwater;
            DRIFT.pO2_corrected.inwater = Work.surface.argo.pO2_water .* DRIFT.fitRegression_airwater;
            Work.PPOX_DRIFTCORR_inair = DRIFT.pO2_corrected.inair;
            Work.PPOX_DRIFTCORR_inwater = DRIFT.pO2_corrected.inwater;
        end
        
        
        if ~isfield(Work,'isnotREFcorr')
            %Modifie output of the function
            DRIFT.psat_corrected(DRIFT.psat_corrected<=0) = NaN;
            DRIFT.doxy_corrected(DRIFT.doxy_corrected<=0) = NaN;
            DRIFT.ppox_corrected(DRIFT.ppox_corrected<=0) = NaN;
            Work.DOXY_DRIFTCORR = DRIFT.doxy_corrected;
            Work.PSAT_DRIFTCORR = DRIFT.psat_corrected;
            Work.PPOX_DRIFTCORR = DRIFT.ppox_corrected;
            
            % -------------- SUMMARY FILE -----------------
            %if ~(isfield(Work,'fitResult_WOA') && strcmp(whichCorr,'INAIR')) %condition to avoid writing in double the drift info
            fprintf(Work.log_summ,'Time drift correction applied : \n');
            if Work.drift_spec == 1
                fprintf(Work.log_summ,sprintf('Time drift computed on %s and estimated by applying a fitting proposed by the user. The fitting equation is %s\n',drift_on,formula(Work.drift_fittype)));
                fprintf(Work.log_summ,DRIFT.explanation1);
                fprintf(Work.log_summ,DRIFT.explanation2);
            elseif Work.drift_fitPolynomialDegree == 1
                fprintf(Work.log_summ,sprintf('Time drift computed on %s and estimated by applying a linear regression.\n',drift_on));
                fprintf(Work.log_summ,DRIFT.explanation0);
                fprintf(Work.log_summ,DRIFT.explanation1);
                fprintf(Work.log_summ,DRIFT.explanation2);
            elseif Work.drift_fitPolynomialDegree == 2
                fprintf(Work.log_summ, sprintf('Time drift computed on %s and estimated by applying a polynomial regression of order two.\n',drift_on));
                fprintf(Work.log_summ,DRIFT.explanation0);
                fprintf(Work.log_summ,DRIFT.explanation1);
                fprintf(Work.log_summ,DRIFT.explanation2);
            elseif Work.drift_fitPolynomialDegree == 3
                fprintf(Work.log_summ, sprintf('Time drift computed on %s and estimated by applying a polynomial regression of order three.\n',drift_on));
                fprintf(Work.log_summ,DRIFT.explanation0);
                fprintf(Work.log_summ,DRIFT.explanation1);
                fprintf(Work.log_summ,DRIFT.explanation2);
            end
            
            if isfield(Work,'ind_drift_stop')
                fprintf(Work.log_summ,sprintf('A constant drift was applied from day number %d\n',Work.ind_drift_stop(2)));
            end
            fprintf(Work.log_summ,'\n');
            fprintf(Work.log_summ,'-----------------------------------------------------------\n');
        end
        %end
        
        
        
    case 'NO'
        fprintf('\t Time drift is not applied \n')
        Work.DODRIFT = false;
        DRIFT.psat_corrected = Work.PSAT;
        DRIFT.doxy_corrected = Work.DOXY;
        DRIFT.ppox_corrected = Work.PPOX;
        
        % -------------- SUMMARY FILE -----------------
        fprintf(Work.log_summ,'Time drift correction not applied. \n');
        fprintf(Work.log_summ,'\n');
        fprintf(Work.log_summ,'-----------------------------------------------------------\n');
        
end
