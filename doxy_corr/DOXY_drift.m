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
%            09/04/2021         Thierry Reynaud
%                               Cosmetic Changes "Drift" ==> "Time Drift"
%            12/04/2021         Thierry Reynaud
%                               Drift calculation Depth Variables renamed
%            09/07/2021         Thierry Reynaud fitting degree 2 find
%                               extremum day.
%            11.04.2022         Bug corrected for Work.launchdate (T.Reynaud)
%            31.10.2023         Online questions versus Dialog Box choice added (T. Reynaud and C. Kermabon)


function [Work, DRIFT] = DOXY_drift(Work, argoWork, argo, argoTrajWork)
% =========================================================================
%% Initialisation
% =========================================================================
if ~isfield(Work,'refCyc')
    okdim = strcmp(argoWork.pres_adjusted.dim,'N_PROF');
    dim = size(argoWork.pres_adjusted.data);
    refCyc = 1:dim(okdim);
else
    refCyc = Work.refCyc;
end

fprintf('\t Drift correction is done with PRES\n');

% Define if the drift computation for WOA and REF method is done on deep
% levels or at the sea surface

if    strcmp(Work.whichDrift,'WOA')
    if Work.onlineq
        a=input('Do you want to compute the Time Drift on deep levels (No ==> Time Drift will be computed at the sea surface ? Yes/No','s');% Cathy Kermabon
    else
        a=questdlg('Do you want to compute the drift on deep levels (No means that the drift will be computed at the sea surface ?',sprintf('Levels of WOA_based dirft correction - %d',Work.wmo),'Yes','No','No');
    end

    if strcmp(a,'Yes')
        Work.driftondeeplevels = 1;
    else
        Work.driftondeeplevels = 0;
    end
end

%nbr_lin_reg= Minimum number of points to compute a drift. Default value :10.
nbr_lin_reg=10;

% for saving plot
if Work.savePlot == 1
    Work.savePlot = 0;
    mem = 1;
else
    mem = 0;
end

% depth at which to compute the drift and unit
if strcmp(Work.whichDrift,'WOA') %drift computed with WOA data
    % Work with PPOX to put in work field
    unit = 'millibar';  
    
    if Work.driftondeeplevels
      Zq=Work.min_drift_depth_deep:Work.step_drift_depth_deep:Work.max_drift_depth_deep;
    else
        %VT 23/11/2020 + TR 12.04/2021
        %Work.min_drift_depth=0;% Moved in locodox_config
        %Work.max_drift_depth=25;% Moved in locodox_config
        %Work.step_drift_depth=5;% Moved in locodox_config
        Zq=Work.min_drift_depth_surf:Work.step_drift_depth_surf:Work.max_drift_depth_surf;
    end
elseif strcmp(Work.whichDrift,'NCEP') %drift computed with NCEP data
    Zq = 1;
    unit = argoTrajWork.ppox_doxy_adjusted.units;
else
    disp('The drift computation method has not been defined')
end

% initialize
tmpTab = NaN(length(refCyc),length(Zq));

DRIFT.ppox = tmpTab;
DRIFT.ppoxref = tmpTab;
DRIFT.gainppox = tmpTab;
DRIFT.gainppox_Mean = NaN(length(refCyc),1);

DRIFT.applyDriftC='NO';

% =========================================================================
%% Prepare data for drift computation
% =========================================================================
% -------------------------------------------------------------------------
% For drift compute on WOA
% Interpolate data at deep depth, every 100m
% If no data under min(Zq), no interpolation (NaN)
% -------------------------------------------------------------------------
%if strcmp(whichCorr,'WOA') || strcmp(whichCorr,'REF')
if strcmp(Work.whichDrift,'WOA')
    
    % Nb of day after deployment
    %daydiff = (argo.juld.data - argo.juld.data(1));
    % Nb of day since launch date
    daydiff = (argo.juld.data - Work.launchdate);
    
    
    
    if length(refCyc)<nbr_lin_reg
        fprintf('\t No Drift Computation: not enougth data for this float');
        waitfor(warndlg(sprintf('\t No Drift Computation: not enougth data for this float (%d cycle) \n', length(refCyc)) ,...
            'DOXY_drift.m')); %marine 11/06/19
        Work.DODRIFT = false;
        
        if Work.makePlot
            
            DOXY_PLOT_drift(DRIFT, Work, daydiff,unit,1);
            if mem == 1  %marine 11/06/19
                Work.savePlot = 1;
            end
            return
        end
        
    else
        for i = refCyc
            if Work.driftondeeplevels
                isok = find(~isnan(Work.PRES(i,:)) & Work.PRES(i,:) >= Work.min_drift_depth_deep & Work.PRES(i,:) <= Work.max_drift_depth_deep);
            else
                isok = find(~isnan(Work.PRES(i,:)) & Work.PRES(i,:) >= Work.min_drift_depth_surf & Work.PRES(i,:) <= Work.max_drift_depth_surf);
            end  
            Zar = Work.PRES(i,isok);
            
            % At least two points to allow interpolate
            if  length(Zar) < 2
                DRIFT.ppox(i,:) = NaN;
                DRIFT.ppoxref(i,:) = NaN;
            else
                Var = Work.PPOX(i,isok);
                Vref = Work.PPOX_WOA_interpP(i,isok);
                
                DRIFT.ppox(i,:) = interp1(Zar,Var,Zq,'linear');
                DRIFT.ppoxref(i,:) = interp1(Zar,Vref,Zq,'linear');
            end
            DRIFT.gainppox(i,:) = (DRIFT.ppoxref(i,:)./DRIFT.ppox(i,:))';
            DRIFT.gainppox_Mean(i,:) = nanmean(DRIFT.gainppox(i,:));
        end
        
        if all(all(isnan(DRIFT.ppox)))
            fprintf('\t No Drift Computation: no data below %d m\n', min(Zq));
            waitfor(warndlg(sprintf('\t No Drift Computation: no data below %d m\n --> The minimum depth for drift calculation can be changed in the configuration file\n', min(Zq)),...
                'DOXY_drift.m')); %marine 11/06/19
            Work.DODRIFT = false;
            if Work.makePlot
                DOXY_PLOT_drift(DRIFT, Work, daydiff,unit,1);
                if mem == 1  %marine 11/06/19
                    Work.savePlot = 1;
                end
                return
            end
        end
    end
end

% -------------------------------------------------------------------------
% For INAIR correction
% -------------------------------------------------------------------------

%if strcmp(whichCorr,'INAIR')
if strcmp(Work.whichDrift,'NCEP')
    DRIFT.ppox = Work.surface.argo.pO2_air;
    DRIFT.ppoxref = Work.surface.ncep.pO2;
    DRIFT.gainppox = (DRIFT.ppoxref./DRIFT.ppox)';
    DRIFT.gainppox_Mean = nanmean(DRIFT.gainppox);
    
    % Nb of day after deployment
    %   daydiff = (argoTrajWork.juld.data - min(argoTrajWork.juld.data))';
    %--EB : daymin was set to argoTrajWork.juld.data{1}(1), whereas it
    % could be null or empty
    if Work.NSIAfloat
        %daymin = Work.surface.argo.juld_air(1);
        % Added by T. Reynaud 11.04.2022
        if exist('Work.launchdate','var')
            daymin=Work.launchdate;
        else
            tmp=load('tmp.mat');
            daymin=tmp.Work.launchdate;
            delete('tmp.mat');
        end

        daydiff = Work.surface.argo.juld_air - daymin;
        nbr_meas=Work.nbrMeas.inAir;
    else
        %daymin = argoTrajWork.juld.data(1);
        daymin=Work.launchdate;
        daydiff = argoTrajWork.juld.data - daymin;
    end
end

%--EB - 20180418
% Uniformize the data format (no single and double in the same field)
if iscell(DRIFT.ppox)
    driftFields = fieldnames(DRIFT);
    nbfield = numel(driftFields);
    for f = 1:nbfield
        tmpfield = driftFields{f};
        DRIFT.(tmpfield) = cellfun(@double,DRIFT.(tmpfield),'UniformOutput',false);
    end
end
%--EB - 20180418

% -------------------------------------------------------------------------
%% Control Plot
% Plot interpolated doxy from argo and doxy from WOA respect to the number
% of day after deployement, and plot the difference.
% Or plot the ppox from argo inflated and from NCEP respect to the number
% of day after deployement, and plot the difference.
% -------------------------------------------------------------------------
if Work.makePlot
    [ax,~] = DOXY_PLOT_drift(DRIFT, Work, daydiff,unit,1);
end

if iscell(daydiff)
    daydiff = cell2mat(daydiff);
end

% =========================================================================
%% Drift computation
%  Compute the linear regression off the doxy difference timeserie
%  (argo - woa, or argo - NCEP).
%  Based on Takeshita
%VT
%  Compute the linear regression off the ppox ratio
% VT
% =========================================================================

if Work.drift_spec == 1
    v =  ver;
    if all(cellfun('isempty',strfind({v.Name},'Statistics')))
        inputval = inputdlg({'Statistics Toolbox unavailable : Matlab can''t apply your drift equation. The drift will be compute with a polynomial equation. The default degree is 1. Do you wan''t to change it ?'},...
            'DOXY_drift : non statistic toolbox',1,{'1'});
        Work.drift_spec = 0;
        Work.drift_fitPolynomialDegree = str2num(char(inputval)); %#ok<ST2NM>
    end
end
if Work.drift_spec == 0
    coeffFit = NaN(Work.drift_fitPolynomialDegree+1,length(Zq));
else
    coeffFit = NaN(numcoeffs(Work.drift_fittype),length(Zq));
end
cpt = 0;
noDriftComputation = false;

for z=1:length(Zq)
    % linear regression
    %if strcmp(whichCorr,'WOA') || strcmp(whichCorr,'REF')
    if strcmp(Work.whichDrift,'WOA')
        
        isok = ~isnan(daydiff) & ~isnan(DRIFT.gainppox(:,z));
        if sum(isok) >= nbr_lin_reg && any(~isnan(DRIFT.gainppox(:,z)))
            if Work.drift_spec == 0
                fitResult = polyfit(double(daydiff(isok)),double(DRIFT.gainppox(isok,z)),Work.drift_fitPolynomialDegree);
                coeffFit(:,z) = fitResult;
            else
                fitResult = fit(double(daydiff(isok)),double(DRIFT.gainppox(isok,z)),Work.drift_fittype);
                coeffFit(:,z) = coeffvalues(fitResult);
            end
            
            % control plot : check the coefficient and the linear regression at
            % a local depth
            if Work.makePlot
                if Work.drift_spec == 0
                    DRIFT.poliv_c = polyval(coeffFit(:,z), double(daydiff));
                else
                    DRIFT.poliv_c = fitResult(double(daydiff));
                end
                [ax,~] = DOXY_PLOT_drift(DRIFT,Work,daydiff,unit,2,ax);
                DRIFT = rmfield(DRIFT,'poliv_c');
            end
            
        else
            cpt = cpt+1;
            if z==length(Zq) && cpt==length(Zq)
                errorMess = strcat('Less than [', num2str(nbr_lin_reg) ,'] valid data below 1500m for each profile');
                noDriftComputation = true;
            end
        end
        
    %elseif strcmp(whichCorr,'INAIR')
    elseif strcmp(Work.whichDrift,'NCEP')
        data = double(DRIFT.gainppox)';
        isok = ~isnan(data) & ~isnan(daydiff);
        ind_conv=0;
        
%         if isfield(Work,'fitResult_WOA')
%             coeffFit(:,z) = Work.fitResult_WOA;
%             ind_conv=1; %PPOX has to be converted to apply the drift
%         else
        if sum(isok) >= nbr_lin_reg %|| any(~cellfun(@isnan,DRIFT.gainppox,'UniformOutput',false))
            if ~Work.drift_spec
                fitResult = polyfit(daydiff(isok),data(isok),Work.drift_fitPolynomialDegree);
                coeffFit(:,z) = fitResult;
            else
                % fit doesn't accept row vectors.
                if ~iscolumn(daydiff)
                    daydiff = daydiff';
                end
                if ~iscolumn(data)
                    data = data';
                end
                fitResult = fit(daydiff(isok),data(isok),Work.drift_fittype);
                coeffFit(:,z) = coeffvalues(fitResult);
            end
            % control plot : check the coefficient and the linear regression at
            % a local depth
            if Work.makePlot
                if Work.drift_spec == 0
                    DRIFT.poliv_c = polyval(coeffFit(:,z), double(daydiff));
                else
                    DRIFT.poliv_c = fitResult(double(daydiff));
                end
                [ax,~] = DOXY_PLOT_drift(DRIFT,Work,daydiff,unit,2,ax);
                DRIFT = rmfield(DRIFT,'poliv_c');
            end
        else
            cpt = cpt+1;
            if z==length(Zq) && cpt==length(Zq)
                errorMess = 'Less than 10 valid data at surface';
                noDriftComputation = true;
            end
        end
    end
end

% =========================================================================
%Compute drift only if it is possible
% =========================================================================
if ~noDriftComputation
    % -------------------------------------------------------------------------
    % Slope averaging
    % -------------------------------------------------------------------------
    coeffFit_mean = nanmean(coeffFit,2);
    DRIFT.coeffFit_mean=coeffFit_mean;
    
    if Work.drift_fitPolynomialDegree == 2
        first_derivative = polyder(DRIFT.coeffFit_mean);%compute first derivative
        rts = roots(first_derivative);%find roots
        display(strcat('DAY EXTREMUM= ',num2str(round(rts))));
    end  
    
    isok = ~isnan(daydiff);
    % control plot : check the mean coefficient and the linear regression for
    % all depth
    if Work.makePlot
        if mem == 1
            Work.savePlot = 1;
        end
        if Work.drift_spec == 0
            DRIFT.fitRegression = polyval(coeffFit_mean,double(daydiff));
        else
            listCoef = [];
            for i=1:length(coeffFit_mean)
                listCoef = [listCoef,sprintf('coeffFit_mean(%d),',i)];
            end

            eval(['c = cfit(Work.drift_fittype,' listCoef ');']);
            DRIFT.fitRegression = feval(c,double(daydiff));
        end
        [~,Work]=DOXY_PLOT_drift(DRIFT,Work,daydiff(isok),unit,3, ax);
    end
    
    %Recompute fitRegression with offset = 0 now
    %if strcmp(whichCorr,'INAIR') && ind_conv~=1 %Case INAIR correction and drift computed on NCEP
    if strcmp(Work.whichDrift,'NCEP')
        %Change of time we use profile time and no more in air or in water
        %daydiff = (argo.juld.data - argo.juld.data(1));
        daydiff = (argo.juld.data - argo.launchdate);
        
        if Work.drift_spec == 0
            DRIFT.fitRegression = polyval(coeffFit_mean,double(daydiff));
            DRIFT.daydiff=daydiff;% Added by T. Reynaud 16/07/2021
        else
            listCoef = [];
            for i=1:length(coeffFit_mean)
                listCoef = [listCoef,sprintf('coeffFit_mean(%d),',i)];
            end
            eval(['c = cfit(Work.drift_fittype,' listCoef ');']);
            DRIFT.fitRegression = feval(c,double(daydiff));
            DRIFT.daydiff=daydiff;% Added by T. Reynaud 16/07/2021
        end
        if size(DRIFT.fitRegression,2) == size(Work.DOXY,1)
            DRIFT.fitRegression = DRIFT.fitRegression';
            DRIFT.daydiff=daydiff;% Added by T. Reynaud 16/07/2021
        end
    else
        if Work.drift_spec == 0
            DRIFT.fitRegression = polyval(coeffFit_mean,double(daydiff));
            DRIFT.daydiff=daydiff;% Added by T. Reynaud 16/07/2021
        else
            listCoef = [];
            for i=1:length(coeffFit_mean)
                listCoef = [listCoef,sprintf('coeffFit_mean(%d),',i)];
            end
            %VT on uitilise l'offset
                %listCoef(end) = '';
            eval(['c = cfit(Work.drift_fittype,' listCoef ');']);
            DRIFT.fitRegression = feval(c,double(daydiff));
            DRIFT.daydiff=daydiff;% Added by T. Reynaud 16/07/2021
        end
    end
    
    % Compute drift for in air and in water data
    % Change time reference for in air and in water data
   % if strcmp(whichCorr,'INAIR') && ~isfield(Work,'isnotINAIRcorr')
    if strcmp(Work.whichCorr,'INAIR') || strcmp(Work.whichDrift,'NCEP')
        
        if Work.NSIAfloat
            %daymin_airwater = min(Work.surface.argo.juld_air);
            daymin_airwater = argo.launchdate;
            daydiff_airwater = Work.surface.argo.juld_air - daymin_airwater;
        else
            %daymin_airwater = min(min(argoTrajWork.juld.data));
            daymin_airwater = argo.launchdate;
            daydiff_airwater = argoTrajWork.juld.data - daymin_airwater;
        end
        %Recompute drift
        if Work.drift_spec == 0
            DRIFT.fitRegression_airwater = polyval(coeffFit_mean,double(daydiff_airwater));
            DRIFT.daydiff_airwater=daydiff_airwater;%Added by T. Reynaud 16/07/2021
        else
            listCoef_airwater = [];
            for i=1:length(coeffFit_mean)
                listCoef_airwater = [listCoef_airwater,sprintf('coeffFit_mean(%d),',i)];
            end
            %VT on uitilise l'offset
            %listCoef_airwater(end) = '';
            eval(['c = cfit(Work.drift_fittype,' listCoef_airwater ');']);
            DRIFT.fitRegression_airwater = feval(c,double(daydiff_airwater));
            DRIFT.daydiff_airwater=daydiff_airwater;%Added by T. Reynaud 16/07/2021
        end
        if size(DRIFT.fitRegression,2) == size(Work.DOXY,1)
            DRIFT.fitRegression_airwater = DRIFT.fitRegression_airwater';
            DRIFT.daydiff_airwater=daydiff_airwater;%Added by T. Reynaud 16/07/2021
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
            DRIFT.daydiff_airwater=daydiff_airwater;%Added by T. Reynaud 16/07/2021
        end
        
        
    end
    
    % -------------------------------------------------------------------------
    % Annual drift
    % -------------------------------------------------------------------------
    %if ~strcmp(whichCorr,'INAIR')
    if strcmp(Work.whichDrift,'WOA')    
        Work.DRIFT = (DRIFT.fitRegression(end) - DRIFT.fitRegression(1)) / ...
            ((nanmax(argo.juld.data)-nanmin(argo.juld.data))/365);
    else
        if iscell(argoTrajWork.juld.data)
            Work.DRIFT = (DRIFT.fitRegression(end) - DRIFT.fitRegression(1)) / ...
                ((nanmax(cell2mat(argoTrajWork.juld.data))-nanmin(cell2mat(argoTrajWork.juld.data)))/365);
        else
            Work.DRIFT = (DRIFT.fitRegression(end) - DRIFT.fitRegression(1)) / ...
                ((nanmax(argoTrajWork.juld.data)-nanmin(argoTrajWork.juld.data))/365);
        end
    end
    
    
    % =========================================================================
    %% Asking user for drift
    % =========================================================================
    if isfield(Work,'ind_drift_stop')
        DRIFT.explanation3=sprintf('A constant drift is applied from day number %d',Work.ind_drift_stop(2));
    else
        DRIFT.explanation3=' ';
    end
    
    aux_year=find((abs(daydiff-365.25))==min(abs(daydiff-365.25))); %% position after a year life
    nyears_floatlife = (argo.juld.data(end)-argo.juld.data(1))/365;
    
    %Display leading coefficient %marine
    if strcmp(Work.whichCorr,'INAIR')
        if Work.drift_fitPolynomialDegree == 1
            %VT
            %DRIFT.explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) + %2.4f * t \n',coeffFit_mean(1));
            %DRIFT.explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) / ( %2.4f * t + %2.4f) \n',coeffFit_mean(1),coeffFit_mean(2));
            % TR 14.04.2021
            DRIFT.explanation0 = sprintf('Time Drift Gain(1  ): %5.4f %s \n',DRIFT.fitRegression(1));
        elseif Work.drift_fitPolynomialDegree == 2
            %VT
            %DRIFT.explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) / ( %2.4f * t^2 + %2.4f * t + %2.4f) \n',coeffFit_mean(1),coeffFit_mean(2),coeffFit_mean(3));
            % TR 14.04.2021
            DRIFT.explanation0 = sprintf('Time Drift Gain(1  ): %5.4f %s \n',DRIFT.fitRegression(1));
        else
            DRIFT.explanation0 = '';
        end
    else
        if Work.drift_fitPolynomialDegree == 1
            %VT
            %DRIFT.explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) + %2.4f * t \n',coeffFit_mean(1));
            %DRIFT.explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) / ( %2.4f * t + %2.4f) \n',coeffFit_mean(1),coeffFit_mean(2));
            % TR 14.04.2021
            DRIFT.explanation0 = sprintf('Time Drift Gain(1  ): %5.4f %s \n',DRIFT.fitRegression(1));

        elseif Work.drift_fitPolynomialDegree == 2
            %VT
            DRIFT.explanation0 = sprintf('The best fitting was: PO2_corr(t) ~= PO2(t) / ( %2.4f * t^2 + %2.4f * t + %2.4f) \n',coeffFit_mean(1),coeffFit_mean(2),coeffFit_mean(3));
            % TR 14.04.2021
            DRIFT.explanation0 = sprintf('The Time Drift Gain(1  ): %5.4f %s \n',DRIFT.fitRegression(1));

        else
            DRIFT.explanation0 = '';
        end
%         %VT 23/11/2020
%         % need to add case compute drift on PSAT
%         if Work.drift_fitPolynomialDegree == 1
%             DRIFT.explanation0 = sprintf('The best fitting was : DOXY_corr(t) = DOXY(t) + %2.4f * t \n',coeffFit_mean(1));
%         elseif Work.drift_fitPolynomialDegree == 2
%             DRIFT.explanation0 = sprintf('The best fitting was : DOXY_corr(t) = DOXY(t) + %2.4f * t + %2.4f * t^2\n',coeffFit_mean(1),coeffFit_mean(2));
%         elseif Work.drift_fitPolynomialDegree == 3
%             DRIFT.explanation0 = sprintf('The best fitting was : DOXY_corr(t) = DOXY(t) + %2.4f * t + %2.4f * t^2 + %2.4f * t^3\n',coeffFit_mean(1),coeffFit_mean(2),coeffFit_mean(3)); '';
%         end
     end
    
    %Display drift after a year and at the end of life
    if nyears_floatlife>=1
        %DRIFT.explanation1 = sprintf('The drift after a year is %4.3f %s \n',DRIFT.fitRegression(aux_year), unit);
        %DRIFT.explanation2 = sprintf('The drift at the end of the float life ( %4.2f years) is %4.3f %s \n', nyears_floatlife, DRIFT.fitRegression(end),unit);
        % TR 14.04.2021
        DRIFT.explanation1 = sprintf('Time Drift Gain(1yr): %5.4f %s \n',DRIFT.fitRegression(aux_year));
        DRIFT.explanation2 = sprintf('Time Drift Gain(end): %5.4f %s \n',DRIFT.fitRegression(end));
    else
        %DRIFT.explanation1 = sprintf('The drift at the end of the float life ( %4.2f years) is %4.3f %s \n', nyears_floatlife, DRIFT.fitRegression(end),unit);
        %DRIFT.explanation2 = sprintf(' \n');
        % TR 14.04.2021
        DRIFT.explanation1 = sprintf('Time Drift Gain(end): %5.4f %s \n',DRIFT.fitRegression(end));
        DRIFT.explanation2 = sprintf(' \n');
    end
    
    %Message of sugestion depending of the availability of data
    deltaO2MeanNok = sum(isnan(DRIFT.gainppox_Mean));
    if daydiff(end) < 365
        % Less than one year of data
        driftCorrSuggestion = sprintf('* Suggestion: NO DRIFT CORRECTION:');
        DRIFT.explanation = sprintf('* State: less than 1 year of data');
    elseif length(argo.cycle_number.data) - deltaO2MeanNok<nbr_lin_reg+1
        % less than 5 profils with depth > 1500m
        driftCorrSuggestion = sprintf('* Suggestion: NO DRIFT CORRECTION:');
        DRIFT.explanation =sprintf( '* State: less than 10 profiles > 1500m  ATTENTION COMMENTAIRE A REVOIR');
    else
        DRIFT.explanation = ' ';
        driftCorrSuggestion= ' ';
    end
    
    %Choice of user
%     if strcmp(whichCorr,'INAIR') && strcmp(Work.whichDrift,'WOA')
%         if isfield(Work,'fitResult_WOA')
%             DRIFT.applyDriftC='YES';
%         else
%             DRIFT.applyDriftC='NO';
%         end
%     else
%        DRIFT.applyDriftC = questdlg({sprintf('Float #%d',Work.wmo);...
%            ' ---------------------------';
%            'The statistical method results in: '; ...
%            ' ---------------------------';DRIFT.explanation0;
%            DRIFT.explanation1; DRIFT.explanation2; ...
%            '---------------------------';...
%            driftCorrSuggestion; DRIFT.explanation; DRIFT.explanation3;...
%            ' ---------------------------';...
%            'Apply time drift correction ?';''},...
%            'TIME DRIFT CORRECTION','YES','NO','YES');
    if Work.onlineq
       DRIFT.applyDriftC = input('Apply the Time Drift correction ? (YES/NO)','s');
    else
       DRIFT.applyDriftC = questdlg('Apply the Time Drift correction ?','Time Drift Gain','YES','NO','YES');
    end
%     end
    
else
    if mem==1
        Work.savePlot=1;
    end
    driftCorrSuggestion = '* DRIFT CORRECTION NOT POSSIBLE : no data';
    DRIFT.explanation = sprintf('* State: %s',errorMess);
    DRIFT.applyDriftC='NO';
    % to small drift
    waitfor(warndlg({driftCorrSuggestion;DRIFT.explanation},'NO DRIFT COMPUTATION')); %marine 12/06/19
end

