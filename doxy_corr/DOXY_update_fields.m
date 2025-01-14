% DOXY_update_fields updates main fields linked to the history of data
% control
%
% SYNTAX
% [argo, DIMD] = DOXY_update_fields(Work,argo,DIM,REF_ARGO,kwrite)
%
% DESCRIPTION
% DOXY_update_fields updates the following fields :
%   date_update
%   data_mode
%   data_state_indicator
%   parameter_data_mode
%   history fields
%   scientific calib fields
%   doxy_adjusted_error
%
% INPUT
%     argo (structure)    Argo float structure (get directly from data).
%                         Example :
%                               argo =
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
%     Work (structure)  Doxy correction working structure, issued and
%                       computed from argo float data.
%                       Example:
%                             Work =
%                                   readme: [1x1 struct]
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
%
%     DIM (structure)    Argo float Dimension structure (get directly
%                        from data)
%                        Example:
%                         DIM =
%                             date_time: [1x1 struct]
%                             string256: [1x1 struct]
%                              string64: [1x1 struct]
%                              string32: [1x1 struct]
%                              string16: [1x1 struct]
%                               string8: [1x1 struct]
%                               string4: [1x1 struct]
%                               string2: [1x1 struct]
%                                n_prof: [1x1 struct]
%                               n_param: [1x1 struct]
%
%   REF_ARGO (structure)    List the WMO of argo float with DOXY sensor,
%                           and useful information : sensor, DAC,
%                           directory, reference profile. Example:
%                            REF_ARGO =
%                                 cycle : [1x1 double]
%                                 refId : [1x1 string]
%                                 capteur: [1x1 string]
%                                 daclist: [1x1 string]
%                                 wmo : [1x1 double]
%                                 replist: [1x1 string]
%
%   kwrite (double)        writing selection number added for PWLF
%
% OUTPUT
%     argo (structure)      Argo float data structure (get directly from data),
%                           updated for metadata linked to history
%                           date_update, data_mode, data_state_indicator,
%                           parameter_data_mode, history fields, scientific
%                           calib fields, doxy_adjusted_error
%
%     DIM (structure)    Argo float Dimension structure (get directly
%                        from data), updated
%
% CALL :
%   check_FirstDimArray, extract_profile_dim, cat_profile_dim,
%
% SEE ALSO
%   DOXY_corr_main

% HISTORY
%   $created: 12/01/2016 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%       v2.1, 25/04/2016, Emilie Brion, Altran Ouest :
%                         start and end pressure taking among not FillValues
%       v3  26/01/2017    Anne Piron, Altran Ouest
%                         argo - ameliorations 2017 phase 1
%           11/10/2018    Marine GALLIAN, Altran Ouest :
%                         Addition of informations for the drift and
%                         the anterrior correction in scientific calibration
%                         fields
%       v3.4 10.04.2020   T.Reynaud temporal ==> time
%       v3.5 09.07.2021   DRIFT/DRIFT2/SLOPE extracted from Work (T.Reynaud)
%                         and written in netcdf files
%       v3.6 11.02.2022   Modified for Multi REF use (T.Reynaud)
%       v3.7 01.09.2022   Modifications added for handling Automatic Pressure corrections in R/D files. 
%       v3.7 05.09.2022   Modifications added for handling Automatic Pressure corrections in R/D files. 
%       v3.8 26.04.2024   Modifications added for handling PWLF:
%                         kwrite = 1 ==> SLOPE(1) and DRIFT(1) all cases+PWLF linear segment 1 
%                         kwrite = 2 ==> SLOPE(2) and DRIFT(2) PWLF linear segment 2
%      


function [argo, DIMD] = DOXY_update_fields(Work,argo,DIM,REF_ARGO,kwrite)

% =========================================================================
%% General
% =========================================================================
argo.date_update.data = datestr(now,'yyyymmddHHMMSS');
[argo.data_mode.data(:)] = 'D';
if isfield(argo,'n_prof')
    n_prof = argo.n_prof;
    argo.data_state_indicator.data = repmat('2C  ',argo.n_prof,1);
else
    n_prof = size(argo.pres.data,1);
    argo.data_state_indicator.data = repmat('2C  ',n_prof,1);
end

% Update parameter_data_mode
for n = 1:n_prof
    % check only the first n_calib dimension
    myparam = cellstr(permute(argo.parameter.data(n,1,:,:),[3,4,2,1]));
    isdoxy = strcmp(myparam,'DOXY');
    argo.parameter_data_mode.data(n,isdoxy) = 'D';
end

% =========================================================================
%% *Scientific calibration fields*
% =========================================================================
% -------------------------------------------------------------------------
%% Updating
% -------------------------------------------------------------------------
idxNewCalib = 1;
allfields = fieldnames(argo);

ii = strfind(allfields,'scientific_calib_');
is_calib = find(~cellfun('isempty',ii));

argo = check_FirstDimArray(argo,'N_CALIB');
DIMD = DIM;
DIMD.n_param.dimlength = size(argo.scientific_calib_comment.data,2);
DIMD.n_calib.dimlength = size(argo.scientific_calib_comment.data,1);

% -------------------------------------------------------------------------
%% Build new calib table with one more dimension if not enough space to write everything

%if (Work.DODRIFT==true || Work.anteriorcorr==1) && DIMD.n_calib.dimlength==1 %marine
if (Work.anteriorcorr==1) && DIMD.n_calib.dimlength==1 %marine
    
    %Fill the added calib line with fillValue
    % -------------------------------------------------------------------------
    if DIM.n_calib.dimlength~=0
        [argo_ex,DIM_ex]=extract_profile_dim(argo,DIM,'N_CALIB',1);
        for ik = is_calib'
            oneChamp = allfields{ik};
            ii=argo_ex.(oneChamp).data~=argo_ex.(oneChamp).FillValue_;
            argo_ex.(oneChamp).data(ii)=argo_ex.(oneChamp).FillValue_;
        end
        [argo,DIMD] = cat_profile_dim(argo,argo_ex,DIMD,DIM_ex,'N_CALIB');
    else
        for ik = is_calib'
            oneChamp = allfields{ik};
            lengthfield = zeros(1,length(argo.(oneChamp).dim));
            lengthfield(1)=1;
            for tk=2:length(argo.(oneChamp).dim)
                lengthfield(tk) = DIM.(lower(argo.(oneChamp).dim{tk})).dimlength;
            end
            argo.(oneChamp).data = repmat(argo_ex.(oneChamp).FillValue_,lengthfield);
            DIMD.n_calib.dimlength = 1;
        end
    end
end


% waiting the files revision at DAC level, modify the scientific_calib
% flields for PRES, MOLAR_DOXY, and TCPHASE, the fields find in the LOPS
% pool of Bio-Argo

% -------------------------------------------------------------------------
%% Scientific calibration fields PRES
%% Modif C.Lagadec (13 septembre 2016) : les valeurs SCIENTIFIC... sont mises
%% a blanc pour PRES
% -------------------------------------------------------------------------
% Equation and comment
equation.PRES    = ' ';
coefficient.PRES = ' ';
date.PRES        = ' ';
comment.PRES     = 'Adjusted values are provided in the core profile file';

% -------------------------------------------------------------------------
%% Scientific calibration fields I-ARGO
% -------------------------------------------------------------------------
% Equation and comment
IARGO = {'MOLAR_DOXY','TPHASE_DOXY','CPHASE_DOXY','C1PHASE_DOXY','C2PHASE_DOXY',...
    'RPHASE_DOXY','TEMP_DOXY','PHASE_DELAY_DOXY','TEMP_DOXY2'};
for ia = 1:length(IARGO)
    equation.(IARGO{ia}) = 'not applicable';
    comment.(IARGO{ia}) = 'not_applicable';
    coefficient.(IARGO{ia}) = 'not applicable';
    date.(IARGO{ia}) = argo.date_update.data;
end

% -------------------------------------------------------------------------
%% Scientific calibration fields B-ARGO
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 1 :  If pressure effect correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Work.anteriorcorr==1 %|| Work.DODRIFT==true
    % Equation and comment
    %     if Work.anteriorcorr==1 && Work.DODRIFT==true
    %         if strcmp(Work.whichDrift,'NCEP')
    %             part1 = 'Pressure effect correction on DOXY; PPOX converted from DOXY; Time drift correction on PPOX (drift computed from NCEP data)';
    %             equation.DOXY = 'DOXY1=DOXY*(1+C*PRES/1000); PPOX = f(DOXY1); PPOX1=PPOX - drift_coef';
    %             name_doxy = 'PPOX1';
    %         else
    %             name_doxy='DOXY2';
    %             equation.DOXY= 'DOXY1=DOXY*(1+C*PRES/1000); DOXY2= DOXY1 - drift_coef';
    %             part1 = 'Pressure effect correction on DOXY; Time drift correction on DOXY (drift computed from  data)';
    %         end
    %     elseif Work.anteriorcorr==1
    name_doxy='DOXY1';
    equation.DOXY= strcat('DOXY1=DOXY*(1+C*PRES/1000)');
    part1 = 'Pressure effect correction on DOXY';
    %     elseif Work.DODRIFT==true
    %         if strcmp(Work.whichDrift,'NCEP')
    %             part1 = sprintf('PPOX converted from DOXY; Time drift correction on PPOX (drift computed from %s data)',Work.whichDrift);
    %             equation.DOXY = 'PPOX = f(DOXY); PPOX1=PPOX - drift_coef';
    %             name_doxy = 'PPOX1';
    %         else
    %             part1 = sprintf('Time drift correction on DOXY (drift computed from %s data)',Work.whichDrift);
    %             equation.DOXY= 'DOXY1=DOXY - drift_coef';
    %             name_doxy='DOXY1';
    %         end
    
    
    % end
    
    comment.DOXY=part1;
    comment.DOXY=[comment.DOXY blanks(256-length(comment.DOXY))];
    equation.DOXY=[equation.DOXY blanks(256-length(equation.DOXY))];
    
    % Coefficient
    %if Work.anteriorcorr==1 && Work.DODRIFT==false
    coefficient.DOXY = sprintf('C=%6.4f',str2double(Work.coeff_corr));
    %     elseif Work.DODRIFT==false
    %         coefficient.DOXY = '';
    %     elseif Work.DODRIFT==true
    %         if Work.anteriorcorr==1
    %            coefficient.DOXY = sprintf('C=%6.4f; drift_coef=%6.2f;',str2double(Work.coeff_corr),Work.drift_val(num_fic));
    %         else
    %            coefficient.DOXY = sprintf('drift_coef=%6.2f; ',Work.drift_val(num_fic));
    %         end
    %     end
    coefficient.DOXY=[coefficient.DOXY blanks(256-length(coefficient.DOXY))];
    date.DOXY = argo.date_update.data;
    
    equation.DOXY_None = 'none';
    coefficient.DOXY_None = 'none';
    comment.DOXY_None = 'Bad data; not adjustable';
    date.DOXY_None = argo.date_update.data;
    
    % Fill the new scientific calibration fields with updated values
    scientifield = {'equation','comment','coefficient','date'};
    VAROK = IARGO;
    VAROK{end+1} = 'DOXY';
    VAROK{end+1} = 'DOXY2';
    VAROK{end+1} = 'PRES';
    for n = 1:n_prof
        myparam = cellstr(permute(argo.parameter.data(1,:,:,n),[2 3 4 1]));
        for p = 1:length(myparam)
            param = myparam{p};
            if strcmp(param,'DOXY2') % !!! 2 DO float !!!
                continue
            elseif isempty(param) % no parameter
                continue
            elseif ~ismember(param,VAROK) % not a DOXY parameter (ex: DOWNWELLING IRRADIANCE)
                continue
            end
            isok = strcmp(myparam,param);
            % if no correction has been done because the whole profile is bad
            if strcmp(param,'DOXY')
                checkAllBad = all(all(Work.DOXY_QC == '4' | Work.DOXY_QC == '9'));
                if checkAllBad
                    param = 'DOXY_None';
                end
            end
            % complete scientific_calib fields
            for f = 1:length(scientifield)
                myStr = eval([scientifield{f} '.' param]);
                %argo.(['scientific_calib_' scientifield{f}]).data(idxNewCalib,isok,1:length(myStr),n) = myStr;%Commented by T.Reynaud 05.09.2022
                % Added by T. Reynaud 05.09.2022
                %if f==4 & ~isempty(deblank(myStr)) %commented T. Reynaud 09.03.2023
                if ~isempty(deblank(myStr))
                  argo.(['scientific_calib_' scientifield{f}]).data(idxNewCalib,isok,1:length(myStr),n) = myStr;
                end

            end
        end
    end
    
    idxNewCalib=idxNewCalib+1;
else
    name_doxy='DOXY';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART 2 : Drift and Usual Scientific calibration fields B-ARGO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Work.whichO2quantity,'PPOX')
    equation.DOXY = strcat('PPOX_DOXY = f(',name_doxy,'); PPOX_DOXY_ADJUSTED=OFFSET+(SLOPE*(1+DRIFT/100.*(profile_date_juld-launch_date_juld)/365+DRIFT2/100.*(profile_date_juld-launch_date_juld)**2/365**2)+INCLINE_T*TEMP)*PPOX_DOXY; DOXY_ADJUSTED=f(PPOX_DOXY_ADJUSTED)');
    %part1 = ['PPOX_DOXY converted from DOXY; Partial pressure corrected %s;', ' DOXY_ADJUSTED converted from PPOX_DOXY_ADJUSTED'];
    part1 = ['Adjustment done on PPOX_DOXY;'];
end

% % Equation and comment
% if ismember(Work.whichCorr,{'WOA','REF'})
%     if strcmp(Work.whichO2quantity,'PSAT')
%         equation.DOXY = strcat('PSAT=f(',name_doxy,'); PSAT_ADJUSTED=A*PSAT+B;DOXY_ADJUSTED=f(PSAT_ADJUSTED)');
%         part1 = ['Percent saturation corrected as a linear function of PSAT;',...
%             ' %s; PSAT converted from DOXY and DOXY_ADJUSTED',...
%             ' converted from PSAT_ADJUSTED'];
%     elseif strcmp(Work.whichO2quantity,'DOXY')
%         equation.DOXY = strcat('DOXY_ADJUSTED=A*',name_doxy,'+B');
%         part1 = ['Oxygen concentration corrected as a linear function of DOXY; %s'];
%     end
% elseif ismember(Work.whichCorr,{'INAIR'})
%     if strcmp(name_doxy,'PPOX1')
%         equation.DOXY = strcat('PPOX_ADJUSTED = A*PPOX1 + B; DOXY_ADJUSTED=f(PPOX_ADJUSTED)');
%         part1 = ['Partial pressure corrected as a linear function of PPOX %s;',...
%              ' DOXY_ADJUSTED converted from PPOX_ADJUSTED'];
%     else
%         equation.DOXY = strcat('PPOX = f(',name_doxy,'); PPOX_ADJUSTED = A*PPOX + B; DOXY_ADJUSTED=f(PPOX_ADJUSTED)');
%         part1 = ['PPOX converted from DOXY; Partial pressure corrected as a linear function of PPOX %s;',...
%              ' DOXY_ADJUSTED converted from PPOX_ADJUSTED'];
%     end
% end
if Work.DODRIFT==true
    if strcmp(Work.whichDrift,'WOA')
        if Work.driftondeeplevels == 1
            part2 = ['Temporal drift estimated from WOA climatology on deep levels'];
        else
            part2 = ['Temporal drift estimated from WOA climatology at surface'];
        end
    elseif strcmp(Work.whichDrift,'NCEP')
        part2 = ['Temporal drift estimated from NCEP data'];
    end
else
    part2='';
end


if strcmp(Work.whichCorr,'WOA')
    part4 = ['Gain based on comparison between float profiles and WOA climatology'];
elseif strcmp(Work.whichCorr,'REF')
    % Modified by T. Reynaud 11.02.2022
    cvar='';
    for k=1:size(REF_ARGO.refId,2)
        cvar=[cvar,REF_ARGO.refId{k},'/'];
    end
    %part4 = ['Gain based on comparison between Argo cycle ' num2str(REF_ARGO.cycle) ' and reference profile ', REF_ARGO.refId];
    part4 = ['Gain based on comparison between Argo cycle ' num2str(REF_ARGO.cycle) ' and reference profile ', cvar];
    clear k cvar;
elseif strcmp(Work.whichCorr,'INAIR')
    if Work.isokC==1
        part4 = 'Gain based on comparison between in-air measurements and NCEP, with carry over';
    else
        part4 = 'Gain based on comparison between in-air measurements and NCEP, no carry over';
    end
else
    part4='';
end
if ~isempty(part2) && ~isempty(part4)
    part3=' - ';
else
    part3='';
end

comment.DOXY =[ part1 part2 part3 part4];
comment.DOXY=[comment.DOXY blanks(256-length(comment.DOXY))];
equation.DOXY=[equation.DOXY blanks(256-length(equation.DOXY))];

% Coefficient
% if ~isfield(Work,['OFFSET_' Work.whichO2quantity])
%     Aformat = '%2.3f'; Bformat = '%2.3f';
%     if Work.(['FIT_slope_' Work.whichO2quantity]) == 0
%         Aformat = '%d';
%     end
%     if Work.(['FIT_intercept_' Work.whichO2quantity]) == 0
%         Bformat = '%d';
%     end
%     coefficient.DOXY = sprintf([['A=' Aformat] ';' [' B=' Bformat]],Work.(['FIT_slope_' Work.whichO2quantity]), Work.(['FIT_intercept_' Work.whichO2quantity]));
% else
%     coefficient.DOXY = sprintf('A=%d; B=%2.3f',0, Work.(['OFFSET_' Work.whichO2quantity]));
% end
if Work.DODRIFT==true
    %Work.drift_val(num_fic)
    if Work.drift_fitPolynomialDegree == 1
        a=Work.PPOX_DRIFTCORR_COEF(1);%at+b
        b=Work.PPOX_DRIFTCORR_COEF(2);
        c=0;%08.07.2021 added by T. Reynaud
    elseif Work.drift_fitPolynomialDegree == 2
        c=Work.PPOX_DRIFTCORR_COEF(1);%ct^2+at+b
        a=Work.PPOX_DRIFTCORR_COEF(2);
        b=Work.PPOX_DRIFTCORR_COEF(3);%08.07.2021 added by T. Reynaud
    else
        disp('Problem, drift equation not in the type ct^2+at+b!');
    end
else
    a=0;
    b=1;
    c=0;%added by T. Reynaud 09.07.2021
end
if ~isfield(Work,['OFFSET_' Work.whichO2quantity])
    G=Work.(['FIT_slope_' Work.whichO2quantity]);
    OFFSET=Work.(['FIT_intercept_' Work.whichO2quantity]);
else
    G=1;
    OFFSET=Work.(['OFFSET_' Work.whichO2quantity]);
end

% Modified TR 09/07/2021
%SLOPE=b*G;
%DRIFT=100*365*G*a/SLOPE;
OFFSET=Work.OFFSET;
INCLINE_T=Work.INCLINE_T;
SLOPE=Work.SLOPE;
DRIFT=Work.DRIFT;
DRIFT2=Work.DRIFT2;
% End TR 09/07/2021

%coefficient.DOXY = sprintf('INCLINE_T=%d; SLOPE=%2.6f, DRIFT=%2.6f, OFFSET=%2.6f',0,SLOPE,DRIFT,OFFSET);
% Modified by T. Reynaud 09.07.2021
ctmp1=sprintf('INCLINE_T=%d, ',INCLINE_T);

ctmp2=sprintf('SLOPE=%2.6f, ',SLOPE(kwrite));
ctmp3=sprintf('DRIFT=%2.6f, ',DRIFT(kwrite));

ctmp4=sprintf('OFFSET=%2.6f, ',OFFSET);
ctmp5=sprintf('DRIFT2=%4.6f',DRIFT2);

% Commented TR 24.04.2024
%coefficient.DOXY = sprintf('INCLINE_T=%d; SLOPE=%2.6f, DRIFT=%2.6f, OFFSET=%2.6f, DRIFT2=%4.6f',0,SLOPE,DRIFT,OFFSET,DRIFT2);
coefficient.DOXY=[ctmp1,ctmp2,ctmp3,ctmp4,ctmp5];
coefficient.DOXY=[coefficient.DOXY blanks(256-length(coefficient.DOXY))];


date.DOXY = argo.date_update.data;

equation.DOXY_None = 'none';
coefficient.DOXY_None = 'none';
comment.DOXY_None = 'Bad data; not adjustable';
date.DOXY_None = argo.date_update.data;

% Fill the new scientific calibration fields with updated values
scientifield = {'equation','comment','coefficient','date'};
VAROK = IARGO;
VAROK{end+1} = 'DOXY';
VAROK{end+1} = 'DOXY2';
VAROK{end+1} = 'PRES';
for n = 1:n_prof
    myparam = cellstr(permute(argo.parameter.data(1,:,:,n),[2 3 4 1]));
    for p = 1:length(myparam)
        param = myparam{p};
        if strcmp(param,'DOXY2') % !!! 2 DO float !!!
            continue
        elseif isempty(param)   % no parameter
            continue
        elseif ~ismember(param,VAROK)  % not a DOXY parameter (ex: DOWNWELLING IRRADIANCE)
            continue
        end
        isok = strcmp(myparam,param);
        % if no correction has been done because the whole profile is bad
        if strcmp(param,'DOXY')
            checkAllBad = all(all(Work.DOXY_QC == '4' | Work.DOXY_QC == '9'));
            if checkAllBad
                param = 'DOXY_None';
            end
        end
        % complete scientific_calib fields
        for f = 1:length(scientifield)
            %Marine
            %             if Work.DODRIFT==true && strcmp(param,'DOXY') && f==3
            %                 for i=1:length(Work.drift_val)
            %                     myStr = eval([scientifield{f} '.' param]);
            %                     argo.(['scientific_calib_' scientifield{f}]).(['data' num2str(i)])(idxNewCalib,isok,1:length(myStr),n) = myStr;
            %                 end
            %             else
            myStr = eval([scientifield{f} '.' param]);
            %argo.(['scientific_calib_' scientifield{f}]).data(idxNewCalib,isok,1:length(myStr),n) =myStr;% Commented TR 01/09/2022
            % Test on isempty myStr added by T. Reynaud 01/09/2022 ==> date only

            %if f==4 & ~isempty(deblank(myStr)) %Commented T. Reynaud 09.03.2023         
            if ~isempty(deblank(myStr))
                argo.(['scientific_calib_' scientifield{f}]).data(idxNewCalib,isok,1:length(myStr),n) = myStr;
            end
            %end
            
        end
    end
end



% =========================================================================
%% History fields
% =========================================================================
idxNewHistory = DIM.n_history.dimlength+1;
allfields = fieldnames(argo);

ii = strfind(allfields,'history_');
is_history = find(~cellfun('isempty',ii));

argo = check_FirstDimArray(argo,'N_HISTORY');
DIMD = DIMD;
DIMD.n_param.dimlength = size(argo.scientific_calib_comment.data,2);
DIMD.n_calib.dimlength = size(argo.scientific_calib_comment.data,1);

% Build new history table with one more dimension. Fill the added history
% line with fillValue
if DIM.n_history.dimlength~=0
    [argo_ex,DIM_ex]=extract_profile_dim(argo,DIMD,'N_HISTORY',1);
    for ik = is_history'
        oneChamp = allfields{ik};
        ii=argo_ex.(oneChamp).data~=argo_ex.(oneChamp).FillValue_;
        argo_ex.(oneChamp).data(ii)=argo_ex.(oneChamp).FillValue_;
    end
    
    [argo,DIMD] = cat_profile_dim(argo,argo_ex,DIMD,DIM_ex,'N_HISTORY');
else
    for ik = is_history'
        oneChamp =allfields{ik};
        lengthfield = zeros(1,length(argo.(oneChamp).dim));
        lengthfield(1)=1;
        for tk=2:length(argo.(oneChamp).dim)
            lengthfield(tk) = DIM.(lower(argo.(oneChamp).dim{tk})).dimlength;
        end
        argo.(oneChamp).data = repmat(argo_ex.(oneChamp).FillValue_,lengthfield);
        DIMD.n_history.dimlength=1;
    end
end

% Fill the new history fields with updated values
for n = 1:n_prof
    %     tmp = cellstr(squeeze(argo.history_institution.data(:,n,:,:)));
    %     idxNewHistory = find(cellfun('isempty',tmp),1,'First');
    myStr = 'IF';
    argo.history_institution.data(idxNewHistory,n,1:length(myStr)) = myStr;
    myStr = 'CV';
    argo.history_action.data(idxNewHistory,n,1:length(myStr)) = myStr;
    myStr = Work.history.history_software_release;
    argo.history_software_release.data(idxNewHistory,n,1:length(myStr)) = myStr;
    myStr = 'ARSQ';
    argo.history_step.data(idxNewHistory,n,1:length(myStr)) = myStr;
    myStr = Work.history.history_software(1:4);
    argo.history_software.data(idxNewHistory,n,1:length(myStr)) = myStr;
    myStr = Work.history.history_reference;
    argo.history_reference.data(idxNewHistory,n,1:length(myStr)) = myStr;
    myStr = 'DOXY';
    argo.history_parameter.data(idxNewHistory,n,1:length(myStr)) = myStr;
    argo.history_date.data(idxNewHistory,n,:) = argo.date_update.data;
    % get first and last non fillvalue pressure
    isok = (argo.pres.data(n,:) ~= argo.pres.FillValue_);
    if any(isok)
        argo.history_start_pres.data(idxNewHistory,n) = argo.pres.data(n,find(isok,1,'first'));
        argo.history_stop_pres.data(idxNewHistory,n) = argo.pres.data(n,find(isok,1,'last'));
    end
end


argo = check_FirstDimArray(argo,'N_PROF');

clear date;

