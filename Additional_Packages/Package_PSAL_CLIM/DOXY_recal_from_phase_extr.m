function [coef,equations,tech]=DOXY_recal_from_phase_extr(cwmo,rep_float)

% Recalculate DOXY
% extraction of equations and coefficients.

dir_meta=strcat(rep_float,'*meta.nc');
file_list=dir(char(dir_meta));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get SENSOR TYPE information:
tmpbis=ncread(strcat(file_list.folder,'/',file_list.name),'SENSOR');
tmpbis=tmpbis';
idx=0;
for im=1:size(tmpbis,1)
    if ~isempty(strfind(tmpbis(im,:),'DOXY'))
        idx=im;
    end
end
clear tmpbis;

tmpter=ncread(strcat(file_list.folder,'/',file_list.name),'SENSOR_MODEL');
tmpter=tmpter';

tech.DOXY_SENSOR=deblank(tmpter(idx,:));

clear tmpter idx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract coefficients
tmp=ncread(strcat(file_list.folder,'/',file_list.name),'PREDEPLOYMENT_CALIB_COEFFICIENT');
tmp=tmp';

% Clean empty lines
tmp_clean=[];
for iline=1:size(tmp,1)
    cline=tmp(iline,:);
    if strfind(cline,'=')
            tmp_clean=[tmp_clean;cline];
    end
end
tmp=tmp_clean;
clear tmp_clean;


for iline=1:size(tmp,1)
    coef_all=[];
    cline=tmp(iline,:);
    cline=strcat(cline,';');
    cline=strrep(cline,';;',';');
    cline=strrep(cline,'not available','NaN');
    cline=strrep(cline,',',';');
    idx=strfind(cline,';');
    if ~isempty(idx)
        ibeg=1;
        for ieq=1:length(idx)
            cf1=blanks(size(tmp,2));
            iend=idx(ieq)-1;
            cf1(1:iend-ibeg+1)=cline(ibeg:iend);
            if strfind(cf1,'=') 
                coef_all=[coef_all;cf1];
            end
            ibeg=iend+3;
        end
    end
    eval(strcat('coef.L',num2str(iline-1),'=coef_all;'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extract equations from tmp
tmp=ncread(strcat(file_list.folder,'/',file_list.name),'PREDEPLOYMENT_CALIB_EQUATION');
tmp=tmp';

% Clean empty lines
tmp_clean=[];
for iline=1:size(tmp,1)
    cline=tmp(iline,:);
    if strfind(cline,'=')
            tmp_clean=[tmp_clean;cline];
    end
end
tmp=tmp_clean;
clear tmp_clean;


% TEMP_DOXY : read from BR or BD files

for iline=1:size(tmp,1)
    equations_all=[];
    cline=tmp(iline,:);
    %cline=strrep(cline,'with',';with');
    cline=strcat(cline,';');
    cline=strrep(cline,';;',';');
    idx=strfind(cline,';');
    cline=strrep(cline,'[','(');
    cline=strrep(cline,']',')');
    if ~isempty(idx)
        ibeg=1;
        for ieq=1:length(idx)
            eq1=blanks(size(tmp,2));
            iend=idx(ieq)-1;
            eq1(1:iend-ibeg+1)=cline(ibeg:iend);
            if ~isempty(strfind(eq1,'=')) && isempty(strfind(eq1,'with')) && isempty(strfind(eq1,'rho'))
                equations_all=[equations_all;eq1];
            end
            ibeg=iend+3;
        end
    end
    eval(strcat('equations.L',num2str(iline-1),'=equations_all;'));
end
return
end
