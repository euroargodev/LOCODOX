function tmp=DOXY_recal_from_phase_comment(file_doxy_full,Work,ip_comment)
% Prepare SCIENTIFIC_CALIB_COMMENT for BR/BD files
% Last date is modified

% ip_comment: profile number 1 (water) or 2 (air)
% file_doxy_full: Complete BD/BR filename including path
% Work:   contains the Climatological information

% Dimensions information
% char PARAMETER(N_PROF, N_CALIB, N_PARAM, STRING64) ;
% char SCIENTIFIC_CALIB_EQUATION(N_PROF, N_CALIB, N_PARAM, STRING256) ;
% char SCIENTIFIC_CALIB_COEFFICIENT(N_PROF, N_CALIB, N_PARAM, STRING256) ;
% char SCIENTIFIC_CALIB_COMMENT(N_PROF, N_CALIB, N_PARAM, STRING256) ;
% char SCIENTIFIC_CALIB_DATE(N_PROF, N_CALIB, N_PARAM, DATE_TIME) ;

%	N_PROF = 2 ;
%	N_CALIB = 2 ;
%	N_PARAM = 3 ;
%	STRING256 = 256 ;
%	STRING64 = 64 ;


tmp.param=ncread(file_doxy_full,'PARAMETER');% 64x3x2x2
tmp.equation=ncread(file_doxy_full,'SCIENTIFIC_CALIB_EQUATION');
tmp.coefficient=ncread(file_doxy_full,'SCIENTIFIC_CALIB_COEFFICIENT');
tmp.comment=ncread(file_doxy_full,'SCIENTIFIC_CALIB_COMMENT');
tmp.date=ncread(file_doxy_full,'SCIENTIFIC_CALIB_DATE');


% 	char PARAMETER(N_PROF, N_CALIB, N_PARAM, STRING64) ;


[STRING64, N_PARAM, N_CALIB, N_PROF]=size(tmp.param);

csearch='DOXY';
for ip=1:N_PARAM
    for ic=1:N_CALIB
        for i_prof=1:N_PROF
            vcar=deblank(squeeze(tmp.param(:,ip,ic,i_prof)'));
            vdate=deblank(squeeze(tmp.date(:,ip,ic,i_prof)'));
            veq=deblank(squeeze(tmp.equation(:,ip,ic,i_prof)'));
            if strcmp(vcar,csearch) && i_prof==ip_comment && ic==N_CALIB
                vcomment=blanks(256);
                tempo=strcat('Recovered BGC profile.',blanks(1),Work.PSAL_REPLACE_CLIM,' product salinity used.');
                vcomment(1:length(tempo))=tempo;
                tmp.comment(:,ip,ic,i_prof)=vcomment;
            end
            vcomment=deblank(squeeze(tmp.comment(:,ip,ic,i_prof)'));
            vcoef=deblank(squeeze(tmp.coefficient(:,ip,ic,i_prof)'));

            %display(vcar);
            if strcmp(vcar,csearch)
                display('=============');
                display(strcat(vcar,'==> ',' nparam=',num2str(ip),' calib=',num2str(ic),' iprof=',num2str(i_prof)));
                display(veq);
                display(vcoef);
                display(vcomment);
                display(vdate);
            end
        end
    end
end

return
end

