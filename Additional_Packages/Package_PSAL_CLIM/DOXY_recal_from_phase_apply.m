function O2=DOXY_recal_from_phase_apply(coef,equations,optode,PRES,TEMP,PSAL,C1PHASE_DOXY,C2PHASE_DOXY,TEMP_DOXY,MOLAR_DOXY,TPHASE_DOXY,ip)

% Attention aux unités: mM/L ... faut diviser par rho.

% DESCRIPTION
% Computes O2 in mM/L
% 
%
% INPUT
%     coef (structure)         coefficient values from Argo files
%     coef.L0: [6×4096 char] ==> for TEMP_DOXY
%     coef.L1: [24×4096 char]==> for water values
%     coef.L2: [17×4096 char]==> for air values
%
%     equatioıs (structure)    equations values from Argo files
%     equations.L0: 'TEMP_DOXY=T0+T1*TEMP_VOLTAGE_DOXY+T2*TEMP_VOLTAGE_DOXY^2+T3*TEMP_VOLTAGE_DOXY^3+T4*TEMP_VOLTAGE_DOXY^4+T5*TEMP_VOLTAGE_DOXY^5
%     equations.L1: [11×4096 char]
%     equations.L2: [11×4096 char]

%     PRES: pressure from CTD
%     TEMP: temperature from CTD
%     PSAL: salinity from CTD
%     C1PHASE_DOXY: from optode
%     C2PHASE_DOXY: from optode
%     TEMP_DOXY: temperature from the optode
%     ip : measures type (0 for TEMP_DOXY, 1 for water and 2 for air)

%
% OUTPUT
%     O2 : Oxygen content ïn mM/Kg

%
% HISTORY
%   $created: 26/05/2015 $author: T. Reynaud
%   $Revision: version $Date: $author:
%   v1 10/01/2024   T. Reynaud
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean variables

if strcmp(optode,'4330-01')
    clear MOLAR_DOXY TPHASE_DOXY;
elseif strcmp(optode,'4330-02')
    clear C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY MOLAR_DOXY;
elseif strcmp(optode,'3830-01')
    clear C1PHASE_DOXY C2PHASE_DOXY TEMP_DOXY TPHASE_DOXY;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read coeficients
eval(strcat('tmp=coef.L',num2str(ip),';'));

for it=1:size(tmp,1)
    eval(strcat(deblank(tmp(it,:)),';'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply equations

eval(strcat('tmp2=equations.L',num2str(ip),';'));
i_exe=ones(size(tmp2,1),1);

while sum(i_exe)
    for ieq=1:size(tmp2,1)
        if i_exe(ieq)
            line=strcat(deblank(tmp2(ieq,:)));

            % Build the function or variable name
            % Split LHS and RHS part of the equations
            inam=strfind(line,'=');
            gauche=line(1:inam-1);
            droite=line(inam+1:end);

            % Determine input arguments on the RHS
            % clean up everything which is not a variable
            droite=strrep(droite,'exp(',' ');
            droite=strrep(droite,'ln(',' ');
            droite=strrep(droite,'-',' ');
            droite=strrep(droite,'+',' ');
            droite=strrep(droite,'/',' ');
            droite=strrep(droite,'*',' ');
            droite=strrep(droite,'(',' ');
            droite=strrep(droite,')',' ');
            droite=strrep(droite,'^2',' ');
            droite=strrep(droite,'^3',' ');
            droite=strrep(droite,'^4',' ');
            droite=strrep(droite,'1000',' ');
            droite=strrep(droite,'1013.25',' ');
            droite=strrep(droite,'100',' ');
            droite=strrep(droite,'298.15',' ');
            droite=strrep(droite,'273.15',' ');
            droite=deblank(droite);
            droite=strrep(droite,'  ',' ');
            droite=strrep(droite,'  ',' ');
            droite=strrep(droite,'  ',' ');
            droite=fliplr(droite);
            droite=deblank(droite);
            droite=fliplr(droite);
            droite=strrep(droite,'  ',' ');
            droite=strrep(droite,'   ',' ');
            droite=strrep(droite,',',' ');

            toto=split(droite);

            %Create function pH2O
            if strcmp(gauche,'pH2O(TEMP,S)')

                cfunc='function ';
                line_return='return';
                line_end='end';

                % Create tmp function directory
                func_rep=strcat(pwd,'/tmp_func/');

                save('D.mat','D0','D1','D2','D3');
                func_name=gauche;
                ibeg=strfind(func_name,'(')-1;
                iend=length(func_name);
                func_name=deblank(func_name(1:ibeg));
                wline1=strcat(cfunc,' out=',func_name,'(TEMP,S)');

                wline2=strcat('load("D.mat")');
                wline3=strcat('out=',line(inam+1:end),';');
                wline3=strrep(wline3,'ln','log');
                wline3=strrep(wline3,'/','./');
                wline4=line_return;
                wline5=line_end;

                filename=fullfile(func_rep,strcat(func_name,'.m'));
                if exist(filename,'file')
                    delete(filename);
                    rmpath(func_rep);
                end


                if exist(func_rep,'dir')
                    rmdir(func_rep);
                end
                mkdir(func_rep);
                eval(fullfile('addpath ',func_rep))
                %

                fid = fopen(filename,'w');

                fprintf(fid,'%s\n',wline1);
                fprintf(fid,'%s\n',wline2);
                fprintf(fid,'%s\n',wline3);
                fprintf(fid,'%s\n',wline4);
                fprintf(fid,'%s\n',wline5);
                fclose(fid);
                i_exe(ieq)=0;
            end

            iok=true;
            for ivar=1:size(toto,1)
                if ~exist(toto{ivar}) && ~strcmp(toto{ivar},'1')
                    iok=false;
                end
            end
            if iok
                line=strrep(line,'^','.^');
                line=strrep(line,'ln','log');
                line=strrep(line,'/','./');
                line=strrep(line,'*','.*');
                display(strcat(line,';'));
                eval(strcat(line,';'));
                i_exe(ieq)=0;% When run set to 0.
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return
end