close all;
clear all;

addpath /Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX/share/seawater/seawater_330_its90_lpo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cwmo='6903083';
% rep_float=strcat('/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_FLOAT_DATA/DMQC_PSAL/coriolis/',cwmo,'/');
%
% [coef,equations]=DOXY_recal_from_phase_extr(cwmo,rep_float);% O2 mM/L
% %L1 DOXY in Water
% %L2 PPOX_DOXY in the air
%
% rep_float=strcat(rep_float,'profiles/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%cwmo='6902882';%4330
%cwmo='6902818';%4330--Bad PSAL
%cwmo='6902806';%4330
%cwmo='6902807';%4330
%cwmo='6902802';%4330
%cwmo='6902686';%4330 ==> Faut lire TPHASE_DOXY. Programme a revoir.
%cwmo='6902800';%4330
%rep_float=strcat('/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_FLOAT_DATA/For_VT/DMQC_EXPORT/coriolis/',cwmo,'/');


cwmo='5902297';%3830
rep_float=strcat('/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_EXTERNAL_FLOAT_DATA/CORIOLIS/coriolis/',cwmo,'/');

[coef,equations,tech]=DOXY_recal_from_phase_extr(cwmo,rep_float);% O2 mM/L
%L1 DOXY in Water
%L2 PPOX_DOXY in the air

rep_float=strcat(rep_float,'profiles/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
dir_files=strcat(rep_float,'B*.nc');
file_list=dir(char(dir_files));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare interpolation parameters

% Work.PSAL_REPLACE_CLIM='ISAS';
% Work.PSAL_REPLACE_DIR='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/ISAS/';
% Work.PSAL_REPLACE_CLIM_file='isas17.mat';

% Work.PSAL_REPLACE_CLIM='WOA';
% Work.PSAL_REPLACE_DIR='/Users/treynaud/IFREMER/MATLAB/LOCODOX/LOCODOX_LOPS_DATA/WOA/';
% Work.PSAL_REPLACE_CLIM_file='woa18.mat';

Work.PSAL_REPLACE_CLIM='ARMOR-3D';
Work.PSAL_REPLACE_DIR='/Users/treynaud/IFREMER/MATLAB/LOCODOX/Tools_treynaud/';
Work.PSAL_REPLACE_CLIM_file='copernicus';

Work.PSAL_REPLACE_plot=true;
Work.PSAL_REPLACE_plot_close=false;
Work.PSAL_REPLACE_plot_dir=fullfile('/Users/treynaud/IFREMER/MATLAB/LOCODOX/Tools_treynaud/plots/',cwmo,'/');
WOA=[];
CLIM=[];

% Create output directory
if ~exist(Work.PSAL_REPLACE_plot_dir,'dir')
    mkdir(Work.PSAL_REPLACE_plot_dir);
end

ptype='primary';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables
global_DOXY_ORI=[];
global_DOXY_REC2=[];
global_jdays=[];
global_mean_DOXY_diff_ori=[];
global_mean_DOXY_diff=[];
global_mean_DOXY_diff_1000=[];
global_mean_DOXY_diff_2000=[];
global_mean_DOXY_diff_6000=[];
global_fcolor=[];
global_num=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1:length(file_list)

imin=1;
imax=length(file_list);
%
% imin=123;
% imax=123;

for i=imin:imax
    %for i=1:1
    file_doxy=file_list(i).name;

    if ~exist('file_first','var')
        file_first=file_doxy
    end
    file_doxy_full=strcat(file_list(i).folder,'/',file_doxy);


    %Initialized variables:
    MOLAR_DOXY=[];
    C1PHASE_DOXY=[];
    C2PHASE_DOXY=[];
    TEMP_DOXY=[];
    TPHASE_DOXY=[];

    toto=ncinfo(file_doxy_full);
    for ivar=1:length(toto.Variables)
        %display(toto.Variables(ivar).Name)
        if ~isempty(strfind(tech.DOXY_SENSOR,'4330')) && strcmp(toto.Variables(ivar).Name,'C1PHASE_DOXY')
            optode='4330-01';
            C1PHASE_DOXY=ncread(file_doxy_full,'C1PHASE_DOXY');
            C2PHASE_DOXY=ncread(file_doxy_full,'C2PHASE_DOXY');
            TEMP_DOXY=ncread(file_doxy_full,'TEMP_DOXY');
            ieq_water=1;
        elseif ~isempty(strfind(tech.DOXY_SENSOR,'4330')) && strcmp(toto.Variables(ivar).Name,'TPHASE_DOXY')
            optode='4330-02';
            TPHASE_DOXY=ncread(file_doxy_full,'TPHASE_DOXY');
            ieq_water=0;
        elseif ~isempty(strfind(tech.DOXY_SENSOR,'3830')) && strcmp(toto.Variables(ivar).Name,'MOLAR_DOXY')
            optode='3830-01';
            MOLAR_DOXY=ncread(file_doxy_full,'MOLAR_DOXY');
            ieq_water=0;
        end
    end

    PRES_DOXY=ncread(file_doxy_full,'PRES');
    DOXY=ncread(file_doxy_full,'DOXY');
    QC_DOXY=ncread(file_doxy_full,'DOXY_QC');

    cycle_float=i;
    Work.PSAL_REPLACE_cycle_beg=cycle_float;
    Work.PSAL_REPLACE_cycle_end=cycle_float;

    lon_float=ncread(file_doxy_full,'LONGITUDE');
    lat_float=ncread(file_doxy_full,'LATITUDE');

    if  ~isempty(strfind(file_doxy,'BR'))
        file_ctd=strrep(file_doxy,'BR','D');
        file_ctd_full=strcat(file_list(i).folder,'/',file_ctd);
        if ~exist(file_ctd_full,'file')
            file_ctd=strrep(file_doxy,'BR','R');
            file_ctd_full=strcat(file_list(i).folder,'/',file_ctd);
        end
    else
        file_ctd=strrep(file_doxy,'BD','D');
        file_ctd_full=strcat(file_list(i).folder,'/',file_ctd);
        if ~exist(file_ctd_full,'file')
            file_ctd=strrep(file_doxy,'BD','R');
            file_ctd_full=strcat(file_list(i).folder,'/',file_ctd);
        end
    end

    TEMP_ctd=ncread(file_ctd_full,'TEMP');
    PSAL_ctd=ncread(file_ctd_full,'PSAL');
    PSAL_adjusted_ctd=ncread(file_ctd_full,'PSAL_ADJUSTED');
    PRES_ctd=ncread(file_ctd_full,'PRES');

    QC_PSAL_ctd=ncread(file_ctd_full,'PSAL_QC');
    QC_PSAL_adjusted_ctd=ncread(file_ctd_full,'PSAL_ADJUSTED_QC');

    juld_float=ncread(file_ctd_full,'JULD');
    scale_units = ncreadatt(file_ctd_full,'JULD','units');
    dateref=sscanf(scale_units,'days since %f-%f-%f %f:%f:%f UTC');
    julref=datenum(dateref');


    %    for ip=1:size(PRES_DOXY,2)
    for ip=1:1
        PRES=PRES_DOXY(:,ip);
        TEMP=TEMP_ctd(:,ip);
        PSAL=PSAL_ctd(:,ip);

        if strcmp(optode,'4330-01')
            C1PHASE=C1PHASE_DOXY(:,ip);
            C2PHASE=C2PHASE_DOXY(:,ip);
            TEMP_DOXY1=TEMP_DOXY(:,ip);
            MOLAR_DOXY1=[];
            TPHASE_DOXY1=[];
        elseif strcmp(optode,'4330-02')
            C1PHASE=[];
            C2PHASE=[];
            TEMP_DOXY1=[];
            MOLAR_DOXY1=[];
            TPHASE_DOXY1=TPHASE_DOXY(:,ip);
        elseif strcmp(optode,'3830-01')
            MOLAR_DOXY1=MOLAR_DOXY(:,ip);
            TEMP_DOXY1=[];
            C1PHASE=[];
            C2PHASE=[];
            TPHASE_DOXY1=[];            
        end

        DOXY_ORI=DOXY(:,ip);

        O2=DOXY_recal_from_phase_apply(coef,equations,optode,PRES,TEMP,PSAL,C1PHASE,C2PHASE,TEMP_DOXY1,MOLAR_DOXY1,TPHASE_DOXY1,ieq_water);

        PR=0;
        PTEMP = sw_ptmp(PSAL,TEMP,PRES,PR);
        dens0 = sw_dens0(PSAL,PTEMP);
        DOXY_REC=O2./(dens0/1000);% Convert units in Mm/Kg

        %dens = sw_dens(PSAL,PTEMP,PRES);
        %DOXY_REC_INSITU=O2./(dens/1000);% Convert units in Mm/Kg

        tmp2=PRES_ctd(:,ip);
        tmp2=tmp2';
        eval(['DATA.',ptype,'.argo.pres.data=tmp2',';']);

        tmp2=PSAL_ctd(:,ip);
        tmp2=tmp2';
        eval(['DATA.',ptype,'.argo.psal.data=tmp2',';']);


        tmp2=PSAL_adjusted_ctd(:,ip);
        tmp2=tmp2';
        eval(['DATA.',ptype,'.argo.psal_adjusted.data=tmp2',';']);


        tmp2=TEMP_ctd(:,ip);
        tmp2=tmp2';
        eval(['DATA.',ptype,'.argo.temp.data=tmp2',';']);

        tmp2=cycle_float;
        tmp2=tmp2';
        eval(['DATA.',ptype,'.argo.cycle_number.data=tmp2',';']);

        tmp2=QC_PSAL_adjusted_ctd(:,ip:ip);
        tmp2=tmp2';
        eval(['DATA.',ptype,'.argo.psal_adjusted_qc.data=tmp2',';']);

        tmp2=QC_PSAL_ctd(:,ip:ip);
        tmp2=tmp2';
        eval(['DATA.',ptype,'.argo.psal_qc.data=tmp2',';']);


        eval(['DATA.',ptype,'.argo.longitude.data=lon_float(ip);']);
        eval(['DATA.',ptype,'.argo.latitude.data=lat_float(ip);']);
        % Float Julian Days
        tmp=datestr(julref,'yyyymmddHHMMSS');
        eval(['DATA.',ptype,'.argo.reference_date_time.data=tmp;']);
        eval(['DATA.',ptype,'.argo.juld.data=juld_float(ip);']);

        % Climatological interpolation
        if ~strcmp(Work.PSAL_REPLACE_CLIM,'ARMOR-3D')
            [DATA,out_std,CLIM]=DOXY_PSAL_clim_replace_main(ptype,DATA,WOA,Work,CLIM);
        else
            [DATA,out_std,CLIM]=DOXY_PSAL_clim_A3D_replace_main(ptype,DATA,WOA,Work,CLIM);
        end


        PRES2=DATA.primary.argo.pres.data';
        TEMP2=DATA.primary.argo.temp.data';
        PSAL2=DATA.primary.argo.psal.data';
        O2=DOXY_recal_from_phase_apply(coef,equations,optode,PRES2,TEMP2,PSAL2,C1PHASE,C2PHASE,TEMP_DOXY1,MOLAR_DOXY1,TPHASE_DOXY1,ieq_water);

        PR=0;
        PTEMP2 = sw_ptmp(PSAL2,TEMP2,PRES2,PR);
        dens0 = sw_dens0(PSAL2,PTEMP2);

        DOXY_REC2=O2./(dens0/1000);% Convert units in Mm/Kg


        if isfield(out_std,'std_clim_psal_float')
            err=out_std.std_clim_psal_float';
        else
            err=1/4;
        end

        % 3 sigma ==> 99%
        err=3*err;

        O2=DOXY_recal_from_phase_apply(coef,equations,optode,PRES2,TEMP2,PSAL2+err,C1PHASE,C2PHASE,TEMP_DOXY1,MOLAR_DOXY1,TPHASE_DOXY1,ieq_water);

        PR=0;
        PTEMP2 = sw_ptmp(PSAL2+err,TEMP2,PRES2,PR);
        dens0 = sw_dens0(PSAL2+err,PTEMP2);

        DOXY_REC2PE=O2./(dens0/1000);% Convert units in Mm/Kg

        err=-err;
        O2=DOXY_recal_from_phase_apply(coef,equations,optode,PRES2,TEMP2,PSAL2+err,C1PHASE,C2PHASE,TEMP_DOXY1,MOLAR_DOXY1,TPHASE_DOXY1,ieq_water);

        PR=0;
        PTEMP2 = sw_ptmp(PSAL2+err,TEMP2,PRES2,PR);
        dens0 = sw_dens0(PSAL2+err,PTEMP2);

        DOXY_REC2ME=O2./(dens0/1000);% Convert units in Mm/Kg

        DOXY_ORI(DOXY_ORI<0)=NaN;
        DOXY_REC2(DOXY_REC2<0)=NaN;
        DOXY_REC2PE(DOXY_REC2PE<0)=NaN;
        DOXY_REC2ME(DOXY_REC2ME<0)=NaN;

        DOXY_CRIT_MAX=600;
        DOXY_ORI(DOXY_ORI>DOXY_CRIT_MAX)=NaN;
        DOXY_REC2(DOXY_REC2>DOXY_CRIT_MAX)=NaN;
        DOXY_REC2PE(DOXY_REC2PE>DOXY_CRIT_MAX)=NaN;
        DOXY_REC2ME(DOXY_REC2ME>DOXY_CRIT_MAX)=NaN;

        if Work.PSAL_REPLACE_plot
            fig=figure(99);
            subplot(1,3,1)
            title(strrep(file_doxy,'_','\_'));
            subplot(1,3,2)
            title(strrep(file_doxy,'_','\_'));

            subplot(1,3,3)

            h1=plot(DOXY_REC2,-PRES2,'r');
            hold on;
            h2=plot(DOXY_ORI,-PRES,'g');

            h3=plot(DOXY_REC2PE,-PRES2,'.m');
            h4=plot(DOXY_REC2ME,-PRES2,'.c');

            %title(strrep(file_doxy,'_','\_'));
            xlabel('DOXY (MM/Kg)');
            ylabel('Pressure (Db)');
            hleg=legend([h1,h2,h3,h4],['Clim ',num2str(nanmean(DOXY_REC2-DOXY_ORI),'%7.2f')],'Float',['3s=> ', num2str(nanmean(abs(err)),'%6.2f')],['3s => ', num2str(nanmean(-abs(err)),'%6.2f')],'Location','SouthEast');
            %psal_float (:,ilat)= squeeze(tmp);
            %pause;

            hAx(1)=gca;
            hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
            hold(hAx(2),'on')
            h1bis=plot(DOXY_REC2-DOXY_ORI, -PRES,'k');

            orient(fig,'portrait');
            fig.Position(3)=700;
            fig.Position(4)=900;
            eval(['print -dpdf ',Work.PSAL_REPLACE_plot_dir,file_doxy(1:end-3),'_',Work.PSAL_REPLACE_CLIM,'.pdf']);
            close(fig);
        end
        clr=[0 1 0];

        % Use climatological recalculated DOXY profile if QC=4 for
        % PSAL_adjusted values
        qc_psal_num=str2num(QC_PSAL_adjusted_ctd(:,ip:ip));

        if all(qc_psal_num==4)
            f=(i-imin)/(imax-imin);
            clr=f*[255,0,0]/255;
        end
        if Work.PSAL_REPLACE_plot
            fig2=figure(100);
            title('Original ',strcat(strrep(file_first,'_','\_'),'-',strrep(file_doxy,'_','\_')));
            hold on;
            h1=plot(DOXY_ORI,-PRES,'Color',clr);
            hold on;
        end

        global_num=[global_num;i.*ones(length(DOXY_ORI),1)];
        global_fcolor=[global_fcolor;clr.*ones(length(DOXY_ORI),3)];

        if all(qc_psal_num==4)
            DOXY(:,ip)=DOXY_REC2;
            f=(i-imin)/(imax-imin);
            %clr=[252,141,14]/255;%orange
            clr=f*[0,0,255]/255;
        end

        if Work.PSAL_REPLACE_plot
            fig3=figure(101);
            title('Clim ',strcat(strrep(file_first,'_','\_'),'-',strrep(file_doxy,'_','\_')));
            hold on;
            h1=plot(DOXY_REC2,-PRES,'Color',clr);
            hold on;
        end

        global_DOXY_ORI=[global_DOXY_ORI;DOXY_ORI];
        global_DOXY_REC2=[global_DOXY_REC2;DOXY_REC2];

        global_jdays=[global_jdays,juld_float(ip)];
        global_mean_DOXY_diff_ori=[global_mean_DOXY_diff_ori,nanmean(DOXY_REC-DOXY_ORI)];
        global_mean_DOXY_diff=[global_mean_DOXY_diff,nanmean(DOXY_REC2-DOXY_ORI)];

        pmin=0;
        pmax=1000;
        idx=find(PRES>=pmin&PRES<=pmax);
        global_mean_DOXY_diff_1000=[global_mean_DOXY_diff_1000,nanmean(DOXY_REC2(idx)-DOXY_ORI(idx))];

        pmin=1000;
        pmax=2000;
        idx=find(PRES>=pmin&PRES<=pmax);
        global_mean_DOXY_diff_2000=[global_mean_DOXY_diff_2000,nanmean(DOXY_REC2(idx)-DOXY_ORI(idx))];

        pmin=2000;
        pmax=6000;
        idx=find(PRES>=pmin&PRES<=pmax);
        global_mean_DOXY_diff_6000=[global_mean_DOXY_diff_6000,nanmean(DOXY_REC2(idx)-DOXY_ORI(idx))];
    end

end
if Work.PSAL_REPLACE_plot
    figure(100);
    xlabel('DOXY (MM/Kg)');
    ylabel('Pressure (Db)');
    orient(fig2,'portrait');
    fig2.Position(3)=700;
    fig2.Position(4)=900;
    grid on;

    vmax=10*ceil(max([global_DOXY_ORI;global_DOXY_REC2])/10);
    vmin=10*floor(min([global_DOXY_ORI;global_DOXY_REC2])/10);

    set(gca,'XLIM',[vmin vmax]);
    eval(['print -dpdf ',Work.PSAL_REPLACE_plot_dir,file_doxy(1:end-3),'_all_ori','_',Work.PSAL_REPLACE_CLIM,'.pdf']);
    %close(fig2);
end

if Work.PSAL_REPLACE_plot
    figure(101);
    xlabel('DOXY (MM/Kg)');
    ylabel('Pressure (Db)');
    orient(fig2,'portrait');
    fig3.Position(3)=700;
    fig3.Position(4)=900;
    grid on;
    set(gca,'XLIM',[vmin vmax]);
    eval(['print -dpdf ',Work.PSAL_REPLACE_plot_dir,file_doxy(1:end-3),'_all_clim','_',Work.PSAL_REPLACE_CLIM,'.pdf']);
    %close(fig3);
end

if Work.PSAL_REPLACE_plot
    fig3=figure(102);

    idx=find(diff(global_num)==1);%Last values
    idx=[idx;length(global_num)];%Adding last cycle


    ibeg=1;
    for ipt=1:length(idx)
        iend=idx(ipt);
        clr=global_fcolor(idx(ipt),:)
        h1=plot(global_DOXY_ORI(ibeg:iend),global_DOXY_REC2(ibeg:iend),'+');
        h1.Color=clr;
    end


    leg1='Float versus Climatological PSAL Values';

    nan_list=find(~isnan(global_DOXY_ORI)&~isnan(global_DOXY_REC2));
    tmp=fitlm(global_DOXY_ORI(nan_list),global_DOXY_REC2(nan_list),'linear','Intercept',false);
    coeff(1)=double(tmp.Coefficients.Estimate);
    coeff(2)=0;
    hold on


    min_doxy=round(min(min(global_DOXY_ORI(nan_list),global_DOXY_REC2(nan_list))));
    min_doxy=5*floor(min_doxy/5);

    max_doxy=round(max(max(global_DOXY_ORI(nan_list),global_DOXY_REC2(nan_list))));
    max_doxy=5*ceil(max_doxy/5);


    h2=plot(min_doxy:max_doxy,coeff(1)*(min_doxy:max_doxy)+coeff(2),'r','Linewidth',2);
    h2.Color='k';
    h2.LineWidth=2;
    leg2=['Regression line : y=' num2str(coeff(1),'%4.3f') '*x + ' num2str(coeff(2),'%4.3f')];

    htl=legend([h1,h2],leg1,leg2,'Location','SouthEast');

    hold on;
    xlabel('DOXY Float (MM/Kg)');
    ylabel('DOXY CLIM (MM/Kg)');
    title('Takeshita ',strcat(strrep(file_first,'_','\_'),'-',strrep(file_doxy,'_','\_')));
    grid on;
    orient(fig3,'portrait');
    fig3.Position(3)=700;
    fig3.Position(4)=900;
    %set(gca,'XLIM',[120 280]);
    eval(['print -dpdf ',Work.PSAL_REPLACE_plot_dir,'Takeshita_',strcat(file_first(1:end-3),'-',file_doxy(1:end-3)),'_',Work.PSAL_REPLACE_CLIM,'.pdf']);
end

if Work.PSAL_REPLACE_plot
    fig4=figure(104);

    clr='m';
    h0=plot(1:length(global_jdays),global_mean_DOXY_diff_ori);
    h0.Color=clr;
    leg0='DOXY REC - ORI: 0 6000';
    hold on;


    clr='k';
    h1=plot(1:length(global_jdays),global_mean_DOXY_diff);
    h1.Color=clr;
    leg1='DOXY CLIM - ORI: 0 6000';


    clr='r';
    h2=plot(1:length(global_jdays),global_mean_DOXY_diff_1000,'+');
    h2.Color=clr;
    leg2='DOXY CLIM - ORI: 0 1000';

    clr='c';
    h3=plot(1:length(global_jdays),global_mean_DOXY_diff_2000,'+');
    h3.Color=clr;
    leg3='DOXY CLIM - ORI: 1000 2000';

    clr='b';
    h4=plot(1:length(global_jdays),global_mean_DOXY_diff_6000,'+');
    h4.Color=clr;
    leg4='DOXY CLIM - ORI: 2000 6000';


    htl=legend([h0,h1,h2,h3,h4],leg0,leg1,leg2,leg3,leg4,'Location','SouthEast');

    xlabel('Cycle Number');
    ylabel('CLIM-FLOAT (MM/Kg)');
    title('DOXY Diff ',strcat(strrep(file_first,'_','\_'),'-',strrep(file_doxy,'_','\_')));
    grid on;
    orient(fig4,'portrait');
    fig4.Position(3)=700;
    fig4.Position(4)=900;
    %set(gca,'XLIM',[120 280]);
    eval(['print -dpdf ',Work.PSAL_REPLACE_plot_dir,'DOXY_DIFF_',strcat(file_first(1:end-3),'-',file_doxy(1:end-3)),'_',Work.PSAL_REPLACE_CLIM,'.pdf']);
end

line=['save global_',cwmo,'_',Work.PSAL_REPLACE_CLIM,'.mat global_mean_DOXY_diff*'];
eval(line);



