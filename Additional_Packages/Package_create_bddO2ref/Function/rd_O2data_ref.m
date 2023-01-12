function [lon,lat,sta_ctd,juld,pres,temp,tpot,psal,sig0,doxy,o2ok] = rd_O2data_ref( campyear,datatype )
%
% Fonction qui lit les données oxygène calibrées (campagne et bouteille)
% pour calibrer les données flotteurs
% Input:
%        campyear = année des données de référence
%                 = ov10th/ov11me/ov11di/gh10md/ov12sg
%        datatype = type de données lues
%                 = 'ctd' ou 'btl'
%

switch campyear
    
    case 'strass'
        o2ok=1;
        % chargement des donn�es r�alis�es au d�ploiement
        strasse=load('/home4/homedir4/perso/vthierry/matlab/ARGO_OVIDE/ARVOR350/CTDeduitup.017');
        strasse=flipud(strasse);
        pres=strasse(:,2);
        temp=strasse(:,3);
        psal=strasse(:,4);
        doxy=strasse(:,5);
        tpot=sw_ptmp(psal,temp,pres,0);
        sig0=sw_pden(psal,temp,pres,0)-1000;
        lon=-35.700;
        lat=25.98; 
        sta_ctd=1;
        juld='NaN';
        
    case 'ov12sg'
        
        if strcmp(datatype,'ctd') == 1
            % Lecture des donnees campagne CTD
            o2ok=1;
            %inpath='/home1/armen/perso/vthierry/PROJETS/OVIDE/CATARINA/DATA/';
            %fname=['cat_PRES.nc'];
            inpath='/home1/penfret/HYDROCEAN/MLT_NC/LPO/OVIDE/';
            fname=['cat12_PRES.nc'];
            %fname=['cat12_DEPH.nc'];
            file_name=[inpath fname];
            [lat,lon,juld,platform,temp,tpot,psal,doxy,pres,sta_ctd]=...
            read_nc_base_hydro(file_name,'TEMP','TPOT','PSAL','OXYK','PRES','STATION_NUMBER');
            temp(temp==-9999)=NaN;
            tpot(tpot==-9999)=NaN;
            psal(psal==-9999)=NaN;
            doxy(doxy==-9999)=NaN;
            pres(pres==-9999)=NaN;
           % [null,sig0]=swstat90(psal,tpot,0);
            sig0=sw_pden(psal,temp,pres,0)-1000;
            sta_ctd=sta_ctd';
            
        elseif strcmp(datatype,'btl') == 1
            % Lecture des donnees campagne bouteilles
            o2ok=0;
            pres=NaN;temp=NaN;tpot=NaN;psal=NaN;sig0=NaN;oxyl=NaN;doxy=NaN;
        end
    
    case 'ov10th'
        
        if strcmp(datatype,'ctd') == 1
            % Lecture des donnees campagne CTD
            o2ok=1;
            %inpath='/home1/kereon/OVIDE2010/ctd/';
            inpath='/home1/penfret/HYDROCEAN/MLT_NC/LPO/OVIDE/';
            fname=['ov10_PRES.nc'];
            file_name=[inpath fname];
            [lat,lon,juld,platform,temp,tpot,psal,oxyl,sig0,pres,sta_ctd]=...
            read_nc_base_hydro(file_name,'TEMP','TPOT','PSAL','OXYL','SIG0','PRES','STATION_NUMBER');
            temp(temp==-9999)=NaN;
            sig0(sig0==-9999)=NaN;
            tpot(tpot==-9999)=NaN;
            psal(psal==-9999)=NaN;
            oxyl(oxyl==-9999)=NaN;
            doxy=convert_oxygen(oxyl,'mL/L','mumol/kg',sig0);
            sta_ctd=sta_ctd';
        
        elseif strcmp(datatype,'btl') == 1
            % Lecture des donnees campagne bouteilles
            o2ok=1;
            fname='DATA/ovid10.sal';
            fid=fopen(fname,'r');
            ent=fgetl(fid);
            ind=1;
            while ind < 2448
                li=fgetl(fid);
                virg=findstr(li,',');
                stat_btl(ind)=str2num(li(1:virg(1)-1));
                pres(ind)=str2num(li(virg(2)+1:virg(3)-1));
                temp(ind)=str2num(li(virg(3)+1:virg(4)-1));
                tmp2_btl(ind)=str2num(li(virg(6)+1:virg(7)-1));
                clear tmp
                tmp=str2num(li(virg(11)+1:virg(12)-1));
                if isempty(tmp) == 1
                    psal(ind)=NaN;
                else
                    psal(ind)=str2num(li(virg(11)+1:virg(12)-1));
                end

                clear tmp
                tmp=str2num(li(virg(12)+1:virg(13)-1));
                if isempty(tmp) == 1
                    oxyl(ind)=NaN;
                else
                    oxyl(ind)=str2num(li(virg(12)+1:virg(13)-1));
                end
                ind=ind+1;
            end
            psal=psal';temp=temp';pres=pres';oxyl=oxyl';
            %tpot=pottemp(psal,temp,pres,0);
            %[null,sig0]=swstat90(psal,tpot,0);
            tpot=sw_ptmp(psal,temp,pres,0);
            sig0=sw_pden(psal,temp,pres,0)-1000;
            doxy=convert_oxygen(oxyl,'mL/L','mumol/kg',sig0);
            fclose(fid)  
        end
        
       case 'ov08th'
        
        if strcmp(datatype,'ctd') == 1
            % Lecture des donnees campagne CTD
            o2ok=1;
            %inpath='/home1/kereon/OVIDE2008/ctd/';
            inpath='/home1/penfret/HYDROCEAN/MLT_NC/LPO/OVIDE/';
            fname=['ovid08_prs.nc'];
            file_name=[inpath fname];
            [lat,lon,juld,platform,temp,tpot,psal,oxyl,sig0,pres,sta_ctd]=...
            read_nc_base_hydro(file_name,'TEMP','TPOT','PSAL','OXYL','SIG0','PRES','STATION_NUMBER');
            temp(temp==-9999)=NaN;
            sig0(sig0==-9999)=NaN;
            tpot(tpot==-9999)=NaN;
            psal(psal==-9999)=NaN;
            oxyl(oxyl==-9999)=NaN;
            doxy=convert_oxygen(oxyl,'mL/L','mumol/kg',sig0);
            
            sta_ctd=sta_ctd';
        elseif strcmp(datatype,'btl') == 1
             % Lecture des donnees campagne bouteilles
            o2ok=0;
            pres=NaN;temp=NaN;tpot=NaN;psal=NaN;sig0=NaN;oxyl=NaN;doxy=NaN;
        end
        
    case 'ov11me'

        if strcmp(datatype,'ctd') == 1
            % Lecture des donnees campagne CTD
            o2ok=1;
            
            tabnum=[1901209;
                1901210;
                1901211;
                1901212;
                1901213;
                1901214;
                1901215;
                1901217;
                1901218];

            
            tabo2file=['m852_003.ctd';
                       'm852_027.ctd';
                       'm851_014.dat';
                       'm851_018.dat';
                       'm851_026.dat';
                       'm851_079.dat';
                       'm851_082.dat';
                       'm851_095.dat';
                       'm851_106.dat'];
            [nfile,null]=size(tabo2file);
            inpathlpo='/home4/homedir4/perso/vthierry/';
            inpathmac='/Users/vthierry/';
            if exist(inpathlpo) == 7
                inpath=[inpathlpo 'PROJETS/OXYGEN/ACHATS/2010/OXYBTL4CALIB/'];
            elseif exist(inpathmac) == 7
                inpath=['/Users/vthierry/PROJETS/OXYGEN/ACHATS/2010/OXYBTL4CALIB/'];
            else
                disp('Pas de donnees de reference');
                stop
            end
            lat=NaN*ones(nfile,1);
            lon=NaN*ones(nfile,1);
            juld=NaN*ones(nfile,1);
            pres=NaN*ones(nfile,2100);
            %deph=NaN*ones(nfile,2100);
            temp=NaN*ones(nfile,2100);
            psal=NaN*ones(nfile,2100);
            doxy=NaN*ones(nfile,2100);
            sig0=NaN*ones(nfile,2100);
            sta_ctd=[1:nfile]';

            for ifile=1:nfile
                fname=tabo2file(ifile,:);
                filename=[inpath num2str(tabnum(ifile)) '/' fname]
                switch fname(1:4)
                    case 'm852'
                        [latf,lonf,juldf,presf,dephf,tempf,psalf,doxyovf,sig0f] = rd_m852(filename);
                        lat(ifile)=latf;
                        lon(ifile)=lonf;
                        juld(ifile)=juldf;
                        pres(ifile,1:min(length(presf),2100))=presf(1:min(length(presf),2100));
                        %deph(ifile,1:min(length(presf),2100))=dephf(1:min(length(presf),2100));
                        temp(ifile,1:min(length(presf),2100))=tempf(1:min(length(presf),2100));
                        psal(ifile,1:min(length(presf),2100))=psalf(1:min(length(presf),2100));
                        doxy(ifile,1:min(length(presf),2100))=doxyovf(1:min(length(presf),2100));
                        sig0(ifile,1:min(length(presf),2100))=sig0f(1:min(length(presf),2100));
                        tpot=pottemp(psal,temp,pres,0);
          

                    case 'm851'
                        [latf,lonf,juldf,presf,tempf,tpotf,psalf,doxyovf,sig0f] = rd_m851(filename);
                        lat(ifile)=latf;
                        lon(ifile)=lonf;
                        juld(ifile)=juldf;
                        pres(ifile,1:min(length(presf),2100))=presf(1:min(length(presf),2100));
                        temp(ifile,1:min(length(presf),2100))=tempf(1:min(length(presf),2100));
                        psal(ifile,1:min(length(presf),2100))=psalf(1:min(length(presf),2100));
                        doxy(ifile,1:min(length(presf),2100))=doxyovf(1:min(length(presf),2100));
                        sig0(ifile,1:min(length(presf),2100))=sig0f(1:min(length(presf),2100));
                        tpot=pottemp(psal,temp,pres,0);
                    case 'disc'
                end
            end

        elseif strcmp(datatype,'btl') == 1
            % Lecture des donnees campagne bouteilles
            o2ok=0;
            pres=NaN;temp=NaN;tpot=NaN;psal=NaN;sig0=NaN;oxyl=NaN;doxy=NaN;
        end
        
    case 'ov11di'
        if strcmp(datatype,'ctd') == 1
            o2ok=1;
            % Lecture des donnees campagne CTD
            inpath='/home4/homedir4/perso/vthierry/PROJETS/OXYGEN/ACHATS/2010/DISCOVERY/DATA/';
            pres=NaN*ones(29,6000);
            presqc=NaN*ones(29,6000);
            temp=NaN*ones(29,6000);
            tempqc=NaN*ones(29,6000);
            psal=NaN*ones(29,6000);
            psalqc=NaN*ones(29,6000);
            doxy=NaN*ones(29,6000);
            doxyqc=NaN*ones(29,6000);
            sta_ctd=NaN*ones(29,1);
            lat=NaN*ones(29,1);
            lon=NaN*ones(29,1);
            juld=NaN*ones(29,1);
            for iprf=1:29
                fname=sprintf('a16n2011_000%2.2d_00001_ct1.csv',iprf)
                fid=fopen([inpath fname],'r');
                sta_ctd(iprf)=iprf;
                for il=1:8
                    vline=fgetl(fid);
                end
                vline=fgetl(fid);
                ivirg=findstr(vline,'=');
                vdate=deblank(vline(ivirg+1:end))
                juld=datenum(str2double(vdate(1:5)),str2double(vdate(6:7)),str2double(vdate(8:9)));
                datevec(juld)
                vline=fgetl(fid);
                ivirg=findstr(vline,'=');
                vtime=deblank(vline(ivirg+1:end));
                vline=fgetl(fid);
                ivirg=findstr(vline,'=');
                lat(iprf)=str2num(vline(ivirg+1:end));
                vline=fgetl(fid);
                ivirg=findstr(vline,'=');
                lon(iprf)=str2num(vline(ivirg+1:end));
                
                for il=1:3
                    vline=fgetl(fid);
                end
                vline=fgetl(fid);
                ik=1;
                while strcmp(vline,'END_DATA') ~= 1
                    valline=str2num(vline);
                    pres(iprf,ik)=valline(1);
                    presqc(iprf,ik)=valline(2);
                    temp(iprf,ik)=valline(3);
                    tempqc(iprf,ik)=valline(4);
                    psal(iprf,ik)=valline(5);
                    psalqc(iprf,ik)=valline(6);
                    doxy(iprf,ik)=valline(7);
                    doxyqc(iprf,ik)=valline(8);
                    vline=fgetl(fid);
                    ik=ik+1;
                end
            end
            pres(find(pres==-999.0 | presqc == 9))=NaN;
            temp(find(temp==-999.0 | presqc == 9))=NaN;
            psal(find(psal==-999.0 | presqc == 9))=NaN;
            doxy(find(doxy==-999.0 | presqc == 9))=NaN;
            tpot = sw_ptmp(psal,temp,pres,0);
            sig0 = sw_pden(psal,temp,pres,0);

           
                 
        elseif strcmp(datatype,'btl') == 1
            % Lecture des donnees campagne bouteilles
            o2ok=0;
            pres=NaN;temp=NaN;tpot=NaN;psal=NaN;sig0=NaN;oxyl=NaN;doxy=NaN;
        end
end


