% DOXY_PSAL_clim_replace_apply replaces float PSAL values by climatological ones.
%
% SYNTAX
% [[DATA]=DOXY_PSAL_clim_replace_apply(ptype,DATA,Work,timeclim, depthclim, latclim, lonclim,psalclim,cycle_beg,cycle_end)
%
% DESCRIPTION
% DOXY_PSAL_clim_replace replaces float PSAL values by climatological ones
% from cycle_beg to cycle_end
%
% INPUT
%     timeclim            climatology: time (days) relative to the current year: 15,45, .....
%     depthclim           climatology: depth (m)
%     latclim:            climatology: latitude (degrees) -90 to 90
%     lonclim:            climatology: longitude (degrees) 0 to 360
%     psalclim:           climatology field: 4D variables [time,vertical, latitude,longitude]]
%     ptype               profil type 'primary','nearsurf','seconday', ...
%     cycle_beg           replacement occures from cycle_beg
%     cycle_end           and ends at cycle_end
%     DATA (structure)    argo float "multiprofile" structure gathering the
%                         three kind of vertical sampling scheme. Subfields
%                         are argo (the data) and Dim (the dimension).
%                         Example :
%                             DATA =
%                                   primary: [1x1 struct]
%                                  nearsurf: [1x1 struct]
%                                 secondary: [1x1 struct]
%                                     other: [1x1 struct]
%
%                             DATA.primary =
%                                 argo: [1x1 struct]
%                                  Dim: [1x1 struct]
%
%                             DATA.primary.argo =
%                                   dimorder: 'C'
%                                  data_type: [1x1 struct]
%                                          ...
%                                       pres: [1x1 struct]
%                                       doxy: [1x1 struct]
%                                    doxy_qc: [1x1 struct]
%                              doxy_adjusted: [1x1 struct]
%                           doxy_adjusted_qc: [1x1 struct]
%                        doxy_adjusted_error: [1x1 struct]
%                                           ...
%    Work (structure)     Doxy correction working structure, issued and
%                         computed from argo float data.
%                         Example:
%                             Work =
%                                   readme: [1x1 struct]
%                                     unit: 'mumol/kg'
%                                      wmo: 1901205
%                                 DOXY_RAW: [85x120 single]
%                                  CAPTEUR: {'Optode'}
%                                  DOXY_QC: [85x120 char]
% OUTPUT
% DATA
% out_std                 structures containing float climatological STD
%                         interpolated values for PSAL and TEMP.
% HISTORY
%   $created: 14/12/2023 $author: Thierry Reynaud
%   $Revision: version $Date: $author:
%       v1.0 12/12/2023   Thierry Reynaud, LOPS
%                         initial
%
%       v2.0 18/01/2024   Thierry Reynaud, LOPS
%                         adding figures details and std values

function [DATA,out_std]=DOXY_PSAL_clim_replace_apply(ptype,DATA,Work,CLIM,cycle_beg,cycle_end)


timeclim=CLIM.timeclim;
depthclim=CLIM.depthclim;
latclim=CLIM.latclim;
lonclim=CLIM.lonclim;
psalclim=CLIM.psalclim;
tempclim=CLIM.tempclim;

if isfield(CLIM,'std_psalclim')
    std_psalclim=CLIM.std_psalclim;
end

if isfield(CLIM,'std_tempclim')
    std_tempclim=CLIM.std_tempclim;
end
clear CLIM

% Float P,T,S values
eval(['pres_float=DATA.',ptype,'.argo.pres.data;']);
eval(['psal_float=DATA.',ptype,'.argo.psal.data;']);
eval(['psal_float=DATA.',ptype,'.argo.psal.data;']);
eval(['psal_adjusted_float=DATA.',ptype,'.argo.psal_adjusted.data;']);
eval(['temp_float=DATA.',ptype,'.argo.temp.data;']);
eval(['cycle_float=DATA.',ptype,'.argo.cycle_number.data;']);
eval(['qc_psal_adjusted_float=DATA.',ptype,'.argo.psal_adjusted_qc.data;']);
eval(['qc_psal_float=DATA.',ptype,'.argo.psal_qc.data;']);


pres_float=double(pres_float);
psal_float=double(psal_float);
temp_float=double(temp_float);

if exist('std_psalclim','var')
    std_clim_psal_float=nan*zeros(size(psal_float));
end

if exist('std_tempclim','var')
    std_clim_temp_float=nan*zeros(size(temp_float));
end

% Float Positions
eval(['lon_float = DATA.',ptype,'.argo.longitude.data;']);
eval(['lat_float = DATA.',ptype,'.argo.latitude.data;']);

% Float Julian Days
eval(['tmp=DATA.',ptype,'.argo.reference_date_time.data;']);
jref=datenum(str2num(tmp(1:4)),str2num(tmp(5:6)),str2num(tmp(7:8)),str2num(tmp(9:10)),str2num(tmp(11:12)),str2num(tmp(13:14)));
eval(['juld_float = DATA.',ptype,'.argo.juld.data + jref;']);

tabyear = datestr(juld_float,10);
running_day_float= juld_float - datenum(str2num(tabyear),1,1)+1;

idx=find(lon_float<0);
lon_float_tmp = lon_float;
lon_float_tmp(idx) = lon_float(idx) + 360;

for ilat=1:length(lat_float)
    if cycle_float(ilat)>=cycle_beg && cycle_float(ilat)<= cycle_end

        display(num2str(ilat));
        depth_float = sw_dpth(pres_float(ilat,:),lat_float(ilat));

        if sum(~isnan(depth_float))
            % Interpolation clim.
            % interpolation on Z levels:

            % INTERPOLATION: CLIMATOLGY PSAL

            clear tmp2 deep kdeep deep2 tmp
            tmp = interpn(timeclim, depthclim, latclim, lonclim,psalclim,running_day_float(ilat),depth_float,lat_float(ilat),lon_float_tmp(ilat));
            if sum(~isnan(tmp))<2
                tmp(:)=NaN;
            end

            % Find the deepest profile around
            [deep,kdeep]=get_clim_index(psalclim,timeclim,lonclim,latclim,running_day_float(ilat),lat_float(ilat),lon_float_tmp(ilat));
            deep2=interp1(depthclim,deep,depth_float,'linear');

            %Fill NaN values
            if all(~isnan(tmp))
                %tmp without Nan values
                tmp2=tmp;
            else
                % tmp with NaN values
                ilast=find(isnan(tmp)==1,1,'first')-1;% Last real value

                % fill NaN values within the profile
                if ilast>0
                    dz=diff(tmp);
                    idz=find(dz==0);
                    tmp(idz)=NaN;
                    tmp2=interp1(depth_float(~isnan(tmp)),tmp(~isnan(tmp)),depth_float,'linear');
                else
                    tmp2=tmp;
                end

                % Fill deepest values
                tmp2(ilast+1:end)=deep2(ilast+1:end);
            end

            if Work.PSAL_REPLACE_plot
                fig=figure(99);
                %Plot PSAL
                subplot(1,3,1);
                h1=plot(tmp2,-pres_float(ilat,:),'+r');
                hold on;
                h2=plot(psal_float(ilat,:),-pres_float(ilat,:),'+g');
                xlabel('Salinity (PSU)');
                ylabel('Pressure (Db)');
                hleg=legend([h1,h2],'Clim','Float','Location','SouthEast');
                %psal_float (:,ilat)= squeeze(tmp);
            end

            psal_float(ilat,:)=tmp2;
            psal_adjusted_float(ilat,:)=tmp2;
            idx=find(~isnan(tmp2));
            qc_psal_adjusted_float(ilat,idx)=1;%Switch QCs
            qc_psal_float(ilat,idx)=1;%Switch QCs

            % Additionnal interpolations
            % INTERPOLATION: CLIMATOLGY PSAL STD
            if exist('std_psalclim','var')
                clear tmp2 deep kdeep deep2 tmp
                tmp = interpn(depthclim, latclim, lonclim,std_psalclim,depth_float,lat_float(ilat),lon_float_tmp(ilat));
                if sum(~isnan(tmp))<2
                    tmp(:)=NaN;
                end
                % Find the deepest profile around
                [deep,kdeep]=get_clim_index(std_psalclim,timeclim,lonclim,latclim,running_day_float(ilat),lat_float(ilat),lon_float_tmp(ilat));
                deep2=interp1(depthclim,deep,depth_float,'linear');

                %Fill NaN values
                if all(~isnan(tmp))
                    %tmp without Nan values
                    tmp2=tmp;
                else
                    % tmp with NaN values
                    ilast=find(isnan(tmp)==1,1,'first')-1;% Last real value

                    % fill NaN values within the profile
                    if ilast>0
                        dz=diff(tmp);
                        idz=find(dz==0);
                        tmp(idz)=NaN;
                        tmp2=interp1(depth_float(~isnan(tmp)),tmp(~isnan(tmp)),depth_float,'linear');
                    else
                        tmp2=tmp;
                    end

                    % Fill deepest values
                    tmp2(ilast+1:end)=deep2(ilast+1:end);
                end

                std_clim_psal_float(ilat,:)=tmp2;
            end

            % INTERPOLATION: CLIMATOLGY TEMP
            clear tmp2 deep kdeep deep2 tmp
            tmp = interpn(timeclim, depthclim, latclim, lonclim,tempclim,running_day_float(ilat),depth_float,lat_float(ilat),lon_float_tmp(ilat));
            if sum(~isnan(tmp))<2
                tmp(:)=NaN;
            end
            % Find the deepest profile around
            [deep,kdeep]=get_clim_index(tempclim,timeclim,lonclim,latclim,running_day_float(ilat),lat_float(ilat),lon_float_tmp(ilat));
            deep2=interp1(depthclim,deep,depth_float,'linear');

            %Fill NaN values
            if all(~isnan(tmp))
                %tmp without Nan values
                tmp2=tmp;
            else
                % tmp with NaN values
                ilast=find(isnan(tmp)==1,1,'first')-1;% Last real value

                % fill NaN values within the profile
                if ilast>0
                    dz=diff(tmp);
                    idz=find(dz==0);
                    tmp(idz)=NaN;
                    tmp2=interp1(depth_float(~isnan(tmp)),tmp(~isnan(tmp)),depth_float,'linear');
                else
                    tmp2=tmp;
                end

                % Fill deepest values
                tmp2(ilast+1:end)=deep2(ilast+1:end);
            end


            if Work.PSAL_REPLACE_plot
                fig=figure(99);
                %Plot TEMP
                subplot(1,3,2);
                h1=plot(tmp2,-pres_float(ilat,:),'+r');
                hold on;
                h2=plot(temp_float(ilat,:),-pres_float(ilat,:),'+g');
                xlabel('Temperature (PSU)');
                ylabel('Pressure (Db)');
                hleg=legend([h1,h2],'Clim','Float','Location','SouthEast');
                %psal_float (:,ilat)= squeeze(tmp);
                if Work.PSAL_REPLACE_plot_close
                    close(fig);
                end
            end

            % INTERPOLATION: CLIMATOLGY PSAL STD
            if exist('std_tempclim','var')

                clear tmp2 deep kdeep deep2 tmp
                tmp = interpn(depthclim, latclim, lonclim,std_tempclim,depth_float,lat_float(ilat),lon_float_tmp(ilat));
                if sum(~isnan(tmp))<2
                    tmp(:)=NaN;
                end

                % Find the deepest profile around
                [deep,kdeep]=get_clim_index(std_tempclim,timeclim,lonclim,latclim,running_day_float(ilat),lat_float(ilat),lon_float_tmp(ilat));
                deep2=interp1(depthclim,deep,depth_float,'linear');

                %Fill NaN values
                if all(~isnan(tmp))
                    %tmp without Nan values
                    tmp2=tmp;
                else
                    % tmp with NaN values
                    ilast=find(isnan(tmp)==1,1,'first')-1;% Last real value

                    % fill NaN values within the profile
                    if ilast>0
                        dz=diff(tmp);
                        idz=find(dz==0);
                        tmp(idz)=NaN;
                        tmp2=interp1(depth_float(~isnan(tmp)),tmp(~isnan(tmp)),depth_float,'linear');
                    else
                        tmp2=tmp;
                    end

                    % Fill deepest values
                    tmp2(ilast+1:end)=deep2(ilast+1:end);
                end

                std_clim_temp_float(ilat,:)=tmp2;
            end


        end

    end
end

out_std=[];
if exist('std_psalclim','var')
    out_std.std_clim_psal_float=std_clim_psal_float;
end

if exist('std_tempclim','var')
    out_std.std_clim_temp_float=std_clim_temp_float;
end

eval(['DATA.',ptype,'.argo.psal.data=psal_float;']);
eval(['DATA.',ptype,'.argo.psal_adjusted.data=psal_adjusted_float;']);
eval(['DATA.',ptype,'.argo.psal_adjusted_qc.data=qc_psal_adjusted_float;']);
eval(['DATA.',ptype,'.argo.psal_qc.data=qc_psal_float;']);


return
end
