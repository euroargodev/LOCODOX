function [deep,kdeep]=get_clim_index(clim,timeclim,lonclim,latclim,running_day_float,lat_float,lon_float_tmp)

iii(2)=find(lonclim>=lon_float_tmp,1,'First');
iii(1)=iii(2)-1;

jjj(2)=find(latclim>=lat_float,1,'First');
jjj(1)=jjj(2)-1;

ideep=0;
jdeep=0;
kdeep=0;

ttt(2)=find(timeclim>=running_day_float,1,'First');
ttt(1)=ttt(2)-1;



if length(size(clim))==4


    for j=1:length(jjj)
        for i=1:length(iii)
            ii=iii(i);
            jj=jjj(j);
            ilast=find(isnan(squeeze(clim(ttt(1),:,jj,ii)))==1,1,'first');
            if ilast>kdeep
                kdeep=ilast;% First isnan level
                ideep=ii;
                jdeep=jj;
            end

        end
    end

    s2=squeeze(clim(ttt(2),:,jdeep,ideep));
    s1=squeeze(clim(ttt(1),:,jdeep,ideep));

    deltat=timeclim(ttt(2))-timeclim(ttt(1));
    f2=(running_day_float-timeclim(ttt(1)))/deltat;
    f1=(timeclim(ttt(2))-running_day_float)/deltat;

    deep=f1.*s1+f2.*s2;

else

    for j=1:length(jjj)
        for i=1:length(iii)
            ii=iii(i);
            jj=jjj(j);
            ilast=find(isnan(squeeze(clim(:,jj,ii)))==1,1,'first');
            if ilast>kdeep
                kdeep=ilast;% First isnan level
                ideep=ii;
                jdeep=jj;
            end

        end
    end

    f1=1;
    s1=squeeze(clim(:,jdeep,ideep));

    deep=f1.*s1;


end

% Fill up the profile bottom part with the last value.
if kdeep>1
    deep(kdeep:end)=deep(kdeep-1);
end

%if ~isempty(ilast)&&ilast>1
%    tmp2(ilast+1:end)=tmp2(ilast);
%end
return
end