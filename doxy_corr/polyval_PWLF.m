function YI=polyval_PWLF(coef,XI,x)
% HISTORY
%   $created: 26.04.2024 $author: Thierry Reynaud
%   $Revision: version $Date: $author:
%            v1 26.04.2024     Thierry Reynaud LOPS

% This function is a multi linear segment version of the matlab function polyval
% It returns interpolated values upon linear segment.
% 
% INPUT
%           Linear coefficents for each segment:
%           y1=a1*x + b1
%           y2=a2*x + b2
%
%           coef=  a1     a2
%                  b1     b2
%           XI linear segment X limits
%           x values for interpolation
% OUTPUT
%           YI interpolated values
%           

YI=nan*ones(size(x));
for i=1:size(coef,2)
    ibeg=i;
    iend=i+1;

    xmin=XI(ibeg);
    xmax=XI(iend);
    if i<size(coef,2)
     idx=find(x>=xmin & x<xmax);
    else
     idx=find(x>=xmin & x<=xmax);
    end

    cf=coef(:,i);
    y=polyval(cf,x(idx));
    YI(idx)=y;
end

return
end