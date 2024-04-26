function coef=polyfit_PWLF(x,y,XI,n)
% HISTORY
%   $created: 26.04.2024 $author: Thierry Reynaud
%   $Revision: version $Date: $author:
%            v1 26.04.2024     Thierry Reynaud LOPS

% This function returns the linear coefficents for each segment:
% y1=a1*x + b1
% y2=a2*x + b2
%
% INPUT
%           x,y data
%           XI linear segment X limits
%           n number of linear segments
% OUTPUT
%           coef=  a1     a2
%                  b1     b2

% function lsq_lut_piecewise.m dowloaded from the matlab website:
% Guido Albertin (2024). Piecewise linear least square fit (https://www.mathworks.com/matlabcentral/fileexchange/40913-piecewise-linear-least-square-fit), MATLAB Central File Exchange. 

% obtain vector of 1-D look-up table "y" points
if size(x,1)==1
    YI = lsq_lut_piecewise( x', y', XI );
else
    YI = lsq_lut_piecewise( x, y, XI );
end

coef=NaN*ones(2,n);% [ 2, number of line segments]
for i=1:n
    ibeg=i;
    iend=i+1;
    p=polyfit(XI(ibeg:iend),YI(ibeg:iend),1);
    coef(:,i)=p;
end

return
end
