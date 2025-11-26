function out=pH2O(TEMP,S)
load("D.mat")
out=1013.25*exp(D0+D1*(100./(TEMP+273.15))+D2*log((TEMP+273.15)./100)+D3*S);
return
end
