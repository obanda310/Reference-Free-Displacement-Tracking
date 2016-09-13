function [pF] = pillarFit2(x,S,dS,dP,P,Pg,Mi,Mf)
pF = Mi./(1+exp(S.*x-dS))+(Mf.*sin(((2*pi)./(P+Pg.*x)).*x+dP)).*(1./(1+exp((S).*x-dS)));