function [pF] = pillarFit(x,S,iS,iP,P,Pg,Mo,Mf)
% Independent variable 'x' is z-frame
% Dependant variable is Intensity
% All other variables are coefficients

% Some reasonable ICs for the coefficients:

% S=0.1985; Steepness of Decay - Should Relate to Material Density Profile of
% Hydrogel in Z-dimension, which I approximate as a sigmoidal function.

% iS=7.792; Decay Shift - Should relate to pattern position relative to surface

% iP=1.759; Period Shift - Should relate to pattern position relative to surface

% P=10.87; Base Period - Should relate to center-center distance of
% ellipsoids without swelling

% Pg=0.03287; Period Elongation - Should relate to Hydrogel Swelling

% Mi=28270; Maximum feature intensity with oscillations subtracted

% Mo=5761; Magnitude of oscillations - Should relate to PSF shape from lithography

% IMPORTANT: Mi and Mf will scale based on how you process the intensity
% information. Here they were obtained using TrackMate, which outputs the
% average intensity within features it detects.


Decay = 1./(1+exp(S.*x-iS)); %General Logistic Decay Formula - Should 
% describe material density profile in z.

Period = (2*pi)./(P+Pg.*x); %Should describe Z-Spacing information with a
% swelling factor: 'Pg*x'. There is probably a more appropriate way to
% describe swelling, but this seems to work for now.

Oscillation = sin(Period.*x+iP); %General Oscillation - Should describe 
% PSF locations. Local maxima should line up with ellipsoid centers

pF = (Mo + Mf.*Oscillation).*Decay; %Scaled Product of Osc. and Decay
end

% Here is how the function is implemented in my code:
%
% xVals = linspace(1,totalNumFrames,totalNumFrames);
% fitWeights = zeros(totalNumFrames,1);
% fitWeights(:,1) = xVals(1,:) + totalNumFrames;
% s = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[4000,25000,0,0,0,0,0],...
%                'Upper',[6000,35000,15,.1,1,3.14,15],...
%                'StartPoint',[5500,29000,11.5,.015,0.1731,2.2,6.712]);
% f = fittype('pillarFit2(x,S,iS,iP,P,Pg,Mo,Mf)','options',s);
% [fitCoef,fitQual] = fit(xVals',currentAverages,f,'weight',fitWeights)
% fitCoef2 = coeffvalues(fitCoef);
% yVals = pillarFit(xVals,fitCoef2(1,5),fitCoef2(1,7),fitCoef2(1,6),fitCoef2(1,3),fitCoef2(1,4),fitCoef2(1,2),fitCoef2(1,1));
% plot(xVals,yVals,'.','MarkerSize',20,'Color','bl')