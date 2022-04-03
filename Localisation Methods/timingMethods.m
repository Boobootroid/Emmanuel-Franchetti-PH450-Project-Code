clear all; close all;  
warning off

% Set Parameters for Least Squares fit. Worse parameters will increase
% computation time
imgSigma = 5;
noiseVariance = 100;
photonNumber = 1e3;
QE = 0.9;
    photonNumber = photonNumber*QE;
location = '1000 Image Data\';
imgName = strcat(location,'IMG50_2000.tif');
image = imread(imgName);
[X,Y] = size(image);

numSamples = 100; %Number of times timeit is ran for each method

%% Establish guess parameters for LSq
% N
lowerPhotonNumber=photonNumber-200; startPhotonNumber=photonNumber; upperPhotonNumber=photonNumber+200;
% xo and yo
locDispersion = 2; 
    startX = X/2; lowerX = startX-locDispersion; upperX = startX + locDispersion;
    startY = Y/2; lowerY = startY-locDispersion; upperY = startY + locDispersion;
% Sigma
sigDisp = 0.5;
    lowerSigma = (imgSigma-sigDisp); startSigma = imgSigma; upperSigma = (imgSigma+sigDisp);
%Establish final guess parameter for fit, A
minIm = min(min(image));
    lowerOffset=minIm; startOffset=minIm+2*sqrt(noiseVariance); upperOffset=minIm+4*sqrt(noiseVariance);

guess.Lower = [lowerOffset,lowerPhotonNumber,lowerX,lowerY,lowerSigma]; % lower bounds of models parameters
guess.StartPoint = [startOffset,startPhotonNumber,startX,startY,startSigma]; % starting parameter guess
guess.Upper = [upperOffset,upperPhotonNumber,upperX,upperY,upperSigma]; % upper bounds

%% Least Squares
l = @() lsq2DGaussian(image,guess,0);
%% Weighted Centroiding
% Threshold image for centroiding
threshold = minIm+4*sqrt(noiseVariance);
BW = image > threshold;
% Perform weighted centroid
w = @()  regionprops(BW,image,{'Centroid','WeightedCentroid'});

%% Timing

for i = 1:numSamples
    LSTime(i)=timeit(l);
     WCTime(i)=timeit(w);
end


LSMean = mean(LSTime)
LSStd = std(LSTime)

WCMean = mean(WCTime)
WCStd = std(WCTime)

