%% Simulates PSF one photon at a time, then adds additive gaussian white noise.
function [imgDat, realCoordinates] = DripModelSpotSimulation(spotParams,imageParams)

    %%% PASS IN SPOT PARAMETERS
    numberOfSpots = spotParams.numberOfSpots;
    spotSpreadRange =spotParams.spotSpreadRange;% Sigma value in PSF
    photonNumber = spotParams.avgSpotPhotonNumber;% Photons/Spot/shutterTime

    %%% PASS IN IMAGE PARAMETERS
    quantumEfficiency = imageParams.quantumEfficiency;
    imageResolution = imageParams.imageResolution;
    pixelScale = imageParams.pixelScale;
    imageSize = imageParams.imageSize;
    noiseMean = imageParams.noiseMean;
    noiseVariance = imageParams.noiseVariance;


%%% RANDOMLY DISTRIBUTE SPOTS ACROSS REAL SPACE,
%%% INTEGRATE OVER SPOTS FOR EACH PIXEL FOR PDF
photonCount = zeros(imageResolution(1),imageResolution(2));
photoelectronCount = zeros(imageResolution(1),imageResolution(2));
for n= 1:numberOfSpots
    
    % Generating random realspace corrdinates within 1 pixel of the centre
    % of the image
    upperLimX = imageSize(1)*0.5+pixelScale; lowerLimX = imageSize(1)*0.5-pixelScale;  
    upperLimY = imageSize(2)*0.5+pixelScale; lowerLimY = imageSize(2)*0.5-pixelScale;
    
    realCoordinates=[(upperLimX-lowerLimX)*rand+lowerLimX,(upperLimY-lowerLimY)*rand+lowerLimY];
    
    %Generate random spread within set boundaries
    spotSpread = (spotSpreadRange(2)-spotSpreadRange(1))*rand+spotSpreadRange(1);
    
    %%% GENERATE GAUSSIAN SPREAD FOR 1 SPOT
    [PSF,~] = (gaussianPSF(realCoordinates(1),realCoordinates(2),spotSpread));
    PSFunc = matlabFunction(PSF);       
    
    %%% 2D INTEGRATION TO FIND INTENSITY PER PIXEL
    expectationArray = integratePSF(PSFunc, imageResolution, pixelScale, spotSpread, realCoordinates(1), realCoordinates(2));

 
%     figure(2)
%     pcolor(1:30, 1:30,expectationArray); xlabel( 'X' ); ylabel( 'Y' );
%     colormap(gray)
%     set(gca, 'YDir','reverse')
%     hold on;
%     plot(realCoordinates(1)+1,realCoordinates(2)+1, 'ro');
%     pause(1)
     
    %%% FLATTEN 2D ARRAY AND STORE ORIGINAL COORDINATES, CREATE NEW ARRAY
    %%% IN WHICH EACH ELEMENT IS DEFINED AS THE SUM OF ALL PREVIOUS
    %%% ELEMENTS. THIS WILL RANGE FROM 0 TO 1
    clear flatExpArray flatSum(i)
    i=1;
    for x = 1:imageResolution(1)
        for y = 1:imageResolution(2)
            flatExpArray(i) = expectationArray(x,y);
            flatSum(i) = sum(flatExpArray);
            flatSumPos(1,i) = x;
            flatSumPos(2,i) = y;
            i=i+1;
        end               
    end    
    %%% GENERATES RND NUM & USES FLATTENED SUM BETWEEN 0&1 IN CONJUCTION
    %%% WITH THE STORED COORDINATES OF THE 2D ARRAY TO DETERMINE WHERE IN
    %%% THE CDF THE PHOTON LANDS
    for pn = 1:photonNumber       
        [~,minPos] = min(abs(flatSum-rand));
        xPos = flatSumPos(1,minPos);
        yPos = flatSumPos(2,minPos);
        photonCount(xPos,yPos) = photonCount(xPos,yPos)+1;
    end    
    %%% QUANTUM EFFICIENCY 
    photoelectronCount = measurePhotolectrons(photonCount, quantumEfficiency, imageResolution);
end 

%%% GAUSSIAN NOISE
imageNoise = gaussianNoise(imageResolution,noiseMean,noiseVariance);

%%% IMAGE CREATION
image = imageNoise + photoelectronCount;
imgDat = image;

% figure(3)
% pcolor(1:30, 1:30,image); xlabel( 'X' ); ylabel( 'Y' );
% colormap(gray)
% set(gca, 'YDir','reverse')
% hold on;
% plot(realCoordinates(1)+1,realCoordinates(2)+1, 'ro');
% pause(1)

% %%% DISPLAY IMAGE
% % SCALED TO FIT MAXIMUM BRIGHTNESS
% 
% figure(3)
% imshow(photoelectronCount, [0 max(max(photonCount))]);
% title('Photoelectron Count');
% truesize(figure(3), [300,300]);
% 
% figure(4)
% imshow(image, [0 max(max(image))]);
% title('Photoelectron Count w/ Noise');
% truesize(figure(4), [300,300]);
% 
% %%% EXPORT IMAGE IN TIF FORMAT
% name='IMG1.tif'; location = 'imgDump';
% f = fullfile(location,name);
% image=mat2gray(image, [0 max(max(image))]);
% imwrite(image,f,'tif','Compression','none');



