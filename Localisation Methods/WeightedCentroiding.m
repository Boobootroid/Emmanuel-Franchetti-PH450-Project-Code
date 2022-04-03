clear all;close all;clc;
% Performs weighted centroiding on stack of spots, then shows analysis

%% Set Parameters
% parameters used in spot simulation
imgSigmas =[1,1.2,1.4,1.6,1.8,2]; % width of PSF
variance = 50; % mean and variance in poissonian background noise
    noiseMean = variance;
numFilters = 0; % number of processed images
QE=0.9; % quantum efficiency
photonNumber = 2e3; % number of photons hitting photodetector

procType = 'Gaussian'; % if anything other than 'Gaussian' will not modify threshold

% Pre load Arrays
unProcMeanMiss = zeros(1, length(imgSigmas));
unProcError = zeros(1, length(imgSigmas));
procMeanMiss = zeros(length(imgSigmas),numFilters);
procError = zeros(length(imgSigmas),numFilters);

%% Analyse Unprocessed Image
for var =1:length(imgSigmas)
    %noiseMean = imgSigmas(var);
    imgSigma = imgSigmas(var);
    nameNum = imgSigmas(var)*10;
    imgName = strcat('IMG50_',num2str(photonNumber),'_r',num2str(nameNum),'.tif'); imgName2 = erase(imgName,'.tif');
 
    threshold = noiseMean+4*sqrt(variance);
    imageinf =imfinfo(imgName); numImages = length(imageinf);
    
    % Pre load arrays
    unProcMiss = zeros(1,numImages);
    unProcMissPerc = zeros(1,numImages);

    for i = 1:numImages
        image = imread(imgName,i);

        % Threshold image for centroiding
        BW = image > threshold;
        % Perform weighted centroid
        rp = regionprops(BW,image,{'Centroid','WeightedCentroid'});
        cent = rp.WeightedCentroid;

        %Extract Ground Truth from image metadata
        desc = imageinf(i).ImageDescription; pat=digitsPattern; metadata = extract(desc,pat);  
        gt.xo(i) = str2num(metadata{2})+str2num(metadata{3})*10^(-length(metadata{3}));
        gt.yo(i) = str2num(metadata{4})+str2num(metadata{5})*10^(-length(metadata{5}));        
        gt.scale = str2num(metadata{6});
        %gt.xo(s) = 15; gt.yo(s)=15; gt.scale =1;
        
        %Calculate how much the fit misses the ground truth 
        unProcMiss(i) = norm([cent(1) cent(2)]-[gt.xo(i)/gt.scale gt.yo(i)/gt.scale]);
        unProcMissPerc(i) = 100*unProcMiss(i)/imgSigma;

%         %Option to Display Images
        if i==1&& var ==1
            figure(1)
            pcolor(image); 
            xlabel( 'X' ); ylabel( 'Y' ); title(strcat('Image ',num2str(i)));
            set(gca, 'YDir','reverse')
            hold on;  colormap(gray)
            plot(cent(1)+1,cent(2)+1,'ro'); plot((gt.xo(i)+1)/gt.scale,(gt.yo(i)+1)/gt.scale,'bo');
            fitstring = strcat('Weighted Centroid, xo:',num2str(cent(1)),',yo:',num2str(cent(2)));
            gtstring = strcat('Ground Truth, xo=',num2str(gt.xo(i)/gt.scale),',yo=',num2str(gt.yo(i)/gt.scale));
            legend('Photoelectrons',fitstring ,gtstring);
            pause(0.5)         
            
            figure(2)
            pcolor(BW); 
            xlabel( 'X' ); ylabel( 'Y' ); title(strcat('Image ',num2str(i), ', Threshold:', num2str(threshold)));
            set(gca, 'YDir','reverse')
            hold on;  colormap(gray)
            plot(cent(1)+1,cent(2)+1,'ro'); plot((gt.xo(i)+1)/gt.scale,(gt.yo(i)+1)/gt.scale,'bo');
            fitstring = strcat('Weighted Centroid, xo:',num2str(cent(1)),',yo:',num2str(cent(2)));
            gtstring = strcat('Ground Truth, xo=',num2str(gt.xo(i)/gt.scale),',yo=',num2str(gt.yo(i)/gt.scale));
            legend('',fitstring ,gtstring);
            
            pause(0.5)
        end
    end 
    %Calculate mean miss over stack for one threshhold 
    unProcMeanMiss(var) = mean(unProcMissPerc);  
    unProcError(var) = std(unProcMissPerc)/sqrt(numImages);
    
    %% Analyse Processed Images

    % Pre load arrays   
    procMiss = zeros(numImages, numFilters);
    procMissPerc = zeros(numImages, numFilters);

    % Threshold initially 4 times noise variance plus mean
    threshold = noiseMean+(4*sqrt(imgSigmas(var)));
    for n = 1:numFilters
        % Image filename
        str = strcat(imgName2,{'x'},num2str(n),'.tif');
        % Checks if image has applied enough gaussian blur to reduce
        % threshhold
        if (QE*photonNumber)/(2*pi*(imgSigma+0.1*n)^2)+noiseMean < threshold && strcmp(procType,'Gaussian') == 1 
            threshold =  QE*photonNumber*(2*pi*(imgSigma+0.1*n)^2)^-1+noiseMean
        end
            for i = 1:numImages           
                image = imread(str{1},i);
                % Threshold image for centroiding
                BW = image > threshold;
                % Perform weighted centroid
                rp = regionprops(BW,image,{'Centroid','WeightedCentroid'});
                cent = rp.WeightedCentroid;
                
                % Option to display image
%                 if s==1
%                     figure(2)
%                     pcolor(BW); 
%                     xlabel( 'X' ); ylabel( 'Y' ); title(strcat('Image ',num2str(s), ', Filter:', num2str(n)));
%                     set(gca, 'YDir','reverse')
%                     hold on;  colormap(gray)
%                     plot(cent(1)+1,cent(2)+1,'ro'); plot((gt.xo(s)+1)/gt.scale,(gt.yo(s)+1)/gt.scale,'bo');
%                     fitstring = strcat('Weighted Centroid, xo:',num2str(cent(1)),',yo:',num2str(cent(2)));
%                     gtstring = strcat('Ground Truth, xo=',num2str(gt.xo(s)/gt.scale),',yo=',num2str(gt.yo(s)/gt.scale));
%                     legend('',fitstring ,gtstring);
%                 end

                %Calculate miss from ground Truth
                procMiss(i) = norm([cent(1) cent(2)]-[gt.xo(i)/gt.scale gt.yo(i)/gt.scale]);
                procMissPerc(i) = 100*procMiss(i)/imgSigma;
            end


        % take mean and std of percentage misses from fit
        procMeanMiss(var,n) = mean(procMissPerc);    
        procError(var,n) = std(procMissPerc)/sqrt(numImages);
    end

end

%% Plot processed image effect 
if numFilters >0
    figure(3)
    hold on;
    for var=1:length(imgSigmas)   
        plot(1:numFilters,unProcMeanMiss(var)-procMeanMiss(var,:))
    end
    legend(num2str(imgSigmas(1)),num2str(imgSigmas(2)),num2str(imgSigmas(3)),num2str(imgSigmas(4)),num2str(imgSigmas(5)))
    if strcmp(procType,'Gaussian') == 1 
        xlabel('Gaussian blur sigma')
    else 
        xlabel('Number of MNRs applied')
    end
    ylabel( 'Mean Increase in accuracy over 25 images [% of PSF Radius]');
    title(strcat('Mean Effect of ',procType,' on Centroiding Accuracy'));
end

%% Plot Unprocessed performance with noise
figure(4)
errorbar(imgSigmas, 100-unProcMeanMiss,unProcError)
xlabel('PSF Radius');ylabel('Acuracy as a % of PSF Radius'); title('Weighted Centroiding Accuracy with Noise');

% Save data to file
save(strcat('WCmeanMiss_4k_r',num2str(imgSigmas(var)),'.mat'),'unProcMeanMiss');
save(strcat('WCErrors_4k_r',num2str(imgSigmas(var)),'.mat'),'unProcError');


