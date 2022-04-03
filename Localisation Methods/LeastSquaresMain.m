clear all;close all;clc;
warning('off')
display(strcat('Progress: 0%')) % Progress counter

%% Set Parameters
% Calls Least Squares Fitting on stack of spots, performs analysis
% parameters used in spot generation
variances =[0,25,50,75,100,150,200,250,300,350,400]; % noise variance in spot stacks
imgSigma = 1.2; % PSF sigma
numFilters = 0; % Number of processed stacks
QE = 0.9; % Quantum Efficiency

procType = 'MNR'; % if anything other than 'Gaussian' will not modify guess sigma
    
% Pre load Arrays
unProcMeanMiss = zeros(1, length(variances));
unProcError = zeros(1, length(variances));
procMeanMiss = zeros(length(variances),numFilters);
procError = zeros(length(variances),numFilters);

% Loop through each stack
for var = 1:length(variances)
    photonNumber = 2e3; % Photon number set in spot simulation
    noiseVariance = variances(var);
    %% analyse unprocessed stack
    nameNum=variances(var)*10;
    imgName = strcat('IMG',num2str(noiseVariance),'_',num2str(photonNumber),'.tif'); % Image title based on variance & photon number
    %noiseVariance = imgSigmas(var);
    photonNumber = QE*photonNumber; %quantum efficiency reducing total number of photoelectrons
    imageinf =imfinfo(imgName); numImages = length(imageinf); %finding length of stack
    %Get size of image 
    image = imread(imgName,1);
    [X,Y] = size(image);
    
    %Establish guess parameters for fit
    % N
    lowerPhotonNumber=photonNumber-200; startPhotonNumber=photonNumber; upperPhotonNumber=photonNumber+200;
    % xo and yo
    locDispersion = 2; 
        startX = X/2; lowerX = startX-locDispersion; upperX = startX + locDispersion;
        startY = Y/2; lowerY = startY-locDispersion; upperY = startY + locDispersion;
    % Sigma
    sigDisp = 0.2;
        lowerSigma = (imgSigma-sigDisp); startSigma = imgSigma; upperSigma = (imgSigma+sigDisp);
    
    % Pre load arrays
    lsqxo =  zeros(1,numImages);
    lsqyo = zeros(1,numImages);
    unProcMiss = zeros(1,numImages);
    unProcMissPerc = zeros(1,numImages);
    
    % Loop through each image in unprocessed stack
    for i = 1:numImages    
        image = imread(imgName,i);
        
        %Establish final guess parameter for fit, A
        minIm = min(min(image));
            lowerOffset=minIm; startOffset=minIm+2*sqrt(noiseVariance); upperOffset=minIm+4*sqrt(noiseVariance);
        
        %Guesses are: Offset,Photon Number, X0, Y0, Sigma
        guess.Lower = [lowerOffset,lowerPhotonNumber,lowerX,lowerY,lowerSigma]; % lower bounds of models parameters
        guess.StartPoint = [startOffset,startPhotonNumber,startX,startY,startSigma]; % starting parameter guess
        guess.Upper = [upperOffset,upperPhotonNumber,upperX,upperY,upperSigma]; % upper bounds
    
        %Perform Fit
        [lsqFit, gof, fout] = lsq2DGaussian(image, guess,0);
        %Store centroid
        lsqxo(i) = lsqFit.xo; lsqyo(i) = lsqFit.yo;
        %Extract Ground Truth from image metadata
        desc = imageinf(i).ImageDescription; pat=digitsPattern; metadata = extract(desc,pat);  
        gt.xo(i) = str2num(metadata{2})+str2num(metadata{3})*10^(-length(metadata{3}));
        gt.yo(i) = str2num(metadata{4})+str2num(metadata{5})*10^(-length(metadata{5}));
        gt.scale = str2num(metadata{6});
        
        %Distance from centroid to ground truth    
        unProcMiss(i) = norm([lsqxo(i) lsqyo(i)]-[gt.xo(i)/gt.scale gt.yo(i)/gt.scale]);
        %As a percentage of sigma
        unProcMissPerc(i) = 100*unProcMiss(i)/imgSigma;

%         %Option to show unprocessed image
        if i==1
            figure(1)
            pcolor(image); xlabel( 'X' ); ylabel( 'Y' ); title(strcat('Var:',num2str(variances(var)),', Image ',num2str(i)));
            set(gca, 'YDir','reverse')
            hold on;
            plot(lsqxo(i)+1,lsqyo(i)+1,'ro'); plot((gt.xo(i)+1)/gt.scale,(gt.yo(i)+1)/gt.scale,'bo');
    
            fitstring = strcat('LSQ Centre, xo:',num2str(lsqxo(i)),',yo:',num2str(lsqyo(i)));
            gtstring = strcat('Ground Truth, xo=',num2str(gt.xo(i)/gt.scale),',yo=',num2str(gt.yo(i)/gt.scale));
            legend('Pixel Count',fitstring ,gtstring);
            colormap(gray)
        end
   
    end

    % Mean miss of unprocessed stack 
    unProcMeanMiss(var) = mean(unProcMissPerc);
    unProcError(var) = std(unProcMissPerc)/sqrt(numImages);
    %% Analyse Processed stacks

    % Pre load arrays   
    procMiss = zeros(numImages, numFilters);
    procMissPerc = zeros(numImages, numFilters);
    
    % Loop through each processed stack
    for n=1:numFilters
        n
        if strcmp(procType,'Gaussian') == 1 
            str = strcat(erase(imgName,'.tif'),{'x'},num2str(n),'.tif');
        else
            str = strcat(erase(imgName,'.tif'),{'_MNR'},num2str(n),'.tif');
        end
        % Pre load arrays
        lsqxo = zeros(1,numImages);
        lsqyo = zeros(1,numImages);
        
        %Loop through each image in the stack
        for i = 1:numImages
            image = imread(imgName,i); %original image to compare with blurred image
            image2 = imread(str{1},i); %New blurred image
            minIm = min(min(image2)); 

            [X,Y] = size(image2);
            %Establish guess parameters for fit
            % A
            lowerOffset=minIm; startOffset=minIm+2*sqrt(noiseVariance); upperOffset=minIm+4*sqrt(noiseVariance);
     
            %Sigma increases as image blurs, n is 10*blur sigma, image not
            %affected by sigma less than 0.3
            if strcmp(procType,'Gaussian') == 1 && any(image2-image,'all') %Checks that image has been changed by filter
                % Modified sigma is Sum of two sigmas in quadrature
                lowerSigma = ((imgSigma)^2+(0.1*n)^2)^0.5-sigDisp; startSigma = (imgSigma^2+(0.1*n)^2)^0.5; upperSigma = ((imgSigma)^2+(0.1*n)^2)^0.5+sigDisp;
            end
            
            
            %Guesses are: Offset,Photon Number, X0, Y0, Sigma
            guess.Lower = [lowerOffset,lowerPhotonNumber,lowerX,lowerY,lowerSigma]; % lower bounds of models parameters
            guess.StartPoint = [startOffset,startPhotonNumber,startX,startY,startSigma]; % starting parameter guess
            guess.Upper = [upperOffset,upperPhotonNumber,upperX,upperY,upperSigma]; % upper bounds
    
            %Perform Fit         
            [lsqFit, gof, fout] = lsq2DGaussian(image2, guess,0);
            
            %Store centroid
            lsqxo(i) = lsqFit.xo; lsqyo(i) = lsqFit.yo;

            %Distance from centroid to ground truth
            procMiss(i,n) = norm([lsqxo(i) lsqyo(i)]-[gt.xo(i)/gt.scale gt.yo(i)/gt.scale]);
            %As a percentage of PSF Sigma
            procMissPerc(i,n) = 100*procMiss(i,n)/imgSigma;

            % Optional Display of Processed image
            if n == 1 && i==1  
                figure(2)
                pcolor(image); xlabel( 'X' ); ylabel( 'Y' ); title(strcat('Image: ',num2str(i), ', Gaussian Blur: 2'));
                set(gca, 'YDir','reverse')
                hold on;
                plot(lsqxo(i)+1,lsqyo(i)+1,'ro'); plot((gt.xo(i)+1)/gt.scale,(gt.yo(i)+1)/gt.scale,'bo');
        
                fitstring = strcat('LSQ Centre, xo:',num2str(lsqxo(i)),',yo:',num2str(lsqyo(i)));
                gtstring = strcat('Ground Truth, xo=',num2str(gt.xo(i)/gt.scale),',yo=',num2str(gt.yo(i)/gt.scale));
                legend('',fitstring ,gtstring);
                colormap(gray)
                pause(0.5)
            end            
        
        end    
        %find mean difference in fit-to-ground truth of processed vs
        %unprocessed images. Averaged over every image in stack
        procMeanMiss(var,n) = mean(unProcMissPerc-rot90(procMissPerc(:,n)));
        %Standard error in mean of above
        procError(var,n) = std(unProcMissPerc-rot90(procMissPerc(:,n)))/sqrt(numImages);

        % Optional break if change in accuracy becomes negative
%         if strcmp(procType,'Gaussian') == 1   
%             if meanPerc(var,n) <0 && n>20% If change in accuracy is -ve 
%                 display(strcat('Break at n=',num2str(n)))
%                 break % break out of loop
%             end
%         end
    end       
    display(strcat('Progress: ',num2str(100*var/length(variances)),'%')) % Progress counter
end

%% Plot Unprocessed performance with noise
% figure(4)
% errorbar(variances, 100-unProcMeanMiss,unProcError); hold on;
% xlabel('PSF Radius');ylabel('Acuracy as a % of PSF radius'); title('Unprocessed LSQ Accuracy with PSF radius, 2000 photons');
% 
% figure(11)
% legend('Noise Variance: 0','25','50','75','100','200','400')

%% Export results with standard errors
save(strcat('LSQmeanMiss_Unprocessed_',num2str(variances(var)),'.mat'),'unProcMeanMiss');
save(strcat('LSQmeanMiss_Unprocessed_Errors_',num2str(variances(var)),'.mat'),'unProcError');
if numFilters >0
    save(strcat('LSQmeanMiss_MNR_',num2str(variances(var)),'.mat'),'procMeanMiss');
    save(strcat('LSQmeanMiss_MNR_Errors_',num2str(variances(var)),'.mat'),'procError');
end




