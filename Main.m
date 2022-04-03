close all; clear all; clc;

display(strcat('Completion: 0%'))

variances =[50];
for var =1:length(variances)

    %%% SET PARAMETERS
    
        %%% SPOT PARAMS
        SpotParam.numberOfSpots =1; % Number of spots in image  
        SpotParam.spotSpreadRange = [1.2,1.2]; % Sigma value in PSF
        SpotParam.avgSpotPhotonNumber = 2e3;% Number of photons per spot in some shutter time
    
        %%% DEFINE IMAGE PARAMETERS
        ImageParam.quantumEfficiency = 0.9; % as a percentage
        ImageParam.imageResolution = [30,30]; % in pixels [height,width] 
        ImageParam.pixelScale = 1; % in unit distance/pixel [ie pixel scale (Astronomy), Pixel size (Microscopy)]
            ImageParam.imageSize = [ImageParam.pixelScale*ImageParam.imageResolution(1), ImageParam.pixelScale*ImageParam.imageResolution(2)]; %in units of distance
        ImageParam.noiseMean =variances(var); % Mean of gaussian noise distributed through image
        ImageParam.noiseVariance = variances(var); % Variance in gaussian noise, kept as mean to approximate poissonian distribution       
         
        numImages = 1000; % Number of images in stack
        
    %%% PASS PARAMETERS INTO SIMULATION
    
    %Get 10 seperate images
    for i = 1:numImages
	    [imageData{i},spotCoords{i}] = DripModelSpotSimulation(SpotParam,ImageParam);        
    end
    display(strcat('Completion: ',num2str(100*var/length(variances)),'%'))
    
    %%% EXPORT IMAGE IN TIF FORMAT
    nameNum=variances(var);
    name=strcat('IMG',num2str(nameNum),'_',num2str(SpotParam.avgSpotPhotonNumber),'_r10.tif'); location = 'ImageFolder';
    %name=strcat('TEST',num2str(nameNum),'_',num2str(SpotParam.avgSpotPhotonNumber),'.tif'); location = 'imgDump';
    f = fullfile(location,name);
    delete(f); % deletes file from previous run
    t = Tiff(f,'a'); % creating new Tiff file
        
        tagstruct.BitsPerSample = 16; % 16 bit image
        tagstruct.Compression = Tiff.Compression.None; %no compression
        tagstruct.Photometric = 1; %BlackIsZero. For bilevel and grayscale images: 0 is imaged as black. 2^BitsPerSample-1 is imaged as white.
        tagstruct.SamplesPerPixel = 1;% SamplesPerPixel  is  1  for bilevel, grayscale, and palette color images.  SamplesPerPixel is 3 for RGB images.
        tagstruct.RowsPerStrip = 1; %The image data is organized into strips for  fast access  to individual  rows  when  the  data  is compressed - though this field is valid even  if the  data is not compressed.
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % irrelevant for greyscale, but needs to be set
        tagstruct.Software = 'MATLAB';
        
    for i = 1:numImages
        img = im2uint16(mat2gray(imageData{i},[0 65535])); %converting to 16 bit grayscale, 1 bit corresponds to 1 photoelectron
        tagstruct.ImageLength = size(img,1);
        tagstruct.ImageWidth = size(img,2);
        
        desc = append(' IMG',num2str(i),',xo:',num2str(spotCoords{1,i}(2),6),' yo:',num2str(spotCoords{1,i}(1),6),' Scale:',num2str(ImageParam.pixelScale));
        tagstruct.ImageDescription = desc;
    
        setTag(t,tagstruct); % set Tiff tags
        write(t,img); % save image to Tiff file
        writeDirectory(t); % incriment image directory to next position in the stack
    % 
    %     figure(i+2)
    %     imshow(img, [0 max(max(img))]);
    %     truesize(figure(i+2), [300,300]);
    
    end
    close(t);
end

