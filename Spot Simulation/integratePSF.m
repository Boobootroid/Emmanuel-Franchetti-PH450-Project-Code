function intensityArray = integratePSF(PSF, imageResolution, pixelScale, sigma, Xpos,Ypos)


%%% CALCULATE MEASURED INTENSITY AT EACH PIXEL FOR EACH SPOT
    
    intensityArray = zeros(imageResolution(1),imageResolution(2));

    for x = 0:imageResolution(1)-1
        scaledX = x*pixelScale; %Converting into real coords for integration limits
            for y = 0:imageResolution(2)-1
                scaledY = y*pixelScale;
                    %upper limit set by real coord plus minRes, ie the size of one 
                    %pixel in nm
                    int = integral2(PSF,scaledX,scaledX+pixelScale,scaledY,scaledY+pixelScale);  
                    intensityArray(y+1, x+1) = int; 
            end    
    end

    return  
end
