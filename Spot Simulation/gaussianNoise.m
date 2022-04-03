%% Distributes additive gaussian white noise

function noise = gaussianNoise(dim,mean,var)
    sig = sqrt(var); % define std deviation
    
    % create probability distribution
    pd = makedist('Normal','mu',mean,'sigma',sig);
    % create range of distribution (mean-3sigma to mean+3sigma)
    noiseLevel = round(mean-3*sig):round(mean+3*sig);
    % create cumulative probability distribution
    noiseProb = cdf(pd,noiseLevel);    
    
    for x = 1:dim(1)
       for y = 1:dim(2)
           %Use cdf to randomly select noise level
           [~,noiseLevelPos] = min(abs(noiseProb-rand));
           % will only accept minimum of 0
           if noiseLevel(noiseLevelPos) > 0
                noise(x,y) = noiseLevel(noiseLevelPos);                 
           else
                noise(x,y)=0;
           end
       end
    end
end
