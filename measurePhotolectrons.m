%% Simulates Quantum efficiency in pixel array

function photoElectronCount = measurePhotolectrons(photonCount,quantumEfficiency,imageResolution)

    %%% QUANTUM EFFICIENCY 
    photoElectronCount = zeros(imageResolution(1),imageResolution(2));
    for x = 1:imageResolution(1)
        for y = 1:imageResolution(2)
            for pn = 1:photonCount(x,y)
               rnd = rand;
               if rnd <= quantumEfficiency
                   photoElectronCount(x,y) = photoElectronCount(x,y)+1;
               end
            end
        end               
    end
end