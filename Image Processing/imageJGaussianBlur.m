clear all; close all; clc;
addpath 'D:\Users\Mani\fiji-win64\Fiji.app\scripts' % Update for your ImageJ2 (or Fiji) installation as appropriate
ImageJ;
% Applies increasingly strong gaussian blurs

variances =[75];
photonNumber = 2e3;
for i = 1:length(variances)
    display(strcat('Progress: ',num2str(100*i/length(variances)),'%'))
    nameNum=variances(i);
    src = strcat('IMG',num2str(nameNum),'_',num2str(photonNumber),'.tif'); 
    % Loop
    for j=1:40
        %Open Image in IM
        imp = ij.IJ.openImage(src);
        imp.show();
        % Blur image
        sigma= 0.1*j;
        str = strcat("sigma=",num2str(sigma)," stack");
        ij.IJ.run("Gaussian Blur...", str);
        str = strcat("IMG",num2str(nameNum),'_',num2str(photonNumber),"x",num2str(j),".tif");
        ij.IJ.saveAs("Tiff", str);
        imp.close();
    end
end