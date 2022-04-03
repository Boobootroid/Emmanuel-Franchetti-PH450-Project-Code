clear all; close all; clc;
addpath 'D:\Users\Mani\fiji-win64\Fiji.app\scripts' % Update for your ImageJ2 (or Fiji) installation as appropriate
ImageJ;
% repeats MNR multiple times on same image

variances =[0,25,50,75,100,150,200,250,300,350,400];
photonNumber = 2e3;
for i = 1:length(variances)
    nameNum=variances(i);
    src = strcat('IMG',num2str(nameNum),'_',num2str(photonNumber),'.tif'); 
    % Loop
    for j=1:1 % number of MNR applied
        j
        %Open Image in IM
        if j==1
            imp = ij.IJ.openImage(src);
        else
            imp = ij.IJ.openImage(str);
        end
        imp.show();
        % Blur image
        r= 0.5; % blur radius in pixels
        str = strcat("radius=",num2str(r)," stack");
        ij.IJ.run("Median...", str);
    %     imp.show()
        str = strcat("IMG",num2str(nameNum),'_',num2str(photonNumber),"_MNR",num2str(j),".tif");
        ij.IJ.saveAs("Tiff", str);
        imp.close();
    end
end