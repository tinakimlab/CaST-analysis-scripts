%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================

%% change filepath, filefolder, nFOV base on the data, create 'results' folder

% -- FIRST, combine all FOVs horizonally in ImageJ/Fiji -- %
function GFP_SA647_bgsubtract_cell(filename_base)

%clear all
%close all

% load the image tif files data
filepath='/Users/rz/Library/Matlab/Data/';
filefolder='230105_CaST_IRES_repeat/';
filename_GFP=strcat(filename_base,'_GFP.tif');

imsegm=imread(strcat(filepath,filefolder,filename_GFP));
imsegm=double(imsegm);
filename_SA647=strcat(filename_base,'.tif');
imsegm_SA647=imread(strcat(filepath,filefolder,filename_SA647));
imsegm_SA647=double(imsegm_SA647);


% This assumes FOVs are combined in 1 row x # of columns.
nFOV=10;
nCON=1;
imx=size(imsegm,2)/nFOV; % width
imy=size(imsegm,1)/nCON; % height
FOV_ratios=zeros(nCON,nFOV);


% Segmentation by adaptive thresholding
prm.method = 'adth';
prm.adth.th = 0.01; %adjust thresholding to improve results

[cellbw1,wat,imsegmout,prmout] = cellsegm.segmct(imsegm,0.01,500,'prm',prm);
cellsegm.show(imsegmout,1);title('Raw image');axis off;
cellsegm.show(cellbw1,1);title('Cell segmentation by ADTH');axis off;

% improving the results by splitting of cells
splitth = 1;
plane = 1;
% cells above this threshold are split (all cells here)
n = 100;
h = [0.5 0.5 1.5];
cellbw2 = cellsegm.splitcells(cellbw1,splitth,n,h);
cellsegm.show(cellbw2,2);title('Cell segmentation by ADTH with splitting');axis off;
%cellbw2=cellbw1;

% create cellmasks
[L,num_cells]=bwlabel(cellbw2);
segcentroids=zeros(num_cells,2);
meanfluo=zeros(num_cells,1);
meanfluo_SA647=zeros(num_cells,1);
meanratio=zeros(num_cells,1);
cell_FOV_ind = [];

for a=1:num_cells
    [tempx,tempy]=find(L==a);
    fluo=[]; % calculate each cell's GFP fluorescence
    fluo_SA647=[]; % calculate each cell's SA647 fluroescence
    cellmasktemp=zeros(size(cellbw2,1),size(cellbw2,2));
    for b=1:length(tempx)
        cellmasktemp(tempx(b),tempy(b))=1;
        fluo(b)=imsegm(tempx(b),tempy(b));
        fluo_SA647(b)=imsegm_SA647(tempx(b),tempy(b));
    end
    temp=regionprops(cellmasktemp);
    segcentroids(a,:)=temp.Centroid;
    meanfluo(a)=mean(fluo);
    meanfluo_SA647(a)=mean(fluo_SA647);
    meanratio(a)=mean(fluo_SA647)./mean(fluo); % calculate average mCh/GFP ratio
    cell_FOV_ind(a)=ceil(temp.Centroid(1)/imx);
end

no_ca_inds=find(segcentroids(:,1)<=imx*nFOV/2);
ca_inds=find(segcentroids(:,1)>imx*nFOV/2);

% save data 
savename=strcat(filename_base,'_bgsub_ratios.mat');
save(strcat(filepath,filefolder,'results/',savename),'num_cells','meanfluo','meanfluo_SA647','meanratio','no_ca_inds','ca_inds','cell_FOV_ind');

%% Calculate the average ratios per FOV

for z=1:nCON % number of rows in combined image
    temp=find(segcentroids(:,2)>=(z-1)*imy+1 & segcentroids(:,2)<(z)*imy); %finds cells corresponding to current row
    tempx=segcentroids(temp,1); %gets x position of each cell in current row
    
    imsegm_row=imsegm_SA647((z-1)*imy+1:(z)*imy,:);
    
    Biotin_NoCa=[];
    Biotin_Ca=[];
    
    for a=1:nFOV % number of columns in combined image
        
        tempa=find(tempx>=(a-1)*imx+1 & tempx<(a)*imx); %finds cells in current row for current column
        currmask=cellbw2((z-1)*imy+1:z*imy,(a-1)*imx+1:a*imx);
        currim=imsegm_row(:,(a-1)*imx+1:a*imx);
        %currim_bg=prctile(currim(~currmask),8); %bg subtract base on percentile
        currim_bg=mean(currim(~currmask)); %bg subtract base on mean
        %currim_bg=0; %no bg subtract (when bg is subtracted in previous steps)
        
        GFP=meanfluo(temp(tempa));
        SA647=meanfluo_SA647(temp(tempa))-currim_bg;
        
        if a<(nFOV/2)+1
            Biotin_NoCa=[Biotin_NoCa; SA647./GFP];
        else
            Biotin_Ca=[Biotin_Ca; SA647./GFP];
        end
        
        FOV_values_SA647(z,a)=mean(SA647);
        FOV_values_GFP(z,a)=mean(GFP);
        FOV_ratios(z,a)=mean(SA647)/mean(GFP);
        %FOV_cellcount(z,a)=length(tempa);
        
    end
end

savename=strcat(filename_base,'_FOV_bgsub_ratios.mat');
save(strcat(filepath,filefolder,'results/',savename),'FOV_ratios','Biotin_Ca','Biotin_NoCa','FOV_values_SA647','FOV_values_GFP');

display('done!')

end
