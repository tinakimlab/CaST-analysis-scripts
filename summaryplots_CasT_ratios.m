%% Change the filefolder, filename_base, adjust xlim for fig3
clear all
close all

filepath='/Users/rz/Library/Matlab/Data/';
filefolder='230105_CaST_IRES_repeat_copy/results/';
filename_base='IRES_repeat_ds';
filename=strcat(filepath,filefolder,filename_base,'_FOV_bgsub_ratios.mat');
load(filename);


nFOV=size(FOV_ratios,2);
nCon=size(FOV_ratios,1);

% 96 pixels per 2.54cm
% 1 pixel = 0.265 mm
% 480x360 pixel images
%

% figure1
figure;
hold on;
for a=1:nCon
    tempNoCa=FOV_ratios(a,1:nFOV/2);
    tempCa=FOV_ratios(a,nFOV/2+1:end);

    bar([a-0.15 a+0.15],[mean(tempNoCa) mean(tempCa)],'facecolor', 'k');
    errorbar([a-0.15 a+0.15],[mean(tempNoCa) mean(tempCa)],[std(tempNoCa)./sqrt(nFOV) std(tempCa)./sqrt(nFOV)],'k','linestyle','none');
    ratios(a)=mean(tempCa)/mean(tempNoCa);
    pvals(a)=ranksum(tempNoCa,tempCa);
    plot(a-0.15,tempNoCa,'ok');
    plot(a+0.15,tempCa,'or');
end

box off

%%
% figure2

filename2=strcat(filepath,filefolder,filename_base,'_bgsub_ratios.mat');
load(filename2);

figure;
ratio_noca=meanfluo_SA647(no_ca_inds)./meanfluo(no_ca_inds);
ratio_ca=meanfluo_SA647(ca_inds)./meanfluo(ca_inds);
bar([1 2],[mean(ratio_noca) mean(ratio_ca)]);
hold on
errorbar([1 2],[mean(ratio_noca) mean(ratio_ca)], [std(ratio_noca)./sqrt(length(ratio_noca)) std(ratio_ca)./sqrt(length(ratio_ca))]);
plot(1,meanratio(no_ca_inds),'ok');
plot(2,meanratio(ca_inds),'or');

%%

% figure3
figure;
hold on;
plot(meanfluo(ca_inds),meanfluo_SA647(ca_inds),'.r');

plot(meanfluo(no_ca_inds),meanfluo_SA647(no_ca_inds),'.k');
xlim([5 22])

%%

%figure4 5
data=[meanfluo_SA647./meanfluo];
labels=zeros(1,length(num_cells));
labels(ca_inds)=1;
figure;
vs=violinplot(data,labels);
% TEST CASE 4
figure;
%C = colororder;
C=[1 0 0; 0 0 0];
vs4 = violinplot({meanfluo_SA647(ca_inds)./meanfluo(ca_inds),meanfluo_SA647(no_ca_inds)./meanfluo(no_ca_inds)},labels,'ViolinColor',{C,C(2,:)},'ViolinAlpha',{0.3 0.3}, 'ShowMean', false);
hold on
plot([0.7 1.3],[mean(meanfluo_SA647(no_ca_inds)./meanfluo(no_ca_inds)) mean(meanfluo_SA647(no_ca_inds)./meanfluo(no_ca_inds))],'-k');
plot([0.7 1.3],[mean(meanfluo_SA647(ca_inds)./meanfluo(ca_inds)) mean(meanfluo_SA647(ca_inds)./meanfluo(ca_inds))],'-r');

%%

% figure6
thresh=prctile(meanfluo_SA647(no_ca_inds),90);
FOV_thresh_noca=[];
FOV_thresh_ca=[];

nFOV=size(FOV_ratios,2);

for b=1:nFOV
    curr_SA647=meanfluo_SA647(cell_FOV_ind==b);
    curr_above=length(find(curr_SA647>thresh))/length(curr_SA647);
    
    if b<(nFOV/2)+1
        FOV_thresh_noca=[FOV_thresh_noca; curr_above];
    else %Ca_inds
        FOV_thresh_ca=[FOV_thresh_ca; curr_above];
    end
end

figure;
swarmchart(ones(1,length(FOV_thresh_noca)),FOV_thresh_noca);
hold on
plot([0.75 1.25],[mean(FOV_thresh_noca) mean(FOV_thresh_noca)]);
swarmchart(ones(1,length(FOV_thresh_noca))*2,FOV_thresh_ca);
plot([1.75 2.25],[mean(FOV_thresh_ca) mean(FOV_thresh_ca)]);

% figure7
figure;
bar([1 2],[mean(FOV_thresh_noca) mean(FOV_thresh_ca)]);
hold on
errorbar([1 2],[mean(FOV_thresh_noca) mean(FOV_thresh_ca)], [std(FOV_thresh_noca)./sqrt(length(FOV_thresh_noca)) std(FOV_thresh_ca)./sqrt(length(FOV_thresh_ca))]);
plot(1,FOV_thresh_noca,'ok');
plot(2,FOV_thresh_ca,'or');
