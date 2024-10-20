% spotsQCFromManual
% clear
% makeGloblZombieSpots
matlabtempfiles='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles'
load(fullfile(matlabtempfiles,'merfishMask.mat'))

matlabzombie='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles/'
mypath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/data'
nucleipath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/nuclei'
% load(fullfile(matlabzombie,'otherglobalzombiespotsregionSmallRotation.mat'));%loads globalzombiespots
% % load((fullfile(matlabzombie,'binaryDapiregion0.mat')));%loads a dapimosaic variable
% load(fullfile(matlabzombie,'otherglobalzombiespotsregionSmallRotationFIXEDfromHW2hhmi.mat'));%loads globalzombiespots

% %% transform mask
% load(fullfile(matlabzombie,'tranformMtoZReal.mat'))%loads tform
% mymask=bwlabel(mymask>0);
% mymask=flipdim(mymask,1);
% rotatedmask = imwarp_same(mymask,tformlarge);
% imwrite(reduceImage(reduceImage(reduceImage(rotatedmask))),fullfile(matlabzombie,'shrunkrotatedmask.tif'))
% save(fullfile(matlabzombie,strcat('rotatedmask.mat')),'rotatedmask')
load(fullfile(matlabzombie,strcat('rotatedmask.mat')))
disp('loaded mask')
%%
numberofcells=max(max(rotatedmask));
% globalzombiespots=globalzombiespots(globalzombiespots(:,end-3)>0&globalzombiespots(:,end-2)>0&globalzombiespots(:,end-3)<size(rotatedmask,1)&globalzombiespots(:,end-2)<size(rotatedmask,2),:);
% size(globalzombiespots)
% ids=sub2ind([size(rotatedmask,1),size(rotatedmask,2)],globalzombiespots(:,end-3),globalzombiespots(:,end-2));
% globalzombiespots(:,11)=ids;
% save('/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles/globalZspotswithcoordinatesFIXEDHW2UGER.mat','globalzombiespots','-v7.3')
load('/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles/globalZspotswithcoordinatesFIXEDHW2UGER.mat')
disp('loaded spots')
% zombieimage=zeros(size(rotatedmask));
% for i=1:size(ids,1)
%     zombieimage(ids(i))=globalzombiespots(i,4);
% end
% imwrite(reduceImage(reduceImage(reduceImage(imdilate(zombieimage,strel('disk',30))))),fullfile(matlabzombie,'mappedzomiespotsR0.tif'))
% size(zombieimage)
% size(rotatedmask)

% cellstats=regionprops(rotatedmask>0,rotatedmask+0,'Area','PixelIdxList','MeanIntensity');
for i=1:numberofcells
    if rem(i,300)==0
        disp(strcat('current cell is ',num2str(i), 'out of the ', num2str(size(numberofcells,1))))%progress report
    end
    if rem(i,20000)==0
        save(strcat('/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles/globalZspotswithIDFIXEDHW2UGER_',num2str(i),'.mat'),'globalzombiespots','-v7.3')
    end
    mypixels=find(rotatedmask==i);
    positions=ismember(globalzombiespots(:,11),mypixels);
    globalzombiespots(positions,12)=i;%add a cell ID to each transcript
end
disp('start saving')
save('/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles/globalZspotswithIDFIXEDHW2UGER.mat','globalzombiespots','-v7.3')
% load('/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles/globalZspotswithIDFIXEDHW2UGER.mat')
disp('start counttable')
countable=zeros(numberofcells,2+max(globalzombiespots(:,2)));
%% new code
zombietableR0=zeros(numberofcells,max(globalzombiespots(:,4)));
for i=1:size(globalzombiespots,1)
    if rem(i,300)==0
        disp(strcat('current cell is ',num2str(i), 'out of spots ', num2str(size(globalzombiespots,1))))%progress report
    end
    if (globalzombiespots(i,12)>0)
        zombietableR0(globalzombiespots(i,12),globalzombiespots(i,4))=1;
    end
end
%% end new code

writematrix(zombietableR0,fullfile(matlabzombie,'otherzombiecounttableNEWR0FIXEDHW2.csv'))
% clear
% disp('finding neighbors')
% findtouchingcellstumor
