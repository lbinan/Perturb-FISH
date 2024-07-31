matlabmerfish='/broad/clearylab/Users/Loic/ASDMerfishRevs4/matlabtempfiles'
matlabzombie='/broad/clearylab/Users/Loic/ASdRevsZombiedish4/matlabtempfiles'
outputs='/broad/clearylab/Users/Loic/ASdRevsZombiedish4/outputs'
load(fullfile(matlabmerfish,'mymaskregion0V2.mat'));%loads  mymask
load(fullfile(matlabzombie,'globalzombiespotsregion0.mat'));%loads globalzombiespots
load((fullfile(matlabzombie,'binaryDapiregion0.mat')));%loads a dapimosaic variable

%% transform mask
load(fullfile(outputs,'tformMerfishstraighttozombiesample4.mat'))%loads tform
mymask=bwlabel(mymask>0);
mymask=flipdim(mymask,1);
rotatedmask = imwarp_same(mymask,tform);
imwrite(reduceImage(reduceImage(reduceImage(rotatedmask))),fullfile(matlabzombie,'rotatedmaskR0.tif'))
imwrite(reduceImage(reduceImage(reduceImage(dapimosaic(1:size(rotatedmask,1),1:size(rotatedmask,2))))),fullfile(matlabzombie,'matchingzombienuclei.tif'))
% cropzombie=zombieimage(5000:15000,5000:15000);
% croprotatedmask=rotatedmask(5000:15000,5000:15000);
% imwrite(imdilate(cropzombie>0,strel('disk',10)),fullfile(matlabzombie,'cropzombieregion0.tif'))
% imwrite(croprotatedmask,fullfile(matlabzombie,'croprotatedmaskregion0.tif'))

%%
numberofcells=max(max(rotatedmask));
croppedspots=globalzombiespots(globalzombiespots(:,end-3)<size(rotatedmask,1)&globalzombiespots(:,end-2)<size(rotatedmask,2),:);
ids=sub2ind([size(rotatedmask,1),size(rotatedmask,2)],croppedspots(:,end-3),croppedspots(:,end-2));
croppedspots(:,11)=ids;
zombieimage=zeros(size(rotatedmask));
for i=1:size(ids,1)
    zombieimage(ids(i))=croppedspots(i,4);
end
imwrite(reduceImage(reduceImage(reduceImage(imdilate(zombieimage,strel('disk',30))))),fullfile(matlabzombie,'mappedzomiespotsR0.tif'))
size(zombieimage)
size(rotatedmask)
globalzombiespots=croppedspots;

% rotatedmask=rotatedmask(1:min(size(rotatedmask,1),size(zombieimage,1)),1:min(size(rotatedmask,2),size(zombieimage,2)));
% zombieimage=zombieimage(1:min(size(rotatedmask,1),size(zombieimage,1)),1:min(size(rotatedmask,2),size(zombieimage,2)));

% cellstats=regionprops(rotatedmask>0,rotatedmask+0,'Area','PixelIdxList','MeanIntensity');
for i=1:numberofcells
    mypixels=find(rotatedmask==i);
    positions=ismember(globalzombiespots(:,11),mypixels);
    globalzombiespots(positions,12)=i;%add a cell ID to each transcript
end
%% new code
zombietableR0=zeros(numberofcells,max(globalzombiespots(:,4)));
for i=1:size(globalzombiespots,1)
    if (globalzombiespots(i,12)>0)
        zombietableR0(globalzombiespots(i,12),globalzombiespots(i,4))=1;
    end
end
%% end new code

writematrix(zombietableR0,fullfile(matlabzombie,'zombiecounttableNEWsample4.csv'))
%%