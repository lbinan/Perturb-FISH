mypath='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/data'
matlabtempfiles='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles'
nucleipath='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/nuclei'
load(fullfile(matlabtempfiles,'binaryDapiregion0optiNoOverlapAGAINSmallRotation2.mat'));%loads dapimosaic
dilated=imdilate(dapimosaic>0,strel('disk',20));%dilate to include cytoplasm
disp('dilated')
load(fullfile(matlabtempfiles,'binaryDapiregion0optiNoOverlapAGAINSmallRotation2Erodednocontour.mat'));%loads dapimosaic
bigseeds=imerode(dapimosaic>0,strel('disk',10));%shrinks the nuceli to serve as seeds, also seperates nuclei that got merged at the stiching interface
clear dapimosaic
load(fullfile(matlabtempfiles,'globalmerfishspotsregion0optiSmallRotation2.mat'))%loads globalmerfishspots


D = bwdist(dilated);
I3 = imhmin(D,5);
gmag2 = imimposemin(I3 , bigseeds);
mymask = watershed(gmag2);
mymask(~dilated)=0;
mymask=mymask>0;
mymask=bwareaopen(mymask,2000);

figure, imshow(mymask);
save(fullfile(matlabtempfiles,'merfishMask.mat'),'mymask','-v7.3')
disp('starting')
ids=sub2ind([size(mymask,1),size(mymask,2)],globalmerfishspots(:,end-3),globalmerfishspots(:,end-2));
globalmerfishspots(:,13)=ids;
cellstats=regionprops(mymask,'Area','Centroid','PixelIdxList');
for i=1:size(cellstats,1)
    mypixels=cellstats(i).PixelIdxList;
    positions=ismember(globalmerfishspots(:,13),mypixels);
    globalmerfishspots(positions,14)=i;%add a cell ID to each transcript
end
countable=zeros(size(cellstats,1),4+max(globalmerfishspots(:,2)));
for thiscell=1:size(cellstats,1)%build the count table
    if rem(thiscell,200)==0
        disp(strcat('current cell is ',num2str(thiscell), 'out of ', num2str(size(countable,1))))%progress report
    end
    if cellstats(thiscell).Area>0 
        countable(thiscell,1)=thiscell;%#index of that cell
        countable(thiscell,2)=cellstats(thiscell).Area;%#area of that cell
        countable(thiscell,3)=sub2ind([size(mymask,1),size(mymask,2)],cellstats(thiscell).Centroid(2),cellstats(thiscell).Centroid(1));%#position of that cell
        mytable2=globalmerfishspots(globalmerfishspots(:,14)==thiscell,:);
        for gene=1:max(globalmerfishspots(:,2))+1
            cellcountable(1,gene)=sum(mytable2(:,2)==gene-1);
        end
        countable(thiscell,4:end)=cellcountable;
    end
end
writematrix(countable,fullfile(matlabtempfiles,'merfishcounttable.csv'))

clear
tumorsZombiecounttable