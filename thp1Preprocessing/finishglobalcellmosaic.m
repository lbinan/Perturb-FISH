seedsObject = matfile(fullfile(matlabtempfiles,'seedsmosaicAGAIN.mat'));%load seeds for watershed
seedsObject.Properties.Writable = true;
seeds = double(seedsObject.seedsmosaic(64000:114200,81000:111500)); %watershed reguide a lot of memory. image is segmented one half at a time

blurredObject = matfile(fullfile(matlabtempfiles,'blurred.mat'));%load blurred image (generated from transcript density. defins basins for watershed
blurredObject.Properties.Writable = true;
blurred = blurredObject.blurred(64000:114200,81000:111500); 

blurred=imcomplement(blurred)+imcomplement(30*imgaussfilt(seeds,30));%blurred seeds are added back to the basins, to get smoother segmentation around seeds
gmag2 = imimposemin(blurred, seeds);
clear blurred seeds
save(fullfile(matlabtempfiles,'gmag2part2.mat'),'gmag2','-v7.3')

mymosaic=watershed(gmag2);
clear gmag2
save(fullfile(matlabtempfiles,'watershedmymosaicpart2temp.mat'),'mymosaic','-v7.3')
mymosaic(mymosaic>0)=1;
cellsObject = matfile(fullfile(matlabtempfiles,'cellsmosaic.mat'));
cellsObject.Properties.Writable = true;
cellsmosaic = cellsObject.cellsmosaic(64000:114200,81000:111500);
mymosaic(cellsmosaic==0)=0;
clear cellsmosaic
mymosaic=imbinarize(mymosaic);
mymosaicpart2=bwareaopen(mymosaic,3000);
save(fullfile(matlabtempfiles,'watershedmymosaicpart2.mat'),'mymosaicpart2','-v7.3')
%%
seedsObject = matfile(fullfile(matlabtempfiles,'seedsmosaicAGAINcleanedup.mat'));
seedsObject.Properties.Writable = true;
seeds = double(seedsObject.seedsmosaic(65800:116700,1:47000)); 

blurredObject = matfile(fullfile(matlabtempfiles,'blurred.mat'));
blurredObject.Properties.Writable = true;
blurred = blurredObject.blurred(65800:116700,1:47000); 

blurred=imcomplement(blurred)+imcomplement(30*imgaussfilt(seeds,30));
gmag2 = imimposemin(blurred, seeds);
clear blurred seeds
save(fullfile(matlabtempfiles,'gmag2part1.mat'),'gmag2','-v7.3')

mymosaic=watershed(gmag2);
clear gmag2
save(fullfile(matlabtempfiles,'watershedmymosaicpart1temp.mat'),'mymosaic','-v7.3')
mymosaic(mymosaic>0)=1;
cellsObject = matfile(fullfile(matlabtempfiles,'cellsmosaic.mat'));
cellsObject.Properties.Writable = true;
cellsmosaic = cellsObject.cellsmosaic(65800:116700,1:47000);
mymosaic(cellsmosaic==0)=0;
clear cellsmosaic
mymosaic=imbinarize(mymosaic);
mymosaicpart1=bwareaopen(mymosaic,3000);
save(fullfile(matlabtempfiles,'watershedmymosaicpart1.mat'),'mymosaicpart1','-v7.3')
%% ok
load(fullfile(matlabtempfiles,'seeds.mat'))
merlinmask=zeros(size(seeds));
load(fullfile(matlabtempfiles,'globalmerfishspots.mat'))

merlinmask(globalmerfishspots(:,10))=1;
save(fullfile(matlabtempfiles,'spotsimage.mat'),'merlinmask','-v7.3')

merlinmask=imdilate(imbinarize(merlinmask),strel('disk',14));
finalmask=zeros(size(seeds));
finalmask(65800:116700,1:47000)=mymosaicpart1;
finalmask(64000:114200,81000:111500)=mymosaicpart2;%pull back both ahalves together
finalmask(~merlinmask)=0;
finalmask=finalmask>0;
save(fullfile(matlabtempfiles,'finalmask.mat'),'finalmask','-v7.3')

clear globalmerfishspots merlinmask
cellsstats=regionprops(finalmask,'Centroid','PixelIdxList','Eccentricity','Area','Circularity','Extent');
statswithimages=regionprops(finalmask,"Image");
clear mymosaic
load(fullfile(matlabtempfiles,'globalmerfishspots.mat'))
load(fullfile(matlabtempfiles,'globalzombiespots.mat'))
countable=zeros(size(cellsstats,1),8+max(globalmerfishspots(:,1)));
countableZombie=zeros(size(cellsstats,1),9+78);

disp(strcat(num2str(size(cellsstats)),' objects detected'))
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
for i=1:size(cellsstats,1)%for each cell, find transcripts and guides that are part of it
    mypixels=cellsstats(i).PixelIdxList;
    positions=ismember(globalmerfishspots(:,10),mypixels);
    positionszombie=ismember(globalzombiespots(:,9),mypixels);
    %sum(double(positions));
    globalmerfishspots(positions,11)=i;%add a cell ID to each transcript
    globalzombiespots(positionszombie,10)=i;%add a cell ID to each guide
end
writematrix(globalzombiespots,fullfile(matlabtempfiles,'globalzombiespotswithcellIDbacktomax.csv'))
writematrix(globalmerfishspots,fullfile(matlabtempfiles,'globalmerfishspotswithcellIDbacktomax.csv'))


for thiscell=1:size(cellsstats,1)%build the 2 counts matrices
    if cellsstats(thiscell).Area>0
        countable(thiscell,1)=1;%thisFOV;%#index of that cell
        countable(thiscell,2)=thiscell;%#index of that cell
        countable(thiscell,3)=cellsstats(thiscell).Area;
        countable(thiscell,4)=cellsstats(thiscell).Eccentricity;
        countable(thiscell,5)=cellsstats(thiscell).Circularity;
        countable(thiscell,6)=cellsstats(thiscell).PixelIdxList(1);
        countable(thiscell,7)=cellsstats(thiscell).Extent(1);
        mytable=globalmerfishspots(globalmerfishspots(:,11)==thiscell,:);
        cellcountable=zeros(1,130);
        if size(mytable,1)>0
            for gene=1:max(globalmerfishspots(:,1))
                cellcountable(gene)=sum(mytable(:,1)==gene);
            end
            countable(thiscell,9:end)=cellcountable;
        end
    end
    if cellsstats(thiscell).Area>0 
        countableZombie(thiscell,1)=1;%thisFOV;%#index of that cell
        countableZombie(thiscell,2)=thiscell;%#index of that cell
        countableZombie(thiscell,3)=cellsstats(thiscell).Area;
        countableZombie(thiscell,4)=cellsstats(thiscell).Eccentricity;
        countableZombie(thiscell,5)=cellsstats(thiscell).Circularity;
        countableZombie(thiscell,6)=cellsstats(thiscell).PixelIdxList(1);
        countableZombie(thiscell,7)=cellsstats(thiscell).Extent(1);
        mytable2=globalzombiespots(globalzombiespots(:,10)==thiscell,:);
        if size(mytable2,1)>0
             countableZombie(thiscell,8)=size(mytable2,1);
             for j=1:size(mytable2,1)
                 countableZombie(thiscell,8+mytable2(j,4))=1;
             end
        end
    end
end
writematrix(countableZombie,fullfile(matlabtempfiles,'zombietableonmosaicbacktomax.csv'))
writematrix(countable,fullfile(matlabtempfiles,'merfishtableonmosaicbacktomax.csv'))

for i=1:70
    name=strcat('guide',num2str(i));
    mkdir(strcat('/broad/clearylab/Users/Loic/thp1homemadezombie_1/cellimages/',name))
    myindeces=countableZombie(:,8+i)>0;
    croppedstats=statswithimages(myindeces);
    for j=1:min(sum(myindeces),40)
        imwrite(croppedstats(j).Image,strcat('/broad/clearylab/Users/Loic/thp1homemadezombie_1/cellimages/',name,'/cell',num2str(j),'.png'));
    end
end

