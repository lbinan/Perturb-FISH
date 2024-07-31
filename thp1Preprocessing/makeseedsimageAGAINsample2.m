localpath='\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\Data'
globalpath='/broad/clearylab/Users/Loic/thp1homedish2_zombie/Data'
savepath='/broad/clearylab/Users/Loic/thp1homedish2_zombie/preprocessed/'
matlabtempfiles='/broad/clearylab/Users/Loic/thp1homedish2_zombie/matlabtempfilesrealigned/'
maxprojpath='/broad/clearylab/Users/Loic/thp1homedish2_zombie/maxproj/'
filteredpath='/broad/clearylab/Users/Loic/thp1homedish2_zombie/filteredimages/'
merlinpath='/broad/clearylab/Users/Loic/thp1homedish2_merfish/merlin/merlin/Data/GenerateMosaic/images'
maskedagain='/broad/clearylab/Users/Loic/thp1homedish2_zombie/masked'
filtered='/broad/clearylab/Users/Loic/thp1homedish2_zombie/filtered'
mypath='/broad/clearylab/Users/Loic/thp1homedish2_zombie/Data'
localpath=globalpath;
minX=0;
minY=0;
maxX=0;
maxY=0;
myfovs=[0:1252];%change for current dataset
% myfovs=[0:150];
% myfovs=[0:2];
%load decoded transcripts, from merlin
merfishbarcodes=table2array(readtable(fullfile(merlinpath,'barcodes.csv')));
%load QC filtered guide spots
zombiebarcodes=table2array(readtable(fullfile(matlabtempfiles,'filtereddecodedGuidesmanualmanualadjust.csv')));
globalmerfishspots=merfishbarcodes;
globalzombiespots=zombiebarcodes;
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalmerfishspots(:,1)=globalmerfishspots(:,1)+1;
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
rotatedtemplate=imrotate(thisimage,-90);

thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
reverserotatedtemplate=imrotate(thisimage,90);%get image coordinates to facilitate rotations later

% for fov=1:size(myfovs,2)%if varialbles minX min Y etc don't exist, run
% this first to find mosaic size in number of images and its center
% thisFOV=myfovs(fov);
% disp(strcat('doing FOV ', num2str(thisFOV)))
% if thisFOV<10
%     thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_1.xml')));
% elseif thisFOV<100
%     thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_1.xml')));
% elseif thisFOV<1000
%     thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_1.xml')));
% else
%     thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_1.xml')));
% end
% thisposition=thisinfo.settings.acquisition.stage_position.Text;
% C = textscan(thisposition, '%f,%f');
% thisX=C{1}/200;
% thisY=C{2}/200;
% minX=min(minX,thisX);
% minY=min(minY,thisY);
% maxX=max(maxX,thisX);
% maxY=max(maxY,thisY);
% end
% mymosaic=imread(fullfile(matlabtempfiles,'testmosaic4.dax'));
% save(fullfile(matlabtempfiles,'minY.mat'),'minY')
% save(fullfile(matlabtempfiles,'minX.mat'),'minX')
% save(fullfile(matlabtempfiles,'maxY.mat'),'maxY')
% save(fullfile(matlabtempfiles,'maxX.mat'),'maxX')
load(fullfile(matlabtempfiles,'minY.mat'))
load(fullfile(matlabtempfiles,'minX.mat'))
load(fullfile(matlabtempfiles,'maxY.mat'))
load(fullfile(matlabtempfiles,'maxX.mat'))
milieuX=round((abs(minX)+abs(maxX))/2)+1;
milieuY=round((abs(minY)+abs(maxY))/2)+1;
seedsmosaic=zeros(floor(2048+(maxY-minY)*2048),floor(2048+(maxX-minX)*2048));
seedsmosaic=seedsmosaic>0;
% grid=seedsmosaic;
% gridimage=zeros(2048,2048);
% gridimage(1:100,:)=1;
% gridimage(end-100:end,:)=1;
% gridimage(:,1:100)=1;
% gridimage(:,end-100:end)=1;
% gridimage=gridimage>0;
%% 
for fov=1:size(myfovs,2)%this generates a mosaic of the seeds and converts coordinates of rna and guide spots from
    %FOV cooredinates to global mosaic coordinates
    thisFOV=myfovs(fov);
    if rem(thisFOV,100)==0
        disp(strcat('current FOV is ',num2str(thisFOV)))
    end
    thisfilesname=fullfile(maskedagain,strcat('cellmask', num2str(thisFOV),'.tif'));
    seedimage=imread(thisfilesname);
    if thisFOV<10
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_1.xml')));
    elseif thisFOV<100
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_1.xml')));
    elseif thisFOV<1000
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_1.xml')));
    else
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_1.xml')));
    end
    thisposition=thisinfo.settings.acquisition.stage_position.Text;
    %%
    % thisimage=flipdim(thisimage,1);
    % thisimage=flipdim(thisimage,2);
    % thisimage=imrotate(thisimage,-90);
    seedimage=flipdim(seedimage,1);%coorect camera orientation
    seedimage=flipdim(seedimage,2);
    seedimage=imrotate(seedimage,-90);
    seedimage=bwareaopen(seedimage,800);
    cells=regionprops(seedimage,'Centroid');
    X=zeros(2048,2048);
    for i=1:size(cells,1)
        X(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
    end
    X=X>0;
    X=imdilate(X,strel('disk',25));
    thisposition=thisinfo.settings.acquisition.stage_position.Text;
    C = textscan(thisposition, '%f,%f');
    thisX=C{1}/200+abs(minX)+1;
    thisY=C{2}/200+abs(minY)+1;
    seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),X);
    % grid(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(grid(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),gridimage);
    %below, we convert coorediantes of spots for the orientation of the
    %camera, like we did twith the actual image
    zombiecoo=floor(globalzombiespots(globalzombiespots(:,3)==(thisFOV),1:2));
    zombiecoo(:,1)=2049-zombiecoo(:,1);
    zombiecoo(:,2)=2049-zombiecoo(:,2);
    zombieindeces=sub2ind([2048 2048],zombiecoo(:,1),zombiecoo(:,2));
    rotatedindeces=reverserotatedtemplate(zombieindeces);
    [zombiecoo(:,1),zombiecoo(:,2)]=ind2sub([2048 2048],rotatedindeces);
    zombiecoo(:,1)=floor((thisY-1)*1843.2+1+zombiecoo(:,1));
    zombiecoo(:,2)=floor((thisX-1)*1843.2+1+zombiecoo(:,2));
    globalzombiespots(globalzombiespots(:,3)==(thisFOV),end-2:end)=[zombiecoo,rotatedindeces];
    merfishcoo=floor(globalmerfishspots(globalmerfishspots(:,7)==thisFOV,5:6));
    merfishcoo(:,2)=2049-merfishcoo(:,2);
    merfishcoo(:,2)=merfishcoo(:,2);
    merfishcoo(:,1)=merfishcoo(:,1);
    merfishindeces=sub2ind([2048 2048],merfishcoo(:,1),merfishcoo(:,2));
    rotatedindecesmerfish=rotatedtemplate(merfishindeces);
    [merfishcoo(:,1),merfishcoo(:,2)]=ind2sub([2048 2048],rotatedindecesmerfish);
    merfishcoo(:,1)=floor((thisY-1)*1843.2+1+merfishcoo(:,1));
    merfishcoo(:,2)=floor((thisX-1)*1843.2+1+merfishcoo(:,2));
    globalmerfishspots(globalmerfishspots(:,7)==thisFOV,end-2:end)=[merfishcoo,rotatedindecesmerfish];
end
%%

save(fullfile(matlabtempfiles,'seedsmosaicAGAIN.mat'),'seedsmosaic','-v7.3')
save(fullfile(matlabtempfiles,'globalmerfishspots.mat'),'globalmerfishspots','-v7.3')
save(fullfile(matlabtempfiles,'globalzombiespots.mat'),'globalzombiespots','-v7.3')

% clear globalmerfishspots globalzombiespots
% cells=regionprops(seedsmosaic,'Centroid');
% X=zeros(size(seedsmosaic));
% for i=1:size(cells,1)
%     X(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
% end
% X=X>0;
% seedsmosaic=imdilate(X,strel('disk',25));
% clear X
% imwrite(seedsmosaic(100000:110000,40000:50000),fullfile(matlabtempfiles,'sampleseedsmosaic.tif'))
% save(fullfile(matlabtempfiles,'seedsmosaicAGAINcleanedup.mat'),'seedsmosaic','-v7.3')
% 
merfishspots=zeros(size(seedsmosaic));%generate a mosaic with all decoded RNAs show as 1
clear seedsmosaic
load(fullfile(matlabtempfiles,'globalmerfishspots.mat'))
for i=1:size(globalmerfishspots,1)
    merfishspots(globalmerfishspots(i,end-2),globalmerfishspots(i,end-1))=1;
end
% imwrite(merfishspots(100000:110000,40000:50000),fullfile(matlabtempfiles,'samplemerfishspots.tif'))
% save(fullfile(matlabtempfiles,'merfishspotsimage.mat'),'merfishspots','-v7.3')
blurred=imgaussfilt(80*merfishspots,35);%generate a transcript density map. this is used as basins during watershed segmentation
save(fullfile(matlabtempfiles,'blurred.mat'),'blurred','-v7.3')
globalmerfishspots(:,10)=sub2ind([size(merfishspots,1) size(merfishspots,2)],globalmerfishspots(:,8),globalmerfishspots(:,9));
% merfishspots=zeros(size(merfishspots));
% merfishspots(globalmerfishspots(:,10))=1;
save(fullfile(matlabtempfiles,'globalmerfishspots.mat'),'globalmerfishspots','-v7.3')
% load(fullfile(matlabtempfiles,'globalzombiespots.mat'))
load('/broad/clearylab/Users/Loic/thp1homedish2_zombie/matlabtempfiles/finalmask.mat')
globalzombiespots(:,9)=sub2ind([size(finalmask,1) size(finalmask,2)],globalzombiespots(:,7),globalzombiespots(:,8));
save(fullfile(matlabtempfiles,'globalzombiespots.mat'),'globalzombiespots','-v7.3')
save(fullfile(matlabtempfiles,'blurred.mat'),'blurred','-v7.3')

finishglobalcellmosaic