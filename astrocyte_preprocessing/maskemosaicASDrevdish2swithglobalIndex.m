mypath='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/data'
matlabtempfiles='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/matlabtempfiles'
nucleipath='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/nuclei'
barcodes=readtable('/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/202402142122_loicasdRevsDish2_VMSC11302/ExportBarcodes/region_0/barcodes.csv');
barcodes=table2array(barcodes(:,1:8));
minX=0;
minY=0;
maxX=0;
maxY=0;
myfovs=[min(barcodes(:,8)):max(barcodes(:,8))];
globalmerfishspots=barcodes;
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
rotatedtemplate=imrotate(thisimage,-90);%generate templates of pixel IDS to facilitate rotations later
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
reverserotatedtemplate=imrotate(thisimage,90);
for fov=1:size(myfovs,2)
    thisFOV=myfovs(1,fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
    else
        thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
    end
    fileID = fopen(thispath,'r');
    formatSpec = '%s';
    A = fscanf(fileID,formatSpec);
    position=strfind(A,'position');
    nextposition=strfind(A(position:end),',');
    X=A(position+10:position+nextposition-2);
    Y=A(position+nextposition(1,1):position+nextposition(1,2)-2);
    X=sscanf(X, '%d')
    Y=sscanf(Y, '%d')
    fclose(fileID)
    thisX=X/200
    thisY=Y/200
    minX=min(minX,thisX);
    minY=min(minY,thisY);
    maxX=max(maxX,thisX);
    maxY=max(maxY,thisY);
end
save(fullfile(matlabtempfiles,'minYRegion0.mat'),'minY')
save(fullfile(matlabtempfiles,'minXRegion0.mat'),'minX')
save(fullfile(matlabtempfiles,'maxYRegion0.mat'),'maxY')
save(fullfile(matlabtempfiles,'maxXRegion0.mat'),'maxX')
milieuX=round((abs(minX)+abs(maxX))/2)+1;
milieuY=round((abs(minY)+abs(maxY))/2)+1;
seedsmosaic=double(zeros(floor(2048+(maxY-minY)*2048),floor(2048+(maxX-minX)*2048)));
dapimosaic=seedsmosaic;
% seedsmosaic=seedsmosaic>0;
disp('seedsmosaic size')

for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_000',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_00',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_0',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    else
        thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    end
    fileID = fopen(thispath,'r');
    formatSpec = '%s';
    A = fscanf(fileID,formatSpec);
    position=strfind(A,'position');
    nextposition=strfind(A(position:end),',');
    X=A(position+10:position+nextposition-2);
    Y=A(position+nextposition(1,1):position+nextposition(1,2)-2);
    X = sscanf(X, '%d');
    Y= sscanf(Y, '%d');
    fclose(fileID)
    thisX=X/200+abs(minX)+1;
    thisY=Y/200+abs(minY)+1;
    nuclei=imbinarize(thisImage,'adaptive');
    nuclei=bwareaopen(nuclei,7000);
    nuclei=imfill(nuclei,'holes');
    nuclei=imdilate(nuclei,strel('disk',5));
    split=bwareaopen(nuclei,35000);
    keep=(nuclei-split)>0;
    split=split-bwareaopen(split,100000);
    D = bwdist(~split);
    D = -D;
    I3 = imhmin(D,2);
    disp('got to I3')
    resplitnuclei = watershed(I3);
    disp('finish watershed')
    resplitnuclei(~split)=0;
    mymask=bwareaopen(((resplitnuclei>0)+(keep>0))>0,2000);
    dapi=flipdim(mymask,2);
    dapi=imrotate(dapi,-90);
    imwrite(mymask,fullfile(nucleipath,strcat('nuclei',num2str(thisFOV),'.tif')));
    cells=regionprops(dapi,'Centroid');
    X=zeros(2048,2048);
    for i=1:size(cells,1)
        X(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
    end
    X=X>0;
%     X=imdilate(X,strel('disk',25));
    seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),X);
    dapimosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(dapimosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),dapi);
    merfishcoo=floor(globalmerfishspots(globalmerfishspots(:,8)==thisFOV,6:7));
%     merfishcoo(:,2)=2049-merfishcoo(:,2);
%     merfishcoo(:,2)=merfishcoo(:,2);
%     merfishcoo(:,1)=merfishcoo(:,1);
    merfishindeces=sub2ind([2048 2048],1+merfishcoo(:,2),1+merfishcoo(:,1));
%     rotatedindecesmerfish=rotatedtemplate(merfishindeces);
    [merfishcoo(:,1),merfishcoo(:,2)]=ind2sub([2048 2048],merfishindeces);
    merfishcoo(:,1)=floor((thisY-1)*1843.2+1+merfishcoo(:,1));
    merfishcoo(:,2)=floor((thisX-1)*1843.2+1+merfishcoo(:,2));
     globalindex=sub2ind([size(seedsmosaic,1) size(seedsmosaic,2)],1+merfishcoo(:,1),1+merfishcoo(:,2));

    %converting coordinates of RNA transcript spot from coordinates inside
    %the FOV to coordinates in the mosaic
    globalmerfishspots(globalmerfishspots(:,8)==thisFOV,end-3:end)=[merfishcoo,merfishindeces,globalindex];

end
save(fullfile(matlabtempfiles,'globalmerfishspotsregion0.mat'),'globalmerfishspots','-v7.3')
save(fullfile(matlabtempfiles,'seedsregion0.mat'),'seedsmosaic','-v7.3')
save(fullfile(matlabtempfiles,'binaryDapiregion0.mat'),'dapimosaic','-v7.3')
% Se=reduceImage(reduceImage(seedsmosaic));
% save(fullfile(matlabtempfiles,'Searea0.mat'),'Se','-v7.3')
% imwrite(reduceImage(Se),fullfile(matlabtempfiles,'reducedSearea0.tif'));
% clear Se
% SD=reduceImage(reduceImage(dapimosaic));
% save(fullfile(matlabtempfiles,'SDarea0.mat'),'SD','-v7.3')
% imwrite(reduceImage(SD),fullfile(matlabtempfiles,'reducedSDarea0.tif'));
% clear SD
% testimage=zeros(size(seedsmosaic));
% clear seedsmosaic dapimosaic
% 
% % for i=1:size(globalmerfishspots,1)
% %     testimage(floor(globalmerfishspots(i,end-2)),floor(globalmerfishspots(i,end-3)))=1;
% % end
% % testimage=imdilate(testimage,strel('disk',15));
% % testimage=reduceImage(reduceImage(reduceImage(testimage)));
% % imwrite(testimage,fullfile(matlabtempfiles,'spotsarea0.tif'));
% % testimage=zeros(size(testimage));
% for i=1:size(globalmerfishspots,1)
%     testimage(floor(globalmerfishspots(i,end-3)),floor(globalmerfishspots(i,end-2)))=1;
% end
% testimage=imdilate(testimage,strel('disk',15));
% testimage=reduceImage(reduceImage(reduceImage(testimage)));
% imwrite(testimage,fullfile(matlabtempfiles,'spotsbisarea0.tif'));



%%
barcodes=readtable('/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/202402142122_loicasdRevsDish2_VMSC11302/ExportBarcodes/region_1/barcodes.csv');
barcodes=table2array(barcodes(:,1:8));
minX=0;
minY=0;
maxX=0;
maxY=0;
myfovs=[min(barcodes(:,8)):max(barcodes(:,8))];
globalmerfishspots=barcodes;
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
rotatedtemplate=imrotate(thisimage,-90);%generate templates of pixel IDS to facilitate rotations later
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
reverserotatedtemplate=imrotate(thisimage,90);
for fov=1:size(myfovs,2)
    thisFOV=myfovs(1,fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
    else
        thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
    end
    fileID = fopen(thispath,'r');
    formatSpec = '%s';
    A = fscanf(fileID,formatSpec);
    position=strfind(A,'position');
    nextposition=strfind(A(position:end),',');
    X=A(position+10:position+nextposition-2);
    Y=A(position+nextposition(1,1):position+nextposition(1,2)-2);
    X=sscanf(X, '%d')
    Y=sscanf(Y, '%d')
    fclose(fileID)
    thisX=X/200
    thisY=Y/200
    minX=min(minX,thisX);
    minY=min(minY,thisY);
    maxX=max(maxX,thisX);
    maxY=max(maxY,thisY);
end
save(fullfile(matlabtempfiles,'minYRegion1.mat'),'minY')
save(fullfile(matlabtempfiles,'minXRegion1.mat'),'minX')
save(fullfile(matlabtempfiles,'maxYRegion1.mat'),'maxY')
save(fullfile(matlabtempfiles,'maxXRegion1.mat'),'maxX')
milieuX=round((abs(minX)+abs(maxX))/2)+1;
milieuY=round((abs(minY)+abs(maxY))/2)+1;
seedsmosaic=double(zeros(floor(2048+(maxY-minY)*2048),floor(2048+(maxX-minX)*2048)));
dapimosaic=seedsmosaic;
% seedsmosaic=seedsmosaic>0;
disp('seedsmosaic size')
myfovs=[min(barcodes(:,8)):max(barcodes(:,8))];

for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_000',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_00',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_0',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    else
        thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_',num2str(thisFOV),'.dax')));
        thisImage=double(thisImage(:,:,24))/750;
    end
    fileID = fopen(thispath,'r');
    formatSpec = '%s';
    A = fscanf(fileID,formatSpec);
    position=strfind(A,'position');
    nextposition=strfind(A(position:end),',');
    X=A(position+10:position+nextposition-2);
    Y=A(position+nextposition(1,1):position+nextposition(1,2)-2);
    X = sscanf(X, '%d');
    Y= sscanf(Y, '%d');
    fclose(fileID)
    thisX=X/200+abs(minX)+1;
    thisY=Y/200+abs(minY)+1;
    nuclei=imbinarize(thisImage,'adaptive');
    nuclei=bwareaopen(nuclei,7000);
    nuclei=imfill(nuclei,'holes');
    nuclei=imdilate(nuclei,strel('disk',5));
    split=bwareaopen(nuclei,35000);
    keep=(nuclei-split)>0;
    split=split-bwareaopen(split,100000);
    D = bwdist(~split);
    D = -D;
    I3 = imhmin(D,2);
    disp('got to I3')
    resplitnuclei = watershed(I3);
    disp('finish watershed')
    resplitnuclei(~split)=0;
    mymask=bwareaopen(((resplitnuclei>0)+(keep>0))>0,2000);
    dapi=flipdim(mymask,2);
    dapi=imrotate(dapi,-90);
    imwrite(mymask,fullfile(nucleipath,strcat('nuclei',num2str(thisFOV),'.tif')));
    cells=regionprops(dapi,'Centroid');
    X=zeros(2048,2048);
    for i=1:size(cells,1)
        X(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
    end
    X=X>0;
%     X=imdilate(X,strel('disk',25));
    seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),X);
    dapimosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(dapimosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),dapi);
    merfishcoo=floor(globalmerfishspots(globalmerfishspots(:,8)==thisFOV,6:7));
%     merfishcoo(:,2)=2049-merfishcoo(:,2);
%     merfishcoo(:,2)=merfishcoo(:,2);
%     merfishcoo(:,1)=merfishcoo(:,1);
    merfishindeces=sub2ind([2048 2048],1+merfishcoo(:,2),1+merfishcoo(:,1));
%     rotatedindecesmerfish=rotatedtemplate(merfishindeces);
    [merfishcoo(:,1),merfishcoo(:,2)]=ind2sub([2048 2048],merfishindeces);
    merfishcoo(:,1)=floor((thisY-1)*1843.2+1+merfishcoo(:,1));
    merfishcoo(:,2)=floor((thisX-1)*1843.2+1+merfishcoo(:,2));
    globalindex=sub2ind([size(seedsmosaic,1) size(seedsmosaic,2)],1+merfishcoo(:,1),1+merfishcoo(:,2));
    %converting coordinates of RNA transcript spot from coordinates inside
    %the FOV to coordinates in the mosaic
    globalmerfishspots(globalmerfishspots(:,8)==thisFOV,end-3:end)=[merfishcoo,merfishindeces,globalindex];

end
save(fullfile(matlabtempfiles,'globalmerfishspotsregion1.mat'),'globalmerfishspots','-v7.3')
save(fullfile(matlabtempfiles,'seedsregion1.mat'),'seedsmosaic','-v7.3')
save(fullfile(matlabtempfiles,'binaryDapiregion1.mat'),'dapimosaic','-v7.3')

% Se=reduceImage(reduceImage(seedsmosaic));
% save(fullfile(matlabtempfiles,'Searea1.mat'),'Se','-v7.3')
% imwrite(reduceImage(Se),fullfile(matlabtempfiles,'reducedSearea1.tif'));
% clear Se
% SD=reduceImage(reduceImage(dapimosaic));
% save(fullfile(matlabtempfiles,'SDarea0.mat'),'SD','-v7.3')
% imwrite(reduceImage(SD),fullfile(matlabtempfiles,'reducedSDarea1.tif'));
% clear SD
% testimage=zeros(size(seedsmosaic));
% clear seedsmosaic dapimosaic

% for i=1:size(globalmerfishspots,1)
%     testimage(floor(globalmerfishspots(i,end-2)),floor(globalmerfishspots(i,end-3)))=1;
% end
% testimage=imdilate(testimage,strel('disk',15));
% testimage=reduceImage(reduceImage(reduceImage(testimage)));
% imwrite(testimage,fullfile(matlabtempfiles,'spotsarea1.tif'));
% testimage=zeros(size(testimage));
% for i=1:size(globalmerfishspots,1)
%     testimage(floor(globalmerfishspots(i,end-3)),floor(globalmerfishspots(i,end-2)))=1;
% end
% testimage=imdilate(testimage,strel('disk',15));
% testimage=reduceImage(reduceImage(reduceImage(testimage)));
% imwrite(testimage,fullfile(matlabtempfiles,'spotsbisarea1.tif'));

clear
maskemosaicASDrevdish3swithglobalIndex
