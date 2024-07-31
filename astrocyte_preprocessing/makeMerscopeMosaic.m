mypath='/broad/clearylab/Users/Loic/revisionsdish1merfish/Data'
matlabtempfiles='/broad/clearylab/Users/Loic/revisionsdish1merfish/matlabtempfiles'
barcodes=readtable('/broad/clearylab/Users/Loic/revisionsdish1merfish/ExportBarcodes/region_0/barcodes.csv');
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
% seedsmosaic=seedsmosaic>0;
disp('seedsmosaic size')
for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_000',num2str(thisFOV),'.dax')));
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_00',num2str(thisFOV),'.dax')));
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_0',num2str(thisFOV),'.dax')));
    else
        thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_',num2str(thisFOV),'.dax')));
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
    dapi=double(thisImage(:,:,24));
%     dapi=flipdim(dapi,1);%compensate for camera orientation on the scope
    dapi=flipdim(dapi,2);
    dapi=imrotate(dapi,-90);
    seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),dapi);
    merfishcoo=floor(globalmerfishspots(globalmerfishspots(:,8)==thisFOV,6:7));
    merfishcoo(:,2)=2049-merfishcoo(:,2);
    merfishcoo(:,2)=merfishcoo(:,2);
    merfishcoo(:,1)=merfishcoo(:,1);
    merfishindeces=sub2ind([2048 2048],merfishcoo(:,1),merfishcoo(:,2));
    rotatedindecesmerfish=rotatedtemplate(merfishindeces);
    [merfishcoo(:,1),merfishcoo(:,2)]=ind2sub([2048 2048],rotatedindecesmerfish);
    merfishcoo(:,1)=floor((thisY-1)*1843.2+1+merfishcoo(:,1));
    merfishcoo(:,2)=floor((thisX-1)*1843.2+1+merfishcoo(:,2));
    %converting coordinates of RNA transcript spot from coordinates inside
    %the FOV to coordinates in the mosaic
    globalmerfishspots(globalmerfishspots(:,8)==thisFOV,end-2:end)=[merfishcoo,rotatedindecesmerfish];
    
end
save(fullfile(matlabtempfiles,'dapiregion0.mat'),'seedsmosaic','-v7.3')
save(fullfile(matlabtempfiles,'globalmerfishspotsregion0.mat'),'globalmerfishspots','-v7.3')

   
barcodes=readtable('/broad/clearylab/Users/Loic/revisionsdish1merfish/ExportBarcodes/region_1/barcodes.csv');
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
    thisFOV=myfovs(fov);
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
    X = sscanf(X, '%d');
    Y= sscanf(Y, '%d');
    fclose(fileID)
    thisX=X/200;
    thisY=Y/200;
    minX=min(minX,thisX);
    minY=min(minY,thisY);
    maxX=max(maxX,thisX);
    maxY=max(maxY,thisY);
end
save(fullfile(matlabtempfiles,'minYregion1.mat'),'minY')
save(fullfile(matlabtempfiles,'minXregion1.mat'),'minX')
save(fullfile(matlabtempfiles,'maxYregion1.mat'),'maxY')
save(fullfile(matlabtempfiles,'maxXregion1.mat'),'maxX')
milieuX=round((abs(minX)+abs(maxX))/2)+1;
milieuY=round((abs(minY)+abs(maxY))/2)+1;
seedsmosaic=double(zeros(floor(2048+(maxY-minY)*2048),floor(2048+(maxX-minX)*2048)));
% seedsmosaic=seedsmosaic>0;
disp('seedsmosaic size')
for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_000',num2str(thisFOV),'.dax')));
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_00',num2str(thisFOV),'.dax')));
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_0',num2str(thisFOV),'.dax')));
    else
        thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_prestain_',num2str(thisFOV),'.dax')));
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
    dapi=double(thisImage(:,:,24));
%     dapi=flipdim(dapi,1);%compensate for camera orientation on the scope
    dapi=flipdim(dapi,2);
    dapi=imrotate(dapi,-90);
    seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048)=max(seedsmosaic(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+2048,floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+2048),dapi);
    merfishcoo=floor(globalmerfishspots(globalmerfishspots(:,8)==thisFOV,6:7));
    merfishcoo(:,2)=2049-merfishcoo(:,2);
    merfishcoo(:,2)=merfishcoo(:,2);
    merfishcoo(:,1)=merfishcoo(:,1);
    merfishindeces=sub2ind([2048 2048],merfishcoo(:,1),merfishcoo(:,2));
    rotatedindecesmerfish=rotatedtemplate(merfishindeces);
    [merfishcoo(:,1),merfishcoo(:,2)]=ind2sub([2048 2048],rotatedindecesmerfish);
    merfishcoo(:,1)=floor((thisY-1)*1843.2+1+merfishcoo(:,1));
    merfishcoo(:,2)=floor((thisX-1)*1843.2+1+merfishcoo(:,2));
    %converting coordinates of RNA transcript spot from coordinates inside
    %the FOV to coordinates in the mosaic
    globalmerfishspots(globalmerfishspots(:,8)==thisFOV,end-2:end)=[merfishcoo,rotatedindecesmerfish];
    
end
save(fullfile(matlabtempfiles,'dapiregion1.mat'),'seedsmosaic','-v7.3')
save(fullfile(matlabtempfiles,'globalmerfishspotsregion1.mat'),'globalmerfishspots','-v7.3')
shrunkregion1=reduceImage(seedsmosaic);
shrunkregion1=reduceImage(shrunkregion1);
save(fullfile(matlabtempfiles,'dapiregion1shrunk.mat'),'shrunkregion1','-v7.3')
