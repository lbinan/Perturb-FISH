%need to finish after data transfer
mypath='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/data'
matlabtempfiles='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles'
nucleipath='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/nuclei'
barcodes=readtable('/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/barcodes.csv');
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
contour=zeros(2048,2048);
contour(1:98,:)=1;
contour(:,1:98)=1;
contour(1949:2048,:)=1;
contour(:,1949:2048)=1;
contour=contour>0;
for i=1:2048*2048
    thisimage(i)=i;
end
rotatedtemplate=imrotate(thisimage,-91);%generate templates of pixel IDS to facilitate rotations later
% thisimage=zeros(2048,2048);
% for i=1:2048*2048
%     thisimage(i)=i;
% end
% reverserotatedtemplate=imrotate(thisimage,90);
% % % for fov=1:size(myfovs,2)
% % %     thisFOV=myfovs(1,fov);
% % %     disp(strcat('doing FOV ', num2str(thisFOV)))
% % %     if thisFOV<10
% % %         thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
% % %     elseif thisFOV<100
% % %         thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
% % %     elseif thisFOV<1000
% % %         thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
% % %     else
% % %         thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
% % %     end
% % %     fileID = fopen(thispath,'r');
% % %     formatSpec = '%s';
% % %     A = fscanf(fileID,formatSpec);
% % %     position=strfind(A,'position');
% % %     nextposition=strfind(A(position:end),',');
% % %     X=A(position+10:position+nextposition-2);
% % %     Y=A(position+nextposition(1,1):position+nextposition(1,2)-2);
% % %     X=sscanf(X, '%d')
% % %     Y=sscanf(Y, '%d')
% % %     fclose(fileID)
% % %     thisX=X/200
% % %     thisY=Y/200
% % %     minX=min(minX,thisX);
% % %     minY=min(minY,thisY);
% % %     maxX=max(maxX,thisX);
% % %     maxY=max(maxY,thisY);
% % % end
% % % save(fullfile(matlabtempfiles,'minYRegion0.mat'),'minY')
% % % save(fullfile(matlabtempfiles,'minXRegion0.mat'),'minX')
% % % save(fullfile(matlabtempfiles,'maxYRegion0.mat'),'maxY')
% % % save(fullfile(matlabtempfiles,'maxXRegion0.mat'),'maxX')

load(fullfile(matlabtempfiles,'minYRegion0.mat'))
load(fullfile(matlabtempfiles,'minXRegion0.mat'))
load(fullfile(matlabtempfiles,'maxYRegion0.mat'))
load(fullfile(matlabtempfiles,'maxXRegion0.mat'))

milieuX=round((abs(minX)+abs(maxX))/2)+1;
milieuY=round((abs(minY)+abs(maxY))/2)+1;
registrationimage=double(zeros(floor(2048+(maxY-minY)*2048),floor(2048+(maxX-minX)*2048)));
% dapimosaic=seedsmosaic;
% seedsmosaic=seedsmosaic>0;
disp('seedsmosaic size')

for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_0_000',num2str(thisFOV),'.dax')));
         thisImage=double(thisImage(:,:,15))/3000;
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_0_00',num2str(thisFOV),'.dax')));
         thisImage=double(thisImage(:,:,15))/3000;
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_0_0',num2str(thisFOV),'.dax')));
         thisImage=double(thisImage(:,:,15))/3000;
    else
        thispath=(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.inf')));
        thisImage=myOwnReadDax(fullfile(mypath,strcat('stack_0_',num2str(thisFOV),'.dax')));
         thisImage=double(thisImage(:,:,15))/3000;
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
    
    
    dapi=flipdim(thisImage,2);
    dapi=imrotate(dapi,-91);
    
    
    
%     X=imdilate(X,strel('disk',25));
        registrationimage(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+size(rotatedtemplate,1),floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+size(rotatedtemplate,1))=max(registrationimage(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+size(rotatedtemplate,1),floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+size(rotatedtemplate,1)),dapi);
    %     merfishcoo(:,2)=2049-merfishcoo(:,2);
%     merfishcoo(:,2)=merfishcoo(:,2);
%     merfishcoo(:,1)=merfishcoo(:,1);
%     merfishindeces=sub2ind([2048 2048],1+merfishcoo(:,2),1+merfishcoo(:,1));
%     rotatedindecesmerfish=rotatedtemplate(merfishindeces);
%     [merfishcoo(:,1),merfishcoo(:,2)]=ind2sub([2048 2048],merfishindeces);
%     merfishcoo(:,1)=floor((thisY-1)*1843.2+1+merfishcoo(:,1));
%     merfishcoo(:,2)=floor((thisX-1)*1843.2+1+merfishcoo(:,2));
%      globalindex=sub2ind([size(seedsmosaic,1) size(seedsmosaic,2)],1+merfishcoo(:,1),1+merfishcoo(:,2));

    %converting coordinates of RNA transcript spot from coordinates inside
    %the FOV to coordinates in the mosaic
%     globalmerfishspots(globalmerfishspots(:,8)==thisFOV,end-3:end)=[merfishcoo,merfishindeces,globalindex];

end
% save(fullfile(matlabtempfiles,'globalmerfishspotsregion0optiSmallRotation2.mat'),'globalmerfishspots','-v7.3')
save(fullfile(matlabtempfiles,'registrationimageMERFISH.mat'),'registrationimage','-v7.3')
% save(fullfile(matlabtempfiles,'binaryDapiregion0optiNoOverlapAGAINSmallRotation2.mat'),'dapimosaic','-v7.3')
Se=reduceImage(reduceImage(registrationimage));
save(fullfile(matlabtempfiles,'Smallregistrationimage.mat'),'Se','-v7.3')
