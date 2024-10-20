clear

mypath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/Data'
filtered='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/filtered'%output path
matlabtempfiles='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles'
nucleipath2='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/nuclei2'
nucleipathRotated='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/nucleiRotated'
barcodes=readtable('/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles/round1/manualfilteredGREAT.csv');
barcodes=table2array(barcodes);
minX=0;
minY=0;
maxX=0;
maxY=0;
myfovs=[0:1169];
globalzombiespots=barcodes;
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
rotatedtemplate=imrotate(thisimage,-87.7);%generate templates of pixel IDS to facilitate rotations later
% thisimage=zeros(2048,2048);
% for i=1:2048*2048
%     thisimage(i)=i;
% end
% reverserotatedtemplate=imrotate(thisimage,-2.3);
% for fov=1:size(myfovs,2)
%     thisFOV=myfovs(1,fov);
%     disp(strcat('doing FOV ', num2str(thisFOV)))
%     if thisFOV<10
%         thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_000',num2str(thisFOV),'.xml')));
%     elseif thisFOV<100
%         thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_00',num2str(thisFOV),'.xml')));
%     elseif thisFOV<1000
%         thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_0',num2str(thisFOV),'.xml')));
%     else
%         thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_',num2str(thisFOV),'.xml')));
%     end
%     thisinfo=xml2struct(thispath);
%     thisposition=thisinfo.settings.acquisition.stage_position.Text;
%     C = textscan(thisposition, '%f,%f');
%     thisX=C{1}/200;
%     thisY=C{2}/200;
%     minX=min(minX,thisX);
%     minY=min(minY,thisY);
%     maxX=max(maxX,thisX);
%     maxY=max(maxY,thisY);
% end
% save(fullfile(matlabtempfiles,'minYRegion0.mat'),'minY')
% save(fullfile(matlabtempfiles,'minXRegion0.mat'),'minX')
% save(fullfile(matlabtempfiles,'maxYRegion0.mat'),'maxY')
% save(fullfile(matlabtempfiles,'maxXRegion0.mat'),'maxX')
load(fullfile(matlabtempfiles,'minYRegion0.mat'))
load(fullfile(matlabtempfiles,'minXRegion0.mat'))
load(fullfile(matlabtempfiles,'maxYRegion0.mat'))
load(fullfile(matlabtempfiles,'maxXRegion0.mat'))
milieuX=round((abs(minX)+abs(maxX))/2)+1;
milieuY=round((abs(minY)+abs(maxY))/2)+1;
registrationimage=double(zeros(floor(2048+(maxY-minY)*2048),floor(2048+(maxX-minX)*2048)));
disp('seedsmosaic size')

for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
    if thisFOV<10
        thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_000',num2str(thisFOV),'.xml')));
        thisImage=ReadDax(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_0.dax')));
        thisImage=double(thisImage(:,:,16))/2500;
    elseif thisFOV<100
        thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_00',num2str(thisFOV),'.xml')));
        thisImage=ReadDax(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_0.dax')));
        thisImage=double(thisImage(:,:,16))/2500;
    elseif thisFOV<1000
        thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_0',num2str(thisFOV),'.xml')));
        thisImage=ReadDax(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_0.dax')));
        thisImage=double(thisImage(:,:,16))/2500;
    else
        thispath=(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_',num2str(thisFOV),'.xml')));
        thisImage=ReadDax(fullfile(mypath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_0.dax')));
        thisImage=double(thisImage(:,:,16))/2500;
    end
    thisinfo=xml2struct(thispath);
    thisposition=thisinfo.settings.acquisition.stage_position.Text;
    C = textscan(thisposition, '%f,%f');
    thisX=C{1}/200+abs(minX)+1;
    thisY=C{2}/200+abs(minY)+1;
    
    dapi=flipdim(thisImage,1);
    dapi=flipdim(dapi,2);
    dapi=imrotate(dapi,-87.7);
    
    registrationimage(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+size(rotatedtemplate,1),floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+size(rotatedtemplate,1))=max(registrationimage(floor((thisY-1)*1843.2)+1:floor((thisY-1)*1843.2)+size(rotatedtemplate,1),floor((thisX-1)*1843.2)+1:floor((thisX-1)*1843.2)+size(rotatedtemplate,1)),dapi);
    
end
save(fullfile(matlabtempfiles,'registrationimageZOMBIE.mat'),'registrationimage','-v7.3')
Se=reduceImage(reduceImage(registrationimage));
save(fullfile(matlabtempfiles,'ReducedregistrationimageZombie.mat'),'Se','-v7.3')
clear
maskeMERFISHmosaictumor1forregistration

