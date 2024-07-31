% mypath='\\helium\broad_clearylab\Users\Loic\Zombiewithscaff\maxproj';
% getvaluesforthresholdfilteredimagesclearedthp1 %creates the variable mythresholds
% globalFoundSpots=cell(size(myfovs,2),1);
% globalFoundSpotsfromBW=cell(size(myfovs,2),1);
mypath='/broad/clearylab/Users/Loic/thp1homedish2_zombie/Data'
filtered='/broad/clearylab/Users/Loic/thp1homedish2_zombie/filtered'
matlabtempfiles='/broad/clearylab/Users/Loic/thp1homedish2_zombie/matlabtempfiles'
datapath='/broad/clearylab/Users/Loic/thp1homedish2_zombie/Data'
masked='/broad/clearylab/Users/Loic/thp1homedish2_zombie/maskedagain'
lowestthresholds=[100,90,120,150,140,100,120,60,90,100,220,140,210,120,120];
% lowestthresholds=[30,30,20,30,30,30,140,70,15,140,90,40,110,90,75];
mymask=zeros(2048,2048);
mymask(521,2035)=1;%remove dead pixels fom camera
mymask(657,2003)=1;
mymask(718,1990)=1;
mymask=imdilate(mymask,strel('disk',15));
mymask(1:25,:)=1;%remove edges from picture (overlap with next picture)
mymask(:,(end-25):end)=1;
mymask=imbinarize(mymask);
% lowestthresholds=table2array(readtable(fullfile(matlabtempfiles,'scalingfactorrestart1.csv')));
% lowestthresholds=ma;
delete(gcp('nocreate'))
thatparpool=parpool(22)  
myfovs=[0:1252];
parfor (fov=1:size(myfovs,2),22)
%     background=zeros(2048,2048,5);
%     for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    if rem(thisFOV,100)==0
        disp(strcat('current FOV is ',num2str(thisFOV)))%progress report
    end
    thisimage=zeros(2048,2048,15);
    for thisbit=1:15%load the image into 15 stack
        if thisFOV<10
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_');
        elseif thisFOV<100
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_');
        elseif thisFOV<1000
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_');
        else
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_');
        end
        tempimage=double(imread(fullfile(filtered,strcat(thisfilesname,num2str(thisbit),'.tif'))));
        tempimage(mymask)=0;
        thisimage(:,:,thisbit)=tempimage;
    end
%     backgroundimage=imread(fullfile(filteredpath,strcat('orangebackground_FOV_',num2str(thisFOV),'.tif')));
    %% continue
    original=zeros(2048,2048);
    bw=zeros(2048,2048,15);
%     bwbackground=zeros(2048,2048,5);
    minHW=0;
    countescapes=0;
    oldthresholds=zeros(1,15);
%     if sum(sum(nuclei))>0
%         for j=1:size(mythresholds,2)
%                 bwbackgroundimage=backgroundimage>mybackgroundthresholds(1,0.13);
%                 bwbackgroundimage=backgroundimage>mybackgroundthresholds(1,6);
%                 bwbackgroundimage=imdilate(bwbackgroundimage,strel('disk',6));
%                 countours=bwareaopen(imbinarize(backgroundimage,'adaptive'),2500);
    for frame=1:15
        thisbw=zeros(2048,2048);
        tempbw=thisimage(:,:,frame)>lowestthresholds(frame);%binarize image with round specific threshold
        bw2=imfill(tempbw,'holes');%fill holes
        bw3=bw2-bwareaopen(bw2,300);%remove spots that are too big
        bw3=bwareaopen(bw3,7);%remove spots that are too small
        bw3=bwpropfilt(bw3,"Eccentricity",[0 0.8]);%remove spots with wrong shape
        bw4=bw3+thisbw;
        bw4(bw4>0)=1;
        thisbw=bw4;
        bw(:,:,frame)=thisbw;
    end
    disp(strcat("done ",num2str(j)))
    spotsbw=decodeThisBinaryImagedontcheckrepeatsbutlookfurther(bw);%call the function that actually decodes
    if size(spotsbw,1)~=0
         spotsbw(:,end+1)=thisFOV;
         globalFoundSpotsfromBW{fov}=spotsbw;%ATTENTION : spots in position x are from image x-1
    end
end
delete(thatparpool)
save(fullfile(matlabtempfiles,strcat('globalFoundSpotsmanualV2Norepeatfurther.mat')),'globalFoundSpotsfromBW')
% save(fullfile(matlabtempfiles,'globalFoundSpots.mat'),'globalFoundSpots')

% disp('starting global preprocess')
% clear
% filterpreprocess
getgeneIDfromBWmanualthresholdsclaedlowestdontcheckrepeats
% getaveragespotIntensity
% getgeneIDopti
%% delete(gcp('nocreate'))
