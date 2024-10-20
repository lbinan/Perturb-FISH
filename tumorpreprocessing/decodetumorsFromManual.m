% mypath='\\helium\broad_clearylab\Users\Loic\Zombiewithscaff\maxproj';
% getvaluesforthresholdfilteredimagesclearedthp1 %creates the variable mythresholds
% globalFoundSpots=cell(size(myfovs,2),1);
% globalFoundSpotsfromBW=cell(size(myfovs,2),1);
mypath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/Data'
filtered='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/filtered'
matlabtempfiles='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles'
datapath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/tumorZombiedish1Donor1FFPE\Data'
% masked='/broad/clearylab/Users/Loic/tumorszombie2/maskedagain'
lowestthresholds=[100,90,120,150,140,100,120,60,90,100,220,140,210,120,120];
lowestthresholds=[120,110,140,170,160,120,140,80,110,120,240,160,230,140,140];
lowestthresholds=[227,245,216,255,255,24,171,173,207,255,240,140,148,177,190];
lowestthresholds=[30,30,20,30,30,30,140,70,15,140,90,40,110,90,75];
lowestthresholds=[30,80,80,30,80,80,80,60,80,80,90,70,70,70,220];
lowestthresholds=[45,70,90,45,70,80,70,50,70,70,70,80,70,60,200];
lowestthresholds=[40,60,85,45,60,80,60,40,60,65,80,80,70,45,115];
lowestthresholds=[30,50,50,30,50,45,45,40,30,50,50,40,50,45,35];
% lowestthresholds=[30,85,65,40,60,85,50,20,30,50,90,85,70,40,90];
lowestthresholds=[17,71,99,19,52,101,52,16,33,23,79,108,52,35,89];
lowestthresholds=[27,80,80,20,40,60,52,20,25,25,65,65,35,20,89];
lowestthresholds=[40,71,60,80,50,80,60,30,40,35,65,60,60,35,70];
lowestthresholds=[30,60,50,60,50,70,50,25,30,45,75,80,50,25,90];
lowestthresholds=[25,65,55,50,50,75,50,25,35,45,75,70,60,35,100];
% lowestthresholds=[27,60,60,20,40,50,52,20,35,25,50,55,55,55,89];
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
thatparpool=parpool(20)  
myfovs=[0:1169];
parfor (fov=1:size(myfovs,2),20)
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
%         imwrite(bw2,fullfile(testpath,strcat('testbw2',num2str(frame),'.png')))
        bw3=bw2-bwareaopen(bw2,250);%remove spots that are too big
%         imwrite(bw3,fullfile(testpath,strcat('testnobig',num2str(frame),'.png')))
        bw3=bwareaopen(bw3,8);%remove spots that are too small was 7
        bw3=bwpropfilt(bw3,"Eccentricity",[0 0.8]);%remove spots with wrong shape
%         imwrite(bw3,fullfile(testpath,strcat('testfiltered',num2str(frame),'.png')))
        bw4=bw3+thisbw;
        bw4(bw4>0)=1;
        thisbw=bw4;
        bw(:,:,frame)=thisbw;
%         imwrite(bw(:,:,frame),fullfile(testpath,strcat('testfinal',num2str(frame),'.png')))
    end
%     disp(strcat("done ",num2str(j)))
    spotsbw=decodeThisBinaryImageTUMOR(bw);%call the function that actually decodes
    altspots=decodeThisBinaryImageTUMORHW2(bw);
    if size(spotsbw,1)~=0
         spotsbw(:,end+1)=thisFOV;
         globalFoundSpotsfromBW{fov}=spotsbw;%ATTENTION : spots in position x are from image x-1
    end
    if size(altspots,1)~=0
         altspots(:,end+1)=thisFOV;
         globalFoundSpotsfromBWHW2{fov}=altspots;%ATTENTION : spots in position x are from image x-1
    end
end
delete(thatparpool)
delete(gcp('nocreate'))
save(fullfile(matlabtempfiles,strcat('globalFoundFromManual.mat')),'globalFoundSpotsfromBW')
save(fullfile(matlabtempfiles,strcat('globalFoundFromManualHW2.mat')),'globalFoundSpotsfromBWHW2')
% save(fullfile(matlabtempfiles,'globalFoundSpots.mat'),'globalFoundSpots')

% disp('starting global preprocess')
% clear
% filterpreprocess
getgeneIDFromManual
getgeneIDFromManualHW2
disp('finished getgenes')

% getaveragespotIntensity
% getgeneIDopti
%% delete(gcp('nocreate'))
spotsQCFromManual
spotsQCFromManualHW2
disp('finished spotsQCtumors')
delete(gcp('nocreate'))
%%
