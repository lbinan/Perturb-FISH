thiscycle=1%runs through a sample of all images to set thresholds for intensity for spot detection (final thresholds were further adjusted manually)
numberofimagestoinitialiseon=600
thisparpool=parpool(22)  
mymask=zeros(2048,2048);
mymask(521,2035)=1;%remove dead pixels on camera
mymask(657,2003)=1;
mymask(718,1990)=1;
mymask=imdilate(mymask,strel('disk',15));
mymask(1:25,:)=1;%also remove edges
mymask(:,(end-25):end)=1;
mymask=imbinarize(mymask);
parfor (thisFOV=1:numberofimagestoinitialiseon,22) %look at ~40 images per channel to find the right threshold
    if thisFOV<10
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_');
    elseif thisFOV<100
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_');
    elseif thisFOV<1000
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_');
    else
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_');
    end
    thisimage=zeros(2048,2048,3);
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'1.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,1)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'2.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,2)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'3.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,3)=tempimage;
    imsbit1(:,:,thisFOV)=thisimage(:,:,1);
    imsbit2(:,:,thisFOV)=thisimage(:,:,2);
    imsbit3(:,:,thisFOV)=thisimage(:,:,3);
    
end
% imsbit1(imsbit1==1)=0;
% imsbit2(imsbit2==1)=0;
% imsbit3(imsbit3==1)=0;
h=histogram(imsbit1);%this approach is taken from "Spatially resolved epigenomic profiling of single
%cells in complex tissues"
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end  
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(1,:)=thisthresh;
h=histogram(imsbit2);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(2,:)=thisthresh;
h=histogram(imsbit3);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(3,:)=thisthresh;
% save(fullfile(matlabtempfiles,'imsbit1.mat'),'imsbit1')
% save(fullfile(matlabtempfiles,'imsbit2.mat'),'imsbit2')
% save(fullfile(matlabtempfiles,'imsbit3.mat'),'imsbit3')
clear imsbit1 imsbit2 imsbit3
%%
thiscycle=2;
imsbit1nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
imsbit2nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
imsbit3nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
parfor (thisFOV=1:numberofimagestoinitialiseon,22) %look at ~40 images per channel to find the right threshold
    if thisFOV<10
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_');
    elseif thisFOV<100
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_');
    elseif thisFOV<1000
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_');
    else
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_');
    end
        thisimage=zeros(2048,2048,3);
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'4.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,1)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'5.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,2)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'6.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,3)=tempimage;
    imsbit1nextcycles(:,:,thisFOV)=thisimage(:,:,1);
    imsbit2nextcycles(:,:,thisFOV)=thisimage(:,:,2);
    imsbit3nextcycles(:,:,thisFOV)=thisimage(:,:,3);
end

% imsbit1nextcycles(imsbit1nextcycles==1)=0;
% imsbit2nextcycles(imsbit2nextcycles==1)=0;
% imsbit3nextcycles(imsbit3nextcycles==1)=0;
h=histogram(imsbit1nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(4,:)=thisthresh;
h=histogram(imsbit2nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(5,:)=thisthresh;
h=histogram(imsbit3nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(6,:)=thisthresh;
% save(fullfile(matlabtempfiles,'imsbit1nextcyclescluster.mat'),'imsbit1nextcycles')
% save(fullfile(matlabtempfiles,'imsbit2nextcyclescluster.mat'),'imsbit2nextcycles')
% save(fullfile(matlabtempfiles,'imsbit3nextcyclescluster.mat'),'imsbit3nextcycles')
clear imsbit1nextcycles imsbit2nextcycles imsbit3nextcycles imsbackground

%% 
thiscycle=3;
imsbit7nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
imsbit8nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
imsbit9nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
parfor (thisFOV=1:numberofimagestoinitialiseon,22) %look at ~40 images per channel to find the right threshold
    if thisFOV<10
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_');
    elseif thisFOV<100
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_');
    elseif thisFOV<1000
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_');
    else
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_');
    end
        thisimage=zeros(2048,2048,3);
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'7.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,1)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'8.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,2)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'9.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,3)=tempimage;
    imsbit7nextcycles(:,:,thisFOV)=thisimage(:,:,1);
    imsbit8nextcycles(:,:,thisFOV)=thisimage(:,:,2);
    imsbit9nextcycles(:,:,thisFOV)=thisimage(:,:,3);
end

% imsbit7nextcycles(imsbit7nextcycles==1)=0;
% imsbit8nextcycles(imsbit8nextcycles==1)=0;
% imsbit9nextcycles(imsbit9nextcycles==1)=0;
h=histogram(imsbit7nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(7,:)=thisthresh;
h=histogram(imsbit8nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(8,:)=thisthresh;
h=histogram(imsbit9nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(9,:)=thisthresh;
% save(fullfile(matlabtempfiles,'imsbit7nextcyclescluster.mat'),'imsbit7nextcycles')
% save(fullfile(matlabtempfiles,'imsbit8nextcyclescluster.mat'),'imsbit8nextcycles')
% save(fullfile(matlabtempfiles,'imsbit9nextcyclescluster.mat'),'imsbit9nextcycles')
clear imsbit7nextcycles imsbit8nextcycles imsbit9nextcycles imsbackground
%% 
thiscycle=4;
% imsbit10nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
% imsbit11nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
% imsbit12nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
parfor (thisFOV=1:numberofimagestoinitialiseon,22) %look at ~40 images per channel to find the right threshold
    if thisFOV<10
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_');
    elseif thisFOV<100
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_');
    elseif thisFOV<1000
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_');
    else
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_');
    end
    thisimage=zeros(2048,2048,3);
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'10.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,1)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'11.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,2)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'12.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,3)=tempimage;
    imsbit10nextcycles(:,:,thisFOV)=thisimage(:,:,1);
    imsbit11nextcycles(:,:,thisFOV)=thisimage(:,:,2);
    imsbit12nextcycles(:,:,thisFOV)=thisimage(:,:,3);
end

% imsbit10nextcycles(imsbit10nextcycles==1)=0;
% imsbit11nextcycles(imsbit11nextcycles==1)=0;
% imsbit12nextcycles(imsbit12nextcycles==1)=0;
h=histogram(imsbit10nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(10,:)=thisthresh;
h=histogram(imsbit11nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(11,:)=thisthresh;
h=histogram(imsbit12nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(12,:)=thisthresh;
% save(fullfile(matlabtempfiles,'imsbit10nextcyclescluster.mat'),'imsbit10nextcycles')
% save(fullfile(matlabtempfiles,'imsbit11nextcyclescluster.mat'),'imsbit11nextcycles')
% save(fullfile(matlabtempfiles,'imsbit12nextcyclescluster.mat'),'imsbit12nextcycles')
clear imsbit10nextcycles imsbit11nextcycles imsbit12nextcycles imsbackground
%% 
thiscycle=5;
imsbit13nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
imsbit14nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
imsbit15nextcycles=zeros(2048,2048,numberofimagestoinitialiseon-19);
parfor (thisFOV=1:numberofimagestoinitialiseon,22) %look at ~40 images per channel to find the right threshold
    if thisFOV<10
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_');
    elseif thisFOV<100
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_');
    elseif thisFOV<1000
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_');
    else
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_');
    end  
    thisimage=zeros(2048,2048,3);
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'13.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,1)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'14.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,2)=tempimage;
    tempimage=imread(fullfile(filtered,strcat(thisfilesname,'15.tif')));
    tempimage(mymask)=0;
    thisimage(:,:,3)=tempimage;
    imsbit13nextcycles(:,:,thisFOV)=thisimage(:,:,1);
    imsbit14nextcycles(:,:,thisFOV)=thisimage(:,:,2);
    imsbit15nextcycles(:,:,thisFOV)=thisimage(:,:,3);
end

% imsbit13nextcycles(imsbit13nextcycles==1)=0;
% imsbit14nextcycles(imsbit14nextcycles==1)=0;
% imsbit15nextcycles(imsbit15nextcycles==1)=0;
h=histogram(imsbit13nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(13,:)=thisthresh;
h=histogram(imsbit14nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(14,:)=thisthresh;
h=histogram(imsbit15nextcycles);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
    ccum(i)=sum(h.BinCounts(1:i));
end   
idx=find(ccum>0.99999*ccum(end));  
thisthresh=h.BinEdges(idx(1)) ;
for threshold=0.99995:-0.0001:0.98
    idx=find(ccum>threshold*ccum(end));  
    thisthresh=[thisthresh;h.BinEdges(idx(1))];
end
mythresholds(15,:)=thisthresh;
% save(fullfile(matlabtempfiles,'imsbit13nextcyclescluster.mat'),'imsbit13nextcycles')
% save(fullfile(matlabtempfiles,'imsbit14nextcyclescluster.mat'),'imsbit14nextcycles')
% save(fullfile(matlabtempfiles,'imsbit15nextcyclescluster.mat'),'imsbit15nextcycles')
clear imsbit13nextcycles imsbit14nextcycles imsbit15nextcycles imsbackground
delete(thisparpool)
save(fullfile(matlabtempfiles,'mythresholds.mat'),'mythresholds')

