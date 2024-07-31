delete(gcp('nocreate'))

clear mystats
decoded=readtable(fullfile(matlabtempfiles,'decodedGuidesmanualV2dontcheckrepeatslookfurther.csv'));%load output of decoding script
filtered='/broad/clearylab/Users/Loic/thp1homedish2_zombie/filteredagain'%path to gaussian filtered images

decoded=table2array(decoded);
thatparpool=parpool(22)  
parfor (thisspot=1:size(decoded,1),22)%here, for each decoded guide, we load the actual intensity values from the images
    thisY=decoded(thisspot,2);
    thisX=decoded(thisspot,1);
    thisFOV=decoded(thisspot,3);
    thisimage=zeros(2048,2048,15);
    spotimage=zeros(11,11,15);
    thisbw=zeros(11,11,15);
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
        thisimage(:,:,thisbit)=tempimage;
    end
    spotimage=thisimage((thisX-5):(thisX+5),(thisY-5):(thisY+5),1:15);
    for i=1:15%only temporary for manuual thresh
        thisbw(:,:,i)=spotimage(:,:,i)>110;
    end
    bw=sum(thisbw,3);
    thisarea=sum(sum(bw==max(max(bw))));
    thismean=mean(spotimage(6,6,:));
    thisstd=std(spotimage(6,6,:));
    thismax=max(spotimage(6,6,:));
    thismin=min(spotimage(6,6,:));
    thisseq=zeros(1,15);
    for j=1:15
        thisseq(j)=spotimage(6,6,j)
    end
    mystats(thisspot,:)=[decoded(thisspot,4),decoded(thisspot,3),thisarea,thismean,thisstd,thismax,thismin,sum(thisseq),thisseq];
end
delete(thatparpool)
% writematrix(mystats,fullfile(matlabtempfiles,strcat('qcStatsthresh',num2str(whichtoQC),'.csv')))
writematrix(mystats,fullfile(matlabtempfiles,strcat('qcStatsthresh','QClowestmanualfiltereddecreasing.csv')))
delete(gcp('nocreate'))


%% filter wrong calls
figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.') %scatter guides using STD and average intensity
%identify the cloud that contains blank barcodes, get the equation of the
%line separating the 2 clouds, and filter out the wrong calls
hold on
scatter(qcStatsthresh7(qcStatsthresh7(:,1)>75,5),qcStatsthresh7(qcStatsthresh7(:,1)>75,4),'r')
% scatter(qcStatsthresh7(qcStatsthresh7(:,8)==5,5),qcStatsthresh7(qcStatsthresh7(:,8)==5,4),'o')
% figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.')
% hold on
% scatter(qcStatsthresh7(qcStatsthresh7(:,1)>75,5),qcStatsthresh7(qcStatsthresh7(:,1)>75,4),'r')
% 
x = linspace(0,120,2);
y = 2.1*x-97.582;
plot(x,y)
% 
% title('example filter, circles = blanks')
qcStatsthresh7=qcStatsthresh7(:,1:23);
for i=1:size(qcStatsthresh7,1)
    qcStatsthresh7(i,24)=2.1*qcStatsthresh7(i,5)-97.582-qcStatsthresh7(i,4);
end
figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.')
hold on
scatter(qcStatsthresh7(qcStatsthresh7(:,24)>0,5),qcStatsthresh7(qcStatsthresh7(:,24)>0,4),'r')
title('kept spots (10510 out of 16844)')
sum(qcStatsthresh7(:,24)>0)
size(qcStatsthresh7)
filtered=qcStatsthresh7(qcStatsthresh7(:,24)>0,:);%save filtered as input used to build count tables
% 
% figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.')
% hold on
% scatter(qcStatsthresh7(qcStatsthresh7(:,8)==4,5),qcStatsthresh7(qcStatsthresh7(:,8)==4,4),'r')
% 

