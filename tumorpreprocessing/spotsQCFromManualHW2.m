% matlabtempfiles='\\helium\broad_clearylab\Users\Loic\thp1ubi_encodetogether_gel_clear_bleachAt4CaWhile\matlabtempfiles'
% filtered='\\helium\broad_clearylab\Users\Loic\thp1ubi_encodetogether_gel_clear_bleachAt4CaWhile\filtered'
% whichtoQC=6;
% matlabtempfiles='/broad/clearylab/Users/Loic/asddish1_zombie/asddish1/matlabtempfiles'
delete(gcp('nocreate'))

clear mystats
% decoded=readtable(fullfile(matlabtempfiles,strcat('decodedGuidesdeocdedonBWmoresensitive',num2str(whichtoQC),'.csv')));
mypath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/Data'
filtered='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/filtered'
matlabtempfiles='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles'
datapath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/tumorZombiedish1Donor1FFPE\Data'
decoded=readtable(fullfile(matlabtempfiles,'decodedGuidesFromManualHW2.csv'));
% filtered='/broad/clearylab/Users/Loic/asddish1_zombie/asddish1/filtered'

decoded=table2array(decoded);
thatparpool=parpool(20)  
parfor (thisspot=1:size(decoded,1),20)
thisY=decoded(thisspot,2);
thisX=decoded(thisspot,1);

thisFOV=decoded(thisspot,3);
    thisimage=zeros(2048,2048,15);
    spotimage=zeros(17,17,15);
    thisbw=zeros(17,17,15);
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
%     tempimage(mymask)=0;
    thisimage(:,:,thisbit)=tempimage;
    end

spotimage=thisimage((thisX-8):(thisX+8),(thisY-8):(thisY+8),1:15);
% for i=1:15
%     thisbw(:,:,i)=spotimage(:,:,i)>mythresholds(i,whichtoQC);
% end
for i=1:15%only temporary for manuual thresh
    thisbw(:,:,i)=spotimage(:,:,i)>120;
end
bw=sum(thisbw,3);
thisarea=sum(sum(bw==max(max(bw))));
thismean=mean(spotimage(9,9,:));
thisstd=std(spotimage(9,9,:));
thismax=max(spotimage(9,9,:));
thismin=min(spotimage(9,9,:));
thisseq=zeros(1,21);
for j=1:15
thisseq(j)=spotimage(9,9,j)
end
mystats(thisspot,:)=[decoded(thisspot,4),decoded(thisspot,3),thisarea,thismean,thisstd,thismax,thismin,sum(thisseq),thisseq];
end
delete(thatparpool)
% writematrix(mystats,fullfile(matlabtempfiles,strcat('qcStatsthresh',num2str(whichtoQC),'.csv')))
writematrix(mystats,fullfile(matlabtempfiles,strcat('qcStatsManualAdjust','QCFromManualHW2.csv')))

delete(gcp('nocreate'))


%% play with data
% figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.')
% hold on
% scatter(qcStatsthresh7(qcStatsthresh7(:,1)>75,5),qcStatsthresh7(qcStatsthresh7(:,1)>75,4),'r')
% scatter(qcStatsthresh7(qcStatsthresh7(:,8)==5,5),qcStatsthresh7(qcStatsthresh7(:,8)==5,4),'o')
% figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.')
% hold on
% scatter(qcStatsthresh7(qcStatsthresh7(:,1)>75,5),qcStatsthresh7(qcStatsthresh7(:,1)>75,4),'r')
% 
% x = linspace(0,120,2);
% y = 2.1*x-97.582;
% plot(x,y)
% 
% title('example filter, circles = blanks')
% qcStatsthresh7=qcStatsthresh7(:,1:23);
% for i=1:size(qcStatsthresh7,1)
%     qcStatsthresh7(i,24)=2.1*qcStatsthresh7(i,5)-97.582-qcStatsthresh7(i,4);
% end
% figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.')
% hold on
% scatter(qcStatsthresh7(qcStatsthresh7(:,24)>0,5),qcStatsthresh7(qcStatsthresh7(:,24)>0,4),'r')
% title('kept spots (10510 out of 16844)')
% sum(qcStatsthresh7(:,24)>0)
% size(qcStatsthresh7)
% filtered=qcStatsthresh7(qcStatsthresh7(:,24)>0,:);
% 
% figure, scatter(qcStatsthresh7(:,5),qcStatsthresh7(:,4),'.')
% hold on
% scatter(qcStatsthresh7(qcStatsthresh7(:,8)==4,5),qcStatsthresh7(qcStatsthresh7(:,8)==4,4),'r')
% 

