% matlabtempfiles='\\helium\broad_clearylab\Users\Loic\thp1ubi_encodetogether_gel_clear_bleachAt4CaWhile\matlabtempfiles'
% filtered='\\helium\broad_clearylab\Users\Loic\thp1ubi_encodetogether_gel_clear_bleachAt4CaWhile\filtered'
whichtoQC=6;
% matlabtempfiles='/broad/clearylab/Users/Loic/asddish1_zombie/asddish1/matlabtempfiles'
delete(gcp('nocreate'))

clear mystats
% decoded=readtable(fullfile(matlabtempfiles,strcat('decodedGuidesdeocdedonBWmoresensitive',num2str(whichtoQC),'.csv')));
mypath='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/Data'%path to raw images (dax files)
filtered='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/filtered'%output path
matlabtempfiles='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/matlabtempfiles'
decoded=readtable(fullfile(matlabtempfiles,'decodedGuidesManualAdjust.csv'));
% filtered='/broad/clearylab/Users/Loic/asddish1_zombie/asddish1/filtered'

decoded=table2array(decoded);
thatparpool=parpool(22)  
parfor (thisspot=1:size(decoded,1),22)
thisY=decoded(thisspot,2);
thisX=decoded(thisspot,1);

thisFOV=decoded(thisspot,3);
    thisimage=zeros(2048,2048,21);
    spotimage=zeros(11,11,21);
    thisbw=zeros(11,11,21);
    for thisbit=1:21%load the image into 15 stack
      if thisFOV<10
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_000',num2str(thisFOV),'_');
        elseif thisFOV<100
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_00',num2str(thisFOV),'_');
        elseif thisFOV<1000
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_0',num2str(thisFOV),'_');
        else
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_',num2str(thisFOV),'_');
      end
    tempimage=double(imread(fullfile(filtered,strcat(thisfilesname,num2str(thisbit),'.tif'))));
%     tempimage(mymask)=0;
    thisimage(:,:,thisbit)=tempimage;
    end

spotimage=thisimage((thisX-5):(thisX+5),(thisY-5):(thisY+5),1:21);
% for i=1:15
%     thisbw(:,:,i)=spotimage(:,:,i)>mythresholds(i,whichtoQC);
% end
for i=1:21%only temporary for manuual thresh
    thisbw(:,:,i)=spotimage(:,:,i)>120;
end
bw=sum(thisbw,3);
thisarea=sum(sum(bw==max(max(bw))));
thismean=mean(spotimage(6,6,:));
thisstd=std(spotimage(6,6,:));
thismax=max(spotimage(6,6,:));
thismin=min(spotimage(6,6,:));
thisseq=zeros(1,21);
for j=1:21
thisseq(j)=spotimage(6,6,j)
end
mystats(thisspot,:)=[decoded(thisspot,4),decoded(thisspot,3),thisarea,thismean,thisstd,thismax,thismin,sum(thisseq),thisseq];
end
delete(thatparpool)
% writematrix(mystats,fullfile(matlabtempfiles,strcat('qcStatsthresh',num2str(whichtoQC),'.csv')))
writematrix(mystats,fullfile(matlabtempfiles,strcat('qcStatsManualAdjust','QC1manualAdjustlarger.csv')))

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

