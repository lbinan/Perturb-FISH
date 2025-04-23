%this script takes in the tables coming from image analysis, runs
%filtering, and formats the tables as needed for FR-Perturb and python
%plots
outputpath='\\helium\broad_clearylab\Users\Loic\tumorTablesMoreCells'
%load MERFISH data
uiopen('\\helium\broad_clearylab\Users\Loic\tumorMERFISHdish1Donor1FFPE\202408041743_tumor1DonorMerfishFFPE_VMSC11302\matlabtempfiles\merfishcounttable.csv',1)
%load gene names
uiopen('\\helium\broad_clearylab\Users\Loic\tumorMERFISHdish1Donor1FFPE\202408041743_tumor1DonorMerfishFFPE_VMSC11302\codebook_0_ImmunOncology_0.csv',1)
genenames=genenames';
%size(mymask)=169973    102389

%% just look at the data 
[X,Y]=ind2sub([102389,169973],merfishcounttable(:,3));
counts=merfishcounttable(:,4:end);
counts=counts(:,1:500);
area=merfishcounttable(:,2);
totalcounts=sum(counts,2);
bar(totalcounts)
hist(area,1000)
writematrix(counts,fullfile(outputpath,'countTableMerfish.csv'))
writematrix(area,fullfile(outputpath,'areaMerfish.csv'))
writematrix(totalcounts,fullfile(outputpath,'totalcountsMerfish.csv'))

positions=1:size(X,1);
positions=positions';
positions=[positions,X,Y];
% genefreq=sum(counts,1);
writematrix(positions,fullfile(outputpath,'positions.csv'))

%% perturbations
%filter cells to remove bad cells, mouse cells and doublets
indeces=totalcounts>40 & totalcounts<800 & area>2000 & area<30000;
zombiefiltered=zombie(indeces,:);
areafiltered=area(indeces,:);
merfishfiltered=counts(indeces,:);
corrected=zombiefiltered;
corrected(:,7)=0;corrected(:,30)=0;corrected(:,68)=0;%decoding of barcodes 7, 30 and 68 was noisy and resulted in many more decoded spots than expected, so these Ids are removed from the data

corrected(:,71)=corrected(:,71)+corrected(:,72)+corrected(:,73)+corrected(:,74);
corrected=corrected(:,1:71);
% pool together 2 guides with the same target
for i=1:35
    pooled(:,i)=corrected(:,2*i-1)+corrected(:,2*i);
end
pooled(:,36)=corrected(:,71);
writematrix(pooled,fullfile(outputpath,'pooledZombieAllcells.csv'))
writematrix(corrected,fullfile(outputpath,'notpooledZombieAllcells.csv'))
writematrix(areafiltered,fullfile(outputpath,'areafiltered.csv'))
% filter to keep only cells with identified perturbation
perturbedpooled=pooled(sum(pooled,2)>0,:);
perturbed=corrected(sum(corrected,2)>0,:);
writematrix(perturbedpooled',fullfile(outputpath,'perturbedpooledZombie.csv'))
writematrix(perturbed',fullfile(outputpath,'perturbedZombie.csv'))
perturbedMerfish=merfishfiltered(sum(pooled,2)>0,:);
writematrix(perturbedMerfish,fullfile(outputpath,'perturbedMerfish.csv'))
writematrix(merfishfiltered,fullfile(outputpath,'merfishgoodCells.csv'))
uiopen('C:\Users\lbinan\Desktop\tumor\coordinates.csv',1)
writematrix(coordinates(indeces,:),fullfile(outputpath,'positionsgoodcells.csv'))
% prepare a file with the list of perturbation in each cell, to load in
% scanpy

for i=1:size(merfishfiltered,1)
    perturbation(i)="ND";
    if sum(pooled(i,:))>0
        idx=find(pooled(i,:)>0);
        perturbation(i)=perturbedNames(idx(1));
    end
end

writematrix(perturbation',fullfile(outputpath,'perturbationName.csv'))

for i=1:size(merfishfiltered,1)
    perturbation(i)=0;
    if sum(pooled(i,:))>0
        idx=find(pooled(i,:)>0);
        perturbation(i)=idx(1);
    end
end

writematrix(perturbation',fullfile(outputpath,'perturbationIDx.csv'))
% load cell identity identify in scanpy, using leiden clustering (see notebook)
uiopen('C:\Users\lbinan\Desktop\tumor\morecells\pythonAllclusters.csv',1)
%% neighbors
% load and filter neighbor design matrix
load('\\helium\broad_clearylab\Users\Loic\tumorMERFISHdish1Donor1FFPE\202408041743_tumor1DonorMerfishFFPE_VMSC11302\matlabtempfiles\touchingnumberofneighborsoverlapUGERAGAINsmaller.mat')
%load('\\helium\broad_clearylab\Users\Loic\tumorMERFISHdish1Donor1FFPE\202408041743_tumor1DonorMerfishFFPE_VMSC11302\matlabtempfiles\touchingneighboroverlapUGERfinally.mat')
filteredneighbors=smallertable(indeces,:);
neighborsofperturbedCells=filteredneighbors(sum(pooled,2)>0,:);
cellID=1:size(filteredneighbors,1);cellID=cellID';
%immunecellIDs=cellID(pythonAllclusters(:,2)==5|pythonAllclusters(:,2)==7,:);%turns
%out 7 are doublets

%find tumor cells with a neighboring tcell
immunecellIDs=cellID(pythonAllclusters(:,2)==5,:);

immuneneighbor=zeros(size(perturbedMerfish,1),1);
for i=1:size(neighborsofperturbedCells,1)
    neighbors=neighborsofperturbedCells(i,neighborsofperturbedCells(i,:)>0);
    if size(neighbors,2)>0
        for thisneighbor=1:size(neighbors,2)
            if sum(immunecellIDs==neighbors(thisneighbor))>0
                immuneneighbor(i)=1;
            end
        end
    end
end

notimmuneneighborzombie=perturbedpooled(immuneneighbor==0,:);
notimmuneneighborMerfish=perturbedMerfish(immuneneighbor==0,:);
immuneneighborzombie=perturbedpooled(immuneneighbor==1,:);
immuneneighborMerfish=perturbedMerfish(immuneneighbor==1,:);
writematrix(notimmuneneighborzombie',fullfile(outputpath,'notimmuneneighborzombie.csv'))
writematrix(notimmuneneighborMerfish,fullfile(outputpath,'notimmuneneighborMerfish.csv'))
writematrix(immuneneighborzombie',fullfile(outputpath,'immuneneighborzombie.csv'))
writematrix(immuneneighborMerfish,fullfile(outputpath,'immuneneighborMerfish.csv'))

writematrix(notimmuneneighborzombie',fullfile(outputpath,'notimmuneneighborzombieonly5.csv'))
writematrix(notimmuneneighborMerfish,fullfile(outputpath,'notimmuneneighborMerfishonly5.csv'))
writematrix(immuneneighborzombie',fullfile(outputpath,'immuneneighborzombieonly5.csv'))
writematrix(immuneneighborMerfish,fullfile(outputpath,'immuneneighborMerfishonly5.csv'))

%% center on Tcells, look for perturbaed tumor cell in the neighbohood
% merfishTcells=merfishfiltered(pythonAllclusters(:,2)==5|pythonAllclusters(:,2)==7,:);
% immunecellsneighbors=filteredneighbors(pythonAllclusters(:,2)==5|pythonAllclusters(:,2)==7,:);
uiopen('C:\Users\lbinan\Desktop\tumor\morecells\pythonimmunebis12.csv',1)
% merfishTcells=merfishfiltered(pythonAllclusters(:,2)==5|pythonAllclusters(:,2)==7,:);
merfishTcells=merfishfiltered(pythonAllclusters(:,2)==5,:);%cluster 7 in fact contains a lot of doublets
merfishTcells=merfishTcells(pythonimmunebis1(:,2)==1|pythonimmunebis1(:,2)==3|pythonimmunebis1(:,2)==4,:);
immunecellsneighbors=filteredneighbors(pythonAllclusters(:,2)==5,:);
cellID=1:size(smallertable,1);cellID=cellID';
cellIDfiltered=cellID(indeces,:);
PerturbedcellIDs=cellIDfiltered(sum(pooled,2)>0,:);

cancerneighbortoTcell=zeros(size(merfishTcells,1),36);
%build the "perturbation" matrix: for each t cell (row) we have a 1 in the
%column of the perturbation received by neighboring tumor cells
for i=1:size(merfishTcells,1)
    neighbors=immunecellsneighbors(i,immunecellsneighbors(i,:)>0);
    if size(neighbors,2)>0
        for thisneighbor=1:size(neighbors,2)
            currentneighbor=neighbors(thisneighbor);
            if sum(PerturbedcellIDs==currentneighbor)>0
                cancerneighbortoTcell(i,:)=cancerneighbortoTcell(i,:)+pooled(find(PerturbedcellIDs==currentneighbor),:);
            end
        end
    end
end
tcellnextoperturbed=cancerneighbortoTcell(sum(cancerneighbortoTcell,2)>0,:);

%% all neighbors

correctedraw=zombie;
correctedraw(:,7)=0;correctedraw(:,30)=0;correctedraw(:,68)=0;

correctedraw(:,71)=correctedraw(:,71)+correctedraw(:,72)+correctedraw(:,73)+correctedraw(:,74);
correctedraw=correctedraw(:,1:71);

for i=1:35
    pooledraw(:,i)=correctedraw(:,2*i-1)+correctedraw(:,2*i);
end
pooledraw(:,36)=correctedraw(:,71);

neighborpertubations=zeros(size(pooledraw,1),36);

for i=1:size(pooledraw,1)
    neighbors=smallertable(i,smallertable(i,:)>0);
    if size(neighbors,2)>0
        for thisneighbor=1:size(neighbors,2)
            currentneighbor=neighbors(thisneighbor);
                neighborpertubations(i,:)=neighborpertubations(i,:)+pooledraw(currentneighbor,:);
        end
    end
end
neighborperturbationfiltered=neighborpertubations(indeces,:);
% tcellneighborsdesign=neighborperturbationfiltered(pythonAllclusters(:,2)==5|pythonAllclusters(:,2)==7,:);
tcellneighborsdesign=neighborperturbationfiltered(pythonAllclusters(:,2)==5,:);
tcellneighborsdesign=tcellneighborsdesign(~badcells,:);
perturbedtcellszombie=tcellneighborsdesign(sum(tcellneighborsdesign,2)>0,:);
perturbedtcellszombieonlycontrol=perturbedtcellszombie;
perturbedtcellszombieonlycontrol(sum(perturbedtcellszombie(:,1:35),2)>0,36)=0;
perturbedtcellsmerfish=merfishTcells(sum(tcellneighborsdesign,2)>0,:);
otherTcellsmerfish=merfishTcells(sum(tcellneighborsdesign,2)==0,:);
writematrix(perturbedtcellszombieonlycontrol',fullfile(outputpath,'perturbedTcellsZombieonly5.csv'))
writematrix(perturbedtcellsmerfish,fullfile(outputpath,'perturbedtcellsmerfishonly5.csv'))


distribution=sum(perturbedtcellszombieonlycontrol,1);
perturbedtcellszombieonlycontrolenoughguides=perturbedtcellszombieonlycontrol(:,distribution>40);
keepguides=distribution>40;%keep only guides with enough detections
bar(distribution)
expressionimmune=mean(merfishTcells,1);
sum(expressionimmune>0.3)%keep genes expressed in tcelsl only

writematrix([perturbedtcellsmerfish(sum(perturbedtcellszombieonlycontrolenoughguides,2)>0,expressionimmune>0.3);otherTcellsmerfish(:,expressionimmune>0.3)],fullfile(outputpath,'perturbedtcellsmerfishenoughguidesonly5NobadcellsEVENSTRICTER0917includeallcells.csv'))
extracells=zeros(size(otherTcellsmerfish,1),sum(keepguides));
extracells(:,sum(keepguides))=1;
writematrix([(perturbedtcellszombieonlycontrolenoughguides(sum(perturbedtcellszombieonlycontrolenoughguides,2)>0,:)>0);extracells]',fullfile(outputpath,'perturbedTcellsZombieenoughguidesonly5NobadcellsEVENSTRICTER0917includeallcells.csv'))

guides(keepguides)
genenames(expressionimmune>0.3)
save(fullfile(outputpath,'workspace.mat'))

