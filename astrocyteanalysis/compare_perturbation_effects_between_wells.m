mypath='\\helium\broad_clearylab\Users\Loic\AsdRevstables'
thistable=(readtable(fullfile(mypath,'goodgenes_well1_LFCs.csv')));
effectsWell1=table2array(thistable(:,2:end));
thistable=readtable(fullfile(mypath,'goodgenes_well1_qvals.csv'));
qWell1=table2array(thistable(:,2:end));

thistable=(readtable(fullfile(mypath,'goodgenes_well2_LFCs.csv')));
effectsWell2=table2array(thistable(:,2:end));
thistable=readtable(fullfile(mypath,'goodgenes_well2_qvals.csv'));
qWell2=table2array(thistable(:,2:end));

thistable=(readtable(fullfile(mypath,'goodgenes_well3_LFCs.csv')));
effectsWell3=table2array(thistable(:,2:end));
thistable=readtable(fullfile(mypath,'goodgenes_well3_qvals.csv'));
qWell3=table2array(thistable(:,2:end));
otherpath='\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes'
thistable=(readtable(fullfile(otherpath,'all_LFCs.csv')));
effectsGlobal=table2array(thistable(:,2:end));
thistable=readtable(fullfile(otherpath,'all_qvals.csv'));
qGlobal=table2array(thistable(:,2:end));


flateffectsWell1=reshape(effectsWell1,[277*129,1]);
flatqWell1=reshape(qWell1,[277*129,1]);

flateffectsWell2=reshape(effectsWell2,[277*129,1]);
flatqWell2=reshape(qWell2,[277*129,1]);

flateffectsWell3=reshape(effectsWell3,[277*129,1]);
flatqWell3=reshape(qWell3,[277*129,1]);

flateffectsGlobal=reshape(effectsGlobal,[277*129,1]);
flatqGlobal=reshape(qGlobal,[277*129,1]);


subplot(3,3,1)
scatter(flateffectsWell1,flateffectsWell2,'.')
title('All effects')
xlabel('Well 1')
ylabel('Well 2')
C=corrcoef(flateffectsWell1,flateffectsWell2);
annotation('textbox',[0.15 0.92 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(3,3,2)
scatter(flateffectsWell1,flateffectsWell3,'.')
title('All effects')
xlabel('Well 1')
ylabel('Well 3')
C=corrcoef(flateffectsWell1,flateffectsWell3);
annotation('textbox',[0.45 0.92 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(3,3,3)
scatter(flateffectsWell2,flateffectsWell3,'.')
title('All effects')
xlabel('Well 2')
ylabel('Well 3')
C=corrcoef(flateffectsWell2,flateffectsWell3);
annotation('textbox',[0.72 0.92 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     


subplot(3,3,4)
scatter(flateffectsWell1(flatqWell1<0.1 & flatqWell2< 0.1),flateffectsWell2(flatqWell1<0.1 & flatqWell2< 0.1),'.')
title('Effects significant in both datasets')
xlabel('Well 1')
ylabel('Well 2')
C=corrcoef(flateffectsWell1(flatqWell1<0.1 & flatqWell2< 0.1),flateffectsWell2(flatqWell1<0.1 & flatqWell2< 0.1));
annotation('textbox',[0.15 0.62 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(3,3,5)
scatter(flateffectsWell1(flatqWell1<0.1 & flatqWell3< 0.1),flateffectsWell3(flatqWell1<0.1 & flatqWell3< 0.1),'.')
title('Effects significant in both datasets')
xlabel('Well 1')
ylabel('Well 3')
C=corrcoef(flateffectsWell1(flatqWell1<0.1 & flatqWell3< 0.1),flateffectsWell3(flatqWell1<0.1 & flatqWell3< 0.1));
annotation('textbox',[0.45 0.62 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(3,3,6)
scatter(flateffectsWell2(flatqWell2<0.1 & flatqWell3< 0.1),flateffectsWell3(flatqWell2<0.1 & flatqWell3< 0.1),'.')
title('Effects significant in both datasets')
xlabel('Well 2')
ylabel('Well 3')
C=corrcoef(flateffectsWell2(flatqWell2<0.1 & flatqWell3< 0.1),flateffectsWell3(flatqWell2<0.1 & flatqWell3< 0.1));
annotation('textbox',[0.72 0.62 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     


subplot(3,3,7)
scatter(flateffectsWell1(flatqWell1<0.1 | flatqWell2<0.1),flateffectsWell2(flatqWell1<0.1 | flatqWell2<0.1),'.')
title('Effects significant in either dataset')
xlabel('Well 1')
ylabel('Well 2')
C=corrcoef(flateffectsWell1(flatqWell1<0.1 | flatqWell2<0.1),flateffectsWell2(flatqWell1<0.1 | flatqWell2<0.1));
annotation('textbox',[0.15 0.32 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(3,3,8)
scatter(flateffectsWell1(flatqWell1<0.1 | flatqWell3<0.1),flateffectsWell3(flatqWell1<0.1 | flatqWell3<0.1),'.')
title('Effects significant in either dataset')
xlabel('Well 1')
ylabel('Well 3')
C=corrcoef(flateffectsWell1(flatqWell1<0.1 | flatqWell3<0.1),flateffectsWell3(flatqWell1<0.1 | flatqWell3<0.1));
annotation('textbox',[0.45 0.32 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(3,3,9)
scatter(flateffectsWell2(flatqWell2<0.1 | flatqWell3<0.1),flateffectsWell3(flatqWell2<0.1 | flatqWell3<0.1),'.')
title('Effects significant in either dataset')
xlabel('Well 2')
ylabel('Well 3')
C=corrcoef(flateffectsWell2(flatqWell2<0.1 | flatqWell3<0.1),flateffectsWell3(flatqWell2<0.1 | flatqWell3<0.1));
annotation('textbox',[0.72 0.32 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

set(gcf,'color','w');



oldpath='\\helium\broad_clearylab\Users\Loic\AsdRevstables'
thistable=(readtable(fullfile(mypath,'old_data_LFCs.csv')));
effectsOldWell1=table2array(thistable(:,2:end));
thistable=readtable(fullfile(mypath,'old_data_qvals.csv'));
qOldWell1=table2array(thistable(:,2:end));

oldpath='\\helium\broad_clearylab\Users\Loic\AsdRevstables'
thistable=(readtable(fullfile(mypath,'old_datadish3_LFCs.csv')));
effectsOldWell3=table2array(thistable(:,2:end));
thistable=readtable(fullfile(mypath,'old_datadish3_qvals.csv'));
qOldWell3=table2array(thistable(:,2:end));

oldpath='\\helium\broad_clearylab\Users\Loic\AsdRevstables'
thistable=(readtable(fullfile(mypath,'all_old_data_LFCs.csv')));
effectsPooledOld=table2array(thistable(:,2:end));
thistable=readtable(fullfile(mypath,'all_old_data_qvals.csv'));
qOldPooledOld=table2array(thistable(:,2:end));


flateffectsOldWell1=reshape(effectsOldWell1,[277*129,1]);
flatqOldWell1=reshape(qOldWell1,[277*129,1]);
flateffectsOldWell3=reshape(effectsOldWell3,[277*129,1]);
flatqOldWell3=reshape(qOldWell3,[277*129,1]);
flateffectsPooledOld=reshape(effectsPooledOld,[277*129,1]);
flatqOldPooledOld=reshape(qOldPooledOld,[277*129,1]);

scatter(flateffectsGlobal(flatqGlobal<0.1),flateffectsOldWell1(flatqGlobal<0.1),'.')
title('Effects significant in new data')
xlabel('Pooled new data')
ylabel('Old well 1')

scatter(flateffectsGlobal(flatqGlobal<0.1),flateffectsOldWell3(flatqGlobal<0.1),'.')
title('Effects significant in new data')
xlabel('Pooled new data')
ylabel('Old well 3')

scatter(flateffectsGlobal(flatqGlobal<0.05 & flatqOldPooledOld<0.05),flateffectsPooledOld(flatqGlobal<0.05 & flatqOldPooledOld<0.05),'.')
title('Effects significant in both')
xlabel('Pooled new data')
ylabel('Pooled old data')

% looking at self targeting
selfIDs=zeros(277, 129);
for i=1:129
    selfIDs(:,i)=strcmp(genes,guides(i));
end
flatselfids=reshape(selfIDs,[277*129,1]);

flatsplit1_LFC=reshape(split1LFCs,[232*128,1]);
flatsplit2_LFC=reshape(split2LFCs,[232*128,1]);
% selfIDs=zeros(232, 128);
% for i=1:128
%     selfIDs(:,i)=strcmp(genes,guides(i));
% end
% flatselfids=reshape(selfIDs,[232*128,1]);

% scatter(flatsplit1_LFC,flatsplit2_LFC,'.b')
% hold on
% scatter(flatsplit1_LFC(flatselfids>0),flatsplit2_LFC(flatselfids>0),'.r')
% set(gcf,'color','w');
% title('All vs self effects')
% xlabel('1st split')
% ylabel('2nd split')
scatter(flateffectsWell1(flatselfids>0),flateffectsGlobal(flatselfids>0),'.r')

title('Self effects')
xlabel('Well 1')
ylabel('Pooled')
scatter(flateffectsWell1,flateffectsWell2,'.')

for i=1:128
    if sum(selfIDs(:,i))>0
figure,
scatter(split1LFCs(:,i),split2LFCs(:,i),'.b')
hold on
scatter(split1LFCs(selfIDs(:,i)>0,i),split2LFCs(selfIDs(:,i)>0,i),'.r')
set(gcf,'color','w');
title(strcat('All vs self effects, guide',guides(i)))
    end
end
%%
split1LFCs1=split1LFCs1(:,1:127);
split2LFCs1=split2LFCs1(:,1:127);
split1qvals=split1qvals(:,1:127);
split2qvals=split2qvals(:,1:127);
flatsplit1_LFC=reshape(split1LFCs1,[277*127,1]);
flatsplit2_LFC=reshape(split2LFCs1,[277*127,1]);
flatsplit1_q=reshape(split1qvals,[277*127,1]);
flatsplit2_q=reshape(split2qvals,[277*127,1]);
selfIDs=zeros(277, 127);
for i=1:127
    selfIDs(:,i)=strcmp(genes,guides(i));
end
flatselfids=reshape(selfIDs,[277*127,1]);
subplot(2,3,1)

scatter(flatsplit1_LFC,flatsplit2_LFC,'.')
title('All effects')
xlabel('1st split')
ylabel('2nd split')
subplot(2,3,2)

scatter(flatsplit1_LFC(flatsplit1_q<0.1 | flatsplit2_q<0.1),flatsplit2_LFC(flatsplit1_q<0.1 | flatsplit2_q<0.1),'.')
title('Effects significant in at least one half')
xlabel('1st split')
ylabel('2nd split')
subplot(2,3,3)

scatter(flatsplit1_LFC(flatsplit1_q<0.1 & flatsplit2_q<0.1),flatsplit2_LFC(flatsplit1_q<0.1 & flatsplit2_q<0.1),'.')
title('Effects significant in both halves')
xlabel('1st split')
ylabel('2nd split')

subplot(2,3,4)
scatter(flatsplit1_LFC,flatsplit2_LFC,'.')
hold on
scatter(flatsplit1_LFC(flatselfids>0),flatsplit2_LFC(flatselfids>0),'.r')
title('All vs self effects, significant in atl east one half')
xlabel('1st split')
ylabel('2nd split')

subplot(2,3,5)
scatter(flatsplit1_LFC(flatsplit1_q<0.1 | flatsplit2_q<0.1),flatsplit2_LFC(flatsplit1_q<0.1 | flatsplit2_q<0.1),'.')
hold on
scatter(flatsplit1_LFC((flatsplit1_q<0.1 | flatsplit2_q<0.1) & flatselfids>0),flatsplit2_LFC((flatsplit1_q<0.1 | flatsplit2_q<0.1) &  flatselfids>0),'.r')
title('All vs self effects, significant in atl east one half')
xlabel('1st split')
ylabel('2nd split')

set(gcf,'color','w');
%%
sigID=qGlobal<0.1;
sigguides=sum(sigID,1)>1;
sigguides=sigguides(1:end-1);
siggenes=sum(sigID,2)>1;
significanteffects=effectsGlobal(siggenes,sigguides)';
heatmap(genes(siggenes),guides(sigguides),significanteffects,'Colormap',redblue(256))
ylabel("Perturbations")
xlabel("Genes")
title("Heatmap of effects")
caxis([-0.3 0.3])
set(gcf,'color','w');
%% clustering
keptguides=guides(sigguides);
keptgenes=genes(siggenes);
% keptguides=keptguides(1:end-1);
% significanteffects=significanteffects(1:end-1,:);

opts = statset('Display','final');
[idxguides,C] = kmeans(significanteffects,4,'Distance','cosine','Replicates',1000,'MaxIter',1000,'Options',opts);
sorted=significanteffects(idxguides==1,:);
sortedguides=keptguides(idxguides==1);
for i=2:max(idxguides)
sorted=[sorted;significanteffects(idxguides==i,:)];
sortedguides=[sortedguides,keptguides(idxguides==i)];
end

[idxgenes,C] = kmeans(significanteffects',5,'Distance','cosine','Replicates',1000,'MaxIter',1000,'Options',opts);
sortedagain=sorted(:,idxgenes==1);
sortedgenes=keptgenes(idxgenes==1);
for i=2:max(idxgenes)
sortedagain=[sortedagain,sorted(:,idxgenes==i)];
sortedgenes=[sortedgenes;keptgenes(idxgenes==i)];
end
figure,
heatmap(sortedgenes,sortedguides,sortedagain,'Colormap',redblue(256))
ylabel("Perturbations")
xlabel("Genes")
title("Heatmap of effects")
caxis([-0.2 0.2])
set(gcf,'color','w');
writematrix(sortedagain,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtoperturbationMatching\sortedeffects0424.csv')
writematrix(sortedguides,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtoperturbationMatching\sortedguides0424.csv')
writematrix(sortedgenes,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtoperturbationMatching\sortedgenes0424.csv')
writematrix(idxgenes,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtoperturbationMatching\idxgenes0424.csv')
writematrix(idxguides,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtoperturbationMatching\idxguides0424.csv')

%% with calcium
perturbcalcium=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\clusteredCalciumperturbcells.csv'));
perturbations=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\perturbationsInCalcium.csv'));

perturbcalciumlessclusters=perturbcalcium;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==7,3)=5;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==6,3)=2;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==8,3)=6;

guidespercluster=zeros(6,128);
for cluster=1:max(perturbcalciumlessclusters(:,3))
    theseperturbations=perturbations(perturbcalciumlessclusters(:,3)==cluster,:);
    mylist=theseperturbations(find(theseperturbations));
    for i=1:128
        guidespercluster(cluster,i)=sum(mylist==i);
    end
end
cellspercluster=[]
for i=1:max(perturbcalciumlessclusters(:,3))
cellspercluster=[cellspercluster;sum(perturbcalciumlessclusters(:,3)==i)];
end
expectedfreq=zeros(6,128);
for j=1:128
    expectedfreq(:,j)=sum(guidespercluster(:,j))*cellspercluster'/sum(cellspercluster);
end
FC=log(guidespercluster./expectedfreq);
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

M=size(perturbations,1);%sample size
K=sum(guidespercluster(1,:));%number of good items in pop
N=sum(guidespercluster(:,1));%number drawn
pval=1-hygecdf(guidespercluster(:,1)-1,M,K,N);
p=zeros(max(perturbcalciumlessclusters(:,3)),size(guidespercluster,2));
for i=1:size(guidespercluster,2)
    for j=1:size(guidespercluster,1)
        K=sum(guidespercluster(j,:));
        N=sum(guidespercluster(:,i));
        pval=1-hygecdf(guidespercluster(j,i)-1,M,K,N);
        p(j,i)=pval;
    end
end
guides=guides(1:128);
heatmap(guides(sum(p<0.1,1)>0),clusternames,FC(:,sum(p<0.1,1)>0),'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])
%% sort calcium enrichment
keptFC=FC(:,sigguides(1:128));
keptp=p(:,sigguides(1:128));
sortedFC=keptFC(:,idxguides==1);
sortedp=keptp(:,idxguides==1);

for i=2:max(idxguides)
sortedFC=[sortedFC,keptFC(:,idxguides==i)];
sortedp=[sortedp,keptFC(:,idxguides==i)];
end
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
heatmap(clusternames,sortedguides,sortedFC','Colormap',redblue(256))
xlabel("Cluster")
ylabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

sigOnly=sortedFC;
sigOnly(sortedp>0.1)=NaN;
heatmap(clusternames,sortedguides,sigOnly','Colormap',redblue(256))
xlabel("Cluster")
ylabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

%% calcium + genes
calciumclusterallidxitcalciumZscoreONEatAtime=table2array(readtable('\\helium\broad_clearylab\Users\Loic\pooledrevision\singleclusters\calciumclusterall_idx_it_calciumZscoreONEatAtime.csv'));
pooledmerfish=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\pooledmerfishmatchedtocalcium.csv'));

pooledmerfish=pooledmerfish(:,3:end);
totalcounts=sum(pooledmerfish,2);
badcells=totalcounts<200;
badcells=badcells|(totalcounts>1900);
pooledmerfish=pooledmerfish(~badcells,:);
calcium=calciumclusterallidxitcalciumZscoreONEatAtime(~badcells,:);
pooledmerfish=pooledmerfish(:,1:485);
zscored=pooledmerfish;
zscored(calcium(:,1)==21,:)=(pooledmerfish(calcium(:,1)==21,:)-mean(pooledmerfish(calcium(:,1)==21,:),1))./std(pooledmerfish(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(pooledmerfish(calcium(:,1)==22,:)-mean(pooledmerfish(calcium(:,1)==22,:),1))./std(pooledmerfish(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(pooledmerfish(calcium(:,1)==31,:)-mean(pooledmerfish(calcium(:,1)==31,:),1))./std(pooledmerfish(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(pooledmerfish(calcium(:,1)==32,:)-mean(pooledmerfish(calcium(:,1)==32,:),1))./std(pooledmerfish(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(pooledmerfish(calcium(:,1)==4,:)-mean(pooledmerfish(calcium(:,1)==4,:),1))./std(pooledmerfish(calcium(:,1)==4,:),1);


for j=1:max(merfishedcalcium(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(calcium(:,3)==j,:));
    controlsignature(j,:)=median(zscored(calcium(:,3)~=j,:));
end
% figure, heatmap(Zscoredgenesignature,'Colormap',redblue(256))

%build geneNames by keeping the good genes from FR, and reslice
%Zscoredgenesignature with that. might need to get those signatures on data
%normed to total counts
for i=1:485
    if sum(strcmp(genes,allgenes(i))==1)
        goodGenesIdx(i)=1;
    end
end
figure, heatmap(Zscoredgenesignature(:,goodGenesIdx),'Colormap',redblue(256))
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
figure, heatmap(genes,clusternames,Zscoredgenesignature(:,goodGenesIdx>0),'Colormap',redblue(256))

S=Zscoredgenesignature(:,goodGenesIdx>0);
sortedZscore=S(:,idxgenes==1);
for i=2:max(idxgenes)
sortedZscore=[sortedZscore,S(:,idxgenes==i)];
end
figure, heatmap(sortedgenes,clusternames,sortedZscore,'Colormap',redblue(256))
set(gcf,'color','w');
caxis([-0.3 0.3])


%% calcium + genes normed to total counts
calciumclusterallidxitcalciumZscoreONEatAtime=table2array(readtable('\\helium\broad_clearylab\Users\Loic\pooledrevision\singleclusters\calciumclusterall_idx_it_calciumZscoreONEatAtime.csv'));
pooledmerfish=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\pooledmerfishmatchedtocalcium.csv'));

pooledmerfish=pooledmerfish(:,3:end);
totalcounts=sum(pooledmerfish,2);
badcells=totalcounts<200;
badcells=badcells|(totalcounts>1900);
pooledmerfish=pooledmerfish(~badcells,:);
calcium=calciumclusterallidxitcalciumZscoreONEatAtime(~badcells,:);
pooledmerfish=pooledmerfish(:,1:485);
pooledmerfish=pooledmerfish./sum(pooledmerfish,2);
zscored=pooledmerfish;
zscored(calcium(:,1)==21,:)=(pooledmerfish(calcium(:,1)==21,:)-mean(pooledmerfish(calcium(:,1)==21,:),1))./std(pooledmerfish(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(pooledmerfish(calcium(:,1)==22,:)-mean(pooledmerfish(calcium(:,1)==22,:),1))./std(pooledmerfish(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(pooledmerfish(calcium(:,1)==31,:)-mean(pooledmerfish(calcium(:,1)==31,:),1))./std(pooledmerfish(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(pooledmerfish(calcium(:,1)==32,:)-mean(pooledmerfish(calcium(:,1)==32,:),1))./std(pooledmerfish(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(pooledmerfish(calcium(:,1)==4,:)-mean(pooledmerfish(calcium(:,1)==4,:),1))./std(pooledmerfish(calcium(:,1)==4,:),1);


for j=1:max(merfishedcalcium(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(calcium(:,3)==j,:));
    controlsignature(j,:)=median(zscored(calcium(:,3)~=j,:));
end
% figure, heatmap(Zscoredgenesignature,'Colormap',redblue(256))

%build geneNames by keeping the good genes from FR, and reslice
%Zscoredgenesignature with that. might need to get those signatures on data
%normed to total counts
for i=1:485
    if sum(strcmp(genes,allgenes(i))==1)
        goodGenesIdx(i)=1;
    end
end
% figure, heatmap(Zscoredgenesignature(:,goodGenesIdx),'Colormap',redblue(256))
% clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
% figure, heatmap(genes,clusternames,Zscoredgenesignature(:,goodGenesIdx>0),'Colormap',redblue(256))

S=Zscoredgenesignature(:,goodGenesIdx>0);
sortedZscore=S(:,idxgenes==1);
for i=2:max(idxgenes)
sortedZscore=[sortedZscore,S(:,idxgenes==i)];
end
figure, heatmap(sortedgenes,clusternames,sortedZscore,'Colormap',redblue(256))
set(gcf,'color','w');
caxis([-0.3 0.3])


%% repeat all to plot genes and guides with 2 signigficant effects at least
sigID=qGlobal<0.1;
sigguides=sum(sigID,1)>2;
sigguides=sigguides(1:end-1);
siggenes=sum(sigID,2)>2;
significanteffects=effectsGlobal(siggenes,sigguides)';
heatmap(genes(siggenes),guides(sigguides),significanteffects,'Colormap',redblue(256))
ylabel("Perturbations")
xlabel("Genes")
title("Heatmap of effects")
caxis([-0.3 0.3])
set(gcf,'color','w');
%% clustering
keptguides=guides(sigguides);
keptgenes=genes(siggenes);
keptguides=keptguides(1:end-1);
significanteffects=significanteffects(1:end-1,:);

opts = statset('Display','final');
[idxguides,C] = kmeans(significanteffects,3,'Distance','correlation','Replicates',60,'MaxIter',1000,'Options',opts);
sorted=significanteffects(idxguides==1,:);
sortedguides=keptguides(idxguides==1);
for i=2:max(idxguides)
sorted=[sorted;significanteffects(idxguides==i,:)];
sortedguides=[sortedguides,keptguides(idxguides==i)];
end

[idxgenes,C] = kmeans(significanteffects',3,'Distance','correlation','Replicates',60,'MaxIter',1000,'Options',opts);
sortedagain=sorted(:,idxgenes==1);
sortedgenes=keptgenes(idxgenes==1);
for i=2:max(idxgenes)
sortedagain=[sortedagain,sorted(:,idxgenes==i)];
sortedgenes=[sortedgenes;keptgenes(idxgenes==i)];
end
figure,
heatmap(sortedgenes,sortedguides,sortedagain,'Colormap',redblue(256))
ylabel("Perturbations")
xlabel("Genes")
title("Heatmap of effects")
caxis([-0.3 0.3])
set(gcf,'color','w');
%% with calcium
perturbcalcium=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\clusteredCalciumperturbcells.csv'));
perturbations=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\perturbationsInCalcium.csv'));

perturbcalciumlessclusters=perturbcalcium;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==7,3)=5;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==6,3)=2;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==8,3)=6;

guidespercluster=zeros(6,128);
for cluster=1:max(perturbcalciumlessclusters(:,3))
    theseperturbations=perturbations(perturbcalciumlessclusters(:,3)==cluster,:);
    mylist=theseperturbations(find(theseperturbations));
    for i=1:128
        guidespercluster(cluster,i)=sum(mylist==i);
    end
end
cellspercluster=[]
for i=1:max(perturbcalciumlessclusters(:,3))
cellspercluster=[cellspercluster;sum(perturbcalciumlessclusters(:,3)==i)];
end
expectedfreq=zeros(6,128);
for j=1:128
    expectedfreq(:,j)=sum(guidespercluster(:,j))*cellspercluster'/sum(cellspercluster);
end
FC=log(guidespercluster./expectedfreq);
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

M=size(perturbations,1);%sample size
K=sum(guidespercluster(1,:));%number of good items in pop
N=sum(guidespercluster(:,1));%number drawn
pval=1-hygecdf(guidespercluster(:,1)-1,M,K,N);
p=zeros(max(perturbcalciumlessclusters(:,3)),size(guidespercluster,2));
for i=1:size(guidespercluster,2)
    for j=1:size(guidespercluster,1)
        K=sum(guidespercluster(j,:));
        N=sum(guidespercluster(:,i));
        pval=1-hygecdf(guidespercluster(j,i)-1,M,K,N);
        p(j,i)=pval;
    end
end
guides=guides(1:128);
heatmap(guides(sum(p<0.1,1)>0),clusternames,FC(:,sum(p<0.1,1)>0),'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])
%% sort calcium enrichment
keptFC=FC(:,sigguides(1:128));
keptp=p(:,sigguides(1:128));
sortedFC=keptFC(:,idxguides==1);
sortedp=keptp(:,idxguides==1);

for i=2:max(idxguides)
sortedFC=[sortedFC,keptFC(:,idxguides==i)];
sortedp=[sortedp,keptFC(:,idxguides==i)];
end
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
heatmap(clusternames,sortedguides,sortedFC','Colormap',redblue(256))
xlabel("Cluster")
ylabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

sigOnly=sortedFC;
sigOnly(sortedp>0.1)=NaN;
figure,heatmap(clusternames,sortedguides,sigOnly','Colormap',redblue(256))
xlabel("Cluster")
ylabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1 1])

%% calcium + genes
calciumclusterallidxitcalciumZscoreONEatAtime=table2array(readtable('\\helium\broad_clearylab\Users\Loic\pooledrevision\singleclusters\calciumclusterall_idx_it_calciumZscoreONEatAtime.csv'));
pooledmerfish=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\pooledmerfishmatchedtocalcium.csv'));

pooledmerfish=pooledmerfish(:,3:end);
totalcounts=sum(pooledmerfish,2);
badcells=totalcounts<200;
badcells=badcells|(totalcounts>1900);
pooledmerfish=pooledmerfish(~badcells,:);
calcium=calciumclusterallidxitcalciumZscoreONEatAtime(~badcells,:);
pooledmerfish=pooledmerfish(:,1:485);
zscored=pooledmerfish;
zscored(calcium(:,1)==21,:)=(pooledmerfish(calcium(:,1)==21,:)-mean(pooledmerfish(calcium(:,1)==21,:),1))./std(pooledmerfish(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(pooledmerfish(calcium(:,1)==22,:)-mean(pooledmerfish(calcium(:,1)==22,:),1))./std(pooledmerfish(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(pooledmerfish(calcium(:,1)==31,:)-mean(pooledmerfish(calcium(:,1)==31,:),1))./std(pooledmerfish(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(pooledmerfish(calcium(:,1)==32,:)-mean(pooledmerfish(calcium(:,1)==32,:),1))./std(pooledmerfish(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(pooledmerfish(calcium(:,1)==4,:)-mean(pooledmerfish(calcium(:,1)==4,:),1))./std(pooledmerfish(calcium(:,1)==4,:),1);


for j=1:max(merfishedcalcium(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(calcium(:,3)==j,:));
    controlsignature(j,:)=median(zscored(calcium(:,3)~=j,:));
end
% figure, heatmap(Zscoredgenesignature,'Colormap',redblue(256))

%build geneNames by keeping the good genes from FR, and reslice
%Zscoredgenesignature with that. might need to get those signatures on data
%normed to total counts
for i=1:485
    if sum(strcmp(genes,allgenes(i))==1)
        goodGenesIdx(i)=1;
    end
end
figure, heatmap(Zscoredgenesignature(:,goodGenesIdx),'Colormap',redblue(256))
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
figure, heatmap(genes,clusternames,Zscoredgenesignature(:,goodGenesIdx>0),'Colormap',redblue(256))

S=Zscoredgenesignature(:,goodGenesIdx>0);
sortedZscore=S(:,idxgenes==1);
for i=2:max(idxgenes)
sortedZscore=[sortedZscore,S(:,idxgenes==i)];
end
figure, heatmap(sortedgenes,clusternames,sortedZscore,'Colormap',redblue(256))
set(gcf,'color','w');
caxis([-0.3 0.3])


%% calcium + genes normed to total counts
calciumclusterallidxitcalciumZscoreONEatAtime=table2array(readtable('\\helium\broad_clearylab\Users\Loic\pooledrevision\singleclusters\calciumclusterall_idx_it_calciumZscoreONEatAtime.csv'));
pooledmerfish=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\pooledmerfishmatchedtocalcium.csv'));

pooledmerfish=pooledmerfish(:,3:end);
totalcounts=sum(pooledmerfish,2);
badcells=totalcounts<200;
badcells=badcells|(totalcounts>1900);
pooledmerfish=pooledmerfish(~badcells,:);
calcium=calciumclusterallidxitcalciumZscoreONEatAtime(~badcells,:);
pooledmerfish=pooledmerfish(:,1:485);
pooledmerfish=pooledmerfish./sum(pooledmerfish,2);
zscored=pooledmerfish;
zscored(calcium(:,1)==21,:)=(pooledmerfish(calcium(:,1)==21,:)-mean(pooledmerfish(calcium(:,1)==21,:),1))./std(pooledmerfish(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(pooledmerfish(calcium(:,1)==22,:)-mean(pooledmerfish(calcium(:,1)==22,:),1))./std(pooledmerfish(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(pooledmerfish(calcium(:,1)==31,:)-mean(pooledmerfish(calcium(:,1)==31,:),1))./std(pooledmerfish(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(pooledmerfish(calcium(:,1)==32,:)-mean(pooledmerfish(calcium(:,1)==32,:),1))./std(pooledmerfish(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(pooledmerfish(calcium(:,1)==4,:)-mean(pooledmerfish(calcium(:,1)==4,:),1))./std(pooledmerfish(calcium(:,1)==4,:),1);


for j=1:max(merfishedcalcium(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(calcium(:,3)==j,:));
    controlsignature(j,:)=median(zscored(calcium(:,3)~=j,:));
end
% figure, heatmap(Zscoredgenesignature,'Colormap',redblue(256))

%build geneNames by keeping the good genes from FR, and reslice
%Zscoredgenesignature with that. might need to get those signatures on data
%normed to total counts
for i=1:485
    if sum(strcmp(genes,allgenes(i))==1)
        goodGenesIdx(i)=1;
    end
end
% figure, heatmap(Zscoredgenesignature(:,goodGenesIdx),'Colormap',redblue(256))
% clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
% figure, heatmap(genes,clusternames,Zscoredgenesignature(:,goodGenesIdx>0),'Colormap',redblue(256))

S=Zscoredgenesignature(:,goodGenesIdx>0);
sortedZscore=S(:,idxgenes==1);
for i=2:max(idxgenes)
sortedZscore=[sortedZscore,S(:,idxgenes==i)];
end
figure, heatmap(sortedgenes,clusternames,sortedZscore,'Colormap',redblue(256))
set(gcf,'color','w');
caxis([-0.3 0.3])

