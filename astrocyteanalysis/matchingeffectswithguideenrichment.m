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
% heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
% ylabel("Cluster")
% xlabel("Perturbation")
% title("Perturbation enrichment per cluster")
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
        pval=min(1-hygecdf(guidespercluster(j,i)-1,M,K,N),hygecdf(guidespercluster(j,i)-1,M,K,N));
        p(j,i)=pval;
    end
end
globalguideFC=FC;
globalguidep=p;
guides=guides(1:128);
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
figure,
heatmap(clusternames,guides,globalguideFC','Colormap',redblue(256))

globalSIGguideFC=globalguideFC(:,sigguides);
globalSIGguidep=globalguidep(:,sigguides);
cropguides=guides(sigguides);

sortedglobalguideFC=globalSIGguideFC(:,idxguides==1);
sortedglobalguidep=globalSIGguidep(:,idxguides==1);
for i=2:max(idxguides)
sortedglobalguideFC=[sortedglobalguideFC,globalSIGguideFC(:,idxguides==i)];
sortedglobalguidep=[sortedglobalguidep,globalSIGguidep(:,idxguides==i)];
end
sortedglobalguideFCNANs=sortedglobalguideFC;
sortedglobalguideFCNANs(sortedglobalguidep>0.1)=NaN;

figure,
heatmap(clusternames,sortedguides,sortedglobalguideFC','Colormap',redblue(256),'CellLabelColor','none')
xlabel("Calcium phenotype")
ylabel("Gene Knock down")
title("Guide enrichment")
caxis([-0.35 0.35])
set(gcf,'color','w');


figure,
heatmap(clusternames,sortedguides,sortedglobalguideFCNANs','Colormap',redblue(256),'CellLabelColor','none')
xlabel("Calcium phenotype")
ylabel("Gene Knock down")
title("Guide enrichment")
caxis([-0.35 0.35])
set(gcf,'color','w');
globalguideFC(:,strcmp(guides,"TRAF7"))

%% pooling guides per cluster
sigguidespercluster=guidespercluster(:,sigguides);
for i=1:max(idxguides)
    pooledguidespercluster(:,i)=sum(sigguidespercluster(:,idxguides==i),2);
end

cellspercluster=sum(pooledguidespercluster,2);
expectedfreq=zeros(6,max(idxguides));
for j=1:max(idxguides)
    expectedfreq(:,j)=sum(pooledguidespercluster(:,j))*cellspercluster'/sum(cellspercluster);
end
FC=log(pooledguidespercluster./expectedfreq);
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
% heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
% ylabel("Cluster")
% xlabel("Perturbation")
% title("Perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

M=sum(sum(pooledguidespercluster));%sample size
K=sum(pooledguidespercluster(1,:));%number of good items in pop
N=sum(pooledguidespercluster(:,1));%number drawn
pval=1-hygecdf(pooledguidespercluster(:,1)-1,M,K,N);
p=zeros(max(perturbcalciumlessclusters(:,3)),size(pooledguidespercluster,2));
for i=1:size(pooledguidespercluster,2)
    for j=1:size(pooledguidespercluster,1)
        K=sum(pooledguidespercluster(j,:));
        N=sum(pooledguidespercluster(:,i));
        pval=min(1-hygecdf(pooledguidespercluster(j,i)-1,M,K,N),hygecdf(pooledguidespercluster(j,i)-1,M,K,N));
        p(j,i)=pval;
    end
end
SigguideFC=FC;
Sigguidep=p;

figure,
heatmap(clusternames,["Cluster 1","Cluster 2","Cluster 3","Cluster 4"],SigguideFC','Colormap',redblue(256),'CellLabelColor','none')
xlabel("Calcium phenotype")
ylabel("Perturbation cluster")
title("Guide enrichment")
caxis([-0.35 0.35])
set(gcf,'color','w');
