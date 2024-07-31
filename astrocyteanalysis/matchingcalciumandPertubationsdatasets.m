perturbedcalcium=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==1,:);
perturbedcalcium=[perturbedcalcium;calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==2,:)];
perturbedcalcium=[perturbedcalcium;calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==4,:)];

well1=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==21,:);
Zwell1=goodcellsallindex(goodcellsallindex(:,1)==1,:);
perturbcalcium=[];
for i=1:size(well1,1)
    if ismember(well1(i,2),Zwell1(:,2));
        perturbcalcium=[perturbcalcium;well1(i,:)];
    end
end

well2=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==22,:);
Zwell2=goodcellsallindex(goodcellsallindex(:,1)==2,:);
for i=1:size(well2,1)
    if ismember(well2(i,2),Zwell2(:,2));
        perturbcalcium=[perturbcalcium;well2(i,:)];
    end
end


well4=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==4,:);
Zwell4=goodcellsallindex(goodcellsallindex(:,1)==3,:);
for i=1:size(well4,1)
    if ismember(well4(i,2),Zwell4(:,2));
        perturbcalcium=[perturbcalcium;well4(i,:)];
    end
end
writematrix(perturbcalcium,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\clusteredCalciumperturbcells.csv')
% goodcellsallzombie=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\goodcellsallzombie.csv'));
zombie=goodcellsallzombie';
zombie(:,128)=zombie(:,128)+zombie(:,129);
zombie=zombie(:,1:128);
zombie=zombie>0;
zombie=[1:128].*zombie;
perturbations=[];
p1idx=goodcellsallindex(goodcellsallindex(:,1)==1,:);
p1=zombie(goodcellsallindex(:,1)==1,:);
for i=1:size(p1,1)
    if ismember(p1idx(i,2),well1(:,2));
    perturbations=[perturbations;p1(i,:)];
    end
end

p2idx=goodcellsallindex(goodcellsallindex(:,1)==2,:);
p2=zombie(goodcellsallindex(:,1)==2,:);
for i=1:size(p2,1)
    if ismember(p2idx(i,2),well2(:,2));
    perturbations=[perturbations;p2(i,:)];
    end
end

   
p4idx=goodcellsallindex(goodcellsallindex(:,1)==3,:);
p4=zombie(goodcellsallindex(:,1)==3,:);
for i=1:size(p4,1)
    if ismember(p4idx(i,2),well4(:,2));
    perturbations=[perturbations;p4(i,:)];
    end
end
writematrix(perturbations,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\perturbationsInCalcium.csv')

%% cluster enrichment in perturbations
guidespercluster=zeros(8,128);
for cluster=1:max(perturbcalcium(:,3))
    theseperturbations=perturbations(perturbcalcium(:,3)==cluster,:);
    mylist=theseperturbations(find(theseperturbations));
    for i=1:128
        guidespercluster(cluster,i)=sum(mylist==i);
    end
end
cellspercluster=[]
for i=1:max(perturbcalcium(:,3))
cellspercluster=[cellspercluster;sum(perturbcalcium(:,3)==i)];
end
expectedfreq=zeros(8,128);
for j=1:128
    expectedfreq(:,j)=sum(guidespercluster(:,j))*cellspercluster'/sum(cellspercluster);
end
FC=log(guidespercluster./expectedfreq);
clusternames=["large peak","inactive","large early transient","small peak","step","oscillating","step2","delayed"];
heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

%% build pvalue

M=size(perturbations,1);%sample size
K=sum(guidespercluster(1,:));%number of good items in pop
N=sum(guidespercluster(:,1));%number drawn
pval=1-hygecdf(guidespercluster(:,1)-1,M,K,N);
p=zeros(max(perturbcalcium(:,3)),size(guidespercluster,2));
for i=1:size(guidespercluster,2)
    for j=1:size(guidespercluster,1)
        K=sum(guidespercluster(j,:));
        N=sum(guidespercluster(:,i));
        pval=1-hygecdf(guidespercluster(j,i)-1,M,K,N);
        p(j,i)=pval;
    end
end
guides=guides(1:128);
heatmap(guides(sum(p<0.2,1)>0),clusternames,FC(:,sum(p<0.2,1)>0),'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

%% lessclusters
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
figure, heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

%% build pvalue

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
guides=guides(1:128);
heatmap(guides(sum(p<0.1,1)>0),clusternames,FC(:,sum(p<0.1,1)>0),'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("Significant perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

heatmap(guides(sum(p<0.1,1)>0),clusternames,p(:,sum(p<0.1,1)>0),'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("p values per cluster")
set(gcf,'color','w');caxis([0 1])

heatmap(guides,clusternames,p,'Colormap',redblue(256))
ylabel("Cluster")
xlabel("Perturbation")
title("p values per cluster")
set(gcf,'color','w');caxis([0 1])