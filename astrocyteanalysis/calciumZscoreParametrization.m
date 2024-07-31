allcalciumidx=table2array(readtable('/broad/clearylab/Users/Loic/pooledrevision/allcalciumidx_0329.csv'));
allcalciumtracesnormed=table2array(readtable('/broad/clearylab/Users/Loic/pooledrevision/allcalciumtracesnormed_0329.csv'));
calciumsummary=table2array(readtable('/broad/clearylab/Users/Loic/pooledrevision/calciumsummary_0329.csv'));

zscored=[];
for sample=1:5
        if sample==1
            thissampleID=21;
        elseif sample==2
            thissampleID=22;
        elseif sample==3
            thissampleID=31;          
        elseif sample==4
            thissampleID=32;           
        elseif sample==5
            thissampleID=4;
        end
    zscored=[zscored;(calciumsummary(allcalciumidx(:,1)==thissampleID,:)-mean(calciumsummary(allcalciumidx(:,1)==thissampleID,:),1))./std(calciumsummary(allcalciumidx(:,1)==thissampleID,:),1)];
end
opts = statset('Display','final');
[idx,C] = kmeans(zscored,8,'Distance','correlation','Replicates',1000,'MaxIter',1000,'Options',opts);
calciumclusters=[allcalciumidx,idx];
writematrix(calciumclusters,strcat('/broad/clearylab/Users/Loic/pooledrevision/singleclusters/calciumclusterall_idx_it_calciumZscoreONEatAtime.csv'))
writematrix(C,strcat('/broad/clearylab/Users/Loic/pooledrevision/singleclusters/calciumclustersall_C_itcalciumZscoreONEatAtime.csv'))
calcium=allcalciumtracesnormed;
time=[1:200];
for sample=1:5
        if sample==1
            thissampleID=21;
        elseif sample==2
            thissampleID=22;
        elseif sample==3
            thissampleID=31;          
        elseif sample==4
            thissampleID=32;           
        elseif sample==5
            thissampleID=4;
        end
        fig1 = figure;
filename=strcat('/broad/clearylab/Users/Loic/pooledrevision/singleclusters/traceZscoredONEatAtimesample_',num2str(sample),'.png');

for i=1:8
    subtraces=calcium(idx==i & allcalciumidx(:,1)==thissampleID,:);
    for j=10:20
    subplot(2,4,i)
    plot(time,subtraces(j,:)+(j-10)*10,'LineWidth',2)
%     title(clusternames(i))
    ylim([0 350])
    hold on
    end
end
set(gcf,'color','w');
saveas(fig1, filename);
end
