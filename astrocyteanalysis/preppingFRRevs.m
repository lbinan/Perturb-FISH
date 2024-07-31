zombieS2=[zombiecounttableNEWR0;zombiecounttableNEWR1];
pooledtable=zeros(31879,129);
for i=1:129
    pooledtable(:,i)=zombieS2(:,(2*i-1))+zombieS2(:,(2*i));
end
writematrix(pooledtable,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\pooledzombieS2.csv')

merfishS2=[merfishcounttableregion0rightIndex;merfishcounttableregion1rightIndex];
s2index=[[ones(size(merfishcounttableregion0rightIndex,1),1);2*ones(size(merfishcounttableregion1rightIndex,1),1)],merfishS2(:,1)];
writematrix(merfishS2(:,3:end),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\rawmerfishS2.csv')
writematrix(s2index,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\rawindecesS2.csv')

allzombie=[pooledtable;pooledzombieS4];
allindex=[s2index;[3*ones(size(pooledzombieS4,1),1),merfishcounttablerightIndex(:,1)]];
writematrix(allindex,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\allindex.csv')
writematrix(allzombie,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\allzombie.csv')
allmerfish=[merfishS2;merfishcounttablerightIndex];
writematrix(allmerfish(:,3:end),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\allmerfish.csv')

totalcounts=sum(allmerfish(:,3:end),2);
badcells=totalcounts<200;
badcells=badcells|totalcounts>1900;
badgennes=sum(allmerfish(:,3:end),1)/size(allmerfish,1)<0.5;
writematrix(allmerfish(~badcells,3:end),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodcellsallmerfish.csv')
writematrix(allzombie(~badcells,:),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodcellsallzombie.csv')
writematrix(allindex(~badcells,:),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodcellsallindex.csv')
badgennes=[1,1,badgennes];
writematrix(allmerfish(~badcells,~badgennes),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodcellgoodgenessallmerfish.csv')
goodgenenames=genenames(~(badgennes(3:end)>0));
writematrix(goodgenenames','\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenenames.csv')
writematrix(genenames','\\helium\broad_clearylab\Users\Loic\AsdRevstables\genenames.csv')
goodcellsmerfish=allmerfish(~badcells,3:end);
goodcellsallzombie=allzombie(~badcells,:);
goodcellsallindex=allindex(~badcells,:);
perturbedcells=sum(goodcellsallzombie,2)>0&sum(goodcellsallzombie,2)<3;

writematrix(goodcellsmerfish(perturbedcells,:),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodperturbedcellsallmerfish.csv')
writematrix(goodcellsallzombie(perturbedcells,:)','\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodcellsallzombie.csv')
writematrix(goodcellsallindex(perturbedcells,:),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodcellsallindex.csv')

writematrix(goodcellsmerfish(perturbedcells,:),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodperturbedcellsAllGenesmerfish.csv')
writematrix(goodcellsmerfish(perturbedcells,~badgennes(3:end)),'\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodperturbedcellsGoodGenesmerfish.csv')

well1goodcellszombie=goodcellsallzombie(perturbedcells & goodcellsallindex(:,1)==1,:)';
well2goodcellszombie=goodcellsallzombie(perturbedcells & goodcellsallindex(:,1)==2,:)';
well3goodcellszombie=goodcellsallzombie(perturbedcells & goodcellsallindex(:,1)==3,:)';

well1goodcellsmerfish=goodcellsmerfish(perturbedcells & goodcellsallindex(:,1)==1,~badgennes(3:end));
well2goodcellsmerfish=goodcellsmerfish(perturbedcells & goodcellsallindex(:,1)==2,~badgennes(3:end));
well3goodcellsmerfish=goodcellsmerfish(perturbedcells & goodcellsallindex(:,1)==3,~badgennes(3:end));

writematrix(well1goodcellsmerfish,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\well1goodcellsgoodgenesmerfish.csv')
writematrix(well1goodcellszombie,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\well1goodcellsgoodgeneszombie.csv')
writematrix(well2goodcellsmerfish,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\well2goodcellsgoodgenesmerfish.csv')
writematrix(well2goodcellszombie,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\well2goodcellsgoodgeneszombie.csv')
writematrix(well3goodcellsmerfish,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\well3goodcellsgoodgenesmerfish.csv')
writematrix(well3goodcellszombie,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\well3goodcellsgoodgeneszombie.csv')
