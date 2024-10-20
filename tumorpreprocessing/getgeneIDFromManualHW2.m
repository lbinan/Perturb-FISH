zombiecodebookHD4reordered=readtable(fullfile(matlabtempfiles,'zombiecodebookTHp12023reorderedGOOD.csv'));
zombiecodebookHD4reordered=table2array(zombiecodebookHD4reordered(1:end,2:end));
globalguides=cell(size(myfovs,2),1);

thatparpool=parpool(20)

parfor (fov=1:size(myfovs,2),20)
%     for fov=1:size(myfovs,2)
thisFOV=myfovs(fov);
spots=globalFoundSpotsfromBWHW2{thisFOV+1};
thisFOVsGUides=[];
%% cheating for different ilumination in round 1
% if size(spots,2)==18
% spots(:,3:5)=spots(:,3:5)*1.5-0.5;
% spots(spots<0)=0;
for ii=1:size(spots,1)
    thisguide=[];
    thisspot=spots(ii,1:end-1);     
    % thisspot(3:end)=4*thisspot(3:end)/sum(thisspot(3:end));
    distances=pdist2(thisspot(3:end),zombiecodebookHD4reordered);
    thiserror=min(distances);
    signalsum=sum(thisspot(3:end),2);
    guideID=find(distances==thiserror);
    nextup=min(distances(distances~=min(distances)));
    confidence=sum(thisspot(3:end));
    if  size(guideID,2)==1 & thiserror<1.2
        %it goes x | y | FOV | ID |error | error to next up
        thisguide=thisspot(1,1:2);
        thisguide(1,3)=thisFOV;
        thisguide(1,4)=guideID;
        thisguide(1,5)=thiserror;
        thisguide(1,6)=confidence;
    end
        thisFOVsGUides=[thisFOVsGUides;thisguide];
    end
% end
globalguides{fov}=thisFOVsGUides;
disp(strcat('did ', num2str(fov)));
end
delete(thatparpool)
decodedGuides=[];
for i=1:size(globalguides,1)
    decodedGuides=[decodedGuides;globalguides{i}];
end
writematrix(decodedGuides,fullfile(matlabtempfiles,strcat('decodedGuidesFromManualHW2.csv')))
foundguides=zeros(1,size(zombiecodebookHD4reordered,1));
for i=1:size(foundguides,2)
    foundguides(i)=sum(decodedGuides(:,4)==i);
end
writematrix(foundguides,fullfile(matlabtempfiles,strcat('foundguidesFromManualHW2.csv')))