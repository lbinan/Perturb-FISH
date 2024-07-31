zombiecodebookHD4reordered=readtable(fullfile(matlabtempfiles,'zombiecodebookTHp12023reorderedGOOD.csv'));
zombiecodebookHD4reordered=table2array(zombiecodebookHD4reordered(1:end,2:end));
globalguides=cell(size(myfovs,2),1);
thatparpool=parpool(22)
parfor (fov=1:size(myfovs,2),22)
    thisFOV=myfovs(fov);
    spots=globalFoundSpotsfromBW{thisFOV+1};
    thisFOVsGUides=[];
    if size(spots,2)==18
        for ii=1:size(spots,1)
            thisguide=[];
            thisspot=spots(ii,1:17);     
            distances=pdist2(thisspot(3:end),zombiecodebookHD4reordered);%compute distance to all words from codebook, to find closest one
            thiserror=min(distances);
            guideID=find(distances==thiserror);%find the guide ID with smallest distance
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
    end
    globalguides{fov}=thisFOVsGUides;
    disp(strcat('did ', num2str(fov)));
end
delete(thatparpool)
decodedGuides=[];
for i=1:size(globalguides,1)
    decodedGuides=[decodedGuides;globalguides{i}];
end
writematrix(decodedGuides,fullfile(matlabtempfiles,strcat('decodedGuidesmanualV2dontcheckrepeatslookfurther.csv')))
foundguides=zeros(1,size(zombiecodebookHD4reordered,1));
for i=1:size(foundguides,2)
    foundguides(i)=sum(decodedGuides(:,4)==i);%get deteceted guide distribution
end
writematrix(foundguides,fullfile(matlabtempfiles,strcat('foundguidesmanualV2dontcheckrepeatslookfurther.csv')))