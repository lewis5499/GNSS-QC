function [DI_matGps,DI_matBds] = GetFrqIntergrity(cell_allGps,cell_allBds)
DI_matGps=zeros(length(cell_allGps),2);
DI_matBds=zeros(length(cell_allBds),2);

len=length(cell_allGps);
for i=1:len
    num=0;
    nums=length(cell_allGps{i,1}.SatelliteID);
    PRN=cell_allGps{i,1}.SatelliteID(1);
    for j=1:nums
        %if cell_allGps{i,1}.C1C(j)~=0 && cell_allGps{i,1}.C2W(j)~=0 && cell_allGps{i,1}.C5X(j)~=0 && cell_allGps{i,1}.L1C(j)~=0 && cell_allGps{i,1}.L2W(j)~=0 && cell_allGps{i,1}.L5X(j)~=0
        if cell_allGps{i,1}.C1C(j)~=0 && cell_allGps{i,1}.C2W(j)~=0 && cell_allGps{i,1}.L1C(j)~=0 && cell_allGps{i,1}.L2W(j)~=0
           num=num+1;
        end
    end
    [DI_matGps(i,1),DI_matGps(i,2)]=deal(PRN,num/nums);
end

len=length(cell_allBds);
for i=1:len
    num=0;
    nums=length(cell_allBds{i,1}.SatelliteID);
    PRN=cell_allBds{i,1}.SatelliteID(1);
    for j=1:nums
        if cell_allBds{i,1}.C2I(j)~=0 && cell_allBds{i,1}.C6I(j)~=0 && cell_allBds{i,1}.L2I(j)~=0 && cell_allBds{i,1}.L6I(j)~=0
        %if cell_allBds{i,1}.C2I(j)~=0 && cell_allBds{i,1}.C7I(j)~=0 && cell_allBds{i,1}.C6I(j)~=0 && cell_allBds{i,1}.L2I(j)~=0 && cell_allBds{i,1}.L7I(j)~=0 && cell_allBds{i,1}.L6I(j)~=0
           num=num+1;
        end
    end
    [DI_matBds(i,1),DI_matBds(i,2)]=deal(PRN,num/nums);
end
end

