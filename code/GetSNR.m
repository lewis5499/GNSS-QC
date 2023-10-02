function [SNR_Sys_Gps,SNR_Sys_Bds,SNR_matGps,SNR_matBds] = GetSNR(cell_allGps,cell_allBds)
SNR_matGps=zeros(length(cell_allGps),4);%['PRN','SNR_L1','SNR_L2','SNR_L5']
SNR_matBds=zeros(length(cell_allBds),4);%['PRN','SNR_C2','SNR_C7','SNR_C6']
SNR_Sys_Gps=zeros(4,1);%['AllFreq_meanSNR','meanSNR_L1','meanSNR_L2','meanSNR_L5']
SNR_Sys_Bds=zeros(4,1);%['AllFreq_meanSNR','meanSNR_C2','meanSNR_C7','meanSNR_C6']
sumSNR_Sys_Gps=0;
sumSNR1_Sys_Gps=0;
sumSNR2_Sys_Gps=0;
sumSNR5_Sys_Gps=0;
sumnums=0;
sumnums1=0;
sumnums2=0;
sumnums5=0;

len=length(cell_allGps);
for i=1:len
    sumSNR1=0;
    sumSNR2=0;
    sumSNR5=0;
    nums=length(cell_allGps{i,1}.SatelliteID);
    nums1=length(cell_allGps{i,1}.S1C(cell_allGps{i,1}.S1C~=0));
    nums2=length(cell_allGps{i,1}.S2W(cell_allGps{i,1}.S2W~=0));
    nums5=length(cell_allGps{i,1}.S5X(cell_allGps{i,1}.S5X~=0));
    PRN=cell_allGps{i,1}.SatelliteID(1);
    for j=1:nums
        sumSNR1=sumSNR1+cell_allGps{i,1}.S1C(j);
        sumSNR2=sumSNR2+cell_allGps{i,1}.S2W(j);
        sumSNR5=sumSNR5+cell_allGps{i,1}.S5X(j);
    end
    [SNR_matGps(i,1),SNR_matGps(i,2),SNR_matGps(i,3),SNR_matGps(i,4)]=deal(PRN,sumSNR1/nums1,sumSNR2/nums2,sumSNR5/nums5);
    sumSNR_Sys_Gps=sumSNR_Sys_Gps+sumSNR1+sumSNR2+sumSNR5;
    sumSNR1_Sys_Gps=sumSNR1_Sys_Gps+sumSNR1;
    sumSNR2_Sys_Gps=sumSNR2_Sys_Gps+sumSNR2;
    sumSNR5_Sys_Gps=sumSNR5_Sys_Gps+sumSNR5;
    sumnums=sumnums+nums1+nums2+nums5;
    sumnums1=sumnums1+nums1;
    sumnums2=sumnums2+nums2;
    sumnums5=sumnums5+nums5;
end
[SNR_Sys_Gps(1,1),SNR_Sys_Gps(2,1),SNR_Sys_Gps(3,1),SNR_Sys_Gps(4,1)]=deal(sumSNR_Sys_Gps/sumnums,sumSNR1_Sys_Gps/sumnums1,sumSNR2_Sys_Gps/sumnums2,sumSNR5_Sys_Gps/sumnums5);

sumSNR_Sys_Bds=0;
sumSNR2_Sys_Bds=0;
sumSNR7_Sys_Bds=0;
sumSNR6_Sys_Bds=0;
sumnums=0;
sumnums2=0;
sumnums7=0;
sumnums6=0;
len=length(cell_allBds);
for i=1:len
    sumSNR2=0;
    sumSNR7=0;
    sumSNR6=0;
    nums=length(cell_allBds{i,1}.SatelliteID);
    nums2=length(cell_allBds{i,1}.S2I(cell_allBds{i,1}.S2I~=0));
    nums7=length(cell_allBds{i,1}.S7I(cell_allBds{i,1}.S7I~=0));
    nums6=length(cell_allBds{i,1}.S6I(cell_allBds{i,1}.S6I~=0));
    PRN=cell_allBds{i,1}.SatelliteID(1);
    for j=1:nums
        sumSNR2=sumSNR2+cell_allBds{i,1}.S2I(j);
        sumSNR7=sumSNR7+cell_allBds{i,1}.S7I(j);
        sumSNR6=sumSNR6+cell_allBds{i,1}.S6I(j);
    end
    [SNR_matBds(i,1),SNR_matBds(i,2),SNR_matBds(i,3),SNR_matBds(i,4)]=deal(PRN,sumSNR2/nums2,sumSNR7/nums7,sumSNR6/nums6);
    sumSNR_Sys_Bds=sumSNR_Sys_Bds+sumSNR2+sumSNR7+sumSNR6;
    sumSNR2_Sys_Bds=sumSNR2_Sys_Bds+sumSNR2;
    sumSNR7_Sys_Bds=sumSNR7_Sys_Bds+sumSNR7;
    sumSNR6_Sys_Bds=sumSNR6_Sys_Bds+sumSNR6;
    sumnums=sumnums+nums2+nums7+nums6;
    sumnums2=sumnums2+nums2;
    sumnums7=sumnums7+nums7;
    sumnums6=sumnums6+nums6;
end
[SNR_Sys_Bds(1,1),SNR_Sys_Bds(2,1),SNR_Sys_Bds(3,1),SNR_Sys_Bds(4,1)]=deal(sumSNR_Sys_Bds/sumnums,sumSNR2_Sys_Bds/sumnums2,sumSNR7_Sys_Bds/sumnums7,sumSNR6_Sys_Bds/sumnums6);
end





