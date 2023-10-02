function [CycleSlip_Gps,CycleSlip_Bds,CycleSlipRatio_Gps,CycleSlipRatio_Bds] = GetCycleSlip(cell_allGps,cell_allBds)
CycleSlip_Gps=cell(length(cell_allGps),1);
CycleSlip_Bds=cell(length(cell_allBds),1);
CycleSlipRatio_Gps=cell(length(cell_allGps),2);
CycleSlipRatio_Bds=cell(length(cell_allBds),2);
c=299792458.0;
%Gps:L1/L2->L1C/L2W, Bds:B1/B3->L2I/L6I
[L1,L2,L5,B1,B2,B3]=deal(1575420000,1227600000,1176450000,1561098000,1207140000,1268520000);
[L1_wavlen,L2_wavlen,~,B1_wavlen,~,B3_wavlen]=deal(c/L1,c/L2,c/L5,c/B1,c/B2,c/B3);
%Undifferenced pseudorange accuracy
SigmaC=0.3;

%GPS
%----------------------------------------------------------------------------------------------------------------------------------------------------
wavlen_mw=c/(L1-L2);
len=length(cell_allGps);
for i=1:len
    %disp('i=');disp(i);
    nums=length(cell_allGps{i,1}.SatelliteID);
    PRN=cell_allGps{i,1}.SatelliteID;
    TimeSequence=cell_allGps{i,1}.Time;
    CycleSlipMark=zeros(nums,1);%CycleSlipMark: 0:healthy 1:cycle slip 2:outlier -1:lock-loss
    Nmw_list=[];
    Lmw_list=[];
    LGF_list=[];
    PGF_list=[];
    rms0_Lmw=SigmaC/(L1+L2)*sqrt(L1^2+L2^2);
    % mark the Time-Interruption and calculate the MW value for every epoch
    for j=1:nums
        %judge the lock-loss
        if cell_allGps{i,1}.L1C(j)==0 || cell_allGps{i,1}.L2W(j)==0 || cell_allGps{i,1}.C1C(j)==0 || cell_allGps{i,1}.C2W(j)==0
            CycleSlipMark(j)=-1; 
        end
        %save the MW value
        Nmw=(cell_allGps{i,1}.L1C(j)-cell_allGps{i,1}.L2W(j))-1/((L1+L2)*wavlen_mw)*(L1*cell_allGps{i,1}.C1C(j)+L2*cell_allGps{i,1}.C2W(j));
        Nmw_list=[Nmw_list;Nmw];
        Lmw_list=wavlen_mw.*Nmw_list;
        %save the GF value
        LGF=L1_wavlen*cell_allGps{i,1}.L1C(j)-L2_wavlen*cell_allGps{i,1}.L2W(j);
        PGF=cell_allGps{i,1}.C1C(j)-cell_allGps{i,1}.C2W(j);
        LGF_list=[LGF_list;LGF];
        PGF_list=[PGF_list;PGF];
    end
    %Get dLGF
    dLGF_list=zeros(length(LGF_list)-1,1);
    for m=2:length(LGF_list)
        dLGF_list(m-1)=LGF_list(m)-LGF_list(m-1);
    end
    %Get ddLGF
    ddLGF_list=zeros(length(LGF_list)-2,1);
    for m=2:length(dLGF_list)
        ddLGF_list(m-1)=dLGF_list(m)-dLGF_list(m-1);
    end

    k=1;%k th MW observation as the fisrt MW value for a continuous time span
    n=1;%n th continuous time span
    while true
        if CycleSlipMark(k)==0 && k~=1 
            k=k+1;
            if k==nums || k==nums+1
                break;
            end
        elseif CycleSlipMark(k)==-1
            %at least 3 continuous obs epoches for a time span for detection
            while CycleSlipMark(k)==-1 || CycleSlipMark(k+1)==-1 || CycleSlipMark(k+2)==-1 
                k=k+1;
                if k==nums-2
                    break;
                end
            end
        end  
        if k==nums-2 || k==nums-1
            break;
        end
        %disp('k=');disp(k);
        % especially, j=k: initialize rms
        rms_Lmw=rms0_Lmw;
        mean_Lmw=Lmw_list(k);

        rms_Lmw_list=[rms_Lmw];%initialize list
        mean_Lmw_list=[mean_Lmw];

        % j=k+1
        rms_Lmw=sqrt((2-1)/2*rms_Lmw^2+1/2*(Lmw_list(k+1)-mean_Lmw)^2);
        mean_Lmw=(2-1)/2*mean_Lmw+1/2*Lmw_list(k+1);
        
        rms_Lmw_list=[rms_Lmw_list;rms_Lmw];
        mean_Lmw_list=[mean_Lmw_list;mean_Lmw];

        width=2;
        for j=(k+2):nums
            width=width+1;
            if CycleSlipMark(j)==-1
                break;
            end
            rms_Lmw=sqrt((j-k)/(j-k+1)*rms_Lmw^2+1/(j-k+1)*(Lmw_list(j)-mean_Lmw)^2);
            mean_Lmw=(j-k)/(j-k+1)*mean_Lmw+1/(j-k+1)*Lmw_list(j);

            % 1:cycle slip 2:outlier
            % we start from 'k+1'th MW value, (meanly j-1,j=k+2)
            judge1=abs(Lmw_list(j-1)-mean_Lmw_list(j-k+1-2));
            threshold1=4*rms_Lmw_list(j-k+1-2);
            judge2=abs(Lmw_list(j)-Lmw_list(j-1));
            threshold2=1;
            judge3=abs(Lmw_list(j)-mean_Lmw_list(j-k+1-1));
            threshold3=4*rms_Lmw_list(j-k+1-1);
            
            if width==3
                mm=0;
                nn=0;
            else
                mm=ddLGF_list(j-3)+ddLGF_list(j-2);
                nn=ddLGF_list(j-3);
            end
            
            if (judge1>=threshold1 && judge3>=threshold3 && judge2>threshold2) || (judge1>=threshold1 && judge3<threshold3)
                CycleSlipMark(j-1)=2;
                break;
            elseif judge1>=threshold1 && judge3>=threshold3 && judge2<=threshold2
                CycleSlipMark(j-1)=1;
                break;
            elseif abs(nn)>=0.02 && abs(mm)<0.006 && width>3 
                CycleSlipMark(j-1)=1;
                break;
            else
                rms_Lmw_list=[rms_Lmw_list;rms_Lmw];
                mean_Lmw_list=[mean_Lmw_list;mean_Lmw];
            end            
        end
        % save the data for each continuous time span, aligned with Timesequence
        k=j;
        %disp('j=');disp(k);       
        if(j==nums)
            break;
        end
        n=n+1;
    end  
    
    ratio=tabulate(CycleSlipMark);
    CycleSlipRatio_Gps(i,:)={{PRN(1)},{ratio}};

    data=[PRN,TimeSequence,CycleSlipMark,Nmw_list,Lmw_list,cell_allGps{i,1}.C1C,cell_allGps{i,1}.L1C,cell_allGps{i,1}.C2W,cell_allGps{i,1}.L2W];
    CycleSlip_Gps(i,1)={data};
end

%BDS
%----------------------------------------------------------------------------------------------------------------------------------------------------
wavlen_mw=c/(B1-B3);
len=length(cell_allBds);
for i=1:len
    %disp('i=');disp(i);
    nums=length(cell_allBds{i,1}.SatelliteID);
    PRN=cell_allBds{i,1}.SatelliteID;
    TimeSequence=cell_allBds{i,1}.Time;
    CycleSlipMark=zeros(nums,1);%CycleSlipMark: 0:healthy 1:cycle slip 2:outlier -1:lock-loss
    Nmw_list=[];
    Lmw_list=[];
    LGF_list=[];
    PGF_list=[];
    rms0_Lmw=SigmaC/(B1+B3)*sqrt(B1^2+B3^2);
    % mark the Time-Interruption and calculate the MW value for every epoch
    for j=1:nums
        %judge the lock-loss
        if cell_allBds{i,1}.L2I(j)==0 || cell_allBds{i,1}.L6I(j)==0 || cell_allBds{i,1}.C2I(j)==0 || cell_allBds{i,1}.C6I(j)==0
            CycleSlipMark(j)=-1; 
        end
        %save the MW value
        Nmw=(cell_allBds{i,1}.L2I(j)-cell_allBds{i,1}.L6I(j))-1/((B1+B3)*wavlen_mw)*(B1*cell_allBds{i,1}.C2I(j)+B3*cell_allBds{i,1}.C6I(j));
        Nmw_list=[Nmw_list;Nmw];
        Lmw_list=wavlen_mw.*Nmw_list;
        %save the GF value
        LGF=B1_wavlen*cell_allBds{i,1}.L2I(j)-B3_wavlen*cell_allBds{i,1}.L6I(j);
        PGF=cell_allBds{i,1}.C2I(j)-cell_allBds{i,1}.C6I(j);
        LGF_list=[LGF_list;LGF];
        PGF_list=[PGF_list;PGF];
    end
    %Get dLGF
    dLGF_list=zeros(length(LGF_list)-1,1);
    for m=2:length(LGF_list)
        dLGF_list(m-1)=LGF_list(m)-LGF_list(m-1);
    end
    %Get ddLGF
    ddLGF_list=zeros(length(LGF_list)-2,1);
    for m=2:length(dLGF_list)
        ddLGF_list(m-1)=dLGF_list(m)-dLGF_list(m-1);
    end

    k=1;%k th MW observation as the fisrt MW value for a continuous time span
    n=1;%n th continuous time span
    while true
        if CycleSlipMark(k)==0 && k~=1 
            k=k+1;
            if k==nums || k==nums+1
                break;
            end
        elseif CycleSlipMark(k)==-1
            %at least 3 continuous obs epoches for a time span for detection
            while CycleSlipMark(k)==-1 || CycleSlipMark(k+1)==-1 || CycleSlipMark(k+2)==-1 
                k=k+1;
                if k==nums-2
                    break;
                end
            end
        end  
        if k==nums-2 || k==nums-1
            break;
        end
        %disp('k=');disp(k);
        % especially, j=k: initialize rms
        rms_Lmw=rms0_Lmw;
        mean_Lmw=Lmw_list(k);

        rms_Lmw_list=[rms_Lmw];%initialize list
        mean_Lmw_list=[mean_Lmw];

        % j=k+1
        rms_Lmw=sqrt((2-1)/2*rms_Lmw^2+1/2*(Lmw_list(k+1)-mean_Lmw)^2);
        mean_Lmw=(2-1)/2*mean_Lmw+1/2*Lmw_list(k+1);
        
        rms_Lmw_list=[rms_Lmw_list;rms_Lmw];
        mean_Lmw_list=[mean_Lmw_list;mean_Lmw];

        width=2;
        for j=(k+2):nums
            width=width+1;
            if CycleSlipMark(j)==-1
                break;
            end
            rms_Lmw=sqrt((j-k)/(j-k+1)*rms_Lmw^2+1/(j-k+1)*(Lmw_list(j)-mean_Lmw)^2);
            mean_Lmw=(j-k)/(j-k+1)*mean_Lmw+1/(j-k+1)*Lmw_list(j);

            % 1:cycle slip 2:outlier
            % we start from 'k+1'th MW value, (meanly j-1,j=k+2)
            judge1=abs(Lmw_list(j-1)-mean_Lmw_list(j-k+1-2));
            threshold1=4*rms_Lmw_list(j-k+1-2);
            judge2=abs(Lmw_list(j)-Lmw_list(j-1));
            threshold2=1;
            judge3=abs(Lmw_list(j)-mean_Lmw_list(j-k+1-1));
            threshold3=4*rms_Lmw_list(j-k+1-1);
            
            if width==3
                mm=0;
                nn=0;
            else
                mm=ddLGF_list(j-3)+ddLGF_list(j-2);
                nn=ddLGF_list(j-3);
            end
            
            if (judge1>=threshold1 && judge3>=threshold3 && judge2>threshold2) || (judge1>=threshold1 && judge3<threshold3)
                CycleSlipMark(j-1)=2;
                break;
            elseif judge1>=threshold1 && judge3>=threshold3 && judge2<=threshold2
                CycleSlipMark(j-1)=1;
                break;
            elseif abs(nn)>=0.02 && abs(mm)<0.006 && width>3 
                CycleSlipMark(j-1)=1;
                break;
            else
                rms_Lmw_list=[rms_Lmw_list;rms_Lmw];
                mean_Lmw_list=[mean_Lmw_list;mean_Lmw];
            end            
        end
        % save the data for each continuous time span, aligned with Timesequence
        k=j;
        %disp('j=');disp(k);       
        if(j==nums)
            break;
        end
        n=n+1;
    end   

    ratio=tabulate(CycleSlipMark);
    CycleSlipRatio_Bds(i,:)={{PRN(1)},{ratio}};

    data=[PRN,TimeSequence,CycleSlipMark,Nmw_list,Lmw_list,cell_allBds{i,1}.C2I,cell_allBds{i,1}.L2I,cell_allBds{i,1}.C6I,cell_allBds{i,1}.L6I];
    CycleSlip_Bds(i,1)={data};
end

end

