function [MP_Gps,MP_Bds,MP_SingleGpsSequence,MP_SingleBdsSequence] = GetMultipath(Slipdata_Gps,Slipdata_Bds)
% Gps:L1/L2->L1C/L2W, Bds:B1/B3->L2I/L6I
c=299792458.0;
[L1,L2,L5,B1,B2,B3]=deal(1575420000,1227600000,1176450000,1561098000,1207140000,1268520000);
[L1_wavlen,L2_wavlen,~,B1_wavlen,~,B3_wavlen]=deal(c/L1,c/L2,c/L5,c/B1,c/B2,c/B3);
MP_Gps=zeros(length(Slipdata_Gps),3);
MP_Bds=zeros(length(Slipdata_Bds),3);
MP_SingleGpsSequence=cell(length(Slipdata_Gps),1);
MP_SingleBdsSequence=cell(length(Slipdata_Bds),1);

% GPS
%------------------------------------------------------------------------------------------------
AMB=[];
len=length(Slipdata_Gps);
for i=1:len
    num=0;
    nums=length(Slipdata_Gps{i,1});
    PRN=Slipdata_Gps{i,1}(1,1);
    Time=Slipdata_Gps{i,1}(:,2);
    CycleSlipMark=Slipdata_Gps{i,1}(:,3);
    C1C=Slipdata_Gps{i,1}(:,6);
    L1C=Slipdata_Gps{i,1}(:,7);
    C2W=Slipdata_Gps{i,1}(:,8);
    L2W=Slipdata_Gps{i,1}(:,9);
    sequence_mark=[];
    delta=Time(1)-1;
    k=1;

    % find ALL continuous time sequence without outliers or cycle slips
    while k<nums
        if k==nums
            break;
        end
        if CycleSlipMark(k)~=0
            k=k+1;
            continue;
        end
        if CycleSlipMark(k)==0
            num_start=k;
            for p=(k+1):nums
                if p==nums
                    num_end=p;
                    sequence_mark=[sequence_mark;num_start,num_end];
                    break;
                elseif CycleSlipMark(p)~=0
                    k=p;
                    num_end=p-1;
                    sequence_mark=[sequence_mark;num_start,num_end];
                    break;
                end
            end
            if p==nums
                break;
            end
        end
        k=k+1;
    end

    % find BEST continuous time sequence without outliers or cycle slips
    if length(sequence_mark(:,1))==1
        m=1;
    else
        m=1;
        max=sequence_mark(1,2)-sequence_mark(1,1);
        for w=1:length(sequence_mark)-1
            d=sequence_mark(w,2)-sequence_mark(w,1);
            d1=sequence_mark(w+1,2)-sequence_mark(w+1,1);
            if (d1>d && d1>max)
                m=w+1;
                max=d1;
            end
        end
    end

    % get span-time
    span=sequence_mark(m,:);
    spantime=zeros(span(2)-span(1)+1,1);
    t=span(1):1:span(2);
    prn=zeros(span(2)-span(1)+1,1);
    for q=1:span(2)-span(1)+1
        spantime(q)=t(q)+delta;
        prn(q)=PRN(1);
    end

    % get Single-MP-TimeSequence
    N=span(1,2)-span(1,1)+1;
    MP_Sequence=zeros(N,4);
    for j=span(1,1):span(1,2)
        MP1=C1C(j)-(L1^2+L2^2)/(L1^2-L2^2)*L1_wavlen*L1C(j)+2*L2^2/(L1^2-L2^2)*L2_wavlen*L2W(j);
        MP2=C2W(j)-2*L1^2/(L1^2-L2^2)*L1_wavlen*L1C(j)+(L1^2+L2^2)/(L1^2-L2^2)*L2_wavlen*L2W(j);
        MP_Sequence(j,3)=MP1;
        MP_Sequence(j,4)=MP2;
    end
    MP_Sequence(all(MP_Sequence==0,2),:)=[];
    MP_Sequence(:,1)=prn;
    MP_Sequence(:,2)=spantime;

    % Estimate AMBIGUITY
    mMP1=round(mean(MP_Sequence(:,3)));
    mMP2=round(mean(MP_Sequence(:,4)));
    AMB=[AMB;PRN,mMP1,mMP2];
    for u=1:length(MP_Sequence)
        MP_Sequence(u,3)=MP_Sequence(u,3)-mMP1;
        MP_Sequence(u,4)=MP_Sequence(u,4)-mMP2;
    end

    % MP_Sequence: [prn, spantime, L1_multi, L2_multi]
    MP_SingleGpsSequence(i)={MP_Sequence};

    % get bar-MP for each PRN
    b1=0;
    b2=0;
    for u=1:length(MP_Sequence)
        b1=b1+(MP_Sequence(u,3)-mMP1)^2;
        b2=b2+(MP_Sequence(u,4)-mMP2)^2;
    end
    barMP1=sqrt((1/(N-1)*b1));
    barMP2=sqrt((1/(N-1)*b2));

    MP_Gps(i,1)=PRN;
    MP_Gps(i,2)=barMP1-mMP1;
    MP_Gps(i,3)=barMP2-mMP2;

    % output correction
    if MP_Gps(i,2)>2 || MP_Gps(i,3)>2 % failure of fixing AMB 
        MP_Gps(i,2)=mean(MP_Sequence(:,3));
        MP_Gps(i,3)=mean(MP_Sequence(:,4));
    end
end
% BDS
%------------------------------------------------------------------------------------------------
AMB=[];
len=length(Slipdata_Bds);
for i=1:len
    num=0;
    nums=length(Slipdata_Bds{i,1});
    PRN=Slipdata_Bds{i,1}(1,1);
    Time=Slipdata_Bds{i,1}(:,2);
    CycleSlipMark=Slipdata_Bds{i,1}(:,3);
    C2I=Slipdata_Bds{i,1}(:,6);
    L2I=Slipdata_Bds{i,1}(:,7);
    C6I=Slipdata_Bds{i,1}(:,8);
    L6I=Slipdata_Bds{i,1}(:,9);
    sequence_mark=[];
    delta=Time(1)-1;
    k=1;

    % find ALL continuous time sequence without outliers or cycle slips
    while k<nums
        if k==nums
            break;
        end
        if CycleSlipMark(k)~=0
            k=k+1;
            continue;
        end
        if CycleSlipMark(k)==0
            num_start=k;
            for p=(k+1):nums
                if p==nums
                    num_end=p;
                    sequence_mark=[sequence_mark;num_start,num_end];
                    break;
                elseif CycleSlipMark(p)~=0
                    k=p;
                    num_end=p-1;
                    sequence_mark=[sequence_mark;num_start,num_end];
                    break;
                end
            end
            if p==nums
                break;
            end
        end
        k=k+1;
    end

    % find BEST continuous time sequence without outliers or cycle slips
    if length(sequence_mark(:,1))==1
        m=1;
    else
        m=1;
        max=sequence_mark(1,2)-sequence_mark(1,1);
        for w=1:length(sequence_mark)-1
            d=sequence_mark(w,2)-sequence_mark(w,1);
            d1=sequence_mark(w+1,2)-sequence_mark(w+1,1);
            if (d1>d && d1>max)
                m=w+1;
                max=d1;
            end
        end
    end

    % get span-time
    span=sequence_mark(m,:);
    spantime=zeros(span(2)-span(1)+1,1);
    t=span(1):1:span(2);
    prn=zeros(span(2)-span(1)+1,1);
    for q=1:span(2)-span(1)+1
        spantime(q)=t(q)+delta;
        prn(q)=PRN(1);
    end

    % get Single-MP-TimeSequence
    N=span(1,2)-span(1,1)+1;
    MP_Sequence=zeros(N,4);
    for j=span(1,1):span(1,2)
        MP1=C2I(j)-(B1^2+B3^2)/(B1^2-B3^2)*B1_wavlen*L2I(j)+2*B3^2/(B1^2-B3^2)*B3_wavlen*L6I(j);
        MP2=C6I(j)-2*B1^2/(B1^2-B3^2)*B1_wavlen*L2I(j)+(B1^2+B3^2)/(B1^2-B3^2)*B3_wavlen*L6I(j);
        MP_Sequence(j,3)=MP1;
        MP_Sequence(j,4)=MP2;
    end
    MP_Sequence(all(MP_Sequence==0,2),:)=[];
    MP_Sequence(:,1)=prn;
    MP_Sequence(:,2)=spantime;

    % Estimate AMBIGUITY
    mMP1=round(mean(MP_Sequence(:,3)));
    mMP2=round(mean(MP_Sequence(:,4)));
    AMB=[AMB;PRN,mMP1,mMP2];
    for u=1:length(MP_Sequence)
        MP_Sequence(u,3)=MP_Sequence(u,3)-mMP1;
        MP_Sequence(u,4)=MP_Sequence(u,4)-mMP2;
    end

    % MP_Sequence: [prn, spantime, L1_multi, L2_multi]
    MP_SingleBdsSequence(i)={MP_Sequence};

    % get bar-MP for each PRN
    b1=0;
    b2=0;
    for u=1:length(MP_Sequence)
        b1=b1+(MP_Sequence(u,3)-mMP1)^2;
        b2=b2+(MP_Sequence(u,4)-mMP2)^2;
    end
    barMP1=sqrt((1/(N-1)*b1));
    barMP2=sqrt((1/(N-1)*b2));

    MP_Bds(i,1)=PRN;
    MP_Bds(i,2)=barMP1-mMP1;
    MP_Bds(i,3)=barMP2-mMP2;

    % output correction
    if MP_Bds(i,2)>2 || MP_Bds(i,3)>2 % failure of fixing AMB 
        MP_Bds(i,2)=mean(MP_Sequence(:,3));
        MP_Bds(i,3)=mean(MP_Sequence(:,4));
    end

end

end

