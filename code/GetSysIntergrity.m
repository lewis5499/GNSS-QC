function [DI_Sys_Gps,DI_Sys_Bds] = GetSysIntergrity(struct_SysGpsData,struct_SysBdsData)
len=length(struct_SysGpsData.SatelliteID);
num=0;
for i=1:len
    if struct_SysGpsData.C1C(i)~=0 && struct_SysGpsData.C2W(i)~=0 && struct_SysGpsData.L1C(i)~=0 && struct_SysGpsData.L2W(i)~=0
    %if struct_SysGpsData.C1C(i)~=0 && struct_SysGpsData.C2W(i)~=0 && struct_SysGpsData.C5X(i)~=0 && struct_SysGpsData.L1C(i)~=0 && struct_SysGpsData.L2W(i)~=0 && struct_SysGpsData.L5X(i)~=0
        num=num+1;
    end
end
DI_Sys_Gps=num/len;

len=length(struct_SysBdsData.SatelliteID);
num=0;
for i=1:len
    if struct_SysBdsData.C2I(i)~=0 && struct_SysBdsData.C6I(i)~=0 && struct_SysBdsData.L2I(i)~=0 && struct_SysBdsData.L6I(i)~=0
    %if struct_SysBdsData.C2I(i)~=0 && struct_SysBdsData.C7I(i)~=0 && struct_SysBdsData.C6I(i)~=0 && struct_SysBdsData.L2I(i)~=0 && struct_SysBdsData.L7I(i)~=0 && struct_SysBdsData.L6I(i)~=0
        num=num+1;
    end
end
DI_Sys_Bds=num/len;

end

