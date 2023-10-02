function [GpsDateSec,BdsDateSec] = GetSatelliteTime(dataGps,dataBds)
GpsDateString = datestr(dataGps.Time);
GpsDateSec=zeros(length(GpsDateString),1);
for i=1:length(GpsDateString)
    strh=GpsDateString(i,13:14);
    strm=GpsDateString(i,16:17);
    strs=GpsDateString(i,19:20);
    sec=str2double(strh)*3600+str2double(strm)*60+str2double(strs);
    GpsDateSec(i,1)=sec;
end
BdsDateString = datestr(dataBds.Time);
BdsDateSec=zeros(length(BdsDateString),1);
for i=1:length(BdsDateString)
    strh=BdsDateString(i,13:14);
    strm=BdsDateString(i,16:17);
    strs=BdsDateString(i,19:20);
    sec=str2double(strh)*3600+str2double(strm)*60+str2double(strs);
    BdsDateSec(i,1)=sec;
end
end

