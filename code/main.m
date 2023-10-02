clc,clear,close all; %designed by hzLiu, July 9th,2023
%% Load data
filename = "../renix/D005.23o"; 
data = rinexread(filename);
structGps=table2struct(data.GPS,"ToScalar",true);
structBds=table2struct(data.BeiDou,"ToScalar",true);
% calculate time in day
[GpsDateSec,BdsDateSec]=GetSatelliteTime(data.GPS,data.BeiDou);
% add field in struct: Time
structGps.Time=GpsDateSec;
structBds.Time=BdsDateSec;
% get target type: table
tableGps=struct2table(structGps);
tableBds=struct2table(structBds);
clear GpsDateSec BdsDateSec GpsDateString BdsDateString i sec strh strm strs filename;
disp("Loading is done");
% carrier frequency / wavelength / c : Gps Bds
c=299792458.0;
[L1,L2,L5,B1,B2,B3]=deal(1575420000,1227600000,1176450000,1561098000,1207140000,1268520000);
[L1_wavlen,L2_wavlen,L5_wavlen,B1_wavlen,B2_wavlen,B3_wavlen]=deal(c/L1,c/L2,c/L5,c/B1,c/B2,c/B3);

%% Get all single PRN data: Gps/Bds
allGpsPRNData=GetPRNData(tableGps);
allBdsPRNData=GetPRNData(tableBds);
disp("Have already got all single PRN data: Gps/Bds");
%% Intergrity rate of observations
[DI_Sys_Gps,DI_Sys_Bds]=GetSysIntergrity(structGps,structBds);
[DI_matGps,DI_matBds]=GetFrqIntergrity(allGpsPRNData,allBdsPRNData);
disp("Have already got results of Intergrity rate of observations: Gps/Bds");
%% SNR data
[SNR_Sys_Gps,SNR_Sys_Bds,SNR_matGps,SNR_matBds]=GetSNR(allGpsPRNData,allBdsPRNData);
disp("Have already got results of SNR data: Gps/Bds");
%% Carrier Phase Cycle Slip Detection
[CycleSlip_Gps,CycleSlip_Bds,CycleSlipRatio_Gps,CycleSlipRatio_Bds]=GetCycleSlip(allGpsPRNData,allBdsPRNData);
disp("Have already got results of Cycle Slip data: Gps/Bds");
%% Pseudo-range multipath effects
[MP_Gps,MP_Bds,MP_SingleGpsSequence,MP_SingleBdsSequence]=GetMultipath(CycleSlip_Gps,CycleSlip_Bds);
disp("Have already got results of Pseudo-range multipath effects: Gps/Bds");

%% Visualization
Visualization(CycleSlip_Gps,CycleSlip_Bds,MP_SingleGpsSequence,MP_SingleBdsSequence,allGpsPRNData,allBdsPRNData);
disp("Visualization is done");













