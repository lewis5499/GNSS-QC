function c = GetPRNData(tableData)
%matData:['Time','allPRNnum',(2)
%         'C1C','C1C_SSI','C2W','C2W_SSI','C5X','C5X_SSI',(8)
%         'L1C','L1C_LLI','L1C_SSI','L2W','L2W_LLI','L2W_SSI','L5X','L5X_LLI','L5X_SSI',(17)
%         'S1C','S2W','S5X'(20)]
c={};
for PRNnum=1:60
    rows = (tableData.EpochFlag==0 & tableData.SatelliteID==PRNnum);
    if all(rows==0)
        continue;
    end
    PRNtableData=tableData(rows,:);
    PRNtableData.Properties.Description=int2str(PRNnum);
    sPRN=table2struct(PRNtableData,"ToScalar",true);
    c=[c;{sPRN}];
end
end
