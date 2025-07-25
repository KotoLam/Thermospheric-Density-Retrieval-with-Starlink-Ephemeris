function eopdata = ReadEOPData(processDate)

fid = fopen('EOP-Last5Years.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s
%  ----------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k ~= 0)
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for cnt=1:numrecsobs
            eopdata(:,cnt) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for cnt=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);
        for cnt=numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,cnt) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

% 
dateVec = datevec(processDate);
mjd = mjuliandate(dateVec);
ind = find(mjd==eopdata(4,:),1,'first');
eopdata = eopdata(:,ind-5:ind+5);

end