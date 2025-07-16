function Ephemeris2DeltaE(processDate, inputPath, outputPath)

if isempty(regexp(processDate, '^\d{4}-\d{2}-\d{2}$', 'once'))
    error('处理日期格式应为 "yyyy-mm-dd"');
end

foldersAll = dir(inputPath);
isDir = [foldersAll.isdir];
folderNames = {foldersAll(isDir).name};
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));
folderNames = sort(folderNames);
processFolders = folderNames(~cellfun(@isempty, regexp(folderNames, processDate)));

if isempty(processFolders)
    error('没有找到符合日期 %s 的文件夹', processDate);
end

eopdata = ReadEOPData(processDate);
load("data/DE440_Coeff.mat","DE440Coeff")
load("data/EGM2008_Coeff.mat","Cnm","Snm");

if isunix && ~ismac
    poolobj = gcp("nocreate");
    if isempty(poolobj)
        parpool(16)
    end
end

files = dir(fullfile(inputPath,processFolders{cnt_folder},"*.mat"));

MAX_COUNT_FILE = length(files);
fprintf("files counts: %d\n",MAX_COUNT_FILE)
for cnt_file = 1:MAX_COUNT_FILE
    tic;
    fname = fullfile(files(cnt_file).folder,files(cnt_file).name);
    s = load(fname);
    Datevec_UTC = s.Datevec_UTC;
    Pos_ECI = s.Pos_ECI;
    Vel_ECI = s.Vel_ECI;
    Cov = s.Cov;

    
    MAX_COUNT_TIME = 2*24*60+1; % 2天

    Datevec_UTC = Datevec_UTC(1:MAX_COUNT_TIME,:);
    Datetime_UTC = datetime(Datevec_UTC);
    r_ECI = Pos_ECI(1:MAX_COUNT_TIME,:)*1e3; % km -> m
    v_ECI = Vel_ECI(1:MAX_COUNT_TIME,:)*1e3; % km/s -> m/s
    Cov = Cov(1:MAX_COUNT_TIME,:);

    Mjd_UTC = mjuliandate(Datevec_UTC);
    fprintf("%d/%d start %s\n",cnt_file,MAX_COUNT_FILE,fname)
    degree = 20;
    r_ECEF = nan(size(r_ECI));
    v_ECEF = nan(size(v_ECI));
    OE_Osc = nan(MAX_COUNT_TIME,6);
    OE_Mean = nan(MAX_COUNT_TIME,6);

    DELTA_T = 60;
    STEP   = 60;   % [s] integration step size
    Delta_E = nan(MAX_COUNT_TIME,1);
    for cnt_data=1:MAX_COUNT_TIME-1
        Mjd_UTC_Cur = Mjd_UTC(cnt_data);
        r_ECI_Cur = r_ECI(cnt_data,:);
        v_ECI_Cur = v_ECI(cnt_data,:);
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_Cur);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        MJD_TT  = Mjd_UTC_Cur + TT_UTC/86400;
        MJD_TDB = Mjday_TDB(MJD_TT);
        JD_TDB = MJD_TDB+2400000.5;
        [r_Moon_Cur,r_Sun_Cur,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
        %
        Mjd_UTC_Next = Mjd_UTC(cnt_data+1);
        r_ECI_Next = r_ECI(cnt_data+1,:);
        v_ECI_Next = v_ECI(cnt_data+1,:);
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_Next);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        MJD_TT  = Mjd_UTC_Next + TT_UTC/86400;
        MJD_TDB = Mjday_TDB(MJD_TT);
        JD_TDB = MJD_TDB+2400000.5;
        [r_Moon_Next,r_Sun_Next,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);

        V_K_Cur = 0.5*norm(v_ECI_Cur)^2;
        V_K_Next = 0.5*norm(v_ECI_Next)^2;
        DeltaE_K = V_K_Next - V_K_Cur;

        V_Moon_Cur = Potential_Moon(r_ECI_Cur,r_Moon_Cur);
        V_Moon_Next = Potential_Moon(r_ECI_Next,r_Moon_Next);
        DeltaE_Moon = V_Moon_Next - V_Moon_Cur;

        V_Sun_Cur = Potential_Sun(r_ECI_Cur,r_Sun_Cur);
        V_Sun_Next = Potential_Sun(r_ECI_Next,r_Sun_Next);
        DeltaE_Sun = V_Sun_Next - V_Sun_Cur;

        h = STEP; % Step-size of integration [s]
        span = 0:STEP:DELTA_T;
        num = length(span);
        Work = 0;
        for ii = 1:num-1
            Mjd_UTC_Left = Mjd_UTC_Cur + span(ii)/86400;
            cnt_t = find(Mjd_UTC_Left>=Mjd_UTC,1,"last");
            if cnt_t < 5
                ind = 1:8;
            elseif cnt_t > MAX_COUNT_TIME-5
                ind = MAX_COUNT_TIME-7:MAX_COUNT_TIME;
            else
                ind = cnt_t-3:cnt_t+4;
            end
            Mjd_UTC_Ind = Mjd_UTC(ind);
            r_ECI_Ind = r_ECI(ind,:);
            v_ECI_Ind = v_ECI(ind,:);

            r_ECI_Ind = r_ECI_Ind';
            v_ECI_Ind = v_ECI_Ind';

            a1 = 0;
            Mjd_UTC_a1 = Mjd_UTC_Left + a1*h/86400;
            r_ECI_a1=[LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(1,:),Mjd_UTC_a1);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(2,:),Mjd_UTC_a1);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(3,:),Mjd_UTC_a1);];
            v_ECI_a1=[LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(1,:),Mjd_UTC_a1);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(2,:),Mjd_UTC_a1);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(3,:),Mjd_UTC_a1);];
            [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a1);
            [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT  = Mjd_UTC_a1 + TT_UTC/86400;
            MJD_TDB = Mjday_TDB(MJD_TT);
            JD_TDB = MJD_TDB+2400000.5;
            [r_Moon_a1,r_Sun_a1,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
            acc1 = Accel_EarthGravity_Tides(Mjd_UTC_a1,r_Sun_a1,r_Moon_a1,...
                r_ECI_a1,Cnm,Snm,eopdata);

            a11 = 1;
            Mjd_UTC_a11 = Mjd_UTC_Left + a11*h/86400;
            r_ECI_a11=[LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(1,:),Mjd_UTC_a11);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(2,:),Mjd_UTC_a11);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(3,:),Mjd_UTC_a11);];
            v_ECI_a11=[LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(1,:),Mjd_UTC_a11);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(2,:),Mjd_UTC_a11);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(3,:),Mjd_UTC_a11);];
            [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a11);
            [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT  = Mjd_UTC_a11 + TT_UTC/86400;
            MJD_TDB = Mjday_TDB(MJD_TT);
            JD_TDB = MJD_TDB+2400000.5;
            [r_Moon_a11,r_Sun_a11,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
            acc11 = Accel_EarthGravity_Tides(Mjd_UTC_a11,r_Sun_a11,r_Moon_a11,...
                r_ECI_a11,Cnm,Snm,eopdata);

            a6 = 1.0 / 2.0;
            Mjd_UTC_a6 = Mjd_UTC_Left + a6*h/86400;
            r_ECI_a6=[LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(1,:),Mjd_UTC_a6);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(2,:),Mjd_UTC_a6);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(3,:),Mjd_UTC_a6);];
            v_ECI_a6=[LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(1,:),Mjd_UTC_a6);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(2,:),Mjd_UTC_a6);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(3,:),Mjd_UTC_a6);];
            [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a6);
            [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT  = Mjd_UTC_a6 + TT_UTC/86400;
            MJD_TDB = Mjday_TDB(MJD_TT);
            JD_TDB = MJD_TDB+2400000.5;
            [r_Moon_a6,r_Sun_a6,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
            acc6 = Accel_EarthGravity_Tides(Mjd_UTC_a6,r_Sun_a6,r_Moon_a6,...
                r_ECI_a6,Cnm,Snm,eopdata);

            a7 = 5.0 / 6.0;
            Mjd_UTC_a7 = Mjd_UTC_Left + a7*h/86400;
            r_ECI_a7=[LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(1,:),Mjd_UTC_a7);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(2,:),Mjd_UTC_a7);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(3,:),Mjd_UTC_a7);];
            v_ECI_a7=[LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(1,:),Mjd_UTC_a7);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(2,:),Mjd_UTC_a7);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(3,:),Mjd_UTC_a7);];
            [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a7);
            [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT  = Mjd_UTC_a7 + TT_UTC/86400;
            MJD_TDB = Mjday_TDB(MJD_TT);
            JD_TDB = MJD_TDB+2400000.5;
            [r_Moon_a7,r_Sun_a7,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
            acc7 = Accel_EarthGravity_Tides(Mjd_UTC_a7,r_Sun_a7,r_Moon_a7,...
                r_ECI_a7,Cnm,Snm,eopdata);

            a8 = 1.0 / 6.0;
            Mjd_UTC_a8 = Mjd_UTC_Left + a8*h/86400;
            r_ECI_a8=[LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(1,:),Mjd_UTC_a8);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(2,:),Mjd_UTC_a8);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(3,:),Mjd_UTC_a8);];
            v_ECI_a8=[LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(1,:),Mjd_UTC_a8);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(2,:),Mjd_UTC_a8);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(3,:),Mjd_UTC_a8);];
            [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a8);
            [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT  = Mjd_UTC_a8 + TT_UTC/86400;
            MJD_TDB = Mjday_TDB(MJD_TT);
            JD_TDB = MJD_TDB+2400000.5;
            [r_Moon_a8,r_Sun_a8,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
            acc8 = Accel_EarthGravity_Tides(Mjd_UTC_a8,r_Sun_a8,r_Moon_a8,...
                r_ECI_a8,Cnm,Snm,eopdata);

            a9 = 2.0 / 3.0;
            Mjd_UTC_a9 = Mjd_UTC_Left + a9*h/86400;
            r_ECI_a9=[LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(1,:),Mjd_UTC_a9);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(2,:),Mjd_UTC_a9);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(3,:),Mjd_UTC_a9);];
            v_ECI_a9=[LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(1,:),Mjd_UTC_a9);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(2,:),Mjd_UTC_a9);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(3,:),Mjd_UTC_a9);];
            [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a9);
            [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT  = Mjd_UTC_a9 + TT_UTC/86400;
            MJD_TDB = Mjday_TDB(MJD_TT);
            JD_TDB = MJD_TDB+2400000.5;
            [r_Moon_a9,r_Sun_a9,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
            acc9 = Accel_EarthGravity_Tides(Mjd_UTC_a9,r_Sun_a9,r_Moon_a9,...
                r_ECI_a9,Cnm,Snm,eopdata);

            a10 = 1.0 / 3.0;
            Mjd_UTC_a10 = Mjd_UTC_Left + a10*h/86400;
            r_ECI_a10=[LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(1,:),Mjd_UTC_a10);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(2,:),Mjd_UTC_a10);
                LagrangeInterpolation(Mjd_UTC_Ind,r_ECI_Ind(3,:),Mjd_UTC_a10);];
            v_ECI_a10=[LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(1,:),Mjd_UTC_a10);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(2,:),Mjd_UTC_a10);
                LagrangeInterpolation(Mjd_UTC_Ind,v_ECI_Ind(3,:),Mjd_UTC_a10);];
            [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a10);
            [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT  = Mjd_UTC_a10 + TT_UTC/86400;
            MJD_TDB = Mjday_TDB(MJD_TT);
            JD_TDB = MJD_TDB+2400000.5;
            [r_Moon_a10,r_Sun_a10,~] = JPL_Eph_DE440(JD_TDB,DE440Coeff);
            acc10 = Accel_EarthGravity_Tides(Mjd_UTC_a10,r_Sun_a10,r_Moon_a10,...
                r_ECI_a10,Cnm,Snm,eopdata);

            c_1_11 = 41.0 / 840.0;
            c6 = 34.0 / 105.0;
            c_7_8= 9.0 / 35.0;
            c_9_10 = 9.0 / 280.0;
            Work_c = h * (c_1_11 * (dot(acc1,v_ECI_a1) + dot(acc11,v_ECI_a11)) +...
                c6 * dot(acc6,v_ECI_a6) + c_7_8 * (dot(acc7,v_ECI_a7) + dot(acc8,v_ECI_a8))...
                + c_9_10 * (dot(acc9,v_ECI_a9) + dot(acc10,v_ECI_a10)));
            Work = Work + Work_c;

        end
        Delta_E_Cur = (DeltaE_K - DeltaE_Sun - DeltaE_Moon - Work);
        Delta_E(cnt_data) = -Delta_E_Cur;

    end
    LLA = ecef2lla(r_ECEF);
    LST = mod(hours(Datetime_UTC+hours(LLA(:,2)/15) ...
        - dateshift(Datetime_UTC+hours(LLA(:,2)/15),"start","day")),24);


    s = struct("Datevec_UTC",Datevec_UTC,...
        "OE_Mean",OE_Mean,"OE_Osc",OE_Osc,"r_ECI",r_ECI,...
        "v_ECI",v_ECI,"r_ECEF",r_ECEF,"v_ECEF",v_ECEF, ...
        "LLA",LLA,"LST",LST,"Delta_E",Delta_E,"Cov",Cov);
    outputDir = replace(files(cnt_file).folder,inputPath,outputPath);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputname = regexp(files(cnt_file).name, '(\d+_STARLINK-\d+)', 'tokens');
    if isempty(outputname)
        continue
    end

    save(fullfile(outputDir,sprintf("%s.mat",outputname{1}{1})),"-fromstruct",s);
    fprintf("%d/%d %s time:%.2fs\n",cnt_file,MAX_COUNT_FILE,outputname{1}{1},toc)
end

end



function q=LagrangeInterpolation(rs,ry,x) %rs stand for the x-node, ry for the y-nodes
rs = rs(:)';ry = ry(:)';
mlocx=rs'*ones(1,length(rs));
msave=mlocx;
mloci=mlocx;
mlocx=-mlocx+x;
mlocx=mlocx-diag(diag(mlocx))+diag(ones(1,length(rs)));
mloci=-mloci+msave';
mloci=mloci-diag(diag(mloci))+diag(ones(1,length(rs)));
px=prod(mlocx);
pi=prod(mloci);
polyvect=px./pi;
q=dot(ry,polyvect);
end

function result = interpolateRows(time_nodes, data_rows, target_time)
num_rows = size(data_rows, 1);
result = zeros(num_rows, 1);
for i = 1:num_rows
    result(i) = LagrangeInterpolation(time_nodes, data_rows(i,:), target_time);
end
end