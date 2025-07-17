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

files = dir(fullfile(inputPath,sprintf("%s*",processDate),"*.mat"));

MAX_COUNT_FILE = length(files);
fprintf("files counts: %d\n",MAX_COUNT_FILE)
for cnt_file = 1:MAX_COUNT_FILE
    tic;
    fname = fullfile(files(cnt_file).folder,files(cnt_file).name);
    s = load(fname);
    datevec_UTC = s.datevec_UTC;
    pos_ECI = s.pos_ECI;
    vel_ECI = s.vel_ECI;
    cov = s.cov;

    MAX_COUNT_TIME = 2*24*60+1; % 2天

    datevec_UTC = datevec_UTC(1:MAX_COUNT_TIME,:);
    datetime_UTC = datetime(datevec_UTC);
    mjd_UTC = mjuliandate(datevec_UTC);
    r_ECI = pos_ECI(1:MAX_COUNT_TIME,:)*1e3; % km -> m
    v_ECI = vel_ECI(1:MAX_COUNT_TIME,:)*1e3; % km/s -> m/s
    cov = cov(1:MAX_COUNT_TIME,:);

    r_ECEF = nan(size(r_ECI));
    v_ECEF = nan(size(v_ECI));
    deltaE = nan(MAX_COUNT_TIME,1);

    for cnt_data=1:MAX_COUNT_TIME-1

        TT_TAI  = +32.184;          % TT-TAI time difference [s]
        %
        mjd_Cur = mjd_UTC(cnt_data);
        r_ECI_Cur = r_ECI(cnt_data,:);
        v_ECI_Cur = v_ECI(cnt_data,:);

        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,mjd_Cur);


        TT_UTC  = TT_TAI+TAI_UTC;   % TT-UTC time difference [s]
        mjd_TT  = mjd_Cur + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_Cur,r_Sun_Cur] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
        %
        mjd_Next = mjd_UTC(cnt_data+1);
        r_ECI_Next = r_ECI(cnt_data+1,:);
        v_ECI_Next = v_ECI(cnt_data+1,:);
        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,mjd_Next);
        TT_UTC  = TT_TAI+TAI_UTC;   % TT-UTC time difference [s]
        mjd_TT  = mjd_Next + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_Next,r_Sun_Next] = JPL_Eph_DE440(jd_TDB,DE440Coeff);

        V_K_Cur = 0.5*norm(v_ECI_Cur)^2;
        V_K_Next = 0.5*norm(v_ECI_Next)^2;
        DeltaE_K = V_K_Next - V_K_Cur;

        V_Moon_Cur = Potential_Moon(r_ECI_Cur,r_Moon_Cur);
        V_Moon_Next = Potential_Moon(r_ECI_Next,r_Moon_Next);
        DeltaE_Moon = V_Moon_Next - V_Moon_Cur;

        V_Sun_Cur = Potential_Sun(r_ECI_Cur,r_Sun_Cur);
        V_Sun_Next = Potential_Sun(r_ECI_Next,r_Sun_Next);
        DeltaE_Sun = V_Sun_Next - V_Sun_Cur;

        if cnt_data < 5
            idx_Originn = 1:8;
        elseif cnt_data > MAX_COUNT_TIME-5
            idx_Originn = MAX_COUNT_TIME-7:MAX_COUNT_TIME;
        else
            idx_Originn = cnt_data-3:cnt_data+4;
        end
        mjd_Ind = mjd_UTC(idx_Originn);
        r_ECI_Ind = r_ECI(idx_Originn,:)';
        v_ECI_Ind = v_ECI(idx_Originn,:)';
        
        h = mjd_Next - mjd_Cur;
        a_i = [0, 1, 1.0/2.0, 5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0];
        mjd_i = mjd_Cur + a_i.*h;
        acc_i = nan(3,length(a_i));
        r_ECI_i = nan(3,length(a_i));
        v_ECI_i = nan(3,length(a_i));
        
        for cnt_Point = 1:length(a_i)
            mjd_Point = mjd_i(cnt_Point);
            r_ECI_Point = InterpolateRows(mjd_Ind, r_ECI_Ind, mjd_Point);
            v_ECI_Point = InterpolateRows(mjd_Ind, v_ECI_Ind, mjd_Point);
            idx = find(floor(mjd_Point)==eopdata(4,:),1,'first');
            preeop = eopdata(:,idx);
            TAI_UTC = preeop(13);         % TAI-UTC time difference [s]
            TT_UTC  = TT_TAI+TAI_UTC;   % TT-UTC time difference [s]
            mjd_TT  = mjd_Point + TT_UTC/86400;
            mjd_TDB = Mjday_TDB(mjd_TT);
            jd_TDB = mjd_TDB+2400000.5;
            [r_Moon_Point,r_Sun_Point] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
            acc_i(:,cnt_Point) = Accel_EarthGravity_Tides(mjd_Point,r_Sun_Point,r_Moon_Point,...
            r_ECI_Point,Cnm,Snm,eopdata);
        end

        r_ECI_a1=[LagrangeInterpolation(mjd_Ind,r_ECI_Ind(1,:),Mjd_UTC_a1);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(2,:),Mjd_UTC_a1);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(3,:),Mjd_UTC_a1);];
        v_ECI_a1=[LagrangeInterpolation(mjd_Ind,v_ECI_Ind(1,:),Mjd_UTC_a1);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(2,:),Mjd_UTC_a1);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(3,:),Mjd_UTC_a1);];
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a1);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        mjd_TT  = Mjd_UTC_a1 + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_a1,r_Sun_a1,~] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
        acc1 = Accel_EarthGravity_Tides(Mjd_UTC_a1,r_Sun_a1,r_Moon_a1,...
            r_ECI_a1,Cnm,Snm,eopdata);

        a11 = 1;
        Mjd_UTC_a11 = Mjd_UTC_Left + a11*h/86400;
        r_ECI_a11=[LagrangeInterpolation(mjd_Ind,r_ECI_Ind(1,:),Mjd_UTC_a11);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(2,:),Mjd_UTC_a11);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(3,:),Mjd_UTC_a11);];
        v_ECI_a11=[LagrangeInterpolation(mjd_Ind,v_ECI_Ind(1,:),Mjd_UTC_a11);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(2,:),Mjd_UTC_a11);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(3,:),Mjd_UTC_a11);];
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a11);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        mjd_TT  = Mjd_UTC_a11 + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_a11,r_Sun_a11,~] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
        acc11 = Accel_EarthGravity_Tides(Mjd_UTC_a11,r_Sun_a11,r_Moon_a11,...
            r_ECI_a11,Cnm,Snm,eopdata);

        a6 = 1.0 / 2.0;
        Mjd_UTC_a6 = Mjd_UTC_Left + a6*h/86400;
        r_ECI_a6=[LagrangeInterpolation(mjd_Ind,r_ECI_Ind(1,:),Mjd_UTC_a6);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(2,:),Mjd_UTC_a6);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(3,:),Mjd_UTC_a6);];
        v_ECI_a6=[LagrangeInterpolation(mjd_Ind,v_ECI_Ind(1,:),Mjd_UTC_a6);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(2,:),Mjd_UTC_a6);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(3,:),Mjd_UTC_a6);];
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a6);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        mjd_TT  = Mjd_UTC_a6 + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_a6,r_Sun_a6,~] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
        acc6 = Accel_EarthGravity_Tides(Mjd_UTC_a6,r_Sun_a6,r_Moon_a6,...
            r_ECI_a6,Cnm,Snm,eopdata);

        a7 = 5.0 / 6.0;
        Mjd_UTC_a7 = Mjd_UTC_Left + a7*h/86400;
        r_ECI_a7=[LagrangeInterpolation(mjd_Ind,r_ECI_Ind(1,:),Mjd_UTC_a7);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(2,:),Mjd_UTC_a7);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(3,:),Mjd_UTC_a7);];
        v_ECI_a7=[LagrangeInterpolation(mjd_Ind,v_ECI_Ind(1,:),Mjd_UTC_a7);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(2,:),Mjd_UTC_a7);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(3,:),Mjd_UTC_a7);];
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a7);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        mjd_TT  = Mjd_UTC_a7 + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_a7,r_Sun_a7,~] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
        acc7 = Accel_EarthGravity_Tides(Mjd_UTC_a7,r_Sun_a7,r_Moon_a7,...
            r_ECI_a7,Cnm,Snm,eopdata);

        a8 = 1.0 / 6.0;
        Mjd_UTC_a8 = Mjd_UTC_Left + a8*h/86400;
        r_ECI_a8=[LagrangeInterpolation(mjd_Ind,r_ECI_Ind(1,:),Mjd_UTC_a8);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(2,:),Mjd_UTC_a8);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(3,:),Mjd_UTC_a8);];
        v_ECI_a8=[LagrangeInterpolation(mjd_Ind,v_ECI_Ind(1,:),Mjd_UTC_a8);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(2,:),Mjd_UTC_a8);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(3,:),Mjd_UTC_a8);];
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a8);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        mjd_TT  = Mjd_UTC_a8 + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_a8,r_Sun_a8,~] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
        acc8 = Accel_EarthGravity_Tides(Mjd_UTC_a8,r_Sun_a8,r_Moon_a8,...
            r_ECI_a8,Cnm,Snm,eopdata);

        a9 = 2.0 / 3.0;
        Mjd_UTC_a9 = Mjd_UTC_Left + a9*h/86400;
        r_ECI_a9=[LagrangeInterpolation(mjd_Ind,r_ECI_Ind(1,:),Mjd_UTC_a9);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(2,:),Mjd_UTC_a9);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(3,:),Mjd_UTC_a9);];
        v_ECI_a9=[LagrangeInterpolation(mjd_Ind,v_ECI_Ind(1,:),Mjd_UTC_a9);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(2,:),Mjd_UTC_a9);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(3,:),Mjd_UTC_a9);];
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a9);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        mjd_TT  = Mjd_UTC_a9 + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_a9,r_Sun_a9,~] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
        acc9 = Accel_EarthGravity_Tides(Mjd_UTC_a9,r_Sun_a9,r_Moon_a9,...
            r_ECI_a9,Cnm,Snm,eopdata);

        a10 = 1.0 / 3.0;
        Mjd_UTC_a10 = Mjd_UTC_Left + a10*h/86400;
        r_ECI_a10=[LagrangeInterpolation(mjd_Ind,r_ECI_Ind(1,:),Mjd_UTC_a10);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(2,:),Mjd_UTC_a10);
            LagrangeInterpolation(mjd_Ind,r_ECI_Ind(3,:),Mjd_UTC_a10);];
        v_ECI_a10=[LagrangeInterpolation(mjd_Ind,v_ECI_Ind(1,:),Mjd_UTC_a10);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(2,:),Mjd_UTC_a10);
            LagrangeInterpolation(mjd_Ind,v_ECI_Ind(3,:),Mjd_UTC_a10);];
        [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(eopdata,Mjd_UTC_a10);
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        mjd_TT  = Mjd_UTC_a10 + TT_UTC/86400;
        mjd_TDB = Mjday_TDB(mjd_TT);
        jd_TDB = mjd_TDB+2400000.5;
        [r_Moon_a10,r_Sun_a10,~] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
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

        Delta_E_Cur = (DeltaE_K - DeltaE_Sun - DeltaE_Moon - Work);
        deltaE(cnt_data) = -Delta_E_Cur;

    end
    LLA = ecef2lla(r_ECEF);
    LST = mod(hours(datetime_UTC+hours(LLA(:,2)/15) ...
        - dateshift(datetime_UTC+hours(LLA(:,2)/15),"start","day")),24);


    s = struct("Datevec_UTC",datevec_UTC,...
        "OE_Mean",OE_Mean,"OE_Osc",OE_Osc,"r_ECI",r_ECI,...
        "v_ECI",v_ECI,"r_ECEF",r_ECEF,"v_ECEF",v_ECEF, ...
        "LLA",LLA,"LST",LST,"Delta_E",deltaE,"Cov",cov);
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



function q=LagrangeInterpolation(rs,ry,x) 
% rs stand for the x-node, ry for the y-nodes
% 
% refer to Vincent Naudot under the answer of 
% https://ww2.mathworks.cn/matlabcentral/answers/305169-what-is-the-code-for-lagrange-interpolating-polynomial-for-a-set-of-given-data#answer_500110
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

function result = InterpolateRows(time_nodes, data_rows, target_time)
num_rows = size(data_rows, 1);
result = zeros(num_rows, 1);
for i = 1:num_rows
    result(i) = LagrangeInterpolation(time_nodes, data_rows(i,:), target_time);
end
end