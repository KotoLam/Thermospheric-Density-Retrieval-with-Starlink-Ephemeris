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

addpath("Data")
addpath("Function")

eopdata = ReadEOPData(processDate);
load("DE440_Coeff.mat","DE440Coeff")
load("EGM2008_Coeff.mat","Cnm","Snm");

if isunix && ~ismac
    poolobj = gcp("nocreate");
    if isempty(poolobj)
        parpool(16)
    end
end

files = dir(fullfile(inputPath,sprintf("%s*",processDate),"*.mat"));

MAX_COUNT_FILE = length(files);
fprintf("Files Counts: %d\n",MAX_COUNT_FILE)
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
    mjd_UTC = mjuliandate(datevec_UTC);
    r_ECI = pos_ECI(1:MAX_COUNT_TIME,:)*1e3; % km -> m
    v_ECI = vel_ECI(1:MAX_COUNT_TIME,:)*1e3; % km/s -> m/s
    cov = cov(1:MAX_COUNT_TIME,:);

    deltaE = nan(MAX_COUNT_TIME,1);
    parfor cnt_data=1:MAX_COUNT_TIME-1

        if cnt_data < 5
            idx_Origin = 1:8;
        elseif cnt_data > MAX_COUNT_TIME-5
            idx_Origin = MAX_COUNT_TIME-7:MAX_COUNT_TIME;
        else
            idx_Origin = cnt_data-3:cnt_data+4;
        end
        mjd_Ind = mjd_UTC(idx_Origin);
        r_ECI_Ind = r_ECI(idx_Origin,:)';
        v_ECI_Ind = v_ECI(idx_Origin,:)';

        h = 60; % [s]
        a_i = [0, 1, 1.0/2.0, 5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0];
        mjd_i = mjd_UTC(cnt_data) + a_i.*h/86400;
        acc_i = nan(3,length(a_i));
        r_ECI_i = nan(3,length(a_i));
        v_ECI_i = nan(3,length(a_i));
        r_Moon_i = nan(3,length(a_i));
        r_Sun_i = nan(3,length(a_i));

        for cnt_Point = 1:length(a_i)
            mjd_Point = mjd_i(cnt_Point);
            r_ECI_Point = InterpolateRows(mjd_Ind, r_ECI_Ind, mjd_Point);
            v_ECI_Point = InterpolateRows(mjd_Ind, v_ECI_Ind, mjd_Point);
            r_ECI_i(:,cnt_Point) = r_ECI_Point;
            v_ECI_i(:,cnt_Point) = v_ECI_Point;
            idx = find(floor(mjd_Point)==eopdata(4,:),1,'first');
            preeop = eopdata(:,idx);
            TT_TAI  = +32.184;          % TT-TAI time difference [s]
            TAI_UTC = preeop(13);       % TAI-UTC time difference [s]
            TT_UTC  = TT_TAI+TAI_UTC;   % TT-UTC time difference [s]
            mjd_TT  = mjd_Point + TT_UTC/86400;
            mjd_TDB = Mjday_TDB(mjd_TT);
            jd_TDB = mjd_TDB+2400000.5;
            [r_Moon_Point,r_Sun_Point] = JPL_Eph_DE440(jd_TDB,DE440Coeff);
            r_Moon_i(:,cnt_Point) = r_Moon_Point;
            r_Sun_i(:,cnt_Point) = r_Sun_Point;
            acc_i(:,cnt_Point) = Accel_EarthGravity_Tides(mjd_Point,r_Sun_Point,r_Moon_Point,...
                r_ECI_Point,Cnm,Snm,eopdata);
        end

        c_1_11 = 41.0 / 840.0;
        c6 = 34.0 / 105.0;
        c_7_8= 9.0 / 35.0;
        c_9_10 = 9.0 / 280.0;

        Work = h * (c_1_11 * (dot(acc_i(:,1),v_ECI_i(:,1)) + dot(acc_i(:,2),v_ECI_i(:,2))) +...
            c6 * dot(acc_i(:,3),v_ECI_i(:,3)) + c_7_8 * (dot(acc_i(:,4),v_ECI_i(:,4)) + dot(acc_i(:,5),v_ECI_i(:,5)))...
            + c_9_10 * (dot(acc_i(:,6),v_ECI_i(:,6)) + dot(acc_i(:,7),v_ECI_i(:,7))));

        deltaE(cnt_data) = -(...
            ( 0.5*norm(v_ECI_i(:,2))^2 - (0.5*norm(v_ECI_i(:,1))^2) ) -...
            ( Potential_Moon(r_ECI_i(:,2),r_Moon_i(:,2)) - Potential_Moon(r_ECI_i(:,1),r_Moon_i(:,1)) ) -...
            ( Potential_Sun(r_ECI_i(:,2),r_Sun_i(:,2)) - Potential_Sun(r_ECI_i(:,1),r_Sun_i(:,1)) ) -...
            Work );
    end

    s = struct("Datevec_UTC",datevec_UTC,...
        "r_ECI",r_ECI,"v_ECI",v_ECI,"deltaE",deltaE,"Cov",cov);
    outputDir = replace(files(cnt_file).folder,inputPath,outputPath);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputname = regexp(files(cnt_file).name, '(\d+_STARLINK-\d+)', 'tokens');
    if isempty(outputname)
        continue
    end
    % save(fullfile(outputDir,sprintf("%s.mat",outputname{1}{1})),"-fromstruct",s);
    save(fullfile(outputDir,files(cnt_file).name),"-fromstruct",s);
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
