function Ephemeris_txt2mat(processDate, inputPath, outputPath)
fprintf("Start processing %s\n",processDate)

files = dir(fullfile(inputPath,sprintf("%s*",processDate),"*.txt"));

MAX_COUNT_FILE = length(files);
fprintf("%s Files Counts: %d\n",processDate,MAX_COUNT_FILE)
parfor cnt_file = 1:MAX_COUNT_FILE
    % tic;
    fname = fullfile(files(cnt_file).folder,files(cnt_file).name);
    fid = fopen(fname);
    MAX_COUNT_TIME = 3*24*60+1; % 3天
    date_str = strings(MAX_COUNT_TIME,1);
    pos_ECI = nan(MAX_COUNT_TIME,3);
    vel_ECI = nan(MAX_COUNT_TIME,3);
    cov = nan(MAX_COUNT_TIME,21);
    epochCount=1;
    flag = 0;
    while(~feof(fid))
        line_t = fgetl(fid);
        if ~ischar(line_t), break, end
        if flag == 1
            A = strsplit(line_t);
            date_str(epochCount) = A{1};
            pos_ECI(epochCount,:) = str2double(A(2:4));
            vel_ECI(epochCount,:) = str2double(A(5:7));
            line_t = fgetl(fid);
            C1 = strsplit(line_t);
            cov(epochCount,1:7) = str2double(C1(1:7));
            line_t = fgetl(fid);
            C2 = strsplit(line_t);
            cov(epochCount,8:14) = str2double(C2(1:7));
            line_t = fgetl(fid);
            C3 = strsplit(line_t);
            cov(epochCount,15:21) = str2double(C3(1:7));

            epochCount = epochCount + 1;
        elseif strcmp(line_t(1:3),'UVW')
            flag = 1;
        else
            continue
        end
    end
    fclose(fid);
    % 由于数据源采用年积日doy，所以用datetime方便转换
    datetime_UTC = datetime(date_str,"Format","uuuuDDDHHmmss.SSS");
    datevec_UTC = datevec(datetime_UTC);
    s = struct("datevec_UTC",datevec_UTC,"pos_ECI",pos_ECI,"vel_ECI",vel_ECI,"cov",cov);
    
    outputDir = replace(files(cnt_file).folder,inputPath,outputPath);
    outputname = replace(files(cnt_file).name, ".txt", ".mat");
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir,outputname),"-fromstruct",s);
    % fprintf("%d/%d %s time:%.2fs\n",cnt_file,MAX_COUNT_FILE,outputname,toc)

end
end

% Ephemeris_txt2mat("2025-06-29", "/home/data2/iono/Starlink/Ephemeris", "/home/data2/iono/Starlink/Ephemeris_mat")