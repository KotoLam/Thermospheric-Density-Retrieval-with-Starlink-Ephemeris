addpath("Function")
% 设置起始和结束日期
start_date = datetime('2024-10-14');
end_date = datetime('2024-12-31');

% 输入输出目录
input_dir = "/home/data2/iono/Starlink/Ephemeris";
output_dir = "/home/data2/iono/Starlink/Ephemeris_mat";

% 循环处理每一天
current_date = start_date;
while current_date <= end_date
    % 将日期格式化为'YYYY-MM-DD'字符串
    date_str = datestr(current_date, 'yyyy-mm-dd');
    
    % 显示正在处理的日期
    disp(['Processing date: ', date_str]);
    
    % 调用你的函数
    Ephemeris_txt2mat(date_str, input_dir, output_dir);
    
    % 增加一天
    current_date = current_date + days(1);
end

disp('All dates processed successfully!');