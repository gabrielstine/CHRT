function psych_struct = psychIntervalCalc(psych_struct)
% non_rt_val & non_rt_val include inter-trial intervals (~2sec feedback)

subject = psych_struct.subject;
epoch = psych_struct.epoch;

% file_name = [psych_struct.root_dir,'/SavedData/',subject,'Data',epoch];
% load(file_name);

data = psych_struct.data;


% data = [result', response', coh', rt', start_t', end_t', pts', max_dot_dur']

% R(strfind(str_result,'WRONG'))    = 0;
% R(strfind(str_result,'CORRECT'))  = 1;
% R(strfind(str_result,'NOFIX'))    = 2;
% R(strfind(str_result,'FIXBREAK')) = 3;
% R(strfind(str_result,'NOCHOICE')) = 4;
% R(strfind(str_result,'TIMEOUT'))  = 5;

% Calculated in "sec"

val_ind = logical(data(:,1)==0 | data(:,1)==1);
time_lim_ind = logical(data(:,1)==5);
trial_length = data(:,6) - data(:,5);
rt = data(:,4)/1000;
max_dot_dur = data(:,8);
interval = data(2:end,5) - data(1:(end-1),6);
interval = interval(interval>1);

non_rt_val = trial_length(val_ind) - rt(val_ind);
non_rt_inv = trial_length(time_lim_ind) - max_dot_dur(time_lim_ind);
psych_struct.non_rt_val = mean(non_rt_val)+mean(interval);
psych_struct.non_rt_inv = mean(non_rt_inv)+mean(interval);
