function [good_start, good_dur, n_good] = bad2good_segments(bad_start, bad_dur, EEG)
% bad2good_segments converts start points and durations of bad data segments to start points and durations of the corresponding good segments
%   Inputs: start_points_bad_segments, duration_bad_segments, EEG_dataset
%   Outputs: start_points_good_segments, duration_good_segments, number_of_good_epochs
% Author: Tilman Stephani, May 2018

    if size(bad_start,1)>0  %if there are any bad segments
        % number of "good" segments
        if (bad_start(end) + bad_dur(end) == size(EEG.data, 2) && bad_start(1)~=1) ...
                || (bad_start(end) + bad_dur(end) ~= size(EEG.data, 2) && bad_start(1)==1)
            n_good = length(bad_start);
        elseif bad_start(end) + bad_dur(end) == size(EEG.data, 2) && bad_start(1)==1
            n_good = length(bad_start) -1;
        else
            n_good = length(bad_start) +1;
        end

        % get start points and durations
        tmp_good_start = 1; % temporary variable
        for i = 1:n_good
            good_start(i) = tmp_good_start;
            if i < n_good
                good_dur(i) = bad_start(i)-1 - tmp_good_start;
                tmp_good_start = bad_start(i) + bad_dur(i) +1;
            elseif i == n_good
                good_dur(i) = size(EEG.data, 2) - tmp_good_start; % end of the data
            end
        end
        
    else %no bad segments
        good_start = 1;
        good_dur = size(EEG.data, 2);
        n_good = 1;
    end        
end

