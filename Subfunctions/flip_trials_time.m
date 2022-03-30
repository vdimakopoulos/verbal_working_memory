function flipped_data = flip_trials_time(data)

nTrials = length(data.trial)'
flipped_data = data;
for i = 1:nTrials
    flipped_data.trial{i} = fliplr(data.trial{i})
end

end