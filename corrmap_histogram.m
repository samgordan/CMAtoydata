function [handles, histogram] = corrmap_histogram(corrmap, threshold, normalization)
%[handles, histogram] = corrmap_histogram(corrmap, threshold, normalization)
%
%
% Compute the histogram of the correlation map using the specified
% threshold for correlation strength.  Threshold must be a real number 0 < threshold < 1
% The normalization option can be either 'length','hits', or false.

% Written by Samantha Gordon Danner (2018), based on original script by Rob
% Fuhrman (2014)

%% INPUT VALIDATION
if ~isa(corrmap,'CorrelationMap'); 
    error('Input corrmap must be an object of type: CorrelationMap'); 
end;

if nargin < 3
    error('This function requires three inputs: corrmap, threshold, and normalization');
end

if ~(strcmpi(normalization,'hits') || strcmpi(normalization,'length') || (normalization == false))
    error('Normalization must either be false or one of the options ''hits'' or ''length''');
end

%% MAIN
delta = corrmap.delta;
corr_map = corrmap.corr_map;
fs = corrmap.rate;

%N corrmap samples for normalization
N = length(corr_map);

pos_map = corr_map;
neg_map = corr_map;

%Map thresholding

%positive correlations
pos_map(pos_map > threshold) = 1;
pos_map(pos_map < threshold) = 0;

%negative correlations
neg_map(neg_map > -threshold) = 0;
neg_map(neg_map < -threshold) = 1;

%sums
pos_sum = sum(pos_map,2);
neg_sum = sum(neg_map,2);

%Nhits for normalization
Npos = sum(pos_sum);
Nneg = sum(neg_sum);


%vector of offsets for plotting
offset_vector = [-delta/fs:1/fs:delta/fs];
offset_vector = fliplr(offset_vector);

%create an vector of attenuation lengths for each offset based on delta
attenuation = [[delta:-1:1]';0;[1:delta]']*2;


%switch-case for normalization option
switch normalization
    
    case 'hits'
        p_hist = pos_sum / Npos; 
        n_hist = neg_sum / Nneg;
        histogram.pos = p_hist;
        histogram.neg = n_hist;
        
    case 'length'
        lengths = N * ones(size(corr_map,1),1) - attenuation;
        p_hist = bsxfun(@rdivide,pos_sum,lengths);
        n_hist = bsxfun(@rdivide,neg_sum,lengths);
        histogram.pos = p_hist;
        histogram.neg = n_hist;
        
    case false
        histogram.pos = pos_sum;
        histogram.neg = neg_sum;      
        
end

histogram.plot_axis = offset_vector;

%compute the mean offset as an index
mh_pos = sum([1:2*delta+1]' .* histogram.pos) / sum(histogram.pos);
mh_neg = sum([1:2*delta+1]' .* histogram.neg) / sum(histogram.neg);

%compute the mode of the histogram as an index
[~,mode_pos] = max(histogram.pos);
[~,mode_neg] = max(histogram.neg);


histogram.pos_mean = round(mh_pos); %round() ensures that these are integers.
histogram.pos_mode = mode_pos;

histogram.neg_mean = round(mh_neg);
histogram.neg_mode = mode_neg;


%% PLOTTING
handles.pos_handle = figure; set(handles.pos_handle, 'Visible', 'off'); bar(offset_vector, histogram.pos,'r'); %added set(handles.pos_handle, 'Visible', 'off') 5/7/17 to prevent figure display when using function in loop
hold on; plot(histogram.plot_axis(histogram.pos_mode),histogram.pos(histogram.pos_mode),'kx','MarkerSize',15,'LineWidth',2);
xlabel('Offset');



handles.neg_handle = figure;set(handles.neg_handle, 'Visible', 'off'); bar(offset_vector, histogram.neg,'b'); %added set(handles.neg_handle, 'Visible', 'off') 5/7/17 to prevent figure display when using function in loop
hold on; plot(histogram.plot_axis(histogram.neg_mode),histogram.neg(histogram.neg_mode),'kx','MarkerSize',15,'LineWidth',2);
xlabel('Offset');

end

