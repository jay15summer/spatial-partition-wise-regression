%% greedy merging of 2D spatial segmented piecewise linear regression
%% (partition wise linear regression)
%%=========================================================================
function [partition_all, partiaion_slt] = spatial_partition_reg(S, X, y, h, v, T)
%%
% Input: -S (N-by-2)is the spatial location of the data
%        -X (N-by-m)are the independent variables at location S
%        -y (N-by-1)is the response at location S
%        -h initializes the segments in horizontal direction
%        -v initializes the segments in vertical direction
%        (initial partition is h*v rectangle grids )# can be improved in
%        future release, i.e., more flexible partition
%        -T is a threshold (stopping criteria) applied to the percentage
%        increase of merging errors, e.g., (SSE^i-SSE^{i-1})/SSE^{i-1}; T = 0.05.
% Output:-partition_all is a cell recording the segmentations of every
%        iterations (partition_all.segments) and their fitting errors
%        (partition_all.metric).
%        -partition_slt is a map representing the selected partitions/segmentations;
%        the keys of the map can be obtained via 'cell2mat(keys(partition_slt))';
%        the data [S X y] of every partition can be obtained via 'partition_slt(keys).data';
%        the neighbor of the selected segmentation is 'partition_slt(keys).neighbor'.
%        -two plots: 1. percentage-increase of SSE 2. selected final partitioning.
% For every segmentions, the relationship of X-y is different,
% i.e., the coefficient \beta of y = \beta*X + noise is piecewise.
%%
% 1. initialize the segmentation (this program uses map to represent partitioning)
% ini_partition is a map with keys representing different segmentations
% ini_partition(keys).data is the [S X y] data for every segmentation
% ini_partition(keys).neighbor is keys' neighboring segments/partitions (neighboring keys)
ini_partition = initialize_segment(S, X, y, h, v);
% 2. calculate the metric: SSE for the initial segmentation
ini_metric = fit_evaluate_all(ini_partition);
% 3. iterative merge until only one segmentation is left or meeting the
%% stopping criteria
partition = struct();
partition.segments = ini_partition; % segments is the map of the partitioning/segments
partition.metric = ini_metric; % metric is the SSE for the partitioning
partition_all = cell(h*v, 1);
partition_all{1} = partition;

increase_percent = 0; ite = 1;
while ite <= h*v-1 && increase_percent <= T
    partition = struct();
    [partition.segments, metric_increase] = best_merge(partition_all{ite}.segments);
    partition.metric = partition_all{ite}.metric + metric_increase;
    partition_all{ite+1} = partition;
    increase_percent = metric_increase / partition_all{ite}.metric;
    fprintf('Merge iteration %d of possible maximum iteration %d ... \n', ite, h*v-1)
    ite = ite + 1;
end

% 4. select partition/segmentation, and plot the segmentation
if ite == h*v && increase_percent <= T
    slt = ite; % only one segmentation
    fprintf('Merge finished! Not exceeding threshold T! 1 segmentation is found!\n')
elseif ite == h*v && increase_percent > T
    slt = ite - 1; % two segmentations
    fprintf('Merge finished! Exceed threshold T! 2 segmentations are found!\n')
else
    slt = ite - 1; % select the one before exceeding threshold T
    fprintf('Merge finished early! Exceed threshold T! %d segmentations are found!\n', h*v - slt + 1)
end
partiaion_slt = partition_all{slt}.segments;
Seg= []; color = 1;
for i = cell2mat(keys(partiaion_slt))
    S_tmp = partiaion_slt(i).data(:,1:2);
    c_value_tmp = color*ones(size(S_tmp, 1), 1);
    S_tmp = [S_tmp c_value_tmp];
    Seg = [Seg; S_tmp];
    color = color + 1;
end
figure;
scatter(Seg(:,1), Seg(:,2),3, Seg(:,3),'square','filled'); axis equal; colorbar;
title('spatial partitioning/segmentation');
xlabel('S1'); ylabel('S2');

% 5. plot the percentage of increase
idx_no_empty = find(~cellfun(@isempty, partition_all));
percent = zeros(length(idx_no_empty)-1, 1);
for i = 2:length(idx_no_empty)
    increase = partition_all{i}.metric - partition_all{i-1}.metric;
    percent(i-1) = increase/partition_all{i-1}.metric*100;
end
figure;
plot(percent)
title('percentage of the increase');
xlabel('merge iteration number'); ylabel('increase-percentage');

end

%%========================sub-functions====================================
function [new_data_seg, metric_best] = best_merge(old_data_seg)
%%
% go through all the neighbor combinations in old_data_seg and select the
% combination with the best metric, e.g., least SSE,
% merge the selected combination, return new_data_seg.
%%
metric_increase = NaN + zeros(2*length(old_data_seg), 1);
merge_track_tmp = NaN + zeros(2*length(old_data_seg), 2);
ite = 1;
for idx_merge1 = cell2mat(keys(old_data_seg))
    for j_merge = 1:length(old_data_seg(idx_merge1).neighbor)
        idx_merge2 = old_data_seg(idx_merge1).neighbor(j_merge);
        data1 = old_data_seg(idx_merge1).data;
        data2 = old_data_seg(idx_merge2).data;
        data_comb = [data1; data2];
        % a pair of neighbors fitting independently error
        metric_ind = fit_evaluate_part(data1(:, 3:end-1), data1(:, end)) + fit_evaluate_part(data2(:, 3:end-1), data2(:, end));
        % a pair of neighbors merged fitting error
        metric_merge = fit_evaluate_part(data_comb(:, 3:end-1), data_comb(:, end));
        % calculate metric increase
        metric_increase(ite) = metric_merge - metric_ind;

        merge_track_tmp(ite,:) = [idx_merge1 idx_merge2];
        ite = ite + 1;
    end
end
% find the merge combination with the best metric evaluation
% metric_best = max(metric)
metric_best = min(metric_increase);
merged = merge_track_tmp(metric_increase == metric_best, :);
merge_track = merged(1, :);
% return the revised segmentation map with updated segmentation, data and
% neighbors
new_data_seg = copymap(old_data_seg);
struct_tmp = struct();
struct_tmp.data = [old_data_seg(merge_track(1)).data; old_data_seg(merge_track(2)).data];
nb_tmp = [old_data_seg(merge_track(1)).neighbor old_data_seg(merge_track(2)).neighbor];
nb_tmp = nb_tmp(nb_tmp ~= merge_track(1)& nb_tmp ~= merge_track(2));
nb_tmp = unique(nb_tmp);
struct_tmp.neighbor = nb_tmp;
% keep key merge_track(2) and remove key merge_track(1)
new_data_seg(merge_track(2)) = struct_tmp;
remove(new_data_seg, merge_track(1));
% update related segmentation's neighbor index due to the merge,
% i.e., change neighbor index of idx_merge1 to idx_merge2
for i_nb = cell2mat(keys(new_data_seg))
    nb_check = new_data_seg(i_nb).neighbor;
    if ismember(merge_track(1), nb_check)
        nb_check(nb_check == merge_track(1)) = merge_track(2);
        nb_check = unique(nb_check);
        struct_tmp = struct();
        struct_tmp.data = new_data_seg(i_nb).data;
        struct_tmp.neighbor = nb_check;
        new_data_seg(i_nb) = struct_tmp;
    end
end
% calculate the new metric for all the segmentations after merging
% metric_new = fit_evaluate_all(new_data_seg);
end



function ini_data_seg = initialize_segment(S, X, y, h, v)
%%
% make initial #h*v number of segmentations
%%
S1_low = min(S(:,1));
S1_high = max(S(:,1));
S2_low = min(S(:,2));
S2_high = max(S(:,2));
% segmentation is initialized by assuming a uniform partition:
% axial direction: horizontal: h segments; vertical: v segments
% so the organization of the index for the segmentations are e.g.,
% h = 4; v = 3;
% 9 10 11 12
% 5  6  7  8
% 1  2  3  4
h_segment = (S1_high-S1_low)/h*(1:h-1) + S1_low; h_segment = [S1_low h_segment S1_high];
v_segment = (S2_high-S2_low)/v*(1:v-1) + S2_low; v_segment = [S2_low v_segment S2_high];
% 2. allocate [S, X, y] to each segmentations
index_seg = 1:1:h*v;
data_by_index = cell(1, h*v);
for i = 1:h*v
    if mod(i, h)~= 0 && ceil(i/h) ~= v
        idx_tmp = find(S(:,1)>=h_segment(mod(i, h))&S(:,1)<h_segment(mod(i, h)+1)&S(:,2)>=v_segment(ceil(i/h))&S(:,2)<v_segment(ceil(i/h)+1));
        S_tmp = S(idx_tmp,:);
        X_tmp = X(idx_tmp,:);
        y_tmp = y(idx_tmp);
    elseif mod(i, h)~= 0 && ceil(i/h) == v
        idx_tmp = find(S(:,1)>=h_segment(mod(i, h))&S(:,1)<h_segment(mod(i, h)+1)&S(:,2)>=v_segment(ceil(i/h))&S(:,2)<=v_segment(ceil(i/h)+1));
        S_tmp = S(idx_tmp,:);
        X_tmp = X(idx_tmp,:);
        y_tmp = y(idx_tmp);
    elseif ceil(i/h) ~= v
        idx_tmp = find(S(:,1)>=h_segment(h)&S(:,1)<=h_segment(h+1)&S(:,2)>=v_segment(i/h)&S(:,2)<v_segment(i/h+1));
        S_tmp = S(idx_tmp,:);
        X_tmp = X(idx_tmp,:);
        y_tmp = y(idx_tmp);
    else
        idx_tmp = find(S(:,1)>=h_segment(h)&S(:,1)<=h_segment(h+1)&S(:,2)>=v_segment(i/h)&S(:,2)<=v_segment(i/h+1));
        S_tmp = S(idx_tmp,:);
        X_tmp = X(idx_tmp,:);
        y_tmp = y(idx_tmp);
    end

    if i == h*v
        neighbor_tmp = [];
    elseif mod(i,h) == 0
        neighbor_tmp = i+h;
    elseif ceil(i/h) == v
        neighbor_tmp = i+1;
    else
        neighbor_tmp = [i+1, i+h];
    end

    segment_data_tmp = struct();
    segment_data_tmp.data = [S_tmp X_tmp y_tmp];
    segment_data_tmp.neighbor = neighbor_tmp;
    data_by_index{i} = segment_data_tmp;
end
ini_data_seg = containers.Map(index_seg, data_by_index);
end


function eval_metric = fit_evaluate_all(data_seg_map)
% use sum of squared error as metric
% initial SSE
SSE = 0.;
for seg = cell2mat(keys(data_seg_map))
    % select [X, y] data in this segmentation;
    % exclude S
    X_seg = data_seg_map(seg).data(:, 3:end-1);
    y_seg = data_seg_map(seg).data(:, end);
    regf=@(XTRAIN,ytrain,XTEST)([XTEST ones(size(XTEST, 1), 1)] * ([XTRAIN ones(size(XTRAIN, 1), 1)]\ytrain));
    SSE = sum((y_seg - regf(X_seg, y_seg, X_seg)).^2) + SSE;

end
eval_metric = SSE;
end


function eval_metric = fit_evaluate_part(X_seg, y_seg)
% use sum of squared error as metric
regf=@(XTRAIN,ytrain,XTEST)([XTEST ones(size(XTEST, 1), 1)] * ([XTRAIN ones(size(XTRAIN, 1), 1)]\ytrain));
SSE = sum((y_seg - regf(X_seg, y_seg, X_seg)).^2);
eval_metric = SSE;
end


function copiedmap = copymap(map)
%%
% for copying the segmentation data map
%%
copy  = cell(1, length(map));
k = 1;
for i = cell2mat(keys(map))
    tmp_struct = struct();
    tmp_struct.data = map(i).data;
    tmp_struct.neighbor = map(i).neighbor;
    copy{k} = tmp_struct;
    k = k + 1;
end
copiedmap = containers.Map(cell2mat(keys(map)), copy);
end
