function t = make_feasible(t, node)
% From a possibly non-feasible tree below node, iterate until feasible
% A tree is feasible when
%   1. groups have rates equal the sum of child rates
%   2. wells/groups have rates \leq individual limits and if
%      available for group control, \leq group its group target  
%   3. all nodes are categorized with 1 of 4 'mode_cnt' with the following
%      requirements:
%      INDIVIDUAL: 
%      Node is a well or group with a rate equal to individual limits
%      GROUP: 
%      Node is a well or group available for group control and target mode 
%      rate equal to group target
%      NONE: 
%      Node is group available for group control where none of its children 
%      are labeled with GROUP and all rates less than individual limits 
%      and group target
%      UNDETERMINED: 
%      Node is a group *not* available for group-control with all rates 
%      less than individual limits 
%                 
tol = 1e-3;
if nargin < 2
    node = 1; % FIELD-node
end
done = false;
count = 0;
% set group target of top node to nan -> starts with limits
target_first = true;
if any(isnan(t(node).rates))
    t(node).rates = [0 0 0];
else
    target_first = true;
end
[switched, rates_changed] = deal(false);
while ~done
    rates = t(node).rates;
    if target_first
        [t, switched] = cap_individual_and_sum(t, node, tol);
        rates_changed = any(abs(t(node).rates - rates) > tol*abs(rates));
        if any(t(node).rates < 0)
            % initial guess too far from feasible, revert to approx zero-rate
            % state
            t = initializeWithTinyRates(t, node, node);
        end
    else
        target_first = false;
    end
    [t, target_changed] = set_and_update_targets(t, node, tol);
    done = ~switched && ~rates_changed && ~target_changed;
    count = count +1;
end
if count == 1
    fprintf('make_feasible: no changes required from node %d\n', node)
else
    fprintf('make_feasible used %d iterations from node %d\n', count, node)
end
%plotGroupTree(t, 'rate');
end

%--------------------------------------------------------------------------
function t = initializeWithTinyRates(t, node, top_node)
% introduce small well rates such that all wells/groups
% operates under group-control (assuming limits are reasonnable)
if ~t(node).allow_group_control && node ~= top_node 
    % don't touch sub-trees not part of current balancing
    return
end
if strcmp(t(node).type, 'well')
    f_tmp = max(sqrt(eps), t(node).fractions);
    fac = 1e-6*min(t(node).limits./f_tmp);
    t(node).rates = fac*t(node).fractions;
else
    c = t(node).children;
    for k = 1:numel(c)
        t = initializeWithTinyRates(t, c(k), top_node);
    end
    t(node).rates = sum(vertcat(t(c).rates), 1);
end
end
%--------------------------------------------------------------------------
function [t, switched] = cap_individual_and_sum(t, node, tol)
% From well and up, cap at individual limits and sum
c = t(node).children;
any_group_controlled_children = false;
switched = false;
if ~isempty(c)
    [rates, group_rates] = deal(zeros(size(t(node).rates)));
    for k = 1:numel(c)
        [t, switched_child] = cap_individual_and_sum(t, c(k), tol);
        any_group_controlled_children = any_group_controlled_children || strcmp(t(c(k)).mode_cnt, 'GRUP');
        switched = switched || switched_child;  
        rates = rates + t(c(k)).rates;
        if ~strcmp(t(c(k)).mode_cnt, 'INDIVIDUAL')
           group_rates = group_rates + t(c(k)).group_rates;
        end
    end
    t(node).rates = rates;
    t(node).group_rates = group_rates;
else
    if strcmp(t(node).mode_cnt, 'INDIVIDUAL')
        t(node).group_rates = [0 0 0];
    else
        t(node).group_rates = t(node).rates;
    end
end
mode_cnt = t(node).mode_cnt;    
[r, t(node).mode_ind] = max(t(node).rates ./ t(node).limits);
if r > 1 - tol
    % any broken individual limits
    t(node).rates = t(node).rates./r;
    t(node).mode_cnt = 'INDIVIDUAL';
elseif ~isempty(t(node).children) && ~any_group_controlled_children
    % all children at individual limit
    t(node).mode_cnt = 'NONE';
elseif ~t(node).allow_group_control
    t(node).mode_cnt = 'UNDETERMINED';
else
    t(node).mode_cnt = 'GRUP';
end

% update fractions for groups
if strcmp(t(node).type, 'GROUP') && ~all(t(node).rates == 0)
    t(node).fractions = t(node).rates/sum(t(node).rates);
end

switch_local = ~strcmp(mode_cnt, t(node).mode_cnt);
switched = switched || switch_local;
end

%--------------------------------------------------------------------------
function [t, target_changed] = set_and_update_targets(t, node, tol)
    target_changed = false;
    %cmp = t(node).guideComp;
    % switch to group if group-target is broken
    gt = t(node).group_target;
    if t(node).allow_group_control && t(node).rates(gt.mode) > gt.value*(1 + tol)
        t(node).mode_cnt = 'GRUP';
    end
    % scale rate to meet target
    if strcmp(t(node).mode_cnt, 'GRUP')
        mode_rate = t(node).rates(gt.mode);
        if mode_rate > 0
            r = gt.value/mode_rate;
            t(node).rates = t(node).rates.*r;
        else
            mode_frac = t(node).fractions(gt.mode);
            r = gt.value/mode_frac;
            t(node).rates = t(node).fractions.*r;
        end
    end
    if isempty(t(node).children)
        return
    end
    % distribute to children according to current groups guide-mode
    switch t(node).mode_cnt
        case 'INDIVIDUAL'
            mode = t(node).mode_ind;
            group = node;
        case 'GRUP'
            mode = gt.mode;
            group = gt.group;
        case 'UNDETERMINED'
            mode = t(node).mode_pref;
            group = node;
        case 'NONE'
            mode = gt.mode;
            group = gt.group;
    end
    target = t(node).rates(mode);
    c = t(node).children;
    % get guide-rates corresponding to mode
    guide_rates = arrayfun(@(n)n.guide_rates(mode), t(c));
    % should be guide_rates correponding to mode (not done here)
    %guide_rates = arrayfun(@(n)n.group_target.guiderate, t(c));
    is_grp = cellfun(@(s)strcmp('GRUP', s), {t(c).mode_cnt}).';
    guide_sum = sum(guide_rates(is_grp));
    child_rates_cmp = vertcat(t(c).rates);
    child_rates_cmp = child_rates_cmp(:, mode);
    target_sum = target - sum(child_rates_cmp(~is_grp));
    node_guide_rate_mode = t(node).guide_rates(mode);
    for k = 1:numel(c)
        t(c(k)).group_target.group = group;
        t(c(k)).group_target.mode = mode;
        t(c(k)).group_target.guide_rate = guide_rates(k);
        % ratio used for switching to fall-back, not used here
        t(c(k)).group_target.guide_rate_ratio = guide_rates(k)/node_guide_rate_mode;
        old_target = t(c(k)).group_target.value;
        if is_grp(k)
            new_target = (guide_rates(k)/guide_sum) * target_sum;
        else
            new_target = (guide_rates(k)/(guide_sum + guide_rates(k))) * (target_sum + child_rates_cmp(k));
        end
        t(c(k)).group_target.value = new_target;
        [t, sub_targets_changed] = set_and_update_targets(t, c(k), tol);
        target_changed = target_changed || sub_targets_changed || abs(old_target - new_target) > tol*old_target;
        is_grp_now = strcmp('GRUP', t(c(k)).mode_cnt);
        if is_grp(k) ~= is_grp_now
            if is_grp_now
                % add back to guide_sum/target_sum
                guide_sum = guide_sum + guide_rates(k);
                target_sum = target_sum + child_rates_cmp(k);%t(c(k)).group_target.value;
            else
                % subtract from guide_sum/target_sum
                guide_sum = guide_sum - guide_rates(k);
                target_sum = target_sum - child_rates_cmp(k);%t(c(k)).group_target.value;
            end
        end
    end
end

