function t = balance_group_tree_new(t)
% Fist node t(1) is FIELD-node
% if there are sub-nodes not available for group-controls, corresponding
% sub-trees must be treated first. Get list of sub-tree top-nodes,
% last entry will be global top node (field-node)
topnodes = get_sub_tree_ordering(t);

% loop over all subtrees
for node = topnodes
    % make sure tree is feasible
    t = make_feasible(t, node);
    plotGroupTree(t, 'rate');
    %pause
    t = parametrize_tree(t, node);
    count = 0;
    checkGroupTree(t, node);
    while strcmp(t(node).mode_cnt, 'UNDETERMINED')
        % get valid node with smallest next_alpha
        [next_node, alpha] = getNextLimitNode(t, node);
        % update alpha for all relevant nodes
        t = stepAlpha(t, node, alpha);
        % update node control-status and
        t = updateNode(t, next_node);
        %plotGroupTree(t, 'rate');
        % re-parameterize (this could be done more efficiently)
        t = parametrize_tree(t, node);
        %checkGroupTree(t, node, false);
        count = count + 1;
        %plotGroupTree(t, 'rate');
    end
    % finally update group-target (not fully used/updated during balancing)
    t = set_targets(t);
    checkGroupTree(t, node);
    fprintf('balance_group_tree used %d iterations from node %d\n', count, node)
end
plotGroupTree(t, 'rate');
%pause
end

function t = parametrize_tree(t, node)
    % Assumes a feasible tree!
    % Here we don't need to consider INDIVIDUAL/NONE-nodes since
    % sub-trees of these are already fixed
    if ~t(node).allow_group_control
        t(node).weight = 1;
        mode = t(node).mode_pref;
        t(node).linear_term(mode) = t(node).rates(mode);
    else
        mode = t(node).group_target.mode;
    end
    % update group-rates and set weights for children
    c = t(node).children;
    if ~isempty(c)
        is_grp = cellfun(@(s)strcmp('GRUP', s), {t(c).mode_cnt}).';
        cg = c(is_grp);
        guide_rates = arrayfun(@(n)n.guide_rates(mode), t(cg));
        guide_sum = sum(guide_rates);
        % sum group_rates from children
        linear_term = [0 0 0];
        for k = 1:numel(cg)
            % set linear term
            %weight = t(node).weight * (guide_rates(k)/guide_sum);
            weight =  (guide_rates(k)/guide_sum);
            coeff = weight_transfer_coeff(t, node, cg(k));
            child_mode = t(cg(k)).group_target.mode;
            t(cg(k)).linear_term(child_mode) = t(node).linear_term(mode)*weight*coeff;
            %t(cg(k)).weight = weight*coeff;    
            t = parametrize_tree(t, cg(k));
            linear_term = linear_term + t(cg(k)).linear_term;
        end
        t(node).linear_term = linear_term;  
    else
        % group-controlled well
        t(node).linear_term = (t(node).linear_term(mode)/t(node).rates(mode))*t(node).rates;
    end

    % calculate distance in terms of top-node guided mode
    % avoid division by zero 
    lin_term_tmp = max(eps*t(node).limits, t(node).linear_term);
    [alpha, mode] = min( (t(node).limits - t(node).rates)./lin_term_tmp );
    t(node).next_lim = mode;
    % distance in terms of top-node guide-mode rate
    %guide_mode = t(node).group_target.mode;
    t(node).alpha_lim =alpha;% r*t(node).group_rates(guide_mode)/t(node).weight;
end

function [next_node, alpha] = getNextLimitNode(t, node)
% get node with smallest "distance" to limit from node and downwards
alpha = t(node).alpha_lim;
next_node = node;
c = t(node).children;
if isempty(c)
    return;
end
for k = 1:numel(c)
    % only group controlled children are relevant
    if strcmp(t(c(k)).mode_cnt, 'GRUP')
        [next_node_tmp, dist_tmp] = getNextLimitNode(t, c(k));
        if dist_tmp < alpha
            [next_node, alpha] = deal(next_node_tmp, dist_tmp);
        end
    end
end
end

function  t = stepAlpha(t, node, alpha)
% make a step alpha from node and downwards

%mode_guide = t(node).group_target.mode;
%rate_increase = (t(node).weight*alpha/t(node).group_rates(mode_guide)) * t(node).group_rates;
t(node).rates = t(node).rates + alpha*t(node).linear_term;
t(node).alpha_lim = t(node).alpha_lim - alpha;
%t(node).group_rates = t(node).group_rates + rate_increase;
% calculate distance in terms of top-node guided mode
% avoid division by zero 
%group_rates_tmp = max(eps*t(node).limits, t(node).group_rates);
%[r, mode_lim] = min( (t(node).limits - t(node).rates)./group_rates_tmp );
%t(node).next_lim = mode_lim;
% distance in terms of top-node guide-mode rate
%t(node).dist_to_lim = r*t(node).group_rates(mode_guide)/t(node).weight;
%if t(node).dist_to_lim < 0
%    printf('problem');  
%end
c = t(node).children;
for k = 1:numel(c)
    % only group controlled children are relevant
    if strcmp(t(c(k)).mode_cnt, 'GRUP')
        t = stepAlpha(t, c(k), alpha);
    end
end
end

function t = updateNode(t, node)
% node should have reached it's limit
assert(abs(t(node).alpha_lim) < eps);
t(node).mode_cnt = 'INDIVIDUAL';
t(node).mode_ind = t(node).next_lim;
% if node is not available for group-control, we're done
if ~t(node).allow_group_control
    return
end

% update linear_term for node and all parents
linear_term = t(node).linear_term;
reached_top = false;
cur_node = node;
while ~reached_top
    t(cur_node).linear_term = t(cur_node).linear_term - linear_term;
    t = update_distance_to_lim(t, cur_node);
    if t(cur_node).allow_group_control
        cur_node = t(cur_node).parent;
    else
        reached_top = true;
    end
end

% need to re-parametrize from active parent and down
parent = t(node).parent;
found_active_parent = false;
while ~found_active_parent
    pc = t(parent).children;
    is_grp = cellfun(@(s)strcmp('GRUP', s), {t(pc).mode_cnt}).';
    found_active_parent = any(is_grp);
    if ~found_active_parent
        % parent has no group-controlled children
        % switch to NONE and move up
        t(parent).mode_cnt = 'NONE';
        if t(parent).allow_group_control
            parent = t(parent).parent;
        else
            % we've reached the top -> done
            return
        end
    end
end
% re-parametrize from parent down
t = parametrize_tree(t, parent);
end

function coeff = weight_transfer_coeff(t, parent, child)
% compute transfer coeff in cases guide-mode is different for parent/child
coeff = 1;
from = t(parent).group_target.mode;
to = t(child).group_target.mode;
if strcmp(t(parent).type, 'GROUP') && from ~= to
    coeff = t(parent).linear_term(to)/t(parent).linear_term(from);
end
end

function t = update_distance_to_lim(t, node)
% calculate distance in terms of top-node guided mode
% avoid division by zero
lin_term_tmp = max(eps*t(node).limits, t(node).linear_term);
[alpha, mode] = min( (t(node).limits - t(node).rates)./lin_term_tmp );
t(node).next_lim = mode;
% distance in terms of top-node guide-mode rate
%guide_mode = t(node).group_target.mode;
t(node).alpha_lim =alpha;
end

%--------------------------------------------------------------------------

function t = set_targets(t, node)
if nargin < 2
    node = 1;
end
gt = t(node).group_target;
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
guide_rates = arrayfun(@(n)n.guide_rates(mode), t(c));
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
    if is_grp(k)
        t(c(k)).group_target.value = (guide_rates(k)/guide_sum) * target_sum;
    else
        t(c(k)).group_target.value = (guide_rates(k)/(guide_sum + guide_rates(k))) * (target_sum + child_rates_cmp(k));
    end
    t = set_targets(t, c(k));
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



function topnodes = get_sub_tree_ordering(t, node)
if nargin < 2
    node = 1;
end
topnodes = [];
c = t(node).children;
for k = 1:numel(c)
    topnodes = [get_sub_tree_ordering(t, c(k)), topnodes]; %#ok
end
if ~t(node).allow_group_control
    topnodes = [topnodes, node];
end
if node == 1
    assert(~isempty(topnodes) && topnodes(end)==1)
end
end



