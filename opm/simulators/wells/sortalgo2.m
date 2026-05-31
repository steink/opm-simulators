function t = sortalgo2(t)
tol = 1e-8;
t = cap_individual_and_sum(t, 1, tol);
topnodes = get_sub_tree_ordering(t);
for k = 1:numel(topnodes)
    node = topnodes(k);
    if strcmp(t(node).mode_category, 'INDIVIDUAL')
        mode = t(node).mode;
    else
        mode = t(node).mode_preferred;
    end
    qm = t(node).limits(mode);
    t = balance_group_tree(t, node, mode, qm, tol);
    t = set_targets(t, node);
end
if ~checkGroupTree(t)
    sprintf('Tree not balanced!!')
end
end

%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
function t = balance_group_tree(t, node, target_mode, target_rate, tol)
skip_balancing =  false;%t(node).is_balanced && ~new_target;
% no changes - skip balancing, but add to origin
if ~skip_balancing
    max_resorting_count = 10;
    max_top_switch_count = 2;
    mode = target_mode;
    qm = target_rate;
    % scale rates to target or strictest individual control
    rates = t(node).rates;
    scaled_rates = rates * (qm/rates(mode));
    [r, rmode] = max(scaled_rates ./ t(node).limits);
    if r > 1-tol
        mode = rmode;
        qm = scaled_rates(mode)/r;
    end
    t(node).balancing_count = t(node).balancing_count +1;
    any_group_controlled_children = false;
    if strcmp(t(node).type, 'well')
        % just scale rates
        rates = t(node).rates;
        t(node).rates = rates.*(qm/rates(mode));
        t(node).is_balanced = true;
    else % group
        balanced = false;
        resorting_count = 0;
        top_switch_count = 0;    
        while ~balanced && resorting_count <= max_resorting_count && top_switch_count <= max_top_switch_count
            %resorting_count = resorting_count + 1;
            [c, c_fixed, c_trans] = getLocalTreeDescendants(t, node);
            % reset rate_sums/guide_rates_sums from c and up to node
            t = reset_rates_and_guide_rate_sums(t, c, c_fixed, c_trans, node, mode);
            % add inn rates from fixed sub-nodes
            for k = 1:numel(c_fixed)
                t = increment_parent_rate_sums(t, c_fixed(k), node);
            end
            % Compute limit-to-guide ratios. If fractions for c are
            % sufficiently accurate, this will give us an order where children
            % hitting their individual limits/become 'none' appear before any 
            % group-controlled children (which will balance the current
            % tree). If order ends up not correct (some c ends up non-group and 
            % appears after a group-controlled c), we redo the distribution
            % with hopefully better fractions
            ratios = get_ratios_for_sorting(t, c);
            [ratios, ix] = sort(ratios);
            needs_resorting = false;
            any_group_controlled_children = false;
            any_children_need_fallback = false;
            for k = 1:numel(ix)
                ck = c(ix(k));
                if t(ck).use_fallback
                    % handle fallback after main loop
                    any_children_need_fallback = true;
                    continue;
                end
                % ck needs a free path - some transparent group in between
                % may have become active
                if hasFreePath(t, ck, node)
                    % Check if any transparent group reaches its limit
                    % before ck, and if so update 
                    if ~isempty(c_trans) && ~any_group_controlled_children
                        [t, c_trans_update] = update_transparent_groups(t, node, c_trans, ratios(k), qm, mode, tol);
                        if ~isempty(c_trans_update)
                            % loop through newly updated groups to update
                            % needs_resorting/any_group_controlled_children
                            any_transparent = false;
                            any_individual = false;
                            for kt = 1:numel(c_trans_update)
                                ctk = c_trans_update(kt);
                                if hasFreePath(t, ctk, node)
                                    any_transparent = any_transparent || strcmp(t(ctk).mode_category, 'TRANSPARENT');
                                    any_individual = any_individual || strcmp(t(ctk).mode_category, 'INDIVIDUAL');
                                end
                            end
                            % Need resorting if any individual after distribution of group rates
                            needs_resorting = any_group_controlled_children && any_individual;
                            % updated and still transparent means some groups have been given group rates
                            any_group_controlled_children = any_group_controlled_children || any_transparent;
                            % remove updated from list 
                            c_trans = setdiff(c_trans, c_trans_update);
                            % jump to next it if ck has been blocked from
                            % group control
                            if ~hasFreePath(t, ck, node)
                                continue
                            end
                        end
                    end

                    gsum = t(node).guide_rate_sums(mode);
                    gk = t(ck).guide_rate_sums(mode); 
                        
                    qm_remain = qm - t(node).rate_sums(mode);
                    % For correct fractions, any c with (qm_remain/gsum) > ratios 
                    % will hit it's individual limit.
                    qk = (gk/gsum)*qm_remain; % target for ck
                    % get copy of current guide_rate_sum before balancing
                    guide_rate_sums_orig = t(ck).guide_rate_sums;
                    % balance t(ck):
                    t = balance_group_tree(t, ck, mode, qk, tol);
                    any_group_controlled_children = any_group_controlled_children || strcmp(t(ck).mode_category, 'GRUP');
                    if ~strcmp(t(ck).mode_category, 'GRUP') && any_group_controlled_children
                        % non-group-controlled appear after first
                        % group-controlled: wrong order 
                        needs_resorting = true;
                    end
                    % update rate_sums/guide_rate_sums from 
                    % ck -> [any transparent] -> node 
                    if t(ck).has_guide_rate
                        t = decrement_parent_guide_rate_sums(t, ck, node, guide_rate_sums_orig);
                        t = increment_parent_rate_sums(t, ck, node);
                    end
                end
            end
            % update rates for any remaining transparent groups:
            for k = 1:numel(c_trans)
                t(c_trans(k)).rates = t(c_trans(k)).rate_sums;
                % if any limits are broken, we need to resort 
                r = max(t(c_trans(k)).rates ./ t(c_trans(k)).limits);
                if r > 1 + tol
                    needs_resorting = true;
                end
            end
            if needs_resorting && resorting_count < max_resorting_count
                resorting_count = resorting_count + 1;
                %t(node).rates = rate_sums;
            else
                if any_children_need_fallback
                    t = distribute_fallback_rates(t, node, c, mode, tol);
                    % update transparent groups once again
                    for k = 1:numel(c_trans)
                        t(c_trans(k)).rates = t(c_trans(k)).rate_sums;
                    end
                end
                resorting_count = 0;
                % check if any new limit is violated
                limit_violated = false;
                rate_sums = t(node).rate_sums;
                if top_switch_count < max_top_switch_count
                    limits = t(node).limits;
                    limits(target_mode) = min(limits(target_mode), target_rate);
                    [r, mode_new] = max(rate_sums ./ limits);
                    if r > 1 -tol && mode ~= mode_new
                        % broken individual limit or group target
                        % if we're back at original mode, we make this the
                        % last iteration
                        if mode_new == target_mode
                            top_switch_count = max_top_switch_count;
                        else
                            top_switch_count = top_switch_count + 1;
                        end
                        limit_violated = true;
                        mode = mode_new;
                        qm = rate_sums(mode)/r;
                    elseif r < 1-tol && any_group_controlled_children
                        % did not manage to distribute everything!
                        assert(false)
                        limit_violated = false;
                    end
                else
                    fprintf('Node %d reached max switch count\n', node);
                end
                % set balanced
                balanced = ~limit_violated; % fix !!!
                t(node).is_balanced = balanced;
                t(node).rates = rate_sums;
            end
        end
        if resorting_count > max_resorting_count
            fprintf('Node %d reached max resort count', node);
        end
    end

    % update category
    t(node).mode = mode;
    t(node).visited = true;
    if strcmp(t(node).type, 'GROUP') && ~any_group_controlled_children
        t(node).mode_category = 'NONE';
        t(node).mode = nan;
    elseif abs(t(node).rates(mode) - t(node).limits(mode)) <= tol*t(node).limits(mode)
        t(node).mode_category = 'INDIVIDUAL';
    elseif abs(t(node).rates(mode) - qm) <= tol*abs(qm)
        if t(node).has_guide_rate
            t(node).mode_category = 'GRUP';
        else
            t(node).mode_category = 'TRANSPARENT';
            t(node).mode = nan;
        end
    else
        % should not occur
        warning('Problematic balancing of tree for node %d\n', node);
    end

    if ~t(node).is_balanced
       fprintf('Node %d reached maximum resort count (%d). Continuing with unbalanced sub-tree.\n', node, max_resorting_count)
    end
end
end

%--------------------------------------------------------------------------
function [d, df, dt] = getLocalTreeDescendants(t, node)
% d: nodes with guide-rate available for guide below node
% df: nodes not available for guide-control (fixed)
% dt: nodes without guide-rate in between d and node
c = t(node).children;
[d, df, dt] = deal(cell(numel(c), 1));
for k = 1:numel(c)
    if ~t(c(k)).allow_group_control  
        df{k} = c(k);
    elseif t(c(k)).has_guide_rate 
        d{k} = c(k);
    else
        [d{k}, df{k}, dt_tmp] = getLocalTreeDescendants(t, c(k));
        dt{k} = [dt_tmp; c(k)];
    end
end
d = vertcat(d{:});
df = vertcat(df{:});
dt = vertcat(dt{:});
end

%--------------------------------------------------------------------------
function t = reset_rates_and_guide_rate_sums(t, c, cf, ct, origin, mode)
% reset guide-rate-sums corresponding to all c guide-rates for c
% reset rate-sums to zero
nodes = [origin; c; cf; ct];
[t(nodes).rate_sums] = deal([0.0, 0.0, 0.0]);
[t(nodes).guide_rate_sums] = deal([0.0, 0.0, 0.0]);
[t(nodes).visited] = deal(false);
[t(ct).mode_category] = deal('TRANSPARENT');
% take extra care for cases where active mode is not the preferred
mode_pref = t(origin).mode_preferred;
mode_is_pref = mode == mode_pref;
guide_sum_origin = 0.0;
if ~mode_is_pref
    for k = 1:numel(c)
        guide_sum_origin = guide_sum_origin + t(c(k)).guide_rates(mode);
    end
end
for k = 1:numel(c)
    rates = t(c(k)).rates;
    if ~mode_is_pref
        rate_frac_mode = rates(mode)/sum(rates);
        guide_ratio_mode = t(c(k)).guide_rates(mode)/guide_sum_origin;
        if rate_frac_mode < sqrt(eps) || guide_ratio_mode < sqrt(eps)
            % these will be treated after rates have been distributed
            % to all other children
            t(c(k)).use_fallback = true;
            t(c(k)).guide_rate_sums = [nan, nan, nan];
            continue;
        end
    end
    % set guide-rate sums
    t(c(k)).use_fallback = false;
    guide_rate = t(c(k)).guide_rates(mode);
    guide_rate_sums = (guide_rate/rates(mode)) * rates;
    t(c(k)).guide_rate_sums = guide_rate_sums;
    parent = c(k);
    while parent ~= origin
        parent = t(parent).parent;
        t(parent).guide_rate_sums = t(parent).guide_rate_sums + guide_rate_sums;
    end
end
end

%--------------------------------------------------------------------------
function [r, mode_ind] = get_ratios_for_sorting(t, c)
% get limit-to-guiderate ratios for distribution ordering
nc = numel(c);
[r, mode_ind] = deal(nan(nc, 1));
for k = 1:nc
    if strcmp(t(c(k)).mode_category, 'NONE')
        limits = t(c(k)).rates;
    else
        limits = t(c(k)).limits;
    end
    is_zero = all(t(c(k)).guide_rate_sums./t(c(k)).guide_rates < sqrt(eps));
    if ~is_zero
        [r(k), mode_ind(k)] = min( (limits - t(c(k)).rate_sums)./t(c(k)).guide_rate_sums );
    else
        [r(k), mode_ind(k)] = deal(inf, nan);
    end
end
end

%--------------------------------------------------------------------------
function flag = hasFreePath(t, node, node_control)
% check that path from node up to node_control is "free",
% meaning there is no node in between with active limit, i.e.,
% all nodes in between are transparent
    if node == node_control || t(node).parent == node_control
        flag = true;
    else
        parent = t(node).parent;
        isFreeStep = strcmp(t(parent).mode_category, 'TRANSPARENT') && ~t(parent).visited;
        flag = isFreeStep && hasFreePath(t, parent, node_control);
    end
end

%--------------------------------------------------------------------------
function [t, c_trans_update] = update_transparent_groups(t, node, c_trans, next_ratio, qm, mode, tol)
% check if any of c_trans (transparent groups) have limit-to-guide ratio
% less than next_ratio. Transparent groups guide is sum of group-controlled child
% guide-rates
found_less_that_next = true;
c_trans_update = [];
while found_less_that_next && ~isempty(c_trans)
    ratios = get_ratios_for_sorting(t, c_trans);
    [ratio_min, ix] = min(ratios);
    gsum = t(node).guide_rate_sums(mode);
    qm_remain = qm - t(node).rate_sums(mode);
    if ratio_min < next_ratio && ratio_min < qm_remain/gsum  
        ck = c_trans(ix);
        if hasFreePath(t, ck, node) % may have been blocked
            c_trans_update = [c_trans_update; ck]; %#ok
            % current sum of individual children
            rate_sums = t(ck).rate_sums;
            % current remaining guiderate sum
            guide_rate_sums_orig = t(ck).guide_rate_sums;
            gk = t(ck).guide_rate_sums(mode);
            qk = (gk/gsum) * qm_remain + rate_sums(mode); % target for ck
            t = balance_group_tree(t, ck, mode, qk, tol);
            % update parent rate-sums
            t = increment_parent_rate_sums(t, ck, node, t(ck).rates - rate_sums);
            t = decrement_parent_guide_rate_sums(t, ck, node, guide_rate_sums_orig);
        end
    else
        found_less_that_next = false;
    end
end
end

%--------------------------------------------------------------------------
function t = increment_parent_rate_sums(t, node, origin, rates)
% increment rates/guide-rates from node up to origin
if nargin < 4
    rates = t(node).rates;
end
parent = node;
while parent ~= origin
    parent = t(parent).parent;
    t(parent).rate_sums = t(parent).rate_sums + rates;
end
end

%--------------------------------------------------------------------------
function t = decrement_parent_guide_rate_sums(t, node, origin, guide_rate_sums)
% subtract guide-rates from node up to origin
parent = node;
while parent ~= origin
    parent = t(parent).parent;
    t(parent).guide_rate_sums = t(parent).guide_rate_sums - guide_rate_sums;
end
end

%--------------------------------------------------------------------------
function t = distribute_fallback_rates(t, node, c, mode, tol)
% Some children in c cannot be controlled by mode due to tiny fraction
% They will be given rates for node's preferred mode in proportion
% to ratio of node's rate-sums of preferred mode to node's guiderate-sums
% (sums over group controlled wells). This ensures that all group-controlled
% wells have approx same rate/guiderate for preferred mode.
rate_sums = 0;
guiderate_sums = 0;
mode_pref = t(node).mode_preferred;
c_fallback = [];
for k = 1:numel(c)
    if t(c(k)).use_fallback
        c_fallback = [c_fallback; c(k)]; %#ok
    elseif strcmp(t(c(k)).mode_category, 'GRUP')
        rate_sums = rate_sums + t(c(k)).rates(mode_pref);
        guiderate_sums = guiderate_sums + t(c(k)).guide_rates(mode_pref);
    end
end
% cornercase where no wells are group-controlled should not happen, but 
% deal with it just in case
any_group_controlled_wells = rate_sums > 0.0;
ratio = rate_sums/guiderate_sums;
for k = 1:numel(c_fallback)
    ck = c_fallback(k);
    if any_group_controlled_wells
        qk = t(ck).guide_rates(mode_pref)*ratio;        
    else
        qk = t(ck).limits(mode_pref);
    end
    t = balance_group_tree(t, ck, mode_pref, qk, tol);
    % when we add in rates, set rate for mode to exaclty zero so
    % we don't trigger any re-balancing (we accept a small error)
    rates = t(ck).rates;
    rates(mode) = 0.0;
    t = increment_parent_rate_sums(t, ck, node, rates);
end
end

%--------------------------------------------------------------------------
%- Pre-processing ---------------------------------------------------------

function [t, switched] = cap_individual_and_sum(t, node, tol)
% From well and up, cap at individual limits and sum
c = t(node).children;
t(node).is_balanced = false;
t(node).balancing_count = 0;
any_group_controlled_children = false;
switched = false;
if ~isempty(c)
    rates = zeros(size(t(node).rates));
    for k = 1:numel(c)
        [t, switched_child] = cap_individual_and_sum(t, c(k), tol);
        any_group_controlled_children = any_group_controlled_children || strcmp(t(c(k)).mode_category, 'GRUP');
        switched = switched || switched_child;  
        rates = rates + t(c(k)).rates;
    end
    t(node).rates = rates;
    % use these as explcit rates here:
    t(node).explicit_rates = rates;
end

mode_category = t(node).mode_category;
t = updateCategory(t, node, any_group_controlled_children, nan, tol);

% update fractions for groups
if strcmp(t(node).type, 'GROUP') && ~all(t(node).rates == 0)
    t(node).fractions = t(node).rates/sum(t(node).rates);
end

switch_local = ~strcmp(mode_category, t(node).mode_category);
switched = switched || switch_local;
end

%--------------------------------------------------------------------------
function t = updateCategory(t, node, any_group_controlled_children, target_mode, tol)
% catgorize group/well
[r, mode] = max(t(node).rates ./ t(node).limits);
if r > 1 - tol
    if ~isnan(target_mode) && mode ~= target_mode
        %fprintf('Target mode is %d, but most violating is %d\n', target_mode, mode);
        mode = target_mode;
    end
    % any broken individual limits
    %t(node).rates = t(node).rates./r;
    t(node).mode_category = 'INDIVIDUAL';
    t(node).mode = mode;
elseif ~isempty(t(node).children) && ~any_group_controlled_children
    % all children at individual limit
    t(node).mode_category = 'NONE';
    t(node).mode = nan;
elseif ~t(node).allow_group_control
    t(node).mode_category = 'UNDETERMINED';
    t(node).mode = nan;
elseif t(node).has_guide_rate
    t(node).mode_category = 'GRUP';
    t(node).mode = target_mode;
else
    t(node).mode_category = 'TRANSPARENT';
    t(node).mode = nan;
end
end

%--------------------------------------------------------------------------
%- Set group-targets from balanced tree -----------------------------------
function t = set_targets(t, node)
if nargin < 2
    node = 1;
end
if isempty(t(node).children)
    return
end
c = getDescendants(t, node);
if strcmp(t(node).mode_category, 'NONE')
    % no target for descendants of 'NONE'
    for k = 1:numel(c)
        t(c(k)).group_target.group = node;
        t(c(k)).group_target.mode = nan;
        t(c(k)).group_target.value = nan;
        t = set_targets(t, c(k));
    end
else
    % node is 'INDIVIDUAL' or 'GROUP' (function not called for 'TRANSPARENT')
    %{
    gt = t(node).group_target;
    if strcmp(t(node).mode_category, 'INDIVIDUAL')
        mode = t(node).mode;
        mode_pref = t(node).mode_preferred;
    elseif strcmp(t(node).mode_category, 'GRUP')
        % mode from ancestor
        ancestor = t(node).group_target.group;
        mode = t(ancestor).mode;
        mode_pref = t(ancestor).mode_preferred;
    else
        assert(false)
    end
    %}
    mode = t(node).mode;
    mode_pref = t(node).mode_preferred;
    % check for fallback modes (removed from calc)
    mode_is_pref = mode == mode_pref;
    any_fallback = ~mode_is_pref && any(vertcat(t(c).use_fallback));
    guide_rates = arrayfun(@(n)n.guide_rates(mode), t(c));
    rates = arrayfun(@(n)n.rates(mode), t(c));
    if any_fallback
        % remove from calculation
        fallback_nodes = vertcat(t(c).use_fallback);
        guide_rates(fallback_nodes) = 0.0;
        rates(fallback_nodes) = 0.0;
    end
    is_grp = cellfun(@(s)strcmp('GRUP', s), {t(c).mode_category}).';
    guide_sum = sum(guide_rates);
    guide_sum_is_group = sum(guide_rates(is_grp));
    target = t(node).rates(mode);
    target_sum = target - sum(rates(~is_grp));
    for k = 1:numel(c)
        t(c(k)).group_target.group = node;
        t(c(k)).group_target.mode = mode;
        t(c(k)).group_target.guide_rate = guide_rates(k);
        if ~mode_is_pref
            if t(c(k)).use_fallback
                % forces fallback
                t(c(k)).group_target.guide_rate_ratio = 0;
                t(c(k)).group_target.value = nan;
                t = set_targets(t, c(k));
                continue;
            end
            % ratio used for switching to fall-back if fractions become small
            t(c(k)).group_target.guide_rate_ratio = guide_rates(k)/guide_sum;
        end
        if is_grp(k)
            t(c(k)).group_target.value = (guide_rates(k)/guide_sum_is_group) * target_sum;
        else
            t(c(k)).group_target.value = (guide_rates(k)/(guide_sum_is_group + guide_rates(k))) * (target_sum + rates(k));
        end
        t = set_targets(t, c(k));
    end
    % set fallback targets
    if ~mode_is_pref
        guide_rates_pref = arrayfun(@(n)n.guide_rates(mode_pref), t(c));
        rates_pref = arrayfun(@(n)n.rates(mode_pref), t(c));
        proportion = sum(rates_pref(is_grp))/sum(guide_rates_pref(is_grp));
        for k = 1:numel(c)
            t(c(k)).group_target_fallback.group = node;
            t(c(k)).group_target_fallback.mode = mode_pref;
            t(c(k)).group_target.guide_rate = guide_rates_pref(k);
            if t(c(k)).use_fallback
                % use exact value set
                t(c(k)).group_target_fallback.value = t(c(k)).rates(mode_pref);
            else
                t(c(k)).group_target_fallback.value = guide_rates_pref(k)*proportion;
            end
        end
    end
end
end


function c = getDescendants(t, node)
c = num2cell(t(node).children);
for k = 1:numel(c)
    ck = c{k};
    if strcmp(t(ck).mode_category, 'TRANSPARENT')
        c{k} = getDescendants(t, ck);
    end
end
c = vertcat(c{:});
end


