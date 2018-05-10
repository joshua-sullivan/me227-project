function path = planPath(path)
% This function outputs desired Ux and Ux_dot based on a provided path. The
% format of path must match to that given by the ME 227 course staff. The
% only other configurable specifications are given below, you may need to set
% those a little conservatively as compared to the specifications provided
% in the project prompt. To test this function see testPlanPath.m and
% testPlanPathInSimulation.m
%
% Strategy for the speed profile is as follows: we first calculate the max
% speed allowed in the arcs of the track, based on this speed we calculate
% the max allowable entry and exit speeds of the clothoids. After this we allow 
% the speed in the straight segments to rise as far as possible until braking is
% required to match the next clothoid's entry speed. We use a helper
% function computeZones (see below) to compute the point on the straight
% segments where the transition from acceleration to decelaration happens.

    % max acceleration specifications
    max_ax_pos = 3; % [m/s^2] 
    max_ax_neg = -4; % [m/s^2]
    max_a = 4; % [m/s^2]
    
    % initialize outputs
    Ux_des = zeros(length(path.s_m), 1);
    Ux_dot_des = zeros(length(path.s_m), 1);

    % specify slope factor of the clothoids
    c = sqrt(0.11536 / (2 * 12));

    % specify path coordinates along the clothoids, straights, and arcs
    s_c = [ path.s_m(49) path.s_m(97) ...
            path.s_m(158) path.s_m(206) ...
            path.s_m(301) path.s_m(349) ...
            path.s_m(410) path.s_m(457) ...
            path.s_m(553) path.s_m(600) ...
            path.s_m(661) path.s_m(709) ...
            path.s_m(805) path.s_m(853) ...
            path.s_m(913) path.s_m(961) ...
        ];
    s_s = [ 0 path.s_m(49) ...
            path.s_m(206) path.s_m(301) ...
            path.s_m(457) path.s_m(553) ...
            path.s_m(709) path.s_m(805) ...
            path.s_m(961) path.s_m(1009) ];
    s_a = [ path.s_m(97) path.s_m(158) ...
            path.s_m(349) path.s_m(410) ...
            path.s_m(600) path.s_m(661) ...
            path.s_m(853) path.s_m(913)];

    % determining the speeds for the arcs
    K = 0.1154;
    Ux_arc = sqrt(max_a * (1 / K));
    for i = 1:2:length(s_a)
        % initial and final path coords
        s_0 = s_a(i);
        s_f = s_a(i + 1);

        % find indices for the clothoid
        [~,ind_0] = min(abs(path.s_m - s_0));
        [~,ind_f] = min(abs(path.s_m - s_f));

        % setting constant U desired
        Ux_des(ind_0 : ind_f) = Ux_arc;

        % setting constant Ux dot desired
        Ux_dot_des(ind_0 : ind_f) = 0;

    end

    % determinig the speeds for clothoids that start after a straight segment
    for i = 1:4:length(s_c)

        % initial and final path coords
        s_0 = s_c(i);
        s_f = s_c(i + 1);

        % find indices for the clothoid
        [~,ind_0] = min(abs(path.s_m - s_0));
        [~,ind_f] = min(abs(path.s_m - s_f));

        % set Ux desired at the end
        Ux_des(ind_f) = Ux_arc;

        % integrate backwards
        for j = (ind_f - 1) : -1: ind_0
%             ds = path.s_m(j) - path.s_m(ind_0);
%             dUxds = (1 / Ux_des(j + 1)) * sqrt(max_a^2 - (2 * c^2 * ds * Ux_des(j + 1))^2);
%             Ux_des(j) = Ux_des(j + 1) + dUxds * (path.s_m(j + 1) - path.s_m(j));
            Ux_des(j) = Ux_arc;
            Ux_dot_des(j) = 0;
        end
    end

    % determinig the speeds for clothoids that start after a curve
    for i = 3:4:length(s_c)

        % initial and final path coords
        s_0 = s_c(i);
        s_f = s_c(i + 1);

        % find indices for the clothoid
        [~,ind_0] = min(abs(path.s_m - s_0));
        [~,ind_f] = min(abs(path.s_m - s_f));

        % set Ux desired at the start
        Ux_des(ind_0) = Ux_arc;

        % integrate forwards
        for j = (ind_0 + 1) : 1 : ind_f 
%             ds = path.s_m(j) - path.s_m(ind_0);
%             dUxds = (1 / Ux_des(j - 1)) * sqrt(max_a^2 - (2 * c^2 * ds * Ux_des(j - 1))^2);
%             Ux_des(j) = Ux_des(j - 1) + dUxds * (path.s_m(j) - path.s_m(j - 1));
            Ux_des(j) = Ux_arc;
            Ux_dot_des(j) = 0;
        end
    end

    % setting speeds along straight paths
    for i = 1 : 2 : length(s_s)

        % specify start and end coordinates of the current straight section
        s_0 = s_s(i);
        s_f = s_s(i + 1);

        % find indices for the clothoid
        [~,ind_1] = min(abs(path.s_m - s_0));
        [~,ind_3] = min(abs(path.s_m - s_f));

        % set Ux desired at the start and end of the straight section
        if i == 1
            Ux_des(ind_1) = 0;
            Ux_des(ind_3) = Ux_des(ind_3);
        elseif i == length(s_s) - 1
            % obey the "stop 3m before the end" constraint
            [~, ind_3] =  min(abs(path.s_m - 248));
            Ux_des(ind_1) = Ux_des(ind_1);
            Ux_des(ind_3) = 0;
        else
            Ux_des(ind_1) = Ux_des(ind_1);
            Ux_des(ind_3) = Ux_des(ind_3);
        end

        % compute x2 and ux2
        [ux2, x2] = computeZones(s_0, s_f, Ux_des(ind_1), Ux_des(ind_3), max_ax_pos, max_ax_neg);
        [~, ind_2] =  min(abs(path.s_m - x2));
        Ux_des(ind_2) = ux2;
        
        % if it's the last straight, start braking a little early
        if i == length(s_s) - 1
            ind_2 = ind_2 - 10;
        end

        % integrate forwards in the acceleration zone
        for j = (ind_1 + 1) : 1 : ind_2
            Ux_des(j) = sqrt(Ux_des(j - 1)^2 + 2 * max_ax_pos * (path.s_m(j) - path.s_m(j-1)));
            Ux_dot_des(j) = max_ax_pos;
        end

        % integrate backwards in the braking zone
        for j = (ind_3 - 1) : -1 : ind_2
            Ux_des(j) = sqrt(Ux_des(j + 1)^2 - 2 * max_ax_neg * (path.s_m(j + 1) - path.s_m(j)));
            Ux_dot_des(j) = max_ax_neg;
        end

    end
    
    path.Ux_des_mps = Ux_des;
    path.Ux_dot_des_mps2 = Ux_dot_des;

end




function [ux2, x2] = computeZones(x1, x3, ux1, ux3, max_ax_pos, max_ax_neg)
    syms ux2 x2
    eqns = [ux2 ^ 2 - ux1 ^ 2 - 2 * max_ax_pos * (x2 - x1) == 0, ...
            ux3 ^ 2 - ux2 ^ 2 - 2 * max_ax_neg * (x3 - x2) == 0];
    S = solve(eqns, [ux2 x2]);
    ux2vec = eval(S.ux2);
    x2vec = eval(S.x2);
    
    if ux2vec(1) < 0
        ux2 = ux2vec(2);
    else
        ux2 = ux2vec(1);
    end
    
    if x2vec(1) < 0
        x2 = x2vec(2);
    else
        x2 = x2vec(1);
    end
    
end