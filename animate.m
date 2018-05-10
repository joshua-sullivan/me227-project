function animate(path, veh, dpsi_rad, s_m, e_m, delta_rad)
% ANIMATE animates vehicle simulation data
%   path is a struct describing the desired path, usually from generate_path.m
%   veh is a struct containing vehicle parameters, namely distance from ref point
%       to front axle (a), to rear axle (b), and tire radius (rW)
%   dpsi_rad is the relative heading of the vehicle wrt the path (dpsi = psi_veh - psi_path)
%   dpsi_rad, s_m, and e_m are the states to visualized, these vectors must be the same length
%
% the visualization can be moved forward or backward in time with the arrow keys or a scrollwheel
% animations can be started and paused by pressing 'p'
% videos can be recorded by pressing 'm' to start and stop

% Matt Brown, Vincent Laurense, John Subosits

    arrow_step = 1;                         % stepsize when arrow keys are pressed
    scroll_step = 30;                       % stepsize when mousewheel is scrolled
    timer_step = 5;                         % stepsize when animation is played

    body_color = [.8 0 .2];                 % color of body
    tire_color = [0 0 0];                   % color of tires

    animation_fps = 50;                     % approximate fps of animation
    animation_filename = 'myMovie.avi';     % movie filename
    delta_mag_factor = 1;                   % show animated steer angle as larger than modeled

    % convert s-e back to E-N
    [posE_m posN_m psi_rad] = convert_path_to_global(path, s_m, e_m, dpsi_rad);

    f = figure('Name', 'Simulation Data', 'WindowScrollWheelFcn', @mousescroll_callback, 'WindowKeyPressFcn', @keypress_callback);
    animate_timer = timer('ExecutionMode', 'fixedRate', 'Period', 1/animation_fps);
    set(animate_timer, 'TimerFcn', @timer_callback)
    recording_video = false;
    video_obj = [];

    ax_space = subplot(1,1,1);
    hold on;

    sim.body                = plot(0,0, 'Color', body_color, 'LineWidth', 2.5);
    sim.tire_f              = plot(0,0, 'Color', tire_color, 'LineWidth', 3);
    sim.tire_r              = plot(0,0, 'Color', tire_color, 'LineWidth', 3);
    sim.ref                 = plot(0,0, 'Color', 'k',        'LineWidth', 2);

    plot(ax_space, path.posE_m, path.posN_m, 'k--');

    curI = 1;
    maxI = length(s_m);

    xlabel('E [m]')
    ylabel('N [m]')
    axis equal;
    grid on;
    box on;

    update_axes();

    function keypress_callback(~, event)
        if strcmp(event.Key, 'rightarrow')
            increment_I(arrow_step);
        end
        if strcmp(event.Key, 'leftarrow')
            decrement_I(arrow_step);
        end
        if strcmp(event.Key, 'p')
            if strcmp(animate_timer.Running, 'on')
                stop(animate_timer);
            else
                start(animate_timer);
            end
        end
        if strcmp(event.Key, 'm')
            if recording_video                          % stop recording
                recording_video = false;
                close(video_obj);
                stop(animate_timer);
            elseif strcmp(animate_timer.Running, 'off') % start recording
                recording_video = true;
                video_obj = VideoWriter(animation_filename);
                video_obj.FrameRate = 30;
                video_obj.Quality = 100;
                open(video_obj);
                start(animate_timer);
            end
        end
        update_axes();
    end

    function mousescroll_callback(~,event)
        if event.VerticalScrollCount > 0
            decrement_I(scroll_step);
        else
            increment_I(scroll_step);
        end
        update_axes();
    end

    function timer_callback(~,~)
        if curI == maxI
            stop(animate_timer);
            if recording_video
                recording_video = false;
                close(video_obj);
                return;
            end
        end
        increment_I(timer_step);
        update_axes()
        if recording_video
            cur_frame = getframe(f);
            try % side effect of global recording_video is extra call to timer_callback before recording_video is updated
                writeVideo(video_obj, cur_frame);
            end
        end
    end

    function increment_I(step)
        curI = curI + step;
        if curI > maxI
            curI = maxI;
        end
    end

    function decrement_I(step)
        curI = curI - step;
        if curI < 1
            curI = 1;
        end
    end


    function update_axes()
        update_vehicle(psi_rad(curI), posE_m(curI), posN_m(curI), delta_mag_factor*delta_rad(curI));
        xlim( [ posE_m(curI)-5 posE_m(curI)+5 ]);
        ylim( [ posN_m(curI)-4 posN_m(curI)+4 ]);
        set(f, 'Name', sprintf('Simulation Data, (%d/%d)', curI, maxI));
        drawnow;
    end

    function update_vehicle(psi, posE, posN, delta)
        a = veh.a;
        b = veh.b;
        rW = veh.rW;

        % body
        body_front = [posE - a*sin(psi); posN + a*cos(psi)];
        body_rear = [posE + b*sin(psi); posN - b*cos(psi)];
        body = [body_front body_rear];

        % tires
        tire_f = [body_front(1)-rW*sin(psi+delta) body_front(1)+rW*sin(psi+delta);
                   body_front(2)+rW*cos(psi+delta) body_front(2)-rW*cos(psi+delta)];

        tire_r = [body_rear(1)-rW*sin(psi) body_rear(1)+rW*sin(psi);
                  body_rear(2)+rW*cos(psi) body_rear(2)-rW*cos(psi)];

        % ref point
        ref_rad = .1;
        ang = linspace(0,2*pi,20);
        set(sim.ref, 'Xdata', posE+ref_rad*cos(ang), 'YData', posN+ref_rad*sin(ang));
        set(sim.body, 'XData', body(1,:), 'YData', body(2,:));
        set(sim.tire_f, 'XData', tire_f(1,:), 'YData', tire_f(2,:));
        set(sim.tire_r, 'XData', tire_r(1,:), 'YData', tire_r(2,:));

    end

end
