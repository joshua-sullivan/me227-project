function F_y = computeTireForce(tires, Fz, alpha)

    Ca = tires.C_alpha;
    mu = tires.mu;
    mu_s = tires.mu_s;
    
    % update Fz based on weight dist
    Fz = tires.weightDist * Fz;
    
    alpha_slip = atan2(3 * mu * Fz , Ca);
    

    % Compute F_y
    if abs(alpha) < alpha_slip
        F_y = (-Ca*tan(alpha)) + (((Ca^2)/(3*mu*Fz)*(2-(mu_s/mu)))*abs(tan(alpha))*tan(alpha)) - ...
             (((Ca^3)/(9*mu*mu*Fz*Fz))*(tan(alpha)^3)*(1-(2*mu_s/3/mu)));
    else        
        F_y = -mu_s*Fz*sign(alpha);
    end
    

end