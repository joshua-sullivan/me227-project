function F_y = computeCoupledTireForce(tires, Fz, Fx, alpha)

    Ca = tires.C_alpha;
    mu = tires.mu;
    mu_s = tires.mu_s;
    
    % update Fz based on weight dist
    Fz = tires.weightDist * Fz;
    
    zeta = sqrt((mu * Fz)^2 - Fx^2)/(mu * Fz);
    
    alpha_slip = atan2(3 * zeta * mu * Fz , Ca);
    
    % Compute F_y
    if abs(alpha) < alpha_slip
        F_y = -Ca * tan(alpha) + Ca^2 / (3 * zeta * mu * Fz) * abs(tan(alpha)) * tan(alpha) ...
                -(tan(alpha)^3) * Ca^3 / (27 * zeta^2 * mu^2 * Fz^2);
    else        
        F_y = -zeta * mu * Fz * sign(alpha);
    end
    

end