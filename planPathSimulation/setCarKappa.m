function car = setCarKappa(car, g)
% Computes the understeer gradient

    car.k_rad = ((car.Wf/car.C_alphaf) - (car.Wr/car.C_alphar)) * (1/g); % rad/m/s^2
    car.k_deg = ((car.Wf/(car.C_alphaf * (pi/180))) - (car.Wr/(car.C_alphar * (pi/180)))) ; % deg/g

end