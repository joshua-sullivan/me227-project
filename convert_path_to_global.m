function [E, N, psi] = convert_path_to_global( path, s, e, dpsi )
% CONVERT_PATH_TO_GLOBAL converts s and e vectors along a path defined by
% path into EN global coordinates

% Matt Brown, Vincent Laurense, John Subosits

n = length(s);
E = zeros(n,1);
N = zeros(n,1);

centE = interp1(path.s_m,path.posE_m, s);
centN = interp1(path.s_m,path.posN_m, s);
theta = interp1(path.s_m,path.psi_rad, s);
psi = theta + dpsi;

for i=1:n
    E(i) = centE(i) - e(i)*sin(pi/2-theta(i));
    N(i) = centN(i) - e(i)*cos(pi/2-theta(i));
end

end

