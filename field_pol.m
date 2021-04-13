function [e, Ia, Ib, Qs, r_C3, t_C3] = field_pol(E, plot, C3)
%FIELD_POL Plot the time-harmonic oscillation of the electric field vector
%
% Input:
% E:        3x1 vector field values (complex)
% plot:     boolean flag
% C3:       coherency matrix
%
% Output:
% e:        ellipticity
% Ia, Ib:   intensity along the major and minor axes
% Qs:       real(C3) eigenvectors, sorted by descending eigenvalue order
% r_C3:     rank(C3)
% t_C3:     rank(real(C3))
%
% Author: Lorenzo Pattelli

if ~exist('C3','var')
    [~, C3, ~, r_C3, t_C3] = coherencymatrix3D(E);
end

[Q, D] = eig(real(C3));

[~, ind] = sort(diag(D), 'descend');
Ds = D(ind, ind); % sort the eigenvalues -- see: Gil PRA 90, 043856 (2014)
Qs = norm(E)*Q(:, ind); % normalize by the norm of the field for visualization purposes
I = norm(E)^2;
Ia = I*Ds(1,1); Ib = I*Ds(2,2);
e = sqrt(Ib/Ia); % ellipticity

if plot
    t = linspace(0,2*pi,26); t(end) = [];

    figure
    for ti = numel(t):-1:1
        Et(ti,:) = real(E*exp(-1i*t(ti)));
        quiver3(0,0,0,Et(ti,1),Et(ti,2),Et(ti,3),'k'), hold on
        axis equal
    end

    quiver3(0, 0, 0, sqrt(Ds(1,1))*Qs(1,1), sqrt(Ds(1,1))*Qs(2,1), sqrt(Ds(1,1))*Qs(3,1),'r')
    quiver3(0, 0, 0, sqrt(Ds(2,2))*Qs(1,2), sqrt(Ds(2,2))*Qs(2,2), sqrt(Ds(2,2))*Qs(3,2),'g')
    quiver3(0, 0, 0, Qs(1,3), Qs(2,3), Qs(3,3),'b')
    
    xlabel('E_x'), ylabel('E_y'), zlabel('E_z')
end

end

