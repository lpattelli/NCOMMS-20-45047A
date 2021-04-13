function [Lambda, C3, C3_eig, r_C3, t_C3] = coherencymatrix3D(E)
%COHERENCYMATRIX3D Returns quantities related to the coherency matrix of a
% single (3x1) electric field vector value
%
% Input:
% E:        3x1 vector field values (complex)
%
% Output:
% Lambda:   8x1 array with Stokes parameters from 1 to 8 (Lambda_0 = 1)
% C3:       coherency matrix
% C3_eig:   eigenvalues of the coherency matrix, sorted in descending order
% r_C3:     rank(C3)
% t_C3:     rank(real(C3))
%
% Author: Lorenzo Pattelli

E = E/sqrt(conj(E(1))*E(1) + conj(E(2))*E(2) + conj(E(3))*E(3)); % normalize to Lambda_0 = 1;

% as formulated in PRA 90, 043858 (2014)
Lambda = [3 *     (conj(E(1))*E(2) + conj(E(2))*E(1))/2, ...                  % Lambda_1
          3 * 1i* (conj(E(1))*E(2) - conj(E(2))*E(1))/2, ...                  % Lambda_2
          3 *     (conj(E(1))*E(1) - conj(E(2))*E(2))/2, ...                  % Lambda_3
          3 *     (conj(E(1))*E(3) + conj(E(3))*E(1))/2, ...                  % Lambda_4
          3 * 1i* (conj(E(1))*E(3) - conj(E(3))*E(1))/2, ...                  % Lambda_5
          3 *     (conj(E(2))*E(3) + conj(E(3))*E(2))/2, ...                  % Lambda_6
          3 * 1i* (conj(E(2))*E(3) - conj(E(3))*E(2))/2, ...                  % Lambda_7
          sqrt(3)*(conj(E(1))*E(1) + conj(E(2))*E(2) - 2*conj(E(3))*E(3))/2]; % Lambda_8

C3 = [1+Lambda(3)+Lambda(8)/sqrt(3),  Lambda(1) - 1i*Lambda(2),         Lambda(4) - 1i*Lambda(5); ...
      Lambda(1) + 1i*Lambda(2),       1-Lambda(3)+Lambda(8)/sqrt(3),    Lambda(6) - 1i*Lambda(7); ...
      Lambda(4) + 1i*Lambda(5),       Lambda(6) + 1i*Lambda(7),         1 - 2*Lambda(8)/sqrt(3)] / 3;
  
C3_eig = flipud(eig(C3)); % flipped to have lambda1 >= lambda2 >= lambda3

r_C3 = rank(C3);
t_C3 = rank(real(C3));
end

