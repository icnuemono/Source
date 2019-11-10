function [u_freestream,v_freestream, psi_freestream] = freestream_functions(u_inf,N,Y)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
u_freestream = u_inf * ones(N);
v_freestream = zeros(N);

psi_freestream = u_inf.*Y;
end

