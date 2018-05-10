% A simple script that calls the plan path function and plots the output Ux
% desired and Ux_dot_desired. To see the planPath function used in a
% simulation, see testPlanPathInSimulation.m


clc; clear; close all;

% Load the project data for the path info
load project_data.mat


% Simulation time
t_final = 8;
dT = 0.01;
t_s = 0:dT:t_final;
N = length(t_s);


[Ux_des, Ux_dot_des] = planPath(path);

figure; subplot(2,1,1); plot(Ux_des); grid on; subplot(2,1,2); plot(path.k_1pm); grid on;
figure; subplot(2,1,1); plot(Ux_dot_des); grid on; subplot(2,1,2); plot(path.k_1pm); grid on;