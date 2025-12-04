clear all; close all; clc;
tic
% 创建电力系统优化对象
powerSystem = PowerSystemOptimization();
% 运行优化调度
powerSystem = powerSystem.run_optimization();
% 可视化结果
powerSystem.visualize_results();
toc;