load('..\data\rotsync_cornell.mat')
[R_opt, RR_opt] = alter_rot_sync(I, RR, 1e-2); % use 5e-2 for noisy data
[angleDistCurve_abs] = gt_eval_abs(Rgt, R_opt);
[angleDistCurve_relative_opt] = gt_eval_relative(I, Rgt, RR_opt);
[angleDistCurve_relative_init] = gt_eval_relative(I, Rgt, RR);
plot(angleDistCurve_relative_init, 'r');
hold on;
plot(angleDistCurve_relative_opt, 'b');
plot(angleDistCurve_abs, 'k');