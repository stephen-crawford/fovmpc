cd ../../lib/mpc_cbf/cmake-build-debug
cmake --build ./ --target mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example
date="02272025"
instances="circle"
#instances="formation"
min_r=2
max_r=4
max_exp=5
min_fov=240
max_fov=360
fov_step=120
min_slack_decay=0.2
max_slack_decay=0.3
slack_decay_step=0.1
default_slack_decay=0.2
for (( exp_idx = 0; exp_idx < ${max_exp}; exp_idx++ )); do
  (for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
    for fov in 120 240 360; do
        echo experiment instance type:${instances} num_robot:${num_r}, exp_idx:${exp_idx}
        ./mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example --instance_type ${instances} --num_robots ${num_r} --fov ${fov} --slack_decay ${default_slack_decay} --write_filename "/home/stephencrawford/fovmpc/experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_decay${default_slack_decay}_States_${exp_idx}.json"
    done
  done) &

  (for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
    for fov in 120 240; do
      for slack_decay in 0.1; do
          echo experiment instance type:${instances} num_robot:${num_r}, exp_idx:${exp_idx}
          ./mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example --instance_type ${instances} --num_robots ${num_r} --fov ${fov} --slack_decay ${slack_decay} --write_filename "/home/stephencrawford/fovmpc/experiments/instances/results/log${date}/${instances}${num_r}_fov${fov}_decay${slack_decay}_States_${exp_idx}.json"
      done
    done
  done) &

done
wait # Wait for all background processes to complete