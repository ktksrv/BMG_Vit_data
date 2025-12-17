import STZ
max_alpha = 0.2
max_beta = 0.4
total_disp_steps = 349
lower = -1
upper = 1
steps = 9
total_proc = 32
proc_per_task = 4
current_proc = 1 

for no_atoms in [270, 446, 666, 851, 1013, 1742, 1985, 1499]:
    for iteration in range(1,21):
        STZ.MC_good_cluster_p_red(no_atoms,max_alpha,max_beta,total_disp_steps,lower,upper,steps,iteration,total_proc,proc_per_task,current_proc)



