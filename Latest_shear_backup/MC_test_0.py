import STZ
max_alpha = 0.2
max_beta = 0.4
total_disp_steps = 349
lower = -1
upper = 1
steps = 9
total_proc = 32
proc_per_task = 4 
shift = 0.01 
current_proc = 0 

no_atoms = 1499
for iteration in [6,7,8,9,10,11,12,13,14,15]:
    STZ.MC_good_cluster_p_red(no_atoms,max_alpha,max_beta,total_disp_steps,lower,upper,steps,iteration,total_proc,proc_per_task,current_proc,shift)

# no_atoms = 851
# for iteration in [1,2,3,4,5]:
#     STZ.MC_good_cluster_p_red(no_atoms,max_alpha,max_beta,total_disp_steps,lower,upper,steps,iteration,total_proc,proc_per_task,current_proc,shift)

# no_atoms = 1013
# for iteration in [1,2,3,4,5]:
#     STZ.MC_good_cluster_p_red(no_atoms,max_alpha,max_beta,total_disp_steps,lower,upper,steps,iteration,total_proc,proc_per_task,current_proc,shift)


