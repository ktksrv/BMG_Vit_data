import STZ
max_alpha = 0.2
max_beta =0.4
total_disp_steps = 9
no_atoms = 1013
iteration = 1
total_proc = 90
proc_per_task = 9
current_proc = 0
STZ.Normal_Stress_Matrix_bond_len_multi(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc)