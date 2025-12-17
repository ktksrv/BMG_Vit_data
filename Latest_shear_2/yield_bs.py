import STZ
max_alpha=0.2
max_beta=0.4
total_disp_steps = 349
no_atoms = 666
no_segs = 3
for iteration in range(1,21):
    STZ.yield_strain_multi_bs(max_alpha,max_beta,total_disp_steps,no_atoms,iteration,no_segs) 