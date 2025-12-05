import STZ
max_alpha = 0.2
max_beta =0.4
total_disp_steps = 349

for size_of_cluster in [1742]:
    for i in [1,2,3]:
        STZ.Normal_Stress_Matrix_p(size_of_cluster,max_alpha,max_beta,total_disp_steps,i,90,9,0)