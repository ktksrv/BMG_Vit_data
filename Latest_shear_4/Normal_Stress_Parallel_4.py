import STZ
max_alpha = 0.2
max_beta =0.4
total_disp_steps = 349

for size_of_cluster in [851]:
    for i in [15,16,17,18,19,20]:
        STZ.Normal_Stress_Matrix_p(size_of_cluster,max_alpha,max_beta,total_disp_steps,i,80,10,4)

for size_of_cluster in [1013,1499]:
    for i in [1,2,3,4,5]:
        STZ.Normal_Stress_Matrix_p(size_of_cluster,max_alpha,max_beta,total_disp_steps,i,80,10,4)

for size_of_cluster in [1742,1985]:
    for i in [1,2,3,4,5,6,7,8,9,10]:
        STZ.Normal_Stress_Matrix_p(size_of_cluster,max_alpha,max_beta,total_disp_steps,i,80,10,4)