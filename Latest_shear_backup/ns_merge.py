import numpy as np 
no_task = 8
no_atoms = 1985
total_disp_steps = 159
batch_ns_size = int((total_disp_steps+1)/no_task)
# ns = np.zeros([int(total_disp_steps+1), int(total_disp_steps+1)])
ns = np.zeros([int(total_disp_steps+1), 350])
for iteration in range(1,21):
    for i in range(no_task):
        path = '../Latest_shear_{}/Clusters/{}_cluster/cluster_{}/Normal_stress_{}.txt'.format(i,no_atoms,iteration,i)
        ns_batch = np.loadtxt(path)
        ns[0+batch_ns_size*i:batch_ns_size+batch_ns_size*i,:] = ns_batch
        np.savetxt('../Latest_shear_backup/Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration),ns)











