#Function to extract coordinates from .dat file
from numpy import transpose
from sympy import convex_hull


def extract_dat(path, no_atoms):                                              
    import pandas as pd
    import numpy as np

    # Read the first line to check OVITO version
    with open(path, 'r') as file:
        first_line = file.readline()
    
    # Determine number of rows to skip based on version
    if "OVITO Pro 3.7.7" in first_line:
        skip = 8
    elif "OVITO Basic 3.7.5" in first_line:
        skip = 8
    elif "OVITO Pro 3.12.4" in first_line:
        skip = 10
    else:
        raise ValueError("Unknown OVITO version in .dat file")

    # Read atom data
    df = pd.read_csv(path, sep=" ", skiprows=skip, nrows=no_atoms, header=None)
    df.columns = ['ID', 'type', 'X', 'Y', 'Z']
    
    # Extract columns
    ID = np.linspace(1, no_atoms, no_atoms)
    x = np.array(df.X)
    y = np.array(df.Y)
    z = np.array(df.Z)
    
    # Combine and return
    out = np.column_stack([ID, x, y, z])
    return out

def extract_dat_type(path,no_atoms):
                                                              
    import pandas as pd
    import numpy as np

    # Read the first line to check OVITO version
    with open(path, 'r') as file:
        first_line = file.readline()
    
    # Determine number of rows to skip based on version
    if "OVITO Pro 3.7.7" in first_line:
        skip = 8
    elif "OVITO Basic 3.7.5" in first_line:
        skip = 8
    elif "OVITO Pro 3.12.4" in first_line:
        skip = 10
    else:
        raise ValueError("Unknown OVITO version in .dat file")

    # Read atom data
    df = pd.read_csv(path, sep=" ", skiprows=skip, nrows=no_atoms, header=None)
    df.columns = ['ID', 'type', 'X', 'Y', 'Z']
    
    # Extract columns
    ID = np.linspace(1, no_atoms, no_atoms)
    x = np.array(df.X)
    y = np.array(df.Y)
    z = np.array(df.Z)
    typ = np.array(df.type)
    out=np.column_stack([ID,x,y,z])
    return out,typ

def Simulation_box_dim(path):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=3, nrows=3, header=None)
    df.columns=['lo', 'hi', 'label_1','label_2',]
    axis_lo=np.array(df.lo)
    axis_hi=np.array(df.hi)
    return axis_hi, axis_lo
    
def extract_dump(path,no_atoms):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=9, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z']
    #df.drop(df.columns[[-1,-2]], axis=1, inplace=True)
    ID=np.linspace(1,no_atoms,no_atoms)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    out=np.column_stack([ID,x,y,z])
    return out

def extract_dump_type(path,no_atoms):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=9, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z']
    #df.drop(df.columns[[-1,-2]], axis=1, inplace=True)
    ID=np.linspace(1,no_atoms,no_atoms)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    typ = np.array(df.type)
    out=np.column_stack([ID,x,y,z])
    return out,typ


#Function to extract ID of surface atoms
def surface_atoms(initial,no_atoms):
    from ovito.data import Particles, DataCollection, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.modifiers import ConstructSurfaceModifier
    import numpy as np
    particles = Particles()
    data = DataCollection()
    data.objects.append(particles)
    cell = SimulationCell(pbc = (False, False, False))
    cell[...] = [[200,0,0,0],                               #use cell vectors for given congifuration from ovito
                [0,72.56,0,0],
                [0,0,145.12,0]]
    cell.vis.line_width = 0.1
    data.objects.append(cell)
    pos_prop = particles.create_property('Position', data=initial)
    id=range(1,no_atoms+1)
    pipeline = Pipeline(source = StaticSource(data = data))
    pipeline.add_to_scene()
    pipeline.modifiers.append(ConstructSurfaceModifier(method = ConstructSurfaceModifier.Method.AlphaShape, radius=2.55, select_surface_particles=True))
    data= pipeline.compute()
    selection=np.array(data.particles["Selection"])
    surface_ID_all=selection*id
    surface_ID=surface_ID_all[surface_ID_all !=0]
    
    return  surface_ID

def init_surface_atom_str(no_atoms,iteration):
    import STZ
    initial = STZ.extract_dat('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    d_initial = initial[:,1:]
    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "    

    surface_group='\ngroup surface id {}'.format(surface_str)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "w")
    text_file.write(surface_group)
    text_file.close()

#Function to classify upper and lower groups
def slice_atoms(initial):
    import numpy as np
    z_mean=initial[:,3].mean()
    upper=np.zeros(len(initial[:,3]))
    lower=np.zeros(len(initial[:,3]))

    upper_i=np.zeros([len(initial[:,2]),3])
    lower_i=np.zeros([len(initial[:,2]),3])
    for i in range(len(initial[:,3])):
        if initial[i,3]>z_mean:
            upper[i]=i+1
            upper_i[i,:]=initial[i,1:]
        else:
            lower[i]=i+1
            lower_i[i,:]=initial[i,1:]
    upper_id=upper[upper!=0]
    lower_id=lower[lower!=0]
    upper_i=upper_i[~np.all(upper_i==0,axis=1)]
    lower_i=lower_i[~np.all(lower_i==0,axis=1)]

    return upper_id,lower_id,upper_i,lower_i

def extract_fdump(path,no_atoms):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=9, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z', 'F_x', 'F_y', 'F_z']
    #df.drop(df.columns[[-1,-2]], axis=1, inplace=True)
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    fx=np.array(df.F_x)
    fy=np.array(df.F_y)
    fz=np.array(df.F_z)

    out_c=np.column_stack([x,y,z])                                 
    out_f=np.column_stack([fx,fy,fz]) 

    return out_c, out_f

def stress_calc(coordinates,force):
    import numpy as np
    from scipy.spatial import ConvexHull
    sigma=np.zeros([3,3])
    for i in range(len(force)):
        F = force[i,:]
        F=F.reshape(3,1)
        R = coordinates[i,:]
        sigma = sigma + np.kron(F,R)
    hull= ConvexHull(coordinates) 
    volume=hull.volume
    sigma=-sigma/volume
    sigma= sigma*1.6021766208e+11   #eV/ang^3  to SI
    return sigma[2,2], sigma[0,2]

def stress_config(total_disp_x,total_disp_z,total_disp_steps,no_atoms):
    from lammps import lammps
    import STZ
    initial=STZ.extract_dat('L_100_STZ.dat',no_atoms)
    upper_id,_,_,_,=STZ.slice_atoms(initial)

    for j in  range(total_disp_steps):
        disp_x=total_disp_x/total_disp_steps
        disp_z=total_disp_z/total_disp_steps
        for m in upper_id:
            initial[int(m)-1,3]+=disp_z
            initial[int(m)-1,1]+=disp_x
        lmp=lammps()
        lmp_ovito=lammps()
        lmp_ovito.file('1_3.in')
        create_atoms=['create_atoms 1 single {} {} {}'.format(initial[l,1],initial[l,2],initial[l,3]) for l in range(len(initial[:,0]))]
        lmp_ovito.commands_list(create_atoms)
        lmp_ovito.command('dump dump_1 all custom 1 dump.ovito id type x y z')
        lmp_ovito.command('run 1')
        surface_id=STZ.surface_atoms('dump.ovito',no_atoms)
        lmp.file('1_3.in')
        lmp.commands_list(create_atoms)
        surface_id_str=surface_id.astype('str')
        surface_str=" "
        for m in range(len(surface_id)):
            surface_str=surface_str+surface_id_str[m]
            surface_str+=" "
        surface_group=['group surface id {}'.format(surface_str)]
        lmp.commands_list(surface_group)
        minimization_block='''
        group not_surface subtract all surface
        fix 1 surface setforce 0 0 0 
        minimize 1e-40 1e-40 10000 10000
        unfix 1
        '''
        lmp.commands_string(minimization_block)
        lmp.command('dump dump_1 all custom 1 dump.minimized id type x y z')
        lmp.command('compute force all property/atom fx fy fz')
        lmp.command('dump fcal all custom 1 dump.force id type x y z fx fy fz')    
        lmp.command('run 1')

        initial=STZ.extract_dump('dump.minimized',no_atoms)
        upper_id,_,_,_,=STZ.slice_atoms(initial)
    #lmp.command('minimize 1e-40 1e-40 10000 10000')
    coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
    ns,ss=STZ.stress_calc(coordinates,force)
    return ns, ss

def stress_config_grad(alpha,beta,initial,type_particle,no_atoms,iteration):
    import numpy as np
    import STZ
    from lammps import lammps
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()     
    initialization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001
    region myregion block {} {} {} {} {} {}  units box
    create_box 5 myregion
    mass 1 91.224
    mass 2 63.546
    mass 3 58.693
    mass 4 26.982
    mass 5 47.867 
    pair_style eam/alloy
    pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti
    '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

    create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
    create_atoms_str = '\n'.join(str(e) for e in create_atoms)

    minimization_block='''
    fix freeze surface setforce 0 0 0 
    minimize 0 1e-4 100000 100000
    unfix freeze
    compute force all property/atom fx fy fz
    dump fcal all custom 1 dump.force id type x y z fx fy fz
    run 1 
    '''
    lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
    lmp = lammps()
    lmp.commands_string(lammps_input_script)
    coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
    ns,ss=STZ.stress_calc(coordinates,force)
    return ns, ss


def random_cluster_centroid(x_ul,x_ll,y_ul,y_ll,z_ul,z_ll):
    import numpy as np
    x= x_ll+(x_ul-x_ll)*np.random.random()
    y= y_ll+(y_ul-x_ll)*np.random.random()
    z= z_ll+(z_ul-z_ll)*np.random.random()
    return x, y, z


def random_cluster_generator(size,iteration):

    from ovito.data import Particles, DataCollection, SimulationCell,ParticleType
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.modifiers import ExpressionSelectionModifier, DeleteSelectedModifier
    from lammps import lammps
    from ovito.io import export_file
    import numpy as np
    import STZ

    # Random Cluster of given size
    current_count_diff=10 #flag
    #initial BMG size and parameters which define the usable region for cluster extraction
    no_atoms=32800                                  
    #radius for midpoint distance of the ranges below
    # r1,r2,r3 = 80,170,240   #Cu-Cu    
    # r1,r2,r3 = 85,180,250   #Pd-Ni-P 
    r1,r2,r3 = 105,220,310   #vit-105

    if size>0 and size<500 :
        r_initial=r1
    elif size>500 and size<1000 :
        r_initial=r2
    elif size>1000 and size<1500 :
        r_initial=r3
    r_updated=r_initial
    # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(54,20,54,20,90,54) #Cu-Cu
    # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(45,30,45,30,100,40) #Pd-Ni-P
    x_rand,y_rand,z_rand= STZ.random_cluster_centroid(69,21,69,21,106,74) #vit-105
    initial,type_particle=STZ.extract_dat_type('Initial_BMG.dat',no_atoms)
    iterations=0
    particles = Particles()
    data = DataCollection()
    data.objects.append(particles)
    cell = SimulationCell(pbc = (False, False, False))
    #use cell vectors for given congifuration from ovito
    #Cu-Cu
    # cell[...] = [[72.56,0,0,0],                               
    #             [0,72.56,0,0],
    #             [0,0,145.12,0]]
    # #Pd-Ni-P
    # cell[...] = [[60.2781,0,0,6.14095],                               
    #             [0,60.2781,0,6.14095],
    #             [0,0,120.556,12.2819]]
    #vit
    cell[...] = [[90.51,0,0,0],                               
                [0,90.51,0,0],
                [0,0,181.02,0]]
    cell.vis.line_width = 0.1
    data.objects.append(cell)
    pos_prop = particles.create_property('Position', data=initial[:,1:])
    type_prop = particles.create_property('Particle Type')
    type_prop.types.append(ParticleType(id = 1, name = 'Zr', color = (0.0,1.0,0.0)))
    type_prop.types.append(ParticleType(id = 2, name = 'Cu', color = (1.0,0.0,0.0)))
    type_prop.types.append(ParticleType(id = 3, name = 'Ni', color = (0.0,0.0,1.0)))
    type_prop.types.append(ParticleType(id = 4, name = 'Al', color = (0.0,0.5,0.5)))
    type_prop.types.append(ParticleType(id = 5, name = 'Ti', color = (0.5,0.5,0.5)))

    for ind in range(no_atoms):
        type_prop[ind] = type_particle[ind]
    pipeline = Pipeline(source = StaticSource(data = data))
    pipeline.add_to_scene()
    while (current_count_diff!=0):
        if iterations==1000:
            # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(54,20,54,20,90,54) #Cu-Cu
            # x_rand,y_rand,z_rand= STZ.random_cluster_centroid(45,30,45,30,100,40) #Pd-Ni-P
            x_rand,y_rand,z_rand= STZ.random_cluster_centroid(69,21,69,21,106,74) #vit-105
            iterations=0
        pipeline.modifiers.append(ExpressionSelectionModifier(expression = '(Position.X-{})^2+(Position.Y-{})^2+(Position.Z-{})^2 <{}'.format(x_rand,y_rand,z_rand,r_updated)))
        data= pipeline.compute()
        current_count=data.attributes['ExpressionSelection.count']
        current_count_diff=size-current_count
        del pipeline.modifiers[0]
        iterations+=1
        r_updated= r_updated+0.01*current_count_diff

    pipeline.modifiers.append(ExpressionSelectionModifier(expression = '(Position.X-{})^2+(Position.Y-{})^2+(Position.Z-{})^2 > {}'.format(x_rand,y_rand,z_rand,r_updated)))
    data = pipeline.compute()
    count=data.attributes['ExpressionSelection.count']
    pipeline.modifiers.append(DeleteSelectedModifier())
    export_file(pipeline, "Clusters/{}_cluster/cluster_{}/cluster_{}.dat".format(size,iteration,size), "lammps/data")
    lmp=lammps()
    # Minimize the cluster

    # Cu-Cu
    # minimization_block='''
    # dimension 3
    # units metal
    # boundary s s s
    # atom_style atomic
    # timestep 0.001

    # read_data Clusters/{}_cluster/cluster_{}/cluster_{}.dat
    # mass 1 63.546
    # pair_style eam
    # pair_coeff 1 1 Cu_u6.eam
    # minimize 0 1e-4 100000 100000


    # dump dump_1 all custom 1 Clusters/{}_cluster/cluster_{}/cluster_{}.dat id type x y z
    # run 1

    # '''.format(size,iteration,size,size,iteration,size)
    # lmp.commands_string(minimization_block)


    # #Pd-Ni-P
    # minimization_block='''
    # dimension 3
    # units metal
    # boundary s s s
    # atom_style atomic
    # timestep 0.001

    # read_data Clusters/{}_cluster/cluster_{}/cluster_{}.dat
    # mass 1 106.42 
    # mass 2 63.546
    # mass 3 58.693
    # mass 4 30.974
    # pair_style nep nep.population160.generation1013200.txt
    # pair_coeff * * 

    # minimize 0 1e-4 100000 100000


    # dump dump_1 all custom 1 Clusters/{}_cluster/cluster_{}/cluster_{}.dat id type x y z

    # run 1

    # '''.format(size,iteration,size,size,iteration,size)
    # lmp.commands_string(minimization_block)


    #vit-105
    minimization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001

    read_data Clusters/{}_cluster/cluster_{}/cluster_{}.dat
    mass 1 91.224
    mass 2 63.546
    mass 3 58.693
    mass 4 26.982
    mass 5 47.867 
    pair_style eam/alloy
    pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti 

    minimize 0 1e-4 100000 100000


    dump dump_1 all custom 1 Clusters/{}_cluster/cluster_{}/cluster_{}.dat id type x y z

    run 1

    '''.format(size,iteration,size,size,iteration,size)
    lmp.commands_string(minimization_block)


    # Affine Transformations
    initial,type_particle=STZ.extract_dump_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(size,iteration,size),size)
    particles = Particles()
    data = DataCollection()
    data.objects.append(particles)
    cell = SimulationCell(pbc = (False, False, False))
    #Transformed Simulation box, which should contain all the sheared atoms
    # Cu-Cu
    # cell[...] = [[200,0,0,0],                               
    #             [0,72.56,0,0],
    #             [0,0,145.12,0]]
    #Pd-Ni-P
    # cell[...] = [[180,0,0,6.14095],                               
    #             [0,60.2781,0,6.14095],
    #             [0,0,120.556,12.2819]]
    #vit-105
    cell[...] = [[250,0,0,0],                               
                [0,90.51,0,0],
                [0,0,181.02,0]]

    cell.vis.line_width = 0.1
    data.objects.append(cell)
    pos_prop = particles.create_property('Position', data=initial[:,1:])
    pos_prop = particles.create_property('Position', data=initial[:,1:])
    type_prop = particles.create_property('Particle Type')
    type_prop.types.append(ParticleType(id = 1, name = 'Zr', color = (0.0,1.0,0.0)))
    type_prop.types.append(ParticleType(id = 2, name = 'Cu', color = (1.0,0.0,0.0)))
    type_prop.types.append(ParticleType(id = 3, name = 'Ni', color = (0.0,0.0,1.0)))
    type_prop.types.append(ParticleType(id = 4, name = 'Al', color = (0.0,0.5,0.5)))
    type_prop.types.append(ParticleType(id = 4, name = 'Ti', color = (0.0,0.5,0.5)))

    for ind in range(size):
        type_prop[ind] = type_particle[ind]
    pipeline = Pipeline(source = StaticSource(data = data))
    pipeline.add_to_scene()
    export_file(pipeline, "Clusters/{}_cluster/cluster_{}/cluster_{}.dat".format(size,iteration,size),"lammps/data")



def Normal_Stress_Matrix(no_atoms,max_alpha,max_beta,total_disp_steps,iteration):
    from lammps import lammps
    import numpy as np
    import STZ
    ## NORMAL STRESS MATRIX GENERATION ##
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    ns=np.zeros([len(d_alpha),len(d_beta)])
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    initial,type_particle=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    y=-1
    for i in d_alpha:
        alpha=i
        y=y+1
        for j in range(0,total_disp_steps+1):
            beta=(max_beta/total_disp_steps)*j
            deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
            d_initial=np.zeros([no_atoms,3])
            for k in range(no_atoms):
                d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
            surface_id=STZ.surface_atoms(d_initial,no_atoms)
            surface_id_str=surface_id.astype('str')
            surface_str=" "
            for r in range(len(surface_id)):
                surface_str=surface_str+surface_id_str[r]
                surface_str+=" "     
            lmp=lammps()
            initialization_block='''
            dimension 3
            units metal
            boundary s s s
            atom_style atomic
            timestep 0.001
            region myregion block 6.14095 186.14095 6.14095 66.41905 12.2819 132.8379 units box
            create_box  4 myregion
            mass 1 106.42 
            mass 2 63.546
            mass 3 58.693
            mass 4 30.974
            pair_style nep nep.population160.generation1013200.txt
            pair_coeff * * 
            '''
            lmp.commands_string(initialization_block)
            create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            lmp.commands_list(create_atoms)
            surface_group=['group surface id {}'.format(surface_str)]
            lmp.commands_list(surface_group)
            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            '''
            lmp.commands_string(minimization_block)
            lmp.command('compute force all property/atom fx fy fz')
            lmp.command('dump fcal all custom 1 dump.force id type x y z fx fy fz')    
            lmp.command('run 1')
            coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
            ns[y,j],_=STZ.stress_calc(coordinates,force)
        np.savetxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration),ns)

def MC_test(tolerance,max_alpha,max_beta,total_disp_steps,no_atoms,iteration,lower,upper,steps):

    import numpy as np
    import matplotlib.pyplot as plt
    import STZ
    from scipy import interpolate, stats 

    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    initial=STZ.extract_dat('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    
    #tau_o 
    normal_s = 0
    cs = plt.contour(Beta,Alpha,ns,normal_s)
    dat0=cs.allsegs[1][0]
    ss_lvl=np.zeros(len(dat0))
    s_strain=np.zeros(len(dat0))
    

    # for j in range(len(dat0)):
    #     _,ss_lvl[j]=STZ.stress_config_grad(dat0[j][1],dat0[j][0],initial,no_atoms)
    #     s_strain[j]=dat0[j][0]
    
    # x_vals = np.linspace(0,max_beta,100)
    # splines=interpolate.splrep(s_strain,ss_lvl)
    # y_vals_stress=interpolate.splev(x_vals, splines)
    # y_vals=interpolate.splev(x_vals, splines, der=2)
    # index = np.where(y_vals < tolerance)
    # tau_o=y_vals_stress[index[0][0]]
    
    tau_o = 9.931594491454593658e+09
    #Multiples of tau_o
    factors=np.linspace(lower,upper,steps)
    normal_s = factors*tau_o
    max_shear_lvl = np.zeros(len(normal_s))
    max_yeild_strain = np.zeros(len(normal_s))

    string='Normal_stress_'

    cs = plt.contour(Beta,Alpha,ns,normal_s)
    
    
    for  i in range(0,len(normal_s)):
        dat0=cs.allsegs[i][0]
        ss_lvl=np.zeros(len(dat0))
        s_strain=np.zeros(len(dat0)) 
        ns=np.zeros([len(normal_s),len(dat0)])
        for j in range(len(dat0)):
            ns[i,j],ss_lvl[j]=STZ.stress_config_grad(dat0[j][1],dat0[j][0],initial,no_atoms)
            s_strain[j]=dat0[j][0]
        # x_vals = np.linspace(0,max_beta,100)
        # splines=interpolate.splrep(s_strain,ss_lvl)
        # y_vals_stress=interpolate.splev(x_vals, splines)
        # y_vals=interpolate.splev(x_vals, splines, der=2)
        # index = np.where(y_vals < tolerance)
        # max_shear_lvl[i]=y_vals_stress[index[0][0]]
        # max_yeild_strain[i]=x_vals[index[0][0]]
        shear_lvl_data=np.column_stack((ss_lvl,s_strain))
        string=string+str(round(normal_s[i]/tau_o,3))   
        np.savetxt('Clusters/{}_cluster/cluster_{}/{}.txt'.format(no_atoms,iteration,string),shear_lvl_data)
        # np.savetxt('Clusters/{}_cluster/cluster_{}/ns_{}.txt'.format(no_atoms,iteration,string),ns)
        # plt.plot(x_vals,y_vals_stress)
        # plt.plot(max_yeild_strain[i],max_shear_lvl[i], marker='x', markersize=10, color='r' )
        # plt.title(' Normalised Normal stress={}'.format(round(factors[i],3)))
        # plt.xlabel('Shear strain')
        # plt.ylabel('Shear stress (N/m\u00b2)')
        # plt.savefig('Clusters/{}_cluster/cluster_{}/{}.png'.format(no_atoms,iteration,string), facecolor='white', transparent=False)
        # plt.clf()
        string='Normal_stress_'
    # plt.figure(figsize=[6,4], dpi=300)
    # plt.scatter(factors, max_shear_lvl/tau_o, marker='.',s=150)
    # plt.title(' MC Test (Linear fit)')
    # plt.xlabel("Normalised normal stress")
    # plt.ylabel("Normalised shear stress")
    # slope, intercept, r, _, _ = stats.linregress(factors, max_shear_lvl/tau_o)
    # x=np.linspace(lower,upper,10)
    # y=slope*x+intercept
    # plt.plot(x,y,'r')
    # plt.text(lower,1,'y={}x+{}'.format(round(slope,5),round(intercept,5)))
    # plt.text(lower,0.9,'R\u00b2={}'.format(r**2))
    # plt.savefig('Clusters/{}_cluster/cluster_{}/{}_cluster_MC.png'.format(no_atoms,iteration,no_atoms), facecolor='white', transparent=False)
    # MC_data=np.column_stack((normal_s,max_shear_lvl,max_yeild_strain))
    plt.clf()
    # np.savetxt('Clusters/{}_cluster/cluster_{}/MC_data'.format(no_atoms,iteration),MC_data)
    


def convex_hull_relaxed_config(alpha,beta,initial,no_atoms):
    import numpy as np
    import STZ
    from lammps import lammps
    from scipy.spatial import ConvexHull
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "     
    lmp=lammps()
    initialization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001
    region myregion block 0.0 200.0 0.0 72.56 0.0 145.12  units box
    create_box  1 myregion
    mass 1 63.546

    pair_style eam
    pair_coeff 1 1 Cu_u6.eam '''
    lmp.commands_string(initialization_block)
    create_atoms=['create_atoms 1 single {} {} {}'.format(d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
    lmp.commands_list(create_atoms)
    surface_group=['group surface id {}'.format(surface_str)]
    lmp.commands_list(surface_group)
    
    minimization_block='''
    fix freeze surface setforce 0 0 0 
    minimize 0 1e-4 100000 100000
    unfix freeze
    '''
    lmp.commands_string(minimization_block)
    lmp.command('compute force all property/atom fx fy fz')
    lmp.command('dump fcal all custom 1 dump.force id type x y z fx fy fz')    
    lmp.command('run 1')
    coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
    hull= ConvexHull(coordinates) 
    volume=hull.volume

    return volume


def Convex_hull__volume_zero_ns(max_alpha,max_beta,total_disp_steps,no_atoms,iteration):

    import numpy as np
    import matplotlib.pyplot as plt
    import STZ
    from scipy import interpolate, stats 

    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    initial=STZ.extract_dat('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    
    #tau_o 
    normal_s = 0
    cs = plt.contour(Beta,Alpha,ns,normal_s)
    dat0=cs.allsegs[1][0]
    convex_hull_vol=np.zeros(len(dat0))
    s_strain=np.zeros(len(dat0)) 

    for j in range(len(dat0)):
        convex_hull_vol[j]=STZ.convex_hull_relaxed_config(dat0[j][1],dat0[j][0],initial,no_atoms)
        s_strain[j]=dat0[j][0]
    
    return convex_hull_vol, s_strain


def RDF_plots_normal_stress(initial,no_atoms,beta,iteration,i):
    from ovito.data import Particles, DataCollection, SimulationCell
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.modifiers import CoordinationAnalysisModifier
    import matplotlib.pyplot as plt    
    particles = Particles()
    data = DataCollection()
    data.objects.append(particles)
    cell = SimulationCell(pbc = (False, False, False))
    cell[...] = [[200,0,0,0],                               #use cell vectors for given congifuration from ovito
                [0,72.56,0,0],
                [0,0,145.12,0]]
    cell.vis.line_width = 0.1
    data.objects.append(cell)
    pos_prop = particles.create_property('Position', data=initial)
    id=range(1,no_atoms+1)
    pipeline = Pipeline(source = StaticSource(data = data))
    pipeline.add_to_scene()
    modifier = CoordinationAnalysisModifier(cutoff = 10.0 , number_of_bins=100)
    pipeline.modifiers.append(modifier)
    data = pipeline.compute()
    rdf_alpha_beta=data.tables['coordination-rdf'].xy()
    max_rdf = rdf_alpha_beta[:,1].max()
    plt.plot(rdf_alpha_beta[:,0],rdf_alpha_beta[:,1])
    plt.title('RDF plot for Normal Stress = 0 and shear strain= {}'.format(beta))
    plt.xlabel('Pair Separation Distance')
    plt.ylabel('g(r)')
    plt.savefig('Clusters/{}_cluster/cluster_{}/RDF_step_{}'.format(no_atoms,iteration,i), facecolor='white', transparent=False)
    plt.clf()
    return max_rdf

    
def interatomic_dist_matrix(initial,no_atoms):
    import STZ
    import numpy as np
    import math
    interatomic_dist_mat = np.zeros([no_atoms,no_atoms])
    for i in range(no_atoms):
        for j in range(no_atoms):
            atom_i_coord = initial[i,:]
            atom_j_coord = initial[j,:]
            interatomic_dist_mat[i,j]= math.dist(atom_i_coord,atom_j_coord)
    return interatomic_dist_mat

def bond_matrix_config(alpha,beta,initial,thresh,no_atoms):
    import numpy as np
    import STZ
    from lammps import lammps
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "     
    lmp=lammps()
    initialization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001
    region myregion block 0.0 200.0 0.0 72.56 0.0 145.12  units box
    create_box  1 myregion
    mass 1 63.546

    pair_style eam
    pair_coeff 1 1 Cu_u6.eam '''
    lmp.commands_string(initialization_block)
    create_atoms=['create_atoms 1 single {} {} {}'.format(d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
    lmp.commands_list(create_atoms)
    surface_group=['group surface id {}'.format(surface_str)]
    lmp.commands_list(surface_group)

    minimization_block='''
    fix freeze surface setforce 0 0 0 
    minimize 0 1e-4 100000 100000
    unfix freeze
    '''
    lmp.commands_string(minimization_block)
    lmp.command('compute force all property/atom fx fy fz')
    lmp.command('dump fcal all custom 1 dump.force id type x y z fx fy fz')    
    lmp.command('run 1')
    coordinates,_=STZ.extract_fdump('dump.force',no_atoms)
    interatomic_dist_mat = STZ.interatomic_dist_matrix(coordinates,no_atoms)
    bonded_matrix = np.where(interatomic_dist_mat<thresh,1,0)
    return bonded_matrix


def Normal_Stress_Matrix_p(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc):
    from lammps import lammps
    import numpy as np
    import STZ

    ## NORMAL STRESS MATRIX GENERATION ##
    no_of_task = int(total_proc/proc_per_task)
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    reduced_alpha = d_alpha[140:300]  ## actual use alpha range
    # reduced_alpha = d_alpha ## actual use alpha range
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    batch_alpha = np.split(reduced_alpha,no_of_task)
    ns=np.zeros([int(len(reduced_alpha)/no_of_task),len(d_beta)])
    initial,type_particle=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    rev_arr = np.flip(batch_alpha[current_proc])
    y=-1
    # for i in rev_arr:
    for i in batch_alpha[current_proc]:
        alpha=i
        y=y+1
        for j in range(0,total_disp_steps+1):
            beta=(max_beta/total_disp_steps)*j
            deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
            d_initial=np.zeros([no_atoms,3])
            for k in range(no_atoms):
                d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
            x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)  
            initialization_block='''
            dimension 3
            units metal
            boundary s s s
            atom_style atomic
            timestep 0.001
            region myregion block {} {} {} {} {} {}  units box
            create_box 5 myregion
            mass 1 91.224
            mass 2 63.546
            mass 3 58.693
            mass 4 26.982
            mass 5 47.867 
            pair_style eam/alloy
            pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti
            '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

            create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            create_atoms_str = '\n'.join(str(e) for e in create_atoms)


            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            compute force all property/atom fx fy fz
            dump fcal all custom 1 dump.force id type x y z fx fy fz
            run 1 
            '''
            lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
            lmp = lammps()
            lmp.commands_string(lammps_input_script)

            coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
            ns[y,j],_=STZ.stress_calc(coordinates,force)
        np.savetxt('Clusters/{}_cluster/cluster_{}/Normal_stress_{}.txt'.format(no_atoms,iteration,current_proc),ns)

        

def Normal_Stress_Matrix_one(no_atoms,max_alpha,max_beta,total_disp_steps,iteration):
    from lammps import lammps
    import numpy as np
    import STZ

    ## NORMAL STRESS MATRIX GENERATION ##
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    reduced_alpha = d_alpha[130:300]  ## actual use alpha range
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    ns=np.zeros([int(len(reduced_alpha)),len(d_beta)])
    initial,type_particle=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    y=-1

    for i in reduced_alpha:
        alpha=i
        y=y+1
        for j in range(0,total_disp_steps+1):
            beta=(max_beta/total_disp_steps)*j
            deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
            d_initial=np.zeros([no_atoms,3])
            for k in range(no_atoms):
                d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
            x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)  
            initialization_block='''
            dimension 3
            units metal
            boundary s s s
            atom_style atomic
            timestep 0.001
            region myregion block {} {} {} {} {} {}  units box
            create_box  4 myregion
            mass 1 106.42 
            mass 2 63.546
            mass 3 58.693
            mass 4 30.974
            pair_style nep nep.population160.generation1013200.txt
            pair_coeff * * 
            '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

            create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            create_atoms_str = '\n'.join(str(e) for e in create_atoms)


            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            compute force all property/atom fx fy fz
            dump fcal all custom 1 dump.force id type x y z fx fy fz
            run 1 
            '''
            lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
            lmp = lammps()
            lmp.commands_string(lammps_input_script)

            coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
            ns[y,j],_=STZ.stress_calc(coordinates,force)
        np.savetxt('Clusters/{}_cluster/cluster_{}/Normal_stress_{}.txt'.format(no_atoms,iteration,current_proc),ns)


def box_coordinates (def_initial):
    margin = 2 # based on type of atom and its unit cell
    x_min = def_initial[:,0].min()-margin
    x_max = def_initial[:,0].max()+margin
    y_min = def_initial[:,1].min()-margin
    y_max = def_initial[:,1].max()+margin
    z_min = def_initial[:,2].min()-margin
    z_max = def_initial[:,2].max()+margin
    return x_min,x_max,y_min,y_max,z_min,z_max
def moving_average(a, n):
    import numpy as np
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def segments_fit(X, Y, count):
    from scipy import optimize
    import numpy as np
    import pylab as pl
    xmin = X.min()
    xmax = X.max()

    seg = np.full(count - 1, (xmax - xmin) / count)

    px_init = np.r_[np.r_[xmin, seg].cumsum(), xmax]
    py_init = np.array([Y[np.abs(X - x) < (xmax - xmin) * 0.01].mean() for x in px_init])

    def func(p):
        seg = p[:count - 1]
        py = p[count - 1:]
        px = np.r_[np.r_[xmin, seg].cumsum(), xmax]
        return px, py

    def err(p):
        px, py = func(p)
        Y2 = np.interp(X, px, py)
        return np.mean((Y - Y2)**2)

    r = optimize.minimize(err, x0=np.r_[seg, py_init], method='Nelder-Mead')
    return func(r.x)





def bond_matrix_stress_config(alpha,beta,initial,thresh,no_atoms):
    import numpy as np
    import STZ
    from lammps import lammps
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
  
    initialization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001
    region myregion block {} {} {} {} {} {}  units box
    create_box  1 myregion
    mass 1 63.546

    pair_style eam
    pair_coeff 1 1 Cu_u6.eam
    '''.format(x_min,x_max,y_min,y_max,z_min,z_max)
    
    create_atoms=['create_atoms 1 single {} {} {}'.format(d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
    create_atoms_str = '\n'.join(str(e) for e in create_atoms)

    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "    

    surface_group='\ngroup surface id {}'.format(surface_str)

    minimization_block='''
    fix freeze surface setforce 0 0 0 
    minimize 0 1e-4 100000 100000
    unfix freeze
    compute force all property/atom fx fy fz
    dump fcal all custom 1 dump.force id type x y z fx fy fz
    run 1 
    '''
    lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
    lmp = lammps()
    lmp.commands_string(lammps_input_script)
    coordinates,force=STZ.extract_fdump('dump.force',no_atoms)
    interatomic_dist_mat = STZ.interatomic_dist_matrix(coordinates,no_atoms)
    bonded_matrix = np.where(interatomic_dist_mat<thresh,1,0)
    ns,ss=STZ.stress_calc(coordinates,force)
    return bonded_matrix, ns ,ss

def interatomic_dist_matrix_config(alpha,beta,initial,no_atoms):
    import numpy as np
    import STZ
    from lammps import lammps
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    
    x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "     

    initialization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001
    region myregion block {} {} {} {} {} {}  units box
    create_box  1 myregion
    mass 1 63.546

    pair_style eam
    pair_coeff 1 1 Cu_u6.eam
    '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

    create_atoms=['create_atoms 1 single {} {} {}'.format(d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
    create_atoms_str = '\n'.join(str(e) for e in create_atoms)


    surface_id=STZ.surface_atoms(d_initial,no_atoms)
    surface_id_str=surface_id.astype('str')
    surface_str=" "
    for r in range(len(surface_id)):
        surface_str=surface_str+surface_id_str[r]
        surface_str+=" "    

    surface_group='\ngroup surface id {}'.format(surface_str)

    minimization_block='''
    fix freeze surface setforce 0 0 0 
    minimize 0 1e-4 100000 100000
    unfix freeze
    compute force all property/atom fx fy fz
    dump fcal all custom 1 dump.force id type x y z fx fy fz
    run 1 
    '''
    lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
    lmp = lammps()
    lmp.commands_string(lammps_input_script)
    coordinates,_=STZ.extract_fdump('dump.force',no_atoms)
    interatomic_dist_mat = STZ.interatomic_dist_matrix(coordinates,no_atoms)
    return interatomic_dist_mat




def initial_n_NN_list(initial,n):
    import STZ
    import numpy as np
    no_atoms = initial.shape[0]
    dist_mat = STZ.interatomic_dist_matrix_config(0,0,initial,no_atoms)
    n_NN = np.zeros([no_atoms,n+1])
    for i in range(no_atoms):
        per_atom_dist_sort_arg= dist_mat[i,:].argsort()+1
        n_NN[i,:] = per_atom_dist_sort_arg[0:n+1]
    n_NN = n_NN.astype(int)
    return n_NN



def dist_NN(ID1,ID2,no_atoms):
    import STZ
    import math
    coordinates,_=STZ.extract_fdump('dump.force',no_atoms)
    atom_1_coord = coordinates[ID1-1,:]
    atom_2_coord = coordinates [ID2-1,:]
    interatomic_dist_12= math.dist(atom_1_coord,atom_2_coord)
    return interatomic_dist_12

def Normal_Stress_Matrix_bond_len_multi(no_atoms,max_alpha,max_beta,total_disp_steps,iteration,total_proc,proc_per_task,current_proc):
    from lammps import lammps
    import numpy as np
    import STZ

    no_of_task = int(total_proc/proc_per_task)
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    batch_alpha = np.split(d_alpha,no_of_task)
    no_pass = len(batch_alpha[current_proc])*len(d_beta)
    initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    type_1_count = np.count_nonzero(particle_type==1) # Zr
    type_2_count = np.count_nonzero(particle_type==2) # Cu
    type_3_count = np.count_nonzero(particle_type==3) # Ni
    type_4_count = np.count_nonzero(particle_type==4) # Al
    type_5_count = np.count_nonzero(particle_type==5) # Ti


    dist_mat_1_1=np.zeros([type_1_count*type_1_count, no_pass])
    dist_mat_1_2=np.zeros([type_1_count*type_2_count, no_pass])
    dist_mat_1_3=np.zeros([type_1_count*type_3_count, no_pass])
    dist_mat_1_4=np.zeros([type_1_count*type_4_count, no_pass])
    dist_mat_1_5=np.zeros([type_1_count*type_5_count, no_pass])
    dist_mat_2_2=np.zeros([type_2_count*type_2_count, no_pass])
    dist_mat_2_3=np.zeros([type_2_count*type_3_count, no_pass])
    dist_mat_2_4=np.zeros([type_2_count*type_4_count, no_pass])
    dist_mat_2_5=np.zeros([type_2_count*type_5_count, no_pass])
    dist_mat_3_3=np.zeros([type_3_count*type_3_count, no_pass])
    dist_mat_3_4=np.zeros([type_3_count*type_4_count, no_pass])
    dist_mat_3_5=np.zeros([type_3_count*type_5_count, no_pass])
    dist_mat_4_4=np.zeros([type_4_count*type_4_count, no_pass])
    dist_mat_4_5=np.zeros([type_4_count*type_5_count, no_pass])
    dist_mat_5_5=np.zeros([type_5_count*type_5_count, no_pass])
    y=-1
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    for i in batch_alpha[current_proc]:
        alpha=i
        for j in range(0,total_disp_steps+1):
            y=y+1
            beta=(max_beta/total_disp_steps)*j
            deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
            d_initial=np.zeros([no_atoms,3])
            for k in range(no_atoms):
                d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
            x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
            initialization_block='''
            dimension 3
            units metal
            boundary s s s
            atom_style atomic
            timestep 0.001
            region myregion block {} {} {} {} {} {}  units box
            create_box 5 myregion
            mass 1 91.224
            mass 2 63.546
            mass 3 58.693
            mass 4 26.982
            mass 5 47.867 
            pair_style eam/alloy
            pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti
            '''.format(x_min,x_max,y_min,y_max,z_min,z_max)
            
            create_atoms=['create_atoms {} single {} {} {}'.format(particle_type[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
            create_atoms_str = '\n'.join(str(e) for e in create_atoms)

            minimization_block='''
            fix freeze surface setforce 0 0 0 
            minimize 0 1e-4 100000 100000
            unfix freeze
            compute force all property/atom fx fy fz
            dump fcal all custom 1 dump.force id type x y z fx fy fz
            dump_modify fcal sort id
            run 1 
            '''
            lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
            lmp = lammps()
            lmp.commands_string(lammps_input_script)

            initial_relaxed,_ = STZ.extract_fdump('dump.force',no_atoms)
            type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5 = STZ.init_interatomic_dist_matrix_multi(initial,particle_type,no_atoms)

            #Zr-Zr
            dist_mat_1_1[:,y] = type_1_1.flatten() 
            #Zr-Cu
            dist_mat_1_2[:,y] = type_1_2.flatten() 
            #Zr-Ni
            dist_mat_1_3[:,y] = type_1_3.flatten() 
            #Zr-Al
            dist_mat_1_4[:,y] = type_1_4.flatten() 
            #Zr-Ti
            dist_mat_1_5[:,y] = type_1_5.flatten() 
            #Cu-Cu 
            dist_mat_2_2[:,y] = type_2_2.flatten() 
            #Cu-Ni
            dist_mat_2_3[:,y] = type_2_3.flatten() 
            #Cu-Al
            dist_mat_2_4[:,y] = type_2_4.flatten() 
            #Cu-Ti
            dist_mat_2_5[:,y] = type_2_5.flatten() 
            #Ni-Ni
            dist_mat_3_3[:,y] = type_3_3.flatten() 
            #Ni-Al
            dist_mat_3_4[:,y] = type_3_4.flatten() 
            #Ni-Ti
            dist_mat_3_5[:,y] = type_3_5.flatten() 
            #Al-Al
            dist_mat_4_4[:,y] = type_4_4.flatten() 
            #Al-Ti
            dist_mat_4_5[:,y] = type_4_5.flatten() 
            #Ti-Ti
            dist_mat_5_5[:,y] = type_5_5.flatten() 

    
        
        np.savetxt('bond_len_1-1.txt',dist_mat_1_1)
        np.savetxt('bond_len_1-2.txt',dist_mat_1_2)
        np.savetxt('bond_len_1-3.txt',dist_mat_1_3)
        np.savetxt('bond_len_1-4.txt',dist_mat_1_4)
        np.savetxt('bond_len_1-5.txt',dist_mat_1_5)
        np.savetxt('bond_len_2-2.txt',dist_mat_2_2)
        np.savetxt('bond_len_2-3.txt',dist_mat_2_3)
        np.savetxt('bond_len_2-4.txt',dist_mat_2_4)
        np.savetxt('bond_len_2-5.txt',dist_mat_2_5)
        np.savetxt('bond_len_3-3.txt',dist_mat_3_3)
        np.savetxt('bond_len_3-4.txt',dist_mat_3_4)
        np.savetxt('bond_len_3-5.txt',dist_mat_3_5)
        np.savetxt('bond_len_4-4.txt',dist_mat_4_4)
        np.savetxt('bond_len_4-5.txt',dist_mat_4_5)
        np.savetxt('bond_len_5-5.txt',dist_mat_5_5)



def initial_NN_list_multi(initial,particle_type):
    import STZ
    import numpy as np
    no_atoms = initial.shape[0]
    type_1_count = np.count_nonzero(particle_type==1) # Zr
    type_2_count = np.count_nonzero(particle_type==2) # Cu 
    type_3_count = np.count_nonzero(particle_type==3) # Ni 
    type_4_count = np.count_nonzero(particle_type==4) # Al   
    type_5_count = np.count_nonzero(particle_type==5) # Ti   
    type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,sorted_initial_ID = STZ.init_interatomic_dist_matrix_multi(initial,particle_type,no_atoms)

    #Zr-Zr
    n_NN_1_1 = np.zeros([type_1_count,2])
    n_NN_1_1[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_1[i,:].argsort()
        n_NN_1_1[i,1] = sorted_initial_ID[per_atom_dist_sort_arg[1]]
    n_NN_1_1 = n_NN_1_1.astype(int)

    #Zr-Cu
    n_NN_1_2 = np.zeros([type_1_count,2])
    n_NN_1_2[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_2[i,:].argsort()
        n_NN_1_2[i,1] = sorted_initial_ID[type_1_count+per_atom_dist_sort_arg[1]]
    n_NN_1_2 = n_NN_1_2.astype(int)

    #Zr-Ni
    n_NN_1_3 = np.zeros([type_1_count,2])
    n_NN_1_3[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_3[i,:].argsort()
        n_NN_1_3[i,1] = sorted_initial_ID[type_1_count+type_2_count+per_atom_dist_sort_arg[1]]
    n_NN_1_3 = n_NN_1_3.astype(int)

    #Zr-Al
    n_NN_1_4 = np.zeros([type_1_count,2])
    n_NN_1_4[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_4[i,:].argsort()
        n_NN_1_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_1_4 = n_NN_1_4.astype(int)

    #Zr-Ti
    n_NN_1_5 = np.zeros([type_1_count,2])
    n_NN_1_5[:,0] = np.transpose(sorted_initial_ID[0:type_1_count])
    for i in range(type_1_count):
        per_atom_dist_sort_arg = type_1_5[i,:].argsort()
        n_NN_1_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_1_5 = n_NN_1_5.astype(int)

    #Cu-Cu
    n_NN_2_2 = np.zeros([type_2_count,2])
    n_NN_2_2[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_2[i,:].argsort()
        n_NN_2_2[i,1] = sorted_initial_ID[type_1_count+per_atom_dist_sort_arg[1]]
    n_NN_2_2 = n_NN_2_2.astype(int)

    #Cu-Ni
    n_NN_2_3 = np.zeros([type_2_count,2])
    n_NN_2_3[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_3[i,:].argsort()
        n_NN_2_3[i,1] = sorted_initial_ID[type_1_count+type_2_count+per_atom_dist_sort_arg[1]]
    n_NN_2_3 = n_NN_2_3.astype(int)

    #Cu-Al
    n_NN_2_4 = np.zeros([type_2_count,2])
    n_NN_2_4[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_4[i,:].argsort()
        n_NN_2_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_2_4 = n_NN_2_4.astype(int)

    #Cu-Ti
    n_NN_2_5 = np.zeros([type_2_count,2])
    n_NN_2_5[:,0] = np.transpose(sorted_initial_ID[type_1_count:type_2_count+type_1_count])
    for i in range(type_2_count):
        per_atom_dist_sort_arg = type_2_5[i,:].argsort()
        n_NN_2_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_2_5 = n_NN_2_5.astype(int)
       
    #Ni-Ni
    n_NN_3_3 = np.zeros([type_3_count,2])
    n_NN_3_3[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count:type_3_count+type_2_count+type_1_count])
    for i in range(type_3_count):
        per_atom_dist_sort_arg = type_3_3[i,:].argsort()
        n_NN_3_3[i,1] = sorted_initial_ID[type_1_count+type_2_count+per_atom_dist_sort_arg[1]]
    n_NN_3_3 = n_NN_3_3.astype(int)

    #Ni-Al
    n_NN_3_4 = np.zeros([type_3_count,2])
    n_NN_3_4[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count:type_3_count+type_2_count+type_1_count])
    for i in range(type_3_count):
        per_atom_dist_sort_arg = type_3_4[i,:].argsort()
        n_NN_3_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_3_4 = n_NN_3_4.astype(int)

    #Ni-Ti
    n_NN_3_5 = np.zeros([type_3_count,2])
    n_NN_3_5[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count:type_3_count+type_2_count+type_1_count])
    for i in range(type_3_count):
        per_atom_dist_sort_arg = type_3_5[i,:].argsort()
        n_NN_3_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_3_5 = n_NN_3_5.astype(int)

    #Al-Al
    n_NN_4_4 = np.zeros([type_4_count,2])
    n_NN_4_4[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count+type_3_count:type_4_count+type_3_count+type_2_count+type_1_count])
    for i in range(type_4_count):
        per_atom_dist_sort_arg = type_4_4[i,:].argsort()
        n_NN_4_4[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+per_atom_dist_sort_arg[1]]
    n_NN_4_4 = n_NN_4_4.astype(int)

    #Al-Ti
    n_NN_4_5 = np.zeros([type_4_count,2])
    n_NN_4_5[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count+type_3_count:type_4_count+type_3_count+type_2_count+type_1_count])
    for i in range(type_4_count):
        per_atom_dist_sort_arg = type_4_5[i,:].argsort()
        n_NN_4_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_4_5 = n_NN_4_5.astype(int)


    #Ti-Ti
    n_NN_5_5 = np.zeros([type_5_count,2])
    n_NN_5_5[:,0] = np.transpose(sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count:type_5_count+type_4_count+type_3_count+type_2_count+type_1_count])
    for i in range(type_5_count):
        per_atom_dist_sort_arg = type_5_5[i,:].argsort()
        n_NN_5_5[i,1] = sorted_initial_ID[type_1_count+type_2_count+type_3_count+type_4_count+per_atom_dist_sort_arg[1]]
    n_NN_5_5 = n_NN_5_5.astype(int)


    return n_NN_1_1,n_NN_1_2,n_NN_1_3,n_NN_1_4,n_NN_1_5,n_NN_2_2,n_NN_2_3,n_NN_2_4,n_NN_2_5,n_NN_3_3,n_NN_3_4,n_NN_3_5,n_NN_4_4,n_NN_4_5,n_NN_5_5

def init_interatomic_dist_matrix_multi(initial,particle_type,no_atoms):
    import STZ
    import numpy as np
    import math
    type_1_count = np.count_nonzero(particle_type==1) # Zr
    type_2_count = np.count_nonzero(particle_type==2) # Cu 
    type_3_count = np.count_nonzero(particle_type==3) # Ni 
    type_4_count = np.count_nonzero(particle_type==4) # Al   
    type_5_count = np.count_nonzero(particle_type==5) # Ti   
    type_1_1 = np.zeros([type_1_count,type_1_count])
    type_1_2 = np.zeros([type_1_count,type_2_count])
    type_1_3 = np.zeros([type_1_count,type_3_count])
    type_1_4 = np.zeros([type_1_count,type_4_count])
    type_1_5 = np.zeros([type_1_count,type_5_count])
    type_2_2 = np.zeros([type_2_count,type_2_count])
    type_2_3 = np.zeros([type_2_count,type_3_count])
    type_2_4 = np.zeros([type_2_count,type_4_count])
    type_2_5 = np.zeros([type_2_count,type_5_count])
    type_3_3 = np.zeros([type_3_count,type_3_count])
    type_3_4 = np.zeros([type_3_count,type_4_count])
    type_3_5 = np.zeros([type_3_count,type_5_count])
    type_4_4 = np.zeros([type_4_count,type_4_count])
    type_4_5 = np.zeros([type_4_count,type_5_count])
    type_5_5 = np.zeros([type_5_count,type_5_count])

    sort = np.argsort(particle_type)
    type_sorted_initial = initial[sort]
    type_sorted_initial_coord = type_sorted_initial[:,1:]
    for i in range(type_1_count):
        for j in range(type_1_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j,:]
            type_1_1[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_1_count):
        for j in range(type_2_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count,:]
            type_1_2[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_1_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_1_3[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_1_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_1_4[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_1_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_1_5[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_2_count):
        for j in range(type_2_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count,:]
            type_2_2[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_2_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_2_3[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_2_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_2_4[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_2_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_2_5[i,j]= math.dist(atom_i_coord,atom_j_coord)    
    
    for i in range(type_3_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_3_3[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_3_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_3_4[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_3_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_3_5[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_4_count):
        for j in range(type_4_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count+type_3_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count,:]
            type_4_4[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_4_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count+type_3_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_4_5[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_5_count):
        for j in range(type_5_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count+type_3_count+type_4_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count+type_3_count+type_4_count,:]
            type_5_5[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    return type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,sort

def interatomic_dist_matrix_tertiary(type_sorted_initial_coord,particle_type,no_atoms):
    import STZ
    import numpy as np
    import math
    type_1_count = np.count_nonzero(particle_type==1) # Pd
    type_2_count = np.count_nonzero(particle_type==3) # Ni  LAMMPS type 3
    type_3_count = no_atoms-type_1_count-type_2_count # P   LAMMPS type 4 
    type_1_1 = np.zeros([type_1_count,type_1_count])
    type_1_2 = np.zeros([type_1_count,type_2_count])
    type_1_3 = np.zeros([type_1_count,type_3_count])
    type_2_2 = np.zeros([type_2_count,type_2_count])
    type_2_3 = np.zeros([type_2_count,type_3_count])
    type_3_3 = np.zeros([type_3_count,type_3_count])

    for i in range(type_1_count):
        for j in range(type_1_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j,:]
            type_1_1[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_1_count):
        for j in range(type_2_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count,:]
            type_1_2[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_1_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_1_3[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_2_count):
        for j in range(type_2_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count,:]
            type_2_2[i,j]= math.dist(atom_i_coord,atom_j_coord)

    for i in range(type_2_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_2_3[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    for i in range(type_3_count):
        for j in range(type_3_count):
            atom_i_coord = type_sorted_initial_coord[i+type_1_count+type_2_count,:]
            atom_j_coord = type_sorted_initial_coord[j+type_1_count+type_2_count,:]
            type_3_3[i,j]= math.dist(atom_i_coord,atom_j_coord)
    
    
    return type_1_1,type_1_2,type_1_3,type_2_2,type_2_3,type_3_3


def yield_strain_multi_bs(max_alpha,max_beta,total_disp_steps,no_atoms,iteration):
    import STZ
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import interpolate

    ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    alpha_start, alpha_end = get_alpha_range(ns)
    d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_alpha = d_alpha[alpha_start:alpha_end]
    d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_1_4,initial_interatomic_mat_1_5,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_2_4,initial_interatomic_mat_2_5,initial_interatomic_mat_3_3,initial_interatomic_mat_3_4,initial_interatomic_mat_3_5,initial_interatomic_mat_4_4,initial_interatomic_mat_4_5,initial_interatomic_mat_5_5,sort_order = STZ.init_interatomic_dist_matrix_multi(initial,particle_type,no_atoms)
    initial_bonds = pos_count(initial_interatomic_mat_1_1)+pos_count(initial_interatomic_mat_1_2)+pos_count(initial_interatomic_mat_1_3)+pos_count(initial_interatomic_mat_1_4)+pos_count(initial_interatomic_mat_1_5)+pos_count(initial_interatomic_mat_2_2)+pos_count(initial_interatomic_mat_2_3)+pos_count(initial_interatomic_mat_2_4)+pos_count(initial_interatomic_mat_2_5)+pos_count(initial_interatomic_mat_3_3)+pos_count(initial_interatomic_mat_3_4)+pos_count(initial_interatomic_mat_3_5)+pos_count(initial_interatomic_mat_4_4)+pos_count(initial_interatomic_mat_4_5)
    cs = plt.contour(Beta,Alpha,ns,0)
    dat0=cs.allsegs[1][0]
    ss_lvl=np.zeros(len(dat0))
    s_strain=np.zeros(len(dat0))
    initial_bonded_matrix_1_1,initial_bonded_matrix_1_2,initial_bonded_matrix_1_3,initial_bonded_matrix_1_4,initial_bonded_matrix_1_5,initial_bonded_matrix_2_2,initial_bonded_matrix_2_3,initial_bonded_matrix_2_4,initial_bonded_matrix_2_5,initial_bonded_matrix_3_3,initial_bonded_matrix_3_4,initial_bonded_matrix_3_5,initial_bonded_matrix_4_4,initial_bonded_matrix_4_5,initial_bonded_matrix_5_5 = STZ.bond_mat_from_interatomic_dist_mat_multi(initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_1_4,initial_interatomic_mat_1_5,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_2_4,initial_interatomic_mat_2_5,initial_interatomic_mat_3_3,initial_interatomic_mat_3_4,initial_interatomic_mat_3_5,initial_interatomic_mat_4_4,initial_interatomic_mat_4_5,initial_interatomic_mat_5_5,3.921081,3.891297,3.449271,3.787830,3.767022,3.108552,3.362072,3.535373,3.209539,3.236423,3.587947,3.230980,3.362512,3.507643,3.524257)
    bond_change_frequency = np.zeros(len(dat0))
    dat_transfer = np.zeros(2)
    for j in range(len(dat0)):
        current_bonded_matrix_1_1,current_bonded_matrix_1_2,current_bonded_matrix_1_3,current_bonded_matrix_1_4,current_bonded_matrix_1_5,current_bonded_matrix_2_2,current_bonded_matrix_2_3,current_bonded_matrix_2_4,current_bonded_matrix_2_5,current_bonded_matrix_3_3,current_bonded_matrix_3_4,current_bonded_matrix_3_5,current_bonded_matrix_4_4,current_bonded_matrix_4_5,current_bonded_matrix_5_5,_,ss_lvl[j] = STZ.bond_matrix_stress_config_multi(dat0[j][1],dat0[j][0],initial,particle_type,no_atoms,sort_order,iteration)
        s_strain[j]=dat0[j][0]
        bond_status_change_matrix_1_1 = current_bonded_matrix_1_1-initial_bonded_matrix_1_1
        bond_status_change_matrix_1_2 = current_bonded_matrix_1_2-initial_bonded_matrix_1_2
        bond_status_change_matrix_1_3 = current_bonded_matrix_1_3-initial_bonded_matrix_1_3
        bond_status_change_matrix_1_4 = current_bonded_matrix_1_4-initial_bonded_matrix_1_4
        bond_status_change_matrix_1_5 = current_bonded_matrix_1_5-initial_bonded_matrix_1_5
        bond_status_change_matrix_2_2 = current_bonded_matrix_2_2-initial_bonded_matrix_2_2
        bond_status_change_matrix_2_3 = current_bonded_matrix_2_3-initial_bonded_matrix_2_3
        bond_status_change_matrix_2_4 = current_bonded_matrix_2_4-initial_bonded_matrix_2_4
        bond_status_change_matrix_2_5 = current_bonded_matrix_2_5-initial_bonded_matrix_2_5
        bond_status_change_matrix_3_3 = current_bonded_matrix_3_3-initial_bonded_matrix_3_3
        bond_status_change_matrix_3_4 = current_bonded_matrix_3_4-initial_bonded_matrix_3_4
        bond_status_change_matrix_3_5 = current_bonded_matrix_3_5-initial_bonded_matrix_3_5
        bond_status_change_matrix_4_4 = current_bonded_matrix_4_4-initial_bonded_matrix_4_4
        bond_status_change_matrix_4_5 = current_bonded_matrix_4_5-initial_bonded_matrix_4_5
        # bond_status_change_matrix_5_5 = current_bonded_matrix_5_5-initial_bonded_matrix_5_5

        initial_bonded_matrix_1_1 = current_bonded_matrix_1_1
        initial_bonded_matrix_1_2 = current_bonded_matrix_1_2
        initial_bonded_matrix_1_3 = current_bonded_matrix_1_3
        initial_bonded_matrix_1_4 = current_bonded_matrix_1_4
        initial_bonded_matrix_1_5 = current_bonded_matrix_1_5
        initial_bonded_matrix_2_2 = current_bonded_matrix_2_2
        initial_bonded_matrix_2_3 = current_bonded_matrix_2_3
        initial_bonded_matrix_2_4 = current_bonded_matrix_2_4
        initial_bonded_matrix_2_5 = current_bonded_matrix_2_5
        initial_bonded_matrix_3_3 = current_bonded_matrix_3_3
        initial_bonded_matrix_3_4 = current_bonded_matrix_3_4
        initial_bonded_matrix_3_5 = current_bonded_matrix_3_5
        initial_bonded_matrix_4_4 = current_bonded_matrix_4_4
        initial_bonded_matrix_4_5 = current_bonded_matrix_4_5
        # initial_bonded_matrix_3_3 = current_bonded_matrix_5_5
        bond_change_frequency[j] = STZ.neg_count(bond_status_change_matrix_1_1)+STZ.neg_count(bond_status_change_matrix_1_2)+STZ.neg_count(bond_status_change_matrix_1_3)+STZ.neg_count(bond_status_change_matrix_1_4)+STZ.neg_count(bond_status_change_matrix_1_5)+STZ.neg_count(bond_status_change_matrix_2_2)+STZ.neg_count(bond_status_change_matrix_2_3)+STZ.neg_count(bond_status_change_matrix_2_4)+STZ.neg_count(bond_status_change_matrix_2_5)+STZ.neg_count(bond_status_change_matrix_3_3)+STZ.neg_count(bond_status_change_matrix_3_4)+STZ.neg_count(bond_status_change_matrix_3_5)+STZ.neg_count(bond_status_change_matrix_4_4)+STZ.neg_count(bond_status_change_matrix_4_5)

    bond_breakage_frequency = np.delete(bond_change_frequency,0)
    norm_bond_breakage_freq = bond_breakage_frequency/initial_bonds
    window = 9 # (i-4 to i+4 : 9 strain steps)
    kernel = np.ones(window)
    norm_bond_breakage_density = np.convolve(norm_bond_breakage_freq, kernel, mode='valid')
    y_interp = interpolate.interp1d(s_strain,ss_lvl)
    yield_strain = dat0[np.where(norm_bond_breakage_density>0.1)[0][0]+ 4][0]
    tau_o_zero = y_interp(yield_strain)

    dat_transfer[0] = yield_strain
    dat_transfer[1] = tau_o_zero

    plt.figure(figsize=[6,4], dpi=300)
    plt.title('{} atoms STZ at normalised normal stress = 0.0'.format(no_atoms))
    ax1 = plt.subplot()
    color = 'tab:red'
    ax1.set_xlabel('Shear Strain')
    ax1.set_ylabel('Shear Stress (GPa)', color=color) 
    l1, = ax1.plot(s_strain,ss_lvl/10**9, color=color)
    ax1.plot(yield_strain, tau_o_zero/10**9, marker='x', markersize=10, color='r' )
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('Normalised Bond breaking Status Density', color=color)
    l2, = ax2.plot(s_strain[1:],bond_change_frequency, "-",color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.legend([l2, l1], ["Norm. Bond breaking Density", "Stress vs Strain"])
    plt.savefig('Clusters/{}_cluster/cluster_{}/bond_density_0.0.png'.format(no_atoms,iteration), facecolor='white', transparent=False)
    np.savetxt('Clusters/{}_cluster/cluster_{}/dat_trans.txt'.format(no_atoms,iteration), dat_transfer)
    plt.clf()   

def neg_count(matrix):
    neg_no = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if(matrix[i][j]<0):
                neg_no = neg_no+1
    return neg_no

def pos_count(matrix):
    pos_no = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if(matrix[i][j]>0):
                pos_no = pos_no+1
    return pos_no

def yield_strain_tertiary_bs_red(max_alpha,max_beta,total_disp_steps,no_atoms,iteration,no_segs):
    import STZ
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import interpolate
    import STZ
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import interpolate
    d_alpha = np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_alpha = d_alpha[130:300]
    d_beta = np.linspace(0,max_beta, num=total_disp_steps+1)
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_3_3, sort_order = STZ.init_interatomic_dist_matrix_tertiary(initial,particle_type,no_atoms)
    cs = plt.contour(Beta,Alpha,ns,0)
    dat0=cs.allsegs[1][0]
    ss_lvl=np.zeros(len(dat0))
    s_strain=np.zeros(len(dat0))
    initial_bonded_matrix_1_1,initial_bonded_matrix_1_2,initial_bonded_matrix_1_3,initial_bonded_matrix_2_2,initial_bonded_matrix_2_3,initial_bonded_matrix_3_3 = STZ.bond_mat_from_interatomic_dist_mat(initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_3_3,3.98,3.27,2.96,3.16,2.69,2.83)
    bond_change_frequency = np.zeros(len(dat0))
    dat_transfer = np.zeros(2)
    for j in range(len(dat0)):
        current_bonded_matrix_1_1,current_bonded_matrix_1_2,current_bonded_matrix_1_3,current_bonded_matrix_2_2,current_bonded_matrix_2_3,current_bonded_matrix_3_3,_,ss_lvl[j] = STZ.bond_matrix_stress_config_tertiary(dat0[j][1],dat0[j][0],initial,particle_type,no_atoms,sort_order)
        s_strain[j]=dat0[j][0]
        bond_status_change_matrix_1_1 = current_bonded_matrix_1_1-initial_bonded_matrix_1_1
        bond_status_change_matrix_1_2 = current_bonded_matrix_1_2-initial_bonded_matrix_1_2
        bond_status_change_matrix_1_3 = current_bonded_matrix_1_3-initial_bonded_matrix_1_3
        bond_status_change_matrix_2_2 = current_bonded_matrix_2_2-initial_bonded_matrix_2_2
        bond_status_change_matrix_2_3 = current_bonded_matrix_2_3-initial_bonded_matrix_2_3
        # bond_status_change_matrix_3_3 = current_bonded_matrix_3_3-initial_bonded_matrix_3_3
        initial_bonded_matrix_1_1 = current_bonded_matrix_1_1
        initial_bonded_matrix_1_2 = current_bonded_matrix_1_2
        initial_bonded_matrix_1_3 = current_bonded_matrix_1_3
        initial_bonded_matrix_2_2 = current_bonded_matrix_2_2
        initial_bonded_matrix_2_3 = current_bonded_matrix_2_3
        # initial_bonded_matrix_3_3 = current_bonded_matrix_3_3
        
        bond_change_frequency[j] = np.count_nonzero(bond_status_change_matrix_1_1)+np.count_nonzero(bond_status_change_matrix_1_2)+np.count_nonzero(bond_status_change_matrix_1_3)+np.count_nonzero(bond_status_change_matrix_2_2)+np.count_nonzero(bond_status_change_matrix_2_3)

    bond_change_frequency_trunc = np.delete(bond_change_frequency,0)
    cummulative_sum = 0
    bond_change_cumulative = np.zeros(len(bond_change_frequency_trunc))
    for i in range(len(bond_change_frequency_trunc)):
        cummulative_sum = cummulative_sum+bond_change_frequency[i+1]
        bond_change_cumulative[i] = cummulative_sum

    bond_change_cumulative_averaged= STZ.moving_average(bond_change_cumulative,4)
    bond_status_data =np.column_stack((dat0[0:-4,0],bond_change_cumulative_averaged))                                                               
    X, Y = bond_status_data[0:360,0], bond_status_data[0:360,1]
    px, py = STZ.segments_fit(X, Y, no_segs)
    thresh_zero = px[1]+0.025
    stress_diff = np.zeros(len(ss_lvl))
    for i in range(len(ss_lvl)-1):
        stress_diff[i] = ss_lvl[i+1]-ss_lvl[i]
    drop_strains = s_strain[stress_diff<0]
    drop_strains_useful = drop_strains[drop_strains< thresh_zero]
    if(len(drop_strains_useful)!=0):
        y_interp = interpolate.interp1d(s_strain,ss_lvl)
        tau_o_zero = y_interp(drop_strains_useful).max()
        strain_index = y_interp(drop_strains_useful).argmax()
        yield_strain = drop_strains_useful[strain_index]
    else:
        strains_useful = s_strain[s_strain<thresh_zero]
        y_interp = interpolate.interp1d(s_strain,ss_lvl)
        tau_o_zero = y_interp(strains_useful).max()
        strain_index = y_interp(strains_useful).argmax()
        yield_strain = strains_useful[strain_index]

    dat_transfer[0] = yield_strain
    dat_transfer[1] = tau_o_zero

    plt.figure(figsize=[6,4], dpi=300)
    plt.title('{} atoms STZ at normalised normal stress = 0.0'.format(no_atoms))
    ax1 = plt.subplot()
    color = 'tab:red'
    ax1.set_xlabel('Shear Strain')
    ax1.set_ylabel('Shear Stress (GPa)', color=color) 
    l1, = ax1.plot(s_strain,ss_lvl/10**9, color=color)
    ax1.plot(yield_strain, tau_o_zero/10**9, marker='x', markersize=10, color='r' )
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel('Cummulative Bond Status Changes', color=color)
    l2, = ax2.plot(bond_status_data[0:360,0], bond_status_data[0:360,1], ".",color=color,)
    ax2.plot(px, py, "-or")
    ax2.tick_params(axis='y', labelcolor=color)
    plt.legend([l2, l1], ["Bond status changes", "Stress vs Strain"])
    plt.savefig('Clusters/{}_cluster/cluster_{}/bond_freq_yield_0.0.png'.format(no_atoms,iteration), facecolor='white', transparent=False)
    np.savetxt('Clusters/{}_cluster/cluster_{}/dat_trans.txt'.format(no_atoms,iteration), dat_transfer)
    plt.clf()

    # d_alpha=np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    # d_beta=np.linspace(0,max_beta, num=total_disp_steps+1)
    # Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    # ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    # initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    # initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_3_3, sort_order = STZ.init_interatomic_dist_matrix_tertiary(initial,particle_type,no_atoms)
    # cs = plt.contour(Beta,Alpha,ns,0)
    # dat0=cs.allsegs[1][0]
    # ss_lvl=np.zeros(len(dat0))
    # s_strain=np.zeros(len(dat0))
    # initial_bonded_matrix_1_1,initial_bonded_matrix_1_2,initial_bonded_matrix_1_3,initial_bonded_matrix_2_2,initial_bonded_matrix_2_3,initial_bonded_matrix_3_3 = STZ.bond_mat_from_interatomic_dist_mat(initial_interatomic_mat_1_1,initial_interatomic_mat_1_2,initial_interatomic_mat_1_3,initial_interatomic_mat_2_2,initial_interatomic_mat_2_3,initial_interatomic_mat_3_3,3.98,3.27,2.96,3.16,2.69,2.83)
    # bond_change_frequency = np.zeros(len(dat0))
    # dat_transfer = np.zeros(2)
    # for j in range(len(dat0)):
    #     current_bonded_matrix_1_1,current_bonded_matrix_1_2,current_bonded_matrix_1_3,current_bonded_matrix_2_2,current_bonded_matrix_2_3,current_bonded_matrix_3_3,_,ss_lvl[j] = STZ.bond_matrix_stress_config_tertiary(dat0[j][1],dat0[j][0],initial,particle_type,no_atoms,sort_order)
    #     s_strain[j]=dat0[j][0]
    #     bond_status_change_matrix_1_1 = current_bonded_matrix_1_1-initial_bonded_matrix_1_1
    #     bond_status_change_matrix_1_2 = current_bonded_matrix_1_2-initial_bonded_matrix_1_2
    #     bond_status_change_matrix_1_3 = current_bonded_matrix_1_3-initial_bonded_matrix_1_3
    #     bond_status_change_matrix_2_2 = current_bonded_matrix_2_2-initial_bonded_matrix_2_2
    #     bond_status_change_matrix_2_3 = current_bonded_matrix_2_3-initial_bonded_matrix_2_3
    #     # bond_status_change_matrix_3_3 = current_bonded_matrix_3_3-initial_bonded_matrix_3_3
    #     initial_bonded_matrix_1_1 = current_bonded_matrix_1_1
    #     initial_bonded_matrix_1_2 = current_bonded_matrix_1_2
    #     initial_bonded_matrix_1_3 = current_bonded_matrix_1_3
    #     initial_bonded_matrix_2_2 = current_bonded_matrix_2_2
    #     initial_bonded_matrix_2_3 = current_bonded_matrix_2_3
    #     # initial_bonded_matrix_3_3 = current_bonded_matrix_3_3
        
    #     bond_change_frequency[j] = np.count_nonzero(bond_status_change_matrix_1_1)+np.count_nonzero(bond_status_change_matrix_1_2)+np.count_nonzero(bond_status_change_matrix_1_3)+np.count_nonzero(bond_status_change_matrix_2_2)+np.count_nonzero(bond_status_change_matrix_2_3)

    # bond_change_frequency_trunc = np.delete(bond_change_frequency,0)
    # cummulative_sum = 0
    # bond_change_cumulative = np.zeros(len(bond_change_frequency_trunc))
    # for i in range(len(bond_change_frequency_trunc)):
    #     cummulative_sum = cummulative_sum+bond_change_frequency[i+1]
    #     bond_change_cumulative[i] = cummulative_sum

    # bond_change_cumulative_averaged= STZ.moving_average(bond_change_cumulative,4)
    # bond_status_data =np.column_stack((dat0[0:-4,0],bond_change_cumulative_averaged))                                                               
    # X, Y = bond_status_data[0:360,0], bond_status_data[0:360,1]
    # px, py = STZ.segments_fit(X, Y, no_segs)
    # thresh_zero = px[1]+0.025
    # stress_diff = np.zeros(len(ss_lvl))
    # for i in range(len(ss_lvl)-1):
    #     stress_diff[i] = ss_lvl[i+1]-ss_lvl[i]
    # drop_strains = s_strain[stress_diff<0]
    # drop_strains_useful = drop_strains[drop_strains< thresh_zero]
    # if(len(drop_strains_useful)!=0):
    #     y_interp = interpolate.interp1d(s_strain,ss_lvl)
    #     tau_o_zero = y_interp(drop_strains_useful).max()
    #     strain_index = y_interp(drop_strains_useful).argmax()
    #     yield_strain = drop_strains_useful[strain_index]
    # else:
    #     strains_useful = s_strain[s_strain<thresh_zero]
    #     y_interp = interpolate.interp1d(s_strain,ss_lvl)
    #     tau_o_zero = y_interp(strains_useful).max()
    #     strain_index = y_interp(strains_useful).argmax()
    #     yield_strain = strains_useful[strain_index]

    # dat_transfer[0] = yield_strain
    # dat_transfer[1] = tau_o_zero

    # plt.figure(figsize=[6,4], dpi=300)
    # plt.title('{} atoms STZ at normalised normal stress = 0.0'.format(no_atoms))
    # ax1 = plt.subplot()
    # color = 'tab:red'
    # ax1.set_xlabel('Shear Strain')
    # ax1.set_ylabel('Shear Stress (GPa)', color=color) 
    # l1, = ax1.plot(s_strain,ss_lvl/10**9, color=color)
    # ax1.plot(yield_strain, tau_o_zero/10**9, marker='x', markersize=10, color='r' )
    # ax1.tick_params(axis='y', labelcolor=color)

    # ax2 = ax1.twinx()

    # color = 'tab:blue'
    # ax2.set_ylabel('Cummulative Bond Status Changes', color=color)
    # l2, = ax2.plot(bond_status_data[0:360,0], bond_status_data[0:360,1], ".",color=color,)
    # ax2.plot(px, py, "-or")
    # ax2.tick_params(axis='y', labelcolor=color)
    # plt.legend([l2, l1], ["Bond status changes", "Stress vs Strain"])
    # plt.savefig('Clusters/{}_cluster/cluster_{}/bond_freq_yield_0.0.png'.format(no_atoms,iteration), facecolor='white', transparent=False)
    # np.savetxt('Clusters/{}_cluster/cluster_{}/dat_trans.txt'.format(no_atoms,iteration), dat_transfer)
    # plt.clf()


def bond_matrix_stress_config_multi(alpha,beta,initial,type_particle,no_atoms,sorted_ID,iteration):
    import numpy as np
    import STZ
    from lammps import lammps
    text_file = open('Clusters/{}_cluster/cluster_{}/surface_id.txt'.format(no_atoms,iteration), "r")
    surface_group = text_file.read()
    deformation_grad=np.array([[1,0,beta],[0,1,0],[0,0,1+alpha]])
    d_initial=np.zeros([no_atoms,3])
    for k in range(no_atoms):
            d_initial[k,:]=np.matmul(deformation_grad,initial[k,1:])
    x_min,x_max,y_min,y_max,z_min,z_max = STZ.box_coordinates(d_initial)
    initialization_block='''
    dimension 3
    units metal
    boundary s s s
    atom_style atomic
    timestep 0.001
    region myregion block {} {} {} {} {} {}  units box
    create_box 5 myregion
    mass 1 91.224
    mass 2 63.546
    mass 3 58.693
    mass 4 26.982
    mass 5 47.867 
    pair_style eam/alloy
    pair_coeff * * ZrTiCuNiAl_Zhou04.eam.alloy Zr Cu Ni Al Ti
    '''.format(x_min,x_max,y_min,y_max,z_min,z_max)

    create_atoms=['create_atoms {} single {} {} {}'.format(type_particle[l],d_initial[l,0],d_initial[l,1],d_initial[l,2]) for l in range(len(initial[:,0]))]
    create_atoms_str = '\n'.join(str(e) for e in create_atoms)

    minimization_block='''
    fix freeze surface setforce 0 0 0 
    minimize 0 1e-4 100000 100000
    unfix freeze
    compute force all property/atom fx fy fz
    dump fcal all custom 1 dump.force id type x y z fx fy fz
    dump_modify fcal sort id
    run 1 
    '''

    lammps_input_script = initialization_block+create_atoms_str+surface_group+minimization_block
    lmp = lammps()
    lmp.commands_string(lammps_input_script)
    coordinates,force,particle_type=STZ.extract_fdump_type('dump.force',no_atoms)
    coordinates_sorted = coordinates[sorted_ID]
    type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,_ = STZ.init_interatomic_dist_matrix_multi(coordinates_sorted,particle_type,no_atoms)
    bond_mat_1_1,bond_mat_1_2,bond_mat_1_3,bond_mat_1_4,bond_mat_1_5,bond_mat_2_2,bond_mat_2_3,bond_mat_2_4,bond_mat_2_5,bond_mat_3_3,bond_mat_3_4,bond_mat_3_5,bond_mat_4_4,bond_mat_4_5,bond_mat_5_5 = STZ.bond_mat_from_interatomic_dist_mat_multi(type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,3.35036517,3.55296326,3.29013609,3.20843578,3.38126814,3.21116134,3.09788097,3.74233114,3.50770858,3.29733712,3.84690821,3.60998013,3.64419147,3.70396619,5)
    ns,ss=STZ.stress_calc(coordinates,force)
    return bond_mat_1_1,bond_mat_1_2,bond_mat_1_3,bond_mat_1_4,bond_mat_1_5,bond_mat_2_2,bond_mat_2_3,bond_mat_2_4,bond_mat_2_5,bond_mat_3_3,bond_mat_3_4,bond_mat_3_5,bond_mat_4_4,bond_mat_4_5,bond_mat_5_5, ns ,ss

def bond_mat_from_interatomic_dist_mat(type_1_1,type_1_2,type_1_3,type_2_2,type_2_3,type_3_3,thresh_1_1,thresh_1_2,thresh_1_3,thresh_2_2,thresh_2_3,thresh_3_3):
    import numpy as np
    bond_mat_1_1 = np.where(type_1_1<thresh_1_1,1,0)
    bond_mat_1_2 = np.where(type_1_2<thresh_1_2,1,0)
    bond_mat_1_3 = np.where(type_1_3<thresh_1_3,1,0)
    bond_mat_2_2 = np.where(type_2_2<thresh_2_2,1,0)
    bond_mat_2_3 = np.where(type_2_3<thresh_2_3,1,0)
    bond_mat_3_3 = np.where(type_3_3<thresh_3_3,1,0)

    return bond_mat_1_1,bond_mat_1_2,bond_mat_1_3,bond_mat_2_2,bond_mat_2_3,bond_mat_3_3

def bond_mat_from_interatomic_dist_mat_multi(type_1_1,type_1_2,type_1_3,type_1_4,type_1_5,type_2_2,type_2_3,type_2_4,type_2_5,type_3_3,type_3_4,type_3_5,type_4_4,type_4_5,type_5_5,thresh_1_1,thresh_1_2,thresh_1_3,thresh_1_4,thresh_1_5,thresh_2_2,thresh_2_3,thresh_2_4,thresh_2_5,thresh_3_3,thresh_3_4,thresh_3_5,thresh_4_4,thresh_4_5,thresh_5_5):
    import numpy as np
    bond_mat_1_1 = np.where(type_1_1<thresh_1_1,1,0)
    bond_mat_1_2 = np.where(type_1_2<thresh_1_2,1,0)
    bond_mat_1_3 = np.where(type_1_3<thresh_1_3,1,0)
    bond_mat_1_4 = np.where(type_1_4<thresh_1_4,1,0)
    bond_mat_1_5 = np.where(type_1_5<thresh_1_5,1,0)
    bond_mat_2_2 = np.where(type_2_2<thresh_2_2,1,0)
    bond_mat_2_3 = np.where(type_2_3<thresh_2_3,1,0)
    bond_mat_2_4 = np.where(type_2_4<thresh_2_4,1,0)
    bond_mat_2_5 = np.where(type_2_5<thresh_2_5,1,0)
    bond_mat_3_3 = np.where(type_3_3<thresh_3_3,1,0)
    bond_mat_3_4 = np.where(type_3_4<thresh_3_4,1,0)
    bond_mat_3_5 = np.where(type_3_5<thresh_3_5,1,0)
    bond_mat_4_4 = np.where(type_4_4<thresh_4_4,1,0)
    bond_mat_4_5 = np.where(type_4_5<thresh_4_5,1,0)
    bond_mat_5_5 = np.where(type_5_5<thresh_5_5,1,0)



    return bond_mat_1_1,bond_mat_1_2,bond_mat_1_3,bond_mat_1_4,bond_mat_1_5,bond_mat_2_2,bond_mat_2_3,bond_mat_2_4,bond_mat_2_5,bond_mat_3_3,bond_mat_3_4,bond_mat_3_5,bond_mat_4_4,bond_mat_4_5,bond_mat_5_5



def extract_fdump_type(path,no_atoms):
    import pandas as pd
    import numpy as np
    df=pd.read_csv(path,sep=" ",skiprows=9, nrows=no_atoms, header=None)
    df.columns=['ID', 'type', 'X','Y','Z', 'F_x', 'F_y', 'F_z']
    x=np.array(df.X)
    y=np.array(df.Y)
    z=np.array(df.Z)
    fx=np.array(df.F_x)
    fy=np.array(df.F_y)
    fz=np.array(df.F_z)
    typ = np.array(df.type)
    out_c=np.column_stack([x,y,z])                                 
    out_f=np.column_stack([fx,fy,fz]) 

    return out_c, out_f, typ


    

def MC_good_cluster_p_red(no_atoms,max_alpha,max_beta,total_disp_steps,lower,upper,steps,iteration,total_proc,proc_per_task,current_proc):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import signal, interpolate, stats
    from scipy.optimize import curve_fit
    import STZ
    no_of_task = int(total_proc/proc_per_task)
    dat_trans = np.loadtxt('Clusters/{}_cluster/cluster_{}/dat_trans.txt'.format(no_atoms,iteration))
    tau_o_zero = dat_trans[1]
    bs_change = dat_trans[0]


    factors_ns=np.linspace(lower,upper,steps)
    factors_ns = np.delete(factors_ns,int((steps-1)/2))
    factor_ns_batch = np.split(factors_ns,no_of_task)
    factor = factor_ns_batch[current_proc]
    normal_s = factor*tau_o_zero

    string='Normal_stress_'
    initial,particle_type=STZ.extract_dat_type('Clusters/{}_cluster/cluster_{}/cluster_{}.dat'.format(no_atoms,iteration,no_atoms),no_atoms)
    ns=np.loadtxt('Clusters/{}_cluster/cluster_{}/Normal_stress.txt'.format(no_atoms,iteration))
    alpha_start, alpha_end = get_alpha_range(ns)
    d_alpha = np.linspace(-max_alpha,max_alpha,total_disp_steps+1)
    d_alpha = d_alpha[alpha_start:alpha_end]
    d_beta = np.linspace(0,max_beta, num=total_disp_steps+1)
    Beta,Alpha=np.meshgrid(d_beta,d_alpha)
    cs = plt.contour(Beta,Alpha,ns,normal_s)
    dat0=cs.allsegs[0][0]
    ss_lvl=np.zeros(len(dat0))
    s_strain=np.zeros(len(dat0)) 
    for j in range(len(dat0)):
        _,ss_lvl[j]=STZ.stress_config_grad(dat0[j][1],dat0[j][0],initial,particle_type,no_atoms,iteration)
        s_strain[j]=dat0[j][0]
    string=string+str(round(normal_s[0]/tau_o_zero,4)) 
    ss_data = np.column_stack((s_strain, ss_lvl))


    stress_useful = ss_lvl[s_strain<(bs_change+0.01)]
    max_shear_lvl = stress_useful.max()
    strain_arg = stress_useful.argmax()
    max_yeild_strain = s_strain[strain_arg]

    plt.figure(figsize=[6,4], dpi=300)
    plt.title('{} atoms STZ at normalised normal stress = {}'.format(no_atoms,round(factor[0],4)))
    plt.plot(s_strain,ss_lvl/10**9)
    plt.plot(max_yeild_strain,max_shear_lvl/10**9, marker='x', markersize=10, color='r' )
    plt.xlabel('Shear strain')
    plt.ylabel('Shear Stress (GPa)')
    np.savetxt('Clusters/{}_cluster/cluster_{}/{}.txt'.format(no_atoms,iteration,string),ss_data)
    plt.savefig('Clusters/{}_cluster/cluster_{}/{}.png'.format(no_atoms,iteration,string), facecolor='white', transparent=False)
    plt.clf()
    string='Normal_stress_'

    MC_data=np.column_stack((factor[0],max_shear_lvl/tau_o_zero))
    np.savetxt('Clusters/{}_cluster/cluster_{}/MC_data_{}'.format(no_atoms,iteration,iteration),MC_data)


def get_alpha_range(mat):
    """
    Determine alpha_start and alpha_end based on number of rows in a 2D matrix.

    Parameters
    ----------
    mat : array-like
        2D matrix (NumPy array or similar)

    Returns
    -------
    alpha_start : int
    alpha_end   : int
    """

    n_rows = mat.shape[0]

    alpha_map = {
        160: (140, 300),
        136: (130, 266),
        170: (130, 300),
        100: (136, 236),
    }

    try:
        return alpha_map[n_rows]
    except KeyError:
        raise ValueError(f"Unsupported number of rows: {n_rows}")
