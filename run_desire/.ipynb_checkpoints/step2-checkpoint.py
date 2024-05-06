###############
# IMPORT

import numpy as np               # Needed to perform most of the mathemitcal operations
import matplotlib.pyplot as plt  # We will use it for plotting
from readwrite import readmod    # Needed for reading atmospheres in SIR/DeSIRe format
from readwrite import writemod   # Needed for writing atmospheres in SIR/DeSIRe format
from readwrite import readprof   # Needed for reading profiles in SIR/DeSIRe format
from readwrite import writeprof  # Needed for writing profiles in SIR/DeSIRe format
from readwrite import readrf     # Needed for reading Response Functions in SIR/DeSIRe
import sys                       # We can use it for exiting a program at any time
import os                        # This will allow us to run Linux commands from Python
import shutil
from astropy.io import fits
import time
import multiprocessing as multi

###############
# VARIABLES

# Parent and reference for inputs directory
ref_dir = 'hogehoge/desire/Source_code/desire20221219/run/run_desire/'

# Child directory
head_dir = 'hogehoge/desire/Source_code/desire20221219/run/tmp_'

# 
#nrun = 20
nrun = 10000


# Parameters for parallel computations
num_proc = 20

###############
# Make directories for results

if os.path.isdir(ref_dir+'results/') == False:
    os.mkdir(ref_dir+'results/')
if os.path.isdir(ref_dir+'results/imgs/') == False:
    os.mkdir(ref_dir+'results/imgs/')
if os.path.isdir(ref_dir+'results/models/') == False:
    os.mkdir(ref_dir+'results/models')
if os.path.isdir(ref_dir+'results/profiles/') == False:
    os.mkdir(ref_dir+'results/profiles/')

###############
# FUNCTIONS

def parallel(target, args, num_proc):
    t0 = time.time()
    jobs = []
    proc = 0
    rest = len(args)
    for arg in args:
        if proc < num_proc:
            p = multi.Process(target=target, args=(arg,))
            jobs.append(p)
            p.start()
            proc += 1
            rest -= 1
            if proc == num_proc or rest == 0:
                print('%s process was working, and rest process is %s.'%(num_proc, rest))
                for job in jobs:
                    job.join()
                proc = 0
                jobs = []
    t1 = time.time()
    print('{:.2f}'.format(t1-t0))
    
def make_run_directory(name_dir, ref_dir):
    if os.path.isdir(name_dir):
        print('delite '+name_dir)
        shutil.rmtree(name_dir)

    os.mkdir(name_dir)
    shutil.copy(ref_dir+'desire_inv.dtrol', name_dir)
    shutil.copy(ref_dir+'desire_syn.dtrol', name_dir)
    shutil.copy(ref_dir+'keyword.input', name_dir)
    shutil.copy(ref_dir+'PRD_H_1-0.dat', name_dir)
    shutil.copy(ref_dir+'PRD_H_2-0.dat', name_dir)

    os.mkdir(name_dir+'/__pycache__')
    shutil.copy(ref_dir+'__pycache__/readwrite.cpython-37.pyc', name_dir+'/__pycache__/')

    os.mkdir(name_dir+'/input')
    shutil.copy(ref_dir+'input/ASPLUND', name_dir+'/input/')
    shutil.copy(ref_dir+'input/atoms_backup.input', name_dir+'/input/')
    shutil.copy(ref_dir+'input/atoms.input', name_dir+'/input/')
    shutil.copy(ref_dir+'input/init_atmos.mod', name_dir+'/input/')
    shutil.copy(ref_dir+'input/LINES_LTE', name_dir+'/input/')
    shutil.copy(ref_dir+'input/LINES_NLTE', name_dir+'/input/')
    shutil.copy(ref_dir+'input/malla.grid', name_dir+'/input/')
    shutil.copy(ref_dir+'input/molecules.input', name_dir+'/input/')
    shutil.copy(ref_dir+'input/opacity_fudge.input', name_dir+'/input/')
    shutil.copy(ref_dir+'input/pos.coor', name_dir+'/input/')
    shutil.copy(ref_dir+'input/profiles0.per', name_dir+'/input/')
    shutil.copy(ref_dir+'input/profiles.per', name_dir+'/input/')
    shutil.copy(ref_dir+'input/profiles.wht', name_dir+'/input/')
    shutil.copy(ref_dir+'input/stray.per', name_dir+'/input/')
    shutil.copy(ref_dir+'input/wave.grid', name_dir+'/input/')
    shutil.copy(ref_dir+'input/wave_inv.grid', name_dir+'/input/')
    shutil.copy(ref_dir+'input/wave_syn.grid', name_dir+'/input/')
    shutil.copy(ref_dir+'input/atmos.mod', name_dir+'/input/')
    shutil.copy(ref_dir+'input/atmos2.mod', name_dir+'/input/')

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / stddev) ** 2 / 2)

def make_atmos(ref_dir, name_dir, perturb):

    # Read
    atm0,atm1 = readmod(ref_dir+'/input/atmos.mod')
    atm2,atm3 = readmod(ref_dir+'/input/atmos2.mod')

    # Perturvation
    stddev = 0.3
    tau = [-7,-6,-5,-4,-3,-2,-1,0]
    ntau = len(tau)
    ## Temperature
    dt  = [100.,100.,100.,100.,100.,100.,100.,100.]
    dv  = [1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5]
    db  = [100.,100.,100.,100.,100.,100.,100.,100.]
    for itau in range(0,ntau):
        atm1[:,1] = atm1[:,1] + gaussian(atm1[:,0], dt[itau]*perturb[0,itau],tau[itau], stddev)
        atm3[:,1] = atm3[:,1] + gaussian(atm3[:,0], dt[itau]*perturb[1,itau],tau[itau], stddev)
        atm1[:,5] = atm1[:,5] + gaussian(atm1[:,0], dv[itau]*perturb[2,itau],tau[itau], stddev)
        atm3[:,5] = atm3[:,5] + gaussian(atm3[:,0], dv[itau]*perturb[3,itau],tau[itau], stddev)
        atm1[:,4] = atm1[:,4] + gaussian(atm1[:,0], db[itau]*perturb[4,itau],tau[itau], stddev)
        atm3[:,4] = atm3[:,4] + gaussian(atm3[:,0], db[itau]*perturb[5,itau],tau[itau], stddev)
    
    # Write
    writemod(atm0,atm1,name_dir+'/input/atmos.mod')
    writemod(atm2,atm3,name_dir+'/input/atmos2.mod')
    
def plot_sav(prof_in,prof_out,outfile):
    
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize = (8,8), dpi=80)
    
    i=0
    ax[i].plot(prof_in[:,i+2], marker='.', color='black', linewidth=3, linestyle=' ', alpha=1, label='OBS')
    ax[i].plot(prof_out[:,i+2], color='red', linewidth=1, linestyle='-', alpha=1, label='OBS')

    for i in range(1,4):
        ax[i].plot(prof_in[:,i+2]/prof_in[:,2], marker='.', color='black', linewidth=3, linestyle=' ', alpha=1, label='OBS')
        ax[i].plot(prof_out[:,i+2]/prof_out[:,2], color='red', linewidth=1, linestyle='-', alpha=1, label='OBS')
    #ax[3].set_ylim(-0.003,0.003)
    fig.tight_layout()
    fig.savefig(outfile,dpi=300, bbox_inches = 'tight', pad_inches = 0.1)

def compare_files(file1, file2):
    # ファイル1の内容を読み取る
    with open(file1, 'r') as f1:
        content1 = f1.read()

    # ファイル2の内容を読み取る
    with open(file2, 'r') as f2:
        content2 = f2.read()

    # 内容が一致するかどうかを比較する
    if content1 == content2:
        res = True #print("ファイルの内容は同じです。")
    else:
        res = False #print("ファイルの内容は異なります。")

    return res
    
def fit_desire(args):
    ref_dir  = args['ref_dir']
    name_dir = args['name_dir']
    perturb  = args['perturb']
    
    # Preparations
    make_run_directory(name_dir, ref_dir)
    make_atmos(ref_dir, name_dir, perturb)
    name_dir_tail = name_dir[-6:-1]  #!!!!!
    ncycle=1
    switch_img = False #True    
   
    for i in range(0,9):
        if i > 1:
            np.random.seed(int(name_dir_tail)+i)
            make_atmos(ref_dir, name_dir, (np.random.rand(6, 8)-0.5)*2 )

        # Inversion
        shutil.copy(ref_dir+'input/wave_inv.grid', name_dir+'input/wave.grid')
        shutil.copy(ref_dir+'input/atoms_inv.input', name_dir+'input/atoms.input')
        os.chdir(name_dir)
        cmd = name_dir + '../../bin/desire '+name_dir + 'desire_inv.dtrol'
        os.system(cmd)

        # Synthesize
        shutil.copy(ref_dir+'input/wave_syn.grid', name_dir+'input/wave.grid')
        shutil.copy(ref_dir+'input/atoms_syn.input', name_dir+'input/atoms.input')
        os.chdir(name_dir)
        cmd = name_dir + '../../bin/desire '+name_dir + 'desire_syn.dtrol'
        os.system(cmd)
    
        prof_out = readprof(name_dir+'input/profiles.per')

        if (compare_files(name_dir+'input/profiles.per', name_dir+'input/profiles0.per') == False) & \
                (np.isnan(prof_out[0,2]) == False):
            break


    # Copy results
    shutil.copy(name_dir+'input/profiles.per', \
                ref_dir+'results/profiles/'+name_dir_tail+'.per')
    shutil.copy(name_dir+'input/atmos_'+str(ncycle).zfill(1)+'.mod', \
                ref_dir+'results/models/'+name_dir_tail+'_1.mod')
    shutil.copy(name_dir+'input/atmos2_'+str(ncycle).zfill(1)+'.mod', \
                ref_dir+'results/models/'+name_dir_tail+'_2.mod')


    # Save images
    if switch_img:
        prof_in  = readprof(name_dir+'input/profiles0.per')
        prof_out = readprof(name_dir+'input/profiles.per')
        plot_sav(prof_in,prof_out,ref_dir+'results/imgs/'+name_dir_tail+'.png')
    
    #sys.exit()
    # Remove the temporal run directory
    shutil.rmtree(name_dir) #test!!!!
    
#########################################################################
# Perturvations
ntau = 8
nparam = 6

# シード値の設定
np.random.seed(42)

perturbs = (np.random.rand(nparam, ntau, nrun)-0.5)*2
perturbs[:,:,0] = 0.
#perturbs[:,:,1] = 1. # test


#########################################################################
# Parallel computations
args = [
        {'ref_dir' : ref_dir,  \
         'name_dir': head_dir+str(i).zfill(5)+'/',         \
         'perturb' : perturbs[:,:,i] \
            }
        for i in range(0,nrun)
        ]
parallel(fit_desire, args, num_proc)





