# EB_20220223
Code and data for DKIST/ViSP image processing ( Data was taken on Feb. 23rd, 2022)

## Contents

+ IDL program to remove cross-talks from Stokes I 
+ IDL program to make an input profile for DeSIRe from data provided by the DKIST Data Center Archive
+ Run directory for a DeSIRe inversion

### IDL program to remove cross-talks from Stokes I

The program code is post_calib.pro

#### System requirements

+ Interactive Data Language (IDL) 8.5
+ The Solar Software package, which is available https://www.lmsal.com/solarsoft/
+ linux x86_64 m64

#### Installing guide

+ Download the EB_20220223
+ Typical install time is 10 sec

#### Demo

+ Issue a command '.r post_calib.pro' in the Solar Software IDL. The code remove cross-talk from Stokes I to the other Stokes parameters from the data in the directory 'demodata'.
+ Expected output is a directory for each data set, named 'the data set name'_cal01, under demodata directory.
+ Expected run time is 5 seconds.

#### Instructions for use

+ Download DKIST data in the fits format from the DKIST DATA Archive to the demodata directory, and issue the above command 




### IDL program to make an input profile for DeSIRe from data provided by the DKIST Data Center Archive

The program code file is all_post_calib.pro

#### System requirements

+ Interactive Data Language (IDL) 8.5
+ The Solar Software package, which is available https://www.lmsal.com/solarsoft/
+ linux x86_64 m64

#### Installing guide

+ Download the EB_20220223
+ Typical install time is 10 sec

#### Demo

+ N/A

#### Instructions for use

+ Download DKIST data in the fits format from the DKIST DATA Archive to the demodata directory, and issue a comman '.r all_post_calib'
+ EB_20220223/idl/make_data_08_txt/x025_y030.txt is the profile shown in Fig. 3 and 4



### Run directory for a DeSIRe inversion

#### System requirements

+ linux x86_64 m64
+ DeSIRe, which is available https://github.com/BasilioRuiz 

#### Installing guide

+ Download the EB_20220223
+ Put the directory 'run_desire' under 'run' directory of DeSIRe
+ Typical install time is 1 minute

#### Demo

+ N/A

#### Instructions for use

1. edit dir0 in step1.ipynb as your environment
2. run all cells of step1.ipynb
3. edit ref_dir and head_dir in step2.py as your environment
4. issue a command 'python step2.py'
5. Inversion results are in run_desire/results/
