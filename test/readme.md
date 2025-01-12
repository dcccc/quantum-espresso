In the N2_camb3lyp and ch4_camb3lyp directories, there are files as below:
```
cpmd.out
pw.in
pw.out
pw_div0.in
pw_div0.out
pw_m.in
pw_m.out
pw_m_div0.in
pw_m_div0.out
```

"pw.in" is the example of using cam-b3lyp

"pw_m.in" is the example of using CAM_B3LYP by setting the paramters of "exx_fraction", "exx_fraction_lr", "screening_parameter" with a input_dft of TUNED_CAM_B3LYP

"cpmd.in" is the input file used to calculate with cpmd[https://github.com/CPMD-code/CPMD], and the "cpmd.out" is the output file. Which is used to check the validation of the calculations with quantum espresso.

As the treatments of divergence at |G|=0 are different in cpmd and quantum espresso, to compare their results, the V_exx at |G|=0 are setted to 0. in both codes for the test calculation(the "cpmd.in", the "pw_div0.in", "pw_m_div0.in" and their correspoding output files).


In the N2_wb97mv directory, those are examples to using wb97mrv(id=768), b97mrv(id=769), wb97xrv(id=767) with the modified version libxc[https://github.com/dcccc/libxc]
```
pbe0_b97mrv.in
pbe0_b97mrv.out
pbe0_wb97mrv.in
pbe0_wb97mrv.out
pbe0_wb97xrv.in
pbe0_wb97xrv.out
```
