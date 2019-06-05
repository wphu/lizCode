#lizCode

lizCode is a Particle(粒子) Code, including PIC, MCC and DSMC methods.

This code is developed based on Smilei v2.0:
https://github.com/SmileiPIC/Smilei
http://www.maisondelasimulation.fr/projects/Smilei/html/#

The main difference is that the electromagnetic field is changed to static electric field, and the possion solver is SuperLU for 2d, TDMA for 1d. Besides, some collisions (elastic collision, ionization collision and charge exchange collision) and PSIs (recycling and Backscattering) are included. A smiple plot tool based on Pyqt, matplotlib, hdf5 is developed.

#parallel compile
make -jN (N is number of processes)


data is mainly stored in ElectroMagn (density, temperature, electric potential...), Species (particles information) and Diagnostic (collision energy loss, psi_rate like sputtering rate, reflection rate...). psi_rate is not stored in PSI classes, which is for convenience of IO operation.


Email:  wphu@swip.ac.cn | wphu@qq.com


