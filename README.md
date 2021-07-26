# Extragalactic Microlensing: magnification maps
This is the code to produce the magnification microlensing patterns that you can generate online at https://microlensing.overfitting.es/ . 

There you can find detailed explanations and links to scientific papers describing both the code and applications of the code. A direct link to the last paper about this code is https://arxiv.org/abs/2107.07222

The full list of the files in this project:
js/excanvas.min.js       // JS library to support ancient browsers (earlier than Internet Explorer 6)
js/flotr2.min.js         // JS library to draw
lib/mkl_dfti.f90         // three Fortran modules from Math Kernel Library used by the Poisson Solver
lib/mkl_poisson.f90
lib/mkl_trig_tranforms.f90
ml2.html                // HTML interface 
ml2.f90                 // source code 
style.css               // style file


INSTRUCTIONS:
In order to compile and link the corresponding Fortran files you need to install Intel Fortran Compiler at your computer and make the following steps:
1) Compile the lib/ .f90 files. In a terminal Linux machine would be:
    >  ifort -c mkl_poisson.f90 mkl_dfti.f90 mkl_trig_transforms.f90
2) Copy the main fortran program ml2.f90 in the same directory, compile and link it
    > ifort ml2.f90 -mkl -qopenmp -o ml2.cgi
3) Transfer that statically linked program ml2.cgi in your system cgi-bin/ directory.

The corresponding files in a Linux machine should be distributed as:
/var/www/html/ml2.html
/var/www/cgi-bin/ml2.cgi
/var/www/html/style.css
/var/www/html/js/flotr2.min.js
/var/www/html/js/excanvas.min.js

Please, cite our website and paper if you use this software in your work.


