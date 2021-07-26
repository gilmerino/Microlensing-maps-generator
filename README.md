# Extragalactic Microlensing: magnification maps
This is the code to produce the magnification microlensing patterns that you can generate online at https://microlensing.overfitting.es/ . 

There you can find detailed explanations and links to scientific papers describing both the code and applications of the code. A direct link to the last paper about this code is https://arxiv.org/abs/2107.07222

The full list of the files in this project:<br/>
js/excanvas.min.js      &nbsp // JS library to support ancient browsers (earlier than Internet Explorer 6)<br/>
js/flotr2.min.js        &nbsp // JS library to draw<br/>
lib/mkl_dfti.f90        &nbsp // three Fortran modules from Math Kernel Library used by the Poisson Solver<br/>
lib/mkl_poisson.f90<br/>
lib/mkl_trig_tranforms.f90<br/>
ml2.html                // HTML interface<br/>
ml2.f90                 // source code<br/>
style.css               // style file<br/>


INSTRUCTIONS:
In order to compile and link the corresponding Fortran files you need to install Intel Fortran Compiler at your computer and make the following steps:
1) Compile the lib/ .f90 files. In a terminal Linux machine would be:
    >  ifort -c mkl_poisson.f90 mkl_dfti.f90 mkl_trig_transforms.f90
2) Copy the main fortran program ml2.f90 in the same directory, compile and link it
    > ifort ml2.f90 -mkl -qopenmp -o ml2.cgi
3) Transfer that statically linked program ml2.cgi in your system cgi-bin/ directory.

The corresponding files in a Linux machine should be distributed as:<br/>
/var/www/html/ml2.html<br/>
/var/www/cgi-bin/ml2.cgi<br/>
/var/www/html/style.css<br/>
/var/www/html/js/flotr2.min.js<br/>
/var/www/html/js/excanvas.min.js<br/>

Please, cite our website and paper if you use this software in your work.


