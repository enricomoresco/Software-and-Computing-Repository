## _PYTHON VERSIONS and LIBRARIES_

Is requested the 3.8.5 version of python to run the code (or further updates)
Furthermore are requested the following libraries: 
* numpy
* matplotlib 
* configparser
* hypotesis

## _STRUCTURE OF THE CODE_

The project is divided in different files

* [initial_conditions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/initial_conditions) is designed to contain all the input files
* [IC_case1](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case1) costitute an example case with uniform wind
* [IC_case2](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case2) costitute an example case with convergent wind
* [IC_case3](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case3) costitute an example case with uniform wind field rotor

* [output_figures](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/output_figures) contains all the output figures
* [cfd.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/cfd.py) is the core of the code and it's the file to be launched
* [functions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/functions) contains all the functions previously discussed and other simple tools
* [input_test](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/input_test) is a tool to verify the correctness of the input files
* [testing_functions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/test_functions) contains all the different tests
* [out_plot](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/out_plot) is designed to plot in [output_figures](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/output_figures) the output elevation and current field resulting, and the input wind field from the simulation
* [output_data](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/output_data) contains the output elevation and current field

## _WORK FLOW_


### first step: provide you input data
To run the code you have to provide three text files, in a folder, for example [initial_conditions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/initial_conditions):

* [u_wind.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/initial_conditions/u_wind.txt) : n*m matrix with the values of the **zonal wind velocity**[m/s]
* [v_wind.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/initial_conditions/v_wind.txt) : n*m matrix with the values of the **meridional wind velocity**[m/s]
* [const.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/initial_conditions/const.txt) : contains all the **physical costants**, the grid dimension and the number of cells in each direction in the velocity grid (to verify that data are correctly privided)

The user must provide physically adequate values of variables:

* **Wind velocity** must be below 30[m/s]
* **The grid** itself must count at least three points per dimension 


In order to help the user are already provided three test cases:

* [IC_case1](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case1)
* [IC_case2](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case2)
* [IC_case3](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case3)

these allow you to ignore the first step and go directly to the second one


### second step: verify your input data

To verify the correctnes of the provided data the user should compile [input_test.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/input_test.py).

_python input_test.py_

The program will ask you to insert the name of the IC_folder,

_Enter the name of the IC folder_

you will insert it, for example [IC_case1](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case1) costitute an example case with uniform wind

_Enter the name of the IC folder IC_case1_

when you'll read "ALL TESTS PASSED" you'll be ready to proceed to the third step


### third step: launch the simulation

Then and only then the user will lauch the simulation compiling [cfd.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/cfd.py).

_python cfd.py_

The program will ask you to insert the name of the IC_folder,

_Enter the name of the IC folder_

you will insert it, for example [IC_case1](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/IC_case1) costitute an example case with uniform wind

_Enter the name of the IC folder IC_case1_

and the simulation will run, when it will be complete will tell you

_simulation completed_

then you'll be ready for the fourth step

### fourth step: print your graphics

Finally the user can plot the current velocity field, the wind field and the surface elevation resulted from the simulation, he just has to lauch [out_plot.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/master/out_plot.py).

_python out_plot.py_

### fifth step: enjoy!

enjoy the program, make your own experiments and let me know if you need some more informatiosn than the previously given

Enrico Moresco

mailbox: enrico.moresco@studio.unibo.it


























