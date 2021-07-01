## _PYTHON VERSIONS and LIBRARIES_

Is requested the 3.8.5 version of python to run the code (or further updates)
Furthermore are requested the following libraries: numpy, matplotlib, configparser
## _STRUCTURE OF THE CODE_

The project is divided in different files

* [initial_condition](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_condition) contains all the input files
* [output_figures](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/output_figures) contains all the output figures
* [cfd.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/cfd.py) is the core of the code and it's the file to be launched
* [functions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/functions) contains all the functions previously discussed and other simple tools
* [functions_test](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/functions_test) is a tool to verify the correct behavior of the functions
* [input_test](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/input_test) is a tool to verify the correctness of the input files
* [test_functions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/test_functions) contains all the different tests used in [input_test](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/input_test) and [functions_test](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/functions_test)

## _HOW TO RUN THE CODE_

To run the code you have to provide seven text files, in the folder [initial_conditions](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions):
* [u_top.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/u_top.txt) : n*m matrix with the values of the initial **zonal current velocity of the upper layer**[m/s]
* [u_bot.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/u_bot.txt) : n*m matrix with the values of the initial **zonal current velocity of the lower layer**[m/s]
 * [v_top.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/v_top.txt) : n*m matrix with the values of the initial **meridional current velocity of the upper layer**[m/s]
* [v_bot.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/v_bot.txt) : n*m matrix with the values of the initial **meridional current velocity of the lower layer**[m/s]
* [u_wind.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/u_wind.txt) : n*m matrix with the values of the **zonal wind velocity**[m/s]
* [v_wind.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/v_wind.txt) : n*m matrix with the values of the **meridional wind velocity**[m/s]

* [const.txt](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/initial_conditions/const.txt) : contains all the **physical costants**, the grid dimension and the number of cells in each direction in the velocity grid (to verify that data are correctly privided)

The user must provide physically adequate values of variables:

* **Current velocity** must be below 1[m/s]
* **Wind velocity** must be below 30[m/s]
* **Surface elevation** must be below 2[m] and its mean value must be null
* **The grid** itself must count at least three points per dimension

To verify the correctnes of the provided data the user should compile [input_test.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/input_test.py).

Then and only then he will lauch the simulation compiling [cfd.py](https://github.com/enricomoresco/Software-and-Computing-Repository/blob/main/cfd.py).























