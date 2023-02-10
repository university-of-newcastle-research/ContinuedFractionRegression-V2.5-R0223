# Analytic Continued Fractions for Regression

The Continued Fraction Regression (CFR) program Version is currently in version 2. Further details are in the following papers:

* [cfr-v1](https://ieeexplore.ieee.org/abstract/document/8789889/) - H. Sun and P. Moscato, "A Memetic Algorithm for Symbolic Regression," 2019 IEEE Congress on Evolutionary Computation (CEC), Wellington, New Zealand, 2019, pp. 2167-2174, DOI: 10.1109/CEC.2019.8789889.
* [cfr-v2](https://arxiv.org/abs/2001.00624) - Moscato, P., Sun, H. and Haque, M.N., 2019. Analytic Continued Fractions for Regression: A Memetic Algorithm Approach. arXiv preprint arXiv:2001.00624.

## Prerequisites & Requirements

These prerequisites assume a linux Ubuntu environment. This can be natively, via Windows Subsystem for Linux or via a Docker container. For initial assitance on a windows implementation, see the [Accompanying document](TestData/HowToRun-CFR-Windwos.pdf)

1. ```C++11``` as indicated in the Makefile 
```
CXX = g++
CXXFLAGS = -std=gnu++11 -O2 -Wall -Wextra -pedantic -g -fopenmp
LDFLAGS = -lgomp
```
Later versions not yet tested and dependent on support from the UoN HPC

2. ```gcc``` compiler for building the application
```
$ sudo apt install g++
```
To check the installation for GCC version, run following command in the terminal:
```
$ g++ --version
g++ (Ubuntu 7.2.0-18ubuntu2) 7.2.0
Copyright (C) 2017 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

3. ```make``` to orchestrate compilation
```
sudo apt update
sudo apt install make
```

4. ```Eigen``` library for linear regression and matrix functions
```
sudo apt install libeigen3-dev
```

Installs to ```/usr/include/eigen3``` which is specified in the Makefile with ```-I <Eigen Path>```

## Program Execution

### Build Process

To build the executable of the program we have to create a directory ```bin``` (if does not exist).
Then we have to ```clean``` any previous builds. Now the program can be built.
The following sequence of commands are required to build the program for the first time:
```
mkdir bin
make clean
make
```
Make will generate the executable by compiling the source files into the ```bin``` directory.

### Execution

The program can be run with ```./bin/main```. The latest command-line options can be viewed with the --help flag
```
./bin/main --help
```

The command line arguments as of Sep-2021

```
VERSION:
    Continued Fraction Regression with Matrix Representation.
    Last update: 05-05-2021
    CFR version: 2.5.1

Complete Command-line Arguments:./bin/main --help 
Continued Fraction Regression.

            REQUIRED ARGUMENTS:
        -t --train         Training Data File location.

    OPTIONS:
        -b --base           If the representation for the model is CCL, selects the base used (Integer 0 < b <= 10000)
                            Default = 2

        --cont              Data used to eval model in Local Search is continuous instead of randomly selected (String)
                            [Available Options: true, false]
                            Default = false

        -d --delta          Delta/Penalty in fitness function. (Real Number 0.0 <= d <= 1.0)
                            Default=0.35.

        -f --frac-depth     Depth of the Continued Fraction (Integer 0 <= f <= 20)
                            Default = 4

        -g --num-gen        Number of Generations (Integer 1 <= g <= 10000)
                            Default = 200

        -h --help           Show this screen.

        -io,                Only allow integer coefficients (no argument)
        --interger-only

        -it,                Force the terms f(x) = {gi(x), hi(x)} to resolve to integer values. An integer argument must be provided:
        --integer-terms         0 = integers are not forced
                                1 = floor(f(x)) is take
                                2 = round(f(x)) is taken
                                3 = ceil(f(x)) is taken

                            Local search behaviour may be modified based on this parameter to focus on the exploration of the natural numbers 
                            instead of the reals

        -l --ls-gen         Run Local Search at each l(th) gen (Integer 1 <= l <= g)
                            Default = 1

        -id --inc-depth     Flag (no argument) to enable Incremental Depth Approach where the model from previous depth 
                            will be used to Initialise the Population of current depth.

        -pmf --pr-model-mf  Flag (no argument) to enable the usage of model form previous Depth as a meta-feature in current depth.
                            The `-id` or `--inc-depth` must be enabled to use this feature.

        -ld --ls-data       Percentage of Samples used to eval model in Local Search. (1 <= ld <= 100)
                            Default = 100.

        -loo --loo-cv       Enable Leave-One-Out Cross Validation on the Training Data. No Parameter Values if required. By default it is disabled.

        -m --mut-rate       Mutation Rate (Real Number 0.0 <= m <= 1.0)
                            Default = 0.20

        -n --norm           Data Normalisation Method (Integer 0 <= n <=3)
                            [Available Options: 1=MinMax, 2=z-Score, 3=Decimal Scaling]
                            Default = 0 (No Normalisation).
    
        -nm --nelder-mead   Parameters for Nelder-Mead (NM) (all Integers): <nr>:<gen>:<reset>
                            nr=number of runs of NM (1 to 100), gen=max generations per run of NM (1 to 1000)
                            reset=NM reset after stagnate for consecutive generations (1 to 100).

        -o --obj            Objective Function (String)
                            [Available Options: 
                                mse, 
                                nsme, 
                                pearson, 
                                theil, 
                                med.err, 
                                med.ae, 
                                mase, 
                                ub lb, 
                                mwe, 
                                uwae,
                                cw,
                                mre,
                                aic,
                                mwrae
                                ]
                            Default = mse.

        -or --origin       Force the solution to pass through the origin point for x1,x2,..,xN = 0


        -p --pop-init       Population Initiliastion Method.
                            [Available Options: avg (Initialise linear term of solution with average of the target score),
                                                sol <file-name of previous solution> (Initialise population with a known solution)
                                                linreg <OPTIONAL: file-name of previous solution> (Initialise with linear regression, optionally starting from known solution)]

        -prec --out-prec    Set Precision in Output. (1 <= op <= 18)
                            Default = 6.

        -r --reset-root     Reset Root of population tree after stuck for r generations. (Integer 1 <= r <= g)
                            Default =  5

        -s --seed           Seed for Random Number Generator. (Positive Integer)
                            Default = 13.

        -sw --samp-weight   Flag (no argument) to indicate the last column is sample weight
                            The variable count in the dataset file must be (feature count)+1 when weight supplied e.g.

                            Three features      Three features + sample weight 
                            y,f1,f2,f3          y,f1,f2,f3,w
                            1,0,0,1             1,0,0,1,5
                            2,0,1,0             2,0,1,0,3
                            3,0,1,1             3,0,1,1,5

                            NB  The test data has no purpose for weights, but test and training files are handled the same 
                            Currently, you must include the weight column ```w``` in the test .csv file.
                            The code will ignore this column for the testing data

        -pc --param-count  AIC: path to the file containing the parameter-counter when meta-feature is present.
        -loo --loo-cv      Enable Leave-One-Out Cross Validation on the Training Data. No Parameter Values if required. By default it is disabled.

        -S --is-serial      Run Local Search Optimisation in Serial Mode? (String)
                            [Available Options: true, false]
                            Default = true.
        
        -t --train          Training Data File location. (Valid file Path)

        -T --Test           Testing Data File location. (Valid file Path)
        
        -v --verbose        Verbosity mode of output. (Integer: 0 <= v <=4)
                            Default = 0

        -V --version        Show version.

        -w --window         Set the size of window (in percentage of training data) to use in local search. (1 <= w <= 100)
                            Default = 100.

        -y --type           Type of representation for model output (Integer)
                            [Available Options: 0 for Continued Fraction (CF), 1 for Coefficient Continued Logarithm (CCL)]
                            Default = 0
        
    USAGE:
        .bin/main -t <train-file>
        .bin/main -t <train-file> -T <test-file>
        .bin/main -t <train-file> -T <test-file> -s 107
        .bin/main -t <train-file> -T <test-file> -g 200 -m 0.1 -nm 4:200:5 -v 2
        .bin/main -t <train-file> -T <test-file> -g 200 -m 0.5 -d 0.2 -nm 4:250:10 -f 0 -p avg -v 2
        .bin/main -t <train-file> -T <test-file> -g 200 -m 0.1 -d 0.1 -nm 4:250:10 -f 4 -p sol CFR-v2-Depth-3.sol.txt -v 2
        .bin/main --version

```
An Example of how to run the program (with default 80-20 split of supplied dataset):
```
./bin/main  -t TestData/Argon.csv
```

Another Example: running the program with parameters used in CEC 2020 paper (fitting the function on entire data): [obj=mse, gen=200, muRate=0.10  penalty=0.1 nelder-mead:4:250:10 depth=4]:
```
./bin/main -t TestData/Argon.csv -T TestData/Argon.csv -g 200 -m 0.1 -d 0.1 -nm 4:250:10 -f 4 -v 2
```

## Data Format

The program requires the data to be in Comma Separated Value file (with ```.csv``` extension). The first line of the data file contains the header of variables separated by comma ```,``` character. There are some reserved name to specify some special type of variable/feature in the dataset. They are:
* targetHeader ```y```: Specify the target/dependent variable of the dataset
* weightHeader ```w```: To specify the weight associated for the target of each samples (for weighted learning)
* uncertantyHeader ```dy```: Is required to use Chi-Squared Weighted Error Metric in the Fitness computation.



# Other Information

### Linear Regression Changes

Most of the relevant linear regression code is located in `linear.hpp` and `linear.cpp`, and particular usage is documented in those files.

Currently, linear regression is used in two different ways in the memetic algorithm:

1. Initialization of the population: For every agent, one of the continued fractions (either the pocket or the current, but not both) are initialized using linear regression as follows:

   1.1 The root agent has a fraction initialized completely by linear regression.
   1.2. For non-root agents, a random number `d` is chosen between 0 and the maximum depth of the fraction. Linear regression is used to fit the fraction starting from depth `d` onwards, while the rest of the fraction is initialized as normal. The variables used are randomly chosen in the same way as normal.

   This is enabled through the new initialization option `-p linreg`. Optionally, a solution file may be provided `-p linreg <solution file>`, and the root will have one of its fractions initialized from that solution file. If the fraction depth is greater than the depth of the solution file's fraction, the rest of the fraction is completed by linear regression. The relevant code is in `init.cpp`, in the function `population::initialize_linear_root()` and `population::initialize_linear()`

2. Linear regression before local search: Right before performing nelder-mead, but after mutation and recombination, linear regression will be used to fit the *last* depth of the continued fraction using the variables selected by mutation and recombination. This occurs with a probability `g_linreg_chance`, which is set through the new option `-lr` or `--linreg-prob` followed by a real number between `0.0` and `1.0`, representing the probability. The default chance is `0.0`. It seems that `-lr 0.5` generally works well in practice. The relevant code is in `optimize.cpp`. 

### Run via Docker
A docker container can be built from the Dockerfile within the root directory. This requires Docker Hub is installed on your OS, and within Windows, review the use of Windows Subsystem for Linux (WSL 2).

To build the container, including all required libraries/packages
```
docker build -t mlab/cfr:v2 .
```

To run the container with a terminal (mounting your current directory to root)
```
docker -it -v  $(pwd):/root mlab/cfr:v2
```

From here, navigate to /root and compile as per the build process.
Note that changes to container will only persist if made in /root as the directory is mounted. All other changes (e.g. installation of packages) is refreshed on exit from the terminal container. To persist these changes or configurations, update the Dockerfile



### Increase the Capacity of the Maximum Number of Features \& Maximum Depth
In this matrix-based implementation, we have defined the matrix with a size of ```1000``` features (```N```) and Depth (```M```) of Fraction as ```20```. These values (```M, N```) can be changed at lines 23 and 24 of ```cont_frac.hpp```. If we want to increase the ```depth=30```, then we have to put ```M = 2 * Depth + 1 = 2 * 30 + 1 = 61```. We also have to update the range in ```main.cpp``` with the new value of depth in place of ```20``` in the line 315: ```                if (is_between(g_depth, (size_t) 0, (size_t) 20) == true) {
```. Then we have to ```clean``` and ```make``` the program. 

```
class MatrixContinuedFraction {
    // M x N matrix
    #define M 41    // terms = 2*depth+1, here depth=20
    #define N 1002   // num_vars
```


#### Some dataset to test the program
The ```TestData\``` folder contains some datasets in the required ```.csv``` format:
* **[Argon](TestData\Argon.csv)** : Argon dataset has been used in "Kulakova, L., Arampatzis, G., Angelikopoulos, P. et al. Data driven inference for the repulsive exponent of the Lennard-Jones potential in molecular dynamics simulations. Sci Rep 7, 16576 (2017). https://doi.org/10.1038/s41598-017-16314-4". The version of the Argon data used in the experiments can be found in **Table-2** of the ***Supporting Information*** of "Arthur M. Halpern, Structural and Thermodynamic Properties of the Argon Dimer, Journal of Chemical Education 2010 87 (2), 174-179, DOI: https://doi.org/10.1021/ed800049s."
* **[airfoil](TestData/Airfoil.csv)**
* **[concrete](TestData/Concrete.csv)**
* **[cooling](TestData/Cooling.csv)**
* **[heating](TestData/Heating.csv)**
* **[housing](TestData/Housing.csv)**
* **[yacht](TestData/Yacht.csv)**
* **[PresidentElect](TestData/PresidentElect.csv)**: This data is used in "Moscato, P., Mathieson, L., Mendes, A., & Berretta, R. (2005). The electronic primaries: predicting the U.S. Presidency using feature selection with safe data reduction. In V. Estivill-Castro (Ed.), Computer Science 2005 - 28th Australasian Computer Science Conference, ACSC 2005 (Vol. 38, pp. 371-380)
".
* **[Redshift](TestData/SpaceObj-Redshift)**: The redshift (SpaceObj) dataset is used in the publication of "A. Krone-Martins, E. E. O. Ishida, R. S. de Souza, The first analytical expression to estimate photometric redshifts suggested by a machine, Monthly Notices of the Royal Astronomical Society: Letters, Volume 443, Issue 1, 1 September 2014, Pages L34–L38, https://doi.org/10.1093/mnrasl/slu067". The data folder ```SpaceObj-Redshift``` contains two datasets: ```SpecObj``` is the original dataset with 15 features and ```SpecObjSingleVal``` is the modified version of the Redshift-data which contains only variables with a single character name.
* **[Superconductivity](TestData/Superconductivity)**: The dataset is used in "Hamidieh, Kam, A data-driven statistical model for predicting the critical temperature of a superconductor, Computational Materials Science, Volume 154, November 2018, Pages 346-354", can be accessed from the UCI-Machine Learning Repository at: https://archive.ics.uci.edu/ml/datasets/Superconductivty+Data
* **[BodyFat](TestData/BodyFat)**: This version of the BodyFat dataset has two additional variables (BMI and Waist) than the original dataset at: http://staff.pubhealth.ku.dk/~tag/Teaching/share/data/Bodyfat.html. The value of Body Mass Index (BMI) can be calculated as ``` BMI = Weight/(Height/100)^2```.From Pablo Moscato to Everyone:  10:37 AM
https://people.maths.ox.ac.uk/trefethen/bmi.html
https://www.medicaldaily.com/oxford-mathematician-explains-body-mass-index-flaw-244342


### Error Metrics in Fitness Computation
The program has the facility to use different types of error metric. Those are as follows:
* **MSE**: It uses the Mean Absolute Error (MSE) to compute the fitness score as ```MSE * (1 + used_features() * delta)```. The ```delta``` can be set to the desired value (suppose ```0.1```) by ```-d 0.1```. This is the default error metric, however, can be specified by ```-o mse```.
* **NMSE**: The Normalised MSE can be set by ```-o nmse```. It also uses the ```delta``` to compute the fitness as ```NMSE * (1 + used_features() * g_delta)```
* **Pearson**: The Pearson correlation coefficient can be used as an error metric by ```-o pearson```. It is computed as ```(1 - pearson)```.
* **Median Error**: The median of error score can be used as the error metric by ```-o med.err```.  Here, the error is computed as ```(target-prediction)```.
* **Median of Absolute Error**: The median of absolute error score can be used as the error metric by ```-o med.ae```.  Here, the error is computed as ```abs(target-prediction)```.
* **Mean of Absolute Scaled Error**: The Mean of absolute scaled error (MASE) score can be used as the error metric by ```-o mase```.  Here, the error is computed as described in Wikipedia Article at [https://en.wikipedia.org/wiki/Mean_absolute_scaled_error] for the Non-seasonal time series proposed in Hyndman, R. J. and Koehler A. B. (2006). "Another look at measures of forecast accuracy." International Journal of Forecasting volume 22 issue 4, pages 679-688 [https://doi.org/10.1016%2Fj.ijforecast.2006.03.001].
* **EMC**: The excursion matching criteria (EMC) score can be used as the error metric by ```-o emc```.  Here, it is computed as ```EMC =  \frac{\sum x_i \cdot y_i}{\sum |x_i| \cdot |y_i|}```, where ```x=prediction``` and ```y=target```. It is taken from [https://contextearth.com/2017/10/25/improved-solver-target-error-metric/].

* **Lower Bound**: The lower bound can be computed by specifying ```-o lb```.  If ```(y-f(x)) >= 0 ``` then ```err = 1+ abs(y-f(x))```, else ```err = (1+ (y-f(x)))^8```
* **Upper Bound**: The upper bound can be computed by specifying ```-o ub```.  If ```(y-f(x)) <= 0 ``` then ```err = 1+ abs(y-f(x))```, else ```err = (1+ (y-f(x)))^8```

* **Chi-Squared Weighted**: To use the Ch-squared error metric incorporating the weight/uncertainity we have to use ```-o cw```.   This erro metric is used in Nguyen, H. (2020). Analyzing Pantheon SNeIa data in the context of Barrow's variable speed of light. arXiv preprint arXiv:2010.10292. [https://arxiv.org/abs/2010.10292]. The euqation (61) defines for &sigma;=Uncertainity as:
> $$
> \chi^2=\frac{1}{n} \sum_{i=1}^{n} \frac{1}{\sigma_{i}^2}(y_{obs}^{i} - y_{pred}^{i})^2
> $$

* **Akaike information criterion (AIC)**: To use the AIC error metric we have to use ```-o aic```. To incorporate the meta-feature complexity we have specify the path ```-pc <path to model complexity file>```.  This AIC error metric is computed  as:
> $$
AIC=\begin{cases}n\;ln\left(\frac{RSS}{n}\right)+2K &\text{if } n/K\geq 40\\n\;ln\left(\frac{RSS}{n}\right)+2K+\frac{2K(K+1)}{n-K-1}&\text{otherwise}\end{cases}
> $$
> > Were ```n=number of samples``` and  ```K=number of free parameters in the model```.
</br>
You have to put the free paraemter number for the features as ```0``` and include the number of free parameters in the meta-feature in file supplied by ```-pc``` separated by comma in the sequential order of input data file.


## Publications
* [1]   H. Sun and P. Moscato, "A memetic algorithm for symbolic regression," in 2019 IEEE Congress on Evolutionary Computation (CEC), 2019: IEEE, pp. 2167-2174. [View Article](https://ieeexplore.ieee.org/document/8789889) 
* [2]   P. Moscato, H. Sun, and M. N. Haque, "Analytic Continued Fractions for Regression: A Memetic Algorithm Approach," Expert Systems with Applications, 2021, 115018. [View Article](https://doi.org/10.1016/j.eswa.2021.115018) 
* [3]   P. Moscato, M. N. Haque, K. Huang, J. Sloan, and J. C. de Oliveira, "Learning to extrapolate using continued fractions: Predicting the critical temperature of superconductor materials," arXiv preprint arXiv:2012.03774, 2020. [View Article](https://doi.org/10.48550/arXiv.2012.03774) 
* [4]   P. Moscato, H. Sun, and M. N. Haque, "Analytic Continued Fractions for Regression: Results on 352 datasets from the physical sciences," in 2020 IEEE Congress on Evolutionary Computation (CEC), 2020: IEEE, pp. 1-8. [View Article](https://doi.org/10.1109/CEC48606.2020.9185564)
* [5]   P. Moscato et al., "Multiple regression techniques for modeling dates of first performances of Shakespeare-era plays," arXiv preprint arXiv:2104.05929, 2021. [View Article](https://doi.org/10.1016/j.eswa.2022.116903)


## Program Developers

* **Haoyuan Sun** - [@wuxishy](https://github.com/wuxishy) developed the *cfr-v1* of the program
* **Yoojin Chung** - [@mangoQueen](https://github.com/mangoQueen) took part in the *cfr-v2 program development*
* **Mohammad Nazmul Haque** - [@MohammadNHaque](https://github.com/MohammadNHaque) took part in developing the *cfr-v2 and onwards programs*
* **Andrew Ciezak** - [@andrew-ciezak](https://github.com/andrew-ciezak) took part in developing the *cfr-v2.5 and onwards program*
* **Kevin Huang** - [@KevinHuang8](https://github.com/KevinHuang8) took part in developing the linear regression usage in the *cfr-v2.5* program

## Contributors
* **Pablo Moscato**: Conceptualization, Methodology, Formal analysis, Investigation, Supervision, Project administration.
* **Haoyun Sun**: Methodology, Software, Validation, Formal analysis, Investigation.
* **Mohammad Nazmul Haque**: Methodology, Software, Validation, Formal analysis, Investigation, Visualization.
