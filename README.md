# MOO-SDNBI

# MATLAB-GAMS Optimization Project

This repository contains a MATLAB-based optimization project that interfaces with GAMS for solving optimization problems and uses QConvex for convex hull calculations. The project utilizes the GAMS-MATLAB package for interfacing between MATLAB and GAMS.

## Prerequisites

Before you can run the project, ensure you have the following software installed:

1. **MATLAB**: R2018a or later is recommended.
2. **GAMS**: The General Algebraic Modeling System. Ensure it's installed and properly configured in your system's PATH.
3. **QConvex**: Download the QConvex executable from [QHull](http://www.qhull.org/html/qconvex.htm). (Note that qconvex.exe file is already uploaded here)
4. **GAMS-MATLAB Interface**: Follow the installation instructions from the [GAMS-MATLAB Documentation](https://gams-matlab.readthedocs.io/en/latest/).

## Repository Structure

- `CSx_../`: Directory containing all the necessary MATLAB files for each benchmark problems.
- `main_SDNBI_2obj_RF.m`: The main MATLAB script to start the SDNBI optimization process.
- `CSx_../functions/`: Directory containing all the necessary MATLAB function files.
- `CSx_../GAMS.m`: The GAMS-MATLAB interface class.
- `CSx_../qconvex.exe`: QConvex executable for convex hull calculations.
- `CSx_../data/`: Directory containing example data files (e.g., `PPoints_CS2_oNBI.mat`, `PF_CS2_oNBI.mat`, etc.).
- `README.md`: This readme file.

## Setup Instructions

1. **Clone the Repository**:
    ```bash
    git clone https://github.com/yourusername/yourrepository.git
    cd yourrepository
    ```

2. **Add Paths in MATLAB**:
    Open MATLAB and add the necessary paths to the MATLAB environment:
    ```matlab
    addpath(genpath('functions'));
    addpath(genpath('data'));
    ```

3. **Install GAMS-MATLAB Interface**:
    Follow the [installation guide](https://gams-matlab.readthedocs.io/en/latest/) to set up the GAMS-MATLAB interface. Ensure `GAMS.m` is accessible in your MATLAB path.

4. **Download and Setup QConvex**:
    - Download `qconvex.exe` from [QHull](http://www.qhull.org/html/qconvex.htm).
    - Place `qconvex.exe` in the CSxx directory of this project. 

5. **Verify Data Files**:
    Ensure all required data files (`PPoints_CS2_oNBI.mat`, `PF_CS2_oNBI.mat`, etc.) are available in the `data` directory.

## Running the Project

1. **Run the Main Script**:
    In MATLAB, navigate to the project directory and run the main script:
    ```matlab
    run('main_....m');
    ```

## Main Functions

### main.m

This is the entry point for the project. It initializes parameters, sets up the optimization problem, and iterates through the solution process.

### functions/

- **f_MtoG_WS_2obj.m**: Function to handle specific GAMS-related operations.
- **f_outer_test.m**: Function to perform outer approximation tests.
- **generate_dummy.m**: Helper function to generate dummy points for convex hull calculations.
- **convexHull_R1.m**: Function to compute the convex hull using QConvex.
- **f_clusterParetoPoints.m**: Function to cluster Pareto points and manage subregions.

## Troubleshooting

- **GAMS Errors**: Ensure GAMS is properly installed and the environment variables are correctly set.
- **MATLAB Path Issues**: Verify that all necessary directories are added to the MATLAB path.
- **QConvex Issues**: Ensure `qconvex.exe` is accessible in the system's PATH or located in the project directory.

## Contributing

Feel free to fork this repository and submit pull requests. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
