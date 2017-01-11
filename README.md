# origamiDesignVisualization
Tools for visualization of DNA origami designed using cadnano

To run this code, you must have a copy of cadnano2, and add its location to the path in analysisHelpers.py (sys.path.append).

The base maps use a savistky golay filter, and there may be problems with numpy and nan values. If this occurs, create a virtual environment with python and uninstall MKL. Then reinstall needed packages (https://docs.continuum.io/mkl-optimizations/)
