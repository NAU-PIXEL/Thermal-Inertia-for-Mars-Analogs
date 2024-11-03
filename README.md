# Thermal Inertia for Mars Analogs (TIMA) Model
 TIMA is a finite element numerical model of soil heat transfer that that considers time-variability in thermal conductivity as a function of temperature, moisture, and composition. The detailed parametrization of thermal conductivity is meant to enable direct translation between soil thermophysical properties on Earth and Mars. This repository is maintained by [Ari Koeppel](https://earthsciences.dartmouth.edu/people/ari-koeppel). Please cite Koeppel et al., (2024) when using the model.

# References:
Koeppel, A.H., Edwards, C.S., Edgar, L.A., Nowicki, S., Bennett, K.A., Gullikson, A., Piqueux, S., Eifert, H., Chapline, D. and Rogers, A.D., 2024. A novel surface energy balance method for thermal inertia studies of terrestrial analogs. Earth and Space Science, 11(9), p.e2023EA003259.

# Overview: 
The original implementation of the model was to derive a characteristic soil thermal conductivity under dry conditions at 300K. The procedure achieves this by inputing micrometeorological data into a surface energy balance forward model and adjusting thermal conductivity (along with 5 other modifying parameters) to fit surface skin temperature data (typically obtained through radiometer observation).

# Visualize:

# Input/Output Setup:

## Input data:
Example input files are located in the [`in`](hhttps://github.com/NAU-PIXEL/Thermal-Inertia-for-Mars-Analogs/tree/main/Example%20Data) folder.
