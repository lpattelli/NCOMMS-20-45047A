# NCOMMS-20-45047A

MATLAB scripts used to generate the numerical data for the paper:

Marco Leonetti, Lorenzo Pattelli, Simone De Panfilis, Diederik S. Wiersma, Giancarlo Ruocco. Spatial coherence of light inside three dimensional media. *Nature Communications* (2021).

Simulations performed using [CELES v2.2](https://github.com/disordered-photonics/celes/releases/tag/v2.2)

### Contents

* **spherical_target_120k.mat**: Database containing the coordinates and radii of the scattering particles and of the nanorulers.
* **CELES_gattaquant.m**: Script launching multiple CELES simulations for different randomly oriented plane-wave illumination conditions. For each simulation, two output files will be saved with the field measured along the nanorulers and on a plane cutting through the spherical target of particles.
* **recombine_simulations.m**: Script reading the output files saved by CELES_gattaquant.m and recombining them with random phases. The distribution of integrated intensities at the two ends of each nanoruler is then used to retrieve the number of independent degrees of freedom of the field.
* **plot_intensity.m**: Convenience function plotting an intensity pattern calculated on a plane.
* **CELES_twoparticles.m**: Script launching multiple CELES simulations for different randomly oriented plane-wave illumination conditions for a minimal system comprising only two particles. Field values in the gap between the two particles are computed for several random dephasing combinations to building a distribution of their corresponding propagation direction.
* **cohrencymatrix3D.m**: Function calculating the 3D coherency matrix given a complex-valued electric field vector.
* **field_pol.m**: Function calculating the parameters of the polarization ellipse defined by a complex-valued electric field vector. The propagation direction is perpendicular to the plane of the polarization ellipse.
