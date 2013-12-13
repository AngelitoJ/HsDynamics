HsDynamics

This package contains a tool to perform Semiclassical Molecular Dynamics (NVT)

Given some initial conditions, the usual MD loop uses, some external programs (currently Molcas and Gaussian are supported) capable of calculating the energy gradient, parsings those tools output logs and perform the whole "Semi- Classical part" of the Molecular Dynamics.

Current modules feature these MD components:

	- PES evalutators:
		- abinitio pointwise evalution (via external 3d party packages like Molcas, Gaussian etc..).
		- Fixed Cuadratic interpolation around some point of interest.
		- Dynamic weighted multiple points interpolation (a.k.a Palmeiro interpolator) to more accurate surface interpolation.
	- Symplectic Integrator:
		- Velocity Verlet algorithm.
		- Nosé Hoover thermostate chains for a constant temperature bath.
	- Energy corrections and diverse phenomena:
		- Tully-Hammes-Schiffer hopping algorithm (along with correction of Persico-Granucci).

It also features the possibility to add external forces to the molecule, to simulate constrained conditions that can be found, in example, in a protein binding pocket.

The package exercise multiple haskell facilities, featuring several parsec submodules and making use of Async
in order to manage external programs.

	Thanks: 
		Raul Palmeiro: Fortran coded multinodal interpolator


@2013 Felipe Zapata, Angel Alvarez, Alessio Valentini. The RESMOL group
