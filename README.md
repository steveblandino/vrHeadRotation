# Head Rotation Model for AR/VR System Level Simulations
> A MATLAB open-source model to generate 3-DOF head rotation.


## Table of Contents
* [Installation](#installation)
* [Requirements](#requirements)
* [How to Run](#how-to-run)
* [Contributing](#contributing)
* [Authors](#authors)
* [License](#license)

## Installation
The software does not require any installation procedure: simply download or clone the repository to your local folders.

### Requirements
The codebase is written in MATLAB. It is currently being tested on MATLAB R2021a.
No toolbox are needed.
To reproduce the results obtained in the paper *Steve Blandino, Tanguy Ropitault, Raied Caromi, Jacob Chakareski, Mahmudur Khan, 
  and Nada Golmie. 2021. Head Rotation Model for VR System Level Simulations*, please download the NJIT 6DOF VR Navigation Dataset at [https://www.jakov.org](https://www.jakov.org).

## How to Run 
* To run the model simply run the script `src/mainRotationModel.m`.
* To reproduce the results in the paper *Steve Blandino, Tanguy Ropitault, Raied Caromi, Jacob Chakareski, Mahmudur Khan, 
  and Nada Golmie. 2021. Head Rotation Model for VR System Level Simulations*
	* Download the NJIT 6DOF VR Navigation Dataset at [https://www.jakov.org](https://www.jakov.org)
	* Extract the contentent of Traces_6DOF_NJIT.zip into the folder ./Traces_6DOF_NJIT
	* Run the function `src/mobilityDataProcessing.m`

## Contributing
Feedbacks and additions are more than welcomed! You can directly contact the [authors](#Authors) for any information.


## Authors

| NIST | 

The Q-D Realization software has been developed at NIST by Steve Blandino (steve.blandino@nist.gov)


## License
Please refer to the [LICENSE](LICENSE) file for more information.
