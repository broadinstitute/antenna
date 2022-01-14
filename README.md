# Antenna sgRNA-seq pipeline
Antenna is a cloud pipeline for the identification of sgRNA reads in viral next generation sequencing data.

## System Requirements
The antenna pipeline can be run using any WDL execution engine, such as Cromwell or miniwdl. The pipeline was been tested using miniwdl version 1.4.3. The miniwdl package can be found [here](https://github.com/chanzuckerberg/miniwdl).

To install the python package you will require python 3.8 or later and the required packages listed under REQUIREMENTS.txt

## Installation Guide
Running the antenna pipeline does not require any installation if a suitable WDL run engine is available. In order to run the antenna pipeline export the package from github and run the `wdl/antanna.wdl` file.

```
git clone https://github.com/broadinstitute/antenna.git
miniwdl run --verbose -i input.json ../../wdl/antenna.wdl
```

If you wish to install the antenna python package you can do so using the following commands.
```
git clone https://github.com/broadinstitute/antenna.git
pip install -e antenna/
```

## Demo
An example test script for the pipeline is provided under `testing/testwdl/run.sh`. This script will execute the WDL pipeline using inputs in the input.json file.

To run this example, please edit the file `testing/testwdl/input.json` and replace [LOCATION] with the download location of the repository in your system.

## Instructions for Use
To run the pipeline on your data call the wdl with a custom `input.json` file. For example with miniwdl you would execute:
```
miniwdl run --verbose -i input.json ../../wdl/antenna.wdl
```

# Generating the Docker Image
The docker build script uses the latest main version in github for installation. After updating github, run:
```
cd docker
./build [version]
docker push quay.io/nbarkas_1/antenna:[version]
```

You will need to create a custom quya.io repository if you wish to customize the image.
