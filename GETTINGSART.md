## Getting start
The following instructions explain how to compile **Atif** on your local machine. 

[An example input file]('./Atif/input/input.dat') are also provided to help you get started with your own project using Atif. The meanings of the parameters in the [input file]('./Atif/input/input.dat') have also been explained clearly in it. (If you get confused about the input file, please directly contact [the developer](https://github.com/jiangj-physchem)) 

Note that this guide is written for the Linux working environment, but it should also apply to the Mac OS terminal and Xcode and Visual Studio C++.

### Prerequisites

You should make sure that `cmake` and `gcc` are installed properly on your machine.

### Installation

Download the entire [Atif](./Atif) directory to your local path `your_path_to_atif/`.
```
cd your_path_to_atif/Atif/execute/
make
```
Then you will get an executable file named **PolymerDFT** in the folder `your_path_to_atif/Atif/execute/`. You need to copy this executable file into the same folder with the **input** file. Then you are ready to start your project!
```
cp PolymerDFT your_path_to_atif/Atif/input/
cd your_path_to_atif/Atif/input/
nohup ./PolymerDFT &
```
You will find the output files in the folder as shown in the last line of the [input file]('./Atif/input/input.dat'). If this folder doesn't exit, please create it:
```
cd your_path_to_atif/Atif/
mkdir output
```
