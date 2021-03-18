## Getting start
The following instructions explain how to compile **Atif** on your local machine. 

[An example input file](Atif/example/input_chain.dat) are also provided to help you get started with your own project using Atif. The meanings of the parameters in the [input file](Atif/example/input_chain.dat) have also been explained clearly in it. (If you get confused about the input file, please directly contact [the developer](https://github.com/jiangj-physchem)) 

Note that this guide is written for the Linux working environment, but it should also apply to the Mac OS terminal and Xcode and Visual Studio C++.

### Prerequisites

You should make sure that `cmake` and `gcc` are installed properly on your machine.

### Installation

Download the entire [Atif](Atif) directory to your local path `/Users/jiangj/`.
```
cd /Users/jiangj/Atif/cmake/
make
```
Then you will get an executable file named **ATIFexe** in the folder `/Users/jiangj/Atif/cmake/`. Now you are ready to start your project!
```
nohup ./ATIFexe < /Users/jiangj/Atif/example/input_chain.dat &
```
You will find the output files in the folder as shown in the last line of the [input file](Atif//example/input_chain.dat).
