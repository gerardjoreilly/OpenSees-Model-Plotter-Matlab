# Matlab Plotter for OpenSees Models

This repository contains a Matlab-based plotter for models developed in OpenSees. It is not a GUI, but rather a simple tool that allows the user to plot and examine the model being developed with a single command from Matlab. It works for models developed in either 2D or 3D, although it has been mainly used and tested for 3D models. It is capable of plotting mode shapes also that are generated using the “modalAnalysis.tcl” procedure provided. 

If you find any bugs or have any requested add in, let me know.

## Using the Plotter

1. Build your model. This does not have to be a complete model that can be analysed so you can monitor the gradual progression of your model as you build. 
2. Once you have built you model, use the “print” command to print the entire model to a text file. In the example provided, this is at the end of the file where inserting the command “print model.txt” and then executing the model in OpenSees produces the required file for this plotter.
3. Once you have your model.txt file in your directory, ensure that you have all of the Matlab files provided here on your current Matlab directory. 
4. Now you can use the “model_plot” function. The required input are described in the file but if you just type “help model_plot” of Matlab it will tell you the required inputs.


## Notes:
The following elements will be plotted:
	