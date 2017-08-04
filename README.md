# Matlab-Based Plotter for OpenSees Models

This repository contains a Matlab-based plotter for models developed in OpenSees. 

The only file the user needs to look at is “model_plot.m”, all the others are just functions.

It is not a GUI, but rather a simple tool that allows the user to plot and examine the model being developed with a single command from Matlab. It works for models developed in either 2D or 3D, although it has been mainly used and tested for 3D models. It is capable of plotting mode shapes also that are generated using the “modalAnalysis.tcl” procedure provided. 

There is an example provided with the “Frame_model.tcl”, where by running this file, the necessary input to use this plotter will be provided as a first example of how to implement it.

If you find any bugs or have any requested add in, let me know.

## Using the Plotter

1. Build your model. This does not have to be a complete model that can be analysed; a bunch of floating nodes is fine, for example. This way you can monitor the gradual progression of your model as you build. 
2. Once you have built your model, use the “print” command in OpenSees to print the entire model to a text file. In the example provided, this is at the end of the file where inserting the command “print model.txt” and then executing the file “Frame_model.tcl” in OpenSees produces the required file for this plotter.
3. Once you have your model.txt file in your directory, ensure that you have all of the Matlab files provided here on your current Matlab directory. 
4. Now you can use the “model_plot” function. The required input are described in the file but if you just type “help model_plot” in the Matlab command window it will print them. 
5. Once you execute the “model_plot” function with the relevant inputs in Matlab, it will ask you to select the “model.txt” file that OpenSees produced. Load that and the plotter will plot a 3D figure of your model to examine. Use the camera tools available in Matlab (View>Camera Toolbar) to rotate the model etc.

Additionally, you can use this plotter to examine the mode shapes of your building. 
1. You need to perform a modal analysis of your structure using the “modalAnalysis.tcl” TCL procedure provided. Using this procedure within you model will return the modal values on screen in addition to producing two files: a file containing the eigenvectors and another listing the periods of vibration.
2. Requesting that the modes of vibration be plotted will lead to the plotter also asking you for the eigenvectors and periods file. 
3. Select these and the plotter will plot separate figures for the model and then subsequent modes of vibration with the mode number and period of vibration noted at the top. 

## Available Elements:
The following elements will be plotted:
* ForceBeamColumn (2D and 3D)
* ElasticBeamColumn (2D and 3D)
* ZeroLength
* Truss
* CorotTrussSection

If you use any element that is not listed above, the plotter will not plot it.

If a single node has a mass, it will be plotted in red whereas all other nodes are plotted black.

## Things to Note and Future Developments
Load patterns are not plotted.

Constraints are not plotted.

When plotting mode shapes, element element deformations (i.e. section curvatures) are not plotted.

You can also use this to plot the final state of a model after a pushover or dynamic analysis has been performed. Just put a print command at the point in the TCL model where the analysis is finished and procedure as if it were the “model.txt” file described above.





