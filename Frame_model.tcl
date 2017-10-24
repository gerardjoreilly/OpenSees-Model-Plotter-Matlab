# --------------------------------------
# 3 Storey RC Frame
# --------------------------------------
# Copyright by Gerard J. O'Reilly, 2017
# --------------------------------------
# Gerard J. O'Reilly
# IUSS Pavia
# All units are in kN, m
# --------------------------------------
wipe;
# --------------------------------------
# This model is a simple elastic frame that is used as an example frame with which to use the model_plot Matlab
# tool provided. It is a 3D model and the idea is to be able to plot the model and a number of the mode shapes
# in Matlab. The aim of this plotter is very simple: to be able to see what you are modelling relatively quickly
# before performing and analysis.
# --------------------------------------

# --------------------------------------
# Define some basics
# --------------------------------------
set H	3.; 	# Storey height
set B	5.; 	# Bay width
set M	60.; 	# Mass per node

# Material
set E	20.0e6;
set G	[expr 0.4*$E];

# Column
set hc	0.5;
set Ac	[expr $hc*$hc];
set Iyc	[expr pow($hc,4)/12];
set Izc	[expr pow($hc,4)/12];
set Jc	[expr 2.25*pow($hc,4)];

# Beams
set hb	0.5;
set bb	0.3;
set Ab	[expr $bb*$hb];
set Iyb	[expr $hb*pow($bb,3)/12];
set Izb	[expr $bb*pow($hb,3)/12];
set Jb	[expr 0.2*$hb*pow($bb,3)];

puts [expr $E*$Ac]
puts [expr $E*$Iyc]
puts [expr $E*$Izc]
puts [expr $G*$Jc]

puts [expr $E*$Ab]
puts [expr $E*$Iyb]
puts [expr $E*$Izb]
puts [expr $G*$Jb]

# --------------------------------------
# Define the model
# --------------------------------------
model BasicBuilder -ndm 3 -ndf 6

# --------------------------------------
# Define nodes and fixity
# --------------------------------------
# Note that each of the nodes defined need to have a mass value assigned also

# Ground Level
# node	tag		X		Y		Z		-mass	mX		mY		mZ 		IX 		IY 		IZ
node		110		0.0		0.0		0.0		-mass	0.0		0.0		0.0		0.0		0.0		0.0
node		210		$B		0.0		0.0		-mass	0.0		0.0		0.0		0.0		0.0		0.0
node		310		[expr 2*$B]	0.0		0.0		-mass	0.0		0.0		0.0		0.0		0.0		0.0

node		120		0.0		$B		0.0	-mass	0.0		0.0		0.0		0.0		0.0		0.0
node		220		$B		$B		0.0	-mass	0.0		0.0		0.0		0.0		0.0		0.0
node		320		[expr 2*$B]	$B		0.0	-mass	0.0		0.0		0.0		0.0		0.0		0.0

node		130		0.0		[expr 2*$B]	0.0	-mass	0.0		0.0		0.0		0.0		0.0		0.0
node		230		$B		[expr 2*$B]	0.0	-mass	0.0		0.0		0.0		0.0		0.0		0.0
node		330		[expr 2*$B]	[expr 2*$B]	0.0	-mass	0.0		0.0		0.0		0.0		0.0		0.0

# First Level
# node	tag		X		Y		Z	-mass	mX		mY		mZ 		IX 		IY 		IZ
node		111		0.0		0.0		$H	-mass	$M		$M		0.0		0.0		0.0		0.0
node		211		$B		0.0		$H	-mass	$M		$M		0.0		0.0		0.0		0.0
node		311		[expr 2*$B]	0.0		$H	-mass	$M		$M		0.0		0.0		0.0		0.0

node		121		0.0		$B		$H	-mass	$M		$M		0.0		0.0		0.0		0.0
node		221		$B		$B		$H	-mass	$M		$M		0.0		0.0		0.0		0.0
node		321		[expr 2*$B]	$B		$H	-mass	$M		$M		0.0		0.0		0.0		0.0

node		131		0.0		[expr 2*$B]	$H	-mass	$M		$M		0.0		0.0		0.0		0.0
node		231		$B		[expr 2*$B]	$H	-mass	$M		$M		0.0		0.0		0.0		0.0
node		331		[expr 2*$B]	[expr 2*$B]	$H	-mass	$M		$M		0.0		0.0		0.0		0.0


# Second Level
# node	tag		X		Y		Z			-mass	mX		mY		mZ 		IX 		IY 		IZ
node		112		0.0		0.0		[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		212		$B		0.0		[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		312		[expr 2*$B]	0.0		[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0

node		122		0.0		$B		[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		222		$B		$B		[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		322		[expr 2*$B]	$B		[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0

node		132		0.0		[expr 2*$B]	[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		232		$B		[expr 2*$B]	[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		332		[expr 2*$B]	[expr 2*$B]	[expr 2*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0

# Third Level
# node	tag		X		Y		Z			-mass	mX		mY		mZ 		IX 		IY 		IZ
node		113		0.0		0.0		[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		213		$B		0.0		[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		313		[expr 2*$B]	0.0		[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0

node		123		0.0		$B		[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		223		$B		$B		[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		323		[expr 2*$B]	$B		[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0

node		133		0.0		[expr 2*$B]	[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		233		$B		[expr 2*$B]	[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0
node		333		[expr 2*$B]	[expr 2*$B]	[expr 3*$H]		-mass	$M		$M		0.0		0.0		0.0		0.0

# --------------------------------------
# Define the fixity
# --------------------------------------
# fix 	node	dX 	dY 	dZ 	rX 	rY 	rZ
fix		110		1	1	1	1	1	1
fix		210		1	1	1	1	1	1
fix		310		1	1	1	1	1	1
fix		120		1	1	1	1	1	1
fix		220		1	1	1	1	1	1
fix		320		1	1	1	1	1	1
fix		130		1	1	1	1	1	1
fix		230		1	1	1	1	1	1
fix		330		1	1	1	1	1	1

rigidDiaphragm 3	111 211 311 121 321 131 231 331
rigidDiaphragm 3	112 212 312 122 322 132 232 332
rigidDiaphragm 3	113 213 313 123 323 133 233 333

# --------------------------------------
# Geometric Transformation
# --------------------------------------
geomTransf Linear 	1  0 1 0; # z is in Y
geomTransf Linear 	2  -2 0 0; # z is in -X

# --------------------------------------
# Define the Columns
# --------------------------------------
# Elastic elements
# First Storey
# element elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 	$E 	$G  $J 	$Iy  $Iz 	$transfTag <-mass $massDens> 		<-cMass>
element elasticBeamColumn 	5111 		110 	111		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5211 		210 	211		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5311 		310 	311		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

element elasticBeamColumn 	5121 		120 	121		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5221 		220 	221		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5321 		320 	321		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

element elasticBeamColumn 	5131 		130 	131		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5231 		230 	231		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5331 		330 	331		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

# Second Storey
# element elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 	$E 	$G  $J 	$Iy  $Iz 	$transfTag <-mass $massDens> 		<-cMass>
element elasticBeamColumn 	5112 		111 	112		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5212 		211 	212		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5312 		311 	312		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

element elasticBeamColumn 	5122 		121 	122		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5222 		221 	222		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5322 		321 	322		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

element elasticBeamColumn 	5132 		131 	132		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5232 		231 	232		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5332 		331 	332		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

# Third Storey
# element elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 	$E 	$G  $J 	$Iy  $Iz 	$transfTag <-mass $massDens> 		<-cMass>
element elasticBeamColumn 	5113 		112 	113		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5213 		212 	213		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5313 		312 	313		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

element elasticBeamColumn 	5123 		122 	123		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5223 		222 	223		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5323 		322 	323		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

element elasticBeamColumn 	5133 		132 	133		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5233 		232 	233		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1
element elasticBeamColumn 	5333 		332 	333		$Ac $E 	$G 	$Jc $Iyc $Iyc 	1

# --------------------------------------
# Define the Beams
# --------------------------------------
# Elastic elements
# First Level
# element elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 	$E 	$G  $J 	$Iy  $Iz 	$transfTag <-mass $massDens> 		<-cMass>
element elasticBeamColumn 	7111 		111 	211		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7211 		211 	311		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7121 		121 	221		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7221 		221 	321		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7131 		131 	231		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7231 		231 	331		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1

element elasticBeamColumn 	8111 		111 	121		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8121 		121 	131		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8211 		211 	221		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8221 		221 	231		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8311 		311 	321		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8321 		321 	331		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2

# Second Level
# element elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 	$E 	$G  $J 	$Iy  $Iz 	$transfTag <-mass $massDens> 		<-cMass>
element elasticBeamColumn 	7112 		112 	212		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7212 		212 	312		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7122 		122 	222		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7222 		222 	322		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7132 		132 	232		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7232 		232 	332		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1

element elasticBeamColumn 	8112 		112 	122		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8122 		122 	132		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8212 		212 	222		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8222 		222 	232		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8312 		312 	322		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8322 		322 	332		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2

# Third Level
# element elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 	$E 	$G  $J 	$Iy  $Iz 	$transfTag <-mass $massDens> 		<-cMass>
element elasticBeamColumn 	7113 		113 	213		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7213 		213 	313		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7123 		123 	223		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7223 		223 	323		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7133 		133 	233		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1
element elasticBeamColumn 	7233 		233 	333		$Ab $E 	$G 	$Jb $Iyb $Iyb 	1

element elasticBeamColumn 	8113 		113 	123		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8123 		123 	133		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8213 		213 	223		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8223 		223 	233		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8313 		313 	323		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2
element elasticBeamColumn 	8323 		323 	333		$Ab $E 	$G 	$Jb $Iyb $Iyb 	2

# --------------------------------------
# Modal Analysis
# --------------------------------------
source modalAnalysis.tcl
modalAnalysis 7 1 "_modal"

# --------------------------------------
# Print the model
# --------------------------------------
# Delete the old one first
file delete model.txt
# Print the new one
print model.txt


# Finish
wipe;
