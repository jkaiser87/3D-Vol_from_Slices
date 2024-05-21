**About this Pipeline**
Using FIJI and MATLAB, this pipeline can be used to track injection volumes and translate the location of the injection site into a 3D model.

**Requirements**
- FIJI.app (downloaded from original website)
  - download Fiji folder from here and paste into Fiji.app folder
- folder with single slice tif files of 1 coronal brain, sorted rostral to caudal (!highly recommended to end with padding of slice numbers, *_s001 etc, as Fiji sorts 1, 10, 2 ...)
- MATLAB
	- Go through AP_histology installation as described on their website -	https://github.com/petersaj/AP_histology 
  - install the Curve Fitting Toolbox
	- install the natsortfile add-on (Natural-Order Filename Sort Version 3.4.5 by Stephen23)
	- Download the MATLAB folder from here, put Documents folder into your userfolder “Documents/MATLAB/” (additionally to the AP-histology required files)

**Running the pipeline**
**1. FIJI part**
1.1. Setup
	- Click the setup folder, select the folder that contains the single slice tifs (should not contain any other brain images)

1.2. Pre-process
	- To make the next step a bit easier, go through the pre-processing to rotate the slice upright and potentially flip it to the same side (if necessary)

1.3. Volume tracing
	- Select channel that you want to process (the one that has your signal)
 	- FIJI will automatically isolate that channel. 
		- Click continue if there is no signal
		- If there is signal available in the section, outline it using the pre-selected free selection tool

**2. MATLAB part**
2.1 Once the whole animal has been processed, go into MATLAB and open the code file "AP_1_VOL_SingleAnimal_20240521". 

Adapt setting at beginning as necessary
- rerun_histology should be 1 if you run it the first time.
	- This will guide you to AP_histology, where you can set the corresponding allen atlas section as well as overlap the structures onto your own section
- Choose channels to process
- Choose colors to add
- If you want to additionally highlight certain brain structures from the Atlas, put their name in structure_names. Alternatively comment this line and uncomment the next to make an empty vector

Press run, which should create a plot for you

2.2. The second file in the MATLAB folder ("AP_2_VOL_PlotAllAnimalsInFolder_20240502") can be used to combine several animals into one 3D model.
