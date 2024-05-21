macro "Pre-Process Slices Tool - D16D17D23D26D27D2aD32D33D34D35D36D37D38D39D3aD3bD43D44D45D46D47D48D49D4aD53D54D55D58D59D5aD61D62D63D64D69D6aD6bD6cD71D72D73D74D79D7aD7bD7cD83D84D85D88D89D8aD93D94D95D96D97D98D99D9aDa2Da3Da4Da5Da6Da7Da8Da9DaaDabDb3Db6Db7DbaDc6Dc7CcccD22D24D29D2bD42D4bD56D57D65D68D75D78D86D87D92D9bDb2Db4Db9Dbb"{

var input = call("ij.Prefs.get", "input.x",0);
output = input + "Processed" + File.separator; 

File.makeDirectory(output);

if (!File.exists(output))
exit("Unable to create directory "+output);


//var docoords = call("ij.Prefs.get", "docoords.x", 0);  // Default is "no"

print("\nPre-Processing of Images in folder "+input);

list = getFileList(input);
suffix = ".tif";

//opens all tools
run("Channels Tool...");
run("Brightness/Contrast...");
run("ROI Manager...");

run("Close All");

// Process each subfolder under the input directory
if (docoords == 1) {
createCoordDirectories(subFolderPath);
}

processFolder(input, docoords);

showMessage("Done", "--- All slices pre-processed. \n---- You can now continue to the next step of volume tracing.");
}

}


function processFolder(input, docoords) {
    var list = getFileList(input);
    list = Array.sort(list);
    for (var k = 0; k < list.length; k++) {
        if (endsWith(list[k], ".tif")) {
            processFile(input, docoords, list[k]);
        }
    }
}


function processFile(input, docoords, file) {
   
print(" - Processing file " + file);
open(input + file);

run("Enhance Contrast", "saturated=0.35");

// only turn if higher than wide (coronal slices assumed)
getDimensions(width, height, channels, slices, frames);
if(height>width) {
	run("Rotate 90 Degrees Right");	
}

//rotate to midline
setTool("line");
waitForUser("Action required", "Make a line through the midline \n START AT TOP and go to bottom! \nIf you want to skip image: delete selection (Ctrl+Shift+A)");

if(selectionType()!= -1) {
if (selectionType!=5)
      exit("Straight line selection required");
      getLine(x1, y1, x2, y2, lineWidth);

      angle = getAngle(x1, y1, x2, y2);
	  RotateAngle = angle + 90;
		run("Arbitrarily...", "angle="+RotateAngle+" interpolate");

//save image to save rotation
	
	print(" --- Saving in folder "+output);
	run("Select None");
	saveAs("Tiff",output+file);

if (docoords == 1) {
run("Enhance Contrast", "saturated=0.35");

// get coordinates for midline and outside left and right
ROIs = roiManager("count");
if (ROIs >0) {
roiManager("Deselect");
roiManager("Delete");
}

	//set and rename points
	setTool("point");
	waitForUser("Action required", "Select the left-most point of the left cortex");
	roiManager("Add");
	last = roiManager("count");
	roiManager("select", last-1);
	roiManager("Rename", "LEFT CTX");
	waitForUser("Action required", "Select the Midline right where the tissue meets on the top");
	roiManager("Add");
	last = roiManager("count");
	roiManager("select", last-1);
	roiManager("Rename", "MIDLINE");
	waitForUser("Action required", "Select the right-most point of the right cortex");
	roiManager("Add");
	last = roiManager("count");
	roiManager("select", last-1);
	roiManager("Rename", "RIGHT CTX");
	
	
	saveCoordData(input, file);
	  	
}
}
run("Close All");
}

//---------- FUNCTIONS

function createCoordDirectories(subFolderPath) {
    var subfolders = newArray("ZIP", "ZIP" + File.separator + "COORD", "CSV","CSV" + File.separator + "COORD");
	for (ff = 0; ff < subfolders.length; ff++) {
			animalsub = subFolderPath + File.separator + subfolders[ff];
			File.makeDirectory(animalsub);
 			if (!File.exists(animalsub))
 			 exit("Unable to create directory");
		}
}

function saveCoordData(input, file) {
    
    //save coord roi as zip
	var zipPath = input + "ZIP\\COORD" + File.separator + replace(file, ".tif", ".zip");
	roiManager("save", zipPath);
	//save as csv
	ROIs = roiManager("count");
	array=newArray();
	for(f=0;f<ROIs;f++){
	array[f]=f;
	}

	//remove scale otherwise particle analysis returns in inches
	run("Set Scale...", "distance=0 known=0 unit=pixel");
	
	roiManager("select", array);
	roiManager("Measure");
	var csvPath = input + "CSV\\COORD"+File.separator+replace(file,".tif",".csv"); 
	saveAs("Measurements", csvPath);
	run("Clear Results");
	print(" --- Coordinates saved as ZIP and CSV files");
}

// Returns the angle in degrees between the specified line and the horizontal axis.
  function getAngle(x1, y1, x2, y2) {
      q1=0; q2orq3=2; q4=3; //quadrant
      dx = x2-x1;
      dy = y1-y2;
      if (dx!=0)
          angle = atan(dy/dx);
      else {
          if (dy>=0)
              angle = PI/2;
          else
              angle = -PI/2;
      }
      angle = (180/PI)*angle;
      if (dx>=0 && dy>=0)
           quadrant = q1;
      else if (dx<0)
          quadrant = q2orq3;
      else
          quadrant = q4;
      if (quadrant==q2orq3)
          angle = angle+180.0;
      else if (quadrant==q4)
          angle = angle+360.0;
      return angle;
      
  }