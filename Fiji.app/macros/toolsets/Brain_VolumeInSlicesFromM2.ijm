// Cell Count Analysis tools

macro "Unused Tool-1 - " {}  // leave empty slot

macro "Setup Folder Action Tool - icon:folder.png" {
  
  var defaultPath = call("ij.Prefs.get", "input.x",0); //retrieve last input
  
  Dialog.create("Setup");
  Dialog.addMessage("Choose folder to process \n"
+"(folder in which the TIF files are)\n");
  Dialog.addDirectory("\n", defaultPath)
  Dialog.addMessage("------------------------");
  //Dialog.addCheckbox("Get coordinates? (default false)", false); //this can be useful if you want to make calculations on the width of the coronal slices
  Dialog.show;
  
  var input = Dialog.getString();
  if (!endsWith(input, File.separator)) {
      input += File.separator;
	}
//  var docoords = Dialog.getCheckbox();

	call("ij.Prefs.set", "input.x", input);
	//call("ij.Prefs.set", "docoords.x",docoords );
	
	print("---- Processing folder "+input);
	
} 

macro "Pre-process Slices Action Tool - icon:imageproc.png"{
 runMacro(getDirectory("macros")+"//toolsets//scripts//M2_VOL_SetupSlices.ijm");
 //could change this to NOT select left/right/midline, but rotation is good
}

macro "Draw Volume for All Slices Action Tool - icon:shape.png"{
 runMacro(getDirectory("macros")+"//toolsets//scripts//M2_VOL_channeldep.ijm");
}

macro "Flip image [F]"{
	run("Flip Horizontally");
    print ("--- Selection flipped.");
}