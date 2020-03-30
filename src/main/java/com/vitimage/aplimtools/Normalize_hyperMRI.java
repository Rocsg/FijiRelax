package com.vitimage.aplimtools;
import com.vitimage.common.VitiDialogs;
import com.vitimage.mrutils.HyperMRIT1T2;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.frame.PlugInFrame;

/** Basic example of how to add a function to the menu "AplimTools" in ImageJ > Plugins > Analyze
 * To start a new plugin, right click on the directory com.vitimage.aplimtools, and select "New class". Just choose a name, with an underscore in it (another name than Aplim_Example)
 * When the new code opens, use Aplim_Example to start fast. You can copy these ten lines of code, and changing "Aplim_Example" to "The_name_you_choosed" 
 * Once done, you can run your code using "Run" (CTRL + F11).
 * When  your plugin is ok (and tested), add an entry in src/main/resources/plugins.config to be sure your plugin appears in ImageJ/Fiji. 
 * Then, actualize the code to the last version published (right-click on the project (AplimTools) > Team > Pull)
 * Then push your modifications to the online repository  (right-click on the project (AplimTools) > Team > Commit > Commit and push. 
 * Your plugin will be inspected, then deployed in the next release version
 * Enjoy !

 * @author fernandr
 *
 */


public class Normalize_hyperMRI  extends PlugInFrame{

	private static final long serialVersionUID = 1L;

	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		Normalize_hyperMRI ex=new Normalize_hyperMRI();
		ex.run("");
	}
	
	public Normalize_hyperMRI() {
		super("Normalize hyperMRI");
	}

	public void run(String arg) {
		IJ.log("Starting hypernormalisation");
		HyperMRIT1T2 hyp=new HyperMRIT1T2(VitiDialogs.chooseOneImageUI("Select image","Select your hyperimage to be processed"));
		hyp.showCopy("Initial Hyper MRI");
		hyp.normalizeMRSignal();
		hyp.showCopy("Normalized Hyper MRI");
	}
}
