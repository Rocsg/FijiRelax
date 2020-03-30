package com.vitimage.aplimtools;
import ij.IJ;
import ij.ImageJ;
import ij.plugin.frame.PlugInFrame;
import ij.ImagePlus;

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


public class Aplim_Example  extends PlugInFrame{

	private static final long serialVersionUID = 1L;

	public static void main(String[]args) {
		System.out.println("This is the starting point when testing the code from Eclipse");
		ImageJ ij=new ImageJ();
		Aplim_Example ex=new Aplim_Example();
		ex.run("");
	}
	
	public Aplim_Example() {
		super("Title of the Frame");
	}

	public void run(String arg) {
		System.out.println("This is the starting point when clicking on the plugin in ImageJ menus");
		IJ.showMessage("This button don't do nothing. It is here just for pedagogical purpose. See you !");
	}
}
