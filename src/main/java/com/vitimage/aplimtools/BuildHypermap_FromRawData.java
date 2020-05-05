package com.vitimage.aplimtools;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

import com.vitimage.common.Timer;
import com.vitimage.common.VitiDialogs;
import com.vitimage.mrutils.HyperMRIT1T2;

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


public class BuildHypermap_FromRawData  extends PlugInFrame{

	private static final long serialVersionUID = 1L;
	private String inputDir;
	private String outputDir;
	
	public static void main(String[]args) {
		System.out.println("This is the starting point when testing the code from Eclipse");
		ImageJ ij=new ImageJ();
		boolean testBuildFromBionanoDir=true;
		if(testBuildFromBionanoDir) {
			BuildHypermap_FromRawData ex=new BuildHypermap_FromRawData();
			ex.run(new Object[] {"/home/fernandr/Bureau/Traitements/Sorgho/Donnees_brutes_export_Romain/BM1/BM1_F0_1202_SL_RAW",
					"/home/fernandr/Bureau/Traitements/Sorgho/TestStuff","BM1_1202"});
		}
		else {
			BuildHypermap_FromRawData ex=new BuildHypermap_FromRawData();
			ex.run("");
		}
	}
	
	public BuildHypermap_FromRawData() {
		super("Title of the Frame");
	}

	public void run(String arg) {
		run(new Object[] {});
	}

	public void run(Object[]obj) {
		Timer t=new Timer();
		String nameObservation="";
		//Ask for the directory input, and the directory output
		IJ.log("Opening input and output dirs");
		if(obj != null && obj.length==3){
			inputDir=(String)obj[0];
			outputDir=(String)obj[1];
			nameObservation=(String)obj[2];
		}
		else{
			inputDir=VitiDialogs.chooseDirectoryUI("Localize inputDir (SL_RAW)", "Select this directory");
			outputDir=VitiDialogs.chooseDirectoryUI("Select outputDir", "Select this directory");
		}

		T1T2Seq_Importer t1t2=new T1T2Seq_Importer();
		ImagePlus hyperMap=t1t2.run(inputDir,nameObservation);
		HyperMRIT1T2 hyp=new HyperMRIT1T2(hyperMap);
		
		hyp.computeMapsAgain();
		hyperMap=hyp.getCopy();
		hyperMap.setTitle("Hyperimage");
		String exportFileName=new File(outputDir,"hyperimage_build_time_"+(new SimpleDateFormat("yyyy-MM-dd_hh-mm").format(new Date())+".tif")).getAbsolutePath();
		IJ.saveAsTiff(hyperMap,exportFileName);
	}
}
