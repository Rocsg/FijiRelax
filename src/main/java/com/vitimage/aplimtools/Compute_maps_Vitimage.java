package com.vitimage.aplimtools;
import java.io.File;
import java.util.ArrayList;

import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;

import ij.IJ;
import ij.ImageJ;
import ij.plugin.Duplicator;
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


public class Compute_maps_Vitimage  extends PlugInFrame{

	private static final long serialVersionUID = 1L;

	public static void main(String[]args) {
		System.out.println("This is the starting point when testing the code from Eclipse");
		ImageJ ij=new ImageJ();
		Compute_maps_Vitimage ex=new Compute_maps_Vitimage();
		ex.run("");
	}
	
	public Compute_maps_Vitimage() {
		super("Title of the Frame");
	}

	
	
	
	

	
	
	//TODO : multithread both
	
	public void run(String arg) {
		double part=VitiDialogs.getDoubleUI("Choose a specimen to process", "Answer : \n 0 to process all,\n 1 to process only BM1,\n 2 to process only BM2,\n 3 to process only SSM1, \n 4 to process only SSM2", 0);
		String prefix="";
		boolean onServer=false;
		if(new File("/home/fernandr/").exists())prefix="/home/fernandr/Bureau/Traitements/Bouture6D/Test/";
		else if (new File("/users/bionanonmri/").exists()) {
			prefix="/users/bionanonmri/fernandez/Vitimage/";
			onServer=true;
		}
		else return;
		String sourceDir=prefix+"Donnees_brutes";
		String targetDir=prefix+"Cartes_calculees";
		String[]specs=null;
		specs=new String[] {"B099_PCH"};
		for(String sp : specs) {
			//IJ.log("Traitement des series de l'individu "+sp);
			String spDir=new File(sourceDir,sp).getAbsolutePath();
			String []strTimepoints=new File(spDir).list();
			System.out.println(spDir);
			
			String outputDir=new File(targetDir,sp).getAbsolutePath();
			new File(outputDir).mkdirs();
			String tim=strTimepoints[(int)Math.round(part)];
			IJ.log("Traitement des series de l'individu "+sp+", observation "+tim);
			T1T2Seq_Importer bn=new T1T2Seq_Importer();
			bn.nameObservation="VITIMAGE_"+tim;
			bn.exportFileName=new File(outputDir,sp+"_"+tim.split("_")[2]+".tif").getAbsolutePath();
			System.out.println("Observation name will be "+bn.nameObservation);
			System.out.println("Saving will be in "+bn.exportFileName);
			bn.viewRegistration=false;
			bn.hideMessage=true;
			bn.forgetEarlyReps=false;
			bn.dontShowNothing=true;
			bn.speedupRegistrationForTesting=true;
			bn.forgetEarlyReps=false;
			bn.run(new File(spDir,tim).getAbsolutePath(),bn.exportFileName);
		}
		for(int i=0;i<10;i++)IJ.log("\n");
		IJ.log("Batch processing of Vitimage is now finished. End of process");
	}
	
	
	
	
	
	
	
	
	
	
	
	public void correctifCapillary(String arg) {
		IJ.log("New version in use");
		String prefix="";
		boolean onServer=false;
		if(new File("/home/fernandr/").exists())prefix="/home/fernandr/Bureau/Traitements/Sorgho/";
		else if (new File("/users/bionanonmri/").exists()) {
			prefix="/users/bionanonmri/fernandez/dataBamba/";
			onServer=true;
		}
		else return;
		String sourceDir=prefix+"Donnees_brutes_export_Romain";
		String targetDir=prefix+"Cartes_rangees_par_specimen";
		String[]specs=new String[] {"BM1","BM2","SSM1","SSM2"};
		for(String sp : specs) {
			IJ.log("\nTraitement des series de l'individu "+sp);
			String spDir=new File(sourceDir,sp).getAbsolutePath();
			String []strTimepoints=new File(spDir).list();
			String outputDir=new File(targetDir,sp).getAbsolutePath();
			new File(outputDir).mkdirs();
			for(String tim : strTimepoints) {
				IJ.log("Traitement des series de l'individu "+sp+", observation "+tim);
				String exportFileName=new File(outputDir,sp+"_"+tim.split("_")[2]+".tif").getAbsolutePath();
				System.out.println("Processing "+exportFileName);
				ImagePlus hyper=IJ.openImage(exportFileName);
				int Z=hyper.getNSlices();
				int X=hyper.getWidth();
				int Y=hyper.getHeight();
				int C=hyper.getNChannels();
				ImagePlus imgMaskCap=new Duplicator().run(hyper,5,5,1,Z,1,1);
				float[]capVals;
				int nHits;
				int incr;
				double sum;
				for(int z=0;z<Z;z++) {
					nHits=0;
					ArrayList<int[]> coordsPointsCap=new ArrayList<int[]>();
					capVals=(float[])imgMaskCap.getStack().getProcessor(z+1).getPixels();
					for(int x=0;x<X;x++)for(int y=0;y<Y;y++) if(capVals[VitimageUtils.getPixelIndex(x, X, y)]<-0.5) {nHits++; coordsPointsCap.add(new int[] {x,y});};

					for(int c=0;c<C;c++) {
						double[]vals=new double[nHits];
						incr=0;
						capVals=(float[])hyper.getStack().getProcessor(VitimageUtils.getCorrespondingSliceInHyperImage(hyper, c, z, 0)).getPixels();				
						for(int i=0;i<nHits;i++) {
							vals[incr++]=capVals[VitimageUtils.getPixelIndex(coordsPointsCap.get(i)[0], X, coordsPointsCap.get(i)[1])];
						}
						double[]stats=VitimageUtils.statistics1D(vals);
						hyper.getStack().setSliceLabel(hyper.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyper, c, z, 0))+"_CAP="+VitimageUtils.dou(stats[0])+"+-"+VitimageUtils.dou(stats[1]), VitimageUtils.getCorrespondingSliceInHyperImage(hyper, c, z, 0));
					}
				}	
				IJ.saveAsTiff(hyper, exportFileName);
			}
		}
	}

	
	
	
	
	
}
