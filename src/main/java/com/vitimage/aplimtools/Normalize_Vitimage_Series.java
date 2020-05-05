package com.vitimage.aplimtools;
import java.io.File;
import java.util.ArrayList;

import com.vitimage.common.VitimageUtils;
import com.vitimage.mrutils.MRUtils;

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


public class Normalize_Vitimage_Series  extends PlugInFrame{

	private static final long serialVersionUID = 1L;

	public static void main(String[]args) {
		System.out.println("This is the starting point when testing the code from Eclipse");
		ImageJ ij=new ImageJ();
		Normalize_Vitimage_Series ex=new Normalize_Vitimage_Series();
		ex.run("");
	}
	
	public Normalize_Vitimage_Series() {
		super("Title of the Frame");
	}

	
	
	
	

	
	
	
	
	public void run(String arg) {
		boolean separateSlices=true;
		String prefix="";
		boolean onServer=false;
		if(new File("/home/fernandr/").exists())prefix="/home/fernandr/Bureau/Traitements/Sorgho/";
		else if (new File("/users/bionanonmri/").exists()) {
			prefix="/users/bionanonmri/fernandez/dataBamba/";
			onServer=true;
		}
		else return;
		String sourceDir=prefix+"Cartes_rangees_par_specimen";
		String targetDir=prefix+"Cartes_normalisees_capillaire";
		new File(targetDir).mkdirs();
		
		String[]specs=new String[] {"BM1","BM2","SSM1","SSM2"};
		for(String sp : specs) {
			
			IJ.log("\nTraitement des series de l'individu "+sp);
			String spDir=new File(sourceDir,sp).getAbsolutePath();
			String spOutDir=new File(targetDir,sp).getAbsolutePath();
			new File(spOutDir).mkdirs();
			String []strTimepoints=new File(spDir).list();
			for(String tim : strTimepoints) {
				IJ.log("Traitement des series de l'individu "+sp+", observation "+tim);
				String inputName=new File(spDir,tim).getAbsolutePath();
				String outputName=new File(spOutDir,tim).getAbsolutePath();
				System.out.println("Processing "+inputName+" to set result in "+outputName);
				ImagePlus hyper=IJ.openImage(inputName);
				int Z=hyper.getNSlices();
				int X=hyper.getWidth();
				int Y=hyper.getHeight();
				int C=hyper.getNChannels();
				ImagePlus imgMaskCap=new Duplicator().run(hyper,5,5,1,Z,1,1);
				float[]capVals;
				int nHits;
				int incr;
				double sum=0;
				for(int z=0;z<Z;z++) {
					double []valCap=MRUtils.readCapValuesInSliceLabel(hyper,0,z,0);
					sum+=valCap[0];
				}
				double valMean=sum/Z;
				System.out.println("Valeur du capillaire "+valMean);
				double factorToApply=MRUtils.maxM0ForNormalization/valMean;
				for(int c=0;c<C;c++) {
					if(c>0 && c<5)continue;
					System.out.print(" "+c);
					for(int z=0;z<Z;z++) {
						int indMask=VitimageUtils.getCorrespondingSliceInHyperImage(hyper, 3, z, 0);
						float[]valsMask=(float[])hyper.getStack().getProcessor(indMask).getPixels();
						int indSlice=VitimageUtils.getCorrespondingSliceInHyperImage(hyper, c, z, 0);
						float[]valsImg=(float[])hyper.getStack().getProcessor(indSlice).getPixels();
						for(int x=0;x<X;x++) {
							for(int y=0;y<Y;y++) {
								valsImg[y*X+x]*=factorToApply;
								if(valsMask[y*X+x]==0 && c==0)valsImg[y*X+x]=0;
							}
						}
						hyper.getStack().setSliceLabel(hyper.getStack().getSliceLabel(indSlice).replace("CAP", "CAP-INIT"), indSlice);
					}
				}
				System.out.println();
				IJ.saveAsTiff(hyper, outputName);
			}
		}
	}

	
	
	
	
	
}
