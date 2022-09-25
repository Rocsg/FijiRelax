package io.github.rocsg.fijirelax.test;
import io.github.rocsg.fijiyama.common.VitiDialogs;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;

import java.io.File;
//import trainableSegmentation.*;
public class OLD_estScripts {

	public OLD_estScripts() {
		// TODO Auto-generated constructor stub
	}
	
	
	public static ImagePlus hyperStackingCZT(ImagePlus[][][]img) {
		Concatenator con=new Concatenator();
		con.setIm5D(true);
		int C=img.length;
		int Z=img[0].length;
		int T=img[0][0].length;
		ImagePlus[]tab=new ImagePlus[C*Z*T];
		for(int c=0;c<C;c++)for(int z=0;z<Z;z++)for(int t=0;t<T;t++) {
			System.out.println(c+","+z+","+t+" -> "+(c*(Z*T)+z*T+t)+" / "+(C*Z*T));
			tab[c*(Z*T)+z*T+t]=img[c][z][t];
		}
		ImagePlus hypTemp=con.concatenate(tab,true);
		String codeStacking="xyczt";
		ImagePlus hyperImage=HyperStackConverter.toHyperStack(hypTemp, C, Z,T,codeStacking,"Grayscale");
		return hyperImage;
	}


	
	
	public static void thing() {
		//Input data definition
		String mainDir="/home/fernandr/Bureau/FijiRelax_DOI";
		File f=new File(mainDir,"Experiments_for_figures_and_tables");
		f=new File(f.getAbsolutePath(),"Input_data");
		f=new File(f.getAbsolutePath(),"Images");
		f=new File(f.getAbsolutePath(),"Squared_fantoms");
		File f1=new File(f.getAbsolutePath(),"Simulated_echoes.tif");
		ImagePlus img=IJ.openImage(f1.getAbsolutePath());
		img.show();
		ImagePlus[]tab=new ImagePlus[3];
		double []valsSNR=new double[] {200,500,1000};
		double []valsSigma=new double[valsSNR.length];
		double PD_exp02=1000;

		for(int t=0;t<3;t++) {
			tab[t]=img.duplicate();
			valsSigma[t]=Math.sqrt(PD_exp02*PD_exp02/valsSNR[t]);
			for(int c=0;c<img.getNChannels();c++) {
				for(int z=0;z<img.getNSlices();z++) {
					System.out.println(c+","+z+","+t);
					int sli=VitimageUtils.getCorrespondingSliceInHyperImage(tab[t], c, z, 0);
					tab[t].getStack().setSliceLabel("T2SEQ_Fantom_TR=10000_TE="+(10*(c+1))+"_SIGMARICE="+VitimageUtils.dou(valsSigma[t]), sli);
				}
			}
		}
		ImagePlus hyper=VitimageUtils.hyperStackingFrames(tab);
		hyper.show();
		hyper.setTitle("hyper");
		IJ.saveAsTiff(hyper, f1.getAbsolutePath());
	}
	
	
	public static String getFijiRelaxDOIDir() {
		return VitiDialogs.chooseDirectoryUI("Select FijiRelax DOI directory","Dir that contains Case_0_Brain, Case_1 ...");
	}
	
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		String lut="physics";
//		if(true)return;
		//switch variable to use or not use visual debug information
		boolean useVisualDebugThatIsVeryInterestingButImplySlowerComputation=true;

		//Get FijiRelax archive directory path
		String mainDir="";
		if (new File("/home/fernandr/Bureau/FijiRelax_DOI/debug").exists())mainDir="/home/fernandr/Bureau/FijiRelax_DOI";
		else mainDir=getFijiRelaxDOIDir();

		//Input data definition
		File f=new File(mainDir,"Experiments_for_figures_and_tables");
		f=new File(f.getAbsolutePath(),"Input_data");
		f=new File(f.getAbsolutePath(),"Images");
		f=new File(f.getAbsolutePath(),"Squared_fantoms");
		File f1=new File(f.getAbsolutePath(),"Simulated_echoes.tif");
		String sourcePathSimulatedEchoes=f1.getAbsolutePath();
		f1=new File(f.getAbsolutePath(),"Expected_maps.tif");
		String sourcePathExpectedMaps=f1.getAbsolutePath();

		
		//Output image definition
		f=new File(mainDir,"Experiments_for_figures_and_tables");
		f=new File(f.getAbsolutePath(),"Results_and_evaluation");
		f=new File(f.getAbsolutePath(),"My_results");
		f1=new File(f.getAbsolutePath(),"Images");
		File f2 =new File(f1.getAbsolutePath(),"Result_hypermap_figure_1d_and_S7");
		String targetPathHyperMap=f2.getAbsolutePath();

		f1=new File(f.getAbsolutePath(),"Others");
		f2 =new File(f1.getAbsolutePath(),"Result_table_S2.csv");
		String targetPathCSV=f2.getAbsolutePath();

		ImagePlus imgEchoes=IJ.openImage(sourcePathSimulatedEchoes);
		imgEchoes.show();
		imgEchoes.setTitle("Simulated_echoes");
		

		//Expected fit
		ImagePlus imgExpected=new Duplicator().run(IJ.openImage(sourcePathExpectedMaps),3,3,1,12,1,1);
		imgExpected=new Duplicator().run(imgExpected,3,3,1,12,1,1);
		ImagePlus imgMask=VitimageUtils.thresholdImage(imgExpected, 0.000001, 10E8);
		imgMask.show();
		imgMask.setTitle("Img mask");
		imgMask.setDisplayRange(0, 1);
		imgExpected.show();
		imgExpected.setTitle("Expected FIT");
		imgExpected.setDisplayRange(5, 100);
		ImagePlus imgDIFFExpected=VitimageUtils.makeOperationBetweenTwoImages(imgExpected, imgExpected, 3,true);
		imgDIFFExpected=VitimageUtils.makeOperationBetweenTwoImages(imgDIFFExpected, imgMask, 2,true);
		imgDIFFExpected.setTitle("DIFF Expected");
		IJ.run(imgDIFFExpected, lut, "");
		imgDIFFExpected.setDisplayRange(0.8, 1.2);
		
		//Results without noise handling
		imgEchoes.duplicate().show();
		ImagePlus imgEXP=ValidationExperimentsFijiRelaxPaper.estimateMapsFromEchoesSquaredFantoms(imgEchoes,MRUtils.T2_MONO);
		imgEXP.duplicate().show();
		imgEXP=new Duplicator().run(imgEXP,3,3,1,12,1,1);
		imgEXP.show();
		imgEXP.setDisplayRange(5, 100);
		imgEXP.setTitle("Results_EXP_FIT");
		ImagePlus imgDIFFEXP=VitimageUtils.makeOperationBetweenTwoImages(imgEXP, imgExpected, 3,true);
		imgDIFFEXP=VitimageUtils.makeOperationBetweenTwoImages(imgDIFFEXP, imgMask, 2,true);
		imgDIFFEXP.show();
		imgDIFFEXP.setTitle("DIFF EXP FIT");
		IJ.run(imgDIFFEXP, lut, "");
		imgDIFFEXP.setDisplayRange(0.8, 1.2);
		
		//Results with Offset fit
		ImagePlus imgOFFSET=ValidationExperimentsFijiRelaxPaper.estimateMapsFromEchoesSquaredFantoms(imgEchoes,MRUtils.T2_MONO_BIAS);
		imgOFFSET=new Duplicator().run(imgOFFSET,3,3,1,12,1,1);
		imgOFFSET.show();
		imgOFFSET.setDisplayRange(5, 100);
		imgOFFSET.setTitle("Results_OFFSET_FIT");
		imgOFFSET.setDisplayRange(5, 100);
		ImagePlus imgDIFFOFFSET=VitimageUtils.makeOperationBetweenTwoImages(imgOFFSET, imgExpected, 3,true);
		imgDIFFOFFSET=VitimageUtils.makeOperationBetweenTwoImages(imgDIFFOFFSET, imgMask, 2,true);
		imgDIFFOFFSET.show();
		imgDIFFOFFSET.setTitle("DIFF EXP FIT");
		IJ.run(imgDIFFOFFSET, lut, "");
		imgDIFFOFFSET.setDisplayRange(0.8, 1.2);

		//Results with Rice fit
		ImagePlus imgRICE=ValidationExperimentsFijiRelaxPaper.estimateMapsFromEchoesSquaredFantoms(imgEchoes,MRUtils.T2_MONO_RICE);
		imgRICE=new Duplicator().run(imgRICE,3,3,1,12,1,1);
		imgRICE.show();
		imgRICE.setTitle("Results_RICE_FIT");
//		imgRICE.setC(2);
		imgRICE.setDisplayRange(5, 100);
		ImagePlus imgDIFFRICE=VitimageUtils.makeOperationBetweenTwoImages(imgRICE, imgExpected, 3,true);
		imgDIFFRICE=VitimageUtils.makeOperationBetweenTwoImages(imgDIFFRICE, imgMask, 2,true);
		imgDIFFRICE.show();
		imgDIFFRICE.setTitle("DIFF RICE FIT");
		IJ.run(imgDIFFRICE, "Fire", "");
		imgDIFFRICE.setDisplayRange(0.8, 1.2);
		
		//Here is the figure 1a and S2
		IJ.saveAsTiff(imgExpected,targetPathHyperMap+"_Expected.tif");
		IJ.saveAsTiff(imgDIFFExpected,targetPathHyperMap+"_DIFF_Expected.tif");

		IJ.saveAsTiff(imgEXP,targetPathHyperMap+"_Exponential.tif");
		IJ.saveAsTiff(imgDIFFEXP,targetPathHyperMap+"_DIFF_Exponential.tif");

		IJ.saveAsTiff(imgOFFSET,targetPathHyperMap+"_Offset.tif");
		IJ.saveAsTiff(imgDIFFOFFSET,targetPathHyperMap+"_DIFF_Offset.tif");

		IJ.saveAsTiff(imgRICE,targetPathHyperMap+"_Rice.tif");
		IJ.saveAsTiff(imgDIFFRICE,targetPathHyperMap+"_DIFF_Rice.tif");

		ValidationExperimentsFijiRelaxPaper.processResultsSquaredFantoms(new ImagePlus[] {imgExpected,imgEXP,imgOFFSET,imgRICE},targetPathCSV);
		
		


		IJ.showMessage("Your image results has been saved to \n"+targetPathHyperMap+" \n \n.\n.\nYour table results has been saved to \n"+targetPathCSV+"\nComputation done !");
	}
	
	
	
	
}
