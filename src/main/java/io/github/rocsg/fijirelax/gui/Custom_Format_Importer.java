/*
 * 
 */
package io.github.rocsg.fijirelax.gui;

import java.io.File;
import io.github.rocsg.fijiyama.common.Timer;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.Concatenator;
import ij.plugin.FolderOpener;
import ij.plugin.HyperStackConverter;
import ij.plugin.LutLoader;
import ij.process.ImageConverter;
import io.github.rocsg.fijirelax.mrialgo.HyperMap;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;

/**
 * The Custom_Format_Importer class is in charge of handling importation MRI data
 * from NIFTI, multiple NIFTI or BRUKER formats.
 *
 *  	Copyright (C) 2022  io.github.rocsg
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see https://www.gnu.org/licenses/
 * @author Romain Fernandez (romain.fernandezATcirad.fr)
 * @version 1.2, 15.04.2022
 */
public class Custom_Format_Importer{
	
	/** The input dir. */
	private String inputDir="";
	
	/** The name observation. */
	public String nameObservation="";
	
	/** The voxs. */
	private double[]voxs;
	
	/** The img T 2 seq. */
	ImagePlus []imgT2Seq;
	
	/** The img T 1 seq. */
	ImagePlus []imgT1Seq;
	
	/** The img T 2 seq test. */
	ImagePlus []imgT2SeqTest;
	
	/** The img T 1 seq test. */
	ImagePlus []imgT1SeqTest;
	
	/** The M 0 map. */
	ImagePlus M0map;
	
	/** The M 0 map 2. */
	ImagePlus M0map2;
	
	/** The T 1 map. */
	ImagePlus T1map;
	
	/** The T 1 map 2. */
	ImagePlus T1map2;
	
	/** The T 2 map. */
	ImagePlus T2map;
	
	/** The T 2 map 2. */
	ImagePlus T2map2;
	
	/** The export file name. */
	public String exportFileName="";
	
	/** The force no man. */
	boolean forceNoMan=true;
	
	/** The make crop. */
	public boolean makeCrop=false;
	
	/** The skip reg for testing. */
	boolean skipRegForTesting=false;
	
	/** The view registration. */
	boolean viewRegistration=false;

	/** The normalize before computation. */
	public boolean normalizeBeforeComputation=false;
	
	/** The normalize after computation. */
	public boolean normalizeAfterComputation=false;

	/** The force no bouture trick. */
	public boolean forceNoBoutureTrick=false;
	
	/** The make sorgho trick. */
	public boolean makeSorghoTrick=false;
	
	/** The make bouture trick. */
	public boolean makeBoutureTrick;
	
	/** The forget early reps. */
	public boolean forgetEarlyReps=false;
	
	/** The speedup registration for testing. */
	public boolean speedupRegistrationForTesting=true;//false;
	
	/** The hide message. */
	public boolean hideMessage=true;
	
	/** The dont show nothing. */
	public boolean dontShowNothing=false;
	
	/** The compute multi. */
	public boolean computeMulti=true;
	
	/** The img T 1 T 2. */
	private ImagePlus[][] imgT1T2;
	
	/** The img T 1 T 2 line. */
	private ImagePlus[] imgT1T2Line;
	
	/** The multi M 0. */
	private ImagePlus multiM0;
	
	/** The multi T 1. */
	private ImagePlus multiT1;
	
	/** The name. */
	private String name;
	
	/**
	 * Instantiates a new custom format importer.
	 *
	 * @param inputDir the input dir
	 * @param name the name
	 */
	public Custom_Format_Importer(String inputDir,String name) {
		this.inputDir=inputDir;
		this.name=name;
	}
	
	
	/**
	 * Run the importation, according to defined member fields (name,inputdir)
	 * and return the imported MRI observation as an HyperMap.
	 *
	 * @return The reconstructed HyperMap
	 * @see HyperMap
	 */
	public HyperMap importCustomDataAsHyperMap() {
		Timer t=new Timer();
		this.nameObservation=name;
		this.inputDir=inputDir;

		//Ask for the directory input, and the directory output
		//then ead the voxel sizes and the image sizes		
		IJ.log("Opening input and output dirs");		
		IJ.log("Gathering image informations");
		Object []objs=getFirstImageOfT1T2Serie(inputDir);
		ImagePlus testDcm=(ImagePlus)objs[0];
		this.voxs=VitimageUtils.getVoxelSizes(testDcm);

		t.print("Reading images");
		return readT1T2();
	}
	
	
	/**
	 * Set member fields maps (Images) to null
	 */
	public void makeNullMaps(){
		multiM0=VitimageUtils.nullImage(imgT1T2Line[0]);
		multiT1=VitimageUtils.nullImage(imgT1T2Line[0]);
		T2map=VitimageUtils.nullImage(imgT1T2Line[0]);
	}
	

	/**
	 * Gets the first image of T1-T2 serie, when providing the input dir containing all the data.
	 *
	 * @param inputDir the input dir
	 * @return the first image of T 1 T 2 serie
	 */
	public Object[] getFirstImageOfT1T2Serie(String inputDir) {
		String imgPath=inputDir;
		String[]subFile=null;
		while(!imgPath.substring(imgPath.length()-4, imgPath.length()).equals(".dcm") && !imgPath.substring(imgPath.length()-4, imgPath.length()).equals(".tif")){
			//IJ.log("Diving into subdirs, not yet an img :"+imgPath);
			subFile=VitimageUtils.stringArraySort(new File(imgPath).list());
			imgPath=new File(imgPath,subFile[0]).getAbsolutePath();
		}
		return new Object[] {IJ.openImage(imgPath),subFile.length};
	}
	
	/**
	 * Read T1-T2 data from the given input dir (object member),organize it and return it as an HyperMap.
	 *
	 * @return the hyper map
	 */
	public HyperMap readT1T2() {
		String[]listFiles=VitimageUtils.stringArraySortByTrValue(new File(inputDir).list());
		imgT1T2 =new ImagePlus[listFiles.length][];
		int nbSelected=listFiles.length;
		int indexStart=listFiles.length-nbSelected;
		int numberTot=0;

		ImageConverter ic;
		//Reading successive Tr sequences
		for(int ii=indexStart;ii<listFiles.length;ii++) {
			int i=ii-indexStart;
			String[]strTes=VitimageUtils.stringArraySort(new File(inputDir,listFiles[ii]).list());
			imgT1T2[i]=new ImagePlus[strTes.length];

			//Reading successive Te in each Tr
			for(int j=0;j<strTes.length;j++) {
				numberTot++;
				String str=new File(new File(inputDir,listFiles[ii]),strTes[j]).getAbsolutePath();
				imgT1T2[i][j]=FolderOpener.open(str, "");
				if(i==0 && j==0)detectSorgho(imgT1T2[i][j]);
				if(i==0 && j==0)detectBouture(imgT1T2[i][j]);
				int nbAverage=VitimageUtils.getAveraging(imgT1T2[i][j]);
				ic=new ImageConverter(imgT1T2[i][j]);
				ic.convertToGray32();
//				IJ.run(imgT1T2[i][j],"32-bit","");
				if(nbAverage!=2)imgT1T2[i][j]=VitimageUtils.makeOperationOnOneImage(imgT1T2[i][j], 2, 2.0/nbAverage, false);
			}
		}

		
		//Compute sigma and label the slices
		ImagePlus[]tab=new ImagePlus[numberTot];
		int Z=imgT1T2[0][0].getNSlices();
		int incr=0;
		for(int i=0;i<imgT1T2.length;i++) {			
			String type="";
			if(i==(imgT1T2.length-1) && (imgT1T2[i].length>1))type="T2SEQ_";
			else type="T1SEQ_";
			for(int j=0;j<imgT1T2[i].length;j++) {
				MRUtils.computeTeTrAndRiceSigmaOfEachSliceAndWriteItInTheLabels(imgT1T2[i][j],false,type+nameObservation,makeBoutureTrick);
				tab[incr++]=imgT1T2[i][j];
			}
		}

		//Convert to hyperStack
		double valMax=VitimageUtils.maxOfImage(VitimageUtils.maxOfImageArrayDouble(imgT1T2));
		ImagePlus newHyperImg=Concatenator.run(tab);
		VitimageUtils.printImageResume(newHyperImg,"Hyper");
		newHyperImg=HyperStackConverter.toHyperStack(newHyperImg, tab.length,Z,1,"xyztc","Grayscale");		
		
		//Set display range and lut
		for(int c=0;c<tab.length;c++) {
			newHyperImg.setC(c+1);
			WindowManager.setTempCurrentImage(newHyperImg);
			new LutLoader().run("fire");
			
//			IJ.run(newHyperImg,"Fire","");
			newHyperImg.setDisplayRange(0, valMax);
		}
		
		return new HyperMap(newHyperImg);
	}		
	
	
	/**
	 * Used for detecting specific kind of experiments
	 *
	 * @param img the img
	 * @return true, if successful
	 */
	public boolean detectSorgho(ImagePlus img) {
		String patientName=VitimageUtils.getPatientName(img);
		if( patientName.contains("BM") || patientName.contains("SSM")) {
			/*IJ.log("Detected experience : Sorgho\nSettings : "+
					"\nForget early reps = false (dismiss TR<=1000 ms)"+
					"\nNormalize before computation = true (correct uneven normalization issues before maps computation)"+
					"\nMake sorgho trick = true (set all images in a common geometry to help analysis)\n\n");*/
			makeSorghoTrick=true;
			return true;
		}
		return false;
	}

	/**
	 * Used for detecting specific kind of experiments
	 *
	 * @param img the img
	 * @return true, if successful
	 */
	public boolean detectBouture(ImagePlus img) {
		boolean ret=VitimageUtils.isBouture(img);
		if(!ret && this.nameObservation.contains("BOUT"))ret=true;
		VitimageUtils.waitFor(2000);
		makeBoutureTrick=ret;
		return ret;
	}


	
	
}