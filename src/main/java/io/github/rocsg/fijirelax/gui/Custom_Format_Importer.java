/*
 * 
 */
package io.github.rocsg.fijirelax.gui;

import java.io.File;
import io.github.rocsg.fijiyama.common.Timer;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.FolderOpener;
import ij.plugin.HyperStackConverter;
import io.github.rocsg.fijirelax.mrialgo.HyperMap;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;

// TODO: Auto-generated Javadoc
//TODO : is this class still useful ?

/**
 * The Custom_Format_Importer class is in charge of handling importation MRI data
 * from NIFTI, multiple NIFTI or BRUKER formats.
 *
 * @author Romain Fernandez (romain.fernandez@cirad.fr)
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
	 * Set all data member fields to null.
	 */
	public void makeNullMaps(){
		multiM0=VitimageUtils.nullImage(imgT1T2Line[0]);
		multiT1=VitimageUtils.nullImage(imgT1T2Line[0]);
		T2map=VitimageUtils.nullImage(imgT1T2Line[0]);
	}
	

	/**
	 * Gets the first image of T 1 T 2 serie.
	 *
	 * @param inputDir the input dir
	 * @return the first image of T 1 T 2 serie
	 */
	public Object[] getFirstImageOfT1T2Serie(String inputDir) {
		String imgPath=inputDir;
		String[]subFile=null;
		System.out.println(imgPath);
		while(!imgPath.substring(imgPath.length()-4, imgPath.length()).equals(".dcm") && !imgPath.substring(imgPath.length()-4, imgPath.length()).equals(".tif")){
			IJ.log("Not yet an img :"+imgPath);
			subFile=VitimageUtils.stringArraySort(new File(imgPath).list());
			System.out.println(subFile[0]);
			imgPath=new File(imgPath,subFile[0]).getAbsolutePath();
		}
		return new Object[] {IJ.openImage(imgPath),subFile.length};
	}
	
	/**
	 * Read T 1 T 2.
	 *
	 * @return the hyper map
	 */
	public HyperMap readT1T2() {
		System.out.println(inputDir);
		String[]listFiles=VitimageUtils.stringArraySortByTrValue(new File(inputDir).list());
		imgT1T2 =new ImagePlus[listFiles.length][];
		int nbSelected=listFiles.length;
		int indexStart=listFiles.length-nbSelected;
		int numberTot=0;

		//Reading successive Tr sequences
		for(int ii=indexStart;ii<listFiles.length;ii++) {
			int i=ii-indexStart;
			String[]strTes=VitimageUtils.stringArraySort(new File(inputDir,listFiles[ii]).list());
			imgT1T2[i]=new ImagePlus[strTes.length];

			//Reading successive Te in each Tr
			for(int j=0;j<strTes.length;j++) {
				numberTot++;
				String str=new File(new File(inputDir,listFiles[ii]),strTes[j]).getAbsolutePath();
				System.out.println("Opening MRI Data : "+str);						
				imgT1T2[i][j]=FolderOpener.open(str, "");
				if(i==0 && j==0)detectSorgho(imgT1T2[i][j]);
				if(i==0 && j==0)detectBouture(imgT1T2[i][j]);
				int nbAverage=VitimageUtils.getAveraging(imgT1T2[i][j]);
				IJ.run(imgT1T2[i][j],"32-bit","");
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
		for(int i=0;i<tab.length;i++)VitimageUtils.printImageResume(tab[i],""+i);
		ImagePlus newHyperImg=Concatenator.run(tab);
		VitimageUtils.printImageResume(newHyperImg,"Hyper");
		newHyperImg=HyperStackConverter.toHyperStack(newHyperImg, tab.length,Z,1,"xyztc","Grayscale");		
		
		//Set display range and lut
		for(int c=0;c<tab.length;c++) {
			newHyperImg.setC(c+1);
			IJ.run(newHyperImg,"Fire","");
			newHyperImg.setDisplayRange(0, valMax);
		}
		
		return new HyperMap(newHyperImg);
	}		
	
	
	/**
	 * Detect sorgho.
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
	 * Detect bouture.
	 *
	 * @param img the img
	 * @return true, if successful
	 */
	public boolean detectBouture(ImagePlus img) {
		System.out.println("Investigating if is Bouture");
		boolean ret=VitimageUtils.isBouture(img);
		System.out.println(ret);
		System.out.println();
		if(!ret && this.nameObservation.contains("BOUT"))ret=true;
		VitimageUtils.waitFor(2000);
		makeBoutureTrick=ret;
		return ret;
	}


		
	/**
	 * Register the T1 observations along with T2 observations, and update member fields with the processed data
	 *
	 * 	 */
	//TODO : is this function still needed, as HyperMap have its own facility ?
	//TODO : there is a if(false) there

/*
	public void registerT1T2(){		
		int nT1=this.imgT1T2.length;
		if(nT1==1)return;
		ItkTransform  []trs=new ItkTransform[nT1];
		VitimageUtils.printImageResume(imgMask,"imgMask2");
		VitimageUtils.printImageResume(imgMaskCap,"imgMaskCap2");

		//In case of sorgho data, set them all to the same voxel size
		if(makeSorghoTrick && false) {
			System.out.println("Sorgho detect !");
			double targetVZ=0.5;
			double targetVoxVolume=1E-3;
			double targetVX=Math.sqrt(targetVoxVolume/targetVZ);
			double targetVY=Math.sqrt(targetVoxVolume/targetVZ);
			int[]targetDims=new int[] {512,512,4};
			double[]targetVoxs=new double[] {targetVX,targetVY,targetVZ};		
			for(int j=0;j<imgT1T2[nT1-1].length;j++) {
				imgT1T2[nT1-1][j]=new ItkTransform().transformImageExtensive(targetDims, targetVoxs, imgT1T2[nT1-1][j],false);
			}
			imgMask.getStack().setSliceLabel("", 1);
			imgMask=new ItkTransform().transformImage(targetDims, targetVoxs, imgMask,false);
			imgMask=VitimageUtils.getFloatBinaryMask(imgMask, 1.0, 1E10);
			imgMaskCap=VitimageUtils.getAntiCapillaryMask(imgT1T2[imgT1T2.length-1][0],VitimageUtils.bionanoCapillaryRadius);
		}
		VitimageUtils.printImageResume(imgMask,"imgMask3");
		VitimageUtils.printImageResume(imgMaskCap,"imgMaskCap3");
		imgMask=VitimageUtils.makeOperationBetweenTwoImages(imgMask,imgMaskCap, 1, true);
			
		
		trs[nT1-1]=new ItkTransform();
		ImagePlus imgMov,imgRef;
		ItkTransform tr=new ItkTransform();
		RegistrationAction regAct;
		BlockMatchingRegistration bmRegistration;
		imgMov=VitimageUtils.imageCopy(imgT1T2[nT1-2][0]);
		imgRef=VitimageUtils.imageCopy(imgT1T2[nT1-1][0]);
		
		//If needed, handle a manual registration before
		if(!forceNoMan || (nameObservation.contains("B032") && nameObservation.contains("J35"))) {
			ImagePlus comp=VitimageUtils.compositeOf(imgT1T2[nT1-1][0],imgT1T2[nT1-2][0],"red=T2 first echo , green = T1 last rep time");
			if(!dontShowNothing)comp.show();
			if(VitiDialogs.getYesNoUI("Manual registration first ?", "Manual registration first ?\nNeeded if T1 stack and T2 stack are largely disaligned")) {
				Point3d[][]pts=VitiDialogs.registrationPointsUI(5,imgRef,imgMov,true);
				tr=ItkTransform.estimateBestRigid3D(pts[1],pts[0]); 
			}
			if(!dontShowNothing)comp.close();
		}
		
		//Run automatic registration
		int incr=0;
		imgRef=imgT1T2[nT1-1][0];
		ImagePlus maskReg=new Duplicator().run(imgMask);
		maskReg=VitimageUtils.getFloatBinaryMask(maskReg, 2,10E8);
		maskReg=VitimageUtils.makeOperationOnOneImage(maskReg,1, -1, true);
		maskReg=VitimageUtils.makeOperationOnOneImage(maskReg,2, -1, true);

		for(int iMov=0;iMov<nT1-1;iMov++) {	
			if(skipRegForTesting) {
				trs[iMov]=new ItkTransform();
				continue;
			}
			imgMov=imgT1T2[iMov][0];
			regAct=new RegistrationAction();
			regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
			regAct.typeAutoDisplay=0;
			if(speedupRegistrationForTesting) {
				if(makeBoutureTrick)regAct.setLevelMax(1);
				regAct.strideX=regAct.strideX*3;
				regAct.strideY=regAct.strideY*3;
				regAct.strideZ=regAct.strideZ*3;
			};
			regAct.higherAcc=1;
			if(nameObservation.contains("B098") && nameObservation.contains("J70")) {
				regAct.iterationsBMLin+=10;
				regAct.strideX/=3;
				regAct.strideY/=3;
				regAct.strideZ-=1;
				regAct.levelMaxLinear+=1;
			}
			bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
			bmRegistration.consoleOutputActivated=false;
			bmRegistration.timingMeasurement=true;
			bmRegistration.mask=maskReg.duplicate();
			bmRegistration.setRefRange(new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()});
			bmRegistration.setMovRange(new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()});
			bmRegistration.flagRange=true;
			bmRegistration.minBlockVariance=10;
			bmRegistration.displayRegistration=viewRegistration ? 2 : 0;
			if(nameObservation.contains("B032") && nameObservation.contains("J35"))bmRegistration.displayRegistration=2;
			bmRegistration.displayR2=false;
			bmRegistration.levelMin=-1;
			bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
			trs[iMov]=bmRegistration.runBlockMatching(tr,false);
			bmRegistration.closeLastImages();
	
			regAct=new RegistrationAction();
			regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
			if(speedupRegistrationForTesting) {
				regAct.strideX=regAct.strideX*2;
				regAct.strideY=regAct.strideY*2;
				regAct.strideZ=regAct.strideZ*2;
			};
			regAct.higherAcc=1;
			
			if((!forgetEarlyReps) && (MRUtils.readTrInSliceLabel(imgMov,0,0,0)<=300) || makeBoutureTrick) {}
			else {
				regAct.levelMaxDense=1;
				regAct.typeTrans=Transform3DType.DENSE;
				bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
				bmRegistration.consoleOutputActivated=false;
				bmRegistration.timingMeasurement=true;
				bmRegistration.setRefRange(new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()});
				bmRegistration.setMovRange(new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()});
				bmRegistration.mask=maskReg.duplicate();
				bmRegistration.flagRange=true;
				bmRegistration.minBlockVariance=10;
				bmRegistration.displayRegistration=viewRegistration ? 2 : 0;
				bmRegistration.displayR2=false;
				bmRegistration.levelMin=-1;
				bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
				trs[iMov]=bmRegistration.runBlockMatching(trs[iMov],false);
				bmRegistration.closeLastImages();
			}
			for(int j=0;j<imgT1T2[iMov].length;j++) {
				imgT1T2[iMov][j]=trs[iMov].transformImageExtensive(imgT1T2[nT1-1][0],imgT1T2[iMov][j]);
			}
		}
		IJ.log("Registration finished !");
	}
	*/
	
	
	//TODO : this function is not anymore needed
	/*
	public void convert2DImgArrayTo1DArray() {
		int n=0;
		for(int i=0;i<imgT1T2.length;i++)n+=imgT1T2[i].length;
		imgT1T2Line=new ImagePlus[n];
		n=0;
		for(int i=0;i<imgT1T2.length;i++) {
			for(int j=0;j<imgT1T2[i].length;j++) {
				if( (i==imgT1T2.length-1) && (j==0) )indexFirstT2Seq=n;
				imgT1T2Line[n++]=(!makeCrop ? imgT1T2[i][j] : VitimageUtils.cropImageFloat(imgT1T2[i][j], 50, 50, 0, 100, 100, 4));
			}
		}
		dims=VitimageUtils.getDimensions(imgT1T2Line[0]);
	}
	
	
	public ImagePlus stackSaveAndShowResultsT1T2() {
		System.out.println("STARTING THE STACKING");
		if(makeSorghoTrick)nameObservation="SORGHO_"+nameObservation;
		if(makeBoutureTrick)nameObservation="BOUTURE_"+nameObservation;
		int nImg=imgT1T2Line.length;
		ImagePlus []tempRes=new ImagePlus[nImg+4];
		int curInd=0;
		int nMaps=3;
		tempRes[0]=multiM0.duplicate();		
		tempRes[1]=multiT1.duplicate();
		tempRes[2]=T2map.duplicate();
		tempRes[3]=imgMask;
		for(int z=1;z<=dims[2];z++) {
			tempRes[0].getStack().setSliceLabel("M0MAP_"+nameObservation, z);
			tempRes[1].getStack().setSliceLabel("T1MAP_"+nameObservation, z);
			tempRes[2].getStack().setSliceLabel("T2MAP_"+nameObservation, z);
			tempRes[3].getStack().setSliceLabel("MASKMAP_"+nameObservation, z);
		}
		
		for(int i=0;i<nImg;i++) {
			tempRes[i+4]=imgT1T2Line[i].duplicate();
			for(int z=1;z<=dims[2];z++) {
				tempRes[i+4].getStack().setSliceLabel(tempRes[i+4].getStack().getSliceLabel(z).replace("_T1T2SEQ","T1T2SEQ"), z);
			}
		}
		
	
		//Measurements of capillary value and std
		int Z=imgMask.getNSlices();
		int X=imgMask.getWidth();
		int Y=imgMask.getHeight();
		float[]capVals;
		int nHits;
		int incr;
		double sum;
		for(int z=0;z<Z;z++) {
			ArrayList<int[]> coordsPointsCap=new ArrayList<int[]>();
			capVals=(float[])imgMask.getStack().getProcessor(z+1).getPixels();
			nHits=0;
			for(int x=0;x<X;x++)for(int y=0;y<Y;y++) if(capVals[VitimageUtils.getPixelIndex(x, X, y)]>4.5) {nHits++; coordsPointsCap.add(new int[] {x,y});};

			for(int c=0;c<nImg+4;c++) {
				double[]vals=new double[nHits];
				incr=0;
				System.out.println("Going to do c="+c);
				VitimageUtils.printImageResume(tempRes[c]);
				capVals=(float[])tempRes[c].getStack().getProcessor(z+1).getPixels();				
				for(int i=0;i<nHits;i++) {
					vals[incr++]=capVals[VitimageUtils.getPixelIndex(coordsPointsCap.get(i)[0], X, coordsPointsCap.get(i)[1])];
				}
				double[]stats=VitimageUtils.statistics1D(vals);
				tempRes[c].getStack().setSliceLabel(tempRes[c].getStack().getSliceLabel(z+1)+"_CAP="+VitimageUtils.dou(stats[0])+"+-"+VitimageUtils.dou(stats[1]), z+1);
			}
		}	

		
		for(int c=0;c<tempRes.length;c++)tempRes[c]=VitimageUtils.convertFloatToShortWithoutDynamicChanges(tempRes[c]);			
		ImagePlus hyperImg=Concatenator.run(tempRes);
		hyperImg=HyperStackConverter.toHyperStack(hyperImg, tempRes.length,dims[2],1,"xyztc","Grayscale");

		hyperImg.setC(1);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoM0 );

		hyperImg.setC(2);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoT1 );

		hyperImg.setC(3);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoT2 );

		hyperImg.setC(4);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,5 );

		for(int c=0;c<nImg;c++) {
			hyperImg.setC(5+(computeMulti ? 0 : 0)+c);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoM0 );
		}

		hyperImg.setTitle("hypermap");
		return hyperImg;
	}
	*/

	
	
	
	
}