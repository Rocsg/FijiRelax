package com.vitimage.aplimtools;

import java.io.File;
import java.util.ArrayList;
import com.vitimage.common.Timer;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;
import com.vitimage.fijiyama.RegistrationAction;
import com.vitimage.mrutils.HyperMRIT1T2;
import com.vitimage.registration.BlockMatchingRegistration;
import com.vitimage.registration.ItkTransform;
import com.vitimage.registration.Transform3DType;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.FolderOpener;
import ij.plugin.HyperStackConverter;
import math3d.Point3d;

public class T1T2Seq_Importer{
	private static final long serialVersionUID = 1L;
	private String inputDir="";
	public String nameObservation="";
	private int test;
	private double[]voxs;
	private int[]dims;
	private int []valT1Trs;
	private int []valT2Tes;
	private boolean haveT2Seq;
	ImagePlus []imgT2Seq;
	ImagePlus []imgT1Seq;
	ImagePlus []imgT2SeqTest;
	ImagePlus []imgT1SeqTest;
	ImagePlus M0map;
	ImagePlus M0map2;
	ImagePlus T1map;
	ImagePlus T1map2;
	ImagePlus T2map;
	ImagePlus T2map2;
	public String exportFileName="";
	boolean forceNoMan=true;
	public boolean makeCrop=false;
	boolean skipRegForTesting=false;
	boolean viewRegistration=false;

	public boolean normalizeBeforeComputation=false;
	public boolean normalizeAfterComputation=false;

	public boolean forceNoBoutureTrick=false;
	public boolean makeSorghoTrick=false;
	public boolean makeBoutureTrick;
	public boolean forgetEarlyReps=false;
	public boolean speedupRegistrationForTesting=true;//false;
	
	public boolean hideMessage=true;
	public boolean dontShowNothing=false;
	public boolean computeMulti=true;
	
	private boolean haveT1Seq;
	private ImagePlus[][] imgT1T2;
	private ImagePlus[] imgT1T2Line;
	private ImagePlus M0mapMulti1;
	private ImagePlus M0mapMulti2;
	private ImagePlus T1mapMulti1;
	private ImagePlus T1mapMulti2;
	private ImagePlus T2mapMulti1;
	private ImagePlus T2mapMulti2;
	private ImagePlus imgMask;
	private int indexFirstT2Seq;
	private ImagePlus multiM0;
	private ImagePlus multiT2;
	private ImagePlus multiT1;
	private ImagePlus imgMaskCap;
	private static final int TrValueUsedForT2Seq=10000;

	public static void main (String[]args) {
		
	}
	
	
	public ImagePlus run(String inputDir,String name) {
		Timer t=new Timer();
		this.nameObservation=name;
			this.inputDir=inputDir;
		//Ask for the directory input, and the directory output
		IJ.log("Opening input and output dirs");
		
		//Read the voxel sizes and the image sizes		
		IJ.log("Gathering image informations");
		Object []objs=getFirstImageOfT1T2Serie(inputDir);
		ImagePlus testDcm=(ImagePlus)objs[0];
		this.voxs=VitimageUtils.getVoxelSizes(testDcm);

		t.print("Reading images");
		readT1T2();
		t.print("Starting registration");
		registerT1T2();
		t.print("Finishing registration");
		convert2DImgArrayTo1DArray();
		makeNullMaps();
		t.print("Finishing maps computation");
		IJ.log("Finished !");
		VitimageUtils.waitFor(10);
		return stackSaveAndShowResultsT1T2();
	}
		
	public void makeNullMaps(){
		multiM0=VitimageUtils.nullImage(imgT1T2Line[0]);
		multiT1=VitimageUtils.nullImage(imgT1T2Line[0]);
		T2map=VitimageUtils.nullImage(imgT1T2Line[0]);
	}
	
	public Object[] getFirstImageOfT1T2Serie(String inputDir) {
		String imgPath=inputDir;
		String[]subFile=null;
		while(!imgPath.substring(imgPath.length()-4, imgPath.length()).equals(".dcm")){
			IJ.log("Not yet an img :"+imgPath);
			subFile=VitimageUtils.stringArraySort(new File(imgPath).list());
			imgPath=new File(imgPath,subFile[0]).getAbsolutePath();
		}
		return new Object[] {IJ.openImage(imgPath),subFile.length};
	}
	
	public void readT1T2() {
		System.out.println(inputDir);
		String[]listFiles=VitimageUtils.stringArraySort(new File(inputDir).list());
		imgT1T2 =new ImagePlus[listFiles.length][];
		int nbSelected=listFiles.length;
		int indexStart=listFiles.length-nbSelected;
		for(int ii=indexStart;ii<listFiles.length;ii++) {
			int i=ii-indexStart;
			String[]strTes=VitimageUtils.stringArraySort(new File(inputDir,listFiles[ii]).list());
			imgT1T2[i]=new ImagePlus[strTes.length];
			for(int j=0;j<strTes.length;j++) {
				String str=new File(new File(inputDir,listFiles[ii]),strTes[j]).getAbsolutePath();
				System.out.println("Opening MRI Data : "+str);						
				imgT1T2[i][j]=FolderOpener.open(str, "");
				if(i==0 && j==0)detectSorgho(imgT1T2[i][j]);
				if(i==0 && j==0)detectBouture(imgT1T2[i][j]);
				int nbAverage=VitimageUtils.getAveraging(imgT1T2[i][j]);
				IJ.run(imgT1T2[i][j],"32-bit","");
				if(nbAverage!=2)imgT1T2[i][j]=VitimageUtils.makeOperationOnOneImage(imgT1T2[i][j], 2, 2.0/nbAverage, false);
				MRUtils.computeTeTrAndRiceSigmaOfEachSliceAndWriteItInTheLabels(imgT1T2[i][j],false,"T1T2SEQ_"+nameObservation,makeBoutureTrick);
			}
		}
		int ZT2=VitimageUtils.getDimensions(imgT1T2[imgT1T2.length-1][0])[2];
		System.out.println("ZT2="+ZT2);  
		double[]meanRiceTab=new double[ZT2];
		for(int z=0;z<ZT2;z++) {			meanRiceTab[z]=MRUtils.readSigmaInSliceLabel(imgT1T2[imgT1T2.length-1][0], 0, z, 0); }
		double meanRiceOverT2=VitimageUtils.statistics1D(meanRiceTab)[0];
		if(makeBoutureTrick)imgMask=VitimageUtils.getFloatBinaryMask(VitimageUtils.maxOfImageArrayDouble(imgT1T2),MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*meanRiceOverT2,1E10);
		else imgMask=VitimageUtils.getFloatBinaryMask(imgT1T2[imgT1T2.length-1][0],MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*meanRiceOverT2,1E10);
		imgMaskCap=HyperMRIT1T2.getAntiCapillaryMask(imgT1T2[imgT1T2.length-1][0]);
	}		
	
	
	public boolean detectSorgho(ImagePlus img) {
		String patientName=VitimageUtils.getPatientName(img);
		if( patientName.contains("BM") || patientName.contains("SSM")) {
			IJ.log("Detected experience : Sorgho\nSettings : "+
					"\nForget early reps = false (dismiss TR<=1000 ms)"+
					"\nNormalize before computation = true (correct uneven normalization issues before maps computation)"+
					"\nMake sorgho trick = true (set all images in a common geometry to help analysis)\n\n");
			return true;
		}
		return false;
	}

	public boolean detectBouture(ImagePlus img) {
		return VitimageUtils.isBouture(img);
	}


		
	public void registerT1T2(){		
		int nT1=this.imgT1T2.length;
		ItkTransform  []trs=new ItkTransform[nT1];

		//In case of sorgho data, set them all to the same voxel size
		if(makeSorghoTrick) {
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
			imgMaskCap=HyperMRIT1T2.getAntiCapillaryMask(imgT1T2[imgT1T2.length-1][0]);
		}
		imgMask=VitimageUtils.makeOperationBetweenTwoImages(imgMask,imgMaskCap, 1, true);
			
		
		trs[nT1-1]=new ItkTransform();
		ImagePlus imgMov,imgRef;
		ItkTransform tr=new ItkTransform();
		RegistrationAction regAct;
		BlockMatchingRegistration bmRegistration;
		imgMov=VitimageUtils.imageCopy(imgT1T2[nT1-2][0]);
		imgRef=VitimageUtils.imageCopy(imgT1T2[nT1-1][0]);
		
		//If needed, handle a manual registration before
		if(!forceNoMan) {
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
		for(int iMov=0;iMov<nT1-1;iMov++) {	
			if(skipRegForTesting) {
				trs[iMov]=new ItkTransform();
				continue;
			}
			imgMov=imgT1T2[iMov][0];
			imgRef=imgT1T2[nT1-1][0];
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
			bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
			bmRegistration.consoleOutputActivated=false;
			bmRegistration.timingMeasurement=true;
			bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
			bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
			bmRegistration.flagRange=true;
			bmRegistration.minBlockVariance=10;
			bmRegistration.displayRegistration=viewRegistration ? 2 : 0;
			bmRegistration.displayR2=false;
			bmRegistration.levelMin=-1;
			bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
			trs[iMov]=bmRegistration.runBlockMatching(tr);
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
				if((!forgetEarlyReps) && (MRUtils.readTrInSliceLabel(imgMov,0,0,0)<=400))regAct.levelMaxDense=1;
				regAct.typeTrans=Transform3DType.DENSE;
				bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
				bmRegistration.consoleOutputActivated=false;
				bmRegistration.timingMeasurement=true;
				bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
				bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
				bmRegistration.flagRange=true;
				bmRegistration.minBlockVariance=10;
				bmRegistration.displayRegistration=viewRegistration ? 2 : 0;
				bmRegistration.displayR2=false;
				bmRegistration.levelMin=-1;
				bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
				trs[iMov]=bmRegistration.runBlockMatching(trs[iMov]);
				bmRegistration.closeLastImages();
			}
			for(int j=0;j<imgT1T2[iMov].length;j++) {
				imgT1T2[iMov][j]=trs[iMov].transformImageExtensive(imgT1T2[nT1-1][0],imgT1T2[iMov][j]);
			}
		}
		IJ.log("Registration finished !");
	}
	
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
		if(makeSorghoTrick)nameObservation="SORGHO"+nameObservation;
		if(makeBoutureTrick)nameObservation="BOUTURE"+nameObservation;
		int nImg=imgT1T2Line.length;
		ImagePlus []tempRes=new ImagePlus[nImg+5];
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
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );

		hyperImg.setC(2);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT1 );

		hyperImg.setC(3);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );

		hyperImg.setC(4);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,5 );

		for(int c=0;c<nImg;c++) {
			hyperImg.setC(5+(computeMulti ? 0 : 0)+c);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );
		}

		hyperImg.setTitle("hypermap");
		return hyperImg;
	}

	
	
	
	
}
