package com.vitimage.fijirelax;

import java.io.File;
import java.util.ArrayList;
import com.vitimage.common.Timer;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;
import com.vitimage.fijirelax.HyperMap;
import com.vitimage.fijiyama.RegistrationAction;
import com.vitimage.registration.BlockMatchingRegistration;
import com.vitimage.registration.ItkTransform;
import com.vitimage.registration.Transform3DType;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.FolderOpener;
import ij.plugin.HyperStackConverter;
import math3d.Point3d;



public class Custom_Format_Importer{
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
	private String name;
//	public static void main (String[]args) {
		
	//}
	
	public Custom_Format_Importer(String inputDir,String name) {
		this.inputDir=inputDir;
		this.name=name;
	}
	
	
	public HyperMap importHyperMap() {
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
		return readT1T2();
		//t.print("Starting registration");
		//registerT1T2();
		//t.print("Finishing registration");
		//convert2DImgArrayTo1DArray();
		//makeNullMaps();
		//t.print("Finishing maps computation");
		//IJ.log("Finished !");
		//VitimageUtils.waitFor(10);
		//return new HyperMap(stackSaveAndShowResultsT1T2());
	}
	
	
	public void makeNullMaps(){
		multiM0=VitimageUtils.nullImage(imgT1T2Line[0]);
		multiT1=VitimageUtils.nullImage(imgT1T2Line[0]);
		T2map=VitimageUtils.nullImage(imgT1T2Line[0]);
	}
	
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
	
	
	public boolean detectSorgho(ImagePlus img) {
		String patientName=VitimageUtils.getPatientName(img);
		if( patientName.contains("BM") || patientName.contains("SSM")) {
			IJ.log("Detected experience : Sorgho\nSettings : "+
					"\nForget early reps = false (dismiss TR<=1000 ms)"+
					"\nNormalize before computation = true (correct uneven normalization issues before maps computation)"+
					"\nMake sorgho trick = true (set all images in a common geometry to help analysis)\n\n");
			makeSorghoTrick=true;
			return true;
		}
		return false;
	}

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


		
	public void registerT1T2(){		
		int nT1=this.imgT1T2.length;
		if(nT1==1)return;
		ItkTransform  []trs=new ItkTransform[nT1];
		VitimageUtils.printImageResume(imgMask,"imgMask2");
		VitimageUtils.printImageResume(imgMaskCap,"imgMaskCap2");

		//In case of sorgho data, set them all to the same voxel size
		if(makeSorghoTrick) {
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
			bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
			bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
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
				bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
				bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
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

	
	
	
	
}
