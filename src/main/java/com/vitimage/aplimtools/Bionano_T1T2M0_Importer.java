package com.vitimage.aplimtools;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import com.vitimage.common.Timer;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;
import com.vitimage.fijiyama.RegistrationAction;
import com.vitimage.fijiyama.RegistrationManager;
import com.vitimage.mrutils.HyperMRIT1T2;
import com.vitimage.mrutils.MRUtils;
import com.vitimage.registration.BlockMatchingRegistration;
import com.vitimage.registration.ItkTransform;
import com.vitimage.registration.MetricType;
import com.vitimage.registration.Transform3DType;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.FolderOpener;
import ij.plugin.HyperStackConverter;
import ij.plugin.frame.PlugInFrame;
import math3d.Point3d;

public class Bionano_T1T2M0_Importer extends PlugInFrame{
	private static final long serialVersionUID = 1L;
	private String inputDir="";
	private String outputDir="";
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
	boolean forceNoMan=true;
	public boolean makeCrop=false;
	boolean skipRegForTesting=true;
	boolean viewRegistration=true;

	private boolean normalizeBeforeComputation=false;
	private boolean normalizeAfterComputation=false;

	private boolean forceNoBoutureTrick=false;
	boolean makeSorghoTrick=false;
	private boolean forgetEarlyReps=false;//true   
	private boolean speedupRegistrationForTesting=true;//false
	
	public boolean hideMessage=false;
	private boolean dontShowNothing=false;
	private boolean computeMulti=false;
	
	private boolean haveT1Seq;
	private ImagePlus[][] imgT1T2;
	private ImagePlus[] imgT1T2Line;
//	private double[][] trT1T2;
//	private double[][] teT1T2;
	private String lastDir;
	private ImagePlus M0mapMulti1;
	private ImagePlus M0mapMulti2;
	private ImagePlus T1mapMulti1;
	private ImagePlus T1mapMulti2;
	private ImagePlus T2mapMulti1;
	private ImagePlus T2mapMulti2;
	private ImagePlus imgMask;
	private int indexFirstT2Seq;
	private boolean makeBoutureTrick;
	private static final int TrValueUsedForT2Seq=10000;
	
	public Bionano_T1T2M0_Importer() {
		super("");
	}

	public static void main(String[] args) {
		@SuppressWarnings("unused")
		ImageJ ij=new ImageJ();		
		//runChainSorgho();
		
		ItkTransform tr=new ItkTransform();
		IJ.log("Running T1-T2-M0 import");
		Bionano_T1T2M0_Importer bn=new Bionano_T1T2M0_Importer();
		bn.run(new Object[] {"/home/fernandr/Bureau/Recherche_diff/B099J0/SL_RAW/","/home/fernandr/Bureau/Recherche_diff/B099J0/Imported2"});
	//	bn.run(new Object[] {null});
	}

	
	public void detectSorgho(ImagePlus img) {
		String patientName=VitimageUtils.getPatientName(img);
		if( patientName.contains("BM") || patientName.contains("SSM")) {
			forgetEarlyReps=true;
			normalizeBeforeComputation=true;
			makeSorghoTrick=true;
			System.out.println("Detected experience : Sorgho\nSettings : "+
					"\nForget early reps = true (dismiss TR<=1000 ms)"+
					"\nNormalize before computation = true (correct uneven normalization issues before maps computation)"+
					"\nMake sorgho trick = true (set all images in a common geometry to help analysis)\n\n");
		}
	}
	

	public void detectBouture(ImagePlus img) {
		System.out.println("STARTING DETECTION");
		System.out.println("1 "+makeBoutureTrick+" , "+forceNoBoutureTrick);
		makeBoutureTrick=VitimageUtils.isBouture(img);
		System.out.println("2 "+makeBoutureTrick+" , "+forceNoBoutureTrick);
		if(makeBoutureTrick==false)return;
		System.out.println("3 "+makeBoutureTrick+" , "+forceNoBoutureTrick);
		if(forceNoBoutureTrick)makeBoutureTrick=false;
		System.out.println("4 "+makeBoutureTrick+" , "+forceNoBoutureTrick);
		System.out.println("Detected experience : Bouture\nMake bouture trick activated ? "+makeBoutureTrick);
		System.out.println("5 "+makeBoutureTrick+" , "+forceNoBoutureTrick);
	}
	
	
	
	public static void runChainSorgho() {
		String sourceDir="/home/fernandr/Bureau/Traitements/Sorgho/Donnees_brutes_export_Romain";
		String[]specs=new String[] {"BM1","BM2","SSM1","SSM2",};
		int count=0;
		for(String sp : specs) {
			//if(! sp.equals("SSM1"))continue;
			String spDir=new File(sourceDir,sp).getAbsolutePath();
			System.out.println("Traitement des series de l'individu "+sp);
			String []strTimepoints=new File(spDir).list();
			for(String tim : strTimepoints) {
				count++;
				if(count<17)continue;
				//if(! tim.equals("SSM1_F1_0923_SL_RAW"))continue;
				System.out.println("Traitement des series de l'individu "+sp+", observation "+tim);
				Bionano_T1T2M0_Importer bn=new Bionano_T1T2M0_Importer();
				bn.viewRegistration=false;
				bn.hideMessage=true;
				bn.dontShowNothing=true;
				bn.run(new Object[] {new File(spDir,tim).getAbsolutePath(),"/home/fernandr/Bureau/Traitements/Sorgho/Cartes_calculees_methode_Romain"});
			}
		}
	}
	
	
	public void run(String arg) {
		IJ.log("Starting MRI Importer");
		//run(new Object[] {"/home/fernandr/Bureau/Traitements/Sorgho/export/SSM1/SSM1_F1_0923_SL_RAW","/home/fernandr/Bureau/TestPouet"});
		run(new Object[] {});
	}
	
	public void run(Object[]obj) {
		//Ask for the directory input, and the directory output
		IJ.log("Opening input and output dirs");
		if(obj != null && obj.length==2) {
			inputDir=(String)obj[0];
			outputDir=(String)obj[1];
		}
		else{
			inputDir=VitiDialogs.chooseDirectoryUI("Localize inputDir (SL_RAW)", "Select this directory");
			outputDir=VitiDialogs.chooseDirectoryUI("Select outputDir", "Select this directory");
		}

		lastDir=new File(inputDir).getName();
		
		//Read the voxel sizes and the image sizes		
		IJ.log("Gathering image informations");
		Object []objs=getFirstImageOfT1T2Serie(inputDir);
		ImagePlus testDcm=(ImagePlus)objs[0];
		int nSlices=(int)objs[1];
		this.voxs=VitimageUtils.getVoxelSizes(testDcm);

		readT1T2();
		ImagePlus babyShower1=displayT1T2Sequence(" original position ");
		System.out.println();
		Timer t=new Timer();
		t.print("Starting registration");
		registerT1T2();
		t.print("Finishing registration");
		ImagePlus babyShower2=displayT1T2Sequence(" registered space");
		ImagePlus babyShower4=null;
		if(normalizeBeforeComputation) {
			normalizeBeforeComputation();
			babyShower4=displayT1T2Sequence(" registered space and normalized");
		}
		convert2DImgArrayTo1DArray();
		t=new Timer();
		t.print("Starting maps computation");
		computeT1T2(new Object[] {(makeSorghoTrick ? 0.75 : 0.75)});
//		computeT1T2Multi(new Object[] {5.0});
		t.print("Finishing maps computation");
		if(babyShower1!=null)babyShower1.close();
		if(babyShower2!=null)babyShower2.close();
		if(normalizeBeforeComputation && babyShower4!=null)babyShower4.close();
		ImagePlus babyShower5=stackSaveAndShowResultsT1T2();
		IJ.log("Finished !");
		VitimageUtils.waitFor(10);
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
	
	public ImagePlus displayT1T2Sequence(String type) {
		int n=0;int incr=0;
		for(int i=0;i<imgT1T2.length;i++)n+=imgT1T2[i].length;
		ImagePlus[]imgsTemp=new ImagePlus[n];
		for(int i=0;i<imgT1T2.length;i++)for(int j=0;j<imgT1T2[i].length;j++) {
			imgsTemp[incr++]=imgT1T2[i][j].duplicate();
		}
		int nZ=imgsTemp[0].getNSlices();
		int nF=imgsTemp[0].getNFrames();
		ImagePlus imgNew=Concatenator.run(imgsTemp);
		ImagePlus test=HyperStackConverter.toHyperStack(imgNew, n,nZ,nF,"xyztc","Fire");
		for(int c=1;c<=test.getNChannels();c++) {test.setC(c);IJ.run(test,"Fire","");test.setDisplayRange(0, MRUtils.maxDisplayedM0);}
		test.setC(1);
		test.setTitle("Stacked T1 T2 in "+type);
		if(!dontShowNothing)test.show();
		return test;
	}

	public void readT1T2() {
		String[]listFiles=VitimageUtils.stringArraySort(new File(inputDir).list());
		int nbSelected=0;
		if(forgetEarlyReps) {		
			for(int i=0;i<listFiles.length;i++) {int a=Integer.parseInt(listFiles[i].replace("TR","")); if (a>1000)nbSelected++;}
			imgT1T2 =new ImagePlus[nbSelected][];
//			trT1T2 =new double[nbSelected][];
//			teT1T2 =new double[nbSelected][];
		}
		else{ 
			imgT1T2 =new ImagePlus[listFiles.length][];
//			trT1T2 =new double[listFiles.length][];
//			teT1T2 =new double[listFiles.length][];
			nbSelected=listFiles.length;
		}
		int indexStart=listFiles.length-nbSelected;
		for(int ii=indexStart;ii<listFiles.length;ii++) {
			int i=ii-indexStart;
			String[]strTes=VitimageUtils.stringArraySort(new File(inputDir,listFiles[ii]).list());
			imgT1T2[i]=new ImagePlus[strTes.length];
//			trT1T2[i]=new double[strTes.length];
//			teT1T2[i]=new double[strTes.length];
			for(int j=0;j<strTes.length;j++) {
				String str=new File(new File(inputDir,listFiles[ii]),strTes[j]).getAbsolutePath();
				System.out.println("Opening MRI Data : "+str);						
				imgT1T2[i][j]=FolderOpener.open(str, "");
				if(i==0 && j==0)detectSorgho(imgT1T2[i][j]);
				if(i==0 && j==0)detectBouture(imgT1T2[i][j]);
				int nbAverage=VitimageUtils.getAveraging(imgT1T2[i][j]);
				IJ.run(imgT1T2[i][j],"32-bit","");
				if(nbAverage!=2)imgT1T2[i][j]=VitimageUtils.makeOperationOnOneImage(imgT1T2[i][j], 2, 2.0/nbAverage, false);
				MRUtils.computeTeTrAndRiceSigmaOfEachSliceAndWriteItInTheLabels(imgT1T2[i][j],false,"T1T2SEQ",makeBoutureTrick);
			}
		}
		int ZT2=VitimageUtils.getDimensions(imgT1T2[imgT1T2.length-1][0])[2];
		System.out.println("ZT2="+ZT2);
		double[]meanRiceTab=new double[ZT2];
		for(int z=0;z<ZT2;z++) {			meanRiceTab[z]=MRUtils.readSigmaInSliceLabel(imgT1T2[imgT1T2.length-1][0], 0, z, 0); }
		double meanRiceOverT2=VitimageUtils.statistics1D(meanRiceTab)[0];
		imgMask=VitimageUtils.getFloatBinaryMask(imgT1T2[imgT1T2.length-1][0],MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*meanRiceOverT2,1E10);
	}		
	
	public void registerT1T2(){		
		int nT1=this.imgT1T2.length;
		ItkTransform  []trs=new ItkTransform[nT1];

		if(makeSorghoTrick) {
			double targetVZ=0.5;
			double targetVoxVolume=1E-3;
			double targetVX=Math.sqrt(targetVoxVolume/targetVZ);
			double targetVY=Math.sqrt(targetVoxVolume/targetVZ);
			int[]targetDims=new int[] {512,512,4};
			double[]targetVoxs=new double[] {targetVX,targetVY,targetVZ};
			
			for(int j=0;j<imgT1T2[nT1-1].length;j++) {
//				System.out.println("Making TRANSFORM T2");
//				System.out.println("Before : "+imgT1T2[nT1-1][j].getStack().getSliceLabel(1));
				imgT1T2[nT1-1][j]=new ItkTransform().transformImage(targetDims, targetVoxs, imgT1T2[nT1-1][j],false);
//				System.out.println("After : "+imgT1T2[nT1-1][j].getStack().getSliceLabel(1));
			}
			imgMask.getStack().setSliceLabel("", 1);
			imgMask=new ItkTransform().transformImage(targetDims, targetVoxs, imgMask,false);
			imgMask=VitimageUtils.getFloatBinaryMask(imgMask, 1.0, 1E10);
		}
		
		
		
		
		trs[nT1-1]=new ItkTransform();
		ImagePlus imgMov,imgRef;
		ItkTransform tr=new ItkTransform();
		RegistrationAction regAct;
		BlockMatchingRegistration bmRegistration;
		imgMov=VitimageUtils.imageCopy(imgT1T2[nT1-2][0]);
		imgRef=VitimageUtils.imageCopy(imgT1T2[nT1-1][0]);

		
		
		if(!forceNoMan) {
			ImagePlus comp=VitimageUtils.compositeOf(imgT1T2[nT1-1][0],imgT1T2[nT1-2][0],"red=T2 first echo , green = T1 last rep time");
			if(!dontShowNothing)comp.show();
			if(VitiDialogs.getYesNoUI("Manual registration first ?", "Manual registration first ?\nNeeded if T1 stack and T2 stack are largely disaligned")) {
				Point3d[][]pts=VitiDialogs.registrationPointsUI(5,imgRef,imgMov,true);
				tr=ItkTransform.estimateBestRigid3D(pts[1],pts[0]); 
			}
			if(!dontShowNothing)comp.close();
		}
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
			if(speedupRegistrationForTesting) {
				regAct.strideX=regAct.strideX*2;
				regAct.strideY=regAct.strideY*2;
				regAct.strideZ=regAct.strideZ*2;
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
			if((!forgetEarlyReps) && (iMov<1))regAct.levelMaxDense=1;
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
			for(int j=0;j<imgT1T2[iMov].length;j++) {
				imgT1T2[iMov][j]=trs[iMov].transformImage(imgT1T2[nT1-1][0],imgT1T2[iMov][j]);
			}
		}
		IJ.log("Registration finished !");
	}



	
	public void normalizeBeforeComputation() {
		for(int i=0;i<imgT1T2.length;i++) {
			imgT1T2[i]=MRUtils.normalizeBeforeComputation(imgT1T2[i],imgMask);
		}
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
	}
	
//	for(int i=0;i<imgT1T2.length;i++)		for(int j=0;j<imgT1T2[i].length;j++)imgT1T2[i][j]=VitimageUtils.makeOperationBetweenTwoImages(imgT1T2[i][j], imgMask, 2, true);

	
	
	public void computeT1T2(Object[]params) {
		double sigmaSmoothing=0;
		if(params!=null) sigmaSmoothing=(double)params[0];
		sigmaSmoothing=0.75;
		ImagePlus[] maps=null;
		//if(makeBoutureTrick)MRUtils.computeT1T2MapBoutureTrick(imgT1T2Line, sigmaSmoothing,MRUtils.T1T2_MONO_RICE);
		//else
		maps=MRUtils.computeT1T2Map(imgT1T2Line, sigmaSmoothing,MRUtils.T1T2_MONO_RICE);
		
		ImagePlus antiMask=VitimageUtils.makeOperationOnOneImage(imgMask, 2, -1, true);
		antiMask=VitimageUtils.makeOperationOnOneImage(antiMask, 1, 1,true);

		ImagePlus tmp1=VitimageUtils.makeOperationBetweenTwoImages(imgMask, maps[0], 2, true);
		ImagePlus tmp2=VitimageUtils.makeOperationBetweenTwoImages(antiMask, imgT1T2Line[indexFirstT2Seq], 2, true);
		M0map=VitimageUtils.makeOperationBetweenTwoImages(tmp1, tmp2, 1, true);
		T1map=VitimageUtils.makeOperationBetweenTwoImages(maps[1],imgMask,2,true);
		T2map=VitimageUtils.makeOperationBetweenTwoImages(maps[2],imgMask,2,true);
	}
		
	public void computeT1T2Multi(Object[]params) {
		if(!computeMulti)return;
		double sigmaSmoothing=0.75;
		if(params!=null) sigmaSmoothing=(double)params[0];
		ImagePlus[] maps=MRUtils.computeT1T2Map(imgT1T2Line, sigmaSmoothing,MRUtils.T1T2_MULTI_RICE);
		
		M0mapMulti1=maps[0];
		M0mapMulti2=maps[1];
		T1mapMulti1=maps[2];
		T1mapMulti2=maps[3];
		T2mapMulti1=maps[4];
		T2mapMulti2=maps[5];
	}
		

	
	//TODO : si T1 et T2, verifier que le tout matche en terme de XYZ
	public ImagePlus stackSaveAndShowResultsT1T2() {
		dims=VitimageUtils.getDimensions(imgT1T2Line[0]);
		int nImg=imgT1T2Line.length;
		ImagePlus []tempRes=new ImagePlus[nImg+3+(computeMulti ? 6:0)];
		int curInd=0;
		int nMaps=3+(computeMulti ? 6:0);
		double meanM0Cap=HyperMRIT1T2.measureMeanCapillaryValueAlongZ(M0map.duplicate());
		System.out.println("\nMean M0 measured on the capillary : "+meanM0Cap);
		System.out.println("Target mean M0 : "+MRUtils.maxM0ForNormalization);
		System.out.println("Mult factor : "+MRUtils.maxM0ForNormalization/meanM0Cap);
		tempRes[0]=normalizeAfterComputation ? VitimageUtils.makeOperationOnOneImage(M0map, 2, MRUtils.maxM0ForNormalization/meanM0Cap, true) : M0map;
		tempRes[1]=T1map.duplicate();
		tempRes[2]=T2map.duplicate();
		tempRes[0].setDisplayRange(0,MRUtils.maxDisplayedM0 );
		tempRes[1].setDisplayRange(0,MRUtils.maxDisplayedT1 );
		tempRes[2].setDisplayRange(0,MRUtils.maxDisplayedT2 );
		if(computeMulti) {
			tempRes[3]=M0mapMulti1.duplicate();
			tempRes[4]=M0mapMulti2.duplicate();
			tempRes[5]=T1mapMulti1.duplicate();
			tempRes[6]=T1mapMulti2.duplicate();
			tempRes[7]=T2mapMulti1.duplicate();
			tempRes[8]=T2mapMulti2.duplicate();
			tempRes[3].setDisplayRange(0,MRUtils.maxDisplayedM0 );
			tempRes[4].setDisplayRange(0,MRUtils.maxDisplayedM0 );
			tempRes[5].setDisplayRange(0,MRUtils.maxDisplayedT1 );
			tempRes[6].setDisplayRange(0,MRUtils.maxDisplayedT1 );
			tempRes[7].setDisplayRange(0,MRUtils.maxDisplayedT2 );
			tempRes[8].setDisplayRange(0,MRUtils.maxDisplayedT2 );
		}
		for(int z=1;z<=dims[2];z++) {
			tempRes[0].getStack().setSliceLabel("M0MAP", z);
			tempRes[1].getStack().setSliceLabel("T1MAP", z);
			tempRes[2].getStack().setSliceLabel("T2MAP", z);
			if(computeMulti) {
				tempRes[3].getStack().setSliceLabel("M0MAP_M1", z);
				tempRes[4].getStack().setSliceLabel("M0MAP_M2", z);
				tempRes[5].getStack().setSliceLabel("T1MAP_M1", z);
				tempRes[6].getStack().setSliceLabel("T1MAP_M2", z);
				tempRes[7].getStack().setSliceLabel("T2MAP_M1", z);
				tempRes[8].getStack().setSliceLabel("T2MAP_M2", z);
			}
		}
		for(int c=0;c<nImg;c++)tempRes[3+(computeMulti ? 6 : 0)+c]=
				normalizeAfterComputation ? VitimageUtils.makeOperationOnOneImage(imgT1T2Line[c], 2, MRUtils.maxM0ForNormalization/meanM0Cap, true) : imgT1T2Line[c];

		ImagePlus hyperImg=Concatenator.run(tempRes);
		hyperImg=HyperStackConverter.toHyperStack(hyperImg, nImg+3+(computeMulti ? 6 : 0),dims[2],1,"xyztc","Grayscale");

		hyperImg.setC(1);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );

		hyperImg.setC(2);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT1 );

		hyperImg.setC(3);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );

		if(computeMulti) {
			hyperImg.setC(4);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );
			hyperImg.setC(5);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );

			hyperImg.setC(6);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT1 );
			hyperImg.setC(7);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT1 );

			hyperImg.setC(8);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );
			hyperImg.setC(9);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );

		
		}
		for(int c=0;c<nImg;c++) {
			hyperImg.setC(4+(computeMulti ? 6 : 0)+c);
			IJ.run(hyperImg,"Fire","");			
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );
		}

		hyperImg.setTitle("Hyperimage");
		String name=new File(outputDir,lastDir+"_hyperimage_build_time_"+(new SimpleDateFormat("yyyy-MM-dd_hh-mm").format(new Date())+".tif")).getAbsolutePath();
		IJ.saveAsTiff(hyperImg,name);
		hyperImg.close();
		
		ImagePlus result=IJ.openImage(name);
		if(!hideMessage) {
			result.show();
			IJ.showMessage("Hyper MRI image saved as : "+name);
		}
		return result;
	}



	
	
	
	
	
	
	
	
	
	
	
	
	
/*
	//////////////////////OLD CODES //////////////////////
	if(new File(inputDir,"TR010000").exists()) {
		haveT2Seq=true;
	}

	//Look for a T1 sequence, take the successives Tr
	if(new File(inputDir).list().length>2) {
		haveT1Seq=true;
	}

	if(hasT1Seq
	
	readT1Reps();
	displayT1Sequence();

	readT2Echoes();
		displayT2Sequence();

		//Register T1 seq if needed
		IJ.log("Registering T1 sequence");
		if(haveT1Seq && (true || (debug || VitiDialogs.getYesNoUI("Register T1 repetition images ?","Register T1 repetition images ??")))){
			registerT1Reps();
		}
		displayT1Sequence();
			

		//Show T1 estimation and T1 M0 for multiple sigma, then compute T1 and M0
		IJ.log("Computing T1 map");
		computeT1(new Object[] {1.0});
		T1map.show();
		T1map.setDisplayRange(0, MRUtils.maxDisplayedT1);
		T1map.setTitle("T1 map");


		
		//Show T2 estimation and T2 M0 for multiple sigma, then compute T2 and M0
		IJ.log("Computing T2 map");
		computeT2(new Object[] {1.0});
		T2map.show();
		T2map.setDisplayRange(0, MRUtils.maxDisplayedT2);
		M0map.show();
		M0map.setDisplayRange(0, MRUtils.maxDisplayedM0);
		T2map.setTitle("T2 map");
		M0map.setTitle("M0 map");
		
		
		//Stack M0, T1, T2, T1 reps, T2 echoes and save the result
		IJ.log("Stacking results");
		stackSaveAndShowResults();
	}
	
	else {

	*/
	
	public void displayT2Sequence() {
		ImagePlus[]imgs=new ImagePlus[imgT2Seq.length];int nF=imgs.length;int nZ=imgT2Seq[0].getNSlices();
		for(int i=0;i<imgT2Seq.length;i++) imgs[i]=VitimageUtils.imageCopy(imgT2Seq[i]);	

		ImagePlus test=HyperStackConverter.toHyperStack(Concatenator.run(imgs), 1,nZ,nF,"xyztc","Fire");
		test.setSlice(dims[2]/2);
		test.setTitle("Stacked T2 echoes in original position");
		test.resetDisplayRange();
		test.show();
	}
	
	public void displayT1Sequence() {
		ImagePlus[]imgs=new ImagePlus[imgT1Seq.length];int nF=imgs.length;int nZ=imgT1Seq[0].getNSlices();
		for(int i=0;i<imgT1Seq.length;i++) {
			imgs[i]=VitimageUtils.imageCopy(imgT1Seq[i]);
		}
		ImagePlus imgNew=Concatenator.run(imgs);
		ImagePlus test=HyperStackConverter.toHyperStack(imgNew, 1,nZ,nF,"xyztc","Fire");
		test.setTitle("Stacked T1 rep times in original position");
		test.setSlice(test.getStackSize()-1-dims[2]/2);
		test.resetDisplayRange();
		test.show();
	}


	
	public void readAlreadyImportedData() {
		int n=52;
		imgT1T2Line=new ImagePlus[n];
		for(int i=0;i<n;i++) {
			System.out.println(i);
			imgT1T2Line[i]=IJ.openImage(new File(outputDir,"img_"+(i)+".tif").getAbsolutePath());
		}
		
	}
	

	
	public void readT2Echoes() {
		String[]listFiles=VitimageUtils.stringArraySort(new File(inputDir,"TR010000").list());
		this.valT2Tes=new int[listFiles.length];
		this.imgT2Seq=new ImagePlus[listFiles.length];
		this.imgT2SeqTest=new ImagePlus[listFiles.length];
		for(int i=0;i<this.valT2Tes.length;i++) {
			String str=new File(new File(inputDir,"TR0"+TrValueUsedForT2Seq),listFiles[i]).getAbsolutePath();//TODO : genericity among TRvalue : TrValueUsedForT2Seq is fixed
			imgT2Seq[i] = FolderOpener.open(str, "");
			int nbAverage=VitimageUtils.getAveraging(imgT2Seq[i]);
			this.valT2Tes[i]=Integer.parseInt(listFiles[i].replace("TE",""));
			IJ.run(imgT2Seq[i],"32-bit","");			
			if(nbAverage!=2)imgT2Seq[i]=VitimageUtils.makeOperationOnOneImage(imgT2Seq[i], 2, 2.0/nbAverage, false);
			MRUtils.computeTeTrAndRiceSigmaOfEachSliceAndWriteItInTheLabels(imgT2Seq[i],false,"T2SEQ",makeBoutureTrick);
		}
	}		
			
	public void readT1Reps() {
		String[]listFiles=VitimageUtils.stringArraySort(new File(inputDir).list());
		String[]strTes=VitimageUtils.stringArraySort(new File(inputDir,listFiles[0]).list());
		String strChosenTe=strTes[0];
		int teValue=Integer.parseInt(strChosenTe.replace("TE",""));
		this.valT1Trs=new int[listFiles.length-(haveT2Seq ? 1 : 0)];
		this.imgT1Seq=new ImagePlus[this.valT1Trs.length];
		this.imgT1SeqTest=new ImagePlus[this.valT1Trs.length];
		for(int i=0;i<this.valT1Trs.length;i++) {
			String str=new File(new File(inputDir,listFiles[i]),strChosenTe).getAbsolutePath();
			imgT1Seq[i] = FolderOpener.open(str, "");
			int nbAverage=VitimageUtils.getAveraging(imgT1Seq[i]);
			this.valT1Trs[i]=Integer.parseInt(listFiles[i].replace("TR",""));
			IJ.run(imgT1Seq[i],"32-bit","");
			if(nbAverage!=2)imgT1Seq[i]=VitimageUtils.makeOperationOnOneImage(imgT1Seq[i], 2, 2.0/nbAverage, false);
			MRUtils.computeTeTrAndRiceSigmaOfEachSliceAndWriteItInTheLabels(imgT1Seq[i],false,"T1SEQ",makeBoutureTrick);
		}
	}		


	
	public void computeT1(Object[]params) {
		double sigmaSmoothing=0;
		if(params!=null) sigmaSmoothing=(double)params[0];
		else {
			double[]successiveSigmasSmoothing=new double[] {1};//{0.3,0.5,0.67,1.0,1.5,2.0,3.0,5.0};
			int nSucSig=successiveSigmasSmoothing.length;
			ImagePlus []resultsTab=new ImagePlus[nSucSig*2];
			for(int i=0;i<nSucSig;i++) {
				double valSig=successiveSigmasSmoothing[i];
				ImagePlus[]mapsComputed=MRUtils.computeT1Map(imgT1SeqTest, valSig);
				resultsTab[i*2]=mapsComputed[0];
				resultsTab[i*2].getStack().setSliceLabel("M0 map from T1 seq\nWith sigmaSmooth="+valSig,1);
				resultsTab[i*2+1]=mapsComputed[1];
				resultsTab[i*2+1].getStack().setSliceLabel("T1 map from T1 seq\nWith sigmaSmooth="+valSig,1);
			}
			ImagePlus test=HyperStackConverter.toHyperStack(new Concatenator().concatenate(resultsTab,false), nSucSig, 1,2,"xyztc","Fire");
			test.setTitle("T1 map test");
			test.show();
			int count=30;
			while((count--) >0 && test.isVisible()){
				test.setTitle("T1 seq results, "+count+" seconds before closing");
				VitimageUtils.waitFor(1000);
			}
			test.close();
			sigmaSmoothing=VitiDialogs.getDoubleUI("Choose sigma for smoothing", "",1);
		}
		T1map=MRUtils.computeT1Map(imgT1Seq, sigmaSmoothing)[1];
	}
	
	public void computeT2(Object[]params) {
		double sigmaSmoothing=0;
		if(params!=null) sigmaSmoothing=(double)params[0];
		else {
			double[]successiveSigmasSmoothing=new double[] {1};//{0.3,0.5,0.67,1.0,1.5,2.0,3.0,5.0};
			int nSucSig=successiveSigmasSmoothing.length;
			ImagePlus []resultsTab=new ImagePlus[nSucSig*2];
			for(int i=0;i<nSucSig;i++) {
				double valSig=successiveSigmasSmoothing[i];
				ImagePlus[]mapsComputed=MRUtils.computeT2Map(imgT2SeqTest, valSig);
				resultsTab[i*2]=VitimageUtils.writeWhiteOnImage("M0 map from T2 seq\nWith sigmaSmooth="+valSig, mapsComputed[0], 7,0);
				resultsTab[i*2+1]=VitimageUtils.writeWhiteOnImage("T2 map from T2 seq\nWith sigmaSmooth="+valSig, mapsComputed[1], 7,0);
			}
			ImagePlus test=HyperStackConverter.toHyperStack(new Concatenator().concatenate(resultsTab,false), nSucSig, 1,2,"xyztc","Fire");
			test.show();
			int count=30;
			while((count--) >0 && test.isVisible()){
				test.setTitle("T1 seq results, "+count+" seconds before closing");
				VitimageUtils.waitFor(1000);
			}
			test.close();
			sigmaSmoothing=VitiDialogs.getDoubleUI("Choose sigma for smoothing", "",1);
		}
		ImagePlus []results=null;
		results= MRUtils.computeT2Map(imgT2Seq, sigmaSmoothing);
		T2map=results[1];
		M0map=results[0];
	}
	

	
	
	
	public void registerT1Reps(){		
		boolean crossFit=false;		
		int nT1=this.imgT1Seq.length;
		ItkTransform  []trs=new ItkTransform[nT1];
		trs[nT1-1]=new ItkTransform();
		ImagePlus imgMov,imgRef;
		ItkTransform tr=new ItkTransform();
		RegistrationAction regAct;
		BlockMatchingRegistration bmRegistration;

		//Registration of last repetition time of T1 seq (mov) onto the first echo time of T2 seq (ref)
		if(haveT1Seq && haveT2Seq) {
			ImagePlus comp=VitimageUtils.compositeOf(imgT2Seq[0],imgT1Seq[nT1-1],"red=T2 first echo , green = T1 last rep time");
			comp.show();
			if(true ||  VitiDialogs.getYesNoUI("Register T1 and T2 rigid ?", "Register T1 and T2 rigid ?")){
				comp.close();
				imgMov=VitimageUtils.imageCopy(imgT1Seq[nT1-1]);
				imgRef=VitimageUtils.imageCopy(imgT2Seq[0]);
				if(VitiDialogs.getYesNoUI("Manual registration first ?", "Manual registration first ?\nNeeded if T1 stack and T2 stack are largely disaligned")) {
					Point3d[][]pts=VitiDialogs.registrationPointsUI(5,imgRef,imgMov,true);
					tr=ItkTransform.estimateBestRigid3D(pts[1],pts[0]); 
				}
				else tr=new ItkTransform();
				imgRef.show();
				imgMov.show();
				regAct=new RegistrationAction();
				regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
				regAct.higherAcc=1;
				bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
				bmRegistration.consoleOutputActivated=false;
				bmRegistration.timingMeasurement=true;
				bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
				bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
				bmRegistration.flagRange=true;
				bmRegistration.minBlockVariance=10;
				bmRegistration.displayRegistration=2;
				bmRegistration.displayR2=false;
				bmRegistration.levelMin=-1;
				bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
				tr=bmRegistration.runBlockMatching(tr);
				bmRegistration.closeLastImages();

				//Ask for using the TR10000 to make the T1 fit. In this case, the registration between both should be really accurate, implying potentially
				ImagePlus temp=tr.transformImage( imgT2Seq[0], imgT1Seq[nT1-1], false);
				comp=VitimageUtils.compositeOf(imgT2Seq[0],temp,"red=T2 first echo , green = T1 last rep time");
				comp.show();
				if(true || VitiDialogs.getYesNoUI("Register T1 and T2 dense ? ", "Register T1 and T2 dense ? \n(to obtain a very high accuracy in positioning)")) {
					comp.close();
					regAct=new RegistrationAction();
					regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
					regAct.higherAcc=1;
					regAct.typeTrans=Transform3DType.DENSE;
					bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
					bmRegistration.consoleOutputActivated=false;
					bmRegistration.timingMeasurement=true;
					bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
					bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
					bmRegistration.flagRange=true;
					bmRegistration.minBlockVariance=10;
					bmRegistration.displayRegistration=2;
					bmRegistration.displayR2=false;
					bmRegistration.levelMin=-1;
					bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
					tr=bmRegistration.runBlockMatching(tr);
					bmRegistration.closeLastImages();
				}
			}
			if(comp.isVisible())comp.close();
			ImagePlus temp=tr.transformImage( imgT2Seq[0], imgT1Seq[nT1-1], false);
			comp=VitimageUtils.compositeOf(imgT2Seq[0],temp,"red=T2 first echo , green = T1 last rep time");
			comp.show();
			if(VitiDialogs.getYesNoUI("Use the first T2 echo for T1 seq computation ?", "Use the TR10000 of T2 seq as a T1 seq point")) crossFit=true;	
			comp.close();
		}
	

		
		
		//Registration of T1 successive images, if needed
		ImagePlus comp=VitimageUtils.compositeOf(imgT1Seq[0],imgT1Seq[imgT1Seq.length-1],"T1 seq (red=first rep time , green = last rep time)");
		comp.show();
		if(true || VitiDialogs.getYesNoUI("Register T1 images before T1 map estimation ? ", "Register T1 images before T1 map estimation ? (rigid reg.)")) {
			comp.close();
			imgRef=imgT1Seq[nT1-1];
			for(int iMov=0;iMov<nT1;iMov++) {
				if(iMov<nT1-1) {
					imgMov=imgT1Seq[iMov];
					regAct=new RegistrationAction();
					regAct.defineSettingsSimplyFromTwoImages(imgRef,imgMov);
					regAct.higherAcc=1;
					regAct.typeTrans=Transform3DType.EULER2D;
					bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov,regAct);
					bmRegistration.consoleOutputActivated=false;
					bmRegistration.timingMeasurement=true;
					bmRegistration.refRange=new double[] {imgRef.getDisplayRangeMin(),imgRef.getDisplayRangeMax()};
					bmRegistration.movRange=new double[] {imgMov.getDisplayRangeMin(),imgMov.getDisplayRangeMax()};
					bmRegistration.flagRange=true;
					bmRegistration.minBlockVariance=10;
					bmRegistration.displayRegistration=2;
					bmRegistration.displayR2=false;
					bmRegistration.returnComposedTransformationIncludingTheInitialTransformationGiven=true;
					bmRegistration.levelMin=-1;
					trs[iMov]=bmRegistration.runBlockMatching(null);
					bmRegistration.closeLastImages();
				}
				else trs[iMov]=new ItkTransform();
				trs[iMov].addTransform(tr);
				IJ.log("Registration of T1 times finished !");
				imgT1Seq[iMov]=trs[iMov].transformImage(imgT2Seq[0],imgT1Seq[iMov],false);
			}
		}
		if(comp.isVisible())comp.close();

		
		
		//Adding the first image of T2 seq into T1 seq, as the last repetition time
		if(crossFit) {
			int nOld=nT1;
			ImagePlus[]imSeqNew=new ImagePlus[nT1+1];//TODO : verifier que les T1 succ ont les memes dims, pareil pour les T2
			for(int i=0;i<nOld;i++)imSeqNew[i]=VitimageUtils.imageCopy(imgT1Seq[i]);
			imSeqNew[nOld]=VitimageUtils.imageCopy(imgT2Seq[0]);
			for(int z=0;z<imSeqNew[nOld].getStackSize();z++)imSeqNew[nOld].getStack().setSliceLabel(imSeqNew[nOld].getStack().getSliceLabel(z+1).replace("T2SEQ", "T1SEQ"), z+1);
			imgT1Seq=imSeqNew;			
			int[]cpTab=new int[nOld+1];
			for(int i=0;i<nOld;i++)cpTab[i]=this.valT1Trs[i];
			cpTab[nOld]=TrValueUsedForT2Seq;
			this.valT1Trs=cpTab;
		}
				
		IJ.log("Registration finished !");
	}
	
	

	

	
	
	
	public void stackSaveAndShowResults() {
		ImagePlus []tabIm=new ImagePlus[ (haveT1Seq ? (1+this.valT1Trs.length) : 0) +  (haveT2Seq ? (2+ this.valT2Tes.length) : 0) ];  
		int nZ=0;
		if(haveT2Seq)nZ=imgT2Seq[0].getNSlices();
		else if(haveT1Seq)nZ=imgT1Seq[0].getNSlices();
		int curInd=0;
		int nMaps=3;
		int nT1=this.valT1Trs.length;
		int nT2=this.valT2Tes.length;
		if(haveT2Seq) {//Write M0 map
			tabIm[curInd]=VitimageUtils.imageCopy(M0map);
			tabIm[curInd].setDisplayRange(0,MRUtils.maxDisplayedM0 );
			IJ.run(tabIm[curInd],"32-bit","");
			for(int z=1;z<=nZ;z++) {
				tabIm[curInd].getStack().setSliceLabel("M0MAP", VitimageUtils.getCorrespondingSliceInHyperImage(tabIm[curInd], 0, z-1, 0));
			}
			curInd++;
		}
		if(haveT1Seq) {//Write T1 map
			tabIm[curInd]=VitimageUtils.imageCopy(T1map);
			tabIm[curInd].setDisplayRange(0,MRUtils.maxDisplayedT1 );
			IJ.run(tabIm[curInd],"32-bit","");
			for(int z=1;z<=nZ;z++) {
				tabIm[curInd].getStack().setSliceLabel("T1MAP", VitimageUtils.getCorrespondingSliceInHyperImage(tabIm[curInd], 0, z-1, 0));
			}
			curInd++;
		}
		if(haveT2Seq) {//Write T2 map
			tabIm[curInd]=VitimageUtils.imageCopy(T2map);
			tabIm[curInd].setDisplayRange(0,MRUtils.maxDisplayedT2 );
			IJ.run(tabIm[curInd],"32-bit","");
			for(int z=1;z<=nZ;z++) {
				tabIm[curInd].getStack().setSliceLabel("T2MAP", VitimageUtils.getCorrespondingSliceInHyperImage(tabIm[curInd], 0, z-1, 0));
			}
			curInd++;
		}

		if(haveT1Seq) {//Write T1seq
			for(int i=0;i<this.valT1Trs.length;i++) {
				tabIm[curInd]=VitimageUtils.imageCopy(imgT1Seq[i]);
				tabIm[curInd].setDisplayRange(0,MRUtils.maxDisplayedM0 );
				IJ.run(tabIm[curInd++],"32-bit","");
			}
		}
		if(haveT2Seq) {//Write T2 seq
			for(int i=0;i<this.valT2Tes.length;i++) {
				tabIm[curInd]=VitimageUtils.imageCopy(imgT2Seq[i]);
				tabIm[curInd].setDisplayRange(0,MRUtils.maxDisplayedM0 );
				IJ.run(tabIm[curInd++],"32-bit","");
			}
		}

		int nC=tabIm.length;
		ImagePlus hyperImg=Concatenator.run(tabIm);
		hyperImg=HyperStackConverter.toHyperStack(hyperImg, nC,nZ,1,"xyztc","Grayscale");

		curInd=1;
		if(haveT2Seq) {
			hyperImg.setC(curInd++);
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );
			IJ.run(hyperImg,"Fire","");
		}

		if(haveT1Seq) {
			hyperImg.setC(curInd++);
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT1 );
			IJ.run(hyperImg,"Fire","");
		}

		if(haveT2Seq) {
			hyperImg.setC(curInd++);
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );
			IJ.run(hyperImg,"Fire","");
		}

		for(int c=curInd;c<=nC;c++) {
			hyperImg.setC(c);
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedM0 );				
			IJ.run(hyperImg,"Fire","");
		}
		hyperImg.setTitle("Hyperimage");
		String name=new File(outputDir,"hyperImage_"+(new SimpleDateFormat("yyyy-MM-dd_hh-mm").format(new Date())+".tif")).getAbsolutePath();
		IJ.saveAsTiff(hyperImg,name);
		hyperImg.close();
		
		ImagePlus result=IJ.openImage(name);
		result.show();
		IJ.showMessage("Hyper MRI image saved as : "+name);
	}










}
