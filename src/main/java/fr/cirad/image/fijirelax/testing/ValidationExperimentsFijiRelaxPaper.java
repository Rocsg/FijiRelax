package fr.cirad.image.fijirelax.testing;

/*
t2Init_dif = xData(1) - xData(end-1);
t2Init = t2Init_dif/log(yDat(end-1)/yDat(1));

if t2Init<=0 || isnan(t2Init),
    t2Init=30;
end
*/

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fr.cirad.image.common.Timer;
import fr.cirad.image.common.TransformUtils;
import fr.cirad.image.common.VitimageUtils;
import fr.cirad.image.fijiyama.RegistrationAction;
import fr.cirad.image.registration.BlockMatchingRegistration;
import fr.cirad.image.registration.ItkTransform;

import fr.cirad.image.fijirelax.mrialgo.HyperMap;
import fr.cirad.image.fijirelax.mrialgo.MRDataType;
import fr.cirad.image.fijirelax.mrialgo.MRUtils;
import fr.cirad.image.fijirelax.mrialgo.NoiseManagement;
import fr.cirad.image.fijirelax.mrialgo.RiceEstimator;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;

public class ValidationExperimentsFijiRelaxPaper {

	public ValidationExperimentsFijiRelaxPaper() {}
		
	public static void testRec() {
		ImagePlus imgRef=IJ.openImage("/home/fernandr/Bureau/testRef.tif");
		ImagePlus imgMov=IJ.openImage("/home/fernandr/Bureau/testMov.tif");
		imgRef.setDisplayRange(0,1000);
		imgMov.setDisplayRange(0,1000);
		imgRef.show();
		imgMov.show();
		BlockMatchingRegistration bmRegistration;
		RegistrationAction regAct=new RegistrationAction();
		regAct.defineSettingsSimplyFromTwoImages(imgRef, imgMov);
		regAct.typeAutoDisplay=2;
		regAct.higherAcc=0;
		regAct.levelMaxLinear=2;
		regAct.levelMinLinear=1;
		regAct.strideX*=2;
		regAct.strideY*=2;
		regAct.strideZ*=2;
		System.out.println(regAct.readableString());
		bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(imgRef,imgMov ,regAct);
		bmRegistration.displayRegistration=2;
//		bmRegistration.displayRegistration=0;
		ItkTransform tr=bmRegistration.runBlockMatching(null,false);
		bmRegistration.closeLastImages();
		ImagePlus imgRes=tr.transformImage(imgRef,imgMov);
		VitimageUtils.compositeOf(imgRef, imgMov,"Before").show();;
		VitimageUtils.compositeOf(imgRef, imgRes,"After").show();;

		
	}
	
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		//processExperimentData(true,false,false,false,false);	
		exp_02_generateMaps();
		//exp_03_compute_maps_for_real_data2();
		//generateExperimentData(false,false,true,false,false);
		//exp_03_simulate_data_v2();
	}
	
	///////////////// Definitions for experiments
	static String dirExp01="/home/fernandr/Bureau/FijiRelax_PrepaDOI/Paper_experiments/Exp_01_Registration_T1";
	static String dirExp02="/home/fernandr/Bureau/FijiRelax_PrepaDOI/Paper_experiments/Exp_02_Rice_fit";
	static String dirExp03="/home/fernandr/Bureau/FijiRelax_PrepaDOI/Paper_experiments/Exp_03_Joint_fit";
	static String dirExp04="/home/fernandr/Bureau/FijiRelax_PrepaDOI/Paper_experiments/Exp_04_Computation_time";
	
	//Constants used for simulations of exp01
	public static String imgBaseName_exp_01="SE_TR_";
	public static int[] tr_exp_01=new int[] {600,1200,1800,10000};
	public static int [][]spacingSimulations_exp_01=new int[][] {{80,80},{40,25,55,20,40,60},{35,80,80,106,106,106},{25,10,10,5,5,5},{14,40,66},{14,40,66},{15,20,30},{7}};//rectangle size, centersX, centers Y, radius
	public static float t1Up_exp_01=2000;
	public static float t1Down_exp_01=1000;
	public static float PDUp_exp_01=2000;
	public static float PDDown_exp_01=1000;
	public static double []deltax_exp_01=new double[] {-5.4,-3.8,-2,0};
	public static double []deltay_exp_01=new double[] {-2.1,-1.5,-0.6,0};
	public static int dimZ_exp_01=1;
	public static int [][]indexCategories_exp_01=new int[][] {{0},{1,2},{3,4,5}};
	public static int nCategories_exp_01=3;
	public static double sigmaGauss_exp_01=1;
	
	//Constants used for simulations of exp02
	public static double PD_exp02=1000;
	public static int []spacingSimulations_exp_02=new int[] {8,10,10,10};
	public static double[] valsT2Simulation_exp_02=new double[] {20,40,60,80};
	public static int nRepet_exp_02=4;
	public static int dimZ_exp_02=12;
	public static int nEchoes_exp_02=16;
	public static double teSpacing_exp_02=10; //12.8 chez QMRLAB
	public static int tr_exp_02=10000;
	public static int valT1_exp_02=500;
	
	//Constants used for simulations of exp03
	public static int dimZ_exp_03=1;
	public static int []spacingSimulations_exp_03=new int[] {100,100};
	public static double PD_exp03=1000;
	public static double T1_exp03=1000;
	public static double T2_exp03=60;
	private static int nEchoes_exp_03=16;
	private static int nRepet_exp_03=4;


     
 	/** Generate data*/////////////////////////////////////////////////////////////////////////////////
	public static void generateExperimentData(boolean genExp1,boolean genExp2,boolean genExp3,boolean genExp4,boolean genExp5) {
		if(genExp1) {
			//Exp 01 : Reg T1
				//Build spin echoes simulated data
				exp_01_generateMaps();
				exp_01_generateSpinEchoT1Mono();
				exp_01_registerEchoes();
				exp_01_generateSpinEchoT1MonoRealData();
				exp_01_registerEchoesRealData();
				
				//Take B001 J0 and set the echo spin in directory somewhere. Select only T1 relax + first T2 echoes
				
		}

		if(genExp2) {
			//Exp 02 : Rice mono
				//Build spin echoes simulated data
				//exp_02_generateMaps();
				//exp_02_generateSpinEchoT2Mono();
				//exp_02_estimateMapsFromSimulatedEchoes(MRUtils.T2_MONO_RICE);
				//exp_02_estimateMapsFromSimulatedEchoes(MRUtils.T2_MONO_BIAS);
				//exp_02_estimateMapsFromSimulatedEchoes(MRUtils.T2_MONO);
				
				//Take B001 J105 and set the echo spin in directory somewhere. Select only T2 relax echoes
			exp_02_estimateMapsFromEchoes(true,MRUtils.T2_MONO_RICE);
			exp_02_estimateMapsFromEchoes(true,MRUtils.T2_MONO_BIAS);
			exp_02_estimateMapsFromEchoes(true,MRUtils.T2_MONO);

		}			
		
		if(genExp3) {
			//Exp 03 : Rice biexp
				//Build spin echoes simulated data
				//Take B001 J105 and set the echo spin in directory somewhere. Select only T2 relax echoes
			//	exp_03_generateMaps();
		//		exp_03_generateSpinEchoData();
	//			exp_03_estimateMapsFromEchoes();
//				exp_03_simulate_data();
				exp_03_compute_maps_for_real_data();
		}

		if(genExp4) {
			//Exp 04 : Crossfit
				//Build spin echoes simulated data
				//Take B001 J105 and set the echo spin in directory somewhere. Select only T2 relax echoes
		}

		if(genExp5) {
			//Exp 05 : Computation time
			//Take B001 J105, select 24 slices in the center with object in it, select a crop in it (establish crop size using matlab Qmrlab)
			//HyperMap hyperMap=HyperMap.importHyperMapFromT2RelaxNifti("/home/fernandr/Bureau/FijiRelax/Exp_02_Rice_Mono/Source_data/Data_QMRLAB/SEdata.tif","BRAIN",1000, 12.8);
			//hyperMap.computeMapsAgainAndMask(MRUtils.SIMPLEX,true,NoiseManagement.RICE,false,null,5);
			//hyperMap.hyperImg.show();
			/*HyperMap hyperMap=new HyperMap(IJ.openImage("/home/fernandr/Bureau/FijiRelax/Exp_02_Rice_Mono/Source_data/mono_t2_demo/mono_t2_data/SSM1_0930-2.tif"));
			hyperMap.computeMapsAgainAndMask(true,MRUtils.SIMPLEX,true,false);
			hyperMap.hyperImg.show();*/
			//generateExperimentData(false, false, true, false, false);

		
		}

		
	}
	
	
	
	
	/** Process data*/////////////////////////////////////////////////////////////////////////////////
	public static void processExperimentData(boolean procExp1,boolean procExp2,boolean procExp3,boolean procExp4,boolean procExp5) {
		if(procExp1) {
			//Exp 01 : Reg T1
/*			exp_01_estimateMapsFromSimulatedEchoes();*/
			exp_01_processResults();
			//exp_01_estimateMapsFromSimulatedEchoesRealData();

			//
			
		}

		if(procExp2) {
			//Exp 02 : Rice mono
				//compute diffs and tab of results. Also display python code for simulating curve
				//exp_02_processResults();

		}			
		
		if(procExp3) {
			//Exp 03 : Rice biexp
				//
		
		}

		if(procExp4) {
			//Exp 04 : Crossfit
				//
		}

		if(procExp5) {
			//

		
		}

		
	}
	
	
	public static void generatePythonCurve() {
		double sigma=45;
		double teSpacing=10;
		double Tr=10000;
		double T1=1000;
		double T2=60;
		double PD=1000;
		int factor=5;
		int nEchoes=30;
		int nRepetitions=4;
		int N=50;
		double acc=0;
		int nEcAc=(1+nEchoes)*factor+1;
		double []tes=new double[nEchoes];
		double []tesAcc=new double[nEcAc];
		double []trs=new double[nEchoes];
		double []trsAcc=new double[nEcAc];
		double [][]data=new double[N][nEchoes];
		double [][]dataFittenExp=new double[N][nEcAc];
		double [][]dataFittenBias=new double[N][nEcAc];
		double [][]dataBiasBias=new double[N][nEcAc];
		double []T2Exp=new double[N];
		double []T2Bias=new double[N];
		double []BiasBias=new double[N];
		for(int ite=0;ite<nEchoes;ite++) {
			tes[ite]=teSpacing*(ite+1);
			trs[ite]=Tr;
		}
		for(int ite=0;ite<nEcAc;ite++) {
			tesAcc[ite]=teSpacing*(ite)/factor;
			trsAcc[ite]=Tr;
		}
		for(int n=0;n<N;n++) {
			for(int ite=0;ite<nEchoes;ite++) {
				acc=0;
				for(int r=0;r<nRepetitions;r++) {
					acc+=RiceEstimator.getRandomRiceRealization(
						MRUtils.getFitFunctionValue(trs[ite], tes[ite],new double[] {PD,T1,T2}, 0, MRUtils.T1T2_MONO_RICE),sigma);
				}
				data[n][ite]=(double)(acc/nRepetitions); 
			}			
			double[]params=(double[]) MRUtils.makeFit(trs, tes, data[n], MRUtils.T2_MONO, MRUtils.SIMPLEX,400, sigma, false)[0];
			T2Exp[n]=params[1];
			for(int ite=0;ite<nEcAc;ite++) {
				dataFittenExp[n][ite]=MRUtils.getFitFunctionValue(trsAcc[ite], tesAcc[ite], params, sigma, MRUtils.T2_MONO);
			}
			params=(double[]) MRUtils.makeFit(trs, tes, data[n], MRUtils.T2_MONO_BIAS, MRUtils.SIMPLEX,400, sigma, false)[0];
			T2Bias[n]=params[1];
			BiasBias[n]=params[2];
			for(int ite=0;ite<nEcAc;ite++) {
				dataFittenBias[n][ite]=MRUtils.getFitFunctionValue(trsAcc[ite], tesAcc[ite], params, sigma, MRUtils.T2_MONO_BIAS);
				dataBiasBias[n][ite]=BiasBias[n];
			}
		}
		String s1="te=[";
		String s2="teAcc=[";
		String s3="data=[\n";
		String s4="fitExp=[\n";
		String s5="fitBias=[\n";
		String s52="fitBiasBias=[\n";
		String s6="T2Exp=[\n";
		String s7="T2Bias=[\n";
		String s8="BiasBias=[\n";
		for(int ite=0;ite<nEchoes;ite++) s1+=(tes[ite])+(ite==nEchoes-1 ? "" : ",");
		for(int ite=0;ite<nEcAc;ite++) s2+=(tesAcc[ite])+(ite==nEcAc-1 ? "" : ",");

		for(int n=0;n<N;n++) {
			s3+="[";
			s4+="[";
			s5+="[";
			s52+="[";
			for(int ite=0;ite<nEchoes;ite++)s3+=(data[n][ite])+(ite==nEchoes-1 ? "" : ",");
			for(int ite=0;ite<nEcAc;ite++) {s4+=(dataFittenExp[n][ite])+(ite==nEcAc-1 ? "" : ","); s5+=(dataFittenBias[n][ite])+(ite==nEcAc-1 ? "" : ",");s52+=(dataBiasBias[n][ite])+(ite==nEcAc-1 ? "" : ",");}
			s3+="]"+(n==N-1 ? "" : ",")+"\n";
			s4+="]"+(n==N-1 ? "" : ",")+"\n";
			s5+="]"+(n==N-1 ? "" : ",")+"\n";
			s52+="]"+(n==N-1 ? "" : ",")+"\n";
			s6+=(T2Exp[n])+(n==N-1 ? "" : ",");
			s7+=(T2Bias[n])+(n==N-1 ? "" : ",");
			s8+=(BiasBias[n])+(n==N-1 ? "" : ",");
		}
		s1+="]\n";
		s2+="]\n";
		s3+="]\n";
		s4+="]\n";
		s5+="]\n";
		s52+="]\n";
		s6+="]\n";
		s7+="]\n";
		s8+="]\n";
		String s9=
				"import matplotlib.pyplot as plt\n"+
				"index=1\n"+
				"plt.subplot(121)\n"+
				"plt.plot(te,data[index],'bx')\n"+
				"plt.plot(teAcc,fitExp[index],'r')\n"+
				"plt.subplot(122)\n"+
				"plt.plot(te,data[index],'bx')\n"+
				"plt.plot(teAcc,fitBias[index],'r')\n"+
				"print(T2Exp)\n"+
				"print(T2Bias)\n\n";


		VitimageUtils.writeStringInFile((s1+s2+s3+s4+s5+s52+s6+s7+s8+s9),"/home/fernandr/Bureau/log.txt");
		
	}
	

	
	
	public static void exp_01_generateSpinEchoT1Mono() {
		//Collect simulated maps
		ImagePlus hyper=IJ.openImage(dirExp01+"/Source_data/Simulated_map.tif");
		ImagePlus[]maps=VitimageUtils.stacksFromHyperstackFastBis(hyper);

		int nTr=tr_exp_01.length;
		ImagePlus[]tab=new ImagePlus[nTr];		
		for(int tr=0;tr<nTr;tr++)tab[tr]=VitimageUtils.nullImage(maps[0]);
		float[][]echoVals=new float[nTr][];
		float[][]mapVals=new float[2][];
		int dimZ=tab[0].getNSlices();
		int dimX=tab[0].getWidth();
		int dimY=tab[0].getHeight();
		for(int z=0;z<dimZ;z++) {
			for(int i=0;i<maps.length;i++) {
				mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
			}
			for(int tr=0;tr<nTr;tr++) {
				echoVals[tr]=(float[]) tab[tr].getStack().getPixels(z+1);
				for(int x=0;x<dimX;x++) {
					for(int y=0;y<dimY;y++) {
						int x0=(int) Math.floor(x+deltax_exp_01[tr]);
						int y0=(int) Math.floor(y+deltay_exp_01[tr]);
						if(x0<0)x0=0;if(x0>dimX-2)x0=dimX-2;
						if(y0<0)y0=0;if(y0>dimY-2)y0=dimY-2;
						int x1=x0+1;
						int y1=y0+1;
						double dx=x+deltax_exp_01[tr]-x0;
						double dy=y+deltay_exp_01[tr]-y0;
						double accu=0;
						int offset;
						offset=(y0)*dimX+(x0);
						accu+=(1-dx) * (1-dy) * MRUtils.getFitFunctionValue(tr_exp_01[tr], 0,new double[] {mapVals[0][offset],mapVals[1][offset]}, 0, MRUtils.T1_MONO_RICE);
						offset=((y0)+1)*dimX+(x0);
						accu+=(1-dx) * (dy) * MRUtils.getFitFunctionValue(tr_exp_01[tr], 0,new double[] {mapVals[0][offset],mapVals[1][offset]}, 0, MRUtils.T1_MONO_RICE);
						offset=(y0)*dimX+((x0)+1);
						accu+=(dx) * (1-dy) * MRUtils.getFitFunctionValue(tr_exp_01[tr], 0,new double[] {mapVals[0][offset],mapVals[1][offset]}, 0, MRUtils.T1_MONO_RICE);
						offset=((y0)+1)*dimX+((x0)+1);
						accu+=(dx) * (dy) * MRUtils.getFitFunctionValue(tr_exp_01[tr], 0,new double[] {mapVals[0][offset],mapVals[1][offset]}, 0, MRUtils.T1_MONO_RICE);
						offset=(y)*dimX+(x);
						echoVals[tr][offset]=	(float)accu;
					}
				}
			}
		}
	
		for(int tr=0;tr<nTr;tr++) {
			tab[tr].setTitle(""+tr);
			tab[tr].setDisplayRange(0, PDUp_exp_01*5);
			IJ.run(tab[tr],"Fire","");
		}
		ImagePlus hyper2=VitimageUtils.hyperStackingChannels(tab);
		hyper2.setDisplayRange(0, PDUp_exp_01*5);
		IJ.run(hyper2,"Fire","");
		IJ.save(hyper2,dirExp01+"/Source_data/SimulatedData_echoes.tif");
		
	}
		
	public static void exp_01_registerEchoes() {
		ImagePlus hyper=IJ.openImage(dirExp01+"/Source_data/SimulatedData_echoes.tif");
		ImagePlus []echoes=VitimageUtils.stacksFromHyperstackFastBis(hyper);
		int nTr=echoes.length;

		for(int i=0;i<nTr-1;i++) {
			BlockMatchingRegistration bmRegistration;
			RegistrationAction regAct=new RegistrationAction();
			regAct.defineSettingsSimplyFromTwoImages(echoes[nTr-1], echoes[i]);
			regAct.typeAutoDisplay=2;
			regAct.higherAcc=1;
			System.out.println(regAct.readableString());
			bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(echoes[nTr-1], echoes[i],regAct);
			bmRegistration.minBlockVariance=0.1;
			ItkTransform tr=bmRegistration.runBlockMatching(null,false);
			bmRegistration.closeLastImages();
		    echoes[i]=tr.transformImage(echoes[nTr-1], echoes[i]);
		}

		for(int tr=0;tr<nTr;tr++) {
			echoes[tr].setTitle(""+tr);
			echoes[tr].setDisplayRange(0, PDUp_exp_01);
			IJ.run(echoes[tr],"Fire","");
		}
		hyper=VitimageUtils.hyperStackingChannels(echoes);
		hyper.setDisplayRange(0, PDUp_exp_01);
		IJ.run(hyper,"Fire","");
		IJ.save(hyper,dirExp01+"/Source_data/SimulatedData_echoes_registered.tif");
	}
	
	public static void exp_01_estimateMapsFromSimulatedEchoes() {
		int nTr=tr_exp_01.length;
		ImagePlus hyper=IJ.openImage(dirExp01+"/Source_data/SimulatedData_echoes.tif");
		ImagePlus []echoes=VitimageUtils.stacksFromHyperstackFastBis(hyper);
		ImagePlus hyperReg=IJ.openImage(dirExp01+"/Source_data/SimulatedData_echoes_registered.tif");
		ImagePlus []echoesReg=VitimageUtils.stacksFromHyperstackFastBis(hyperReg);
		ImagePlus []maps=new ImagePlus[2];
		ImagePlus []mapsReg=new ImagePlus[2];
		maps[0]=VitimageUtils.nullImage(echoes[0]);
		maps[1]=VitimageUtils.nullImage(echoes[0]);
		mapsReg[0]=VitimageUtils.nullImage(echoes[0]);
		mapsReg[1]=VitimageUtils.nullImage(echoes[0]);


		double[]tabData=new double[nTr];
		double[]tabDataReg=new double[nTr];
		int dimZ=dimZ_exp_01;
		int[][] spacing=spacingSimulations_exp_01;
		int dimX=spacing[0][0];int dimY=spacing[0][1]; //image size
		double []tabTe=new double[nTr];for(int i=0;i<nTr;i++)tabTe[i]=0;
		double []tabTr=new double[nTr];for(int i=0;i<nTr;i++)tabTr[i]=tr_exp_01[i];
		
		float[][]mapVals=new float[2][];
		float[][]mapValsReg=new float[2][];
		float[][]echoVals=new float[nTr][];
		float[][]echoValsReg=new float[nTr][];
		for(int z=0;z<dimZ;z++) {
			for(int i=0;i<2;i++)mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
			for(int i=0;i<nTr;i++)echoVals[i]=(float[]) echoes[i].getStack().getPixels(z+1);
			for(int i=0;i<2;i++)mapValsReg[i]=(float[]) mapsReg[i].getStack().getPixels(z+1);
			for(int i=0;i<nTr;i++)echoValsReg[i]=(float[]) echoesReg[i].getStack().getPixels(z+1);

			for(int x=0;x<dimX;x++) {
				for(int y=0;y<dimY;y++) {
					int offset=(y*dimX+x);
					for(int tr=0;tr<nTr;tr++) {tabData[tr]=echoVals[tr][offset];}
					double[]params=(double[]) MRUtils.makeFit(tabTr, tabTe, tabData, MRUtils.T1_MONO_RICE, MRUtils.SIMPLEX, 400, 0, false)[0];
					mapVals[0][y*dimX+x]=(float) params[0];
					mapVals[1][y*dimX+x]=(float) params[1];
					for(int tr=0;tr<nTr;tr++)tabDataReg[tr]=echoValsReg[tr][offset];
					params=(double[]) MRUtils.makeFit(tabTr, tabTe, tabDataReg, MRUtils.T1_MONO_RICE, MRUtils.SIMPLEX, 400, 0, false)[0];
					mapValsReg[0][y*dimX+x]=(float) params[0];
					mapValsReg[1][y*dimX+x]=(float) params[1];
				}
				
			}
		}
		ImagePlus map=VitimageUtils.hyperStackingChannels(maps);
		map.setC(1);map.setDisplayRange(0, PDUp_exp_01*3);		IJ.run(map,"Fire","");
		map.setC(2);map.setDisplayRange(0, PDUp_exp_01*3);		IJ.run(map,"Fire","");
		IJ.saveAsTiff(map,dirExp01+"/Results/"+("SimulatedData_maps_estimation.tif"));

		ImagePlus mapReg=VitimageUtils.hyperStackingChannels(mapsReg);
		mapReg.setC(1);mapReg.setDisplayRange(0, PDUp_exp_01*3);		IJ.run(mapReg,"Fire","");
		mapReg.setC(2);mapReg.setDisplayRange(0, PDUp_exp_01*3);		IJ.run(mapReg,"Fire","");
		IJ.saveAsTiff(mapReg,dirExp01+"/Results/"+("SimulatedData_maps_estimation_after_registration.tif"));
	}

	
	
	public static ImagePlus gatherResultsWithVariousSNROnSlices(ImagePlus img) {
		img.show();
		int border1=21;
		int border2=41;
		int C=img.getNChannels();
		int Z=img.getNSlices();
		int T=img.getNFrames();
		int X=img.getWidth();
		int Y=img.getHeight();
		ImagePlus out=new Duplicator().run(img,1,C,1,Z,1,1);
		for(int c=0;c<C;c++)for(int z=0;z<Z;z++) {
			float[]tabOut=(float[])out.getStack().getPixels(VitimageUtils.getCorrespondingSliceInHyperImage(img, c, Z, 0));

			float[]tabIn0=(float[])img.getStack().getPixels(VitimageUtils.getCorrespondingSliceInHyperImage(img, c, Z, 0));
			float[]tabIn1=(float[])img.getStack().getPixels(VitimageUtils.getCorrespondingSliceInHyperImage(img, c, Z, 1));
			float[]tabIn2=(float[])img.getStack().getPixels(VitimageUtils.getCorrespondingSliceInHyperImage(img, c, Z, 2));
			for(int y=0;y<Y;y++) {
				for(int x=0;x<border1;x++) tabOut[y*X+x]=tabIn0[y*X+x];
				for(int x=border1;x<border2;x++) tabOut[y*X+x]=tabIn1[y*X+x];
				for(int x=border2;x<X;x++) tabOut[y*X+x]=tabIn2[y*X+x];
			}
		}
		return out;
	}
	
	
	
	/** Functions for exp01 in chrononological order of calling*/////////////////////////////////////////////////////////////////////////////////	
	public static void exp_01_generateMaps() {
		ImagePlus[]maps=new ImagePlus[2];
		int[][] spacing=spacingSimulations_exp_01;
		int dimX=spacing[0][0];int dimY=spacing[0][1]; //image size
		int nCircles=spacing[1].length; // characteristics of simulated objects (circles)
		int[]centersX=spacing[4];int[]centersY=spacing[5];int[] ratio=spacing[6];int radius=spacing[7][0]; 
		ImagePlus img=IJ.createImage("map", dimX, dimY, dimZ_exp_01, 32);
		for(int i=0;i<maps.length;i++)maps[i]=VitimageUtils.nullImage(img);

		
		float[][]mapVals=new float[2][];
		for(int z=0;z<dimZ_exp_01;z++) {
			for(int i=0;i<maps.length;i++)mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);

			//Set all image to target non-null background value
			for(int x=0;x<dimX;x++){
				for(int y=0;y<dimY;y++){
					int offset=y*dimX+x;
					mapVals[0][offset]=PDDown_exp_01; mapVals[1][offset]=t1Down_exp_01;
				}
			}
			
			//Draw objects
			for(int circX=0;circX<centersX.length;circX++) {
				for(int circY=0;circY<centersY.length;circY++) {
					for(int dx=-radius;dx<=radius;dx++){
						for(int dy=-radius;dy<=radius;dy++){
							double ray=Math.sqrt(dx*dx+dy*dy);
							if(ray<radius) {
								int offset=(centersY[circY]+dy)*dimX+(centersX[circX]+dx);
								mapVals[0][offset]*=(ratio[circX]/10.0);
								mapVals[1][offset]*=(ratio[circY]/10.0);
							}							
						}
					}
				}
			}
		}
		
		for(int i=0;i<maps.length;i++)maps[i]=VitimageUtils.gaussianFilteringIJ(maps[i], sigmaGauss_exp_01, sigmaGauss_exp_01, 0);
		ImagePlus hyper=VitimageUtils.hyperStackingChannels(maps);
		hyper.setC(1);hyper.setDisplayRange(0, PDUp_exp_01*3);		IJ.run(hyper,"Fire","");
		hyper.setC(2);hyper.setDisplayRange(0, PDUp_exp_01*3);		IJ.run(hyper,"Fire","");
		IJ.save(hyper,dirExp01+"/Source_data/SimulatedData_expected_maps.tif");
	}

	
	
	
	public static void exp_01_processResults() {
		ImagePlus []maps=VitimageUtils.stacksFromHyperstackFastBis(IJ.openImage(dirExp01+"/Results/"+("SimulatedData_maps_estimation.tif")));
		ImagePlus []mapsReg=VitimageUtils.stacksFromHyperstackFastBis(IJ.openImage(dirExp01+"/Results/"+("SimulatedData_maps_estimation_after_registration.tif")));
		ImagePlus []mapsInit=VitimageUtils.stacksFromHyperstackFastBis(IJ.openImage(dirExp01+"/Source_data/"+("SimulatedData_expected_maps.tif")));

		int dimZ=dimZ_exp_01;
		int[][] spacing=spacingSimulations_exp_01;
		int dimX=spacing[0][0];int dimY=spacing[0][1]; //image size
		int[]centersX=spacing[4];int[]centersY=spacing[5];int[] ratio=spacing[6];int radius=spacing[7][0]; 
		
		float[][]mapVals=new float[2][];
		float[][]mapValsReg=new float[2][];
		float[][]mapValsInit=new float[2][];

		int nSamples=0;
		for(int dx=-radius;dx<=radius;dx++)for(int dy=-radius;dy<=radius;dy++){
			double ray=Math.sqrt(dx*dx+dy*dy);
			if(ray<radius) nSamples++;
		}
		nSamples*=dimZ_exp_01;
		int[][]incr=new int[centersX.length][centersX.length];
		
		double[][][][]tabValsAggregated=new double[3][centersX.length][centersX.length][nSamples];
		for(int z=0;z<dimZ_exp_01;z++) {
			for(int i=0;i<maps.length;i++) {
				mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
				mapValsReg[i]=(float[]) mapsReg[i].getStack().getPixels(z+1);
				mapValsInit[i]=(float[]) mapsInit[i].getStack().getPixels(z+1);
			}
			
			//Draw objects
			for(int circX=0;circX<centersX.length;circX++) {
				for(int circY=0;circY<centersY.length;circY++) {
					for(int dx=-radius;dx<=radius;dx++){
						for(int dy=-radius;dy<=radius;dy++){
							double ray=Math.sqrt(dx*dx+dy*dy);
							if(ray<radius) {
								int offset=(centersY[circY]+dy)*dimX+(centersX[circX]+dx);
								//System.out.print(circX+" , "+circY+" - Processing point x="+(centersX[circX]+dx)+" y="+(centersY[circY]+dy));
								tabValsAggregated[0][circX][circY][incr[circX][circY]]=mapVals[1][offset];
								tabValsAggregated[1][circX][circY][incr[circX][circY]]=mapValsReg[1][offset];
								tabValsAggregated[2][circX][circY][incr[circX][circY]]=mapValsInit[1][offset];
								//System.out.println("  Collected : "+tabValsAggregated[2][circX][circY][incr[circX][circY]]);
								incr[circX][circY]++;
							}							
						}
					}
					System.out.println("Mean="+VitimageUtils.statistics1D(tabValsAggregated[2][circX][circY])[0]);
				}
			}
		}
		
		String[][]tabResults=new String[7][10];
		tabResults[0]=new String[] {"Without reg. PDx1.5","Without reg. PDx2","Without reg. PDx3", "With reg. PDx1.5","With reg. PDx2","With reg. PDx3", "Expected values PDx1.5","Expected values PDx2","Expected values PDx3"};
		tabResults[1][0]="T1x1.5 - mean+-std";		
		tabResults[2][0]="T1x1.5 - rel %";		
		tabResults[3][0]="T1x2 - mean+-std";		
		tabResults[4][0]="T1x2 - rel %";		
		tabResults[5][0]="T1x3 - mean+-std";		
		tabResults[6][0]="T1x3 - rel %";		
		System.out.println("Nsamples="+nSamples);
		String[]cases=new String[]{"Non reg", "Reg", "Init values"};
		for(int c=0;c<cases.length;c++) {
			System.out.println("\n\nResult "+cases[c]);
			for(int i=0;i<centersX.length;i++) {
				System.out.println("");
				System.out.println("------ T2 ratio : "+ratio[i]+"  ------");
				for(int j=0;j<centersY.length;j++) {
					double[]stats=VitimageUtils.statistics1D(tabValsAggregated[c][i][j]);
					double[]statsInit=VitimageUtils.statistics1D(tabValsAggregated[2][i][j]);
					System.out.print(" ["+doudou(stats[0])+" +-"+doudou(stats[1])+"]" );			
					tabResults[1+i*2+0][1+c*3+j]=""+doudou(stats[0])+" +-"+doudou(stats[1]);
				}
				System.out.println();
				for(int j=0;j<centersY.length;j++) {
					double[]stats=VitimageUtils.statistics1D(tabValsAggregated[c][i][j]);
					double[]statsInit=VitimageUtils.statistics1D(tabValsAggregated[2][i][j]);
					System.out.print(" ["+dou(100*(stats[0]-statsInit[0])/statsInit[0])+"]");			
					tabResults[1+i*2+1][1+c*3+j]=""+dou(100*(stats[0]-statsInit[0])/statsInit[0]);
				}
				System.out.println();
			}
		}
		writeStringTabInExcelFile(tabResults,dirExp01+"/Results/"+("table.csv"));

	}

	public static void processTableS1(ImagePlus expectedMaps,ImagePlus mapsWithout, ImagePlus mapsWith, String pathToCSV) {
		expectedMaps.show();
		expectedMaps.setTitle("Expected");
		ImagePlus []maps=VitimageUtils.stacksFromHyperstackFastBis(mapsWithout);
		ImagePlus []mapsReg=VitimageUtils.stacksFromHyperstackFastBis(mapsWith);
		ImagePlus []mapsInit=VitimageUtils.stacksFromHyperstackFastBis(expectedMaps);

		int dimZ=dimZ_exp_01;
		int[][] spacing=spacingSimulations_exp_01;
		int dimX=spacing[0][0];int dimY=spacing[0][1]; //image size
		int[]centersX=spacing[4];int[]centersY=spacing[5];int[] ratio=spacing[6];int radius=spacing[7][0]; 
		
		float[][]mapVals=new float[2][];
		float[][]mapValsReg=new float[2][];
		float[][]mapValsInit=new float[2][];

		int nSamples=0;
		for(int dx=-radius;dx<=radius;dx++)for(int dy=-radius;dy<=radius;dy++){
			double ray=Math.sqrt(dx*dx+dy*dy);
			if(ray<radius) nSamples++;
		}
		nSamples*=dimZ_exp_01;
		int[][]incr=new int[centersX.length][centersX.length];
		
		double[][][][]tabValsAggregated=new double[3][centersX.length][centersX.length][nSamples];
		for(int z=0;z<dimZ_exp_01;z++) {
			for(int i=0;i<2;i++) {
				mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
				mapValsReg[i]=(float[]) mapsReg[i].getStack().getPixels(z+1);
				mapValsInit[i]=(float[]) mapsInit[i].getStack().getPixels(z+1);
			}
			
			//Draw objects
			for(int circX=0;circX<centersX.length;circX++) {
				for(int circY=0;circY<centersY.length;circY++) {
					for(int dx=-radius;dx<=radius;dx++){
						for(int dy=-radius;dy<=radius;dy++){
							double ray=Math.sqrt(dx*dx+dy*dy);
							if(ray<radius) {
								int offset=(centersY[circY]+dy)*dimX+(centersX[circX]+dx);
								//System.out.print(circX+" , "+circY+" - Processing point x="+(centersX[circX]+dx)+" y="+(centersY[circY]+dy));
								tabValsAggregated[0][circX][circY][incr[circX][circY]]=mapVals[1][offset];
								tabValsAggregated[1][circX][circY][incr[circX][circY]]=mapValsReg[1][offset];
								tabValsAggregated[2][circX][circY][incr[circX][circY]]=mapValsInit[1][offset];
								//System.out.println("  Collected : "+tabValsAggregated[2][circX][circY][incr[circX][circY]]);
								incr[circX][circY]++;
							}							
						}
					}
					System.out.print("Mean NO="+VitimageUtils.dou(VitimageUtils.statistics1D(tabValsAggregated[0][circX][circY])[0]));
					System.out.print("Mean YES="+VitimageUtils.dou(VitimageUtils.statistics1D(tabValsAggregated[1][circX][circY])[0]));
					System.out.println("Mean EXPECT="+VitimageUtils.dou(VitimageUtils.statistics1D(tabValsAggregated[2][circX][circY])[0]));
				}
			}
		}
		
		String[][]tabResults=new String[7][10];
		tabResults[0]=new String[] {"Without reg. PDx1.5","Without reg. PDx2","Without reg. PDx3", "With reg. PDx1.5","With reg. PDx2","With reg. PDx3", "Expected values PDx1.5","Expected values PDx2","Expected values PDx3"};
		tabResults[1][0]="T1x1.5 - mean+-std";		
		tabResults[2][0]="T1x1.5 - rel %";		
		tabResults[3][0]="T1x2 - mean+-std";		
		tabResults[4][0]="T1x2 - rel %";		
		tabResults[5][0]="T1x3 - mean+-std";		
		tabResults[6][0]="T1x3 - rel %";		
		System.out.println("Nsamples="+nSamples);
		String[]cases=new String[]{"Non reg", "Reg", "Init values"};
		for(int c=0;c<cases.length;c++) {
			System.out.println("\n\nResult "+cases[c]);
			for(int i=0;i<centersX.length;i++) {
				System.out.println("");
				System.out.println("------ T2 ratio : "+ratio[i]+"  ------");
				for(int j=0;j<centersY.length;j++) {
					double[]stats=VitimageUtils.statistics1D(tabValsAggregated[c][i][j]);
					double[]statsInit=VitimageUtils.statistics1D(tabValsAggregated[2][i][j]);
					System.out.print(" ["+doudou(stats[0])+" +-"+doudou(stats[1])+"]" );			
					tabResults[1+i*2+0][1+c*3+j]=""+doudou(stats[0])+" +-"+doudou(stats[1]);
				}
				System.out.println();
				for(int j=0;j<centersY.length;j++) {
					double[]stats=VitimageUtils.statistics1D(tabValsAggregated[c][i][j]);
					double[]statsInit=VitimageUtils.statistics1D(tabValsAggregated[2][i][j]);
					System.out.print(" ["+dou(100*(stats[0]-statsInit[0])/statsInit[0])+"]");			
					tabResults[1+i*2+1][1+c*3+j]=""+dou(100*(stats[0]-statsInit[0])/statsInit[0]);
				}
				System.out.println();
			}
		}
		writeStringTabInExcelFile(tabResults,pathToCSV);
	}

	
	
	public static void exp_01_generateSpinEchoT1MonoRealData() {
		//Collect simulated maps

		int nTr=tr_exp_01.length;
		ImagePlus[]tab=new ImagePlus[nTr];		
		for(int tr=0;tr<nTr;tr++) {
			tab[tr]=IJ.openImage(dirExp01+"/Source_data/"+"SE_TR_"+tr_exp_01[tr]+".tif");
			IJ.run(tab[tr],"32-bit","");
		}
		for(int tr=0;tr<nTr;tr++) {
			tab[tr].setTitle(""+tr);
			tab[tr].setDisplayRange(0, PDUp_exp_01*5);
			IJ.run(tab[tr],"Fire","");
		}
		ImagePlus hyper2=VitimageUtils.hyperStackingChannels(tab);
		hyper2.setDisplayRange(0, PDUp_exp_01*5);
		IJ.run(hyper2,"Fire","");
		IJ.save(hyper2,dirExp01+"/Source_data/RealData_echoes.tif");
	}
	
	
	
		

		
	public static void exp_01_registerEchoesRealData() {
		ImagePlus hyper=IJ.openImage(dirExp01+"/Source_data/RealData_echoes.tif");
		ImagePlus []echoes=VitimageUtils.stacksFromHyperstackFastBis(hyper);
		int nTr=echoes.length;
		ItkTransform []trs=new ItkTransform[nTr-1];
		
		for(int i=0;i<nTr-1;i++) {
			BlockMatchingRegistration bmRegistration;
			RegistrationAction regAct=new RegistrationAction();
			regAct.defineSettingsSimplyFromTwoImages(echoes[nTr-1], echoes[i]);
			regAct.typeAutoDisplay=2;
			regAct.higherAcc=0;
			regAct.levelMaxLinear=2;
			regAct.levelMinLinear=-1;
			System.out.println(regAct.readableString());
			bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(echoes[nTr-1], echoes[i],regAct);
			ItkTransform tr=bmRegistration.runBlockMatching(null,false);
			bmRegistration.closeLastImages();
		    echoes[i]=tr.transformImage(echoes[nTr-1], echoes[i]);
		    trs[i]=tr;
		}

		for(int tr=0;tr<nTr;tr++) {
			echoes[tr].setTitle(""+tr);
			echoes[tr].setDisplayRange(0, PDUp_exp_01);
			IJ.run(echoes[tr],"Fire","");
			if(tr<nTr-1)System.out.println("\nTransform "+tr+"\n"+trs[tr].drawableString());
		}
		hyper=VitimageUtils.hyperStackingChannels(echoes);
		hyper.setDisplayRange(0, PDUp_exp_01);
		IJ.run(hyper,"Fire","");
		IJ.save(hyper,dirExp01+"/Source_data/RealData_echoes_registered.tif");
	}

	
	public static void exp_01_estimateMapsFromSimulatedEchoesRealData() {
		int nTr=tr_exp_01.length;
		ImagePlus hyper=IJ.openImage(dirExp01+"/Source_data/RealData_echoes.tif");
		ImagePlus []echoes=VitimageUtils.stacksFromHyperstackFastBis(hyper);
		ImagePlus hyperReg=IJ.openImage(dirExp01+"/Source_data/RealData_echoes_registered.tif");
		ImagePlus []echoesReg=VitimageUtils.stacksFromHyperstackFastBis(hyperReg);
		ImagePlus []maps=new ImagePlus[2];
		ImagePlus []mapsReg=new ImagePlus[2];
		maps[0]=VitimageUtils.nullImage(echoes[0]);
		maps[1]=VitimageUtils.nullImage(echoes[0]);
		mapsReg[0]=VitimageUtils.nullImage(echoes[0]);
		mapsReg[1]=VitimageUtils.nullImage(echoes[0]);


		double[]tabData=new double[nTr];
		double[]tabDataReg=new double[nTr];
		int dimZ=echoesReg[0].getNSlices();
		int dimX=echoesReg[0].getWidth();int dimY=echoesReg[0].getHeight();
		double []tabTe=new double[nTr];for(int i=0;i<nTr;i++)tabTe[i]=0;
		double []tabTr=new double[nTr];for(int i=0;i<nTr;i++)tabTr[i]=tr_exp_01[i];
		
		float[][]mapVals=new float[2][];
		float[][]mapValsReg=new float[2][];
		float[][]echoVals=new float[nTr][];
		float[][]echoValsReg=new float[nTr][];
	//	for(int z=0;z<dimZ;z++) {
		for(int z=0;z<1;z++) {
			for(int i=0;i<2;i++)mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
			for(int i=0;i<nTr;i++)echoVals[i]=(float[]) echoes[i].getStack().getPixels(z+1);
			for(int i=0;i<2;i++)mapValsReg[i]=(float[]) mapsReg[i].getStack().getPixels(z+1);
			for(int i=0;i<nTr;i++)echoValsReg[i]=(float[]) echoesReg[i].getStack().getPixels(z+1);
			double sigma=127;
			System.out.println("Sigma="+sigma);
			//for(int x=0;x<dimX;x++) {
			for(int x=0;x<dimX;x++) {
				System.out.println("zx="+z+","+x);
//				for(int y=0;y<dimY;y++) {
				for(int y=0;y<dimY;y++) {
					int offset=(y*dimX+x);
					for(int tr=0;tr<nTr;tr++) {tabData[tr]=echoVals[tr][offset];}
					if(VitimageUtils.mean(tabData)>3*sigma) {
						double[]params=(double[]) MRUtils.makeFit(tabTr, tabTe, tabData, MRUtils.T1_MONO_RICE, MRUtils.SIMPLEX, 400, 0, false)[0];
						mapVals[0][y*dimX+x]=(float) params[0];
						mapVals[1][y*dimX+x]=(float) params[1];
					}
					for(int tr=0;tr<nTr;tr++)tabDataReg[tr]=echoValsReg[tr][offset];
					if(VitimageUtils.mean(tabDataReg)>3*sigma) {
						double []params=(double[]) MRUtils.makeFit(tabTr, tabTe, tabDataReg, MRUtils.T1_MONO_RICE, MRUtils.SIMPLEX, 400, 0, false)[0];
						mapValsReg[0][y*dimX+x]=(float) params[0];
						mapValsReg[1][y*dimX+x]=(float) params[1];
					}
				}
				
			}
		}
		ImagePlus map=VitimageUtils.hyperStackingChannels(maps);
		map.setC(1);map.setDisplayRange(200, 3500);		IJ.run(map,"Fire","");
		map.setC(2);map.setDisplayRange(200, 3500);		IJ.run(map,"Fire","");
		IJ.saveAsTiff(map,dirExp01+"/Results/"+("RealData_Maps_estimation.tif"));

		ImagePlus mapReg=VitimageUtils.hyperStackingChannels(mapsReg);
		mapReg.setC(1);mapReg.setDisplayRange(200, 3500);		IJ.run(mapReg,"Fire","");
		mapReg.setC(2);mapReg.setDisplayRange(200, 3500);		IJ.run(mapReg,"Fire","");
		IJ.saveAsTiff(mapReg,dirExp01+"/Results/"+("RealData_Maps_estimation_after_registration.tif"));
	}


	
	
	public static int doudou(double d) {return (int)Math.round(d);}

	public static void writeStringTabInExcelFile(String[][]tab,String fileName) {
		System.out.println("Impression de tableau de taille "+tab.length+" x "+tab[0].length);
		try { 
			PrintStream l_out = new PrintStream(new FileOutputStream(fileName)); 
			for(int i=0;i<tab.length;i++) {
				for(int j=0;j<tab[i].length;j++) {
					l_out.print(tab[i][j]+" ;"); 
					System.out.print(tab[i][j]+" ;");
				}
				l_out.println(""); 
				System.out.println();
			}
			l_out.flush(); 
			l_out.close(); 
			l_out=null; 
		} 
		catch(Exception e){System.out.println(e.toString());} 
	}

	
	public static void writeDoubleTabInExcelFile(double[][]tab,String fileName) {
		System.out.println("Impression de tableau de taille "+tab.length+" x "+tab[0].length);
		try { 
			PrintStream l_out = new PrintStream(new FileOutputStream(fileName)); 
			for(int i=0;i<tab.length;i++) {
				for(int j=0;j<tab[i].length;j++) {
					l_out.print(tab[i][j]+" ;"); 
					//System.out.print(tab[i][j]+" ;");
				}
				l_out.println(""); 
				//System.out.println();
			}
			l_out.flush(); 
			l_out.close(); 
			l_out=null; 
		} 
		catch(Exception e){System.out.println(e.toString());} 
	}
	
	
	
	/** Functions for exp02 in chrononological order of calling*/////////////////////////////////////////////////////////////////////////////////
	public static double[]exp_02_getValsSigmaSimulation(){
		double []valsSNR=new double[] {200,500,1000};
		double []valsSigma=new double[valsSNR.length];
		for(int i=0;i<valsSNR.length;i++) {
			valsSigma[i]=Math.sqrt(PD_exp02*PD_exp02/valsSNR[i]);
		}
		return valsSigma;
	}

	public static ImagePlus exp_02_generateMaps() {
		int dimZ=dimZ_exp_02;
		ImagePlus[]maps=new ImagePlus[4];
		ImagePlus[]mapsComputed=new ImagePlus[4];
		int[] spacing=spacingSimulations_exp_02;
		int spaceX=spacing[0];int lx=spacing[1];
		int spaceY=spacing[2];int ly=spacing[3];
		double valT1=valT1_exp_02;
		double []valsT2=new double[] {20,40,60,80};
		double []valsSigma=exp_02_getValsSigmaSimulation();
		int dimX=(spaceX+lx)*(valsSigma.length)+spaceX;
		int dimY=(spaceY+ly)*(valsT2.length)+spaceY;
		ImagePlus img=IJ.createImage("DP", dimX, dimY, dimZ, 32);
		maps[0]=VitimageUtils.nullImage(img);
		maps[1]=VitimageUtils.nullImage(img);
		maps[2]=VitimageUtils.nullImage(img);
		maps[3]=VitimageUtils.nullImage(img);
		mapsComputed[0]=VitimageUtils.nullImage(img);
		mapsComputed[1]=VitimageUtils.nullImage(img);
		mapsComputed[2]=VitimageUtils.nullImage(img);
		mapsComputed[3]=VitimageUtils.nullImage(img);
		float[][]mapVals=new float[4][];
		float[][]mapValsComputed=new float[4][];
		for(int z=0;z<dimZ;z++) {
			mapVals[0]=(float[]) maps[0].getStack().getPixels(z+1);
			mapVals[1]=(float[]) maps[1].getStack().getPixels(z+1);
			mapVals[2]=(float[]) maps[2].getStack().getPixels(z+1);
			mapVals[3]=(float[]) maps[3].getStack().getPixels(z+1);
			mapValsComputed[0]=(float[]) mapsComputed[0].getStack().getPixels(z+1);
			mapValsComputed[1]=(float[]) mapsComputed[1].getStack().getPixels(z+1);
			mapValsComputed[2]=(float[]) mapsComputed[2].getStack().getPixels(z+1);
			mapValsComputed[3]=(float[]) mapsComputed[3].getStack().getPixels(z+1);
			for(int i=0;i<valsT2.length;i++) {
				for(int j=0;j<valsSigma.length;j++) {
					System.out.println("i,j="+i+","+j);
					for(int dx=0;dx<lx;dx++) {
						for(int dy=0;dy<ly;dy++) {
							int x=dx+j*(lx+spaceX)+spaceX;
							int y=dy+i*(ly+spaceY)+spaceY;
							mapVals[0][y*dimX+x]=(float) PD_exp02;
							mapVals[1][y*dimX+x]=(float) valT1;
							mapVals[2][y*dimX+x]=(float) valsT2[i];
							mapVals[3][y*dimX+x]=(float) valsSigma[j];
						}
					}
				}
			}
		}
		ImagePlus hyper=VitimageUtils.hyperStackingChannels(maps);
		hyper.setC(1);hyper.setDisplayRange(0, 1500);		IJ.run(hyper,"Fire","");
		hyper.setC(2);hyper.setDisplayRange(0, 800);		IJ.run(hyper,"Fire","");
		hyper.setC(3);hyper.setDisplayRange(0, 120);		IJ.run(hyper,"Fire","");
		//IJ.save(hyper,dirExp02+"/Source_data/maps.tif");
		//IJ.save(hyper,dirExp02+"/Results/Method_0_Source_values.tif");
		IJ.save(hyper,"home/fernandr/Bureau/FijiRelax_DOI/Experiments_for_figures_and_tables/Input_data/Images/Squared_fantoms/initial_maps.tif");
		hyper.show();
		return hyper;
	}
	
	public static void exp_02_generateSpinEchoT2Mono() {
		ImagePlus hyper=IJ.openImage(dirExp02+"/Source_data/maps.tif");
		int nRepet=nRepet_exp_02;
		int nEchoes=nEchoes_exp_02;
		double Tr=tr_exp_02;
		double TeSpacing=teSpacing_exp_02;
		ImagePlus[]tab=new ImagePlus[nEchoes];		
		ImagePlus[]maps=VitimageUtils.stacksFromHyperstackFastBis(hyper);
		for(int i=0;i<tab.length;i++)tab[i]=VitimageUtils.nullImage(maps[0]);
		float[][]echoVals=new float[nEchoes][];
		float[][]mapVals=new float[4][];
		int dimZ=tab[0].getNSlices();
		int dimX=tab[0].getWidth();
		int dimY=tab[0].getHeight();
		for(int z=0;z<dimZ;z++) {
			for(int i=0;i<maps.length;i++) {
				mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
			}
			for(int i=0;i<tab.length;i++) {
				echoVals[i]=(float[]) tab[i].getStack().getPixels(z+1);
				for(int x=0;x<dimX;x++) {
					for(int y=0;y<dimY;y++) {
						int pix=y*dimX+x;
						double acc=0;
						for(int r=0;r<nRepet;r++)acc+=RiceEstimator.getRandomRiceRealization(
							MRUtils.getFitFunctionValue(Tr, TeSpacing*(i+1),new double[] {mapVals[0][pix],mapVals[1][pix],mapVals[2][pix]}, 0, MRUtils.T1T2_MONO_RICE),mapVals[3][pix]);

						echoVals[i][pix]=(float)(acc/nRepet); 
					}
				}
			}
		}
		for(int i=0;i<tab.length;i++) {
			tab[i].setTitle(""+i);
			tab[i].setDisplayRange(0, 800);
			IJ.run(tab[i],"Fire","");
			//tab[i].show();
		}
		ImagePlus hyper2=VitimageUtils.hyperStackingChannels(tab);
		hyper2.setDisplayRange(0, 800);
		IJ.run(hyper,"Fire","");
		IJ.save(hyper2,dirExp02+"/Source_data/Simulated_echoes.tif");
	}
		
	public static void exp_02_estimateMapsFromEchoes(boolean realData,int fitType) {
		int nEchoes=nEchoes_exp_02;
		double Tr=tr_exp_02;
		double TeSpacing=teSpacing_exp_02;
		double[]tabTrTimes=new double[nEchoes];
		double[]tabTeTimes=new double[nEchoes];

		ImagePlus []imgs=new ImagePlus[nEchoes];
		
		if(realData) {
			for(int ec=0;ec<16;ec++) {
				System.out.println(dirExp02+"/Source_data/SSM_1_J0/TE_000"+(ec<10 ? "0" : "")+(11*(ec+1))+"/MAG_20190930T133721_SL000001.dcm");
				imgs[ec]=VitimageUtils.convertShortToFloatWithoutDynamicChanges(IJ.openImage(dirExp02+"/Source_data/SSM_1_J0/TE000"+
				(ec<9 ? "0" : "")+(11*(ec+1))+"/MAG_20190930T133721_SL000001.dcm"));
				//imgs[ec]=VitimageUtils.convertShortToFloatWithoutDynamicChanges(IJ.openImage(dirExp02+"/Source_data/B031_NP_J0/Echo_"+(ec+1)+".tif"));
				//imgs[ec]=VitimageUtils.cropFloatImage(imgs[ec], 0, 511, 0, 511, 0, 0);
			}
		}
		else {
			ImagePlus img=IJ.openImage(dirExp02+"/Source_data/Simulated_echoes.tif");
			imgs=VitimageUtils.stacksFromHyperstackFastBis(img);
		}
		System.out.println("Here");
		
		ImagePlus []maps=new ImagePlus[3];
		double[]tabData=new double[nEchoes];
		for(int ec=0;ec<nEchoes;ec++) {
			tabTrTimes[ec]=Tr;
			tabTeTimes[ec]=TeSpacing*(ec+1);
		}
		int dimZ=realData ? 1 : dimZ_exp_02;
		int[] spacing=spacingSimulations_exp_02;
		int spaceX=spacing[0];int lx=spacing[1];
		int spaceY=spacing[2];int ly=spacing[3];
		double valT1=valT1_exp_02;
		double []valsT2=valsT2Simulation_exp_02;
		double []valsSigma=null;
		if(realData) valsSigma=new double[] {130};//183.5
		else valsSigma=exp_02_getValsSigmaSimulation();
		int dimX=(spaceX+lx)*(valsSigma.length)+spaceX;
		maps[0]=VitimageUtils.nullImage(imgs[0]);
		maps[1]=VitimageUtils.nullImage(imgs[0]);
		maps[2]=VitimageUtils.nullImage(imgs[0]);
		float[][]mapVals=new float[3][];
		float[][]ecVals=new float[nEchoes][];
		for(int z=0;z<dimZ;z++) {
			for(int i=0;i<3;i++)mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
			for(int i=0;i<nEchoes;i++)ecVals[i]=(float[]) imgs[i].getStack().getPixels(z+1);

			if(realData) {
				dimX=maps[0].getWidth();
				int dimY=maps[0].getHeight();
				for(int x=0;x<dimX;x++) {
					System.out.println("x"+x);
					for(int y=0;y<dimY;y++) {
						int offset=y*dimX+x;
						for(int ec=0;ec<nEchoes;ec++)tabData[ec]=ecVals[ec][offset];
						double[]params=(double[]) MRUtils.makeFit(tabTrTimes, tabTeTimes, tabData, fitType, MRUtils.SIMPLEX, 400, valsSigma[0], false)[0];
						if(VitimageUtils.mean(tabData)>3*valsSigma[0]) {
							mapVals[0][y*dimX+x]=(float) params[0];
							mapVals[1][y*dimX+x]=(float) valT1;
							mapVals[2][y*dimX+x]=(float) params[1];
						}
					
					}
				}
				
			}
			else{
				for(int i=0;i<valsT2.length;i++) {
					System.out.println("i"+i);
					for(int j=0;j<valsSigma.length;j++) {
						for(int dx=0;dx<lx;dx++) {
							for(int dy=0;dy<ly;dy++) {
								int x=dx+j*(lx+spaceX)+spaceX;
								int y=dy+i*(ly+spaceY)+spaceY;
								for(int ec=0;ec<nEchoes;ec++)tabData[ec]=ecVals[ec][y*dimX+x];
								double[]params=(double[]) MRUtils.makeFit(tabTrTimes, tabTeTimes, tabData, fitType, MRUtils.SIMPLEX, 400, valsSigma[j], false)[0];
								mapVals[0][y*dimX+x]=(float) params[0];
								mapVals[1][y*dimX+x]=(float) valT1;
								mapVals[2][y*dimX+x]=(float) params[1];
							}
						}
					}
				}
			}
		}

	
		
		ImagePlus map=VitimageUtils.hyperStackingChannels(maps);
		map.setC(1);map.setDisplayRange(0, realData ? 8000 : 1500);		IJ.run(map,"Fire","");
		map.setC(2);map.setDisplayRange(0, realData ? 3000 : 800);		IJ.run(map,"Fire","");
		map.setC(3);map.setDisplayRange(0, realData ? 140 : 120);		IJ.run(map,"Fire","");
		String dat=realData ? "RealData":"Simulated";
		if (fitType==MRUtils.T2_MONO) 	IJ.saveAsTiff(map,dirExp02+"/Results/"+(dat+"_Method_1_Exp.tif"));
		if (fitType==MRUtils.T2_MONO_BIAS) 	IJ.saveAsTiff(map,dirExp02+"/Results/"+(dat+"_Method_2_Bias.tif"));
		if (fitType==MRUtils.T2_MONO_RICE) 	IJ.saveAsTiff(map,dirExp02+"/Results/"+(dat+"_Method_3_Rice.tif"));
	}

	
	
	
	
	public static ImagePlus estimateMapsFromEchoesSquaredFantoms(ImagePlus imgSimulatedEchoes,int fitType) {
		VitimageUtils.printImageResume(imgSimulatedEchoes);
		int nEchoes=nEchoes_exp_02;
		double Tr=tr_exp_02;
		double TeSpacing=teSpacing_exp_02;
		double[]tabTrTimes=new double[nEchoes];
		double[]tabTeTimes=new double[nEchoes];

		ImagePlus []imgs=new ImagePlus[nEchoes];
		
		ImagePlus img=imgSimulatedEchoes;
		imgs=VitimageUtils.stacksFromHyperstackFastBis(img);
		
		
		ImagePlus []maps=new ImagePlus[3];
		double[]tabData=new double[nEchoes];
		for(int ec=0;ec<nEchoes;ec++) {
			tabTrTimes[ec]=Tr;
			tabTeTimes[ec]=TeSpacing*(ec+1);
		}
		int dimZ= dimZ_exp_02;
		int[] spacing=spacingSimulations_exp_02;
		int spaceX=spacing[0];int lx=spacing[1];
		int spaceY=spacing[2];int ly=spacing[3];
		double valT1=valT1_exp_02;
		double []valsT2=valsT2Simulation_exp_02;
		double []valsSigma=null;
		valsSigma=exp_02_getValsSigmaSimulation();
		int dimX=(spaceX+lx)*(valsSigma.length)+spaceX;
		maps[0]=VitimageUtils.nullImage(imgs[0]);
		maps[1]=VitimageUtils.nullImage(imgs[0]);
		maps[2]=VitimageUtils.nullImage(imgs[0]);
		float[][]mapVals=new float[3][];
		float[][]ecVals=new float[nEchoes][];
		for(int z=0;z<dimZ;z++) {
			for(int i=0;i<3;i++)mapVals[i]=(float[]) maps[i].getStack().getPixels(z+1);
			for(int i=0;i<nEchoes;i++) {
				System.out.println("i="+i);
				ecVals[i]=(float[]) imgs[i].getStack().getPixels(z+1);
			}

			for(int i=0;i<valsT2.length;i++) {
				System.out.println("i"+i);
				for(int j=0;j<valsSigma.length;j++) {
					for(int dx=0;dx<lx;dx++) {
						for(int dy=0;dy<ly;dy++) {
							int x=dx+j*(lx+spaceX)+spaceX;
							int y=dy+i*(ly+spaceY)+spaceY;
							for(int ec=0;ec<nEchoes;ec++)tabData[ec]=ecVals[ec][y*dimX+x];
							double[]params=(double[]) MRUtils.makeFit(tabTrTimes, tabTeTimes, tabData, fitType, MRUtils.SIMPLEX, 400, valsSigma[j], false)[0];
							mapVals[0][y*dimX+x]=(float) params[0];
							mapVals[1][y*dimX+x]=(float) valT1;
							mapVals[2][y*dimX+x]=(float) params[1];
						}
					}
				}
			
			}
		}

	
		
		ImagePlus map=VitimageUtils.hyperStackingChannels(maps);
		map.setC(1);map.setDisplayRange(0,  1500);		IJ.run(map,"Fire","");
		map.setC(2);map.setDisplayRange(0,  800);		IJ.run(map,"Fire","");
		map.setC(3);map.setDisplayRange(0,  120);		IJ.run(map,"Fire","");
		String dat="Simulated";
		return map;
	}


	
	
	
	
	public static void processResultsSquaredFantoms(ImagePlus[]dataToProcess,String pathToCSV) {
		//Get the areas
		String[][]tabResults=new String[5][10];
		tabResults[0]=new String[] {"EXP 200","EXP 500","EXP 1000","OFF 200","OFF 500","OFF 1000","RICE 200","RICE 500","RICE 1000"};
		tabResults[1][0]="T2=20 ms";
		tabResults[2][0]="T2=40 ms";
		tabResults[3][0]="T2=60 ms";
		tabResults[4][0]="T2=80 ms";
		int[] spacing=spacingSimulations_exp_02;
		int spaceX=spacing[0];int lx=spacing[1];
		int spaceY=spacing[2];int ly=spacing[3];
		double []valsT2=valsT2Simulation_exp_02;
		double []valsSigma=exp_02_getValsSigmaSimulation();
		int dimX=(spaceX+lx)*(valsSigma.length)+spaceX;
		int dimY=(spaceY+ly)*(valsT2.length)+spaceY;
		int dimZ=dimZ_exp_02;
		double[][][][] dataCubeReal=new double[2][valsT2.length][valsSigma.length][lx*ly*dimZ];// {M0, T2} ; { 
		String[]testNames=new String[] {"Simulated_Method_0_Source_values.tif","Simulated_Method_1_Exp.tif","Simulated_Method_2_Bias.tif","Simulated_Method_3_Rice.tif"};
		
		for(int te=0;te<testNames.length;te++) {
			System.out.println("\n\n\n\n\nProcessing "+testNames[te]);
			ImagePlus img=dataToProcess[te];
			float[][]valT2=new float[dimZ][];
			for(int z=0;z<dimZ;z++) {
				valT2[z]=(float[]) img.getStack().getPixels(z+1);
			}
			for(int i=0;i<valsT2.length;i++) {
				for(int j=0;j<valsSigma.length;j++) {
					for(int z=0;z<dimZ;z++) {
						for(int dx=0;dx<lx;dx++) {
							for(int dy=0;dy<ly;dy++) {
								int x=dx+j*(lx+spaceX)+spaceX;
								int y=dy+i*(ly+spaceY)+spaceY;
								dataCubeReal[1][i][j][dx*ly+dy+z*lx*ly]=valT2[z][y*dimX+x];
								valT2[z][y*dimX+x]=(float) ((valT2[z][y*dimX+x]/*-valsT2[i]*/)/valsT2[i]);
							}
						}
					}
				}
			}
			for(int i=0;i<valsT2.length;i++) {
				System.out.println("");
				System.out.println("------ "+valsT2[i]+" ms ------");
				for(int j=0;j<valsSigma.length;j++) {
					System.out.println("Stat over "+dataCubeReal[0][i][j].length);
					double[]statsM0=VitimageUtils.statistics1D(dataCubeReal[0][i][j]);
					double[]statsT2=VitimageUtils.statistics1D(dataCubeReal[1][i][j]);
					System.out.print(" ["+dou(100*((statsT2[0]-valsT2[i])/valsT2[i]))+"+-"+dou(100*statsT2[1]/valsT2[i])+"]");
					if(te>0)tabResults[i+1][3*(te-1)+j+1]=""+dou(100*((statsT2[0]-valsT2[i])/valsT2[i]))+"+-"+dou(100*statsT2[1]/valsT2[i]);
				}
			}
		}	
		
		writeStringTabInExcelFile(tabResults,pathToCSV);
	}
	
	
	public static void processResultsSquaredFantomsOld(ImagePlus[]dataToProcess) {
		//Get the areas
		int[] spacing=spacingSimulations_exp_02;
		int spaceX=spacing[0];int lx=spacing[1];
		int spaceY=spacing[2];int ly=spacing[3];
		double []valsT2=valsT2Simulation_exp_02;
		double []valsSigma=exp_02_getValsSigmaSimulation();
		int dimX=(spaceX+lx)*(valsSigma.length)+spaceX;
		int dimY=(spaceY+ly)*(valsT2.length)+spaceY;
		int dimZ=dimZ_exp_02;
		double[][][][] dataCubeReal=new double[2][valsT2.length][valsSigma.length][lx*ly*dimZ];// {M0, T2} ; { 
		String[]testNames=new String[] {"Simulated_Method_0_Source_values.tif","Simulated_Method_1_Exp.tif","Simulated_Method_2_Bias.tif","Simulated_Method_3_Rice.tif"};
		
		for(int te=0;te<testNames.length;te++) {
			System.out.println("\n\n\n\n\nProcessing "+testNames[te]);
			ImagePlus img=dataToProcess[te];
			ImagePlus []imgs=VitimageUtils.stacksFromHyperstackFastBis(img);
			imgs[0].show();
			imgs[2].show();
			float[][]valM0=new float[dimZ][];
			float[][]valT2=new float[dimZ][];
			for(int z=0;z<dimZ;z++) {
				valM0[z]=(float[]) imgs[0].getStack().getPixels(z+1);
				valT2[z]=(float[]) imgs[2].getStack().getPixels(z+1);
			}
			for(int i=0;i<valsT2.length;i++) {
				for(int j=0;j<valsSigma.length;j++) {
					for(int z=0;z<dimZ;z++) {
						for(int dx=0;dx<lx;dx++) {
							for(int dy=0;dy<ly;dy++) {
								int x=dx+j*(lx+spaceX)+spaceX;
								int y=dy+i*(ly+spaceY)+spaceY;
								dataCubeReal[0][i][j][dx*ly+dy+z*lx*ly]=valM0[z][y*dimX+x];
								dataCubeReal[1][i][j][dx*ly+dy+z*lx*ly]=valT2[z][y*dimX+x];
								valM0[z][y*dimX+x]=(float) ((valM0[z][y*dimX+x]/*-PD*/)/PD_exp02);
								valT2[z][y*dimX+x]=(float) ((valT2[z][y*dimX+x]/*-valsT2[i]*/)/valsT2[i]);
							}
						}
					}
				}
			}
			img=VitimageUtils.hyperStackingChannels(imgs);
			img.setC(1);IJ.run(img,"Fire","");img.setDisplayRange(-0.5, 0.5);
			img.setC(3);IJ.run(img,"Fire","");img.setDisplayRange(0.8, 1.2);
			//IJ.saveAsTiff(img, dirExp02+"/Results/DIFF_"+testNames[te]);

			for(int i=0;i<valsT2.length;i++) {
				System.out.println("");
				System.out.println("------ "+valsT2[i]+" ms ------");
				for(int j=0;j<valsSigma.length;j++) {
					double[]statsM0=VitimageUtils.statistics1D(dataCubeReal[0][i][j]);
					double[]statsT2=VitimageUtils.statistics1D(dataCubeReal[1][i][j]);
					System.out.print(" ["+dou(100*((statsM0[0]-PD_exp02)/PD_exp02))+"+-"+dou(100*statsM0[1]/PD_exp02)+" - "+dou(100*((statsT2[0]-valsT2[i])/valsT2[i]))+"+-"+dou(100*statsT2[1]/valsT2[i])+"]");
				}
				System.out.println();
				
			}
		}		
	}
		

	public static ImagePlus[]stackToSlicesTframes(ImagePlus img){
		VitimageUtils.printImageResume(img);
		int X=img.getWidth();int Y=img.getHeight();int Z=img.getNFrames();
		ImagePlus []ret=new ImagePlus[Z];
		for(int z=0;z<Z;z++) {
			ret[z]=new Duplicator().run(img, 1, 1, 1, 1, z+1, z+1);
			VitimageUtils.adjustImageCalibration(ret[z],img);
			ret[z].getStack().setSliceLabel(img.getStack().getSliceLabel(z+1),1);
		}
		return ret;
	}	

	
	/**Functions for exp03*/////////////////////:
	public static double[][]exp_03_getValsSigmaSimulation(){
		double []valsSNR=new double[] {200};
		double []valsSigma=new double[valsSNR.length];
		for(int i=0;i<valsSNR.length;i++) {
			valsSigma[i]=Math.sqrt(PD_exp03*PD_exp03/valsSNR[i]);
			System.out.println("Effectivement :");
			System.out.println("PD="+PD_exp03);
			System.out.println("valsSNR[i]="+valsSNR[i]);
			System.out.println("valsSigma[i]="+valsSigma[i]);
			
		}
		return new double[][] {valsSigma,valsSNR};
	}

	public static void exp_03_generateMaps() {
		int[] spacing=spacingSimulations_exp_03;
		ImagePlus imgPD=IJ.createImage("PD", spacing[0], spacing[1], dimZ_exp_03, 32);
		imgPD=VitimageUtils.makeOperationOnOneImage(imgPD,1, PD_exp03,true);

		ImagePlus imgT1=IJ.createImage("T1", spacing[0], spacing[1], dimZ_exp_03, 32);
		imgT1=VitimageUtils.makeOperationOnOneImage(imgT1,1, T1_exp03,true);
		VitimageUtils.set32bitToValue(imgT1, T1_exp03);

		ImagePlus imgT2=IJ.createImage("T2", spacing[0], spacing[1], dimZ_exp_03, 32);
		imgT2=VitimageUtils.makeOperationOnOneImage(imgT2,1, T2_exp03,true);
		VitimageUtils.set32bitToValue(imgT2, T2_exp03);

		ImagePlus[]maps=new ImagePlus[] {imgPD.duplicate(),imgT1.duplicate(),imgT2.duplicate()};
		ImagePlus hyper=VitimageUtils.hyperStackingChannels(maps);
		IJ.save(hyper,dirExp03+"/Source_data/mapsSimT1T2SEQ.tif");
	}
	
	public static double[][][]exp_03_TrTe(){
		double[][][]tab=new double[3][][];
		tab[0]=new double[][] {{600.0,1200.0,2400.0},{10.0,10.0,10.0}};
		tab[1]=new double[2][nEchoes_exp_03];
		for(int i=0;i<nEchoes_exp_03;i++) {
			tab[1][1][i]=10*(i+1);
			tab[1][0][i]=10000;
		}
		tab[2]=new double[2][nEchoes_exp_03+3];
		for(int i=0;i<2;i++) {
			for(int j=0;j<3;j++)tab[2][i][j]=tab[0][i][j];
			for(int j=3;j<3+16;j++)tab[2][i][j]=tab[1][i][j-3];
		}
		return tab;
	}
	
	public static MRDataType[]exp_03_MRDataTypes(){
		MRDataType[]tab=new MRDataType[nEchoes_exp_03+3];
		for(int i=0;i<tab.length;i++)tab[i]=(i<3 ? MRDataType.T1SEQ : MRDataType.T2SEQ);
		return tab;
	}
	
	public static void exp_03_simulate_data() {
		int nSig=10;
		double valT1=1000;
		double valT2=30;
		double valPD=1000;
		int sigSpacing=10;
		double[]sigma=new double[nSig];
		int N=10000;
		int nRepets=4;
		double [][]tabTrTeT1T2=exp_03_TrTe()[2];
		double [][]tabTrTeT1=exp_03_TrTe()[0];
		double [][]tabTrTeT2=exp_03_TrTe()[1];
		double[][][]resultSeparated=new double[nSig][3][N];
		double[][][]resultJoint=new double[nSig][3][N];
		double[][][]resultSepPy=new double[3][2][nSig];
		double[][][]resultJoiPy=new double[3][2][nSig];
		for(int s=0;s<10;s++) {
			System.out.println(s);
			double SNR=100*(s+1);
			sigma[s]=1000/Math.sqrt(SNR);
		
			for(int n=0;n<N;n++) {
				//Simulate T1 data
				double[]signalT1=new double[3];
				for(int i=0;i<3;i++)signalT1[i]=RiceEstimator.getRandomRiceRealization(MRUtils.getFitFunctionValue(tabTrTeT1[0][i],tabTrTeT1[1][i], new double[] {valPD,valT1,valT2}, 0, MRUtils.T1T2_MONO_RICE),sigma[s],nRepets);
				
				//Simulate T2 data
				double[]signalT2=new double[16];
				for(int i=0;i<16;i++)signalT2[i]=RiceEstimator.getRandomRiceRealization(MRUtils.getFitFunctionValue(tabTrTeT2[0][i],tabTrTeT2[1][i], new double[] {valPD,valT1,valT2}, 0, MRUtils.T1T2_MONO_RICE),sigma[s],nRepets);
		
				//Simulate T1T2 data
				double[]signalT1T2=new double[19];
				for(int i=0;i<3;i++)signalT1T2[i]=signalT1[i];
				for(int i=0;i<16;i++)signalT1T2[i+3]=signalT2[i];
				
				double[]params=(double[])MRUtils.makeFit(tabTrTeT1[0], tabTrTeT1[1], signalT1, MRUtils.T1_MONO_RICE, MRUtils.SIMPLEX,400, sigma[s], false)[0];
				double T1sep=params[1];
				params=(double[])MRUtils.makeFit(tabTrTeT2[0], tabTrTeT2[1], signalT2, MRUtils.T2_MONO_RICE, MRUtils.SIMPLEX,400, sigma[s], false)[0];
				double PDsep=params[0];
				double T2sep=params[1];
		
				params=(double[])MRUtils.makeFit(tabTrTeT1T2[0], tabTrTeT1T2[1], signalT1T2, MRUtils.T1T2_MONO_RICE, MRUtils.SIMPLEX,400, sigma[s], false)[0];
				double PDjoint=params[0];
				double T1joint=params[1];
				double T2joint=params[2];
				
				resultJoint[s][0][n]=PDjoint;
				resultJoint[s][1][n]=T1joint;
				resultJoint[s][2][n]=T2joint;
				resultSeparated[s][0][n]=PDsep;
				resultSeparated[s][1][n]=T1sep;
				resultSeparated[s][2][n]=T2sep;
			}
			for(int pr=0;pr<3;pr++) {
				
				double[]tab=VitimageUtils.statistics1D(resultJoint[s][pr]);
				resultJoiPy[pr][0][s]=tab[0];
				resultJoiPy[pr][1][s]=tab[1];
				tab=VitimageUtils.statistics1D(resultSeparated[s][pr]);
				resultSepPy[pr][0][s]=tab[0];
				resultSepPy[pr][1][s]=tab[1];
			}
		}
		String[] prms=new String[]{"PD","T1","T2"};
		for(int i=0;i<3;i++) {
			writeDoubleTabInExcelFile(resultJoiPy[i],dirExp03+"/Results/Tables/"+"table_Join_"+prms[i]+"_line3.csv");
			writeDoubleTabInExcelFile(resultSepPy[i],dirExp03+"/Results/Tables/"+"table_Sep_"+prms[i]+"_line3.csv");		
		}
		System.exit(0);
	}

	
	public static void exp_03_simulate_data_v2() {
		int nVals=10;
		double valT1=1000;
		double valT2=60;
		double valPD=1000;
		int sigSpacing=10;
		double sigma;
		int N=10000;
		int nRepets=4;
		double [][]tabTrTeT1T2=exp_03_TrTe()[2];
		double [][]tabTrTeT1=exp_03_TrTe()[0];
		double [][]tabTrTeT2=exp_03_TrTe()[1];
		double[][][][]resultSeparated=new double[3][3][nVals][N];
		double[][][][]resultJoint=new double[3][3][nVals][N];
		double[][][][]resultSepPy=new double[3][3][nVals][N];
		double[][][][]resultJoiPy=new double[3][3][nVals][N];
		double[][]T1vals=new double[3][nVals];
		double[][]SNRvals=new double[3][nVals];
		for(int i=0;i<10;i++) {
			T1vals[0][i]=500+200*i;
			SNRvals[0][i]=200;
			T1vals[1][i]=1000;
			SNRvals[1][i]=100*(i+1);
			T1vals[2][i]=1500+200*i;;
			SNRvals[2][i]=100*(i+1);
		}
		for(int exp=0;exp<3;exp++) {
			for(int s=0;s<10;s++) {
				System.out.println(s);
				double SNR=SNRvals[exp][s];
				sigma=1000/Math.sqrt(SNR);
		
				for(int n=0;n<N;n++) {
					//Simulate T1 data
					double[]signalT1=new double[3];
					for(int i=0;i<3;i++)signalT1[i]=RiceEstimator.getRandomRiceRealization(MRUtils.getFitFunctionValue(tabTrTeT1[0][i],tabTrTeT1[1][i], new double[] {valPD,T1vals[exp][s],valT2}, 0, MRUtils.T1T2_MONO_RICE),sigma,nRepets);
					
					//Simulate T2 data
					double[]signalT2=new double[16];
					for(int i=0;i<16;i++)signalT2[i]=RiceEstimator.getRandomRiceRealization(MRUtils.getFitFunctionValue(tabTrTeT2[0][i],tabTrTeT2[1][i], new double[] {valPD,T1vals[exp][s],valT2}, 0, MRUtils.T1T2_MONO_RICE),sigma,nRepets);
			
					//Simulate T1T2 data
					double[]signalT1T2=new double[19];
					for(int i=0;i<3;i++)signalT1T2[i]=signalT1[i];
					for(int i=0;i<16;i++)signalT1T2[i+3]=signalT2[i];
					
					double[]params=(double[])MRUtils.makeFit(tabTrTeT1[0], tabTrTeT1[1], signalT1, MRUtils.T1_MONO_RICE, MRUtils.SIMPLEX,400, sigma, false)[0];
					double T1sep=100*(params[1]-T1vals[exp][s])/T1vals[exp][s];
					params=(double[])MRUtils.makeFit(tabTrTeT2[0], tabTrTeT2[1], signalT2, MRUtils.T2_MONO_RICE, MRUtils.SIMPLEX,400, sigma, false)[0];
					double PDsep=100*(params[0]-valPD)/valPD;
					double T2sep=100*(params[1]-valT2)/valT2;
			
					params=(double[])MRUtils.makeFit(tabTrTeT1T2[0], tabTrTeT1T2[1], signalT1T2, MRUtils.T1T2_MONO_RICE, MRUtils.SIMPLEX,400, sigma, false)[0];
					double PDjoint=100*(params[0]-valPD)/valPD;
					double T1joint=100*(params[1]-T1vals[exp][s])/T1vals[exp][s];
					double T2joint=100*(params[2]-valT2)/valT2;
					
					resultJoiPy[exp][0][s][n]=PDjoint;
					resultJoiPy[exp][1][s][n]=T1joint;
					resultJoiPy[exp][2][s][n]=T2joint;
					resultSepPy[exp][0][s][n]=PDsep;
					resultSepPy[exp][1][s][n]=T1sep;
					resultSepPy[exp][2][s][n]=T2sep;
				}
/*				for(int pr=0;pr<3;pr++) {		
					double[]tab=HyperMap.getQuartiles(resultJoint[exp][s][pr]);
					for(int quar=0;quar<3;quar++)resultJoiPy[exp][pr][quar][s]=tab[quar];
					tab=HyperMap.getQuartiles(resultSeparated[exp][s][pr]);
					for(int quar=0;quar<3;quar++)resultSepPy[exp][pr][quar][s]=tab[quar];
				}
*/
			}
		}
		String[] prms=new String[]{"PD","T1","T2"};
		String[] exps=new String[]{"fig4","fig5","fig3"};
		for(int exp=0;exp<3;exp++) {
			for(int i=0;i<3;i++) {
				writeDoubleTabInExcelFile(resultJoiPy[exp][i],dirExp03+"/Results/Tables/"+"table_Join_"+prms[i]+"_quartile_"+exps[exp]+".csv");
				writeDoubleTabInExcelFile(resultSepPy[exp][i],dirExp03+"/Results/Tables/"+"table_Sep_"+prms[i]+"_quartile_"+exps[exp]+".csv");
			}
		}
		System.exit(0);
	}


	
	public static void exp_03_generateSpinEchoData() {
		MRDataType []type=exp_03_MRDataTypes();
		double[][][]times=exp_03_TrTe();
		System.out.println(TransformUtils.stringVectorN(times[2][0], "Tr"));
		System.out.println(TransformUtils.stringVectorN(times[2][1], "Te"));
		double[][]sigmas=exp_03_getValsSigmaSimulation();
		ImagePlus maps=IJ.openImage(dirExp03+"/Source_data/mapsSim"+(MRDataType.T1T2SEQ)+".tif");
		for(int indSig=0;indSig<sigmas[0].length;indSig++) {
			ImagePlus hyper=MRUtils.simulateEchoesT1T2Relax(maps,times[2][0],times[2][1],nRepet_exp_03,type,sigmas[0][indSig],"Simulation_SNR"+sigmas[1][indSig]);
			IJ.saveAsTiff(hyper,dirExp03+"/Source_data/Simulation_SNR"+sigmas[1][indSig]+".tif");
		}
	}			
			
	public static void exp_03_estimateMapsFromEchoes() {
		double[][][]times=exp_03_TrTe();
		double[][]sigmas=exp_03_getValsSigmaSimulation();

		for(int indSig=0;indSig<sigmas[0].length;indSig++) {

			//Compute separated fits
			HyperMap hyper=new HyperMap(IJ.openImage(dirExp03+"/Source_data/Simulation_SNR"+sigmas[1][indSig]+".tif"));
			hyper.computeMapsAgainAndMask(MRUtils.SIMPLEX,true,NoiseManagement.RICE,false,null,4);
			IJ.saveAsTiff(hyper.hyperImg,dirExp03+"/Results/Simul_SNR"+sigmas[1][indSig]+"_separated.tif");

			//Compute cross fits
			hyper=new HyperMap(IJ.openImage(dirExp03+"/Source_data/Simulation_SNR"+sigmas[1][indSig]+".tif"));
			hyper.computeMapsAgainAndMask(MRUtils.SIMPLEX,false,NoiseManagement.RICE,false,null,4);
			IJ.saveAsTiff(hyper.hyperImg,dirExp03+"/Results/Simul_SNR"+sigmas[1][indSig]+"_crossfit.tif");
		}
	}

	public static void exp_03_processResults() {
		//Get the areas
		int[] spacing=spacingSimulations_exp_02;
		int spaceX=spacing[0];int lx=spacing[1];
		int spaceY=spacing[2];int ly=spacing[3];
		double []valsT2=valsT2Simulation_exp_02;
		double []valsSigma=exp_02_getValsSigmaSimulation();
		int dimX=(spaceX+lx)*(valsSigma.length)+spaceX;
		int dimY=(spaceY+ly)*(valsT2.length)+spaceY;
		int dimZ=dimZ_exp_02;
		double[][][][] dataCubeReal=new double[2][valsT2.length][valsSigma.length][lx*ly*dimZ];// {M0, T2} ; { 
		String[]testNames=new String[] {"Simulated_Method_0_Source_values.tif","Simulated_Method_1_Exp.tif","Simulated_Method_2_Bias.tif","Simulated_Method_3_Rice.tif"};
		
		for(int te=0;te<testNames.length;te++) {
			System.out.println("\n\n\n\n\nProcessing "+testNames[te]);
			ImagePlus img=IJ.openImage(dirExp02+"/Results/"+testNames[te]);
			ImagePlus []imgs=VitimageUtils.stacksFromHyperstackFastBis(img);
			imgs[0].show();
			imgs[2].show();
			float[][]valM0=new float[dimZ][];
			float[][]valT2=new float[dimZ][];
			for(int z=0;z<dimZ;z++) {
				valM0[z]=(float[]) imgs[0].getStack().getPixels(z+1);
				valT2[z]=(float[]) imgs[2].getStack().getPixels(z+1);
			}
			for(int i=0;i<valsT2.length;i++) {
				for(int j=0;j<valsSigma.length;j++) {
					for(int z=0;z<dimZ;z++) {
						for(int dx=0;dx<lx;dx++) {
							for(int dy=0;dy<ly;dy++) {
								int x=dx+j*(lx+spaceX)+spaceX;
								int y=dy+i*(ly+spaceY)+spaceY;
								dataCubeReal[0][i][j][dx*ly+dy+z*lx*ly]=valM0[z][y*dimX+x];
								dataCubeReal[1][i][j][dx*ly+dy+z*lx*ly]=valT2[z][y*dimX+x];
								valM0[z][y*dimX+x]=(float) ((valM0[z][y*dimX+x]/*-PD*/)/PD_exp02);
								valT2[z][y*dimX+x]=(float) ((valT2[z][y*dimX+x]/*-valsT2[i]*/)/valsT2[i]);
							}
						}
					}
				}
			}
			img=VitimageUtils.hyperStackingChannels(imgs);
			img.setC(1);IJ.run(img,"Fire","");img.setDisplayRange(-0.5, 0.5);
			img.setC(3);IJ.run(img,"Fire","");img.setDisplayRange(0.8, 1.2);
			IJ.saveAsTiff(img, dirExp02+"/Results/DIFF_"+testNames[te]);

			for(int i=0;i<valsT2.length;i++) {
				System.out.println("");
				System.out.println("------ "+valsT2[i]+" ms ------");
				for(int j=0;j<valsSigma.length;j++) {
					double[]statsM0=VitimageUtils.statistics1D(dataCubeReal[0][i][j]);
					double[]statsT2=VitimageUtils.statistics1D(dataCubeReal[1][i][j]);
					System.out.print(" ["+dou(100*((statsM0[0]-PD_exp02)/PD_exp02))+"+-"+dou(100*statsM0[1]/PD_exp02)+" - "+dou(100*((statsT2[0]-valsT2[i])/valsT2[i]))+"+-"+dou(100*statsT2[1]/valsT2[i])+"]");
				}
				System.out.println();
				
			}
		}		
	}
		
	public static void exp_03_registerEchoesRealData() {
		ImagePlus hyper=IJ.openImage(dirExp01+"/Source_data/RealData_echoes.tif");
		ImagePlus []echoes=VitimageUtils.stacksFromHyperstackFastBis(hyper);
		int nTr=echoes.length;
		ItkTransform []trs=new ItkTransform[nTr-1];
		
		for(int i=0;i<nTr-1;i++) {
			BlockMatchingRegistration bmRegistration;
			RegistrationAction regAct=new RegistrationAction();
			regAct.defineSettingsSimplyFromTwoImages(echoes[nTr-1], echoes[i]);
			regAct.typeAutoDisplay=2;
			regAct.higherAcc=0;
			regAct.levelMaxLinear=2;
			regAct.levelMinLinear=-1;
			System.out.println(regAct.readableString());
			bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(echoes[nTr-1], echoes[i],regAct);
			ItkTransform tr=bmRegistration.runBlockMatching(null,false);
			bmRegistration.closeLastImages();
		    echoes[i]=tr.transformImage(echoes[nTr-1], echoes[i]);
		    trs[i]=tr;
		}

		for(int tr=0;tr<nTr;tr++) {
			echoes[tr].setTitle(""+tr);
			echoes[tr].setDisplayRange(0, PDUp_exp_01);
			IJ.run(echoes[tr],"Fire","");
			if(tr<nTr-1)System.out.println("\nTransform "+tr+"\n"+trs[tr].drawableString());
		}
		hyper=VitimageUtils.hyperStackingChannels(echoes);
		hyper.setDisplayRange(0, PDUp_exp_01);
		IJ.run(hyper,"Fire","");
		IJ.save(hyper,dirExp01+"/Source_data/RealData_echoes_registered.tif");
	}

	
	public static void exp_03_compute_maps_for_real_data() {
		double[][][]times=exp_03_TrTe();
		double[][]sigmas=exp_03_getValsSigmaSimulation();
		String []sources=new String[]{"/home/fernandr/Bureau/FijiRelax/DATA_FIJIRELAX/Source/Dataset_1_Vine_cutting/J0_crop.tif",
				"/home/fernandr/Bureau/FijiRelax/DATA_FIJIRELAX/Source/Sorghum/SSM1_Day0_Registered/SSM1_0930_echoes_registered.tif"};
		String []names=new String[] {"Vine","Sorgho"};
		for(int s=1;s<sources.length;s++) {
			HyperMap hyper=new HyperMap(IJ.openImage(sources[s]));			
			hyper=hyper.pruneBionanoStuff();
			hyper.computeMapsAgainAndMask(MRUtils.SIMPLEX,true,NoiseManagement.RICE,false,null,4);
			System.out.println(dirExp03+"/Results/Sorgho_separated.tif");
			IJ.saveAsTiff(hyper.hyperImg,dirExp03+"/Results/Sorgho_separated.tif");
		
			hyper=new HyperMap(IJ.openImage(sources[s]));			
			hyper=hyper.pruneBionanoStuff();
			hyper.computeMapsAgainAndMask(MRUtils.SIMPLEX,false,NoiseManagement.RICE,false,null,5);
			IJ.saveAsTiff(hyper.hyperImg,dirExp03+"/Results/Sorgho_joint.tif");
		}
	}
		
	public static void exp_03_compute_maps_for_real_data2() {
		System.out.println("Here1");
		double[][][]times=exp_03_TrTe();
		double[][]sigmas=exp_03_getValsSigmaSimulation();
		System.out.println("Here2");
		String []sources=new String[]{"/home/fernandr/Bureau/FijiRelax/Exp_03_Joint_fit/Source_data/B098_J133.tif"};
		System.out.println("Here3");
		String []names=new String[] {"Vine"};
		System.out.println("Here4");
		for(int s=0;s<sources.length;s++) {
			System.out.println("Here5");
			HyperMap hyper=new HyperMap(IJ.openImage(sources[s]));			
			hyper.computeMapsAgainAndMask(MRUtils.SIMPLEX,false,NoiseManagement.RICE,false,null,4);
			System.out.println(dirExp03+"/Results/Vine_long_joint70.tif");
			IJ.saveAsTiff(hyper.hyperImg,dirExp03+"/Results/Vine_long_joint133.tif");			
		}
	}

	
	public static boolean hasNegativeValues(double[]d) {
		for(int i=0;i<d.length;i++)if(d[i]<0)return true;
		return false;
	}
		
	
	
	public static void testSpinEchoMagnitudesMonoExp(double PD,double T1, double T2,double Tr,double TeSpacing,int nEchoes,double sigmaRice) {
		double[]values=new double[nEchoes];
		for(int i=0;i<nEchoes;i++) {
			values[i]= RiceEstimator.getRandomRiceRealization(MRUtils.getFitFunctionValue(Tr, TeSpacing*(i+1),new double[] {PD,T1,T2}, 0, MRUtils.T1T2_MONO_RICE),sigmaRice);
			System.out.println(values[i]);
		}
	}

	public static void testRandomRice() {
		int N=20000;
		double originalSignal=1000;
		double sigmaRice=100;
		double[]tab=new double[N];
		for(int i=0;i<N;i++) {
			tab[i]=RiceEstimator.getRandomRiceRealization(originalSignal, sigmaRice);
			System.out.println(tab[i]);
		}
		System.out.println("Init, Mean, sigma="+originalSignal+", "+TransformUtils.stringVectorN(VitimageUtils.statistics1D(tab),""));
		System.out.println("GetInitvalue="+MRUtils.besFunkCost(originalSignal, sigmaRice));
	}

	public static double dou(double d) {
		double res=Math.round(d*100)/100.0;
		return res;
	}
	

}
