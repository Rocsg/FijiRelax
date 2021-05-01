package fr.cirad.image.fijirelax.mrialgo;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import fr.cirad.image.common.Timer;
import fr.cirad.image.common.TransformUtils;
import fr.cirad.image.common.VitimageUtils;
import fr.cirad.image.fijiyama.RegistrationAction;
import fr.cirad.image.registration.BlockMatchingRegistration;
import fr.cirad.image.registration.ItkTransform;
import fr.cirad.image.registration.Transform3DType;

import fr.cirad.image.fijirelax.gui.Custom_Format_Importer;
import fr.cirad.image.fijirelax.testing.ValidationExperimentsFijiRelaxPaper;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import no.uib.cipr.matrix.sparse.SSOR;


public class HyperMap {
	public static int defaultOutlierAlgorithm=1;//1=Tukey fences, 2=double-sided MADe
	public static double defaultOutlierStdDev=3;
	public static int defaultOutlierNeighbourXY=3;
	public static int defaultAlgoChoice=0;//0=Simplexe, 1=Levenberg
	public static int defaultNoiseChoice =0;//0=Rice , 1=Offset, 2=nothing
	public static int defaultFirstChoice =0;//0=Yes , 1=No (weird)
	public static int defaultMaskChoice =0;//0=Auto , 1=provided
	public static double defaultStdDevMask =4;//Threshold for significant points is = mean Rice noise + n * sigma rice noise
	public static int defaultJointSchemeChoice =0;//0=Joint , 1=Separated
	public static double factorViewNormalisation=1.2;
	public static double percentageKeepNormalisation=99;

	public int[]dims;
	public double[]voxs;
	public int T;
	public int C;
	public int X;
	public int Y;
	public int Z;
	public int Ctot;
	public int nMaps=0;
	final int m0MapChannel=0;
	final int t1MapChannel=1;
	final int t2mapChannel=2;
	public boolean hasT2sequence=false;
	public boolean hasT1sequence=false;
	public boolean hasT1T2sequence=false;
	public boolean hasMaps=false;
	public ImagePlus hyperImg;
	public ImagePlus hyperMaps;
	public ImagePlus hyperEchoes;
	public double[][][]Te;
	public double[][][]Tr;
	public double [][][]sigmaRice;
	public MRDataType [][] mrDataType;
	public String[]actualDay;
	public double[][]tabSigmasT1Seq;
	public double[][]tabSigmasT2Seq;
	public double[][]tabSigmasT1T2Seq;
	public String name="JOHN-DOE";
	public double sigmaInUse=0;
	private double defaultContrastSat=0.35;
/*	private double TR_ALPHA=0;
	private double TE_ALPHA=0;*/

	

	
	/** Main entry points : Main, Instance constructor and static factory ---------------------------------------------------------*/
	public HyperMap() {}

	public HyperMap(ImagePlus hyperImg) {
		this.hyperImg=hyperImg;
		setup();
		adjustContrast();
	}
	
	
	
	public String getRangeInfos() {
		String s="";
		if(this.hyperMaps!=null) {
			for(int c=0; c<this.hyperMaps.getNChannels();c++) {
				s+="Maps channel "+c+" = "+hyperMaps.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperMaps, c, 0, 0));
				int cTemp=hyperMaps.getC();
				hyperMaps.setC(c+1);
				s+=" Range=["+hyperMaps.getDisplayRangeMin()+" , "+hyperMaps.getDisplayRangeMax()+"]\n";
				hyperMaps.setC(cTemp);
			}
		}
		for(int c=0;(this.hyperEchoes!=null) && ((c<2) && (c<this.hyperEchoes.getNChannels()));c++) {
			int cTemp=hyperEchoes.getC();
			hyperEchoes.setC(c+1);
			//s+="\nEchoes channel "+c+" = "+hyperEchoes.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperEchoes, c, 0, 0))+"\n";
			//s+="Range=["+hyperEchoes.getDisplayRangeMin()+" , "+hyperEchoes.getDisplayRangeMax()+"]\n";
			hyperEchoes.setC(cTemp);
		}
		if(this.hyperMaps!=null) {
			for(int c=0; c<this.hyperMaps.getNChannels()+1;c++) {
				s+="HyperImg channel "+c+" = "+hyperImg.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, 0, 0));
				int cTemp=hyperImg.getC();
				hyperImg.setC(c+1);
				s+=" Range=["+hyperImg.getDisplayRangeMin()+" , "+hyperImg.getDisplayRangeMax()+"]\n";
				hyperImg.setC(cTemp);
			}
		}
		return s;
	}
	
	public static HyperMap hyperMapFactory(HyperMap hyp) {
		return new HyperMap(hyp.getAsImagePlus());
	}
	
	public static HyperMap importHyperMapFromRawDicomData(String inputDir,String nameSpecimen) {
		IJ.log("Starting importation with : \ninputDir="+inputDir+"\nnameSpecimen="+nameSpecimen);
		Custom_Format_Importer t1t2=new Custom_Format_Importer(inputDir,nameSpecimen);
		HyperMap hyp=t1t2.importHyperMap();
		ImagePlus img=hyp.getAsImagePlus();
		IJ.run(img,"Fire","");
		HyperMap hyperNew=new HyperMap(img);
		return hyp;
	}


	
	public static HyperMap importHyperMapFromCustomData(String path,String pattern) {
		int TR=1000;
		int TE=10;
		boolean hasTR=false;
		boolean hasTE=false;
		if( (! pattern.contains("{TE}")) && (!pattern.contains("{TR}")) ){
			IJ.showMessage("Wrong pattern for a T1 and/or T2 sequence : "+pattern);
			return null;
		}
		if (pattern.contains("{TE}"))hasTE=true;
		if (pattern.contains("{TR}"))hasTR=true;		
		System.out.println("Importing from pattern with TE ? "+hasTE+" with TR ?"+hasTR);
		
		//Build an arrayList with fileNames
		String[]files=new File(path).list();
		for(int i=0;i<files.length;i++) System.out.println(files[i]);
		System.out.println("Import from custom : arrayList ok");

		
		
		//Read corresponding TR / TE
		String sep="SEPARATOR";
		ArrayList<double[]>paramList=new ArrayList<double[]>();
		for(int i=0;i<files.length;i++) {
			System.out.print("At "+i+" : add ");
			if(hasTE && hasTR){
				//Guess that 1) File names does not begin or end by a parameter such as 000010imgtruc.tif or imgdsdsd.tif.000010
				//And guess that TR and TE params are not concatenated such as img00010050000.tif
				//And guess that file names does not contain "SEPARATOR"
				//And guess that file names does not repeat TR and TE values
				String pat0=pattern.replace("{TR}",sep).replace("{TE}",sep).split(sep)[0];
				String pat1=pattern.replace("{TR}",sep).replace("{TE}",sep).split(sep)[1];
				String pat2=pattern.replace("{TR}",sep).replace("{TE}",sep).split(sep)[2];
				System.out.println("Pat0="+pat0);
				System.out.println("Pat1="+pat1);
				System.out.println("Pat2="+pat2);
				String str=files[i].replace(pat0,sep).replace(pat1, sep).replace(pat2, sep);
				double val0=Double.parseDouble(str.split(sep)[1]);
				double val1=Double.parseDouble(str.split(sep)[2]);
				if(pattern.indexOf("{TR}")<pattern.indexOf("{TE}")){				
					paramList.add(new double[] {i,val0,val1});
					System.out.println(""+i+","+val0+","+val1);
				}
				else {
					paramList.add(new double[] {i,val1,val0});
					System.out.println(""+i+","+val1+","+val0);
				}
			}			
			else if(hasTE){
				//Guess that 1) File names does not begin or end by a parameter such as 000010imgtruc.tif or imgdsdsd.tif.000010
				//And guess that file names does not contain "SEPARATOR"
				//And guess that file names does not repeat TR and TE values
				String pat0=pattern.replace("{TR}",sep).replace("{TE}",sep).split(sep)[0];
				String pat1=pattern.replace("{TR}",sep).replace("{TE}",sep).split(sep)[1];
				System.out.println("Pat0="+pat0);
				System.out.println("Pat1="+pat1);
				String str=files[i].replace(pat0,sep).replace(pat1, sep);
				double val0=Double.parseDouble(str.split(sep)[1]);
				paramList.add(new double[] {i,TR,val0});
				System.out.println(""+i+","+TR+","+val0);
			}			
			else if(hasTR){
				//Guess that 1) File names does not begin or end by a parameter such as 000010imgtruc.tif or imgdsdsd.tif.000010
				//And guess that file names does not contain "SEPARATOR"
				//And guess that file names does not repeat TR and TE values
				String pat0=pattern.replace("{TR}",sep).replace("{TE}",sep).split(sep)[0];
				String pat1=pattern.replace("{TR}",sep).replace("{TE}",sep).split(sep)[1];
				System.out.println("Pat0="+pat0);
				System.out.println("Pat1="+pat1);
				String str=files[i].replace(pat0,sep).replace(pat1, sep);
				System.out.println(str);
				double val0=Double.parseDouble(str.split(sep)[1]);
				paramList.add(new double[] {i,val0,TE});
				System.out.println(""+i+","+val0+","+TE);
			}			
			else {
				IJ.showMessage("No TR / TE data here : "+files[i]);
			}			
		}
		System.out.println("Import from custom : paramList ok");
		for(int i=0;i<files.length;i++) {
			System.out.println(TransformUtils.stringVector(paramList.get(i)," List["+i+"]"));
		}
		
		//Sort the array
		Collections.sort(paramList, new Comparator<double[]>() {
			@Override
			public int compare(double[] o1, double[] o2) {
				return Double.compare(o1[1]*100000+o1[2],o2[1]*100000+o2[2]);
			}
		});
		System.out.println("Import from custom : sort ok");
		for(int i=0;i<files.length;i++) {
			System.out.println(TransformUtils.stringVector(paramList.get(i)," List["+i+"]"));
		}

		//Stack the data
		ImagePlus[]imgs=new ImagePlus[files.length];
		for(int i=0;i<files.length;i++) {
			imgs[i]=IJ.openImage(new File(path,files[(int)Math.round(paramList.get(i)[0]) ] ).getAbsolutePath() );
			if(! VitimageUtils.dimensionsMatch(imgs[i],imgs[0])){
				IJ.showMessage("Warning : dimensions does not match between "+files[(int)Math.round(paramList.get(0)[0])]+" and "+files[(int)Math.round(paramList.get(0)[0])]);
				return null;
			}
		}
		ImagePlus img=VitimageUtils.hyperStackingChannels(imgs);
		System.out.println("Import from custom : stacking ok "+VitimageUtils.imageResume(img));

		
		//Determine type of data
		boolean hasT1Seq=false;
		boolean hasT2Seq=false;
		boolean hasT1T2Seq=false;
		for(int i=0;i<files.length;i++) {
			if(paramList.get(i)[1]!=paramList.get(0)[1])hasT1Seq=true;
			if(paramList.get(i)[2]!=paramList.get(0)[2])hasT2Seq=true;
		}
		if(hasT1Seq && hasT2Seq)hasT1T2Seq=true;
		String seq=(hasT1T2Seq ? "T1T2SEQ" : (hasT1Seq ? "T1SEQ" : "T2SEQ"));
		System.out.println("Import from custom : type of data ok : "+seq);

		
		//Write seq parameters in metadata
		for(int c=0;c<img.getNChannels();c++) {
			for(int z=0;z<img.getNSlices();z++) {
				int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0);
				ImageProcessor imgP=img.getStack().getProcessor(sli);
				double[]stats=MRUtils.getBackgroundStatsFromProcessorTight(imgP);
				double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
				double tr=paramList.get(c)[1];
				double te=paramList.get(c)[2];
				String chain=seq+"_CUSTOM_TR="+VitimageUtils.dou(tr)+"_TE="+VitimageUtils.dou(te)+"_SIGMARICE="+VitimageUtils.dou(sigmaRice);
				img.getStack().setSliceLabel(chain,VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0) );
			}		
		}
		
		return  new HyperMap(img);
	}
	
	
	//TODO : use more features of nifti metadata (radiologic convention for example) in more handsome functions, in a specific class NiftiImporter
	/** Main entry points : Nifti handling ---------------------------------------------------------*/
	public static HyperMap importHyperMapFromNifti4DT2Sequence(String path,String imageName,double Tr, double TeSpacing) {
		ImagePlus img2=IJ.openImage(path);
		boolean isOldBioImage=false;
		if(img2==null) {
			isOldBioImage=true;
			img2=IJ.getImage();
		}
		IJ.run(img2, "Flip Vertically", "stack");
		
		if(isOldBioImage) {
			String s=img2.getTitle();
			img2.setTitle("Readen Nifti");
			ImagePlus tmp=ij.WindowManager.getImage(new File(path).getName());
			if(tmp != null)tmp.close();
			if(img2!=null)img2.hide();
		}
		
		ImagePlus[]tab=ValidationExperimentsFijiRelaxPaper.stackToSlicesTframes(img2);
		tab[0].setTitle("Readen tab0");
		ImagePlus img=VitimageUtils.hyperStackingChannels(tab);
		img.setTitle("Stacken Nifti");
		
		for(int c=0;c<img.getNChannels();c++) {
			for(int z=0;z<img.getNSlices();z++) {
				int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0);
				ImageProcessor imgP=img.getStack().getProcessor(sli);
				double[]stats=MRUtils.getBackgroundStatsFromProcessorTight(imgP);
				double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
				double tr=Tr;
				double te=TeSpacing*(c+1);
				String chain="T2SEQ_NIFTI_TR="+VitimageUtils.dou(tr)+"_TE="+VitimageUtils.dou(te)+"_SIGMARICE="+VitimageUtils.dou(sigmaRice);
				img.getStack().setSliceLabel(chain,VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0) );
			}		
		}
		
		HyperMap hyper=new HyperMap(img);
		hyper.adjustContrast();
		return hyper;
	}

	public static HyperMap importHyperMapFromNifti4DT2Sequence(String path,String imageName,double Tr, double Tes[]) {
		ImagePlus img2=IJ.openImage(path);
		img2.setTitle("Readen Nifti");
		IJ.run(img2, "Flip Vertically", "stack");
		ImagePlus[]tab=ValidationExperimentsFijiRelaxPaper.stackToSlicesTframes(img2);
		tab[0].setTitle("Readen tab0");
		ImagePlus img=VitimageUtils.hyperStackingChannels(tab);
		img.setTitle("Stacken Nifti");
		
		for(int c=0;c<img.getNChannels();c++) {
			for(int z=0;z<img.getNSlices();z++) {
				int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0);
				ImageProcessor imgP=img.getStack().getProcessor(sli);
				double[]stats=MRUtils.getBackgroundStatsFromProcessorTight(imgP);
				double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
				double tr=Tr;
				double te=Tes[c];
				String chain="T2SEQ_NIFTI_TR="+VitimageUtils.dou(tr)+"_TE="+VitimageUtils.dou(te)+"_SIGMARICE="+VitimageUtils.dou(sigmaRice);
				img.getStack().setSliceLabel(chain,VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0) );
			}		
		}
		
		HyperMap hyper=new HyperMap(img);
		hyper.adjustContrast();
		return hyper;
	}

	public static HyperMap importHyperMapFromNifti4DSequence(String path,String imageName,Object[]niftiInfos) {
		ImagePlus img2=IJ.openImage(path);
		IJ.run(img2, "Flip Vertically", "stack");
		img2.setTitle("Readen Nifti");
		ImagePlus[]tab=ValidationExperimentsFijiRelaxPaper.stackToSlicesTframes(img2);
		tab[0].setTitle("Readen tab0");
		ImagePlus img=VitimageUtils.hyperStackingChannels(tab);
		img.setTitle("Stacken Nifti");
		
		for(int c=0;c<img.getNChannels();c++) {
			for(int z=0;z<img.getNSlices();z++) {
				int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0);
				ImageProcessor imgP=img.getStack().getProcessor(sli);
				double[]stats=MRUtils.getBackgroundStatsFromProcessorTight(imgP);
				double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
				MRDataType mr=((MRDataType[])niftiInfos[0])[c];
				double tr=((double[])(niftiInfos[1]))[c];
				double te=((double[])(niftiInfos[2]))[c];
				String chain=""+mr+"_NIFTI_TR="+VitimageUtils.dou(tr)+"_TE="+VitimageUtils.dou(te)+"_SIGMARICE="+VitimageUtils.dou(sigmaRice);
				img.getStack().setSliceLabel(chain,VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0) );
			}		
		}
		
		HyperMap hyper=new HyperMap(img);
		hyper.adjustContrast();
		return hyper;
	}

	public static HyperMap importHyperMapFromNifti4DT1Sequence(String path,String imageName,double []Tr, double TeSpacing) {
		ImagePlus img2=IJ.openImage(path);
		img2.setTitle("Readen Nifti");
		IJ.run(img2, "Flip Vertically", "stack");
		ImagePlus[]tab=ValidationExperimentsFijiRelaxPaper.stackToSlicesTframes(img2);
		tab[0].setTitle("Readen tab0");
		ImagePlus img=VitimageUtils.hyperStackingChannels(tab);
		img.setTitle("Stacken Nifti");
		if(img.getNChannels()!=Tr.length) {
			IJ.showMessage("Sorry, but the number of channels in Nifti image : "+img.getNChannels()+" differs from the number of provided Tr : "+Tr.length);
			return null;
		}
		
		for(int c=0;c<img.getNChannels();c++) {
			for(int z=0;z<img.getNSlices();z++) {
				int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0);
				ImageProcessor imgP=img.getStack().getProcessor(sli);
				double[]stats=MRUtils.getBackgroundStatsFromProcessorTight(imgP);
				double sigmaRice=RiceEstimator.computeRiceSigmaFromBackgroundValuesStatic(stats[0],stats[1]);
				double tr=Tr[c];
				double te=TeSpacing;
				String chain="T1SEQ_NIFTI_TR="+VitimageUtils.dou(tr)+"_TE="+VitimageUtils.dou(te)+"_SIGMARICE="+VitimageUtils.dou(sigmaRice);
				img.getStack().setSliceLabel(chain,VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, 0) );
			}		
		}
		
		HyperMap hyper=new HyperMap(img);
		hyper.adjustContrast();
		return hyper;
	}

	public static HyperMap importHyperMapFromNifti4DUnknownSequence(String path,String imageName) {
		ImagePlus img2=IJ.openImage(path);
		System.out.println("Opening image from "+path);
		IJ.run(img2, "Flip Vertically", "stack");
		IJ.log(VitimageUtils.imageResume(img2));
		IJ.log("Readen Nifti");
		ImagePlus[]tab=ValidationExperimentsFijiRelaxPaper.stackToSlicesTframes(img2);
		ImagePlus img=VitimageUtils.hyperStackingChannels(tab);

		IJ.log("GetNiftinInfos");
		Object[] objs=getNiftiInfos(img.getNChannels());
		IJ.log("GotNiftinInfos");
		if(objs==null) {
			IJ.log("Niftin infos were null");
			return null;
		}
		return importHyperMapFromNifti4DSequence(path, imageName, objs);
	}		

	public static Object[] getNiftiInfos(int nChannels) {	
		IJ.log("Entering GetNiftinInfos");
		
		//Parameters for automatic registration
		GenericDialog gd= new GenericDialog("Nifti importer");
		String[]choice=new String[] {"It is a series of T1-weighted image","It is a series of T2-weighted image", "It is a series where both recovery and echo times varies"};
		gd.addChoice("Choose data type",choice,choice[1]);
        gd.showDialog();    if (gd.wasCanceled()) return null;	        
		
        int choi=gd.getNextChoiceIndex();
        
        if(choi==1) {//T2 series
    		IJ.log("GetNiftinInfos --> T2 series");
    		GenericDialog gd2= new GenericDialog("Populate channel informations");
    		gd2.addNumericField("Recovery time ", 1000, 0, 4, "ms");
    		gd2.addNumericField("First echo ", 10, 0, 4, "ms");
    		gd2.addNumericField("Echo spacing ", 10, 0, 4, "ms");
            gd2.showDialog();    if (gd2.wasCanceled()) return null;	        
            double tr=gd2.getNextNumber();
            double te0=gd2.getNextNumber();
            double tedelta=gd2.getNextNumber();

            MRDataType[]mr=new MRDataType[nChannels];
        	double[]trs=new double[nChannels];
        	double[]tes=new double[nChannels];
        	
        	for(int c=0;c<nChannels;c++) {
        		mr[c]=MRDataType.T2SEQ;
        		trs[c]=tr;
        		tes[c]=te0+c*tedelta;
        	}
    		IJ.log("GetNiftinInfos --> End T2 series");
        	return new Object[] {mr,trs,tes,new Boolean[]{true,false,false}};
        }
        if(choi==0) {//T1 series
    		GenericDialog gd2= new GenericDialog("Populate channel informations");
    		gd2.addNumericField("Echo time ", 10, 0, 4, "ms");
    		gd2.addStringField("Recovery times (use ; as a separator)=", "600 ; 1200 ; 2400", 100);
            gd.showDialog();    if (gd.wasCanceled()) return null;	        
            double te=gd2.getNextNumber();
            String  stringTr=gd2.getNextString();
            String[]strInt=stringTr.replace(" ","").split(";");
            double[]trints=new double[strInt.length];
        	for(int d=0;d<trints.length;d++) {
                try{
        		trints[d]=Double.parseDouble(strInt[d]);
                }catch(Exception e) {IJ.showMessage("Wrong Double read : "+strInt[d]+" in string "+stringTr);return null;}
        	}

            MRDataType[]mr=new MRDataType[nChannels];
        	double[]trs=new double[nChannels];
        	double[]tes=new double[nChannels];
        	
        	for(int c=0;c<nChannels;c++) {
        		mr[c]=MRDataType.T1SEQ;
        		trs[c]=trints[c];
        		tes[c]=te;
        	}
        	return new Object[] {mr,trs,tes,new Boolean[]{true,false,false}};
        }
        
        if(choi==2) {//Both
        	IJ.showMessage("Sorry, it seems that this functionality is not yet implemented for this kind of data. Please, contact our support : romain.fernandez@cirad.fr");
        	return null;
        }
    	return null;
 	}

	public ImagePlus getAsImagePlus() {
		return hyperImg.duplicate();
	}
	
	public static void main(String[]args) {
   		@SuppressWarnings("unused")
		ImageJ ij=new ImageJ();
   		runDebugTest();
	}

	public static void runDebugTest() {
		ImagePlus imgTest=IJ.openImage("/home/fernandr/Bureau/FijiRelax_PrepaDOI/Tests_Refactoring/3_Computed_Maps/hyper.tif");
		imgTest.show();
		HyperMap hypTest=new HyperMap(imgTest);
		ImagePlus res=hypTest.simulateOutlierRemoval(150,150,250,250,0,0);
		IJ.showMessage("Outlier removal results with various parameters (Z slicer to change parameters)\nParameters are written in the image label (top line)");
		res.show();
		res.setTitle("Outlier removal results");
		res.setPosition(VitimageUtils.getCorrespondingSliceInHyperImage(res, 0, 1, 0));

		VitimageUtils.setImageWindowSizeTo(res,500);
		if(true)return;
		int oneForTukeyFenceTwoForMADe=2;
		double nStdDev=5;
		int blockHalfSize=4;
		ImagePlus result=hypTest.replaceMapsOutliersSlicePerSlice(oneForTukeyFenceTwoForMADe ,nStdDev,blockHalfSize,true);
		result.show();
		if(true)return;
		
		double[]test=new double[1000];
		for(int i=0;i<1000;i++)test[i]=new Random().nextGaussian()+14;

		System.out.println(TransformUtils.stringVectorN(test,""));
		System.out.println("Stats MADe : "+TransformUtils.stringVectorN(MADeStatsDoubleSided(test,null) , ""));  
		System.out.println("Stats Quart : "+TransformUtils.stringVectorN(getQuartiles(test,null) , ""));  
	    
		for(double val=10;val<18;val+=0.2){
			Object[]obj=MADeIsOutlier(val, test, null, 3);
			Object[]obj2=tuckeyIsOutlier(val, test, null, 3);
			System.out.println("Val="+val+"    "+"MADe : "+obj[0]+" -> "+obj[1]+"      "+"Tuckey : "+obj2[0]+" -> "+obj2[1]+"      ");
	   		}
	   	}

	

	
	
	
	/** Setup object---------------------------------------------------------*/	
	public void setup (){
		setupDimensions();
		setupDataType();
		setupImages();
		setupName();
		setupTrTe();
		setupSigmaTabs();
	}
	
	public void setupDimensions() {
		//Dimensions
		this.T=hyperImg.getNFrames();
		this.actualDay=new String[this.T];
		for(int t=0;t<this.T;t++)this.actualDay[t]=""+t;
		this.dims=VitimageUtils.getDimensions(hyperImg);
		this.X=dims[0];
		this.Y=dims[1];
		this.Z=dims[2];
		this.voxs=VitimageUtils.getVoxelSizes(hyperImg);
		this.Ctot=hyperImg.getNChannels();
	}

	public void setupDataType() {
		hasT1sequence=false;
		hasT2sequence=false;
		hasT1T2sequence=false;		
		hasMaps=false;
		mrDataType=new MRDataType[T][Ctot];
		for(int t=0;t<this.T;t++) {
			for(int c=0;c<this.Ctot;c++) {
				mrDataType[t][c]=MRUtils.getDataTypeOfThisMagneticResonanceSlice(hyperImg, c, 0, t);
				//System.out.println("t="+t+" channel="+c+" : "+mrDataType[t][c]);
				if(mrDataType[t][c]==MRDataType.T1SEQ)hasT1sequence=true;
				if(mrDataType[t][c]==MRDataType.T2SEQ)hasT2sequence=true;
				if(mrDataType[t][c]==MRDataType.PDMAP)hasMaps=true;					
				if(mrDataType[t][c]==MRDataType.T1T2SEQ) {//old version of software
					hasT1sequence=true;
					hasT2sequence=true;
					int sli=VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, 0, t);
					String lab=hyperImg.getStack().getSliceLabel(sli);
					if(lab.contains("TR=10000")) mrDataType[t][c]=MRDataType.T2SEQ;
					else mrDataType[t][c]=MRDataType.T1SEQ;
					for(int z=0;z<Z;z++)hyperImg.getStack().setSliceLabel(hyperImg.getStack().getSliceLabel(sli).replace("T1T2SEQ", ""+mrDataType[t][c]),sli);
				}
			}
		}
		
		
		hasT1T2sequence=(hasT1sequence && hasT2sequence);
		nMaps=hasT1T2sequence ? 4 : 3 ;
		if(hasMaps) this.C=this.Ctot-nMaps;
		else {
			C=Ctot;
			Ctot=C+nMaps;
			MRDataType[][]newMrType=new MRDataType[T][Ctot];
			for(int t=0;t<T;t++) {
				for(int c=0;c<C;c++)newMrType[t][c+nMaps]=mrDataType[t][c];
				newMrType[t][0]=MRDataType.PDMAP;
				newMrType[t][nMaps-1]=MRDataType.MASKMAP;
				if(hasT1T2sequence) {
					newMrType[t][1]=MRDataType.T1MAP;
					newMrType[t][2]=MRDataType.T2MAP;
				}
				else if(hasT1sequence) {
					newMrType[t][1]=MRDataType.T1MAP;
				}
				else {
					newMrType[t][1]=MRDataType.T2MAP;		
				}
			}
			mrDataType=newMrType;
		}
	}
		
	public void setupImages() {
		if(hasMaps) {
			hyperMaps=new Duplicator().run(hyperImg,1,nMaps,1,Z,1,T);
			hyperEchoes=new Duplicator().run(hyperImg,nMaps+1,Ctot,1,Z,1,T);
		}
		else {
			hyperMaps=null;
			hyperEchoes=new Duplicator().run(hyperImg,1,C,1,Z,1,T);
		}
	}
	
	public void setupName() {
		String label=this.hyperEchoes.getStack().getSliceLabel(1);
		if(label==null)return;
		String[]labs=label.split("SEQ_");
		if((labs==null) || labs.length<2)return;
		labs=labs[1].split("_TR");
		if((labs==null) || labs.length<1)return;
		this.name=labs[0];
		hyperImg.setTitle("Hypermap_"+name);
	}

	public void setupTrTe(){
		this.Te=new double[this.T][this.Z][this.C];
		this.Tr=new double[this.T][this.Z][this.C];
		this.sigmaRice=new double[this.T][this.Z][this.C];
		for(int t=0;t<this.T;t++) {
			for(int c=0;c<this.C;c++) {
				for(int z=0;z<dims[2];z++) {
					if(mrDataType[t][c+nMaps]!=MRDataType.T1SEQ  && mrDataType[t][c+nMaps]!=MRDataType.T2SEQ  && mrDataType[t][c+nMaps]!=MRDataType.T1T2SEQ) {
						sigmaRice[t][z][c]=Float.NaN;
						Te[t][z][c]=Float.NaN;
						Tr[t][z][c]=Float.NaN;
					}
					else {
						sigmaRice[t][z][c]=MRUtils.readSigmaInSliceLabel(hyperEchoes, c, z,t);
						Te[t][z][c]=MRUtils.readTeInSliceLabel(hyperEchoes, c, z, t);
						Tr[t][z][c]=MRUtils.readTrInSliceLabel(hyperEchoes, c, z, t);	
					}
				}
			}
		}
	}
	
	public void setupSigmaTabs() {
		if(hasT1T2sequence) {
			this.tabSigmasT1T2Seq=new double[this.T][dims[2]];			
			for(int t=0;t<this.T;t++) {
				for(int z=0;z<dims[2];z++) {
					double average=0;
					for(int ind=0;ind<C;ind++)average+=MRUtils.readSigmaInSliceLabel(hyperEchoes, ind, z, t);
					this.tabSigmasT1T2Seq[t][z]=average/(C-nMaps);
				}
			}
		}
	
		if(hasT1sequence) {
			this.tabSigmasT1Seq=new double[this.T][dims[2]];
			for(int t=0;t<this.T;t++) {
				for(int z=0;z<dims[2];z++) {
					double average=0;
					int[]indices=getT1Indices(t);
					for(int ind:indices)average+=MRUtils.readSigmaInSliceLabel(hyperEchoes, ind, z, t);
					this.tabSigmasT1Seq[t][z]=average/indices.length;
				}
			}
		}
		if(hasT2sequence) {
			this.tabSigmasT2Seq=new double[this.T][dims[2]];
			for(int t=0;t<this.T;t++) {
				for(int z=0;z<dims[2];z++) {
					double average=0;
					int[]indices=getT2Indices(t);
					for(int ind:indices)average+=MRUtils.readSigmaInSliceLabel(hyperEchoes, ind, z, t);
					this.tabSigmasT2Seq[t][z]=average/indices.length;
				}
			}
		}
	}
	
	
	
	
	
    
    
	/** Helpers for accessing to data subparts ---------------------------------------------------------*/	
	public ImagePlus getMask() {
		if(!hasMaps)return null;
		else return new Duplicator().run(hyperMaps,this.nMaps,this.nMaps,1,this.Z,1,this.T);
	}
	
	public ImagePlus getEchoesImage() {
		return hyperEchoes.duplicate();
	}
	
	public ImagePlus getMapsImage() {
		return hyperMaps.duplicate();
	}
	
	public String[][][] getEchoesImageText() {
		String[][][]texts=new String[C][dims[2]][T];
		for(int z=0;z<dims[2];z++)for(int t=0;t<T;t++) for(int c=0;c<C;c++) {
		texts[c][z][t]="   "+MRUtils.getDataTypeOfThisMagneticResonanceSlice(hyperEchoes, c, z, t)+"  Day="+actualDay[t]+"   Z="+z+
					"  TR="+(MRUtils.readTrInSliceLabel(hyperEchoes, c, z, t))+"  TE="+(MRUtils.readTeInSliceLabel(hyperEchoes, c, z, t));
		}
		return texts;
	}

	public ImagePlus getT1EchoesImage(int t) {
		int nbT1=getT1SeqNumberReps(t);
		VitimageUtils.printImageResume(hyperEchoes);
		ImagePlus echoes= new Duplicator().run(hyperEchoes,1,nbT1,1,Z,1,T);
		return echoes;
	}

	public ImagePlus getT2EchoesImage(int t) {
		int nbT2=getT2SeqNumberReps(t);		
		if(!hasT1sequence) {
			return new Duplicator().run(hyperEchoes,1,nbT2,1,Z,1,T);
		}
		return new Duplicator().run(hyperEchoes,getT1SeqNumberReps(t)+1,getT1SeqNumberReps(t)+getT2SeqNumberReps(t),1,Z,1,T);
	}

	public ImagePlus getHyperEcho(int c,int t) {
		if (c>=this.C)return getHyperEcho(this.C-1,t);
		if(c<0)return getHyperEcho(0,t);
		if (t>=this.T)return getHyperEcho(c,this.T-1);
		if(t<0)return getHyperEcho(c,0);
		return new Duplicator().run(hyperEchoes,c+1,c+1,1,Z,t,t);				
	}


	public String toString() {
		String s="";
		s+="HyperMap "+this.name+"\n X Y Z = "+this.X+" x "+this.Y+" x "+this.Z+"\n";
		s+=" C Ctot T = "+this.C+" , "+this.Ctot+" , "+this.T;
		s+="\nVoxel sizes"+TransformUtils.stringVector(this.voxs,"");
		s+="\n\nHas T1 data ? "+this.hasT1sequence;
		s+="\nHas T2 data ? "+this.hasT2sequence;
		s+="\nHas T1T2 data ? "+this.hasT1T2sequence;
		s+="\nHas built maps ? "+this.hasMaps+"\n\n";
		
		for(int t=0; t<T ; t++) {
			s+="\n";
			s+="\n";
			for(int c=0; c<mrDataType[t].length;c++) {
				s+="\n";
				s+="MRdataType["+t+"]["+c+"]="+mrDataType[t][c];
			}
			s+="\n";
			for(int c=0; c<Te[t][0].length;c++) {
				s+="\n";
				s+="TrTe["+t+"]["+c+"]="+Tr[t][0][c]+" , "+Te[t][0][c];
			}
		}
		return s;
	}
	


	
	
	
	/** Helpers for accessing to data T1 and T2 specific informations--------------------------------------------------------- */	
	public int getT1T2SeqNumberReps(int t){
		int nb=0;
		for(int c=0;c<this.Ctot;c++)if( (mrDataType[t][c]==MRDataType.T1SEQ) ||  (mrDataType[t][c]==MRDataType.T2SEQ) )nb++;
		return nb;
	}
	
	public int getT1SeqNumberReps(int t){
		int nb=0;
		for(int c=0;c<this.Ctot;c++)if(mrDataType[t][c]==MRDataType.T1SEQ)nb++;
		return nb;
	}
		
	public int getT2SeqNumberReps(int t){
		int nb=0;
		for(int c=0;c<this.Ctot;c++)if(mrDataType[t][c]==MRDataType.T2SEQ)nb++;
		return nb;
	}
			
	public int[]getT1Indices(int t){
		int[]ret=new int[getT1SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.C;c++)if(mrDataType[t][c+nMaps]==MRDataType.T1SEQ)ret[incr++]=c;		
		return ret;
	}
	
	public int[]getT2Indices(int t){
		int[]ret=new int[getT2SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.C;c++)if(mrDataType[t][c+nMaps]==MRDataType.T2SEQ)ret[incr++]=c;		
		return ret;
	}

	public int[]getT1T2Indices(int t){
		int[]ret=new int[getT1T2SeqNumberReps(t)];
		int incr=0;
		for(int c=0;c<this.C;c++)if( (mrDataType[t][c+nMaps]==MRDataType.T1SEQ) || (mrDataType[t][c+nMaps]==MRDataType.T2SEQ) )ret[incr++]=c;		
		return ret;
	}
	
	public int[][]getT1Indices(){
		int[][]ret=new int[this.T][];
		for(int t=0;t<this.T;t++)ret[t]=this.getT1Indices(t);
		return ret;
	}

	public int[][]getT2Indices(){
		int[][]ret=new int[this.T][];
		for(int t=0;t<this.T;t++)ret[t]=this.getT2Indices(t);
		return ret;
	}

	public int[][]getT1T2Indices(){
		int[][]ret=new int[this.T][];
		for(int t=0;t<this.T;t++)ret[t]=this.getT1T2Indices(t);
		return ret;
	}

	

	
	
	
	/** Helpers for accessing to echo and recovery times--------------------------------------------------------- */		
	public double[][][]getT1T2TeTimes(){
		double[][][][]init=getT1T2TrTeTimes();
		double[][][]ret=new double[T][Z][];
		for(int t=0;t<T;t++)for(int z=0;z<Z;z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
				ret[t][z][c]=init[t][z][c][1];
			}
		}
		return ret;
	}

	public double[][][]getT1TeTimes(){
		double[][][][]init=getT1TrTeTimes();
		double[][][]ret=new double[this.T][dims[2]][];
		for(int t=0;t<this.T;t++)for(int z=0;z<dims[2];z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
				ret[t][z][c]=init[t][z][c][1];
			}
		}
		return ret;
	}
	
	public double[][][]getT2TeTimes(){
		double[][][][]init=getT2TrTeTimes();
		double[][][]ret=new double[this.T][dims[2]][];
		for(int t=0;t<this.T;t++)for(int z=0;z<dims[2];z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
				ret[t][z][c]=init[t][z][c][1];
			}
		}
		return ret;
	}

	
	public double[][][]getT1T2TrTimes(){
		double[][][][]init=getT1T2TrTeTimes();
		double[][][]ret=new double[this.T][dims[2]][];
		for(int t=0;t<this.T;t++)for(int z=0;z<dims[2];z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
//				System.out.println("DEB "+init[t][z][c][0]);
				ret[t][z][c]=init[t][z][c][0];
				
			}
		}
		return ret;
	}

	public double[][][]getT1TrTimes(){
		double[][][][]init=getT1TrTeTimes();
		double[][][]ret=new double[this.T][dims[2]][];
		for(int t=0;t<this.T;t++)for(int z=0;z<dims[2];z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
				ret[t][z][c]=init[t][z][c][0];
			}
		}
		return ret;
	}

	public double[][][]getT2TrTimes(){
		double[][][][]init=getT2TrTeTimes();
		double[][][]ret=new double[this.T][dims[2]][];
		for(int t=0;t<this.T;t++)for(int z=0;z<dims[2];z++) {
			ret[t][z]=new double[init[t][z].length];
			for(int c=0;c<init[t][z].length;c++) {
				ret[t][z][c]=init[t][z][c][0];
			}
		}
		return ret;
	}

	
	public double[][][]getT1T2TrTeTimes(int t){
		int incr=0;
		double[][][]ret=new double[dims[2]][getT1T2SeqNumberReps(t)][2];
		for(int c : getT1T2Indices(t)) {
	//	for(int c=0;c<this.Ctot;c++) {
			//if((mrDataType[t][c]==MRDataType.T1SEQ) || (mrDataType[t][c]==MRDataType.T2SEQ)) {
//			for(int z=0;z<dims[2];z++)ret[z][incr]=new double[] {this.Tr[t][z][c-((hasMaps||true) ? nMaps : 0)],this.Te[t][z][c-((hasMaps||true) ? nMaps : 0)]};
				for(int z=0;z<dims[2];z++)ret[z][incr]=new double[] {this.Tr[t][z][c],this.Te[t][z][c]};
				incr++;
			//}
		}
		return ret;
	}
	
	public double[][][]getT1TrTeTimes(int t){
		int incr=0;
		double[][][]ret=new double[dims[2]][getT1SeqNumberReps(t)][2];
		for(int c : getT1Indices(t)) {
//			if(mrDataType[t][c]==MRDataType.T1SEQ) {
				for(int z=0;z<dims[2];z++) {
					ret[z][incr]=new double[] {this.Tr[t][z][c],this.Te[t][z][c]};
					//ret[z][incr]=new double[] {this.Tr[t][z][c-(hasMaps ? nMaps : 0)],this.Te[t][z][c-(hasMaps ? nMaps : 0)]};
				}
				incr++;
	//		}
		}
		return ret;
	}
	
	public double[][][]getT2TrTeTimes(int t){
		int incr=0;
		double[][][]ret=new double[dims[2]][getT2SeqNumberReps(t)][2];
//		for(int c=0;c<this.Ctot;c++) {
		for(int c : getT2Indices(t)) {
//			if(mrDataType[t][c]==MRDataType.T2SEQ) {
						for(int z=0;z<dims[2];z++) ret[z][incr]=new double[] {this.Tr[t][z][c],this.Te[t][z][c]};
	//			for(int z=0;z<dims[2];z++) ret[z][incr]=new double[] {this.Tr[t][z][c-(hasMaps ? nMaps : 0)],this.Te[t][z][c-(hasMaps ? nMaps : 0)]};
				incr++;
	//		}
		}
		return ret;
	}

	
	public double[][][][]getT1T2TrTeTimes(){
		double[][][][]ret=new double[T][][][];
		for(int t=0;t<this.T;t++)ret[t]=getT1T2TrTeTimes(t);
		return ret;
	}

	public double[][][][]getT2TrTeTimes(){
		double[][][][]ret=new double[T][][][];
		for(int t=0;t<this.T;t++)ret[t]=getT2TrTeTimes(t);
		return ret;
	}

	public double[][][][]getT1TrTeTimes(){
		double[][][][]ret=new double[T][][][];
		for(int t=0;t<this.T;t++)ret[t]=getT1TrTeTimes(t);
		return ret;
	}

	
	
	
	
	
	
	
	/** Helpers for dynamic maps computation, give access to magnitude data and coordinates in a neighbourhood around a voxel---------------------------------------------------------*/
	public double[][][][]getFullMRISignalAroundThisVoxelT1T2(int x,int y, int z,int crossWidth,int crossThick,double sigmaXYInVoxels){
		int xm,ym,xM,yM,zm,zM;
		xm=Math.max(0, x-crossWidth);
		xM=Math.min(X-1, x+crossWidth);
		ym=Math.max(0, y-crossWidth);
		yM=Math.min(Y-1, y+crossWidth);
		zm=Math.max(0, z-crossThick);
		zM=Math.min(Z-1, z+crossThick);

		int xm2,ym2,xM2,yM2;
		xm2=Math.max(0, x-crossWidth-5);
		xM2=Math.min(X-1, x+crossWidth+5);
		ym2=Math.max(0, y-crossWidth-5);
		yM2=Math.min(Y-1, y+crossWidth+5);

		
		int[][]t1t2Indices=this.getT1T2Indices();
		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		double[][][][]data= new double[T][nHits][1][];//[times][vox][Seq][echos]
		
		int currentChan;
		int index;
		
		ImagePlus temp=VitimageUtils.cropMultiChannelFloatImage(hyperEchoes,xm2,xM2,ym2,yM2,zm,zM);
		temp=VitimageUtils.gaussianFilteringMultiChannel(temp,sigmaXYInVoxels,sigmaXYInVoxels,0);//It's no error : no "big smoothing" over Z, due to misalignment
	
		for(int t=0;t<this.T;t++) {			
			for(int zz=zm;zz<=zM;zz++) {
				for(int xx=xm;xx<=xM;xx++) {
					for(int yy=ym;yy<=yM;yy++) {
						index=(xx-xm) + (xM-xm+1)*(yy-ym) + (xM-xm+1)*(yM-ym+1)*(zz-zm);
						data[t][index][0]=new double[t1t2Indices[t].length];
						for(int t1trInd=0;t1trInd<t1t2Indices[t].length;t1trInd++) {			
							currentChan=t1t2Indices[t][t1trInd];
							int indexSlice=this.C*(zM-zm+1)*t + this.C*(zz-zm) + currentChan + 1;
							data[t][index][0][t1trInd]=temp.getStack().getProcessor(indexSlice).getPixelValue(xx-xm2, yy-ym2);
						}
					}
				}	
			}
		}
		return data;
	}

	public double[][][][]getFullMRISignalInTheseCoordinatesT1T2(int xCor,int yCor, int zCor,int[][]coordinates,double sigmaXYInVoxels){
		int[][]t1t2Indices=this.getT1T2Indices();
		int nHits=coordinates.length;
		int bboxX0=100000000;
		int bboxXf=0;
		int bboxY0=100000000;
		int bboxYf=0;
		for(int vo=0;vo<nHits;vo++) {
			if(coordinates[vo][0]<bboxX0)bboxX0=coordinates[vo][0];
			if(coordinates[vo][0]>bboxXf)bboxXf=coordinates[vo][0];
			if(coordinates[vo][1]<bboxY0)bboxY0=coordinates[vo][1];
			if(coordinates[vo][1]>bboxYf)bboxYf=coordinates[vo][1];
		}
		bboxX0=Math.max(0, bboxX0-5);
		bboxXf=Math.min(X-1, bboxXf+5);
		bboxY0=Math.max(0, bboxY0-5);
		bboxYf=Math.min(Y-1, bboxYf+5);
	
		ImagePlus temp=VitimageUtils.cropMultiChannelFloatImage(hyperEchoes,bboxX0,bboxXf,bboxY0,bboxYf,zCor,zCor);
		temp=VitimageUtils.gaussianFilteringMultiChannel(temp,sigmaXYInVoxels,sigmaXYInVoxels,0);//It's no error : no "big smootphrahing" over Z, due to misalignment
		double[][][][]data= new double[t1t2Indices.length][nHits][1][];//[times][vox][Seq][echos]		
		int currentChan;
		int xx,yy;	
		for(int t=0;t<this.T;t++) {			
			for(int vo=0;vo<nHits;vo++) {
				xx=coordinates[vo][0]-bboxX0;
				yy=coordinates[vo][1]-bboxY0;
				data[t][vo][0]=new double[t1t2Indices[t].length];
				for(int t1trInd=0;t1trInd<t1t2Indices[t].length;t1trInd++) {			
					currentChan=t1t2Indices[t][t1trInd];
					int indexSlice=this.C*t + currentChan + 1;
					data[t][vo][0][t1trInd]=temp.getStack().getProcessor(indexSlice).getPixelValue(xx, yy);
				}
			}
		}
		return data;
	}

	public int[][]getCoordinatesAroundThisVoxel(int x,int y, int z,int crossWidth,int crossThick){
		int xm,ym,xM,yM,zm,zM;
		xm=Math.max(0, x-crossWidth);
		xM=Math.min(dims[0]-1, x+crossWidth);
		ym=Math.max(0, y-crossWidth);
		yM=Math.min(dims[1]-1, y+crossWidth);
		zm=Math.max(0, z-crossThick);
		zM=Math.min(dims[2]-1, z+crossThick);

		int nHits=(xM-xm+1)*(yM-ym+1)*(zM-zm+1);
		int[][]coords= new int[nHits][3];//[vox][dim]
		for(int zz=zm;zz<=zM;zz++) {
			for(int xx=xm;xx<=xM;xx++) {
				for(int yy=ym;yy<=yM;yy++) {
					int index=(xx-xm) + (xM-xm+1)*(yy-ym) + (xM-xm+1)*(yM-ym+1)*(zz-zm);
					coords[index]=new int[] {xx,yy,zz};
				}
			}
		}
		return coords;
	}	
	
	
/*
	public static double getMaxValueForContrastIntelligent(ImagePlus img,int channel,int time,int z,double percentageKeep,double factor) {
		int nBins=256;
		int cOld=img.getC();
		int tOld=img.getT();
		int zOld=img.getZ();
		img.setC(channel);
		img.setT(time);
		img.setZ(z);
		ImageStatistics stats = img.getStatistics(ImagePlus.AREA+ImagePlus.MEAN+ImagePlus.MODE+ImagePlus.MIN_MAX, nBins);
		double[] y = stats.histogram();
		double wid=stats.binSize;
		double sum=0;for(double yy : y)sum+=yy;
		double cum=0;int indexFin=0;
		while(cum<sum*(percentageKeep/100.0)){cum+=y[indexFin++];if(false && cum>sum*0.95)System.out.println("Trying threshold index "+indexFin+" = "+(indexFin*wid)+" : cum="+cum+"/"+sum);}
		return (factor*(wid*indexFin));
	}
	*/
	public static double getMaxValueForContrastMoreIntelligent(ImagePlus img2,int channel,int time,int z,double percentageKeep,double factor) {
		Timer t=new Timer();
		int nBins=256;
		int X=img2.getWidth();
		int Y=img2.getHeight();
		ImagePlus img=new Duplicator().run(img2,channel+1,channel+1,z+1,z+1,time+1,time+1);
		double []d=VitimageUtils.valuesOfBlock(img, 0, 0, 0, X-1, Y-1, 0);
		Arrays.sort(d);
		int index=(int)Math.round(d.length*percentageKeep/100.0);
		//System.out.println("Setting with index="+index+"/"+d.length+" . Val0="+d[0]+" valindex="+d[index]+" valfinale="+d[d.length-1]);
		t.print("After sorting");
		return (d[index]*factor);
	}



	public void adjustContrast() {
		ImagePlus echoes2=this.getEchoesImage();
		ImagePlus echoes=new Duplicator().run(echoes2,1,echoes2.getNChannels(),1,echoes2.getNSlices(),1,1);
		if(hasMaps) {
			ImagePlus maps2=this.getMapsImage();
			ImagePlus maps=new Duplicator().run(maps2,1,maps2.getNChannels(),1,maps2.getNSlices(),1,1);
			double[][]rangesMaps=new double[maps.getNChannels()][2];
			for(int c=0;c<maps.getNChannels();c++) {
				ImagePlus mapTemp=new Duplicator().run(maps,c+1,c+1,maps.getNSlices()/2+1,maps.getNSlices()/2+1,1,1);
				rangesMaps[c]=new double[] {0,getMaxValueForContrastMoreIntelligent(mapTemp,1,1,1,percentageKeepNormalisation,factorViewNormalisation)};
			}
			maps=this.hyperMaps;
			for(int t=0;t<T;t++) {
				maps.setT(t+1);
				this.hyperImg.setT(t+1);
				for(int c=0;c<nMaps;c++) {
					maps.setC(c+1);
					maps.setDisplayRange(rangesMaps[c][0], rangesMaps[c][1]);
					this.hyperImg.setC(c+1);
					this.hyperImg.setDisplayRange(rangesMaps[c][0], rangesMaps[c][1]);
					//System.out.println("Setting at "+t+" , "+c+" : |"+rangesMaps[c][0]+","+rangesMaps[c][1]+"|");
				}
			}
			maps.updateAndDraw();
		}
		
		
		ImagePlus maxEcho=VitimageUtils.maxOfImageArray(VitimageUtils.stacksFromHyperstackFastBis(echoes));
		maxEcho.setZ(maxEcho.getNSlices()/2+1);
		double[]rangesEchoes=new double[] {0,getMaxValueForContrastMoreIntelligent(maxEcho,1,1,1,percentageKeepNormalisation,factorViewNormalisation)};

		
		echoes=this.hyperEchoes;
		for(int t=0;t<T;t++) {
			echoes.setT(t+1);
			this.hyperImg.setT(t+1);
			for(int c=0;c<echoes.getNChannels();c++) {
				echoes.setC(c+1);
				echoes.setDisplayRange(rangesEchoes[0], rangesEchoes[1]);
				this.hyperImg.setC(c+1+(hasMaps ? nMaps : 0));
				this.hyperImg.setDisplayRange(rangesEchoes[0], rangesEchoes[1]);
				IJ.run(echoes,"Fire","");
				IJ.run(hyperImg,"Fire","");
			}
		}
		echoes.setT(1);
		echoes.setC(1);
		echoes.setZ(1);
		hyperImg.setT(1);
		hyperImg.setC(1);
		hyperImg.setZ(1);
		echoes.updateAndDraw();
		hyperImg.updateAndDraw();
	}


	
	
	/** Maps computation routine---------------------------------------------------------*/	
	public ImagePlus computeMaps() {
		return computeMapsAgainAndMask(MRUtils.SIMPLEX,false,NoiseManagement.RICE,false,null,4);
	}

	public ImagePlus computeMaps(NoiseManagement noise) {
		return computeMapsAgainAndMask(MRUtils.SIMPLEX,false,noise,false,null,4);
	}

	public ImagePlus computeMapsNoJoint() {
		return computeMapsAgainAndMask(MRUtils.SIMPLEX,true,NoiseManagement.RICE,false,null,4);
	}
	
	public ImagePlus computeMapsNoJoint(NoiseManagement noi) {
		return computeMapsAgainAndMask(MRUtils.SIMPLEX,true,noi,false,null,4);
	}


	public ImagePlus computeMapsWithMask(ImagePlus mask) {
		return computeMapsAgainAndMask(MRUtils.SIMPLEX,false,NoiseManagement.RICE,false,mask,4);
	}

	
	public ImagePlus computeMaps(double nStdDevForMask) {
		return computeMapsAgainAndMask(MRUtils.SIMPLEX,false,NoiseManagement.RICE,false,null,nStdDevForMask);
	}

	public ImagePlus computeMaps(double[]paramsGUI,ImagePlus imgMask) {
		int algType=(paramsGUI[0]==0) ? MRUtils.SIMPLEX : MRUtils.LM;
		NoiseManagement noi=null;
		boolean separated=(paramsGUI[5]!=0);
		boolean forgetFistEcho=(paramsGUI[2]!=0);
		if(paramsGUI[1]==0)noi=NoiseManagement.RICE;
		else noi=(paramsGUI[1]==1 ? NoiseManagement.OFFSET : NoiseManagement.NOTHING);
		double nbStd=paramsGUI[4];
		return computeMapsAgainAndMask(algType,separated ,noi,forgetFistEcho,imgMask,nbStd);
	}

	public ImagePlus computeMapsAgainAndMask(int algType,boolean separated,NoiseManagement noise,boolean forgetFirstEcho,ImagePlus imgMaskUser,double nbStdNoiseAroundMeanForSelection) {
		if(T>1) {
			
			ImagePlus []imgTab=new ImagePlus[T];
			for(int t=0;t<T;t++) {
				ImagePlus imgTemp=new Duplicator().run(this.hyperImg,1,hasMaps ? Ctot : C,1,Z,t+1,t+1);
				HyperMap hypTemp=new HyperMap(imgTemp);
				imgTab[t]=hypTemp.computeMapsAgainAndMask(algType, separated, noise, forgetFirstEcho, imgMaskUser, nbStdNoiseAroundMeanForSelection);
			}
			double[][]ranges=new double[Ctot][2];
			for(int c=0;c<Ctot;c++) {
				imgTab[0].setC(c+1);
				ranges[c]=new double[] {imgTab[0].getDisplayRangeMin(),imgTab[0].getDisplayRangeMax()};
			}
			ImagePlus img=VitimageUtils.hyperStackingFrames(imgTab);
			for(int c=0;c<Ctot;c++) {
				img.setC(c+1);
				img.setDisplayRange(ranges[c][0],ranges[c][1]);
			}
			this.hasMaps=true;
			this.updateBothFromNewHyperImg(img);
//			System.out.println("Just before adjust at the end\n"+this.getRangeInfos()+"\n\n");
			adjustContrast();
//			System.out.println("Just after adjust at the end\n"+this.getRangeInfos()+"\n\n");
			return this.hyperImg;
		}

		int bitD=hyperImg.getBitDepth();
		int nVox=Z*T*X*Y;
		if(bitD!=16 && bitD!=32) {IJ.showMessage("Unexepected data type : no 16 nor 32-bits");System.exit(0);}
		ImagePlus imgMask=null;
		boolean maskGiven=false;		

		ImagePlus []maps=new ImagePlus[hasT1T2sequence ? 3 : 2];
		ImagePlus []imgT1T2Line=VitimageUtils.stacksFromHyperstackFastBis(hyperEchoes);
		if(imgT1T2Line[0].getType()==ImagePlus.GRAY16) {
			for(int c=0;c<imgT1T2Line.length;c++) {
				IJ.run(imgT1T2Line[c],"32-bit","");
			}
		}

		//Mask computation
		double meanRice=meanRice(imgT1T2Line);
		if(imgMaskUser!=null) {maskGiven=true;imgMask=imgMaskUser;IJ.run(imgMask,"32-bit","");}
		else {
			int index=0;
			if(hasT2sequence) {
				index=getT2Indices()[0][forgetFirstEcho ? 1 : 0];
			}
			else index=getT1Indices()[0][getT1Indices()[0].length-1];
			ImagePlus max=imgT1T2Line[index];
			imgMask=VitimageUtils.getFloatBinaryMask(max,(1.2+nbStdNoiseAroundMeanForSelection)*meanRice,1E10);
		}

		//Extract sequences
		ImagePlus[]imgT1s=hasT1sequence ? VitimageUtils.stacksFromHyperstackFastBis(getT1EchoesImage(0)) : null;
		ImagePlus[]imgT2s=hasT2sequence ? VitimageUtils.stacksFromHyperstackFastBis(getT2EchoesImage(0)) : null;
		if(imgT1s!=null && imgT1s[0].getType()==ImagePlus.GRAY16) {
			for(int c=0;c<imgT1s.length;c++) {
				IJ.run(imgT1s[c],"32-bit","");
			}
		}
		if(imgT2s!=null && imgT2s[0].getType()==ImagePlus.GRAY16) {
			for(int c=0;c<imgT2s.length;c++) {
				IJ.run(imgT2s[c],"32-bit","");
			}
		}

		//Compute the needed fits
		Timer t;
		t=new Timer();
		if(!hasT1T2sequence) {
			//If only T1 compute PD T1 and stack it into hypermap
			if(hasT1sequence) {
				if(noise==NoiseManagement.OFFSET || noise==NoiseManagement.NOTHING) {
					maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T1_MONO,algType,false,false);		
				}
				else maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T1_MONO_RICE,algType,false,false);	
			}
			//if only T2 compute PD T1 and stack it into hypermap
			if(hasT2sequence) {
				if(noise==NoiseManagement.NOTHING) {
					maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T2_MONO,algType,false,false);		
				}
				else if(noise==NoiseManagement.OFFSET) {
					maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T2_MONO_BIAS,algType,false,false);		
				}
				else if(noise==NoiseManagement.RICE) {
					maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T2_MONO_RICE,algType,false,false);		
				}
				//typeMaps=new MRDataType[] {MRDataType.PDMAP,MRDataType.T2MAP};
			}
		}
		else {
			//typeMaps=new MRDataType[] {MRDataType.PDMAP,MRDataType.T1MAP,MRDataType.T2MAP};
			if(separated) {
				ImagePlus [][]mapsTemp=new ImagePlus[2][];
				if(noise==NoiseManagement.NOTHING) {
					mapsTemp[0]=MRUtils.computeT1T2MapMultiThreadSlices(imgT1s, imgMask,sigmaInUse,MRUtils.T1_MONO,algType,false,false);		
					mapsTemp[1]=MRUtils.computeT1T2MapMultiThreadSlices(imgT2s, imgMask,sigmaInUse,MRUtils.T2_MONO,algType,false,false);		
				}
				else if(noise==NoiseManagement.OFFSET) {
					mapsTemp[0]=MRUtils.computeT1T2MapMultiThreadSlices(imgT1s, imgMask,sigmaInUse,MRUtils.T1_MONO,algType,false,false);		
					mapsTemp[1]=MRUtils.computeT1T2MapMultiThreadSlices(imgT2s, imgMask,sigmaInUse,MRUtils.T2_MONO_BIAS,algType,false,false);		
				}
				else if(noise==NoiseManagement.RICE) {
					mapsTemp[0]=MRUtils.computeT1T2MapMultiThreadSlices(imgT1s, imgMask,sigmaInUse,MRUtils.T1_MONO_RICE,algType,false,false);		
					mapsTemp[1]=MRUtils.computeT1T2MapMultiThreadSlices(imgT2s, imgMask,sigmaInUse,MRUtils.T2_MONO_RICE,algType,false,false);		
				}
				maps=new ImagePlus[] {mapsTemp[1][0],mapsTemp[0][1],mapsTemp[1][1]};

			}
			else {
				if(noise==NoiseManagement.NOTHING) {
					maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T1T2_MONO,algType,false,false);	
				}
				else if(noise==NoiseManagement.OFFSET) {
					maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T1T2_MONO_BIAS,algType,false,false);	
				}
				else if(noise==NoiseManagement.RICE) {
					maps=MRUtils.computeT1T2MapMultiThreadSlices(imgT1T2Line, imgMask,sigmaInUse,MRUtils.T1T2_MONO_RICE,algType,false,false);	
				}
			}
		}
		t.print("End fitting, time per voxel = "+(VitimageUtils.dou(t.getTime()*1000.0/nVox)));
		VitimageUtils.waitFor(20);
		
		//Handle image range for viewing
		double maxPD=VitimageUtils.maxOfImage(VitimageUtils.maxOfImageArray(imgT1T2Line))*MRUtils.maxDisplayedPDratio;
		double maxT1=VitimageUtils.max(getT1TrTimes()[0][0])*MRUtils.maxDisplayedT1ratio;
		if(maxT1<3500)maxT1=MRUtils.minPossibleMaxT1;
		double maxT2=VitimageUtils.max(getT1T2TeTimes()[0][0])*MRUtils.maxDisplayedT2ratio;


		//In case of a computed mask, also apply the threshold routine to the PD results to compose a single excluding mask
		if(!maskGiven) {			
			ImagePlus imgMask2=VitimageUtils.getFloatBinaryMask(maps[0],(1.2+nbStdNoiseAroundMeanForSelection)*meanRice,1E10);			
			imgMask=VitimageUtils.makeOperationBetweenTwoImages(imgMask,imgMask2,2,true);
		}

		
		//Apply mask to result maps and set result names in slice labels
		for(int m=0;m<nMaps-1;m++) {
			maps[m]=VitimageUtils.makeOperationBetweenTwoImages(maps[m],imgMask,2,true);
			for(int z=0;z<Z;z++) {
				maps[m].getStack().setSliceLabel(""+mrDataType[0][m]+"_"+name, z+1);
			}
		}
		for(int z=0;z<Z;z++) {
			imgMask.getStack().setSliceLabel(""+"MASKMAP"+"_"+name, z+1);
			imgMask.setDisplayRange(0, 2);
		}

		ImagePlus []tempRes=new ImagePlus[C+nMaps];
		for(int m=0;m<nMaps-1;m++) {
			tempRes[m]=maps[m];
		}
		tempRes[nMaps-1]=imgMask;
		for(int c=0;c<imgT1T2Line.length;c++) {
			tempRes[c+nMaps]=imgT1T2Line[c];
		}
		
		ImagePlus newHyperImg=Concatenator.run(tempRes);
		VitimageUtils.printImageResume(newHyperImg,"New hyperMap");

		newHyperImg=HyperStackConverter.toHyperStack(newHyperImg, C+nMaps,Z,T,"xyztc","Grayscale");		
		hyperImg=newHyperImg;
		hyperImg.setTitle("Hypermap_"+name);
		hyperMaps=new Duplicator().run(hyperImg,1,nMaps,1,Z,1,T);
		hyperEchoes=new Duplicator().run(hyperImg,nMaps+1,Ctot,1,Z,1,T);

		if(!maskGiven) {
			System.out.println("Mean sigma rice="+meanRice+", object threshold="+(1.2+nbStdNoiseAroundMeanForSelection)*meanRice);
		}
		hasMaps=true;
		adjustContrast();
		return hyperImg;
	}
	
	public void registerEchoes(){
		registerEchoes(this.getDefaultRegistrationSettings());
	}
	
	public void registerEchoes(RegistrationAction regAct){
		double PD=this.hyperEchoes.getDisplayRangeMax();
		if(! hasT1sequence)IJ.log("RegisterEchoes in HyperMap : no T1 sequence, thus nothing to do there");
		ImagePlus []echoes=VitimageUtils.stacksFromHyperstackFastBis(hyperEchoes);
		int nTr=-1;
		int nTe=-1;
		double curTr=-100;
		double[][][][]trteinit=getT1T2TrTeTimes();
		double[][][]ret=new double[this.T][dims[2]][];
		ArrayList<ArrayList<ImagePlus>> images=new ArrayList<ArrayList<ImagePlus>>();

		for(int c=0;c<C;c++) {
			if(curTr!=Tr[0][0][c]) {
				curTr=Tr[0][0][c];
				nTr++;
				nTe=-1;
				images.add(new ArrayList<ImagePlus>());
				images.get(nTr).add(echoes[c]);
			}
			else {
				nTe++;
				images.get(nTr).add(echoes[c]);
			}
		}
		
		nTr++;
		ItkTransform []trs=new ItkTransform[nTr-1];
		
		Timer t=new Timer();
		t.print("Starting Registration");
		for(int i=0;i<nTr-1;i++) {
			BlockMatchingRegistration bmRegistration;
			ItkTransform tr=null;

			//Linear registration
			if(regAct.getIterationsBMLinear()>0) {
				regAct.typeTrans=Transform3DType.RIGID;
				bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(images.get(nTr-1).get(0),images.get(i).get(0) ,regAct);
				tr=bmRegistration.runBlockMatching(null,false);
				bmRegistration.closeLastImages();
			}
			
			if(regAct.getIterationsBMNonLinear()>0) {
				regAct.typeTrans=Transform3DType.DENSE;
				bmRegistration=BlockMatchingRegistration.setupBlockMatchingRegistration(images.get(nTr-1).get(0),images.get(i).get(0) ,regAct);
				tr=bmRegistration.runBlockMatching(tr,false);
				bmRegistration.closeLastImages();
			}
			for(int j=0;j<images.get(i).size();j++) {
				ImagePlus temp=tr.transformImage(images.get(nTr-1).get(0),images.get(i).get(j));
				images.get(i).set(j,temp);
			}
		    trs[i]=tr;
		}
		
		int incr=0;
		for(int i=0;i<images.size();i++)for(int j=0;j<images.get(i).size();j++) {echoes[incr++]=images.get(i).get(j);}
		this.hyperEchoes=VitimageUtils.hyperStackingChannels(echoes);
		if(hasMaps)this.updateHyperImgFromEchoesAndMaps();
		else {
			this.hyperImg=this.hyperEchoes;
			this.setDisplayRange();
		}
		adjustContrast();
	}
		
	public RegistrationAction getDefaultRegistrationSettings() {
		RegistrationAction regAct=new RegistrationAction();
		regAct.defineSettingsSimplyFromTwoImages(this.getHyperEcho(0,0),this.getHyperEcho(0,0));
		regAct.typeAutoDisplay=2;
		regAct.higherAcc=0;
		//In Fijiyama default settings, images could have a very different geometry. Here we have only minor translations / rotations to correct
		regAct.levelMaxLinear=2;
		regAct.levelMinLinear=-1;
		regAct.strideX*=2;
		regAct.strideY*=2;
		regAct.iterationsBMDen=0;
		return regAct;
	}

	
	


	
	/** HyperMap global update when a part of it changed ---------------------------------------------------------*/	
	public void updateHyperImgFromEchoesAndMaps() {
		ImagePlus []tab=new ImagePlus[Ctot*T];
		ImagePlus[]maps=VitimageUtils.stacksFromHyperstackFastBis(hyperMaps);
		ImagePlus[]echoes=VitimageUtils.stacksFromHyperstackFastBis(hyperEchoes);
		for(int n=0;n<maps.length;n++) {
			tab[n]=maps[n];
		}
		for(int n=maps.length;n<maps.length+echoes.length;n++) {
			tab[n]=echoes[n-maps.length];
		}

		Concatenator con=new Concatenator();
		con.setIm5D(true);
		ImagePlus hypTemp=con.concatenate(tab,true);
		String codeStacking="xyztc";
		ImagePlus hyperNew=HyperStackConverter.toHyperStack(hypTemp, Ctot, Z,T,codeStacking,"Grayscale");
		hyperImg=hyperNew;
		for(int c=Ctot-1;c>=0;c--)for(int t=T-1;t>=0;t--) {hyperImg.setC(c+1);hyperImg.setT(t+1);IJ.run(hyperImg,"Fire","");}
		adjustContrast();
	}
	
	/** HyperMap global update when a part of it changed ---------------------------------------------------------*/	
	public void updateHyperImgFromEchoesAndMapsOld() {
		ImagePlus []tab=new ImagePlus[Ctot];
		ImagePlus[]maps=VitimageUtils.stacksFromHyperstackFastBis(hyperMaps);
		ImagePlus[]echoes=VitimageUtils.stacksFromHyperstackFastBis(hyperEchoes);
		for(int c=0;c<nMaps;c++) {
			tab[c]=maps[c];
		}
		for(int c=nMaps;c<Ctot;c++) {
			tab[c]=echoes[c-nMaps];
		}
		hyperImg=VitimageUtils.hyperStackingChannels(tab);
		this.setDisplayRange();
	}
	
	public void updateEchoesAndMaps() {
		if(hasMaps) {
			hyperMaps=new Duplicator().run(hyperImg,1,nMaps,1,Z,1,T);
			hyperEchoes=new Duplicator().run(hyperImg,nMaps+1,Ctot,1,Z,1,T);
			VitimageUtils.printImageResume(hyperMaps,"hypMa after update");
			VitimageUtils.printImageResume(hyperEchoes,"hypEc after update");
		}
		else hyperEchoes=hyperImg.duplicate();
		adjustContrast();
	}
	
	public void updateBothFromNewHyperImg(ImagePlus hyperImgNew) {
		this.hyperImg=hyperImgNew.duplicate();
		updateEchoesAndMaps();
		this.setDisplayRange();
		adjustContrast();
	}
	
	public void setDisplayRange() {
		double maxPD=VitimageUtils.maxOfImage(hyperEchoes)*MRUtils.maxDisplayedPDratio;
		double maxT1=VitimageUtils.max(getT1T2TrTimes()[0][0])*MRUtils.maxDisplayedT1ratio;
		double maxT2=VitimageUtils.max(getT1T2TeTimes()[0][0])*MRUtils.maxDisplayedT2ratio;

		if(hasMaps) {
			for(int c=0;c<this.nMaps;c++) {
				if((mrDataType[0][c]==MRDataType.PDMAP)) {hyperImg.setC(c+1);hyperImg.setDisplayRange(0, maxPD);}
				else if(c<nMaps && (mrDataType[0][c]==MRDataType.T1MAP)) {hyperImg.setC(c+1);hyperImg.setDisplayRange(0, maxT1);}
				else if(c<nMaps && (mrDataType[0][c]==MRDataType.T2MAP)) {hyperImg.setC(c+1);hyperImg.setDisplayRange(0, maxT2);}
				else {hyperImg.setC(c+1);hyperImg.setDisplayRange(0, 2);}//mask
				IJ.run(hyperImg,"Fire","");
			}
			for(int c=nMaps;c<Ctot;c++) {
				hyperImg.setC(c+1);hyperImg.setDisplayRange(0, maxPD);
				IJ.run(hyperImg,"Fire","");
			}
		}
		else {
			for(int c=0;c<C;c++) {
				hyperImg.setC(c+1);hyperImg.setDisplayRange(0, maxPD);
				IJ.run(hyperImg,"Fire","");
			}		
		}
		hyperImg.setC(1);
	}
	
	public HyperMap pruneBionanoStuff() {
		System.out.println("1");
		this.hyperEchoes.show();
		ImagePlus []tab=VitimageUtils.stacksFromHyperstackFastBis(this.hyperEchoes);
		int incr=0;
		int ZZ=tab[0].getNSlices();
		double[]trTab=MRUtils.getTrFrom3DRelaxationImageTab(tab);
		double[]teTab=MRUtils.getTeFrom3DRelaxationImageTab(tab)[0];
		System.out.println("2");
		System.out.println(tab.length+"");
		for(int i=0;i<tab.length;i++) {
			System.out.println(i);
			if( ! (trTab[i]<10000 && teTab[i]>13))incr++;
		}
		ImagePlus []ret=new ImagePlus[incr];
		incr=0;
		System.out.println("3");
		for(int i=0;i<tab.length;i++) {
			if( ! (trTab[i]<10000 && teTab[i]>13))ret[incr++]=tab[i];
		}
		ImagePlus newHyperImg=Concatenator.run(ret);
		VitimageUtils.printImageResume(newHyperImg,"Hyper");
		System.out.println("4");

		newHyperImg=HyperStackConverter.toHyperStack(newHyperImg, incr,ZZ,1,"xyztc","Grayscale");	
		return new HyperMap(newHyperImg);
	}

	
	
	
	
	/** Helpers for outlier detection and removal--------------------------------------------------------- */	
    public static Object[] tuckeyIsOutlier(double val,double []tabIn,double[]mask,double eRatio) {
	    double[]tabStats=getQuartiles(tabIn,mask);
	    double interQuart=tabStats[2]-tabStats[0];
    	if(val<tabStats[1]-eRatio/2.0*interQuart)return new Object[] {true,tabStats[1],tabStats[1]};
    	if(val>tabStats[1]+eRatio/2.0*interQuart)return new Object[] {true,tabStats[1],tabStats[1]};
    	return new Object[]{false,val,tabStats[1]};
    }
    
    public static Object[] MADeIsOutlier(double val,double []tabIn,double[]mask,double eRatio) {
    	double factorB=1.4826;//see reference Leys at al. JESP 2013 - Detecting outliers:
    	double[]tabStats=MADeStatsDoubleSided(tabIn, mask);
    	if(tabStats.length==1)return new Object[] {true,val};
    	double madDown=tabStats[0]-tabStats[1];
    	double madUp=tabStats[2]-tabStats[0];
    	if(val<tabStats[0]-eRatio*factorB*madDown)return new Object[] {true,tabStats[0],tabStats[0]};
    	if(val>tabStats[0]+eRatio*factorB*madUp)return new Object[] {true,tabStats[0],tabStats[0]};
    	return(new Object[] {false,val,tabStats[0]});
    }
	
    public static double[] MADeStatsDoubleSided(double[] tabIn,double []mask) {
		if (tabIn.length==0)return null;
		if(mask==null) {mask=new double[tabIn.length];for(int i=0;i<tabIn.length;i++)mask[i]=1;}
		List<Double> l=new ArrayList<Double>();
		for(int i=0;i<tabIn.length;i++)if(mask[i]>0)l.add(tabIn[i]);
		Double []tab=(Double[])(l.toArray(new Double[l.size()]));
	
		if(tab.length<=5)return new double[] {0,0,0};
		Arrays.sort(tab);
		//System.out.print("[");
	//for(int i=0;i<tab.length;i++)System.out.print(tab[i]+" ,");
	//System.out.println("]");
	    double valMed=tab[tab.length/2];
	    double valMedUp=0;
	    double valMedDown=0;
	    if (tab.length%2==0)valMed=(tab[tab.length/2-1]+tab[tab.length/2])/2.0;
	
	    if (tab.length%4==0) {
	    	valMedDown=(tab[tab.length/4-1]+tab[tab.length/4])/2.0;
	    	valMedUp=(tab[3*tab.length/4-1]+tab[3*tab.length/4])/2.0;
	    }
	    if (tab.length%4==1) {
	    	valMedDown=tab[tab.length/4];
	    	valMedUp=tab[3*tab.length/4];
	    }
	    if (tab.length%4==2) {
	    	valMedDown=tab[tab.length/4];
	    	valMedUp=tab[3*tab.length/4];
	    }
	    if (tab.length%4==3) {
	    	valMedDown=(tab[tab.length/4+1]+tab[tab.length/4])/2.0;
	    	valMedUp=(tab[3*tab.length/4-1]+tab[3*tab.length/4])/2.0;
	    }
	    return new double[] {valMed,valMedDown,valMedUp};
    }

    public static double[] getQuartiles(double[] tabInTmp,double []mask) {
	  	if (tabInTmp.length==0)return null;
		if(mask==null) {mask=new double[tabInTmp.length];for(int i=0;i<tabInTmp.length;i++)mask[i]=1;}
		List<Double> lTmp=new ArrayList<Double>();
		for(int i=0;i<tabInTmp.length;i++)if(mask[i]>0)lTmp.add(tabInTmp[i]);
		Double []tabIn=(Double[])(lTmp.toArray(new Double[lTmp.size()]));
		
		int N=tabIn.length;
		if (tabIn.length==0)return new double[] {0,0,0};
		List<Double> l=new ArrayList<Double>();
		for(int i=0;i<tabIn.length;i++)l.add(tabIn[i]);
		Double []tab=(Double[])(l.toArray(new Double[l.size()]));
	
		if(tab.length<=5)return new double[] {0,0,0};
		Arrays.sort(tab);
		return new double[] {tab[N/4],tab[N/2],tab[(3*N)/4]};
    }

    public ImagePlus replaceMapsOutliersSlicePerSlice(int oneForTukeyFenceTwoForMADe,double nStdDev,int blockHalfSize,boolean informUser) {
	   if(oneForTukeyFenceTwoForMADe<1 || oneForTukeyFenceTwoForMADe>2)return VitimageUtils.nullImage(this.getMask());
	   if(nStdDev<0) {IJ.showMessage("Please choose a positive value for outliers threshold");return VitimageUtils.nullImage(this.getMask());}
   
	   ImagePlus mapsTemp=new Duplicator().run(this.hyperMaps,1,nMaps,1,Z,1,T);
	   ImagePlus mapsTempOut=new Duplicator().run(this.hyperMaps,1,nMaps,1,Z,1,T);
	   ImagePlus maskChange=VitimageUtils.nullImage(new Duplicator().run(this.hyperMaps,1,1,1,Z,1,T));
	   if(VitimageUtils.isNullImage(new Duplicator().run(this.hyperMaps,1,1,1,1,1,1))) {
		   IJ.showMessage("Maps image is not yet computed. Compute it first");
		   return VitimageUtils.nullImage(this.getMask());
	   }
	   
	   Object[][]objs=new Object[nMaps-1][];
	   double[]vals=new double[nMaps-1];
	   double[][]tabs=new double[nMaps-1][];
	   double[][]tabsMask=new double[nMaps-1][];
	   float[][]valMaps=new float[nMaps][];
	   float[][]valMapsOut=new float[nMaps][];
	   float[]valMaskOut=null;
	   int[]counts=new int[nMaps];
	   int countGlob=0;
	   int tot=T*Z;
	   int percent=Math.max(1, tot/10);
	   for(int t=0;t<T;t++) {
		   for(int z=0;z<Z;z++) {
			   if((z*T+t)%percent==0) {
				   if(informUser)IJ.log("Outlier processing : "+(VitimageUtils.dou((z*T+t)/(0.01*tot)))+"%");
				   if(informUser)IJ.showProgress((z*T+t)/(1.0*tot));
			   }
			   for(int map=0;map<nMaps;map++) {
				   valMaps[map]=(float[]) mapsTemp.getStack().getProcessor(VitimageUtils.getCorrespondingSliceInHyperImage(mapsTemp, map, z, t)).getPixels();
				   valMapsOut[map]=(float[]) mapsTempOut.getStack().getProcessor(VitimageUtils.getCorrespondingSliceInHyperImage(mapsTempOut, map, z, t)).getPixels();
				   valMaskOut=(float[]) maskChange.getStack().getProcessor(VitimageUtils.getCorrespondingSliceInHyperImage(maskChange, 0, z, t)).getPixels();
			   }
			   for(int x=0;x<X;x++) {
	    		   for(int y=0;y<Y;y++) {
	    			   if(valMaps[nMaps-1][X*y+x]<=0) continue;
	    			   counts[nMaps-1]++;
	    			   for(int map=0;map<nMaps-1;map++) {
	    				   double val=valMaps[map][X*y+x];
	    				   tabsMask[map]=VitimageUtils.valuesOfBlock(mapsTemp,x-blockHalfSize,y-blockHalfSize,z,x+blockHalfSize,y+blockHalfSize,z,nMaps-1,t);
	    				   tabs[map]=VitimageUtils.valuesOfBlock(mapsTemp,x-blockHalfSize,y-blockHalfSize,z,x+blockHalfSize,y+blockHalfSize,z,map,t);
	    				   vals[map]=val;
	    						   
	    				   objs[map]= (oneForTukeyFenceTwoForMADe==2 ? MADeIsOutlier(val, tabs[map], tabsMask[map], nStdDev) : tuckeyIsOutlier(val, tabs[map], tabsMask[map], nStdDev));
	    			   }
	    			   boolean isOutlier=false;
	    			   boolean []isOutlierTab=new boolean[nMaps-1];
	    			   for(int map=0;map<nMaps-1;map++) {
	    				   isOutlierTab[map]=(boolean) objs[map][0];
	    				   if( isOutlierTab[map]) {
	    					   counts[map]++;
	    					   isOutlier=true;
	    				   }
	    			   }
	    			   if(isOutlier) {
	    				   valMaskOut[X*y+x]=1;
	    				   countGlob++;
	    				   //System.out.println("\nWe got an outlier at x="+x+" y="+y+" z="+z+" t="+t);
	    				   for(int map=0;map<nMaps-1;map++) {
	        				   //System.out.println(" map "+map+"\nData was="+TransformUtils.stringVectorN(tabs[map],"")+"\nMask was="+TransformUtils.stringVectorN(tabsMask[map],""));
	        				   //System.out.println(" val=  "+vals[map]+" is outlier ? "+isOutlierTab[map]+" replaced by "+(double) objs[map][2]);
	        				   valMapsOut[map][X*y+x]=(float) ((double)(objs[map][2]));
	    				   }
	    			   }
	    		   }
			   }
		   }
	   }    	   
	   this.hyperMaps=mapsTempOut;
	   this.updateHyperImgFromEchoesAndMaps();
	   //  	   System.out.println("Bilan");
	   //  	   System.out.println("Countglob="+countGlob);
	   for(int map=0;map<nMaps-1;map++) {
		   //   		   System.out.println("Count out "+map+" = "+counts[map]);
		   }
		   maskChange.setDisplayRange(0, 1);
		   return maskChange;
   }

   public ImagePlus simulateOutlierRemoval(int x0,int y0,int x1,int y1,int z,int t){
	   int nConfigs=2;
	   int []configOutlier=new int[] {1,2};
	   String[]textConfigOutlier=new String[] {"Tukey","MADe"};

   int nStd=3;
   int []nStdDev=new int[] {3,5,7};
   String[]textStdDev=new String[] {"3*sigma","5*sigma","7*sigma"};
   
   int nNeighs=5;
   int []neighXY=new int[] {1,2,3,4,5};
   String[]textneighXY=new String[] {"radius=1","radius=2","radius=3","radius=4","radius=5"};

   int totalConf=nConfigs*nStd*nNeighs;
   int nC=this.nMaps;
   
   ImagePlus imgT1=new Duplicator().run(this.getAsImagePlus(),1,Ctot,z,z,t,t);
   VitimageUtils.copyImageCalibrationAndRange(imgT1,this.getAsImagePlus());
   ImagePlus imgT=VitimageUtils.cropMultiChannelFloatImage(imgT1, x0,x1, y0, y1, 0,0);
   VitimageUtils.copyImageCalibrationAndRange(imgT,this.getAsImagePlus());
   VitimageUtils.printImageResume(imgT,"imgT in HyperMap");
   HyperMap hypShort=new HyperMap(imgT);
   ImagePlus oldMap=  hypShort.getAsImagePlus();
   ImagePlus [][]retTab=new ImagePlus[nC][totalConf+1];
   double[][]ranges=new double[nMaps][2];
   
   for(int conf=0;conf<nConfigs;conf++) {
	   for(int st=0;st<nStd;st++) {
		   for(int nei=0;nei<nNeighs;nei++) {
			   double progress=(conf*nStd*nNeighs+st*nNeighs+nei)/(1.0*totalConf);
			   IJ.log("Outlier processing simulation : "+(VitimageUtils.dou(100*progress)+"%"));
			   IJ.showProgress(progress);

			   HyperMap hypTmp=HyperMap.hyperMapFactory(hypShort);
			   ImagePlus maskChange=hypTmp.replaceMapsOutliersSlicePerSlice(configOutlier[conf], nStdDev[st], neighXY[nei],false);
			   ImagePlus newMap=hypTmp.getAsImagePlus();
			   String text="Result "+textConfigOutlier[conf]+" - "+textStdDev[st]+" - "+textneighXY[nei];
			   for(int m=0;m<this.nMaps-1;m++) {
//						   System.out.println("conf="+conf+" st="+st+" nei="+nei+" m="+m);
				   retTab[m][1+conf*nNeighs*nStd+st*nNeighs+nei]=new Duplicator().run(newMap,m+1,m+1,1,1,1,1);						   
				   retTab[m][1+conf*nNeighs*nStd+st*nNeighs+nei].getStack().setSliceLabel(text, 1);
				   newMap.setC(m+1);
				   ranges[m]=new double[] {newMap.getDisplayRangeMin(),newMap.getDisplayRangeMax()};
			   }
			   ranges[nMaps-1]=new double[] {0,1};
			   retTab[nMaps-1][1+conf*nNeighs*nStd+st*nNeighs+nei]=new Duplicator().run(maskChange,1,1,1,1,1,1);	
			   retTab[nMaps-1][1+conf*nNeighs*nStd+st*nNeighs+nei].getStack().setSliceLabel(text, 1);
		   }
	   }
   }
   String text="Initial maps without outliers removal";
   for(int m=0;m<this.nMaps-1;m++) {
	   retTab[m][0]=new Duplicator().run(oldMap,m+1,m+1,1,1,1,1);
	   retTab[m][0].getStack().setSliceLabel(text, 1);
   }
   retTab[nMaps-1][0]=VitimageUtils.nullImage(retTab[0][0]);
   retTab[nMaps-1][0].getStack().setSliceLabel(text, 1);

   ImagePlus []hypTab=new ImagePlus[nC*(totalConf+1)];
   for(int i=0;i<nC;i++)for(int j=0;j<totalConf+1;j++) {
	   hypTab[i*(totalConf+1)+j]=retTab[i][j];
   }
   Concatenator con=new Concatenator();
	con.setIm5D(true);
	ImagePlus img2=con.concatenate(hypTab,false);
	img2=HyperStackConverter.toHyperStack(img2, nMaps, 1,totalConf+1,"xyztc","Grayscale");
		VitimageUtils.copyImageCalibrationAndRange(img2, imgT1);
		img2.setC(this.nMaps);
		img2.setDisplayRange(0, 1);
		return img2;
   }
   
       
	   
	   
	   
	   
	   
    
	/** Random helpers ---------------------------------------------------------*/
	public static String textAlg(int algType) {return algType==MRUtils.LM ? "LM" : "SIMP";}
	
	public static double meanRice(ImagePlus[]imgT1T2Line) {
		double[]meanRiceTab=new double[imgT1T2Line.length];
		for(int c=0;c<imgT1T2Line.length;c++) {		
			meanRiceTab[c]=MRUtils.readSigmaInSliceLabel(imgT1T2Line[c], 0, imgT1T2Line[c].getNSlices()/2, 0); 
		}
		return VitimageUtils.statistics1D(meanRiceTab)[0];
	}

	public void showCopy(String text) {
		ImagePlus temp=hyperImg.duplicate();
		temp.show();
		temp.setTitle(text);
	}
	
	public ImagePlus getCopy() {
		return hyperImg.duplicate();
	}
	
	
	
	
	
	
	
	
	/** Helpers for capillary measurements and hypermap normalization--------------------------------------------------------- */	
	public static double measureMeanCapillaryValueAlongZ(ImagePlus imgM0) {
		//Find capillary in the central slice
		ImagePlus img=new Duplicator().run(imgM0,1,1,imgM0.getNSlices()/2+1,imgM0.getNSlices()/2+1,1,1);
		int []coordsCentral=VitimageUtils.findCapillaryCenterInSlice(img,VitimageUtils.bionanoCapillaryRadius);
		
		//Find best match for the capillary in other slices, in a neighbourhood around, and measure its value
		double[]vals=VitimageUtils.capillaryValuesAlongZStatic(imgM0,coordsCentral,VitimageUtils.bionanoCapillaryRadius);
		System.out.println("Mean capillary values="+TransformUtils.stringVectorN(vals, ""));
		return VitimageUtils.statistics1D(vals)[0];
	}
		
	public double[]measureCapillaryValuesInM0Map(int time){
		//Find capillary in the central slice
		ImagePlus img=new Duplicator().run(hyperImg,1,1,Z/2+1,Z/2+1,time+1,time+1);
		int []coordsCentral=VitimageUtils.findCapillaryCenterInSlice(img,VitimageUtils.bionanoCapillaryRadius);
		
		//Find best match for the capillary in other slices, in a neighbourhood around, and measure its value
		return capillaryValuesAlongZ(time,coordsCentral,VitimageUtils.bionanoCapillaryRadius);
	}

	public double[]capillaryValuesAlongZ(int time,int[]coords,double capillaryRadius){
		IJ.log("\nCapillary detection at time="+time+".\nStarting from coordinates :"+coords[0]+", "+coords[1]);
		Duplicator dup=new Duplicator();
		int rayPix=(int)Math.round(capillaryRadius/(voxs[0]));
		int semiRayPix=rayPix/2;
		int lookupRadius=(int)Math.round(1.5*rayPix);

		int xMed=coords[0];int yMed=coords[1];int zMed=dims[2]/2;
		double[]capVals=new double[dims[2]];
		ImagePlus imgTemp,img3D;
		int xLast=xMed,yLast=yMed;
		long t0= System.currentTimeMillis();

		img3D=dup.run(hyperMaps,1,1,1,Z,time+1,time+1);
		for(int z=zMed;z>=0;z--) {
			//Extract a patch around the finding of upside
			imgTemp=VitimageUtils.cropImageFloat(img3D,xLast-lookupRadius,yLast-lookupRadius,z, lookupRadius*2,lookupRadius*2,1);
			//Find cap center in it
			int[]coordsNew=VitimageUtils.findCapillaryCenterInSlice(imgTemp,VitimageUtils.bionanoCapillaryRadius);
			
			//Update coordinates of last
			xLast=xLast-lookupRadius+coordsNew[0];
			yLast=yLast-lookupRadius+coordsNew[1];

			//Gather information of M0 value
			double[]stats=VitimageUtils.statistics1D(VitimageUtils.valuesOfBlock(img3D,xLast-semiRayPix, yLast-semiRayPix, z-1,xLast+semiRayPix, yLast+semiRayPix, z+1));
			IJ.log("Capillary detected at z="+z+" at coordinates "+xLast+", "+yLast+" with M0="+stats[0]+" std="+stats[1]);
			capVals[z]=stats[0];
		}
		xLast=xMed;
		yLast=yMed;
		for(int z=zMed;z<dims[2];z++) {
			//Extract a patch around the finding of upside
			imgTemp=VitimageUtils.cropImageFloat(img3D,xLast-lookupRadius,yLast-lookupRadius, z, lookupRadius*2,lookupRadius*2,1);
			
			//Find cap center in it
			int[]coordsNew=VitimageUtils.findCapillaryCenterInSlice(imgTemp,VitimageUtils.bionanoCapillaryRadius);
			
			//Update coordinates of last
			xLast=xLast-lookupRadius+coordsNew[0];
			yLast=yLast-lookupRadius+coordsNew[1];
			
			//Gather information of M0 value
			double[]stats=VitimageUtils.statistics1D(VitimageUtils.valuesOfBlock(img3D,xLast-semiRayPix, yLast-semiRayPix, z-1,xLast+semiRayPix, yLast+semiRayPix, z+1));
			System.out.println("Capillary detected at z="+z+" at coordinates "+xLast+", "+yLast+" with M0="+stats[0]+" std="+stats[1]);
			capVals[z]=stats[0];
		}
		return capVals;
	}



	

	
	
	/** Tidy needed below ---------------------------------------------------------*/
	/*
	public ImagePlus computeMapsAgainAndMaskOld(boolean recomputeMask,int algType,boolean  makeMono,boolean  makeMulti) {
		boolean debug=true;
		boolean makeBoutureTrick=hyperImg.getStack().getSliceLabel(1).contains("BOUTURE");
		if(hyperImg.getNFrames()>1) {VitiDialogs.notYet("ComputeMapsAgain of 5D hypermaps in HyperMap");return null;}
		int Z=hyperImg.getNSlices();
		int C=hyperImg.getNChannels();
		int T=hyperImg.getNFrames();
		int X=hyperImg.getWidth();
		int Y=hyperImg.getHeight();
		int bitD=hyperImg.getBitDepth();
		int nVox=Z*T*X*Y;
		if(bitD!=16 && bitD!=32) {IJ.showMessage("Unexepected data type : no 16 nor 32-bits");System.exit(0);}
		ImagePlus imgMask=null;
		ImagePlus echoes=(ImagePlus)(getEchoesImage(false)[0]);
		ImagePlus []imgT1T2Line=VitimageUtils.stacksFromHyperstackFastBis(echoes);
		double[]meanRiceTab=new double[imgT1T2Line.length];
		for(int c=0;c<imgT1T2Line.length;c++) {		
			meanRiceTab[c]=MRUtils.readSigmaInSliceLabel(imgT1T2Line[c], 0, Z/2, 0); 
		}
		double meanRiceOverT2=VitimageUtils.statistics1D(meanRiceTab)[0];
		if(!recomputeMask)imgMask=new Duplicator().run(hyperImg,4,4,1,Z,1,T);
		else {						//Compute mask
			ImagePlus max=null;
			if(bitD==16)max=VitimageUtils.maxOfImageArrayShort(new ImagePlus[][] {imgT1T2Line});
			else max=VitimageUtils.maxOfImageArrayDouble(new ImagePlus[][] {imgT1T2Line});
			imgMask=VitimageUtils.getFloatBinaryMask(max,MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*(makeBoutureTrick ? MRUtils.THRESHOLD_BOUTURE_FACTOR : 1 )*meanRiceOverT2,1E10);
			ImagePlus imgMaskCap=VitimageUtils.getAntiCapillaryMask(max,VitimageUtils.bionanoCapillaryRadius);			
			imgMask=VitimageUtils.makeOperationBetweenTwoImages(imgMask,imgMaskCap, 1, true);
			if(makeBoutureTrick)IJ.run(imgMask,"Median 3D...", "x=1.2 y=1.2 z=1.2");
		}
		IJ.run(imgMask,"32-bit","");
		for(int i=0;i<imgT1T2Line.length;i++)IJ.run(imgT1T2Line[i],"32-bit","");
	
		Timer t;
		ImagePlus[]	maps;
		if(makeMono) {
			t=new Timer();
			t.print("\nStarting "+algType);
			//Compute T1T2mono
			t.print("Start mono");
			maps=MRUtils.computeT1T2MapMultiThread(imgT1T2Line, sigmaInUse,/*MRUtils.T1T2_MONO_RICEMRUtils.T2_MONO_RICE,algType,debug);		
			t.print("End mono");
			System.out.println("Per voxel = "+(VitimageUtils.dou(t.getTime()*1000.0/nVox))) ;
			ImagePlus M0Mono=maps[0];
			
			ImagePlus T2Mono=maps[1];
			ImagePlus khi2Mono=maps[2];
//			ImagePlus T1Mono=maps[1];
//			String []ss=new String[] {"M0","T1","T2","Khi2"};
			String []ss=new String[] {"M0","T2","Khi2"};
			for(int i=0;i<3;i++) {maps[i].show();maps[i].setTitle(ss[i]+"-"+textAlg(algType));maps[i].getWindow().setSize(200,200);maps[i].getCanvas().fitToWindow();IJ.run(maps[i],"Fire","");}
			khi2Mono.setDisplayRange(-1, 3);khi2Mono.updateAndRepaintWindow();
//			T1Mono.setDisplayRange(0, MRUtils.maxDisplayedBionanoT1);T1Mono.updateAndRepaintWindow();
			T2Mono.setDisplayRange(0, MRUtils.maxDisplayedBionanoT2);T2Mono.updateAndRepaintWindow();
			M0Mono.setDisplayRange(0, MRUtils.maxDisplayedBionanoM0);M0Mono.updateAndRepaintWindow();
			//		ImagePlus deltaMono=maps[4];
		}
		
		
		//Compute T1T2bi
		if(makeMulti) {
			t=new Timer();
			t.print("Start multi");
			/////////if(!makeBoutureTrick)///////////////////////
			maps=MRUtils.computeT1T2MapMultiThread(imgT1T2Line, sigmaInUse,/*MRUtils.T1T2_MULTI_RICEMRUtils.T1T2_MULTI_RICE,algType,debug);
			t.print("End multi.");
			System.out.println("Per voxel = "+(VitimageUtils.dou(t.getTime()*1000.0/nVox))) ;
			ImagePlus M0Multi1=maps[0];
			ImagePlus M0Multi2=maps[1];
			ImagePlus T1Multi=maps[2];
			ImagePlus T2Multi1=maps[3];
			ImagePlus T2Multi2=(!makeBoutureTrick) ? maps[4] : maps[3];
			ImagePlus khi2Multi=(!makeBoutureTrick) ? maps[5] : maps[3];
			String []ss2=new String[] {"M01MULTI","M02MULTI","T1MULTI","T21MULTI","T22MULTI","KHI2MULTI"};
			for(int i=0;i<6;i++) {maps[i].show();maps[i].setTitle(ss2[i]+"-"+textAlg(algType));maps[i].getWindow().setSize(200,200);maps[i].getCanvas().fitToWindow();IJ.run(maps[i],"Fire","");}
			khi2Multi.setDisplayRange(-1, 3);
	
			ImagePlus sum1=VitimageUtils.makeOperationBetweenTwoImages(M0Multi1, T2Multi1, 2,true);
			ImagePlus sum2=VitimageUtils.makeOperationBetweenTwoImages(M0Multi2, T2Multi2, 2,true);
			ImagePlus sum3=VitimageUtils.makeOperationBetweenTwoImages(sum1, sum2, 1, true);
			ImagePlus sum4=VitimageUtils.makeOperationBetweenTwoImages(M0Multi1, M0Multi2, 1, true);
			ImagePlus result=VitimageUtils.makeOperationBetweenTwoImages(sum3, sum4, 3, true);
			result.show();result.setTitle("T2globMulti"+"-"+textAlg(algType));IJ.run(result,"Fire","");
			result.setDisplayRange(0, 150);
			result.getWindow().setSize(200,200);result.getCanvas().fitToWindow();
			sum4.show();sum4.setTitle("M0globMulti"+"-"+textAlg(algType));IJ.run(sum4,"Fire","");
			sum4.getWindow().setSize(200,200);sum4.getCanvas().fitToWindow();
			T1Multi.setDisplayRange(0, MRUtils.maxDisplayedBionanoT1);
			T2Multi1.setDisplayRange(0, MRUtils.maxDisplayedBionanoT2);
			T2Multi1.setDisplayRange(0, MRUtils.maxDisplayedBionanoT2);
			M0Multi1.setDisplayRange(0, MRUtils.maxDisplayedBionanoM0);
			M0Multi2.setDisplayRange(0, MRUtils.maxDisplayedBionanoM0);
			//		ImagePlus deltaMulti=maps[6];
		}

		ImagePlus []tempRes=new ImagePlus[C];
		tempRes[0]=retM0;
		tempRes[1]=retT1;
		tempRes[2]=retT2;
		tempRes[3]=imgMask;
		for(int c=4;c<C;c++) {tempRes[c]=imgT1T2Line[c-4];}
		for(int c=0;c<C;c++)tempRes[c]=VitimageUtils.convertFloatToShortWithoutDynamicChanges(tempRes[c]);			
		imgMask=VitimageUtils.convertFloatToShortWithoutDynamicChanges(imgMask);
		ImagePlus newHyperImg=Concatenator.run(tempRes);
		newHyperImg=HyperStackConverter.toHyperStack(newHyperImg, C,Z,T,"xyztc","Grayscale");

		for(int c=0;c<hyperImg.getNChannels();c++)for(int z=0;z<imgT1T2Line[0].getNSlices();z++) {
			newHyperImg.getStack().setSliceLabel(hyperImg.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, 0)), VitimageUtils.getCorrespondingSliceInHyperImage(newHyperImg, c, z, 0));
		}
		hyperImg=newHyperImg;
	
		//Update capillaryvalues
		short[]capVals;
		int nHits;
		int incr;
		int slice=0;
		double sum;
		for(int t=0;t<T;t++) {
			for(int z=0;z<Z;z++) {
				ArrayList<int[]> coordsPointsCap=new ArrayList<int[]>();
				capVals=(short[])imgMask.getStack().getProcessor(z+1).getPixels();
				nHits=0;
				for(int x=0;x<X;x++)for(int y=0;y<Y;y++) if(( capVals[VitimageUtils.getPixelIndex(x, X, y)] & 0xffff ) >4.5) {nHits++; coordsPointsCap.add(new int[] {x,y});};	
				for(int c=0;c<C;c++) {
					double[]vals=new double[nHits];
					incr=0;
					slice=VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, t);
					capVals=(short[])hyperImg.getStack().getProcessor(slice).getPixels();				
					for(int i=0;i<nHits;i++) {
						vals[incr++]=(int)((capVals[VitimageUtils.getPixelIndex(coordsPointsCap.get(i)[0], X, coordsPointsCap.get(i)[1])]) & 0xffff);
					}
					double[]stats=VitimageUtils.statistics1D(vals);
					String label=hyperImg.getStack().getSliceLabel(slice);
					String []infos=label.split("_");
					label="";
					for(int i=0;i<infos.length;i++) {
						if(infos[i].contains("CAP="))infos[i]=("CAP="+VitimageUtils.dou(stats[0])+"+-"+VitimageUtils.dou(stats[1]));
						label+=(infos[i])+( (i==infos.length-1) ?  "" : "_");
					}
					hyperImg.getStack().setSliceLabel(label,slice);
				}
			}
		}
		
		
		//Update ranges and luts
		hyperImg.setC(1);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoM0 );

		hyperImg.setC(2);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoT1 );

		hyperImg.setC(3);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );

		hyperImg.setC(4);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,5 );

		for(int c=5;c<=C;c++) {
			hyperImg.setC(c);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoM0 );
		}
		return hyperImg;


	/*	
		
		//Combine
		ImagePlus retM0=M0Mono.duplicate();
//		ImagePlus retDelta=deltaMono.duplicate();
		ImagePlus retT2=T2Mono.duplicate();
		ImagePlus retT1=T1Mono.duplicate();
		ImagePlus selectedModel=imgMask.duplicate();
		ImagePlus portionMap=imgMask.duplicate();
		int index;
		for(int z=0;z<dims[2];z++) {
			float[]valsRetM0=(float[]) retM0.getStack().getProcessor(z+1).getPixels();
			float[]valsRetT2=(float[]) retT2.getStack().getProcessor(z+1).getPixels();
			float[]valsRetT1=(float[]) retT1.getStack().getProcessor(z+1).getPixels();
			float[]valsM0=(float[]) M0Mono.getStack().getProcessor(z+1).getPixels();
			float[]valsM0Multi1=(float[]) M0Multi1.getStack().getProcessor(z+1).getPixels();
			float[]valsM0Multi2=(float[]) M0Multi2.getStack().getProcessor(z+1).getPixels();
			float[]valsT1=(float[]) T1Mono.getStack().getProcessor(z+1).getPixels();
			float[]valsT1Multi=(float[]) T1Multi.getStack().getProcessor(z+1).getPixels();
			float[]valsT2=(float[]) T2Mono.getStack().getProcessor(z+1).getPixels();
			float[]valsT2Multi1=(float[]) T2Multi1.getStack().getProcessor(z+1).getPixels();
			float[]valsT2Multi2=(float[]) T2Multi2.getStack().getProcessor(z+1).getPixels();
			float[]valMask=(float[]) imgMask.getStack().getProcessor(z+1).getPixels();
			float[]valKhi2=(float[]) khi2Mono.getStack().getProcessor(z+1).getPixels();
			float[]valKhi2Multi=(float[]) khi2Multi.getStack().getProcessor(z+1).getPixels();
			float[]valSelectedModel=(float[]) selectedModel.getStack().getProcessor(z+1).getPixels();
			float[]valPortionMap=(float[]) portionMap.getStack().getProcessor(z+1).getPixels();
			float[]valDelta=(float[]) retDelta.getStack().getProcessor(z+1).getPixels();
			float[]valDeltaMono=(float[]) deltaMono.getStack().getProcessor(z+1).getPixels();
			float[]valDeltaMulti=(float[]) deltaMulti.getStack().getProcessor(z+1).getPixels();*/
		
		
		/*
			for(int y=0;y<dims[1];y++) {
				for(int x=0;x<dims[0];x++) {
					index=y*dims[0]+x;
					if((int)Math.round(valMask[index])%2==0) {
						valsRetM0[index]=0;
						valsRetT1[index]=0;
						valsRetT2[index]=0;
						//valDelta[index]=0;
						continue;//If out of mask, set to 0  (work already done)
					}
					if(makeBoutureTrick)continue;
					if(valKhi2Multi[index]>=MRUtils.ERROR_KHI2) {
						if(valKhi2[index]<MRUtils.ERROR_KHI2) {
							//if multicomp is out, and monocomp khi < errorkhi, set to monocomp (already done)
						}
						else {
							//else set to 0
							valsRetM0[index]=0;
							valsRetT1[index]=0;
							valsRetT2[index]=0;
							//valDelta[index]=0;
						}
					}
					else {
						if(valKhi2Multi[index]<valKhi2[index]) {
							//if T2 bi less Khi, // eventually : and T2 vals > 10  , < 500 and ratio between components > 0.02, < 50
							valsRetM0[index]=valsM0Multi1[index]+valsM0Multi2[index];
							valsRetT1[index]=valsT1Multi[index];
							//valDelta[index]=valDeltaMulti[index];
//							valsRetT2[index]=(valsT2Multi1[index]*valsM0Multi1[index]+valsT2Multi2[index]*valsM0Multi2[index])/(valsM0Multi1[index]+valsM0Multi2[index]);
							valSelectedModel[index]=2;
							valPortionMap[index]=valsM0Multi1[index]/(valsM0Multi1[index]+valsM0Multi2[index]);
						}
						else { 
							if(valKhi2[index]<MRUtils.ERROR_KHI2) {
								//if T1 khi < errorkhi, set to T1
								valsRetM0[index]=valsM0[index];
								valsRetT1[index]=valsT1[index];
								//valDelta[index]=valDeltaMono[index];
								valsRetT2[index]=valsT2[index];
							}
							else {
								//else set to 0
								valsRetM0[index]=0;
								valsRetT1[index]=0;
								//valDelta[index]=0;
								valsRetT2[index]=0;
							}
						}
					}
				}
			}
		}
		
		//retDelta.show();
		//retDelta.setTitle("Composite delta");
		ImagePlus []tempRes=new ImagePlus[C];
		tempRes[0]=retM0;
		tempRes[1]=retT1;
		tempRes[2]=retT2;
		tempRes[3]=imgMask;
		for(int c=4;c<C;c++) {tempRes[c]=imgT1T2Line[c-4];}
		for(int c=0;c<C;c++)tempRes[c]=VitimageUtils.convertFloatToShortWithoutDynamicChanges(tempRes[c]);			
		imgMask=VitimageUtils.convertFloatToShortWithoutDynamicChanges(imgMask);
		ImagePlus newHyperImg=Concatenator.run(tempRes);
		newHyperImg=HyperStackConverter.toHyperStack(newHyperImg, C,Z,T,"xyztc","Grayscale");

		for(int c=0;c<hyperImg.getNChannels();c++)for(int z=0;z<imgT1T2Line[0].getNSlices();z++) {
			newHyperImg.getStack().setSliceLabel(hyperImg.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, 0)), VitimageUtils.getCorrespondingSliceInHyperImage(newHyperImg, c, z, 0));
		}
		hyperImg=newHyperImg;
	
		//Update capillaryvalues
		short[]capVals;
		int nHits;
		int incr;
		int slice=0;
		double sum;
		for(int t=0;t<T;t++) {
			for(int z=0;z<Z;z++) {
				ArrayList<int[]> coordsPointsCap=new ArrayList<int[]>();
				capVals=(short[])imgMask.getStack().getProcessor(z+1).getPixels();
				nHits=0;
				for(int x=0;x<X;x++)for(int y=0;y<Y;y++) if(( capVals[VitimageUtils.getPixelIndex(x, X, y)] & 0xffff ) >4.5) {nHits++; coordsPointsCap.add(new int[] {x,y});};	
				for(int c=0;c<C;c++) {
					double[]vals=new double[nHits];
					incr=0;
					slice=VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, c, z, t);
					capVals=(short[])hyperImg.getStack().getProcessor(slice).getPixels();				
					for(int i=0;i<nHits;i++) {
						vals[incr++]=(int)((capVals[VitimageUtils.getPixelIndex(coordsPointsCap.get(i)[0], X, coordsPointsCap.get(i)[1])]) & 0xffff);
					}
					double[]stats=VitimageUtils.statistics1D(vals);
					String label=hyperImg.getStack().getSliceLabel(slice);
					String []infos=label.split("_");
					label="";
					for(int i=0;i<infos.length;i++) {
						if(infos[i].contains("CAP="))infos[i]=("CAP="+VitimageUtils.dou(stats[0])+"+-"+VitimageUtils.dou(stats[1]));
						label+=(infos[i])+( (i==infos.length-1) ?  "" : "_");
					}
					hyperImg.getStack().setSliceLabel(label,slice);
				}
			}
		}
		
		
		//Update ranges and luts
		hyperImg.setC(1);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoM0 );

		hyperImg.setC(2);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoT1 );

		hyperImg.setC(3);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,MRUtils.maxDisplayedT2 );

		hyperImg.setC(4);
		IJ.run(hyperImg,"Fire","");
		hyperImg.setDisplayRange(0,5 );

		for(int c=5;c<=C;c++) {
			hyperImg.setC(c);
			IJ.run(hyperImg,"Fire","");
			hyperImg.setDisplayRange(0,MRUtils.maxDisplayedBionanoM0 );
		}
		return hyperImg;
		
		
	}
	
*/
	
	


	
	/* Capillary utilities. The first normalize MR signal along Z and T axis, using the capillary values
	public void normalizeMRSignal() {
		double[][]capillaryValuesInM0Map=new double[T][dims[2]];
		double factor=1;
		//Measure capillary values in T1 and T2 sequences, at each time
		//Adjust image dynamics of M0 map, T1seq, and T2 seq in order that value being 10000;
		for(int t=0;t<this.T;t++) {
			capillaryValuesInM0Map[t]=measureCapillaryValuesInM0Map(t);
			for(int c=0;c<this.C;c++) {
				for(int z=0;z<Z;z++) {
					int slice=VitimageUtils.getCorrespondingSliceInHyperImage(hyperEchoes, c, z, t);
					factor=MRUtils.maxM0BionanoForNormalization/capillaryValuesInM0Map[t][z];
					ImageProcessor ip=hyperImg.getStack().getProcessor(slice);
					ip.multiply(factor);
					if((c)!=m0MapChannel)MRUtils.modifySigma(hyperImg,c,z,t,factor*MRUtils.readSigmaInSliceLabel(hyperImg, c, z, t));
				}
				hyperImg.setC(c+1);
				hyperImg.setDisplayRange(0, MRUtils.maxDisplayedBionanoM0);
			}
		}		
	}
	*/	
	

	
	
/*	public int getNumberOfMaps() {
		int nMaps=0;
		for(int i=0;i<hyperImg.getNChannels();i++) {
			if(hyperImg.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyperImg, i, 0, 0)).contains("MAP")) {
				nMaps=i+1;		
			}
		}
		return nMaps;
	}
}
	public static int getNumberOfMaps(ImagePlus hyper) {
		int nMaps=0;
		for(int i=0;i<hyper.getNChannels();i++) {
			if(hyper.getStack().getSliceLabel(VitimageUtils.getCorrespondingSliceInHyperImage(hyper, i, 0, 0)).contains("MAP")) {
				nMaps=i+1;		
			}
		}
		return nMaps;
	}
	*/

   
	/*	
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();


		String dir="/home/fernandr/Bureau/A_Test/Test_boundaries/Tsmall/";
		String[]tissues=new String[] {"testSmall.tif","Tfull.tif","Tlow.tif","Tout.tif","Tvais.tif"};
		for(int tis=0;tis<1;tis++) {
			System.out.println(dir+tissues[tis]);
			//ImagePlus img=IJ.openImage(dir+tissues[tis]);
						ImagePlus img=IJ.openImage("/home/fernandr/Bureau/A_Test/Test_boundaries/test4.tif");
			HyperMap hyp=new HyperMap(img);
			img.show();
			boolean makeMono=true;
			boolean makeMulti=false;
			boolean makeSimplex=true;
			boolean makeLM=false;
			if(makeSimplex)hyp.computeMapsAgainAndMask(true,MRUtils.SIMPLEX,makeMono,makeMulti);
			if(makeLM)hyp.computeMapsAgainAndMask(true,MRUtils.LM,makeMono,makeMulti);
			
			if(makeMono && makeSimplex && makeLM) {
				String[]namesMono=new String[] {"T1","T2","M0","Khi2"};
				for(int i=0;i<namesMono.length;i++) {
					ImagePlus imgSimp= WindowManager.getImage(namesMono[i]+"-SIMP");
					ImagePlus imgLm=WindowManager.getImage(namesMono[i]+"-LM");
					ImagePlus diff=VitimageUtils.makeOperationBetweenTwoImages(imgSimp, imgLm, 4, true);
					diff.show();diff.setTitle(namesMono[i]+"-DIFF");IJ.run(diff,"Fire","");
					diff.getWindow().setSize(200,600);diff.getCanvas().fitToWindow();
				}
			}
			if(makeMulti && makeSimplex && makeLM) {
				String[]namesMulti=new String[] {"T1MULTI","T21MULTI","T22MULTI","M01MULTI","M02MULTI","M0globMulti","T2globMulti","KHI2MULTI"};
				for(int i=0;i<namesMulti.length;i++) {
					ImagePlus imgSimp=WindowManager.getImage(namesMulti[i]+"-SIMP");
					ImagePlus imgLm=WindowManager.getImage(namesMulti[i]+"-LM");
					ImagePlus diff=VitimageUtils.makeOperationBetweenTwoImages(imgSimp, imgLm, 4, true);
					diff.show();diff.setTitle(namesMulti[i]+"-DIFF");IJ.run(diff,"Fire","");
					diff.getWindow().setSize(200,600);diff.getCanvas().fitToWindow();
				}
			}
			
			
			VitimageUtils.waitFor(4000000);
		}
		System.exit(0);
		//		HyperMap hyp=new HyperMap(IJ.openImage("/home/fernandr/Bureau/Traitements/Sorgho/Test SSM1/Out_input_Fijiyama/Output/Exported_data/Data_combined.tif"));
	/*	hyp.showCopy("Initial Hyper MRI");
		System.out.println("Ok1");
		ImagePlus img=hyp.computeMapsAgainAndMask(true);
		img.show();++
		System.out.println("Ok3");
		hyp.showCopy("Normalized Hyper MRI");
	}
*/
	
	
	
	
	
	
	

}