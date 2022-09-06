package net.imagej.fijirelax.testing;

import java.io.File;


import net.imagej.fijiyama.common.VitimageUtils;
import ij.IJ;
import ij.ImageJ;
import ij.plugin.frame.PlugInFrame;

public class TestingFijiRelax_Package extends PlugInFrame {

	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public TestingFijiRelax_Package() {
		super("");
		IJ.log("Here 1");
	}

	
	public TestingFijiRelax_Package(String title) {
		super(title);
		IJ.log("Here 1");
	}

	public static boolean isInteger(String s) {
		try {
			int a=Integer.parseInt(s);	
		}catch (Exception e) {return false;}
		return true;
	}
	
	public static int safeParsePositiveInt(String s) {
		int ERR_CODE=-999999999;
		int ret=ERR_CODE;
		try {
			ret=Integer.parseInt(s);	
		}catch (Exception e) {int a=1;}
		if(ret!=ERR_CODE)return ret;
		
		ret=0;
		int nbChars=s.length();
		boolean isInteger=false;
		int index=0;
		while((index<(s.length()-1) && !isInteger(""+s.charAt(index)))){index++;}
		while((index<(s.length()-1) && isInteger(""+s.charAt(index)))){
			String s2=""+s.charAt(index);
			ret*=10;
			ret+=Integer.parseInt(s2);
			index++;
		}
		return ret;
	}
	
	public void run(String arg) {
		IJ.log("Here 2");
		testNii();
	}
	
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		TestingFijiRelax_Package test=new TestingFijiRelax_Package("");
		test.testNii();
	}
	
	public void testNii() {
		IJ.log("Here 3");
		String pathTest="/home/rfernandez/Bureau/A_Test/FijiRelax_DOI/Experiments_for_figures_and_tables/Input_data/Images/Human_brain_nifti/SEdata.nii";
	//	ImagePlus imp = IJ.openImage("/home/rfernandez/Bureau/A_Test/FijiRelax_DOI/Experiments_for_figures_and_tables/Input_data/Images/Human_brain_nifti/SEdata.nii");
		
//		ImageJ ij=new ImageJ();
		//Nifti_Reader nif;
		//System.out.println(new File(pathTest).exists());
//		ImagePlus img=IJ.openImage(pathTest);
//		img.show();
		//Image img=SimpleITK.readImage(pathTest);
//		ImageFileReader reader = new ImageFileReader();
//		reader.get
		//reader.setImageIO("NiftiImageIO");
//		  reader.setFileName(pathTest);
//		  Image simg= reader.execute();
		//run("Bio-Formats Importer", "open=[" + id + "] autoscale color_mode=Composite rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");

		  
		//ImagePlus simgIJ=ItkImagePlusInterface.itkImageToImagePlus(img);
		//simgIJ.show();
		IJ.log("Here 4");
	}
	
	public static void testStomata() {
		String s="D:\\Matt\\Stomata project\\Grasses\\Experiment 5\\MGX\\ABA\\ABA 1\\ABA_1_1\\Hv_ABA_1_1_SC.tif";
		String s2="BLABLA\\pouet";
		System.out.println(s2);
		System.out.println(s2.replace("\\", "/"));
		System.exit(0);
		for(int i=0;i<s.length();i++)System.out.println(s.charAt(i));
		System.out.println("["+safeParsePositiveInt("0010")+"]");
		String pathToFjmFile="/home/fernandr/Bureau/A_Test/Matthew/Fijiyama_series_2021-02-08_02-50.fjm";
		//String pathToFjmFile="/home/fernandr/Bureau/A_Test/Matthew/Fijiyama_series_2021-01-28_11-30.fjm";
		
			System.out.println("Starting setup");
			//Verify the configuration file, and detect registration style: serie or two images
			if(pathToFjmFile==null || (!(new File(pathToFjmFile).exists())) ||
					(!(pathToFjmFile.substring(pathToFjmFile.length()-4,pathToFjmFile.length())).equals(".fjm"))) {
				return;
			}
			String []names = new String[] {new File(pathToFjmFile).getParent(),new File(pathToFjmFile).getName()};	
			String name=names[1].substring(0,names[1].length()-4);
			String  serieOutputPath=names[0];

			
			String fjmFileContent=VitimageUtils.readStringFromFile(new File(names[0],names[1]).getAbsolutePath());

			//Get parameters from file
			String[]lines=fjmFileContent.split("\n");
			name=lines[1].split("=")[1];
			//this.serieOutputPath=(lines[3].split("=")[1]);
			String serieInputPath=(lines[4].split("=")[1]);
			int step=safeParsePositiveInt(lines[7].split("=")[1]);
			int nSteps=safeParsePositiveInt(lines[8].split("=")[1]);
			System.out.println(step);
			
/*			if(! isSerie) {
				fijiyamaGui.mode=Fijiyama_GUI.MODE_TWO_IMAGES;
				this.nTimes=1;
				this.nMods=2;
				this.referenceTime=0;this.referenceModality=0;
				this.refTime=0;this.refMod=0;
				this.movTime=0;this.movMod=1;
				this.setupStructures();
				this.paths[0][0]=lines[5].split("=")[1];
				this.paths[0][1]=lines[6].split("=")[1];
				fijiyamaGui.modeWindow=Fijiyama_GUI.WINDOWTWOIMG;
			}
			else {
				fijiyamaGui.mode=Fijiyama_GUI.MODE_SERIE;
				fijiyamaGui.modeWindow=Fijiyama_GUI.WINDOWSERIERUNNING;
				this.referenceTime=Integer.parseInt(lines[5].split("=")[1]);
				this.referenceModality=Integer.parseInt(lines[6].split("=")[1]);
				//Read modalities
				this.nMods=Integer.parseInt(lines[10].split("=")[1]);
				if(this.nMods>1) {
					this.mods=new String[this.nMods];
					for(int i=0;i<this.nMods;i++)this.mods[i]=lines[11+i].split("=")[1];
				}
				
				//Read times and nSteps
				this.nTimes=Integer.parseInt(lines[12+this.nMods].split("=")[1]);
				if(this.nTimes>1) {
					this.times=new String[this.nTimes];
					for(int i=0;i<this.nTimes;i++) {
						this.times[i]=lines[13+this.nMods+i].split("=")[1];
					}
				}
				this.expression=lines[14+this.nMods+this.nTimes].split("=")[1];
				this.pathToMask=lines[15+this.nMods+this.nTimes].split("=")[1];
				if(this.pathToMask.equals("None"))this.maskImage=null;
				else this.maskImage=IJ.openImage(this.pathToMask);

				this.setupStructures();

				//Affecter les paths
				for(int nt=0;nt<this.nTimes;nt++) {
					for(int nm=0;nm<this.nMods;nm++) {
						File f=new File(this.serieInputPath,expression.replace("{Time}", times[nt]).replace("{ModalityName}",mods[nm]));
						IJ.log("Series images lookup : checking existence of image "+f.getAbsolutePath());
						if(f.exists()) {
							IJ.log("       Found.");
							this.paths[nt][nm]=f.getAbsolutePath();
						}
						else IJ.log("      Not found.");
					}
				}
				IJ.log("");
			}
		
			//Check computer capacity and get a working copy of images
			this.checkComputerCapacity(true);
			if(!this.openImagesAndCheckOversizing())return false;
			ItkTransform trTemp = null;
			
			//Read the steps : all RegistrationAction serialized .ser object, and associated transform for those already computed
			File dirReg=new File(this.serieOutputPath,"Registration_files");		
			for(int st=0;st<this.nSteps;st++) {
				File f=new File(dirReg,"RegistrationAction_Step_"+st+".ser");
				RegistrationAction regTemp=RegistrationAction.readFromTxtFile(f.getAbsolutePath());
				if(regTemp.isDone()) {
					f=new File(dirReg,"Transform_Step_"+st+(((regTemp.typeAction!=1) || (regTemp.typeTrans!=Transform3DType.DENSE))?".txt":".transform.tif"));
					IJ.log("Transformation lookup : "+f+" type "+regTemp.typeTrans);
					
					if(regTemp.typeTrans!=Transform3DType.DENSE || regTemp.typeAction!=1)trTemp=ItkTransform.readTransformFromFile(f.getAbsolutePath());
					else trTemp=ItkTransform.readAsDenseField(f.getAbsolutePath());				
				}
				addTransformAndActionBlindlyForBuilding(trTemp,regTemp);
			}
			currentRegAction=regActions.get(step);
			this.unit=this.images[this.referenceTime][this.referenceModality].getCalibration().getUnit();
			System.out.println("Finishing setup");
			return true;
		}
*/	
	}
/*
	public static void testDims() {
		String pathSource="/home/fernandr/Bureau/Traitements/Sorgho/Donnees_brutes_export_Romain";
		String tempStr=pathSource;
		for (String str1 : VitimageUtils.stringArraySort( new File(tempStr).list())) {
			System.out.println("\n\n\nSpecimen="+str1);
			tempStr=pathSource+"/"+str1;
			for (String str2 : VitimageUtils.stringArraySort( new File(tempStr).list())) {
				System.out.println("|\n|--Experiment="+str2);
				tempStr=pathSource+"/"+str1+"/"+str2;
				for (String str3 : VitimageUtils.stringArraySort( new File(tempStr).list())) {
					tempStr=pathSource+"/"+str1+"/"+str2+"/"+str3;
					String str4 = VitimageUtils.stringArraySort( new File(tempStr).list())[0];
					tempStr=pathSource+"/"+str1+"/"+str2+"/"+str3+"/"+str4;
					String str5 = VitimageUtils.stringArraySort( new File(tempStr).list())[0];
					tempStr=pathSource+"/"+str1+"/"+str2+"/"+str3+"/"+str4+"/"+str5;
//					VitimageUtils.printImageResume(IJ.openImage(tempStr));
					double []dims=VitimageUtils.getVoxelSizes(IJ.openImage(tempStr));
					double vol=VitimageUtils.getVoxelVolume(IJ.openImage(tempStr));
					System.out.println("|----"+str2+" "+str3+" "+str4+" VoxelVolume="+VitimageUtils.dou(1000*vol)+" .10^-3 mm3     "+"Voxelsizes = "+dims[0]+" X "+dims[1]+" X "+dims[1]);
				}	
			}			
		}
	}
	
	
	public static void testNorm() {
		
		ImagePlus[]tab=new ImagePlus[4];
		tab[0]=IJ.openImage("/home/fernandr/Bureau/Temp/m0.tif");
		HyperMap.measureMeanCapillaryValueAlongZ(tab[0]);
		VitimageUtils.waitFor(1000000);
		tab[1]=IJ.openImage("/home/fernandr/Bureau/Temp/img22.tif");
		tab[2]=IJ.openImage("/home/fernandr/Bureau/Temp/img33.tif");
		tab[3]=IJ.openImage("/home/fernandr/Bureau/Temp/img44.tif");
		VitimageUtils.waitFor(10000000);
	}
	
	public static void chiasse() {
		String dirIn="/home/fernandr/Bureau/Traitements/Sorgho/Cartes_calculees_methode_Romain";
		boolean makeShort=false;
		String dirOut="";
		dirOut="/home/fernandr/Bureau/Traitements/Sorgho/Cartes_rangees_par_specimen";
		if(makeShort)dirOut="/home/fernandr/Bureau/Traitements/Sorgho/Cartes_short";
		String[]specs=new File(dirIn).list();
		for(String sp : specs) {
			System.out.println("Processing "+sp);
			String[]infoTab=sp.split("_");
			String strNew=infoTab[0]+"_"+infoTab[2]+".tif";
			ImagePlus img=IJ.openImage(new File(dirIn,sp).getAbsolutePath());
			img.setC(1);img.setDisplayRange(0, MRUtils.maxDisplayedBionanoM0);
			img.setC(2);img.setDisplayRange(0, MRUtils.maxDisplayedBionanoT1);
			img.setC(3);img.setDisplayRange(0, MRUtils.maxDisplayedBionanoT2);
			img.setC(4);img.setDisplayRange(0, 1);
			img.setC(5);img.setDisplayRange(-1, 1);
			for(int c=6;c<=img.getNChannels();c++) {img.setC(c);img.setDisplayRange(0, MRUtils.maxDisplayedBionanoM0);}
			
			String fileOut=new File(dirOut,infoTab[0]).getAbsolutePath();
			fileOut=new File(fileOut,strNew).getAbsolutePath();
			System.out.println("Sauvegarde in "+fileOut);
			if(makeShort)img=new Duplicator().run(img,1,3,1,4,1,1);
			IJ.saveAsTiff(img,fileOut);			
		}
		System.exit(0);
	}
	
	
	public static void testKhi2() {
		int[][]vals=new int[][] {
			{5,10 },
			{6,12 },
			{10,20 },
			{10,10 },
			{20,10 },
			{50,10 },
			{100,10 },
			{500,10 }
		};
		for(int i=0;i<vals.length;i++) {		System.out.println(" Khi="+vals[i][0] + " Nprms="+vals[i][1]+"  pval="+MRUtils.getPvalue(vals[i][0],vals[i][1]) );}
	}

	
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		testDims();
		System.exit(0);
		//testKhi2();
		//testNorm();
		chiasse();
		String path="/home/fernandr/Bureau/Test/Ghetto/fdfddfd/";
		ImagePlus imgT1=IJ.openImage(path+"T1short.tif");
		ImagePlus imgT2=IJ.openImage(path+"T2short.tif");
		ImagePlus imgMaps=IJ.openImage(path+"Maps.tif");
		int[]tr=new int[] {600,1200,2400};
		for(int t=1;t<=imgT1.getNFrames();t++) {
			for(int z=1;z<=imgT1.getNSlices();z++) {
				System.out.println("t="+t+" , z="+z);
				for(int c=1;c<=imgMaps.getNChannels();c++)imgMaps.getStack().setSliceLabel("_SIGMARICE=130_MAPS_TR="+tr[c-1]+"_TE=11", VitimageUtils.getCorrespondingSliceInHyperImage(imgMaps, c-1, z-1, t-1));
				for(int c=1;c<=imgT1.getNChannels();c++)imgT1.getStack().setSliceLabel("_SIGMARICE=130_T1SEQ_TR="+tr[c-1]+"_TE=11", VitimageUtils.getCorrespondingSliceInHyperImage(imgT1, c-1, z-1, t-1));
				for(int c=1;c<=imgT2.getNChannels();c++)imgT2.getStack().setSliceLabel("_SIGMARICE=130_T2SEQ_TR=10000_TE="+(11*c), VitimageUtils.getCorrespondingSliceInHyperImage(imgT2, c-1, z-1, t-1));
			}
		}
		ImagePlus []tabT1=stacksFromHyperstackFastBis(imgT1);
		ImagePlus []tabT2=stacksFromHyperstackFastBis(imgT2);
		ImagePlus []tabMaps=stacksFromHyperstackFastBis(imgMaps);
		ImagePlus []imgTab=new ImagePlus[3*5+3*5+16*5];
		for(int i=0;i<3*5;i++)imgTab[i]=tabMaps[i];
		for(int i=0;i<3*5;i++)imgTab[3*5+i]=tabT1[i];
		for(int i=0;i<16*5;i++)imgTab[3*5+3*5+i]=tabT2[i];
		ImagePlus res=Concatenator.run(imgTab);
		VitimageUtils.printImageResume(res);
		System.out.println("22 , " +imgT1.getNSlices()+" , "+imgT1.getNFrames());
		ImagePlus res2=HyperStackConverter.toHyperStack(res, 22,imgT1.getNSlices(),imgT1.getNFrames(),"xyztc","Fire");
		IJ.run(res2,"32-bit","");
		res2.show();
	}
	public static ImagePlus[]stacksFromHyperstackFastBis(ImagePlus hyper){
		int nbZ=hyper.getNSlices();
		int nbT=hyper.getNFrames();
		int nbC=hyper.getNChannels();
		int nb=nbT*nbC;
		ImagePlus []ret=new ImagePlus[nb];
		for(int ic=0;ic<nbC;ic++) {
			for(int it=0;it<nbT;it++) {
				int i=ic*nbT+it;
				System.out.println(ic+"/"+nbC+" ,  "+it+"/"+nbT);
				ret[i] = new Duplicator().run(hyper, 1+ic, 1+ic, 1, nbZ, 1+it, 1+it);
				VitimageUtils.adjustImageCalibration(ret[i],hyper);
				IJ.run(ret[i],"Grays","");
			}
		}
		return ret;
	}

	public TestingMRUtilsPackage() {
		MRI_HyperCurvesExplorer explorer=new MRI_HyperCurvesExplorer();
	}
*/
}
