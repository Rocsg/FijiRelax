package com.vitimage.aplimtools;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.time.Duration;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.Period;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

import org.apache.commons.io.FileUtils;

import com.sun.source.tree.MemberReferenceTree;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;
import com.vitimage.mrutils.HyperMRIT1T2;
import com.vitimage.mrutils.MRDataType;
import com.vitimage.common.TransformUtils.VolumeComparator;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.FolderOpener;
import ij.plugin.frame.PlugInFrame;

public class Lookup_Vitimage_Data  extends PlugInFrame{

	private static final long serialVersionUID = 1L;
	private String specimens[];
	private String daysSpecimens[][];
	private MRDataType sequencesSpecimens[][];
	private Fid[][][] fidBySpecimen;//Specimen , day , sequence
	private LocalDateTime[] referenceStartingDay;	
	String[]dirsVitimage=new String[] {"s_20190902_01","s_20180306_01","s_20191106_01"};
	private String sourceDirForFids;
	private String outputDir;
	private Fid[]fids ;//fidsPath;fidsName;fidsDate;fidsTime;fidsVersion;fidsType;HR 3D TRxxx 
	private boolean[] hasTimeLapse;
	private boolean[] hasT1T2;
	private boolean[] hasT2HR;
	private boolean[] hasGE3D;
	private int[] phaseSpecimens;
	private boolean waitForMatlab=false;
	boolean removeDataBehind=true;
	public Lookup_Vitimage_Data() {
		super("");
	}

	public static void test() {
	}


	public boolean isGe3dSliceInBounds(String name) {
		int a=Integer.parseInt(name.split("_SL")[1].split(".dcm")[0].split(".tif")[0]);
		if(a<256)return false;
		if(a>=768)return false;
		return true;
	}
	

	public boolean isT1T2SliceInBounds(String name) {
		if(true)return true;
		int a=Integer.parseInt(name.split("_SL")[1].split(".dcm")[0].split(".tif")[0]);
		if(a<2)return false;
		if(a>38)return false;
		return true;
	}

	
	@SuppressWarnings("unchecked")
	public void run(String arg) {
		checkInitialPathConfiguration();
		detectAllFids();
		Arrays.sort(fids,new FidComparatorByDate());
		identifySpecimensSequencesAndDays();
		checkWeirdDayNaming(false);
		keepOnlyLastVersionOfEach(false);
		removeUnselectedFids();
		printResumeOfExperiences(true);
		identifyPossibleExperiments(true);
		createInfoForNewExperiences();
		if(!this.waitForMatlab) {
			detectEndOfPhase1ForProcessedFidsAndSwitchThemToPhase2();
			executePhase2();
			executePhase3();
			executePhase4();
		}
	}

	public static void main(String[]args) {
		Lookup_Vitimage_Data look= new Lookup_Vitimage_Data();
//		test();
		look.run("");
	}

	
	/** Tools for gathering computation done*/
	public void executePhase3() {
	}
	

	public void executePhase5() {
		System.out.println("\nSTARTING PHASE 5 : STACKING GE3D");
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			String spec=specimens[indSpec];
			if(!spec.equals("B9001"))continue;
			if(phaseSpecimens[indSpec]!=5)continue;
			if(! hasGE3D[indSpec]) {
				phaseSpecimens[indSpec]=6;
				writeConfigFile(indSpec);
				continue;
			}
			String pathStep=new File(outputDir,"2_Gather_data").getAbsolutePath();
			String pathSpec=new File(pathStep,spec).getAbsolutePath();
			for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
				String day=this.daysSpecimens[indSpec][indDay];
				String pathDay=new File(pathSpec,"GE3D_DATA/"+day).getAbsolutePath();
				String pathOut=new File(pathSpec,"GE3D_DATA/"+day+".tif").getAbsolutePath();
				ImagePlus img=FolderOpener.open(pathDay, "");
				int Z=img.getNSlices();
				IJ.saveAsTiff(img, pathOut);
				//try {FileUtils.deleteDirectory(new File(pathDay));} catch (IOException e) {e.printStackTrace();
			}	
			phaseSpecimens[indSpec]=6;
			writeConfigFile(indSpec);
		}
	}

	public void executePhase6() {
		System.out.println("\nSTARTING PHASE 6 : COMPUTING MAPS");
		int total=0;
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			String spec=specimens[indSpec];
			if(!spec.equals("B9001"))continue;
			if(phaseSpecimens[indSpec]!=6)continue;
			if(! hasGE3D[indSpec]){
				phaseSpecimens[indSpec]=7;
				writeConfigFile(indSpec);
				continue;
			}
			String pathStep=new File(outputDir,"2_Gather_data").getAbsolutePath();
			String pathSpec=new File(pathStep,spec).getAbsolutePath();
			for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
				String day=this.daysSpecimens[indSpec][indDay];
				String pathDay=new File(pathSpec,"GE3D_DATA/"+day).getAbsolutePath();
				String pathOut=new File(pathSpec,"GE3D_DATA/"+day+".tif").getAbsolutePath();
				ImagePlus img=FolderOpener.open(pathDay, "");
				int Z=img.getNSlices();
				IJ.saveAsTiff(img, pathOut);
				//try {FileUtils.deleteDirectory(new File(pathDay));} catch (IOException e) {e.printStackTrace();
			}	
			phaseSpecimens[indSpec]=7;
			writeConfigFile(indSpec);
		}
	}
	
	
	public void executePhase7() {
		int count=0;
		for(int phase : phaseSpecimens)if(phase==7)count++;
		if(count==0)return;
		if(VitiDialogs.getYesNoUI("Compute again maps ?", "All computation are done. Therefore you can run again maps computation for both specimens\n"+
				"This is a very long operation. Do you really want to compute again maps ?")) {
			for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
				String spec=specimens[indSpec];
				if(!spec.equals("B9001"))continue;
				if(phaseSpecimens[indSpec]!=7)continue;
				if(! hasT1T2[indSpec])continue;
				
				String pathStep=new File(outputDir,"3_T1T2Compute_maps").getAbsolutePath();
				String pathSpec=new File(pathStep,spec).getAbsolutePath();
				for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
					String day=this.daysSpecimens[indSpec][indDay];
					String pathDay=new File(pathSpec,"hypermap_"+day+".tif").getAbsolutePath();
					HyperMRIT1T2 hyp=new HyperMRIT1T2(IJ.openImage(pathDay));
					hyp.computeMapsAgain();
					IJ.saveAsTiff(hyp.getCopy(),pathDay+"second.tif");
				}
			}
		}
	}

	
	
	public void executePhase4() {
		System.out.println("\nSTARTING PHASE 4 : CONVERTING DATA");
		int total=0;
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			String spec=this.specimens[indSpec];			
			if(!spec.equals("B001"))continue;
			if(phaseSpecimens[indSpec]!=4)continue;
			if(hasT1T2[indSpec]) {
				String pathStep=new File(outputDir,"2_Gather_data").getAbsolutePath();
				String pathSpec=new File(pathStep,spec).getAbsolutePath();
				for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
					String day=this.daysSpecimens[indSpec][indDay];
					String pathDay=new File(pathSpec,"T1T2_DATA/"+day).getAbsolutePath();
					System.out.print("processing T1T2 spec "+spec+"day"+indDay);
					int a=switchImagesToUnsignedShortAndReduceSize(pathDay,false);
					total+=a;
					System.out.println(", images switched to unsigned short = "+a+" . "+VitimageUtils.dou(a*0.128)+" MB saved");
				}
			}
			if(hasGE3D[indSpec]) {
				String pathStep=new File(outputDir,"2_Gather_data").getAbsolutePath();
				String pathSpec=new File(pathStep,spec).getAbsolutePath();
				for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
					String day=this.daysSpecimens[indSpec][indDay];
					String pathDay=new File(pathSpec,"GE3D_DATA/"+day).getAbsolutePath();
					System.out.print("processing Ge3D spec "+spec+"day"+indDay);
					int a=switchImagesToUnsignedShortAndReduceSize(pathDay,true);
					total+=a;
					System.out.println(", images switched to unsigned short = "+a+" . "+VitimageUtils.dou(a*0.128)+" MB saved");
				}
			}
			phaseSpecimens[indSpec]=5;
			writeConfigFile(indSpec);
		}
		System.out.println(" Total = "+total+" . "+VitimageUtils.dou(total*0.000128)+" GB saved on server");
		System.out.println("\nFINISHING PHASE 4 ");
	}

	public int switchImagesToUnsignedShortAndReduceSize(String path,boolean ge3d) {
		if(! new File(path).isDirectory()) {
			/*if(ge3d && (!isGe3dSliceInBounds(path))) {
				new File(path).delete();
				return 0;
			}*/
			
			ImagePlus img=IJ.openImage(path);
			img=VitimageUtils.convertFloatToShortWithoutDynamicChanges(img);
			img=VitimageUtils.cropImageShort(img, 64, 64, 0, 512-64, 512-64, 1);			
			new File(path).delete();
			IJ.saveAsTiff(img,path);
			return 1;
		}
		else {
			String []sons=new File(path).list();
			int count=0;
			for(String son : sons)count+=switchImagesToUnsignedShortAndReduceSize(new File(path,son).getAbsolutePath(),ge3d);
			return count;
		}
	}
	
	public void executePhase2() {
		System.out.println("\nSTARTING PHASE 2 : GATHERING DATA");
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			if(phaseSpecimens[indSpec]!=2 && phaseSpecimens[indSpec]!=3) continue;
			System.out.println(" Phase2, gathering data. Processing "+specimens[indSpec]);

			String sSpec=this.specimens[indSpec];			
			String pathStep1=new File(outputDir,"1_Process_fids").getAbsolutePath();
			String pathStep2=new File(outputDir,"2_Gather_data").getAbsolutePath();
			String pathSpec=new File(pathStep1,sSpec).getAbsolutePath();
			String pathSpec2=new File(pathStep2,sSpec).getAbsolutePath();
			//Gather T1T2 sequences data (and remove previous traces)
			if(hasT1T2[indSpec]) {
				System.out.println("processing T1T2 data");
				for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
					String day=this.daysSpecimens[indSpec][indDay];
					String pathDay=new File(pathSpec,day).getAbsolutePath();
					String pathDay2=new File(pathSpec2,"T1T2_DATA/"+day).getAbsolutePath();
					new File(pathDay2).mkdirs();
					for(String seq : new String[] {"TR600","TR1200","TR2400","TR10000"}) {
						String pathSeq=new File(pathDay,seq).getAbsolutePath();
						System.out.print("Looking for dicom dir for "+pathSeq+" : ");
						String pathToSourceDcm=findT1T2DicomDir(pathSeq);
						System.out.println(" found "+pathToSourceDcm);
						String pathToTargetDcm=new File(pathDay2,seq).getAbsolutePath();
						new File(pathToTargetDcm).mkdirs();
						try {FileUtils.copyDirectory(new File(pathToSourceDcm),new File(pathToTargetDcm));} catch (IOException e) {e.printStackTrace();}						
					}
				}
			}
			//Gather Ge3D data (and remove previous traces)
			if(hasGE3D[indSpec]) {
				System.out.println("processing Ge3d data");
				for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
					String day=this.daysSpecimens[indSpec][indDay];
					String pathDay=new File(pathSpec,day).getAbsolutePath();
					String pathDay2=new File(pathSpec2,"GE3D_DATA/"+day).getAbsolutePath();
					new File(pathDay2).mkdirs();
					String seq="GE3D";
					String pathSeq=new File(pathDay,seq).getAbsolutePath();
					System.out.print("Looking for dicom dir for "+pathSeq+" : ");
					String pathToSourceDcm=findGE3DDicomDir(pathSeq);
					System.out.println(" found "+pathToSourceDcm);
					String pathToTargetDcm=new File(pathDay2).getAbsolutePath();
					new File(pathToTargetDcm).mkdirs();
					try {FileUtils.copyDirectory(new File(pathToSourceDcm),new File(pathToTargetDcm));} catch (IOException e) {e.printStackTrace();}
				}

			}			
			phaseSpecimens[indSpec]=3;
			writeConfigFile(indSpec);
			if(removeDataBehind)try {FileUtils.deleteDirectory(new File(pathSpec));} catch (IOException e) {e.printStackTrace();}
			phaseSpecimens[indSpec]=4;
			writeConfigFile(indSpec);
		}
		System.out.println("\nFINISHING PHASE 2 ");
	}
		
	public void detectEndOfPhase1ForProcessedFidsAndSwitchThemToPhase2() {
		System.out.println("\nSTARTING PHASE 1 : detect finished matlab computation");
		String out=runCommand("ps auxw | grep matlab | grep fernandez");
		int nProcess=out.split("\n").length-1;
		if(nProcess>0) {IJ.showMessage("There is still some matlab process running. Cannot ensure if all the jobs are finished. Processed fids will be gathered next time");return;}
		String pathStep1=new File(outputDir,"1_Process_fids").getAbsolutePath();
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			if(phaseSpecimens[indSpec]!=1) continue;
			boolean isFinished=true;
			String sSpec=this.specimens[indSpec];
			String pathSpec=new File(pathStep1,sSpec).getAbsolutePath();
			for(String sDay : new File(pathSpec).list()) {
				String pathDay=new File(pathSpec,sDay).getAbsolutePath();
				for(String sSeq : new File(pathDay).list()) {
					String pathSeq=new File(pathDay,sSeq).getAbsolutePath();
					String[]list=new File(pathSeq).list();
					int count=0;
					for(String s : list)if(s.substring(0,4).equals("proc"))count++;
					if(count<2)isFinished=false;
				}
			}
			if(isFinished) {
				phaseSpecimens[indSpec]=2;
				writeConfigFile(indSpec);
			}
			else {
				IJ.showMessage("It seems that a matlab script should have been run to process "+sSpec+"\nbut haven't. Removing configuration file. Will be done on next run");
				String pathConfig=new File(outputDir,"Overview_of_experiments").getAbsolutePath();
				new File(pathConfig,specimens[indSpec]+".cfg").delete();
			}
		}
		System.out.println("FINISHING PHASE 1");
	}
	
	public void writeConfigFile(int indSpec) {
		String pathConfig=new File(outputDir,"Overview_of_experiments").getAbsolutePath();
		String pathSpec=new File(pathConfig,specimens[indSpec]+".cfg").getAbsolutePath();
		writeConfigFile(pathSpec,phaseSpecimens[indSpec], daysSpecimens[indSpec], sequencesSpecimens[indSpec]);
	}
	
	
	/** Management of startup for a specimen */	
	public void createInfoForNewExperiences() {
		String script2D="";
		String script3D="\nsleep 300\n";
		int totT2=0,totT1=0,tot3D=0;
		System.out.println("\nSTARTING PHASE 0 : creatingInfoForNew");
		int count=0;
		phaseSpecimens=new int[this.specimens.length];
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			String spec=this.specimens[indSpec];
			if(indSpec<27 || indSpec>=60)continue;
			//if(!spec.equals("B099"))continue;
			String pathConfig=new File(outputDir,"Overview_of_experiments").getAbsolutePath();
			String pathSpec=new File(pathConfig,spec+".cfg").getAbsolutePath();
			boolean found=false;
			if(new File(pathSpec).exists()) {
				found=true;
				Object []objs=readConfigFile(pathSpec);//Read the file B031. Read the current phase waiting. Read the sequences and days, and compare with the one found.
				phaseSpecimens[indSpec]=(int)objs[0];
				String[]tempDays=(String[])objs[1];
				MRDataType[]tempSeqs=(MRDataType[])objs[2];
				boolean isIdentical=true;
				if(tempDays.length != daysSpecimens.length)isIdentical=false;
				else {
					for(int indD=0;indD<tempDays.length;indD++)if(! tempDays[indD].equals(daysSpecimens[indD]))isIdentical=false;
				}
				if(tempSeqs.length != sequencesSpecimens.length)isIdentical=false;
				else {
					for(int indS=0;indS<tempSeqs.length;indS++)if(! tempSeqs[indS].equals(sequencesSpecimens[indS]))isIdentical=false;
				}
				if(isIdentical==false) {
					if(VitiDialogs.getYesNoUI("Data changed","Retrieving fids data for specimen"+spec+" : data changed since last import\nDo you want to remove all computation for this specimen and start it again ?")) {;
						found=false;
						removeArborescenceOfSpecimen(indSpec);
					}
				}
			}
			
			if(found==false || phaseSpecimens[indSpec]<2) {
				IJ.showMessage(" New data detected : "+specimens[indSpec]);
				phaseSpecimens[indSpec]=1;
				Object[]scrs=buildArborescenceOfSpecimen(indSpec);
				count++;
				int ind3=indSpec%3;
				script2D+=(String)scrs[0];
				script3D+=(String)scrs[1];
				totT1+=(Integer)scrs[2];
				totT2+=(Integer)scrs[3];
				tot3D+=(Integer)scrs[4];
				writeConfigFile(indSpec);
			}
		}
		

		int expectedTime=40*totT2+4*totT1+700*tot3D+300;
		int totalTime=expectedTime;
		int exphours=expectedTime/3600;
		expectedTime=expectedTime-exphours*3600;
		int expmins=(expectedTime/60);
		int expsecs=expectedTime-expmins*60;
		int expMB=512*tot3D+20*totT1+320*totT2;
		double expGB=VitimageUtils.dou(expMB/1024.0);

		
		script2D="#! /usr/bin/bash\n\n"+"date\necho \"Starting processing of "+totT1+" T1 data,  "+totT2+" T2 data,  "+tot3D+" Ge3d data\""+
			"\necho \"Expected computation time "+totalTime+" s = "+exphours+" hours, "+expmins+" minutes and "+expsecs+" seconds\""+script2D;
		VitimageUtils.writeStringInFile(script2D+"\n"+script3D, new File(outputDir,"scriptMatlab.sh").getAbsolutePath());
		new File(outputDir,"scriptMatlab.sh").setExecutable(true);
		

		if(count>0) {
			this.waitForMatlab=true;
			IJ.showMessage("Computation of MRI data is ready to launch using Bionano matlab tools.\n.\n"+
			"The process will generate "+(2*expGB)+" GB during Matlab run,\nthat will be converted into "+(expGB)+" by the following phases of this plugin.\n.\n"+
		" In order to run all the matlab process easily, a script has been generated to help you.\n"+
		"What to do next ? \nAfter this plugin finish, open a terminal and type the following command : \n"+"~/Vitimage/scriptMatlab.sh\n.\nAfter hours, come back here, and run again this ImageJ plugin."+
				"\n It will tell you if matlab computation is finished, and go on the computation.\nGood luck ! Be strong !");
		}
		else {
			new File(outputDir,"scriptMatlab.sh").delete();
		}
	}
	
	public void removeArborescenceOfSpecimen(int indSpec) {
		String pathStep1=new File(outputDir,"1_Process_fids").getAbsolutePath();
		String pathSpec=new File(pathStep1,specimens[indSpec]).getAbsolutePath();
		try {FileUtils.deleteDirectory(new File(pathSpec));	} catch (IOException e) {e.printStackTrace();}
		String pathStep2=new File(outputDir,"3_T1T2Compute_maps").getAbsolutePath();
		pathSpec=new File(pathStep2,specimens[indSpec]).getAbsolutePath();
		try {FileUtils.deleteDirectory(new File(pathSpec));	} catch (IOException e) {e.printStackTrace();}
		String pathStep3=new File(outputDir,"4_T1T2_Time_series").getAbsolutePath();
		pathSpec=new File(pathStep3,specimens[indSpec]).getAbsolutePath();
		try {FileUtils.deleteDirectory(new File(pathSpec));	} catch (IOException e) {e.printStackTrace();}
		String pathStep4=new File(outputDir,"5_Ge3d_Time_series").getAbsolutePath();
		pathSpec=new File(pathStep4,specimens[indSpec]).getAbsolutePath();
		try {FileUtils.deleteDirectory(new File(pathSpec));	} catch (IOException e) {e.printStackTrace();}
	}
				
	public Object[] buildArborescenceOfSpecimen(int indSpec) {
		String scriptRunImport2D="\n";
		String scriptRunImport3D="\n";
		int countT1=0;
		int countT2=0;
		int count3D=0;
		scriptRunImport2D+="#PROCESSING SPECIMEN "+specimens[indSpec]+"\n";
		scriptRunImport2D+="echo \"Starting specimen "+specimens[indSpec]+"\""+"\ndate\n";
		scriptRunImport3D+="#PROCESSING SPECIMEN "+specimens[indSpec]+"\n";
		scriptRunImport3D+="echo \"Starting specimen "+specimens[indSpec]+"\""+"\ndate\n";
		//Build arborescence in process fids, and sherpa files : processing.bnp, listSer.ser
		String pathStep1=new File(outputDir,"1_Process_fids").getAbsolutePath();
		String pathSpec=new File(pathStep1,specimens[indSpec]).getAbsolutePath();
		new File(pathSpec).mkdirs();
		for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
			String day=this.daysSpecimens[indSpec][indDay];
			String pathDay=new File(pathSpec,day).getAbsolutePath();
			new File(pathDay).mkdirs();
			scriptRunImport2D+="date\n";
			scriptRunImport3D+="date\n";
			for(int indSeq=0;indSeq<this.sequencesSpecimens[indSpec].length;indSeq++) {
				String seq=this.sequencesSpecimens[indSpec][indSeq].toString();
				if( (!seq.equals("TR600")) && (!seq.equals("TR1200")) && (!seq.equals("TR2400")) && (!seq.equals("TR10000")) && (!seq.equals("GE3D")) )continue;
				if( ( (seq.equals("TR600")) || (seq.equals("TR1200")) || (seq.equals("TR2400")) || (seq.equals("TR10000")) ) && !hasT1T2[indSpec])continue;
				if( ( (seq.equals("GE3D") ) && !hasGE3D[indSpec]))continue;
				if  (seq.equals("T2SEQHR") )continue;
				String pathSeq=new File(pathDay,seq).getAbsolutePath();
				new File(pathSeq).mkdirs();
				new File(pathSeq,"DATA").mkdirs();
				VitimageUtils.writeStringInFile(getTemplateProcessingBnp(),new File(pathSeq,"processing.bnp").getAbsolutePath());
				VitimageUtils.writeStringInFile(fidBySpecimen[indSpec][indDay][indSeq].path+"/fid", new File(pathSeq,"listSer.ser").getAbsolutePath());
				if( ( (seq.equals("TR600")) || (seq.equals("TR1200")) || (seq.equals("TR2400")) || (seq.equals("TR10000")) ) ) {
					scriptRunImport2D+="cd "+pathSeq+"\n";
					scriptRunImport2D+="echo \"Sending computation for "+specimens[indSpec]+" "+day+" "+seq+"\""+pathSeq+"\n";
					scriptRunImport2D+="nohup matlab -nodesktop -nodisplay -nojvm -r \"addpath('/users/bionanonmri/fernandez/Boulot/BioNanoCODE_2020_03_17'); run Processing_fid2dcm_2D_BIS.m ; quit\" 2>&1 logMatlab.log &\n";
					scriptRunImport2D+="sleep "+(seq.equals("TR10000") ? 40 : 4)+"\n";
					scriptRunImport2D+="\n";
					if(seq.equals("TR10000"))countT2++;
					else countT1++;
				}
				if( seq.equals("GE3D") ) {
					scriptRunImport3D+="cd "+pathSeq+"\n";
					scriptRunImport3D+="echo \"Sending computation for "+specimens[indSpec]+" "+day+" "+seq+"\""+pathSeq+"\n";
					scriptRunImport3D+="nohup matlab -nodesktop -nodisplay -nojvm -r \"addpath('/users/bionanonmri/fernandez/Boulot/BioNanoCODE_2020_03_17'); run Processing_fid2dcm_3D_AXIAL_BIS.m ; quit\" 2>&1 logMatlab.log &\n";
					scriptRunImport3D+="sleep 350\n";
					scriptRunImport3D+="\n";
					count3D++;
				}
			}
		}
			
		//Build arborescence in compute maps, time lapse series and Ge3D
		if(hasT1T2[indSpec]) {
			String pathStep0=new File(outputDir,"2_Gather_data").getAbsolutePath();
			pathSpec=new File(pathStep0,specimens[indSpec]+"/T1T2_DATA").getAbsolutePath();
			new File(pathSpec).mkdirs();

			String pathStep2=new File(outputDir,"3_T1T2Compute_maps").getAbsolutePath();
			pathSpec=new File(pathStep2,specimens[indSpec]).getAbsolutePath();
			new File(pathSpec).mkdirs();
			if(hasTimeLapse[indSpec]) {
				String pathStep3=new File(outputDir,"4_T1T2_Time_series").getAbsolutePath();
				pathSpec=new File(pathStep3,specimens[indSpec]).getAbsolutePath();
				new File(pathSpec).mkdirs();
			}
		}
		if(hasGE3D[indSpec]) {
			String pathStep0=new File(outputDir,"2_Gather_data").getAbsolutePath();
			pathSpec=new File(pathStep0,specimens[indSpec]+"/GE3D_DATA").getAbsolutePath();
			new File(pathSpec).mkdirs();

			String pathStep4=new File(outputDir,"5_Ge3d_Time_series").getAbsolutePath();
			pathSpec=new File(pathStep4,specimens[indSpec]).getAbsolutePath();
			new File(pathSpec).mkdirs();
		}
		return new Object[] {scriptRunImport2D,scriptRunImport3D,countT1,countT2,count3D};
	}
			
	public Object[]readConfigFile(String file){
		String s=VitimageUtils.readStringFromFile(file);
		int phase=Integer.parseInt(s.split("\n")[1].split("=")[1]);
		int ndays =Integer.parseInt(s.split("\n")[2].split("=")[1]);
		String days =(s.split("\n")[3]);
		String[]day=new String[ndays];
		for(int i=0;i<ndays;i++)day[i]=days.split(" ")[i];
		int nseqs =Integer.parseInt(s.split("\n")[4].split("=")[1]);
		String seqs =(s.split("\n")[5]);
		MRDataType[]seq=new MRDataType[nseqs];
		for(int i=0;i<nseqs;i++)seq[i]=getMRDataType(seqs.split(" ")[i]);
		return new Object[] {phase,day,seq};
	}
	
	public void writeConfigFile(String file,int phase,String[]days,MRDataType[]sequences){
		String s=LocalDateTime.now().toString()+"\n";
		s+="PHASE="+phase+"\n";
		s+="N_DAYS="+days.length+"\n";
		for(int i=0;i<days.length-1;i++)s+=days[i]+" ";
		s+=days[days.length-1]+"\n";
		s+="N_SEQS="+sequences.length+"\n";
		for(int i=0;i<sequences.length-1;i++)s+=sequences[i]+" ";
		s+=sequences[sequences.length-1]+"\n";
		VitimageUtils.writeStringInFile(s,file);
	}
	

	/** Utilities for specimen experiment detection*/
	public int[]getIndexes(int indSpec){
		MRDataType[]typesLookup=new MRDataType[] {
				MRDataType.TR600,/* 0 */
				MRDataType.TR1200,/* 1 */
				MRDataType.TR2400,/* 2 */
				MRDataType.TR10000,/* 3 */
				MRDataType.T2SEQHR,/* 4 */
				MRDataType.GE3D     /* 5 */
				};
		int []indices=new int[typesLookup.length];
		for(int i=0;i<indices.length;i++) {
			indices[i]=-1;
			for(int indSeq=0;indSeq<this.sequencesSpecimens[indSpec].length;indSeq++)if(this.sequencesSpecimens[indSpec][indSeq].equals(typesLookup[i]))indices[i]=indSeq;
		}
		return indices;
	}
		
	public void identifyPossibleExperiments(boolean verbose){
		hasTimeLapse=new boolean[this.specimens.length];
		hasT1T2=new boolean[this.specimens.length];
		hasT2HR=new boolean[this.specimens.length];
		hasGE3D=new boolean[this.specimens.length];
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			int[]indices=getIndexes(indSpec);
			boolean hasT1T2Timelapse=true;
			boolean hasT2HRTimelapse=true;
			boolean hasGe3dTimelapse=true;
			for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
				//Look for a T1T2 experiment (TR600 and TR1200 and TR2400 and TR10000 present in both days)
				if(indices[0]<0 || fidBySpecimen[indSpec][indDay][indices[0]]==null)hasT1T2Timelapse=false;
				if(indices[1]<0 || fidBySpecimen[indSpec][indDay][indices[1]]==null)hasT1T2Timelapse=false;
				if(indices[2]<0 || fidBySpecimen[indSpec][indDay][indices[2]]==null)hasT1T2Timelapse=false;
				if(indices[3]<0 || fidBySpecimen[indSpec][indDay][indices[3]]==null)hasT1T2Timelapse=false;

				//Look for a T2HR experiment (TR10000HR present in both days)
				if(indices[4]<0 || fidBySpecimen[indSpec][indDay][indices[4]]==null)hasT2HRTimelapse=false;

				//Look for a Ge3D experiment (TR600 and TR1200 and TR2400 and TR10000 present in both days)
				if(indices[5]<0 || fidBySpecimen[indSpec][indDay][indices[5]]==null)hasGe3dTimelapse=false;
			}
			hasTimeLapse[indSpec]=this.daysSpecimens[indSpec].length>1;
			hasT1T2[indSpec]=hasT1T2Timelapse;  
			hasT2HR[indSpec]=hasT2HRTimelapse;
			hasGE3D[indSpec]=hasGe3dTimelapse;	
			String str=" over "+daysSpecimens[indSpec].length+" days : ";
			for(int i=0;i<daysSpecimens[indSpec].length;i++)str+=" "+daysSpecimens[indSpec][i];
			if(verbose)System.out.println(this.specimens[indSpec]+" has "+(hasT1T2[indSpec] ? " T1T2  ": "")+(hasGE3D[indSpec] ? " GE3D" : "")+(hasT2HR[indSpec] ? " T2HR  ":"")+str);
		}
	}
	
	public void printResumeOfExperiences(boolean isShort) {
		if(!isShort) {
			System.out.println("\nEXHAUSTIVE RESUME OF IDENTIFIED EXPERIMENTS\n-----------------------------\n\n");
			for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
				String spec=this.specimens[indSpec];
				System.out.println("\nSPECIMEN "+spec);
				for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
					String day=this.daysSpecimens[indSpec][indDay];
					System.out.print(" "+(day+(day.length()==2 ? "  " : day.length()==3 ? " " : "")));
					for(int indSeq=0;indSeq<this.sequencesSpecimens[indSpec].length;indSeq++) {
						MRDataType mrd=this.sequencesSpecimens[indSpec][indSeq];
						int incr=0;
						for(Fid f:fids)if(f.day.equals(day) && f.typeAcq.equals(mrd) && f.specimen.equals(spec))incr++;
						if(incr>0)System.out.print("  "+mrd+(incr>1 ?"("+incr+")" : ""));
					}
					System.out.println();
				}		
			}
		}
		else {
			System.out.println("\nSHORTLIST RESUME OF IDENTIFIED EXPERIMENTS\n-----------------------------\n\n");
			for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
				String spec=this.specimens[indSpec];
				System.out.println("\nSPECIMEN "+spec);
				fidBySpecimen[indSpec]=new Fid[this.daysSpecimens[indSpec].length][];
				for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
					String day=this.daysSpecimens[indSpec][indDay];
					System.out.print(" "+(day+(day.length()==2 ? "  " : day.length()==3 ? " " : "")));
					fidBySpecimen[indSpec][indDay]=new Fid[this.sequencesSpecimens[indSpec].length];
					for(int indSeq=0;indSeq<this.sequencesSpecimens[indSpec].length;indSeq++) {
						MRDataType mrd=this.sequencesSpecimens[indSpec][indSeq];
						for(Fid f:fids)if(f.day.equals(day) && f.typeAcq.equals(mrd) && f.specimen.equals(spec)) {
							if(fidBySpecimen[indSpec][indDay][indSeq]==null) {fidBySpecimen[indSpec][indDay][indSeq]=f;}
							else if(f.dateAcq.isAfter(fidBySpecimen[indSpec][indDay][indSeq].dateAcq))fidBySpecimen[indSpec][indDay][indSeq]=f;
						}
						if(fidBySpecimen[indSpec][indDay][indSeq]!=null)System.out.print("  "+mrd);

					}
					System.out.println();
				}
			}
		}
	}
	

	
	/** Utilities for dicom dir detection*/
	public String findT1T2DicomDir(String pathSource) {
		String targetDir="";
		String curDir=pathSource;
		long time=0;

		//In this path, there is a target dir named "DATA", other files, and another dir. The last one of these dirs contains a heavy complexify arborescence. 
		String []dirs=new File(curDir).list();
		for(String dir : dirs) {
			if(dir.substring(0,4).equals("proc") && new File(curDir,dir).lastModified()>time) {
				targetDir=dir;
				time=new File(curDir,dir).lastModified();
			}
		}
		curDir=new File(curDir,targetDir).getAbsolutePath();

		//At the very leaves of this arborescence, there is images. We want to collect them and set them into the DATA dir
		
		//Check directory s_something
		dirs=new File(curDir).list();
		for(String dir : dirs) {
			if(dir.substring(0, 2).equals("s_")) {
				targetDir=dir;
			}
		}		
		curDir=new File(curDir,targetDir).getAbsolutePath();

		//Check acq name
		dirs=new File(curDir).list();
		targetDir=dirs[0];
		curDir=new File(curDir,targetDir).getAbsolutePath();

		//Check SL_RAW
		targetDir="SL_RAW";
		curDir=new File(curDir,targetDir).getAbsolutePath();
		
		//Copy recursively the dir here "TR000600" to target
		dirs=new File(curDir).list();
		targetDir=dirs[0];
		curDir=new File(curDir,targetDir).getAbsolutePath();
		return curDir;		
	}
	
	public String findGE3DDicomDir(String pathSource) {
		String targetDir="";
		String curDir=pathSource;
		long time=0;

		//In this path, there is a target dir named "DATA", other files, and another dir. The last one of these dirs contains a heavy complexify arborescence. 
		String []dirs=new File(curDir).list();
		for(String dir : dirs) {
			if(dir.substring(0,4).equals("proc") && new File(curDir,dir).lastModified()>time) {
				targetDir=dir;
				time=new File(curDir,dir).lastModified();
			}
		}
		curDir=new File(curDir,targetDir).getAbsolutePath();

		//At the very leaves of this arborescence, there is images. We want to collect them and set them into the DATA dir
		
		//Check directory s_something
		dirs=new File(curDir).list();
		for(String dir : dirs) {
			if(dir.substring(0, 2).equals("s_")) {
				targetDir=dir;
			}
		}		
		curDir=new File(curDir,targetDir).getAbsolutePath();

		//Check acq name
		dirs=new File(curDir).list();
		targetDir=dirs[0];
		curDir=new File(curDir,targetDir).getAbsolutePath();

		//Check SL_RAW
		targetDir="SL_RAW";
		curDir=new File(curDir,targetDir).getAbsolutePath();
		
		//Check TRSOMETHING
		targetDir=new File(curDir).list()[0];
		curDir=new File(curDir,targetDir).getAbsolutePath();

		//Check TESOMETHING
		targetDir=new File(curDir).list()[0];
		curDir=new File(curDir,targetDir).getAbsolutePath();

		return curDir;		
	}

	
	
	
	/** Utilities for fid detection and listing*/
	public void removeUnselectedFids() {
		int count=0;
		for(Fid f: fids)if(f.isSelected)count++;
		Fid[]newFids=new Fid[count];
		count=0;
		for(int i=0;i<fids.length;i++) if(fids[i].isSelected)newFids[count++]=fids[i];
		fids=newFids;
	}
	
	public void keepOnlyLastVersionOfEach(boolean verbose) {
		fidBySpecimen=new Fid[this.specimens.length][][];
		for(int indSpec=0;indSpec<this.specimens.length;indSpec++) {
			String spec=this.specimens[indSpec];
			fidBySpecimen[indSpec]=new Fid[this.daysSpecimens[indSpec].length][];
			for(int indDay=0;indDay<this.daysSpecimens[indSpec].length;indDay++) {
				String day=this.daysSpecimens[indSpec][indDay];
				fidBySpecimen[indSpec][indDay]=new Fid[this.sequencesSpecimens[indSpec].length];
				for(int indSeq=0;indSeq<this.sequencesSpecimens[indSpec].length;indSeq++) {
					MRDataType mrd=this.sequencesSpecimens[indSpec][indSeq];
					for(Fid f:fids)if(f.day.equals(day) && f.typeAcq.equals(mrd) && f.specimen.equals(spec)) {
						if(fidBySpecimen[indSpec][indDay][indSeq]==null) {fidBySpecimen[indSpec][indDay][indSeq]=f;f.isSelected=true;}
						else if(f.dateAcq.isAfter(fidBySpecimen[indSpec][indDay][indSeq].dateAcq)) {
							fidBySpecimen[indSpec][indDay][indSeq].isSelected=false;
							f.isSelected=true;
							fidBySpecimen[indSpec][indDay][indSeq]=f;
						}
					}
				}
			}
		}
	}
			
	public void checkWeirdDayNaming(boolean verbose) {
		for(Fid f : fids) {
			int dayStr=Integer.parseInt(f.day.replace("J",""));
			double distance=VitimageUtils.dou(Math.abs(dayStr-f.daysAfterStart));
			if(verbose)System.out.println("Distance="+distance+"   "+f.name);
			if(distance>10) {
				if(VitiDialogs.getYesNoUI("Warning : Weird day of experience","Warning : declared day in experience ("+f.day+") is not going well\nwith the actual day from J0\n.\nDo you want to change this ?")){
					String s=VitiDialogs.getStringUI("Suggested day : J"+(int)Math.round(distance), "Tell the day", "J"+(int)Math.round(distance), false);
					f.day=s;
				}
			}
		}
	}
	
	public void identifySpecimensSequencesAndDays() {
		ArrayList<String>names=new ArrayList<String>();
		ArrayList<Integer>numbers=new ArrayList<Integer>();
		for(Fid f:fids) {
			boolean found=false;
			for(int n=0;n<names.size();n++) {
				if(f.specimen.equals(names.get(n)) ){found=true;numbers.set(n,numbers.get(n)+1);}
			}
			if(!found) {
				names.add(f.specimen);
				numbers.add(1);
			}
		}
		for(int i=0;i<names.size();i++) {
			System.out.println("Specimen identified "+names.get(i)+" has "+ numbers.get(i));
		}
		this.specimens=names.toArray(new String[names.size()]);
		this.daysSpecimens=new String[this.specimens.length][];
		this.sequencesSpecimens=new MRDataType[this.specimens.length][];
		this.referenceStartingDay=new LocalDateTime[names.size()];
		
		System.out.println("\n-------------------------------\n");
		for(int indSpec=0;indSpec<this.specimens.length ;indSpec++) {
			names.clear();
			numbers.clear();
			String spec=this.specimens[indSpec];
			System.out.println("\nIdentification jours du specimen "+spec);
			for(Fid f:fids) {
				if(!f.specimen.equals(spec))continue;
				if(this.referenceStartingDay[indSpec]==null) {
					int delta=Integer.parseInt(f.day.replace("J", ""));
					this.referenceStartingDay[indSpec]=f.dateAcq.minusDays(delta);
					f.daysAfterStart=delta;}
				else {
					f.daysAfterStart=Duration.between(this.referenceStartingDay[indSpec],f.dateAcq).getSeconds()/86400.0;
				}
				boolean found=false;
				for(int n=0;n<names.size();n++) {
					if(f.day.equals(names.get(n)) ){found=true;numbers.set(n,numbers.get(n)+1);}
				}
				if(!found) {
					names.add(f.day);
					numbers.add(1);
				}
			}
			for(int i=0;i<names.size();i++) {
				System.out.println(names.get(i)+" : "+ numbers.get(i)+"  |  ");
			}
			this.daysSpecimens[indSpec]=names.toArray(new String[names.size()]);
		}		

		System.out.println("\n-------------------------------\n");
		for(int indSpec=0;indSpec<this.specimens.length ;indSpec++) {
			ArrayList<MRDataType>namesMR=new ArrayList<MRDataType>();
			namesMR.clear();
			numbers.clear();
			String spec=this.specimens[indSpec];
			System.out.println("\nIdentification sequences du specimen "+spec);
			for(Fid f:fids) {
				if(!f.specimen.equals(spec))continue;
				boolean found=false;
				for(int n=0;n<namesMR.size();n++) {
					if(f.typeAcq.equals(namesMR.get(n)) ){found=true;numbers.set(n,numbers.get(n)+1);}
				}
				if(!found) {
					namesMR.add(f.typeAcq);
					numbers.add(1);
				}
			}
			for(int i=0;i<namesMR.size();i++) {
				System.out.println(namesMR.get(i)+" : "+ numbers.get(i)+"  |  ");
			}
			this.sequencesSpecimens[indSpec]=namesMR.toArray(new MRDataType[namesMR.size()]);
		}		
	}

	public void identifySpecimens() {
		ArrayList<String>names=new ArrayList<String>();
		ArrayList<Integer>numbers=new ArrayList<Integer>();
		for(Fid f:fids) {
			boolean found=false;
			for(int n=0;n<names.size();n++) {
				if(f.specimen.equals(names.get(n)) ){found=true;numbers.set(n,numbers.get(n)+1);}
			}
			if(!found) {
				names.add(f.specimen);
				numbers.add(1);
			}
		}
		for(int i=0;i<names.size();i++) {
			System.out.println("Specimen identified "+names.get(i)+" has "+ numbers.get(i));
		}
		this.specimens=names.toArray(new String[names.size()]);
	}

	public void detectDoublons() {
		for(int i=0;i<fids.length;i++) {
			for(int j=i+1;j<fids.length;j++) {
				if(fids[i].equals(fids[j])) {
					System.out.println("\n\n\n\nACHTUNG : EQUALITY !!");
					System.out.println("\nFIRST--------- i="+i+" : ");
					System.out.println(fids[i].toFullString());
					System.out.println("NEIS= : ");
					for(int k=Math.max(0, i-2);k<Math.min(fids.length,i+3);k++)System.out.println("i_nei="+k+" : "+fids[k].toFullString());
					System.out.println("\nSECOND--------- j="+j+" : ");
					System.out.println(fids[j].toFullString());
					System.out.println("NEIS= : ");
					for(int k=Math.max(0, j-2);k<Math.min(fids.length,j+3);k++)System.out.println("j_nei="+k+" : "+fids[k].toFullString());
				}
			}
		}
	}
	
	public void checkInitialPathConfiguration() {	
		sourceDirForFids="/users/bionanonmri/goze/NMRI/montpellier94t/studies/";
		outputDir=new File(System.getProperty("user.home"),"Vitimage").getAbsolutePath();
		if(! new File("/users/bionanonmri/").exists()) {
			if(! new File("/home/fernandr/").exists()) {IJ.showMessage("Error : this plugin should be run on the bionano server foundry.univ-montp2.fr");return;}
			else { 
				sourceDirForFids="/home/fernandr/Bureau/A_Test/TestFids";
				outputDir="/home/fernandr/Bureau/A_Test/output";
			}
		}
		prepareDirs(outputDir);
	}
	
	public void detectAllFids() {
		ArrayList<String> fidList = new ArrayList<String>();
		for(String s : dirsVitimage) {
			ArrayList<String>tempList=lookupFidDirFromSourcePath(sourceDirForFids+s);
			fidList.addAll(tempList);
		}
		fidList=selectStringFromFids(fidList,"_B");	
		fidList=excludeStringFromFids(fidList,"VB");	
		fidList=excludeStringFromFids(fidList,"BM");	
		fidList=excludeStringFromFids(fidList,"Blank");	
		fidList=excludeStringFromFids(fidList,"Boutures");	
		fidList=excludeStringFromFids(fidList,"diff");	
		fidList=excludeStringFromFids(fidList,"flux");	
		String[]fidArray=VitimageUtils.stringArraySort(fidList.toArray(new String[fidList.size()]));
		fids=new Fid[fidArray.length];
		for(int i=0;i<fidArray.length;i++)fids[i]=new Fid(fidArray[i]);
	}

	public void prepareDirs(String outputDir) {
		if(!new File(outputDir).exists()) {
			new File(outputDir).mkdir();
			new File(outputDir,"1_Process_fids").mkdir();
			new File(outputDir,"2_Gather_data").mkdir();
			new File(outputDir,"3_T1T2Compute_maps").mkdir();
			new File(outputDir,"4_T1T2_Time_series").mkdir();
			new File(outputDir,"5_Ge3d_Time_series").mkdir();
			new File(outputDir,"Config").mkdir();
			new File(outputDir,"Overview_of_experiments").mkdir();
			new File(outputDir,"Config/Logs").mkdir();
			new File(outputDir,"Temp").mkdir();
			new File(outputDir,"Temp/Logs").mkdir();
			new File(outputDir,"Temp/Txts").mkdir();
			new File(outputDir,"Temp/Transformations").mkdir();
			new File(outputDir,"Temp/Imgs").mkdir();
		}
	}
	
	public ArrayList<String> lookupFidDirFromSourcePath(String source){
		ArrayList<String> ar=new ArrayList<String>();
		String []sons=new File(source).list();
		String pathToSon;
		for(String so : sons) {
			pathToSon=new File(source,so).getAbsolutePath();
			if(so.endsWith(".fid")) ar.add(pathToSon);
			else if(new File(source,so).isDirectory()) ar.addAll(lookupFidDirFromSourcePath(pathToSon));
		}
		return ar;
	}
	
	public ArrayList<String> selectStringFromFids(ArrayList<String>in,String str){
		System.out.print("Starting selection of fid with "+str+" in their names... ");
		ArrayList<String>out=new ArrayList<String>();
		for(String s : in) if (s.contains(str))out.add(s);
		System.out.println(", number at start : "+in.size()+" at end : "+out.size());
		return out;
	}
	
	public ArrayList<String> excludeStringFromFids(ArrayList<String>in,String str){
		System.out.print("Starting selection of fid without "+str+" in their names... ");
		ArrayList<String>out=new ArrayList<String>();
		for(String s : in) if (! s.contains(str))out.add(s);
		System.out.println("At start : "+in.size()+" at end : "+out.size());
		return out;
	}

	public void printFids(String s,boolean full){
		System.out.println("\n\nList of fids "+s+"\n---------------------------------\n\n");
		for(int i=0;i<fids.length;i++) {
			System.out.println(i+(i<10 ? "   " : (i<100) ? "  " : i<1000 ? " " : "")+" : "+(full ? fids[i].toFullString() :fids[i]));
		}
	}

	
	/** General purpose utilities*/
	class FidComparatorByDate implements java.util.Comparator {
	   public int compare(Object o1, Object o2) {
	      return ((Fid)o1).dateAcq.compareTo(((Fid)o2).dateAcq);
	   }
	}

	public String getTemplateProcessingBnp() {
		String s="% PROCESSING PARAMETERS FROM DEFAULT VALUES"+
	"\n% BEGIN"+
	"\n% Create RAW DICOM  "+
	"\nyes"+
	"\n% Create CLEAN DICOM  "+
	"\nno"+
	"\n% Create MASK DICOM  "+
	"\nno"+
	"\n% apodization in read direction  <apo_p>   "+
	"\n1"+
	"\n% apodization in phase direction <apo_v>   "+
	"\n1"+
	"\n% apodization in phase 2 direction (available only for 3D) <apo_v2>   "+
	"\n1"+
	"\n% Zero Filling Z? OR final size  in read direction <final_p> %  "+
	"\nZ1"+
	"\n% Zero Filling Z? OR final size in phase direction <final_v> %  "+
	"\nZ1"+
	"\n% Zero Filling Z? OR final size in phase direction (available only for 3D) <final_v2> % "+ 
	"\nZ1"+
	"\n% noise cut off <noise_factor> Masked Images with values below Noise*noise_factor are set at zero %  "+
	"\n3"+
	"\n% cleaning <cleaning_factor> integer  "+
	"\n3"+
	"\n% Cleaning AnimalID  / ident ic1 ic2 %  "+
	"\n1"+
	"\n66"+
	"\n% Find Markers in AnimalID  : Marker and Number of Marker % -1 to force name  %"+
	"\n_"+
	"\n3"+
	"\n% Smooth for Maps <smoothmap> %  "+
	"\n0.0001"+
	"\n% Image Registration "+
	"\n% Number Maximum of Reference Tube"+
	"\n1"+
	"\n% Radius Maximum of Reference Tube in [micron]"+
	"\n600"+
	"\n% END  ";
		return s;
	}
	
	public String runCommand(String command) {
		ProcessBuilder processBuilder = new ProcessBuilder();
		processBuilder.command("bash","-c", command);
		String ret="";
		try {
			Process process = processBuilder.start();
			StringBuilder output = new StringBuilder();
			BufferedReader reader = new BufferedReader(	new InputStreamReader(process.getInputStream()));
			String line;
			while ((line = reader.readLine()) != null) {output.append(line + "\n");	}
			int exitVal = process.waitFor();
			ret=output.toString();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return ret;
	}
	
	public MRDataType getMRDataType(String s) {
		switch(s) {
		case 	"TR600" : return MRDataType.TR600;
		case 	"TR1200" : return MRDataType.TR1200;
		case 	"TR2400" : return MRDataType.TR2400;
		case 	"TR10000" : return MRDataType.TR10000;
		case 	"T2SEQHR" : return MRDataType.T2SEQHR;
		case 	"GE3D" : return MRDataType.GE3D;
		default : return MRDataType.OTHER;
		}
	}
	

	
}
