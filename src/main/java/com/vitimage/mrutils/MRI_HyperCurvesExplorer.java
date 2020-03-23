package com.vitimage.mrutils;


import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.TextField;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.JPanel;

import com.vitimage.common.TransformUtils;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotCanvas;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.Opener;
import ij.plugin.Duplicator;
import ij.plugin.frame.PlugInFrame;
import ij.plugin.frame.RoiManager;

/**
 * 
 * @author Rocsg
 * MRI_HyperCurvesExplorer is a user-friendly tool for exploration of T1 T2 relaxation curves coming from T1 and T2 sequence
 * 
 * Acknowledgements :
 * Maida Cardoso and Christophe Goze for data and explanations about T1 and T2 times
 
 * TODO : 
 ******** Critical bugs ********* 
 * Read T1 values from the sequence
 * Open from a gathered hyperimage
 * 
 ******** Elements to investigate and possible refactoring ********* 
 *(priority=0/3, difficulty=1/3) None
 * 
 ******** Fixes testing ********* 
 *(None)
 * 
 ******** Testing needs ********* 
 * 
 ******** User requests ********* 
 * (priority=2/3, difficulty=1/3) None
 */


public class MRI_HyperCurvesExplorer extends PlugInFrame implements ActionListener,KeyListener, MouseListener {// 

	/** Entry points and startup functions*/
	public static void main(String[]args) {
		runExplorer("/home/fernandr/Bureau/Traitements/Sorgho/Test SSM1/Out_input_Fijiyama/Output/Exported_data/Data_combined.tif");
	}

	private ImagePlus M0map;
	private ImagePlus T1map;
	private ImagePlus T2map;
	private Object ij;

	public MRI_HyperCurvesExplorer() {
		super("Vitimage MRI Water tracker ");
	}
	
	public void run(String arg) {
		System.out.println("RUN CALLED !");
		runExplorer(null);
	}
	
	public static void runExplorer(String imgPath) {
		ImagePlus fullHyp=null;
		if (imgPath!=null)fullHyp=IJ.openImage(imgPath);
		else fullHyp=VitiDialogs.chooseOneImageUI("Open hyperimage", "Open a hyperimage built using T1T2MapImporter");
		
		IJ.log("Opening hyperimage ds s");
		System.out.println("Opening hyperimageS");
		VitimageUtils.printImageResume(fullHyp,"Hyperimage T1 T2 and maps");
		
		MRI_HyperCurvesExplorer mrExplo=new MRI_HyperCurvesExplorer();
		mrExplo.fullHyperImage=fullHyp;
		if(mrExplo.ij==null)mrExplo.ij = new ImageJ();
		else mrExplo.ij = IJ.getInstance();
		IJ.log("Starting MRI Curve Explorer");
		mrExplo.extractMapsAndSequences();
		mrExplo.setupStructures();
		mrExplo.startGui();
	}
	
	public void extractMapsAndSequences() {
		//Extract dimensions, number of time-points and number echoes for each sequence
		dims=VitimageUtils.getDimensions(fullHyperImage);
		voxS=VitimageUtils.getVoxelSizes(fullHyperImage);
		int nbCanTot=fullHyperImage.getNChannels();
		this.numberZ=dims[2];
		this.numberTimes=fullHyperImage.getNFrames();
		this.actualDay=new int[this.numberTimes];

		//Look for maps, extract them, range them, and display them. If not, exit
		this.M0map=new Duplicator().run(fullHyperImage, 1, 1, 1, this.numberZ, 1,this.numberTimes);
		this.T1map=new Duplicator().run(fullHyperImage, 2, 2, 1, this.numberZ, 1,this.numberTimes);
		this.T2map=new Duplicator().run(fullHyperImage, 3, 3, 1, this.numberZ, 1,this.numberTimes);
		this.M0map.setDisplayRange(0, standardCapillaryM0);
		this.T1map.setDisplayRange(minAcceptableT1, maxDisplayedT1);
		this.T2map.setDisplayRange(minAcceptableT2, maxDisplayedT2);

		//Look for sequences
		for(int can=3;can<nbCanTot;can++) {
			if(MRUtils.getDataTypeOfThisMagneticResonanceSlice(fullHyperImage,can,0,0)== MRDataType.T1SEQ)this.nTr++;
			if(MRUtils.getDataTypeOfThisMagneticResonanceSlice(fullHyperImage,can,0,0)== MRDataType.T2SEQ)this.nTe++;
		}
		this.t1Times=new double[nTr];this.nTr=0;
		this.t2Times=new double[nTe];this.nTe=0;
		for(int can=3;can<nbCanTot;can++) {
			if(MRUtils.getDataTypeOfThisMagneticResonanceSlice(fullHyperImage,can,0,0)== MRDataType.T1SEQ)
				this.t1Times[this.nTr++]=MRUtils.readTrInSliceLabel(fullHyperImage, can, 0, 0);
			if(MRUtils.getDataTypeOfThisMagneticResonanceSlice(fullHyperImage,can,0,0)== MRDataType.T2SEQ)
				this.t2Times[this.nTe++]=MRUtils.readTeInSliceLabel(fullHyperImage, can, 0, 0);
		}
		VitimageUtils.printImageResume(fullHyperImage);
		IJ.log("nTr="+this.nTr);
		IJ.log("nTe="+this.nTe);
		this.T1seq=new Duplicator().run(fullHyperImage,4,4+this.nTr-1,1,this.numberZ,1,this.numberTimes);
		this.T2seq=new Duplicator().run(fullHyperImage,4+this.nTr,4+this.nTr+this.nTe-1,1,this.numberZ,1,this.numberTimes);
		for(int i=0;i<this.nTr;i++) {
			this.T1seq.setC(i+1);
			this.T1seq.setDisplayRange(0, standardCapillaryM0);			
		}
		this.T1seq.setC(this.nTr-1);
		for(int i=0;i<this.nTe;i++) {
			this.T2seq.setC(i+1);
			this.T2seq.setDisplayRange(0, standardCapillaryM0);			
		}
		this.T2seq.setC(2);
		ImagePlus temp1=this.T1seq.duplicate();
		temp1.show();
		ImagePlus temp2=this.T2seq.duplicate();
		temp2.show();
		//3+
		
		//Set display parameters at start
		this.curEc1=this.nTr-2;
		this.curEc2=1;
		this.zCor=(int)Math.ceil(this.numberZ/2.0);
		this.tCor=1;

		//Gather informations about the actual day of each sequence in time lapse analysis
		this.actualDay=new int[this.numberTimes];
		for(int t=0;t<this.numberTimes;t++)this.actualDay[t]=t; //TODO : gather labels		
	}
	
	public void setupStructures() {
		initializeScreenConstants();		
		riceEstimator=RiceEstimator.getDefaultRiceEstimatorForNormalizedHyperEchoesT1AndT2Images();
		IJ.log("Starting MRI Curves explorer interface");

		//Data gathered on a voxel, on the neighbourhood of a voxel, and associated measured sigma values
		this.dataTimelapseT1=new double[this.numberTimes][this.nTr];
		this.dataTimelapseT2=new double[this.numberTimes][this.nTe];
		this.dataTimelapseT1Full=new double[this.numberTimes][this.nTr][1];
		this.dataTimelapseT2Full=new double[this.numberTimes][this.nTe][1];
		this.dataTimelapseT1Sigmas=new double[this.numberTimes][this.nTr];
		this.dataTimelapseT2SigmasMono=new double[this.numberTimes][this.nTe];
		this.dataTimelapseT2SigmasBicomp=new double[this.numberTimes][this.nTe];

		//Params estimated from fits
		this.paramsTimelapseT1=new double[this.numberTimes][2];
		this.paramsTimelapseT2=new double[this.numberTimes][12];
		this.paramsTimelapseT1Sigmas=new double[this.numberTimes][2];
		this.paramsTimelapseT2Sigmas=new double[this.numberTimes][12];

		//Corresponding curves and nice curves (for display)
		this.tabFittenT1=new double[this.numberTimes][nTr];
		this.tabFittenT2Mono=new double[this.numberTimes][nTe];
		this.tabFittenT2Bicomp=new double[this.numberTimes][nTe];
		this.tabFittenT2Tricomp=new double[this.numberTimes][nTe];
		this.tabFittenT1Cute=new double[this.numberTimes][nTr];
		this.tabFittenT2MonoCute=new double[this.numberTimes][nTe];
		this.tabFittenT2BicompCute=new double[this.numberTimes][nTe];
		this.tabFittenT2TricompCute=new double[this.numberTimes][nTe];
		this.timesT1Cute=MRUtils.getProportionalTimes(0,maxT1,xStepCuteT1);
		this.timesT2Cute=MRUtils.getProportionalTimes(0,maxT1,xStepCuteT2);

		//Estimation errors
		this.jitterT1=new double[this.numberTimes];
		this.jitterT2Mono=new double[this.numberTimes];
		this.jitterT2Bicomp=new double[this.numberTimes];
		this.jitterT2Tricomp=new double[this.numberTimes];

		//Khi and pvalue for estimation
		this.khi2T1=new double[this.numberTimes];
		this.khi2T2Mono=new double[this.numberTimes];
		this.khi2T2Bicomp=new double[this.numberTimes];
		this.khi2T2Tricomp=new double[this.numberTimes];
		this.pValT1=new double[this.numberTimes];
		this.pValT2Mono=new double[this.numberTimes];
		this.pValT2Bicomp=new double[this.numberTimes];
		this.pValT2Tricomp=new double[this.numberTimes];
		

		//Computed distribution of T1 and T2 values around the neighbourhood of a voxel
		this.pointsCoords=new int[this.numberTimes][1][8];
		this.correspondanceCanvas=new int[1000][2];
		this.valsSpectrum=new double[this.numberTimes][][];

		this.tabSigmasT1Seq=new double[numberTimes][numberZ];
		this.tabSigmasT2Seq=new double[numberTimes][numberZ];
		for(int t=0;t<numberTimes;t++) {
			for(int z=0;z<numberZ;z++) {
				VitimageUtils.printImageResume(T1seq);
				this.tabSigmasT1Seq[t][z]=MRUtils.readSigmaInSliceLabel(T1seq, 0, z, t)*0.5+MRUtils.readSigmaInSliceLabel(T1seq, 0, z, t)*0.5;
				this.tabSigmasT2Seq[t][z]=MRUtils.readSigmaInSliceLabel(T2seq, 0, z, t)*0.5+MRUtils.readSigmaInSliceLabel(T2seq, 1, z, t)*0.5;
			}
		}
		IJ.log(TransformUtils.stringMatrixMN("Tab des sigmas de T1", tabSigmasT1Seq));
		IJ.log(TransformUtils.stringMatrixMN("Tab des sigmas de T2", tabSigmasT2Seq));
	}
	
	public void startGui() {
	    WindowManager.addWindow(this);
	    instance = this;
		this.imgCan1=new ImageCanvas(T1seq);
		this.imgCan2=new ImageCanvas(T2seq);
        startPlotsAndRoi();
		initializeGUI2();
		imgCan2.getImage().setPosition(this.nTe,zCor,tCor);
		imgCan1.getImage().setPosition(this.curEc2+1,zCor,tCor);
		imgCan2.getImage().setPosition(1,zCor,tCor);
		imgCan1.getImage().setPosition(this.nTr,zCor,tCor);
		zoomLevel=imgCan1.getMagnification();
		this.T1seq.setDisplayRange(0, standardCapillaryM0);
		this.T2seq.setDisplayRange(0, standardCapillaryM0);
	}


		

	
	
	
	
	
	/** Gui building helpers*/	
	public void repaintAll() {
		imgCan1.repaint();
		imgCan2.repaint();
		plotCan1.repaint();
		plotCan2.repaint();
		plotCan21.repaint();
		plotCan22.repaint();
		repaint();
	}
		
	public void addComponentToPanel(JPanel panel,Component comp,GridBagConstraints gbc,int gridx, int gridy,int gridwidth,int gridheight,double weightx,double weighty) {
        gbc.gridx = gridx;
        gbc.gridy = gridy;
        gbc.gridwidth = gridwidth;
        gbc.gridheight = gridheight;
        gbc.weightx = weightx;
        gbc.weighty = weighty;
        gbc.fill = GridBagConstraints.BOTH;
        panel.add(comp, gbc);		
	}

	public void initializeScreenConstants() {
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		int screenY=(int) Math.round(screenSize.getHeight());
		int screenX=(int) Math.round(screenSize.getWidth());
		if(screenX>1920)screenX=1920;
		IJ.log("Screen resolution : "+screenX+" X "+screenY);
		if(screenX<1054)IJ.showMessage("Your screen has a very low resolution : "+screenX+" X "+screenY+"\nPlease consider investing in one which have at least 1024 lines.\nPlugin will run in survivor mode, everything can happen");
		int DY_TITLE=30;
		TARGET_HEIGHT_FOR_PLOTS_AND_IMGS=512;
		TARGET_WIDTH_FOR_PLOTS_AND_IMGS=512;
		DX_IMG=TARGET_WIDTH_FOR_PLOTS_AND_IMGS;
		DY_IMG=TARGET_HEIGHT_FOR_PLOTS_AND_IMGS;
		DY_TEXT=30;
		DELTA_X=10;
		DELTA_Y=10;
		DX_PLOT_1=(screenX-DX_IMG-2*DELTA_X)/2;
		DX_PLOT_2=(screenX-DX_IMG-2*DELTA_X)/2;;
		DY_PLOT_1=(2*TARGET_HEIGHT_FOR_PLOTS_AND_IMGS-3*DY_TEXT-4*+DELTA_Y)/2;
		DY_PLOT_2=(2*TARGET_HEIGHT_FOR_PLOTS_AND_IMGS-3*DY_TEXT-4*+DELTA_Y)/2;
		totalSizeX= DX_IMG+DX_PLOT_1+DX_PLOT_2+2*DELTA_X;  
		totalSizeY= 2*DY_IMG+DY_TITLE;  
	}
	
	public void initializeGUI2() {
		ImageJ ij = IJ.getInstance();
		setTitle("MRI Curve explorer V2");
		JPanel panel=new JPanel();
		panel.setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		

		//Image T1
		imgCan1.setMagnification(DX_IMG/dims[0]);
		addComponentToPanel(panel,imgCan1,gbc,0,0,DX_IMG,DY_IMG,0,0);
        imgCan1.removeKeyListener(ij);
        imgCan1.addKeyListener(this);
        imgCan1.removeMouseListener(ij);
        imgCan1.addMouseListener(this);

		//Image T2
        imgCan2.setMagnification(DX_IMG/dims[0]);
        addComponentToPanel(panel,imgCan2,gbc,0,DY_IMG,DX_IMG,DY_IMG,0,0);
		imgCan2.removeKeyListener(ij);
        imgCan2.addKeyListener(this);
        imgCan2.removeMouseListener(ij);
        imgCan2.addMouseListener(this);
		
		
		//Texts
		info1=new TextField("   Click on the image to compute T1",targetNumberCharsInTexts);info1.setEditable(false);info1.setFont(new Font("Helvetica", 0, 12));info1.setBackground(new Color(230,230,230));
		addComponentToPanel(panel,info1,gbc,DX_IMG+DELTA_X,DY_PLOT_1+DELTA_Y,DX_PLOT_1+DELTA_X+DX_PLOT_2,DY_TEXT,0,0);
		info3=new TextField("   Press 'h' stroke to get help",targetNumberCharsInTexts);info3.setEditable(false);info3.setFont(new Font("Helvetica", 0, 18));info3.setBackground(new Color(200,200,200));
		addComponentToPanel(panel,info3,gbc,DX_IMG+DELTA_X,DY_PLOT_1+DY_TEXT+2*DELTA_Y,DX_PLOT_1+DELTA_X+DX_PLOT_2,DY_TEXT,0,0);
		info2=new TextField("   Click on the image to compute T2",targetNumberCharsInTexts);info2.setEditable(false);info2.setFont(new Font("Helvetica", 0, 12));info2.setBackground(new Color(230,230,230));
		addComponentToPanel(panel,info2,gbc,DX_IMG+DELTA_X,DY_PLOT_1+2*DY_TEXT+3*DELTA_Y,DX_PLOT_1+DELTA_X+DX_PLOT_2,DY_TEXT,0,0);


		//PlotT1
		addComponentToPanel(panel,plotCan1,gbc,DX_IMG+DELTA_X,                           0,                      DX_PLOT_1,       DY_PLOT_1,    0,0);
        plotCan1.removeKeyListener(ij);
        plotCan1.addKeyListener(this);
        plotCan1.removeMouseListener(ij);
        plotCan1.addMouseListener(this);

		//PlotT2
		addComponentToPanel(panel,plotCan2,gbc,DX_IMG+DELTA_X,              DY_PLOT_1+3*DY_TEXT+4*DELTA_Y,       DX_PLOT_1,       DY_PLOT_1,   0,0);
        plotCan2.removeKeyListener(ij);
        plotCan2.addKeyListener(this);
        plotCan2.removeMouseListener(ij);
        plotCan2.addMouseListener(this);
        
		//PlotT21
        addComponentToPanel(panel,plotCan21,gbc,DX_IMG+DX_PLOT_1+2*DELTA_X,                0,                    DX_PLOT_2,       DY_PLOT_2,0,0);
        plotCan21.removeKeyListener(ij);
        plotCan21.addKeyListener(this);
        plotCan21.removeMouseListener(ij);
        plotCan21.addMouseListener(this);

		//PlotT22
        addComponentToPanel(panel,plotCan22,gbc,DX_IMG+DX_PLOT_1+2*DELTA_X,   DY_PLOT_1+3*DY_TEXT+4*DELTA_Y,     DX_PLOT_2,       DY_PLOT_2,0,0);
        plotCan22.removeKeyListener(ij);
        plotCan22.addKeyListener(this);
        plotCan22.removeMouseListener(ij);
        plotCan22.addMouseListener(this);

  		add(panel);        
  		setSize(totalSizeX,totalSizeY);
 		pack();
        this.setResizable(true);
		setVisible(true);
		repaint();
	}


	
	
	

	/** Gui updating functions and callbacks*/	
	@Override
	public void keyTyped(KeyEvent e) {
		IJ.log("KEY PRESSED : "+e.getKeyChar());

		
		
		///////ACTIONS TO CHANGE SIZE OF AREA
		if (e.getKeyChar()=='p' && statusRoi!=2) {
			this.autoSizingTimePlots=!this.autoSizingTimePlots;
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='z' && statusRoi!=2) {
			if(this.crossThick<dims[2] && statusRoi!=2)this.crossThick++;
			actualizeCursor();
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='s' && statusRoi!=2) {
			if(this.crossThick>0 && statusRoi!=2)this.crossThick--;
			actualizeCursor();
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='a' && statusRoi!=2) {
			if(this.crossWidth<dims[0])this.crossWidth++;
			actualizeCursor();
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='q' && statusRoi!=2) {
			if(this.crossWidth>0)this.crossWidth--;
			actualizeCursor();
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='e') {
			this.nRepetMonteCarlo+=20;
			IJ.log("Monte carlo go to N="+this.nRepetMonteCarlo+" repetitions");
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='d') {
			this.nRepetMonteCarlo-=20;
			if(this.nRepetMonteCarlo<0)this.nRepetMonteCarlo=0;
			IJ.log("Monte carlo go to N="+this.nRepetMonteCarlo+" repetitions");
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}

		
		
		///////ZOOM IN / OUT		
		if (e.getKeyChar()=='+') {
			if(statusRoi==2) {}//VitiDialogs.notYet("Warning : no zoom in / out during Roi mode");return;
			IJ.log("Zoom in");
			if(currentCanvas==WIN_T1 || currentCanvas==WIN_T2 ) {
				imgCan1.zoomIn(xMouse,yMouse);
				imgCan2.zoomIn(xMouse,yMouse);
				zoomLevel=imgCan1.getMagnification();
				xMouse=dims[0]/2;yMouse=dims[1]/2;
				actualizeCursor();
			}
			if(currentCanvas==WIN_PLOT1 || currentCanvas==WIN_PLOT2) {
				plotCan1.zoomIn(xMouse,yMouse);
				plotCan2.zoomIn(xMouse,yMouse);
			}
		}
		if (e.getKeyChar()=='-') {
			if(statusRoi==2) {}//VitiDialogs.notYet("Warning : no zoom in / out during Roi mode");return;
			IJ.log("Zoom out");
			if( (currentCanvas==WIN_T1 || currentCanvas==WIN_T2) ) {
				imgCan1.setMagnification(imgCan1.getMagnification()/2.0);
				imgCan2.setMagnification(imgCan2.getMagnification()/2.0);
				imgCan1.zoomIn(xMouse,yMouse);
				imgCan2.zoomIn(xMouse,yMouse);
				zoomLevel=imgCan1.getMagnification();
				actualizeCursor();
			}
			if(currentCanvas==WIN_PLOT1 || currentCanvas==WIN_PLOT2) {
				zoomLevel=1;
				plotCan1.setMagnification(zoomLevel);
				plotCan2.setMagnification(zoomLevel);
			}
		}


		
		
		///////LUT	
		if (e.getKeyChar()=='f') {
			isFireLut=!isFireLut;
		}
		///////DISPLACEMENT IN HYPERIMAGE		
		if (e.getKeyChar()=='8' || e.getKeyChar()=='2' || e.getKeyChar()=='4' || e.getKeyChar()=='6' || e.getKeyChar()=='1' || e.getKeyChar()=='3' || e.getKeyChar()=='7' || e.getKeyChar()=='9'  || e.getKeyChar()=='f') {
			//Move t
			if((e.getKeyChar()=='8') && (tCor<numberTimes) )tCor+=1;
			if((e.getKeyChar()=='2') && (tCor>1) )tCor-=1;

			//Move z
			if((e.getKeyChar()=='6') && (zCor<numberZ) )zCor+=1;
			if((e.getKeyChar()=='4') && (zCor>1) )zCor-=1;

			//Move Tr of T1
			if((e.getKeyChar()=='7') && (curEc1<this.nTr-1) )curEc1+=1;
			if((e.getKeyChar()=='1') && (curEc1>0) )curEc1-=1;

			//Move echo of T2
			if((e.getKeyChar()=='9') && (curEc2<this.nTe-1) )curEc2+=1;
			if((e.getKeyChar()=='3') && (curEc2>0) )curEc2-=1;

			
			IJ.log("Changing coordinates to t="+tCor+"/"+this.numberTimes+"  , z="+zCor+"/"+this.numberZ+" Tr="+curEc1+"  Te="+curEc2);
			imgCan2.getImage().setPosition(this.nTr,zCor,tCor);
			imgCan2.getImage().setPosition(1,zCor,tCor);
			imgCan2.getImage().setPosition(this.curEc2+1,zCor,tCor);
			imgCan1.getImage().setPosition(1,zCor,tCor);
			imgCan1.getImage().setPosition(this.nTe,zCor,tCor);
			imgCan1.getImage().setPosition(this.curEc1+1,zCor,tCor);
			if(isFireLut) {
				IJ.run(imgCan1.getImage(),"Fire","");//.setLut();
				IJ.run(imgCan2.getImage(),"Fire","");//.setLut();
			}
			else {
				IJ.run(imgCan1.getImage(),"Grays","");//.setLut();
				IJ.run(imgCan2.getImage(),"Grays","");//.setLut();

			}
			
			if(statusRoi!=2 || (userRoi ==null)) {				
				actualizeCursor();
			}
			else {
				this.imgCan1.setOverlay(new Overlay(userRoi));
				this.imgCan2.setOverlay(new Overlay(userRoi));
			}
			if(e.getKeyChar()=='4' || e.getKeyChar()=='6') {
				actualizeMriObservationsBasedOnData();
				actualizeSpectrumCurves();
				this.actualizeRangingBoundaries();
				this.actualizeRangedData();
				computeFitsAgain();
			}
			if(e.getKeyChar()=='2' || e.getKeyChar()=='8') {
				this.actualizeRangingBoundaries();
				this.actualizeRangedData();
			}
			displayFitsAgain();
		}


		

		
		
		
		///////USER CUSTOM ROI MANAGEMENT
		if (e.getKeyChar()=='r') {
			if(statusRoi==0) {				
				zoomLevel=imgCan1.getMagnification();
				if(zoomLevel>1) {
					while(zoomLevel>1) {
						imgCan1.setMagnification(imgCan1.getMagnification()/2.0);
						imgCan2.setMagnification(imgCan2.getMagnification()/2.0);
						imgCan1.zoomIn(xMouse,yMouse);
						imgCan2.zoomIn(xMouse,yMouse);
						repaint();
						plotCan1.repaint();
						plotCan21.repaint();
						imgCan1.repaint();
						plotCan2.repaint();
						plotCan22.repaint();
						imgCan2.repaint();
						zoomLevel=imgCan1.getMagnification();						
					}
				}
				else if(zoomLevel<1) {
					while(zoomLevel<1) {
						imgCan1.zoomIn(xMouse,yMouse);
						imgCan2.zoomIn(xMouse,yMouse);
						repaint();
						plotCan1.repaint();
						plotCan21.repaint();
						imgCan1.repaint();
						plotCan2.repaint();
						plotCan22.repaint();
						imgCan2.repaint();
						zoomLevel=imgCan1.getMagnification();						
					}
				}					
				IJ.log("Please select a ROI now, add it to the Roi manager, then hit the 'r' strike again");
				String file=VitiDialogs.chooseOneRoiPathUI("Choose a .roi file", "");
				if(file==null)return;
				this.userRoi=new Opener().openRoi(file);
				IJ.selectWindow("MRI Curve explorer V2");
				statusRoi=2;
//				zCor=this.userRoi.getZPosition()+1;
//				tCor=this.userRoi.getTPosition();
				this.userRoi.setPosition(0);
				imgCan2.getImage().setPosition(3,zCor,tCor);
				imgCan1.getImage().setPosition(1,zCor,tCor);
				imgCan2.getImage().setPosition(1,zCor,tCor);
				imgCan1.getImage().setPosition(3,zCor,tCor);
				imgCan1.repaint();
				imgCan2.repaint();
				this.imgCan1.setOverlay(new Overlay(userRoi));
				this.imgCan2.setOverlay(new Overlay(userRoi));
				repaint();
				actualizeMriObservationsBasedOnData();
				actualizeSpectrumCurves();
				computeFitsAgain();
				displayFitsAgain();
				IJ.log("Please hit the 'r' strike again to quit the Roi mode");
			}
			else if(statusRoi==2) {
				statusRoi=0;
				this.userRoi=null;
				this.nPtsRoi=1;
				actualizeCursor();
				actualizeMriObservationsBasedOnData();
				actualizeSpectrumCurves();
				computeFitsAgain();
				displayFitsAgain();
			}
			else {
				actualizeCursor();
				actualizeMriObservationsBasedOnData();
				actualizeSpectrumCurves();
				computeFitsAgain();
				displayFitsAgain();
			}
		}
		
		
		
		///MISCELLANEOUS
		if (e.getKeyChar()==',') {
			//just update
			this.rangingFactor=this.rangingFactor*0.7;
			this.actualizeRangingBoundaries();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()==';') {
			//just update
			this.rangingFactor=this.rangingFactor*1.4;
			this.actualizeRangingBoundaries();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='u') {
			//just update
			this.actualizeRangingBoundaries();
			this.actualizeRangedData();
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='b') {
			//just update
			this.rangeDisplayedFullBlocks=!this.rangeDisplayedFullBlocks;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			actualizeSpectrumCurves();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='c') {
			this.multiThreaded=!this.multiThreaded;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			actualizeSpectrumCurves();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='v') {
			this.separateNormalizationSpectrum=!this.separateNormalizationSpectrum;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			actualizeSpectrumCurves();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='l') {
			if(this.fitAlgorithm==MRUtils.LM) {
				this.fitAlgorithm=MRUtils.SIMPLEX;
				IJ.log("Switching to SIMPLEX");
			}
			else if(this.fitAlgorithm==MRUtils.SIMPLEX) {
				this.fitAlgorithm=MRUtils.LM;
				IJ.log("Switching to Levenberg-Marquardt");
			}
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			actualizeSpectrumCurves();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='t') {
			thickCurve=1-thickCurve;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			actualizeSpectrumCurves();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='g') {
			gaussianSpectrum=(gaussianSpectrum+1)%4;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			actualizeSpectrumCurves();
			this.actualizeRangedData();
			computeFitsAgain();
			displayFitsAgain();
		}
		if (e.getKeyChar()=='h') {
			VitiDialogs.getYesNoUI("Curve explorer help",
			getHelpString());
		}
		if (e.getKeyChar()=='e') {
			IJ.log("\n\nExporting data for imporation in R, Python, or Matlab/Octave");
			IJ.log(" --------- R EXPORT ------");
			IJ.log(rSentence);
			IJ.log(" --------- PYTHON EXPORT ------");
			IJ.log(pySentence);
			IJ.log(" --------- MATLAB/OCTAVE EXPORT ------");
			IJ.log(matSentence);
			IJ.log(" ---------------");
		}
		
		actualizeCursor();
		repaint();
		plotCan1.repaint();
		plotCan21.repaint();
		imgCan1.repaint();
		plotCan2.repaint();
		plotCan22.repaint();
		imgCan2.repaint();
	}
	
	public void actualizeCursor() {
		double xD=0;
		double yD=0;
		if(this.currentCanvas==WIN_T1) {
			xCor=imgCan1.offScreenX(xMouse);
			yCor=imgCan1.offScreenY(yMouse);
			xD=imgCan1.offScreenXD(xMouse);
			yD=imgCan1.offScreenYD(yMouse);
		}
		if(this.currentCanvas==WIN_T2) {
			xCor=imgCan2.offScreenX(xMouse);
			yCor=imgCan2.offScreenY(yMouse);
			xD=imgCan2.offScreenXD(xMouse);
			yD=imgCan2.offScreenYD(yMouse);
		}
		else {
			xD=imgCan1.offScreenXD(xMouse);
			yD=imgCan1.offScreenYD(yMouse);
		}
		int xMouseCenter=(int) Math.round(xMouse+  (   Math.floor(xD)+0.5-xD  )*zoomLevel);//screen location of square center, according to xMouse	
		int yMouseCenter=(int) Math.round(yMouse+  (   Math.floor(yD)+0.5-yD  )*zoomLevel);
		double sizeOnScreen=(0.5+crossWidth)*zoomLevel;//half radius

		PointRoi prT1=new PointRoi(xMouseCenter,yMouseCenter,sizeOfCursor()+" yellow hybrid");
		Overlay overT1=new Overlay(prT1);
		Overlay overT2=new Overlay(prT1);
		if(statusRoi!=2) {
			overT1.add(new Roi(xMouseCenter-sizeOnScreen,yMouseCenter-sizeOnScreen,2*sizeOnScreen,2*sizeOnScreen));
			overT2.add(new Roi(xMouseCenter-sizeOnScreen,yMouseCenter-sizeOnScreen,2*sizeOnScreen,2*sizeOnScreen));
		}
		else if(zoomLevel==1) {
			overT1.add(userRoi);
			overT2.add(userRoi);
		}
		if(this.spectrumRangingMode) {
			for(int pt=0;pt<this.rangeRoiPoints.size();pt++) {
				int dx=rangeRoiPoints.get(pt)[0]-xCor;//Relative coordinates to the cross
				int dy=rangeRoiPoints.get(pt)[1]-yCor;
				Roi r=new Roi(xMouseCenter+(dx-0.5)*zoomLevel,yMouseCenter+(dy-0.5)*zoomLevel,zoomLevel,zoomLevel);
				if(rangeDisplayedFullBlocks)r.setFillColor(new Color(0,0,255));
				else r.setStrokeColor(new Color(0,0,255));
				if(this.rangingBoundaries[1]<500)overT2.add(r);
				else overT1.add(r);
			}
		}
		imgCan1.setOverlay(overT1);
		imgCan2.setOverlay(overT2);
		imgCan1.repaint();
		imgCan2.repaint();
	}
	
	public String sizeOfCursor() {
		return "medium";
	}
		
	public void actualizeRangingBoundaries() {
		double borderLeft=76;
		double borderRight=20;
		double lengPlot=DX_PLOT_2-borderRight-borderLeft;
		double xPos=this.xMouseRange-borderLeft;
		double[]vals=plotT22.getLimits();
		double factMul=vals[1]/vals[0];
		double rangCenter=Math.pow(factMul,xPos*1.0/lengPlot)*vals[0];
		this.rangingBoundaries=new double[] {rangCenter*(1-this.rangingFactor),rangCenter,rangCenter*(1+this.rangingFactor)};
	}

	public void actualizeRangedData() {
		this.rangeRoiPoints=new ArrayList<int[]>();
		if(! this.spectrumRangingMode)return;
		for(int p=0;p<pointsCoords[tCor-1].length;p++) {
			if((this.rangingBoundaries[1]>=500)) {
				if ( ( (this.pointsCoords[tCor-1][p][2]<this.rangingBoundaries[2]) && (this.pointsCoords[tCor-1][p][2]>this.rangingBoundaries[0])) &&  this.pointsCoords[tCor-1][p][3]==1){
					//Ranging T1 values. Current Roi point has its value in the range
					rangeRoiPoints.add(new int[] {this.pointsCoords[tCor-1][p][0],this.pointsCoords[tCor-1][p][1]});
				}
			}
			else{
				if ( (this.pointsCoords[tCor-1][p][5]==1) && ( (this.pointsCoords[tCor-1][p][4]<this.rangingBoundaries[2]) && (this.pointsCoords[tCor-1][p][4]>this.rangingBoundaries[0]))
					|| ( (this.pointsCoords[tCor-1][p][7]==1) && (  (this.pointsCoords[tCor-1][p][6]<this.rangingBoundaries[2]) && (this.pointsCoords[tCor-1][p][6]>this.rangingBoundaries[0]))  ) ){
					//Ranging T2 values. Current Roi point has its value in the range
					rangeRoiPoints.add(new int[] {this.pointsCoords[tCor-1][p][0],this.pointsCoords[tCor-1][p][1]});
				}
			}
		}
	}

	public void startPlotsAndRoi(){
		IJ.log("Starting Plots");
		plotT1 = new Plot("T1 curve explorer","Recovery time","Spin-echo magnitude signal");
		plotT2 = new Plot("T2 curve explorer","Echo time","Spin-echo magnitude signal");
		plotT1.changeFont(new Font("Helvetica", 0, 14));
		plotT2.changeFont(new Font("Helvetica", 0, 14));
		plotT1.setSize(DX_PLOT_1,DY_PLOT_1);
		plotT2.setSize(DX_PLOT_2,DY_PLOT_1);
		plotCan1=new PlotCanvas(plotT1.getImagePlus());
		plotCan2=new PlotCanvas(plotT2.getImagePlus());
		plotCan1.setPlot(plotT1);
		plotCan2.setPlot(plotT2);

		plotT21 = new Plot("T1 timelapse tracker","Estimated T1 (red) and T2 (green)","Observation time",Plot.DEFAULT_FLAGS-Plot.X_GRID-Plot.Y_GRID+Plot.X_LOG_TICKS+Plot.X_LOG_NUMBERS);
		plotT22 = new Plot("T2 timelapse tracker","T1 and T2 spectrum over area (M0-weighted)","Observation time",Plot.DEFAULT_FLAGS-Plot.X_GRID-Plot.Y_GRID+Plot.X_LOG_TICKS+Plot.X_LOG_NUMBERS);
		plotT21.changeFont(new Font("Helvetica", 0, 14));
		plotT22.changeFont(new Font("Helvetica", 0, 14));
		plotT21.setBackgroundColor(new Color(0,0,0));
		plotT22.setBackgroundColor(new Color(0,0,0));
		plotT21.setSize(DX_PLOT_2,DY_PLOT_2);
		plotT22.setSize(DX_PLOT_2,DY_PLOT_2);
		plotCan21=new PlotCanvas(plotT21.getImagePlus());
		plotCan22=new PlotCanvas(plotT22.getImagePlus());
		plotCan21.setPlot(plotT21);
		plotCan22.setPlot(plotT22);
		for(int i =0;i<maxCurves;i++){
			plotT1.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
			plotT2.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
			plotT21.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
			plotT22.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
		}
		plotT1.setLimits(0, maxT1, 0, maxPlotYT1);
		plotT2.setLimits(0, maxT2, 0, maxPlotYT2);
		plotT21.setLimits(t0T1,t1T1, 0, numberTimes);
		plotT22.setLimits(t0T2,t1T2, 0, numberTimes);
	}

	public void closePlotsAndRoi(){
		java.awt.Window win;
		win=WindowManager.getWindow("T1 curve explorer");
		if(win != null){
			IJ.selectWindow("T1 curve explorer");
			WindowManager.getCurrentWindow().close();
		}
		win=WindowManager.getWindow("T2 curve explorer");
		if(win != null){
			IJ.selectWindow("T2 curve explorer");
			WindowManager.getCurrentWindow().close();
		}
		RoiManager rm = RoiManager.getInstance();
		if(rm != null) rm.close();		
		IJ.log("Fin de l'exploration");
	}

	
	
	
	
	

	
	/** Read new data and update estimation and fits*/	
	public void actualizeMriObservationsBasedOnData() {
		for(int tim=1;tim<=this.numberTimes;tim++) {
			if(statusRoi!=2 || userRoi==null) {
				this.nPtsCur=statusRoi==2 ? this.nPtsRoi : (1+2*crossWidth)*(1+2*crossWidth)*(1+2*crossThick);
				this.pointsCoords[tim-1]=new int[nPtsCur][8];
				this.dataTimelapseT1Full[tim-1]=MRUtils.getFullDataForVoxel(this.imgCan1.getImage(),(int)xCor,(int)yCor,tim-1,this.numberTimes,zCor-1,numberZ,this.nTr,this.crossWidth,this.crossThick,this.gaussianWeighting);
				this.dataTimelapseT2Full[tim-1]=MRUtils.getFullDataForVoxel(this.imgCan2.getImage(),(int)xCor,(int)yCor,tim-1,this.numberTimes,zCor-1,numberZ,this.nTe,this.crossWidth,this.crossThick,gaussianWeighting);

				this.pointsCoords[tim-1]=MRUtils.getCoordsOfCorrespondingVoxelsUsedInEstimationAroundThisPoint(this.imgCan1.getImage(),(int)xCor,(int)yCor,tim-1,this.numberTimes,zCor-1,numberZ,this.nTr,this.crossWidth,this.crossThick,this.gaussianWeighting);
			}
			else {
				IJ.log("Actualizing data from Roi, at coordinates Z="+zCor+" , T="+tCor);
				int[][]coordsRoi=VitimageUtils.getRoiAsCoords(this.userRoi);
				this.nPtsRoi=coordsRoi.length;
				this.nPtsCur=this.nPtsRoi;
				double[]accumulatorT1=new double[nTr];
				double[]accumulatorT2=new double[nTe];
				this.dataTimelapseT1[tim-1]=new double[nTr];
				this.dataTimelapseT2[tim-1]=new double[nTe];
				this.dataTimelapseT1Full[tim-1]=new double[nTr][nPtsRoi];
				this.dataTimelapseT2Full[tim-1]=new double[nTe][nPtsRoi];
				this.pointsCoords [tim-1]=new int[nPtsRoi][8];
				for(int pt=0;pt<coordsRoi.length;pt++) {
					accumulatorT1=MRUtils.getDataForVoxel(this.imgCan1.getImage(),coordsRoi[pt][0],coordsRoi[pt][1],tim-1,this.numberTimes,zCor-1,numberZ,this.nTr,0,0,gaussianWeighting);
					accumulatorT2=MRUtils.getDataForVoxel(this.imgCan2.getImage(),coordsRoi[pt][0],coordsRoi[pt][1],tim-1,this.numberTimes,zCor-1,numberZ,this.nTe,0,0,gaussianWeighting);					
					for(int tr=0;tr<nTr;tr++) this.dataTimelapseT1Full[tim-1][tr][pt]=accumulatorT1[tr];
					for(int ec=0;ec<nTe;ec++) this.dataTimelapseT2Full[tim-1][ec][pt]=accumulatorT2[ec];
					this.pointsCoords[tim-1][pt]=new int[] {coordsRoi[pt][0],coordsRoi[pt][1],0,0,0,0,0,0};
				}
			}
			for(int tr=0;tr<nTr;tr++) {
				this.dataTimelapseT1[tim-1][tr]=VitimageUtils.statistics1D(this.dataTimelapseT1Full[tim-1][tr])[0];
				if(this.nPtsCur>1) {
					this.dataTimelapseT1Sigmas[tim-1][tr]=VitimageUtils.statistics1D(this.dataTimelapseT1Full[tim-1][tr])[1];
				}
			}
			for(int ec=0;ec<nTe;ec++) {
				this.dataTimelapseT2[tim-1][ec]=VitimageUtils.statistics1D(this.dataTimelapseT2Full[tim-1][ec])[0];
				if(this.nPtsCur>1) {
					this.dataTimelapseT2SigmasMono[tim-1][ec]=VitimageUtils.statistics1D(this.dataTimelapseT2Full[tim-1][ec])[1];
					this.dataTimelapseT2SigmasBicomp[tim-1][ec]=VitimageUtils.statistics1D(this.dataTimelapseT2Full[tim-1][ec])[1];
				}
			}
		}
	}
		
	public Object[] fitAndEvaluate(double[]tabTimes,double[]tabTimesCute,double[]tabData,double sigmaRice,int fitAlgorithm,int fitCurveType,int nbMonteCarloSimulations,int nbPts) {
		int nParams=(fitCurveType==MRUtils.T1_RECOVERY_RICE) ? 2 : (fitCurveType==MRUtils.T2_RELAX_RICE) ? 2 :(fitCurveType==MRUtils.MULTICOMP_RICE) ? 4 : 5; 
		boolean isT1=(fitCurveType==MRUtils.T1_RECOVERY); 
		double []estimatedParams=new double[nParams];
		double []estimatedSigmas=new double[nParams];
		double [][]estimation=MRUtils.makeFit(tabTimes, tabData,fitCurveType,fitAlgorithm,100,sigmaRice,nbMonteCarloSimulations,nbPts,this.riceEstimator);
		for(int i=0;i<nParams;i++) {estimatedParams[i]=estimation[0][i] ;  estimatedSigmas[i]=estimation[1][i];}
		double []tabFitten=MRUtils.fittenRelaxationCurve(tabTimes,estimatedParams,sigmaRice,fitCurveType);
		double []tabFittenCute=MRUtils.fittenRelaxationCurve(tabTimesCute,estimatedParams,sigmaRice,fitCurveType);
		double[]accs=MRUtils.fittingAccuracies(tabData,tabTimes,sigmaRice,estimatedParams,fitCurveType,false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
		if((nParams==4) && (estimatedParams[1]>estimatedParams[3])) {
			double tmp=estimatedParams[1];estimatedParams[1]=estimatedParams[3];estimatedParams[3]=tmp;
			tmp=estimatedParams[0];estimatedParams[0]=estimatedParams[2];estimatedParams[2]=tmp;
			tmp=estimatedSigmas[1];estimatedSigmas[1]=estimatedSigmas[3];estimatedSigmas[3]=tmp;
			tmp=estimatedSigmas[0];estimatedSigmas[0]=estimatedSigmas[2];estimatedSigmas[2]=tmp;
		}		
		double jitter=0;
		double[]statsRice=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(sigmaRice);
		if( tabData[0]<(statsRice[0]+3*statsRice[1]) ) {jitter=99;accs[1]=100;}
		if(isT1 && (estimatedParams[0]<0 || estimatedParams[1]<minAcceptableT1 || estimatedParams[1]>maxT1 || estimatedParams[0]>maxAcceptableM0T1  ) )  {jitter=99;accs[1]=100;}
		if(!isT1 && (nParams==2) && (estimatedParams[0]<0 || estimatedParams[1]<0 || estimatedParams[0]>maxAcceptableM0T2 || estimatedParams[1]>maxAcceptableT2) ) {jitter=99;accs[1]=100;}
		if(!isT1 && (nParams==4) && (estimatedParams[0]<0 || estimatedParams[1]<minAcceptableT2 || 
				estimatedParams[2]<0 || estimatedParams[3]>maxAcceptableT2 || 
				estimatedParams[0]>maxAcceptableM0T2 ||   estimatedParams[2]>maxAcceptableM0T2) ) {jitter=99;accs[1]=100;}		
		return new Object[] {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accs[0],accs[1],jitter};
	}
	
	public double[]computeAllFitsForOnePoint(double[]dataT1,double[]dataT2,double sigmaT1,double sigmaT2) {
		double[]tabRet=new double[9];//M0T1, T1, active,     M0T21, T21, active    M0T22, T22, active 
 		//Estimer T1 Monocomp
		Object[] obj=fitAndEvaluate(t1Times,this.timesT1Cute,dataT1,sigmaT1,this.fitAlgorithm,MRUtils.T1_RECOVERY_RICE,0,1);
 		tabRet[0]=((double[]) obj[0])[0]; 		tabRet[1]=((double[]) obj[0])[1];
 		tabRet[2]=( (double) obj[6] >=99 ? 0 : 1 );

 		//Estimer T2 Monocomp
 		obj=fitAndEvaluate(t2Times,this.timesT2Cute,dataT2,sigmaT2,this.fitAlgorithm,MRUtils.T2_RELAX_RICE,0,1);
 		double M0mono=((double[]) obj[0])[0]; 		double T2mono=((double[]) obj[0])[1];
 		double khimono=((double) obj[4]);  double jimono=((double) obj[6]);
	
 		//Estimer T2 Bicomp
 		obj=fitAndEvaluate(t2Times,this.timesT2Cute,dataT2,sigmaT2,this.fitAlgorithm,MRUtils.MULTICOMP_RICE,0,1);
 		double M01bi=((double[]) obj[0])[0]; 		double T21bi=((double[]) obj[0])[1];
 		double M02bi=((double[]) obj[0])[2]; 		double T22bi=((double[]) obj[0])[3];
 		double khibi=((double) obj[4]);  double jibi=((double) obj[6]);
 		if(khimono>khibi) { 
 			if(jibi<99) { 	
 				tabRet[3]=M01bi; tabRet[4]=T21bi ; tabRet[5]=1; tabRet[6]=M02bi; tabRet[7]=T22bi ; tabRet[8]=1; 	
			}
 			else if(jimono<99) { 			tabRet[3]=M0mono; tabRet[4]=T2mono ; tabRet[5]=1;	} 	
		}
 		else {	if(jimono<99) { 				tabRet[3]=M0mono; tabRet[4]=T2mono ; tabRet[5]=1;	} 		}
 		return tabRet;
	}		
	
	public void actualizeEstimations2() {
		for(int tim=1;tim<= this.numberTimes;tim++) {
		if(tim==tCor)IJ.log("  --> Estimations for current time "+tim+" / "+this.numberTimes+" at z="+zCor);
	 		//Estimer T1
	 		if(tim==tCor)IJ.log("    T1 MRI times = "+TransformUtils.stringVectorNDou(t1Times, "")+" with sigma rice = "+VitimageUtils.dou(this.tabSigmasT1Seq[tim-1][zCor-1])+"\n    T1 MRI data = "+TransformUtils.stringVectorNDou(dataTimelapseT1[tim-1], ""));
	 		Object[] obj=fitAndEvaluate(t1Times,this.timesT1Cute,dataTimelapseT1[tim-1],this.tabSigmasT1Seq[tim-1][zCor-1],this.fitAlgorithm,MRUtils.T1_RECOVERY_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
	 		paramsTimelapseT1[tim-1]=(double[]) obj[0];
	 		paramsTimelapseT1Sigmas[tim-1]=(double[]) obj[1];
	 		tabFittenT1[tim-1]=(double[]) obj[2];
	 		tabFittenT1Cute[tim-1]=(double[]) obj[3];
	 		khi2T1[tim-1]=(double) obj[4];
	 		pValT1[tim-1]=(double) obj[5];
	 		jitterT1[tim-1]=(double) obj[6];
	 		if(this.nPtsCur==1)this.dataTimelapseT1Sigmas[tim-1]=riceEstimator.estimateSigmas(this.tabSigmasT1Seq[tim-1][zCor-1], this.tabFittenT1[tim-1]);

	 		//Estimer T2 Monocomp
	 		if(tim==tCor)IJ.log("    T2 MRI times = "+TransformUtils.stringVectorNDou(t2Times, "")+" with sigma rice = "+VitimageUtils.dou(this.tabSigmasT2Seq[tim-1][zCor-1])+"\n    T2 MRI data = "+TransformUtils.stringVectorNDou(dataTimelapseT2[tim-1], ""));
	 		obj=fitAndEvaluate(t2Times,this.timesT2Cute,dataTimelapseT2[tim-1],this.tabSigmasT2Seq[tim-1][zCor-1],this.fitAlgorithm,MRUtils.T2_RELAX_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
	 		for(int prm=0;prm<2;prm++) {
	 			paramsTimelapseT2[tim-1][prm]=((double[]) obj[0])[prm];
	 			paramsTimelapseT2Sigmas[tim-1][prm]=((double[]) obj[1])[prm];
	 		}
	 		tabFittenT2Mono[tim-1]=(double[]) obj[2];
	 		tabFittenT2MonoCute[tim-1]=(double[]) obj[3];
	 		khi2T2Mono[tim-1]=(double) obj[4];
	 		pValT2Mono[tim-1]=(double) obj[5];
	 		jitterT2Mono[tim-1]=(double) obj[6];
	 		if(this.nPtsCur==1)this.dataTimelapseT2SigmasMono[tim-1]=riceEstimator.estimateSigmas(this.tabSigmasT2Seq[tim-1][zCor-1], this.tabFittenT2Mono[tim-1]);
		
	 		//Estimer T2 Bicomp
	 		obj=fitAndEvaluate(t2Times,this.timesT2Cute,dataTimelapseT2[tim-1],this.tabSigmasT2Seq[tim-1][zCor-1],this.fitAlgorithm,MRUtils.MULTICOMP_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
	 		for(int prm=2;prm<6;prm++) {
	 			paramsTimelapseT2[tim-1][prm]=((double[]) obj[0])[prm-2];
	 			paramsTimelapseT2Sigmas[tim-1][prm]=((double[]) obj[1])[prm-2];
	 		}
	 		tabFittenT2Bicomp[tim-1]=(double[]) obj[2];
	 		tabFittenT2BicompCute[tim-1]=(double[]) obj[3];
	 		khi2T2Bicomp[tim-1]=(double) obj[4];
	 		pValT2Bicomp[tim-1]=(double) obj[5];
	 		jitterT2Bicomp[tim-1]=(double) obj[6];
	 		if(this.nPtsCur==1)this.dataTimelapseT2SigmasBicomp[tim-1]=riceEstimator.estimateSigmas(this.tabSigmasT2Seq[tim-1][zCor-1], this.tabFittenT2Bicomp[tim-1]);
		}		
	}

	public double[][] getDisplayableSpectrumCurveMultiThread(int time) {
		int numberBins=150 /(1+gaussianSpectrum) ;//over 3 decades it makes every 5% or every 10 %
		final double[]binInf=new double[numberBins];
		final double[]binSup=new double[numberBins];
		final double[]binMed=new double[numberBins];
		double[]histoM0T1=new double[numberBins];
		double[]histoM0T2=new double[numberBins];
		double [][]output=new double[3][numberBins];
		final double[][]inputsAllT1=new double[this.nPtsCur][this.nTr];
		final double[][]inputsAllT2=new double[this.nPtsCur][this.nTe];
		final double[]sigmaT1=new double[this.nPtsCur];
		final double[]sigmaT2=new double[this.nPtsCur];
		final double[][]valsFitAll=new double[this.nPtsCur][9];
		double minBin=t0T2;
		double maxBin=t1T1;
		double multFromMinToMax=maxBin/minBin;
		for(int i=0;i<numberBins;i++) {
			binInf[i]=minBin*Math.pow(multFromMinToMax, i*1.0/numberBins);
			binSup[i]=minBin*Math.pow(multFromMinToMax, (i+1)*1.0/numberBins);
			binMed[i]=binSup[i]*0.5+binInf[i]*0.5;
		}

		
		final int stepPt=Math.max(this.nPtsCur/20, 1);

		//Prepare data
		for(int p=0;p<this.nPtsCur;p++) {
			for(int tr=0;tr<this.nTr;tr++)inputsAllT1[p][tr]=this.dataTimelapseT1Full[time][tr][p];
			for(int te=0;te<this.nTe;te++)inputsAllT2[p][te]=this.dataTimelapseT2Full[time][te][p];
			sigmaT1[p]=this.tabSigmasT1Seq[time][zCor-1];
			sigmaT2[p]=this.tabSigmasT2Seq[time][zCor-1];
		}
		int nProc=VitimageUtils.getNbCores()-1;
		final int[][] listProcs=VitimageUtils.listForThreads(this.nPtsCur,nProc);
		final int nbTot=this.nPtsCur;
		AtomicInteger atomNumPt=new AtomicInteger();
		AtomicInteger atomNumThread=new AtomicInteger();
		final Thread[] threads = VitimageUtils.newThreadArray(nProc);  
		for (int ithread = 0; ithread < nProc; ithread++) {  
			
			threads[ithread] = new Thread() {  { setPriority(Thread.NORM_PRIORITY); }  
			public void run() {
				int thread=atomNumThread.getAndIncrement();
				int nPtThread=listProcs[thread].length;
				for(int pt=0;pt<nPtThread;pt++) {
					int number=listProcs[thread][pt];
					valsFitAll[number]=computeAllFitsForOnePoint(inputsAllT1[number],inputsAllT2[number],sigmaT1[number],sigmaT2[number]);
					int nt=atomNumPt.getAndIncrement();
					if(nt%stepPt==0 && nPtThread>1) {IJ.log((100*nt)/nbTot + " %");IJ.showProgress(nt/nbTot);}					
					}
				} //fin run
				};
			}

		VitimageUtils.startAndJoin(threads);  

			
		//Collect results
		for(int p=0;p<this.nPtsCur;p++) {
			//Fit T1
			if(valsFitAll[p][2]>=1) {
				for(int bin=0;bin<numberBins;bin++)if( (valsFitAll[p][1]<binSup[bin]) && (valsFitAll[p][1]>=binInf[bin]) ) {
					histoM0T1[bin]+=valsFitAll[p][0];
					this.pointsCoords[time][p][2]=(int)Math.round(binMed[bin]);this.pointsCoords[time][p][3]=1;
				}
			}

			//Fit T2 Mono et bicomp
			if(valsFitAll[p][5]>=1) {
				for(int bin=0;bin<numberBins;bin++)if( (valsFitAll[p][4]<binSup[bin]) && (valsFitAll[p][4]>=binInf[bin]) ) {
					histoM0T2[bin]+=valsFitAll[p][3];
					this.pointsCoords[time][p][4]=(int)Math.round(binMed[bin]);this.pointsCoords[time][p][5]=1;
				}
			}
			if(valsFitAll[p][8]>=1) {
				for(int bin=0;bin<numberBins;bin++)if( (valsFitAll[p][7]<binSup[bin]) && (valsFitAll[p][7]>=binInf[bin]) ) {
					histoM0T2[bin]+=valsFitAll[p][6];
					this.pointsCoords[time][p][6]=(int)Math.round(binMed[bin]);this.pointsCoords[time][p][7]=1;
				}
			}
		}
				
		//Smooth this histogram with a factor to be defined, maybe depending on the estimation error on parameters
//		output[0]=smoothHisto(histoM0T1);
//		output[1]=smoothHisto(histoM0T2);
		output[0]=histoM0T1;
		output[1]=histoM0T2;
		output[2]=binMed;

		
		//Normalize it to 0.9		//Add time to the Y values
//		double maxValT1=VitimageUtils.max(output[0]);
		//		double maxValT2=VitimageUtils.max(output[1]);
		//for(int i=0;i<numberBins;i++)output[0][i]=time+output[0][i]/maxValT1;
		//for(int i=0;i<numberBins;i++)output[1][i]=time+output[1][i]/maxValT2;
		return output;
	}
		
	public double[]smoothHisto(double[]histo){
		//Test avec sigma=3 bins;
		int sig=7;
		int nbBins=histo.length;
		double[]histSmooth=new double[nbBins];
		for(int i=sig;i<nbBins-sig;i++)histSmooth[i]=0.3*histo[i]+0.26*histo[i-1]+0.26*histo[i+1]+0.23*histo[i-2]+0.23*histo[i+2]+0.15*histo[i+3]+0.15*histo[i-3]+0.08*histo[i+4]+0.08*histo[i-4]+0.04*histo[i+5]+0.04*histo[i-5]+0.02*histo[i+6]+0.02*histo[i-6]+0.01*histo[i+7]+0.01*histo[i-7];
		return histSmooth;  
	}
	
	public void actualizeDisplayedNumbers() {
		String infoT1 = String.format("  M0=%5.3f  |  T1=%6.1f ms  |  Khi=%5.4f  p=%2.0f%c",
				paramsTimelapseT1[this.tCor-1][0],paramsTimelapseT1[this.tCor-1][1],this.khi2T1[this.tCor-1],this.pValT1[this.tCor-1],'%');
		info1.setText(infoT1);
		String infoT2 = String.format("  M0=%5.3f | M0Bi=%5.3f , %5.3f ||  T2=%5.1f ms | T2Bi=%5.1f ms , %5.1f ms || K1=%4.3f p1=%2.0f%c | K2=%4.3f p2=%2.0f%c",
				paramsTimelapseT2[this.tCor-1][0],paramsTimelapseT2[this.tCor-1][2],paramsTimelapseT2[this.tCor-1][4],
				paramsTimelapseT2[this.tCor-1][1],paramsTimelapseT2[this.tCor-1][3],(paramsTimelapseT2[this.tCor-1][5]>10000 ? 9999 : paramsTimelapseT2[this.tCor-1][5]),
				this.khi2T2Mono[this.tCor-1],this.pValT2Mono[this.tCor-1],'%',this.khi2T2Bicomp[this.tCor-1],this.pValT2Bicomp[this.tCor-1],'%');
		info2.setText(infoT2);
		String infoT3 = String.format("  T="+(this.tCor)+"/"+this.numberTimes+"  | Z="+(this.zCor)+"/"+this.numberZ+" | Area=%2dx%2dx%2d-%s | #Points=%d . Running computation = %d%c",
				1+2*this.crossWidth,1+2*this.crossWidth,1+2*this.crossThick,
				(this.gaussianWeighting?"" : ""),this.nPtsCur,this.computationRate,'%');
		info3.setText(infoT3);
	}

	public void computeFitsAgain() {
		actualizeEstimations2();
		actualizeExportSentence();
	}

	public void displayFitsAgain() {
		actualizeFirstPlots();
		actualizeSecondPlots();
		actualizeDisplayedNumbers();
	}
	






	/** Update graphs using new computed results*/		
	public void actualizeFirstPlots() {
		//HANDLE T1 CURVE
		Color crossColor=new Color (60,60,255);
		Color t1Moy=new Color (255,50,50);
		String []tabNamesT1=new String[]{"-Fit with noise estimation"};
		int incrT1=0;
		int deltaXt2=3;
		if(this.autoSizingTimePlots) {
			maxPlotYT1=VitimageUtils.max(dataTimelapseT1[tCor-1])*1.3;
			maxPlotYT2=VitimageUtils.max(dataTimelapseT2[tCor-1])*1.3;
		}
	
        //NOISE
        plotT1.setLineWidth(1);
		double[]statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(tabSigmasT1Seq[tCor-1][zCor-1]);
		double nSum=statusRoi==2 ? Math.sqrt(this.nPtsRoi) : (1+crossWidth);

		this.meanNoiseT1Cur=statsNoise[0];
		this.sigmaNoiseT1Cur=statsNoise[1];
		plotT1.setColor(new Color(210,210 ,210) );
		plotT1.setColor(new Color(0,0 ,0) );
        plotT1.setLineWidth(2);
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1},new double[]{statsNoise[0],statsNoise[0]});//Afficher le mean+sigma
        plotT1.setLineWidth(1);
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1},new double[]{statsNoise[0]+statsNoise[1]/nSum,statsNoise[0]+statsNoise[1]/nSum});//Afficher le mean+sigma
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1},new double[]{statsNoise[0]-statsNoise[1]/nSum,statsNoise[0]-statsNoise[1]/nSum});//Afficher le mean-sigma

		
        //OBSERVATIONS IRM
		plotT1.setLineWidth(3+thickCurve);
		plotT1.setColor(crossColor );
		plotT1.replace(incrT1++,"x",t1Times,dataTimelapseT1[tCor-1]);//Afficher les points IRM
        plotT1.setLineWidth(2);
        

        //FIT T1
        if(jitterT1[tCor-1]>=99)        		plotT1.setColor(Color.gray);
        else plotT1.setColor(t1Moy);
        
        plotT1.setLineWidth(1+thickCurve);
        plotT1.replace(incrT1++, "line",timesT1Cute,tabFittenT1Cute[tCor-1]);//Afficher la courbe 
        
        //PARAM T1 MONO
        if(jitterT1[tCor-1]>=99)        		plotT2.setColor(Color.gray);
        else         plotT1.setColor(t1Moy);
        plotT1.setLineWidth(3);
		plotT1.replace(incrT1++, "line",new double[]{0,(double)(paramsTimelapseT1[tCor-1][1])},new double[]{maxPlotYT1*0.9,maxPlotYT1*0.9});//Afficher le T1
        plotT1.setColor(crossColor);
        plotT1.setLineWidth(1);
		for(int t=0;t<this.nTr;t++)plotT1.replace(incrT1++, "line",new double[]{this.t1Times[t]-deltaXt2,this.t1Times[t]-deltaXt2},new double[]{dataTimelapseT1[tCor-1][t]-dataTimelapseT1Sigmas[tCor-1][t],dataTimelapseT1[tCor-1][t]+dataTimelapseT1Sigmas[tCor-1][t]});//Afficher le T2

        //LEGENDE T1
        plotT1.setLineWidth(1);
		String strLegendT1="";
		for(int hh=0;hh<1;hh++) strLegendT1+="\n";
		strLegendT1+="Noise +-sigma\n\n"+"MRI T1 relaxation"+"\n"+tabNamesT1[0];
		plotT1.setLimits(0, maxT1, 0,maxPlotYT1);
		plotT1.setColor(new Color(150,150 ,150) );
		plotT1.addLegend(strLegendT1,"bottom-right");
		plotCan1.setPlot(plotT1);	

		
		
		
		
	
		//HANDLE T2 CURVE
		Color t2Moy=new Color (255,50,50);
		Color t2Bi=new Color (0,175,0);
		String[] tabNamesT2=new String[]{"One T2","Two T2","Three T2"};
		int incrT2=0;

		//NOISE LEVEL
		plotT2.setLineWidth(1);
		statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(tabSigmasT2Seq[tCor-1][zCor-1]);
		this.meanNoiseT2Cur=statsNoise[0];
		this.sigmaNoiseT2Cur=statsNoise[1];
		plotT2.setColor(new Color(210,210 ,210) );
		nSum=Math.sqrt(this.nPtsCur);
		plotT2.setColor(new Color(0,0 ,0) );
		plotT2.setLineWidth(2);
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2},new double[]{statsNoise[0],statsNoise[0]});//Afficher le mean+sigma
		plotT2.setLineWidth(1);
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2},new double[]{statsNoise[0]+statsNoise[1]/nSum,statsNoise[0]+statsNoise[1]/nSum});//Afficher le mean+sigma
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2},new double[]{statsNoise[0]-statsNoise[1]/nSum,statsNoise[0]-statsNoise[1]/nSum});//Afficher le mean-sigma

        //DONNEES IRM
		plotT2.setLineWidth(3+thickCurve);
		plotT2.setColor(crossColor);
		plotT2.replace(incrT2++,"x",t2Times,dataTimelapseT2[tCor-1]);//Afficher les points IRM	//	plotT2.replace(incrT2++, "line",new double[]{0,0},new double[]{0,0});
        plotT2.setLineWidth(2);
        
        //COURBES FITTED
        if(jitterT2Mono[tCor-1]>=99)        		plotT2.setColor(Color.gray);
        else plotT2.setColor(t2Moy);
        plotT2.setLineWidth(1+thickCurve);
        plotT2.replace(incrT2++, "line",timesT2Cute,tabFittenT2MonoCute[tCor-1]);//Afficher la courbe monocomp
        if(jitterT2Bicomp[tCor-1]>=99)        		plotT2.setColor(Color.gray);
        else plotT2.setColor(t2Bi);
        plotT2.setLineWidth(1+thickCurve);
        plotT2.replace(incrT2++, "line",timesT2Cute,tabFittenT2BicompCute[tCor-1]);//Afficher la courbe bicomp
        
        //PARAM T2 MONO ET STD MONO
        if(jitterT2Mono[tCor-1]>=99)        		plotT2.setColor(Color.gray);
        else         plotT2.setColor(t2Moy);
        plotT2.setLineWidth(3);
		plotT2.replace(incrT2++, "line",new double[]{0,(double)(paramsTimelapseT2[tCor-1][1])},new double[]{maxPlotYT2*0.9,maxPlotYT2*0.9});//Afficher le T2
        plotT2.setLineWidth(1);
        plotT2.setColor(crossColor);
		for(int t=0;t<this.nTe;t++)plotT2.replace(incrT2++, "line",new double[]{this.t2Times[t],this.t2Times[t]},new double[]{dataTimelapseT2[tCor-1][t]-dataTimelapseT2SigmasMono[tCor-1][t],dataTimelapseT2[tCor-1][t]+dataTimelapseT2SigmasMono[tCor-1][t]});//Afficher le T2

		
        //PARAM T2 BI
        if(jitterT2Bicomp[tCor-1]>=99)        		plotT2.setColor(Color.gray);
        else         plotT2.setColor(t2Bi);
        plotT2.setLineWidth(3);
		plotT2.replace(incrT2++, "line",new double[]{0,paramsTimelapseT2[tCor-1][3]},new double[]{maxPlotYT2*0.83,maxPlotYT2*0.83});//Afficher le T2

        if(jitterT2Bicomp[tCor-1]>=99)        		plotT2.setColor(Color.gray);
        else         plotT2.setColor(t2Bi);
        plotT2.setLineWidth(3);
		plotT2.replace(incrT2++, "line",new double[]{0,paramsTimelapseT2[tCor-1][5]},new double[]{maxPlotYT2*0.79,maxPlotYT2*0.79});//Afficher le T2
        plotT2.setLineWidth(1);
 
        //LEGENDE
		plotT2.setLineWidth(1);
		String strLegendT2="";
		for(int hh=0;hh<1;hh++) strLegendT2+="\n";
		strLegendT2+="\nNoise +-sigma\nMRI T2 relaxation"+"\n"+tabNamesT2[0]+"\n"+tabNamesT2[1];
		if(this.autoSizingTimePlots) {
			maxPlotYT2=VitimageUtils.max(dataTimelapseT2[tCor-1]);
			maxPlotYT2*=1.3;
		}
		plotT2.setLimits(0, maxT2, 0,maxPlotYT2);
		plotT2.setColor(new Color(150,150 ,150) );
		plotT2.addLegend(strLegendT2,"bottom-right");
		plotCan2.setPlot(plotT2);
	}
		
	public void actualizeSecondPlots() {
		//washTheCurves
		for(int c=0;c<maxCurves;c++) {
			plotT21.replace(c, "line",new double[] {-200,-199},new double[] {-100,-100});
			plotT22.replace(c, "line",new double[] {-200,-199},new double[] {-100,-100});
		}
		double [][]compColMin=new double[][] {{255,0,0},{0,215,0},{0,0,255}};
		int [][]compColMinCentral=new int[][] {{90,0,0},{0,90,0},{0,0,90}};
		
		int incrT1=0;
		int incrT2=0;
		double y0T1;		double y1T1;				double x0T2;				double y0T2;		double y1T2;	double x1T1; 

		
        //HORIZONTAL GRID T21 ET T22
   		for(int tim=1;tim<=this.numberTimes;tim++) {
			plotT21.setLineWidth(2);
			plotT21.setColor(new Color(70,70,70));
			plotT21.replace(incrT1++, "line",new double[] {t0T2,t1T1},new double[] {tim,tim});
			plotT22.setLineWidth(2);
			plotT22.setColor(new Color(70,70,70));
			plotT22.replace(incrT2++, "line",new double[] {t0T2,t1T2},new double[] {tim,tim});
   		}

        //VERTICAL GRID T21
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT21.setLineWidth(1);        
			plotT21.setColor(new Color(150,150,150));
        	for(int mul=2;mul<=9;mul++)plotT21.replace(incrT1++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.numberTimes});
	        plotT21.setLineWidth(1);        
			plotT21.setColor(new Color(250,250,250));
        	plotT21.replace(incrT1++, "line",new double[] {tt,tt},new double[] {0,this.numberTimes});
        	plotT21.replace(incrT1++, "line",new double[] {tt*3,tt*3},new double[] {0,this.numberTimes});
        }

        //VERTICAL GRID T22
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT22.setLineWidth(1);        
			plotT22.setColor(new Color(150,150,150));
			for(int mul=2;mul<=9;mul++)plotT22.replace(incrT2++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.numberTimes});
	        plotT22.setLineWidth(1);        
			plotT22.setColor(new Color(250,250,250));
        	plotT22.replace(incrT2++, "line",new double[] {tt,tt},new double[] {0,this.numberTimes});
        	plotT22.replace(incrT2++, "line",new double[] {tt*3,tt*3},new double[] {0,this.numberTimes});
        }
        //Draw rectangle of T position
				x1T1=t1T1*0.95;
		y0T1=tCor-0.95;		y1T1=tCor-0.05;		
		x0T2=t0T2*1.05;		
		y0T2=tCor-0.95;		y1T2=tCor-0.05;		
		plotT21.setLineWidth(4);
		plotT21.setColor(new Color(220,220,0));
		plotT21.replace(incrT1++, "line",new double[] {x0T2,x1T1},new double[] {y0T1,y0T1});
        plotT21.replace(incrT1++, "line",new double[] {x0T2,x1T1},new double[] {y1T1,y1T1});
        plotT21.replace(incrT1++, "line",new double[] {x0T2,x0T2},new double[] {y0T1,y1T1});
        plotT21.replace(incrT1++, "line",new double[] {x1T1,x1T1},new double[] {y0T1,y1T1});

        plotT22.setLineWidth(4);        
		plotT22.setColor(new Color(220,220,0));
        plotT22.replace(incrT2++, "line",new double[] {x0T2,x1T1},new double[] {y0T2,y0T2});
        plotT22.replace(incrT2++, "line",new double[] {x0T2,x1T1},new double[] {y1T2,y1T2});
        plotT22.replace(incrT2++, "line",new double[] {x0T2,x0T2},new double[] {y0T2,y1T2});
        plotT22.replace(incrT2++, "line",new double[] {x1T1,x1T1},new double[] {y0T2,y1T2});
        
        
   		for(int tim=1;tim<=this.numberTimes;tim++) {
 			double maxValT1=standardCapillaryM0;
			double maxValT2Mono=standardCapillaryM0;
	 		double quota;
	 		double radius=0;
	 		boolean flagGray;
	 		int rad;
	 		double factRad=0.005;

			compColMin=new double[][] {{255,0,0},{0,215,0},{0,0,255}};
			compColMinCentral=new int[][] {{0,220,0},{220,0,0},{0,0,90}};	 		
	 		quota=1;flagGray=false;
	 		radius=0.05+0.5*paramsTimelapseT1[tim-1][0]/maxValT1;if(radius>1)radius=1;if(radius<0.05)radius=0.05;	 		rad=1+(int)(Math.round((radius-0.05)/0.015));plotT21.setLineWidth(rad);
	 		if(jitterT1[tim-1]>30) flagGray=true;
			if(!flagGray)plotT21.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
	 		else plotT21.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));
	 		plotT21.setLineWidth((rad));
			plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim-1][1],paramsTimelapseT1[tim-1][1]},new double[] {tim-1+0.02+0.01*rad,tim-0.01*rad});//Afficher la courbe
	 		plotT21.setLineWidth(1);
	 		plotT21.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
			plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim-1][1],paramsTimelapseT1[tim-1][1]},new double[] {tim-1-0.05,tim-1+1.05});//Afficher la courbe
			if(this.nRepetMonteCarlo>0) {
				plotT21.setLineWidth(2);
		 		plotT21.setColor(new Color(255,255,255));
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim-1][1]-paramsTimelapseT1Sigmas[tim-1][1],paramsTimelapseT1[tim-1][1]+paramsTimelapseT1Sigmas[tim-1][1]},new double[] {tim-0.5,tim-0.5});//Afficher la courbe
			}
			
			compColMin=new double[][] {{0,215,0},{255,0,0},{0,0,255}};
			compColMinCentral=new int[][] {{220,0,0},{90,0,0},{0,0,90}};
	 		quota=1;flagGray=false;
			if(khi2T2Mono[tim-1]<khi2T2Bicomp[tim-1]) {
		 		if(jitterT2Mono[tim-1]>30) flagGray=true;
				if(!flagGray)plotT21.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
		 		else plotT21.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));

				radius=0.05+0.5*paramsTimelapseT2[tim-1][0]/maxValT2Mono;if(radius>1.0)radius=1.0;if(radius<0.05)radius=0.05;	 		rad=2+(int)(Math.round((radius-0.05)/0.015));plotT22.setLineWidth(rad);
		 		plotT21.setLineWidth((rad));
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][1],paramsTimelapseT2[tim-1][1]},new double[] {tim-1+0.1+Math.min(0.5,factRad*rad),tim-0.1-Math.min(0.5,factRad*rad)});//Afficher la courbe
		 		plotT21.setLineWidth(1);
		 		plotT21.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][1],paramsTimelapseT2[tim-1][1]},new double[] {tim-1-0.05,tim+0.05});//Afficher la courbe
				if(this.nRepetMonteCarlo>0) {
			 		plotT21.setLineWidth(2);
			 		plotT21.setColor(new Color(255,255,255));
					plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][1]-paramsTimelapseT2Sigmas[tim-1][1],paramsTimelapseT2[tim-1][1]+paramsTimelapseT2Sigmas[tim-1][1]},new double[] {tim-0.5,tim-0.5});//Afficher la courbe
				}
			}
			else {
		 		if(jitterT2Bicomp[tim-1]>30) flagGray=true;
				if(!flagGray)plotT21.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
		 		else plotT21.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));

				radius=0.05+0.5*paramsTimelapseT2[tim-1][2]/maxValT2Mono;if(radius>1.0)radius=1.0;if(radius<0.05)radius=0.05;	 		rad=2+(int)(Math.round((radius-0.05)/0.015));plotT22.setLineWidth(rad);
		 		plotT21.setLineWidth((rad));
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][3],paramsTimelapseT2[tim-1][3]},new double[] {tim-1+0.1+Math.min(0.5,factRad*rad),tim-0.1-Math.min(0.5,factRad*rad)});//Afficher la courbe
		 		plotT21.setLineWidth(1);
		 		plotT21.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][3],paramsTimelapseT2[tim-1][3]},new double[] {tim-1-0.05,tim+0.05});//Afficher la courbe
				if(this.nRepetMonteCarlo>0) {
			 		plotT21.setLineWidth(2);
			 		plotT21.setColor(new Color(255,255,255));
					plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][3]-paramsTimelapseT2Sigmas[tim-1][3],paramsTimelapseT2[tim-1][3]+paramsTimelapseT2Sigmas[tim-1][3]},new double[] {tim-0.5,tim-0.5});//Afficher la courbe
				}
		 		quota=1;flagGray=false;
		 		radius=0.05+0.5*paramsTimelapseT2[tim-1][4]/maxValT2Mono;if(radius>1.0)radius=1.0;if(radius<0.05)radius=0.05;	 		rad=2+(int)(Math.round((radius-0.05)/0.015));plotT22.setLineWidth(rad);
		 		if(jitterT2Bicomp[tim-1]>30) flagGray=true;
				if(!flagGray)plotT21.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
		 		else plotT21.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));
		 		plotT21.setLineWidth((rad));
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][5],paramsTimelapseT2[tim-1][5]},new double[] {tim-1+0.1+Math.min(0.5,factRad*rad),tim-0.1-Math.min(0.5,factRad*rad)});//Afficher la courbe
		 		plotT21.setLineWidth(1);
		 		plotT21.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][5],paramsTimelapseT2[tim-1][5]},new double[] {tim-1-0.05,tim+0.05});//Afficher la courbe
				if(this.nRepetMonteCarlo>0) {
			 		plotT21.setLineWidth(2);
			 		plotT21.setColor(new Color(255,255,255));
					plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim-1][5]-paramsTimelapseT2Sigmas[tim-1][5],paramsTimelapseT2[tim-1][5]+paramsTimelapseT2Sigmas[tim-1][5]},new double[] {tim-0.25,tim-0.25});//Afficher la courbe
				}
			}		
	 		plotT22.setLineWidth(2);
	 		plotT22.setColor(new Color(255,0,0));
	 		plotT22.replace(incrT2++, "line", valsSpectrum[tim-1][2], valsSpectrum[tim-1][0]);
	 		plotT22.setColor(new Color(0,255,0));
	 		plotT22.replace(incrT2++, "line", valsSpectrum[tim-1][2], valsSpectrum[tim-1][1]);
   			//if needed, draw ranging interval
	 		if(this.spectrumRangingMode) {
	 			plotT22.setColor(new Color(0,0,255));	 		
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[0]}, new double[] {tim-0.7,tim-0.3});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[2],this.rangingBoundaries[2]}, new double[] {tim-0.7,tim-0.3});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[1],this.rangingBoundaries[1]}, new double[] {tim-0.8,tim-0.2});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[2]}, new double[] {tim-0.5,tim-0.5});
	 		}   			
   		}				
		plotT21.setLimits(t0T2,t1T1, 0, numberTimes);
		plotT22.setLimits(t0T2,t1T1, 0, numberTimes);
	}				

	public void actualizeSpectrumCurves() {
		for(int tim=1;tim<=this.numberTimes;tim++) {
			this.valsSpectrum[tim-1]=getDisplayableSpectrumCurveMultiThread(tim-1);
		}
		if(!this.separateNormalizationSpectrum) {
			for(int tim=1;tim<=this.numberTimes;tim++) {
				double maxValT1=VitimageUtils.max(this.valsSpectrum[tim-1][0]);
				double maxValT2=VitimageUtils.max(this.valsSpectrum[tim-1][1]);
				for(int i=0;i<this.valsSpectrum[tim-1][0].length;i++)this.valsSpectrum[tim-1][0][i]=tim-1+this.valsSpectrum[tim-1][0][i]/maxValT1;
				for(int i=0;i<this.valsSpectrum[tim-1][1].length;i++)this.valsSpectrum[tim-1][1][i]=tim-1+this.valsSpectrum[tim-1][1][i]/maxValT2;
			}
		}
		else {
			double maxTotT1=0;
			double maxTotT2=0;
			for(int tim=1;tim<=this.numberTimes;tim++) {
				double maxValT1=VitimageUtils.max(this.valsSpectrum[tim-1][0]);
				double maxValT2=VitimageUtils.max(this.valsSpectrum[tim-1][1]);
				if(maxTotT1<maxValT1)maxTotT1=maxValT1;
				if(maxTotT2<maxValT2)maxTotT2=maxValT2;
			}
			for(int tim=1;tim<=this.numberTimes;tim++) {
				for(int i=0;i<this.valsSpectrum[tim-1][0].length;i++)this.valsSpectrum[tim-1][0][i]=tim-1+this.valsSpectrum[tim-1][0][i]/maxTotT1;
				for(int i=0;i<this.valsSpectrum[tim-1][1].length;i++)this.valsSpectrum[tim-1][1][i]=tim-1+this.valsSpectrum[tim-1][1][i]/maxTotT2;
			}
		}
		actualizeRangedData();
	}

	public void actualizeExportSentence() {
		rSentence="\n\n#R CODE\n";
		rSentence+="# * t1TimesAndObservations : \n#      first line = recovery times used for observations\n#      following lines are actual echo spin observations (one line per time-point in case of timelapse)\n";
		rSentence+="#     Data are normalised : the values have been divided by the last echo value of the capillary in the central slice)\n";
		rSentence+="# * t2TimesAndObservations : \n#      first line = echo times used for observations\n#      following lines are actual echo spin observations (one line per time-point in case of timelapse)\n";
		rSentence+="#     Data are normalised : the values have been divided by the estimation of M0 of the capillary in the central slice)\n";
		rSentence+="# * estimatedParameters : \n#      each line gives estimated parameters (M,T1,T2) obtained for each time point (in case of timelapse)\n";
		rSentence+="#      all parameters are computed using exponential fit, after rice noise estimation and correction\n";
		rSentence+="#      Params in each line : [ ObservationTime  ,  M0(of T1 relax curve)  ,  T1 ,  M0(of T2 relax curve monocomponent) , \n#                                           T2(mono-exponential), M0-1(first part of bi-exponential) , T2-1 (first T2) ,  M0-2 (second part) , T2-2 (second T2) , \n#                                            Jitter T1 (% error in fit), Jitter T2 mono , Jitter T2 bicomp\n";
		rSentence+="# code begins here\n#\n";
		rSentence+="t1TimesAndObservations <- data.frame(acquisition_number=c(0,";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(t+1)+",";
		rSentence+=""+numberTimes+")";
		rSentence+=", acquisition_day=c(0,";
		for(int t=1;t<numberTimes-1;t++)rSentence+=""+this.actualDay[t]+",";
		rSentence+=""+this.actualDay[numberTimes-1]+")";
		for(int ec=0;ec<t1Times.length;ec++) {
			rSentence+=", Echo"+ec+"=c("+t1Times[ec]+",";
			for(int t=0;t<numberTimes-1;t++)rSentence+=""+(dataTimelapseT1[t][ec])+",";
			rSentence+=""+(dataTimelapseT1[numberTimes-1][ec])+")";
		}
		rSentence+=")\n";
		rSentence+="t2TimesAndObservations <- data.frame(acquisition_number=c(0,";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(t+1)+",";
		rSentence+=""+numberTimes+")";
		rSentence+=", acquisition_day=c(0,";
		for(int t=1;t<numberTimes-1;t++)rSentence+=""+this.actualDay[t]+",";
		rSentence+=""+this.actualDay[numberTimes-1]+")";
		for(int ec=0;ec<t2Times.length;ec++) {
			rSentence+=", Echo"+ec+"=c("+t2Times[ec]+",";
			for(int t=0;t<numberTimes-1;t++)rSentence+=""+(dataTimelapseT2[t][ec])+",";
			rSentence+=""+(dataTimelapseT2[numberTimes-1][ec])+")";
		}
		rSentence+=")\n";
	
	
		
		rSentence+="estimatedParameters <- data.frame(acquisition_number=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(t+1)+",";
		rSentence+=""+numberTimes+"),acquisition_day=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+actualDay[t]+",";
		rSentence+=""+actualDay[numberTimes-1]+"), M0T1=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT1[t][0])+",";
		rSentence+=""+""+(paramsTimelapseT1[numberTimes-1][0])+"), T1=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT1[t][1])+",";
		rSentence+=""+""+(paramsTimelapseT1[numberTimes-1][1])+"), M0T2=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT2[t][0])+",";
		rSentence+=""+""+(paramsTimelapseT2[numberTimes-1][0])+"), T2=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT2[t][1])+",";
		rSentence+=""+""+(paramsTimelapseT2[numberTimes-1][1])+"), M0T21=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT2[t][2])+",";
		rSentence+=""+""+(paramsTimelapseT2[numberTimes-1][2])+"), T21=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT2[t][3])+",";
		rSentence+=""+""+(paramsTimelapseT2[numberTimes-1][3])+"), M0T22=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT2[t][4])+",";
		rSentence+=""+""+(paramsTimelapseT2[numberTimes-1][4])+"), T22=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(paramsTimelapseT2[t][5])+",";
		rSentence+=""+""+(paramsTimelapseT2[numberTimes-1][5])+"), JitterT1=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(jitterT1[t])+",";
		rSentence+=""+""+(jitterT1[numberTimes-1])+"), JitterT2MonoExp=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(jitterT2Mono[t])+",";
		rSentence+=""+""+(jitterT2Mono[numberTimes-1])+"), JitterT2BiExp=c(";
		for(int t=0;t<numberTimes-1;t++)rSentence+=""+(jitterT2Bicomp[t])+",";
		rSentence+=""+""+(jitterT2Bicomp[numberTimes-1])+"))\n\n";
	
	
		
		pySentence="\n\n#CODE PYTHON\n";
		pySentence+="# * t1TimesAndObservations : \n#      first line = recovery times used for observations\n#      following lines are actual echo spin observations (one line per time-point in case of timelapse)\n";
		pySentence+="#     Data are normalised : the values have been divided by the last echo value of the capillary in the central slice)\n";
		pySentence+="# * t2TimesAndObservations : \n#      first line = echo times used for observations\n#      following lines are actual echo spin observations (one line per time-point in case of timelapse)\n";
		pySentence+="#     Data are normalised : the values have been divided by the estimation of M0 of the capillary in the central slice)\n";
		pySentence+="# * estimatedParameters : \n#      each line gives estimated parameters (M,T1,T2) obtained for each time point (in case of timelapse)\n";
		pySentence+="#      all parameters are computed using exponential fit, after rice noise estimation and correction\n";
		pySentence+="#      Params in each line : [ ObservationTime  ,  M0(of T1 relax curve)  ,  T1 ,  M0(of T2 relax curve monocomponent) , \n#                                  T2(mono-exponential), M0-1(first part of bi-exponential) , T2-1 (first T2) ,  M0-2 (second part) , T2-2 (second T2) , \n#                                   Jitter T1 (% error in fit), Jitter T2 mono , Jitter T2 bicomp\n";
		pySentence+="# code begins here\n#\nimport numpy as np\n";
		pySentence+="t1TimesAndObservations=np.array([[0,0";
		for(int ec=0;ec<t1Times.length;ec++)pySentence+=","+t1Times[ec];
		pySentence+="]";					
		for(int t=1;t<=numberTimes;t++) {
			pySentence+=",\n["+t+","+this.actualDay[t-1];
			for(int ec=0;ec<t1Times.length;ec++) pySentence+=","+dataTimelapseT1[t-1][ec];
			pySentence+="]";
		}
		pySentence+="\n])\n";
	
		pySentence+="t2TimesAndObservations=np.array([[0,0";
		for(int ec=0;ec<t2Times.length;ec++)pySentence+=","+t2Times[ec];
		pySentence+="]";					
		for(int t=1;t<=numberTimes;t++) {
			pySentence+=",\n["+t+","+actualDay[t-1];
			for(int ec=0;ec<t2Times.length;ec++) pySentence+=","+dataTimelapseT2[t-1][ec];
			pySentence+="]";
		}
		pySentence+="\n])\n";
	
		
		pySentence+="estimatedParameters=np.array([\n";
		for(int t=0;t<numberTimes;t++) {
			pySentence+="["+(t+1)+","+actualDay[t];		
			for(int val=0;val<2;val++)pySentence+=","+paramsTimelapseT1[t][val];
			for(int val=0;val<6;val++)pySentence+=","+paramsTimelapseT2[t][val];
			pySentence+=","+jitterT1[t]+","+jitterT2Mono[t]+","+jitterT2Bicomp[t];
			pySentence+="]"+(t==numberTimes-1 ? "\n" : ",\n");
		}
		pySentence+="])\n";
	
		matSentence="\n\n%MATLAB\n";
		matSentence+="% * t1TimesAndObservations : \n%      first line = recovery times used for observations\n%      following lines are actual echo spin observations (one line per time-point in case of timelapse)\n";
		matSentence+="%     Data are normalised : the values have been divided by the last echo value of the capillary in the central slice)\n";
		matSentence+="% * t2TimesAndObservations : \n%      first line = echo times used for observations\n%      following lines are actual echo spin observations (one line per time-point in case of timelapse)\n";
		matSentence+="%     Data are normalised : the values have been divided by the estimation of M0 of the capillary in the central slice)\n";
		matSentence+="% * estimatedParameters : \n%      each line gives estimated parameters (M,T1,T2) obtained for each time point (in case of timelapse)\n";
		matSentence+="%      all parameters are computed using exponential fit, after rice noise estimation and correction\n";
		matSentence+="%      Params in each line : [ ObservationTime  ,  M0(of T1 relax curve)  ,  T1 ,  M0(of T2 relax curve monocomponent) , \n%                                           T2(mono-exponential), M0-1(first part of bi-exponential) , T2-1 (first T2) ,  M0-2 (second part) , T2-2 (second T2) , \n%                                            Jitter T1 (% error in fit), Jitter T2 mono , Jitter T2 bicomp\n";
		matSentence+="% code begins here\n%\n";
	
		
		matSentence+="t1TimesAndObservations=[\n[0,0";
		for(int ec=0;ec<t1Times.length;ec++)matSentence+=","+t1Times[ec];
		matSentence+="]";					
		for(int t=1;t<=numberTimes;t++) {
			matSentence+=";\n["+t+","+actualDay[t-1];
			for(int ec=0;ec<t1Times.length;ec++) matSentence+=","+dataTimelapseT1[t-1][ec];
			matSentence+="]";
		}
		matSentence+="\n]\n";
	
		matSentence+="t2TimesAndObservations=[\n[0,0";
		for(int ec=0;ec<t2Times.length;ec++)matSentence+=","+t2Times[ec];
		matSentence+="]";					
		for(int t=1;t<=numberTimes;t++) {
			matSentence+=";\n["+t+","+actualDay[t-1];
			for(int ec=0;ec<t2Times.length;ec++) matSentence+=","+dataTimelapseT2[t-1][ec];
			matSentence+="]";
		}
		matSentence+="\n]\n";
	
		
		matSentence+="estimatedParameters=[\n";
		for(int t=0;t<numberTimes;t++) {
			matSentence+="["+(t+1)+","+actualDay[t];		
			for(int val=0;val<2;val++)matSentence+=","+paramsTimelapseT1[t][val];
			for(int val=0;val<6;val++)matSentence+=","+paramsTimelapseT2[t][val];
			matSentence+=","+jitterT1[t]+","+jitterT2Mono[t]+","+jitterT2Bicomp[t];
			matSentence+="]"+(t==numberTimes-1 ? "\n" : ";\n");
		}
		matSentence+="]\n";
	}
	
	
		
		
	
	
	

	
	
	
	
	
	
	//Images used and their parameters
	ImagePlus fullHyperImage;
	ImagePlus T1seq;
	ImagePlus T2seq;
	double[]voxS;
	int[]dims;
	int nTe=0;
	int nTr=0;
	double [] t1Times;
	double [] t2Times;
	int[]actualDay;


	//Parameters of graph and distribution estimation
	boolean rangeDisplayedFullBlocks=true;
	public double[]rangingBoundaries=new double[3];
	public double rangingFactor=0.2;
	public boolean spectrumRangingMode=false;
	public int gaussianSpectrum=2;
	public boolean separateNormalizationSpectrum=true;
	boolean multiThreaded=true;
	boolean sigmaKhi2WeightedByNumberPoints=false;
	RiceEstimator riceEstimator;
	double[][][]sigmaDataRiceLookupTableT1;
	double[][][]sigmaDataRiceLookupTableT2;
	private int fitAlgorithm=MRUtils.SIMPLEX;
	int computationRate=100;
	int curEc1;
	int curEc2;
	int nRepetMonteCarlo=0;
	boolean gaussianWeighting=false;
	int maxCurves=170;
	double meanNoiseT1Cur=1;
	double sigmaNoiseT1Cur=1;
	double meanNoiseT2Cur=1;
	double sigmaNoiseT2Cur=1;
	int nPtsCur=1;
	int nPtsRoi=1;
	int thickCurve=1;

	
	
	//Gui parameters : Window constants
	public int xMouseRange=0;
	private int TARGET_HEIGHT_FOR_PLOTS_AND_IMGS=512;
	private int TARGET_WIDTH_FOR_PLOTS_AND_IMGS=512;
	private int DX_IMG=TARGET_WIDTH_FOR_PLOTS_AND_IMGS;
	private int DY_IMG=TARGET_HEIGHT_FOR_PLOTS_AND_IMGS;
	private int DY_TEXT=30;
	private int DELTA_X=10;
	private int DELTA_Y=10;
	private int DX_PLOT_1=350;
	private int DX_PLOT_2=350;
	private int DY_PLOT_1=(2*TARGET_HEIGHT_FOR_PLOTS_AND_IMGS-3*DY_TEXT-4*+DELTA_Y)/2;
	private int DY_PLOT_2=(2*TARGET_HEIGHT_FOR_PLOTS_AND_IMGS-3*DY_TEXT-4*+DELTA_Y)/2;
	private int totalSizeX= DX_IMG+DX_PLOT_1+DX_PLOT_2+2*DELTA_X;  
	private int totalSizeY= 2*DY_IMG;  
	boolean autoSizingTimePlots=true;
	boolean isFireLut=false;
	

	//Params for Gui : Objects and texts to display
	int targetNumberCharsInTexts=50;
	boolean sizeBasedOnImage=false;
	boolean debugDisplay=true;
	Roi userRoi;
	String rSentence="";
	String pySentence="";
	String matSentence="";
	int crossWidth=0;
	int crossThick=0;
	int numberTimes=0;
	int numberZ=0;
	int xMouse=0;
	int yMouse=0;
	int xCor=1;
	int yCor=1;
	int tCor=1;
	int zCor=1;
	private Plot plotT1;
	private Plot plotT2;
	private Plot plotT21;
	private Plot plotT22;
	private ImageCanvas imgCan1;
	private ImageCanvas imgCan2;
	private PlotCanvas plotCan1;
	private PlotCanvas plotCan2;
	private PlotCanvas plotCan21;
	private PlotCanvas plotCan22;
	private double zoomLevel=0;
	private double maxPlotYT1= 1.0;
	private double maxPlotYT2= 0.58;
	private int currentCanvas=1;
	TextField info1;
	TextField info2;
	TextField info3;
	static PlugInFrame instance;
	boolean computeMultiComp=true;
	private static final long serialVersionUID = 1L;
	double xStepCuteT1=50;
	double xStepCuteT2=3;
	double []timesT1Cute;
	double []timesT2Cute;

	double[][]tabSigmasT1Seq;
	double[][]tabSigmasT2Seq;

	double[][]dataTimelapseT1;
	double[][]dataTimelapseT2;
	double[][][]dataTimelapseT1Full;
	double[][][]dataTimelapseT2Full;
	double[][]dataTimelapseT1Sigmas;
	double[][]dataTimelapseT2SigmasMono;
	double[][]dataTimelapseT2SigmasBicomp;

	double[][]paramsTimelapseT1;
	double[][]paramsTimelapseT2;
	double[][]paramsTimelapseT1Sigmas;
	double[][]paramsTimelapseT2Sigmas;

	double [][]tabFittenT1;
	double [][]tabFittenT2Mono;
	double [][]tabFittenT2Bicomp;
	double [][]tabFittenT2Tricomp;
	int [][][]pointsCoords;
	double [][]tabFittenT1Cute;
	double [][]tabFittenT2MonoCute;
	double [][]tabFittenT2BicompCute;
	double [][]tabFittenT2TricompCute;

	double[][][]valsSpectrum;
	double[]jitterT1;
	double[]jitterT2Mono;
	double[]jitterT2Bicomp;
	double[]jitterT2Tricomp;

	double[]khi2T1;
	double[]khi2T2Mono;
	double[]khi2T2Bicomp;
	double[]khi2T2Tricomp;
	double[]pValT1;
	double[]pValT2Mono;
	double[]pValT2Bicomp;
	double[]pValT2Tricomp;

	double statusRoi;
	int[][]correspondanceCanvas;
	private ArrayList<int[]> rangeRoiPoints;


	@Override
	public void actionPerformed(ActionEvent e) {
	}

	public String getHelpString() {
		return (
		"Curve explorer is a java plugin dedicated to T1 and T2 caracterization in plant tissue from packed hyperimage 5D (X-Y-Z-Time-Echo).\n"+
		" \n"+
		" \n"+
		"* Click on T1 / T2 image to set the center of the area of interest. \n"+
		"     - The area center is showed as a cross.\n"+
		"     - This area is a 3D cube, centered around the cross. Its size is ( (1+2*radiusXY) x (1+2*radiusXY) x (1+2*radiusZ))\n"+
		"     - Data from the whole cube is used for T1 / T2 estimation : the MRI spin-echo values are averaged over this area. \n"+
		"     - Use '+' or '-' to zoom in / out\n"+
		" \n\n"+
		"* The data can also be gathered from a user-defined Roi.\n"+
		"     - hit 'r' to enter in roi (region of interest) mode.\n"+
		"     - Select your roi using the Opener and load your custom roi._\n"+
		"	  - Hit 'r' again to quit the roi mode and bring back the cross\n"+
		" \n\n"+					
		"* The plots display the current data from the area of interest, and the estimated exponential curves and parameters\n"+
		"     - The parameters are computed using an exponential fit with correction of the rice noise (sigma rice estimation is handled automatically)\n"+
		"     - Default fit optimization use Simplex. Hit 'l' to toggle to / from the Levenberg-Marquard optimization\n"+
		"     - White plots :\n"+
		"         .White left is T1 recovery data from the successive recovery time, and the mono-exponential estimated T1 curve\n"+
		"         .White right is T2 relaxation plot from the successive echoes, with mono-exponential and bi-exponential estimated T2 curve\n"+
		"         .Use 't' to switch between thick-curves-mode and thin-curves-mode to check visually the estimation accuracy\n"+
		"     - Black plots :\n"+
		"         .Black left is the estimated T1 time from the successive timepoints (if any).\n"+
		"         .Black right is the estimated T2 time(s) in mono- and bi-exponential modes from the successive timepoints (if any)\n"+
		"         .The thicker a ray, the greater is the estimated Proton density.\n"+
		"     	  .If estimation parameters confidence is low (jitter > 30 %), rays are shown in grey\n"+
		"  \n\n"+
		" \n"+
		"* Use 'a','q','z','s' 'g' keys to change the size and type of interest area) : \n"+
		"     - 'a' -> remove 1 to radiusXY (area becoming smaller, the MRI signal is noisier and T1/T2 estimation is less significant) \n"+
		"     - 'q' -> add 1 to radiusXY (area becoming bigger, many various tissues are merged during estimation) \n"+
		"     - 'z' -> remove 1 to radiusZ (area becoming smaller, the MRI signal is noisier and T1/T2 estimation is less significant) \n"+
		"     - 's' -> add 1 to radiusZ (area becoming bigger, many various tissues are merged during estimation) \n"+
		"     - 'g' -> Switch to/from uniform averaging to/from gaussian averaging. With gaussian, the central point have more ponderation than the ones in the border\n"+
		" \n"+
		"* Use '4','6','8','2' keys to change current observation time or depth in the stack : \n"+
		"     - '4' -> Change to previous slice (Z-)\n"+
		"     - '6' -> Change to next slice (Z+)\n"+
		"     - '8' -> Change to next time (T+)\n"+
		"     - '2' -> Change to previous time (T+)\n"+
		" \n"+					
		"* Export data and computed parameters : hit 'e' to export current data in matlab/octave/python/r format in the ImageJ log window\n"+
		" \n"+
		" \n"+
		"* The top text fields display continuously precious data :\n"+
		"     - First bar :\n"+
		"	      . M0 : Proton density estimated from the T1 recovery curve\n"+
		"	      . T1 : T1 time estimated from the T1 recovery curve\n"+
		"	      . Jitter : Root mean square error between MRI data and estimated curve, in % of the mean MRI data \n"+
		"	      . Area : dX-dY-dZ-Type display the dimensions of the interest area around the cross, and the type of averaging (Gauss or Uni)\n"+
		"	      . #Points : displays the number of points used for averaging, in cross mode or in roi mode\n"+ 
		"     - Second bar :\n"+
		"	      . M0 : Proton density estimated from the T2 relaxation mono-exponential\n"+
		"	      . M01 : Proton density associated with the first component of the bi-exponential fit\n"+
		"	      . M02 : Proton density associated with the second component of the bi-exponential fit\n"+
		"	      . T2 : T2 time estimated from the T2 relaxation mono-exponential\n"+
		"	      . T21 : T2 time associated with the first component of the bi-exponential fit\n"+
		"	      . T22 : T2 time associated with the second component of the bi-exponential fit\n"+
		"	      . JitMono : Root mean square error between MRI data and estimated mono-exponential curve, in % of the mean MRI data \n"+
		"	      . JitBi : Root mean square error between MRI data and estimated bi-exponential curve, in % of the mean MRI data \n"+					
		" \n");
	}

	
	
	
	public void mouseClicked(MouseEvent e) {
		if(imgCan1.cursorOverImage()) currentCanvas=WIN_T1;
		if(imgCan2.cursorOverImage()) currentCanvas=WIN_T2;
		if(plotCan1.cursorOverImage()) currentCanvas=WIN_PLOT1;
		if(plotCan2.cursorOverImage()) currentCanvas=WIN_PLOT2;
		if(plotCan21.cursorOverImage()) currentCanvas=WIN_PLOT21;
		if(plotCan22.cursorOverImage()) currentCanvas=WIN_PLOT22;

		if(currentCanvas==WIN_T1 || currentCanvas==WIN_T2) {
			IJ.log("Click on MR images | Coordinates=("+ xMouse+","+yMouse+")"+"  |  zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			xMouse=e.getX();
			yMouse=e.getY();
			actualizeCursor();
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			actualizeSpectrumCurves();
			this.actualizeRangedData();
			computeFitsAgain();

			displayFitsAgain();
			actualizeCursor();
			actualizeRangingBoundaries();
			actualizeRangedData();
			this.actualizeRangingBoundaries();
			this.actualizeRangedData();
		
			displayFitsAgain();
			actualizeCursor();
			plotCan1.repaint();
			plotCan21.repaint();
			imgCan1.repaint();
			plotCan2.repaint();
			plotCan22.repaint();
			imgCan2.repaint();
			repaint();
			return;
		}

		else if(currentCanvas==WIN_PLOT22) {
			IJ.log("Click on Plot22 | Coordinates=("+ xMouse+","+yMouse+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			spectrumRangingMode=true;	
			if(this.spectrumRangingMode) {
				this.xMouseRange=e.getX();
				actualizeRangingBoundaries();
				actualizeRangedData();
				this.actualizeRangingBoundaries();
				this.actualizeRangedData();
			
				displayFitsAgain();
				actualizeCursor();
				plotCan1.repaint();
				plotCan21.repaint();
				imgCan1.repaint();
				plotCan2.repaint();
				plotCan22.repaint();
				imgCan2.repaint();
				repaint();
				return;
			}
		}
		else {
			IJ.log("Click on another plot. Ranging mode stop");
			this.spectrumRangingMode=false;	
			actualizeCursor();
			actualizeMriObservationsBasedOnData();
			actualizeSpectrumCurves();
			computeFitsAgain();
			displayFitsAgain();
			repaint();
			plotCan1.repaint();
			plotCan21.repaint();
			imgCan1.repaint();
			plotCan2.repaint();
			plotCan22.repaint();
			imgCan2.repaint();
			return;
		}

	}
	
	@Override
	public void mousePressed(MouseEvent e) {
	}

	@Override
	public void mouseReleased(MouseEvent e) {
	}
	
	public void mouseDragged(MouseEvent e) {
	}

	@Override
	public void mouseEntered(MouseEvent e) {
	}

	@Override
	public void mouseExited(MouseEvent e) {
	}
		
	@Override
	public void keyPressed(KeyEvent e) {
	}

	@Override
	public void keyReleased(KeyEvent e) {
	}

	public void windowClosing(WindowEvent paramWindowEvent)
	  {
	    super.windowClosing(paramWindowEvent);
	    instance = null;
	  }
	

	
	
	private static final double standardCapillaryM0=5000;
	private static final double maxAcceptableM0T2=2*standardCapillaryM0;
	private static final double maxAcceptableM0T1=5*standardCapillaryM0;
	private static final double minAcceptableT2=8;
	private static final double minAcceptableT1=10;
	private static final double maxAcceptableT2=5000;
	private static final double maxT1=10500;
	private static final double maxT2=350;
	private static final double maxDisplayedT1=4000;
	private static final double maxDisplayedT2=120;
	private static final int t0T1=100;
	private static final int t1T1=10000;
	private static final int t0T2=10;
	private static final int t1T2=1000;
	private static final int WIN_T1=1;
	private static final int WIN_T2=2;
	private static final int WIN_PLOT1=3;
	private static final int WIN_PLOT2=4;
	private static final int WIN_PLOT21=6;
	private static final int WIN_PLOT22=7;
	
	
}

