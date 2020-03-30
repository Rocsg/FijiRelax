package com.vitimage.aplimtools;


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

import com.vitimage.common.Timer;
import javax.swing.JPanel;

import com.vitimage.common.TransformUtils;
import com.vitimage.common.VitiDialogs;
import com.vitimage.common.VitimageUtils;
import com.vitimage.mrutils.HyperMRIT1T2;
import com.vitimage.mrutils.MRUtils;
import com.vitimage.mrutils.RiceEstimator;

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
 *(priority=0/3, difficulty=3/3) Introduce multiple T1 estimation
 *(priority=0/3, difficulty=1/3) Wash the indicators of T and Z in order to make it clear
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

	private double xD;
	private double yD;
	private ImagePlus M0map;
	private ImagePlus T1map;
	private ImagePlus T2map;
	private ImageJ ij;
	private HyperMRIT1T2 hyperMRIT1T2;
	private int[] dims;
	private double[] voxs;
	private int nTimepoints;
	private int cEchoes;
	private int cMaps;
	private double[][] dataTimelapseT1;
	private double[][] dataTimelapseT2;
	private double arbitraryLimitBetweenT1andT2=500;
	private final double normValueSpectrum=0.9;
	private ImagePlus mapsImage;
	private ImagePlus echoesImage;
	Color bubbleColor=new Color(220,220,0);
	
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
	
	

	//Parameters of graph and distribution estimation
	boolean rangeDisplayedFullBlocks=true;
	public double[]rangingBoundaries=new double[3];
	public double rangingFactor=0.2;
	public boolean spectrumRangingModeT1=false;
	public boolean spectrumRangingModeT2=false;
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
	boolean isFireLut=true;
	

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
	double [][]t1Times;
	double [][]t2Times;
	double[]timesT1Cute;
	double[]timesT2Cute;

//	double[][]dataTimelapseT1;
//	double[][]dataTimelapseT2;
	double[][][][]dataTimelapseFull;
	int [][]pointsCoords;
	double [][][][]pointsEstimations;
	double[][]dataTimelapseT1Sigmas;
	double[][]dataTimelapseT2Sigmas;
//	double[][]dataTimelapseT2SigmasBicomp;

	double[][]paramsTimelapseT1;
	double[][]paramsTimelapseT2;
//	double[][]paramsTimelapseT1Sigmas;
//	double[][]paramsTimelapseT2Sigmas;


	double[][][]valsSpectrum;

	//Updated just before plot being updated
	double [][]tabFittenT1;
	double [][]tabFittenT2Mono;
	double [][]tabFittenT2Bicomp;
	double [][]tabFittenT2Tricomp;
	double [][]tabFittenT1Cute;
	double [][]tabFittenT2MonoCute;
	double [][]tabFittenT2BicompCute;
	double [][]tabFittenT2TricompCute;
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

	
	
	
	/** Entry points and startup functions*/
	public static void main(String[]args) {
		runExplorer("/home/fernandr/Bureau/Traitements/Sorgho/Test SSM1/Output_normalization/Normalized Hyper MRI.tif");
	}


	
	
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
		IJ.log("Opening hyperimage "+VitimageUtils.imageResume(fullHyp));
		HyperMRIT1T2 hyperMRIT1T2=new HyperMRIT1T2(fullHyp);
		
		MRI_HyperCurvesExplorer mrExplo=new MRI_HyperCurvesExplorer();
		mrExplo.mapsImage=hyperMRIT1T2.getMapsImage(true);
		mrExplo.echoesImage=hyperMRIT1T2.getEchoesImage(true);
		mrExplo.hyperMRIT1T2=hyperMRIT1T2;
		mrExplo.dims=hyperMRIT1T2.dims;
		mrExplo.voxs=hyperMRIT1T2.voxs;
		if(mrExplo.ij==null)mrExplo.ij = new ImageJ();
		else mrExplo.ij = IJ.getInstance();
		IJ.log("Starting MRI Curve Explorer");
		mrExplo.xCor=hyperMRIT1T2.dims[0]/2;
		mrExplo.yCor=hyperMRIT1T2.dims[1]/2;
		mrExplo.zCor=hyperMRIT1T2.dims[2]/2;
		mrExplo.tCor=0;
		mrExplo.nTimepoints=hyperMRIT1T2.nTimepoints;
		
//		mrExplo.extractMaps();
		mrExplo.setupTimesTrTe();
		mrExplo.setupStructures();
		mrExplo.startGui();
	}
	
	public void setupTimesTrTe() {
		this.t1Times=this.hyperMRIT1T2.getT1TrTimes();
		this.t2Times=this.hyperMRIT1T2.getT2TeTimes();
	}
	
	public void setupStructures() {
		initializeScreenConstants();		
		riceEstimator=RiceEstimator.getDefaultRiceEstimatorForNormalizedHyperEchoesT1AndT2Images();

		//Data gathered on a voxel, on the neighbourhood of a voxel, and associated measured sigma values
		this.dataTimelapseT1=new double[this.nTimepoints][];
		this.dataTimelapseT2=new double[this.nTimepoints][];
		this.dataTimelapseFull=new double[this.nTimepoints][][][];
		this.dataTimelapseT1Sigmas=new double[this.nTimepoints][];
		this.dataTimelapseT2Sigmas=new double[this.nTimepoints][];

		//Params estimated from fits
		this.paramsTimelapseT1=new double[this.nTimepoints][2];
		this.paramsTimelapseT2=new double[this.nTimepoints][12];
	//	this.paramsTimelapseT1Sigmas=new double[this.nTimepoints][2];
	//	this.paramsTimelapseT2Sigmas=new double[this.nTimepoints][12];

		//Corresponding curves and nice curves (for display)
		this.tabFittenT1=new double[this.nTimepoints][];
		this.tabFittenT2Mono=new double[this.nTimepoints][];
		this.tabFittenT2Bicomp=new double[this.nTimepoints][];
		this.tabFittenT2Tricomp=new double[this.nTimepoints][];
		this.tabFittenT1Cute=new double[this.nTimepoints][];
		this.tabFittenT2MonoCute=new double[this.nTimepoints][];
		this.tabFittenT2BicompCute=new double[this.nTimepoints][];
		this.tabFittenT2TricompCute=new double[this.nTimepoints][];
		this.timesT1Cute=MRUtils.getProportionalTimes(0,maxT1,xStepCuteT1);
		this.timesT2Cute=MRUtils.getProportionalTimes(0,maxT2,xStepCuteT2);

		//Estimation errors
		this.jitterT1=new double[this.nTimepoints];
		this.jitterT2Mono=new double[this.nTimepoints];
		this.jitterT2Bicomp=new double[this.nTimepoints];
		this.jitterT2Tricomp=new double[this.nTimepoints];

		//Khi and pvalue for estimation
		this.khi2T1=new double[this.nTimepoints];
		this.khi2T2Mono=new double[this.nTimepoints];
		this.khi2T2Bicomp=new double[this.nTimepoints];
		this.khi2T2Tricomp=new double[this.nTimepoints];
		this.pValT1=new double[this.nTimepoints];
		this.pValT2Mono=new double[this.nTimepoints];
		this.pValT2Bicomp=new double[this.nTimepoints];
		this.pValT2Tricomp=new double[this.nTimepoints];
		

		//Computed distribution of T1 and T2 values around the neighbourhood of a voxel
		this.correspondanceCanvas=new int[1000][2];
		this.valsSpectrum=new double[this.nTimepoints][][];
	}
	
	public void startGui() {
	    WindowManager.addWindow(this);
	    instance = this;
		this.imgCan1=new ImageCanvas(this.mapsImage);
		this.imgCan2=new ImageCanvas(this.echoesImage);
        startPlotsAndRoi();
 		initializeGUI();
		imgCan1.getImage().setPosition(1,zCor+1,tCor+1);
		imgCan2.getImage().setPosition(2,zCor+1,tCor+1);
		zoomLevel=imgCan1.getMagnification();
	}


		

	
	
	
	
	
	/** Gui building helpers*/	
	public void repaintAll() {
		imgCan1.repaint();
		imgCan2.repaint();
		plotCan1.repaint();
		plotCan2.repaint();
		plotCan21.repaint();
		plotCan22.repaint();
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
	
	public void initializeGUI() {
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
		info1=new TextField("   Click on the image to compute T1",targetNumberCharsInTexts);info1.setEditable(false);info1.setFont(new Font("Helvetica", 0, 15));info1.setBackground(new Color(230,230,230));
		addComponentToPanel(panel,info1,gbc,
				DX_IMG+DELTA_X               ,   DY_PLOT_1+DELTA_Y,
				DX_PLOT_1+DELTA_X+DX_PLOT_2  ,   DY_TEXT,0,0);

		info2=new TextField("   Press 'h' stroke to get help",targetNumberCharsInTexts);info2.setEditable(false);info2.setFont(new Font("Helvetica", 0, 15));info2.setBackground(new Color(200,200,200));
		addComponentToPanel(panel,info2,gbc,
				DX_IMG+DELTA_X               ,   DY_PLOT_1+DY_TEXT+2*DELTA_Y,
				DX_PLOT_1+DELTA_X+DX_PLOT_2  ,   DY_TEXT,0,0);

		info3=new TextField("   Click on the image to compute T2",targetNumberCharsInTexts);info3.setEditable(false);info3.setFont(new Font("Helvetica", 0, 15));info3.setBackground(new Color(230,230,230));
		addComponentToPanel(panel,info3,gbc,
				DX_IMG+DELTA_X               ,   DY_PLOT_1+2*DY_TEXT+3*DELTA_Y,
				DX_PLOT_1+DELTA_X+DX_PLOT_2  ,   DY_TEXT,0,0);



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
		plotT22 = new Plot("T2 timelapse tracker","T1 and T2 distribution over area (M0-weighted)","Observation time",Plot.DEFAULT_FLAGS-Plot.X_GRID-Plot.Y_GRID+Plot.X_LOG_TICKS+Plot.X_LOG_NUMBERS);
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
		plotT21.setLimits(t0T1,t1T1, 0, nTimepoints);
		plotT22.setLimits(t0T2,t1T2, 0, nTimepoints);
	}

	
	
	

	/** Gui updating functions and callbacks*/	
	@Override
	public void keyTyped(KeyEvent e) {
		IJ.log("KEY PRESSED : "+e.getKeyChar());

		///////ACTIONS TO CHANGE SIZE OF AREA
		if (e.getKeyChar()=='p' && statusRoi!=2) {
			this.autoSizingTimePlots=!this.autoSizingTimePlots;
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='z' && statusRoi!=2) {
			if(this.crossThick<dims[2] && statusRoi!=2)this.crossThick++;
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='s' && statusRoi!=2) {
			if(this.crossThick>0 && statusRoi!=2)this.crossThick--;
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='a' && statusRoi!=2) {
			if(this.crossWidth<dims[0])this.crossWidth++;
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='q' && statusRoi!=2) {
			if(this.crossWidth>0)this.crossWidth--;
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='e') {
			this.nRepetMonteCarlo+=20;
			IJ.log("Monte carlo go to N="+this.nRepetMonteCarlo+" repetitions");
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='d') {
			this.nRepetMonteCarlo-=20;
			if(this.nRepetMonteCarlo<0)this.nRepetMonteCarlo=0;
			IJ.log("Monte carlo go to N="+this.nRepetMonteCarlo+" repetitions");
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}

	
		//Switch multi-thread
		if (e.getKeyChar()=='k') {
			this.multiThreaded=!this.multiThreaded;
			System.out.println("Switching computation mode. New mode = "+(this.multiThreaded ? ("multithreaded over "+VitimageUtils.getNbCores()+" cores.") : "single core."));
			this.computeResultsAgain();
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
			}
			if(currentCanvas==WIN_PLOT1 || currentCanvas==WIN_PLOT2) {
				plotCan1.zoomIn(xMouse,yMouse);
				plotCan2.zoomIn(xMouse,yMouse);
			}
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='-') {
			if(statusRoi==2) {}//VitiDialogs.notYet("Warning : no zoom in / out during Roi mode");return;
			IJ.log("Zoom out");
			if( (currentCanvas==WIN_T1 || currentCanvas==WIN_T2) ) {
				imgCan1.setMagnification(imgCan1.getMagnification()/2.0);
				imgCan1.zoomIn(xMouse,yMouse);
				imgCan2.setMagnification(imgCan2.getMagnification()/2.0);
				imgCan2.zoomIn(xMouse,yMouse);
				zoomLevel=imgCan1.getMagnification();
			}
			if(currentCanvas==WIN_PLOT1 || currentCanvas==WIN_PLOT2) {
				zoomLevel=1;
				plotCan1.setMagnification(zoomLevel);
				plotCan2.setMagnification(zoomLevel);
			}
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}


		
		
		///////LUT	
		if (e.getKeyChar()=='f') {
			isFireLut=!isFireLut;
		}
		///////DISPLACEMENT IN HYPERIMAGE		
		if (e.getKeyChar()=='8' || e.getKeyChar()=='2' || e.getKeyChar()=='4' || e.getKeyChar()=='6' || e.getKeyChar()=='1' || e.getKeyChar()=='3' || e.getKeyChar()=='7' || e.getKeyChar()=='9'  || e.getKeyChar()=='f') {
			//Move t
			if((e.getKeyChar()=='8') && (tCor<nTimepoints-1) )tCor+=1;
			if((e.getKeyChar()=='2') && (tCor>0) )tCor-=1;

			//Move z
			if((e.getKeyChar()=='6') && (zCor<dims[2]-1) )zCor+=1;
			if((e.getKeyChar()=='4') && (zCor>0) )zCor-=1;

			//Move c of maps
			if((e.getKeyChar()=='9') && (cMaps<2) )cMaps+=1;
			if((e.getKeyChar()=='7') && (cMaps>0) )cMaps-=1;

			//Move c of echos
			if((e.getKeyChar()=='3') && (cEchoes<hyperMRIT1T2.nChannels-1) )cEchoes+=1;
			if((e.getKeyChar()=='1') && (cEchoes>0) )cEchoes-=1;

			if(e.getKeyChar()=='4' || e.getKeyChar()=='6') {
				actualizeMriObservationsBasedOnData();
				computeResultsAgain();
			}
			else {
				identifyRangedData();
			}
			displayResultsAgain();
			
			IJ.log("Changing coordinates to t="+tCor+", z="+zCor+" cMaps="+cMaps+" cEchoes="+cEchoes);
			imgCan1.getImage().setPosition(cMaps+1,zCor+1,tCor+1);
			imgCan1.setImageUpdated();
			if(isFireLut) IJ.run(imgCan1.getImage(),"Fire","");
			else IJ.run(imgCan1.getImage(),"Grays","");
			imgCan2.getImage().setPosition(cEchoes+1,zCor+1,tCor+1);
			imgCan2.setImageUpdated();
			if(isFireLut) IJ.run(imgCan2.getImage(),"Fire","");
			else IJ.run(imgCan2.getImage(),"Grays","");
			
			
			if(statusRoi!=2 || (userRoi ==null));
			else {
				imgCan1.setOverlay(new Overlay(userRoi));
				imgCan2.setOverlay(new Overlay(userRoi));
			}
			
			
		}


		

		
		
		
		///////USER CUSTOM ROI MANAGEMENT
		if (e.getKeyChar()=='r') {
			if(statusRoi==0) {				
				zoomLevel=imgCan1.getMagnification();
				if(zoomLevel>1) {
					while(zoomLevel>1) {
						imgCan1.setMagnification(imgCan1.getMagnification()/2.0);
						imgCan1.zoomIn(xMouse,yMouse);
						imgCan2.setMagnification(imgCan2.getMagnification()/2.0);
						imgCan2.zoomIn(xMouse,yMouse);
						repaintAll();
						zoomLevel=imgCan1.getMagnification();						
					}
				}
				else if(zoomLevel<1) {
					while(zoomLevel<1) {
						imgCan1.zoomIn(xMouse,yMouse);
						imgCan2.zoomIn(xMouse,yMouse);
						repaintAll();
						zoomLevel=imgCan1.getMagnification();						
					}
				}					
				IJ.log("Please select a ROI now, add it to the Roi manager, then hit the 'r' strike again");
				String file=VitiDialogs.chooseOneRoiPathUI("Choose a .roi file", "");
				if(file==null)return;
				this.userRoi=new Opener().openRoi(file);
				IJ.selectWindow("MRI Curve explorer V2");
				statusRoi=2;
				this.userRoi.setPosition(0);
				repaintAll();
				this.imgCan1.setOverlay(new Overlay(userRoi));
				this.imgCan2.setOverlay(new Overlay(userRoi));
				computeResultsAgain();
				displayResultsAgain();
				IJ.log("Please hit the 'r' strike again to quit the Roi mode");
			}
			else if(statusRoi==2) {
				statusRoi=0;
				this.userRoi=null;
				this.nPtsCur=1;
				computeResultsAgain();
				displayResultsAgain();
			}
			else {
				computeResultsAgain();
				displayResultsAgain();
			}
		}
		
		
		
		///MISCELLANEOUS
		if (e.getKeyChar()==',') {
			//just update
			this.rangingFactor=this.rangingFactor*0.7;
			this.actualizeRangingBoundaries();
			this.identifyRangedData();
			displayResultsAgain();
		}
		if (e.getKeyChar()==';') {
			//just update
			this.rangingFactor=this.rangingFactor*1.4;
			this.actualizeRangingBoundaries();
			this.identifyRangedData();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='u') {
			//just update
			this.actualizeRangingBoundaries();
			this.identifyRangedData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='b') {
			//just update
			this.rangeDisplayedFullBlocks=!this.rangeDisplayedFullBlocks;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			this.identifyRangedData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='c') {
			this.multiThreaded=!this.multiThreaded;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			this.identifyRangedData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='v') {
			this.separateNormalizationSpectrum=!this.separateNormalizationSpectrum;
			actualizeMriObservationsBasedOnData();
			this.actualizeRangingBoundaries();
			this.identifyRangedData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='l') {
			if(true)return;
			if(this.fitAlgorithm==MRUtils.LM) {
				this.fitAlgorithm=MRUtils.SIMPLEX;
				IJ.log("Switching to SIMPLEX");
			}
			else if(this.fitAlgorithm==MRUtils.SIMPLEX) {
				this.fitAlgorithm=MRUtils.LM;
				IJ.log("Switching to Levenberg-Marquardt");
			}
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='t') {
			thickCurve=1-thickCurve;
			displayResultsAgain();
		}
		if (e.getKeyChar()=='g') {
			gaussianSpectrum=(gaussianSpectrum+1)%4;
			computeResultsAgain();
			displayResultsAgain();
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
		repaintAll();
	}
	

	
	public void actualizeCursor() {
		int xMouseCenter=(int) Math.round(xMouse+  (   Math.floor(xD)+0.5-xD  )*zoomLevel);//screen location of square center, according to xMouse	
		int yMouseCenter=(int) Math.round(yMouse+  (   Math.floor(yD)+0.5-yD  )*zoomLevel);
		double sizeOnScreen=(0.5+crossWidth)*zoomLevel;//half radius

		PointRoi prT1=new PointRoi(xMouseCenter,yMouseCenter,sizeOfCursor()+" yellow hybrid");
		Overlay overT1=new Overlay(prT1);
		Overlay overT2=new Overlay(prT1);

		//If no Roi is defined, draw a square with a cross in the center of it
		if(statusRoi!=2) {
			overT1.add(new Roi(xMouseCenter-sizeOnScreen,yMouseCenter-sizeOnScreen,2*sizeOnScreen,2*sizeOnScreen));
			overT2.add(new Roi(xMouseCenter-sizeOnScreen,yMouseCenter-sizeOnScreen,2*sizeOnScreen,2*sizeOnScreen));
		}
		else if(zoomLevel==1) {
			overT1.add(userRoi);
			overT2.add(userRoi);
		}
		
		//If spectrum ranging is active, draw squares around the pixels which estimate Ttimes falls into the selected interval
		if(this.spectrumRangingModeT1 || this.spectrumRangingModeT2 ) {
			for(int pt=0;pt<this.rangeRoiPoints.size();pt++) {
				int dx=rangeRoiPoints.get(pt)[0]-xCor;//Relative coordinates to the cross
				int dy=rangeRoiPoints.get(pt)[1]-yCor;
				Roi r=new Roi(xMouseCenter+(dx-0.5)*zoomLevel,yMouseCenter+(dy-0.5)*zoomLevel,zoomLevel,zoomLevel);
				if(rangeDisplayedFullBlocks)r.setFillColor(new Color(0,0,255));
				else r.setStrokeColor(new Color(0,0,255));
				overT1.add(r);
				overT2.add(r);
			}
		}
		imgCan1.setOverlay(overT1);
		imgCan2.setOverlay(overT2);
	}
	
	public String sizeOfCursor() {
		return "medium";
	}
		
	public void actualizeRangingBoundaries() {
		double borderLeft=76;
		double borderRight=20;
		double lengPlot=DX_PLOT_2-borderRight-borderLeft;
		double xPos=this.xMouseRange-borderLeft;
		double[]vals=spectrumRangingModeT1 ? plotT21.getLimits() : plotT22.getLimits();
		double factMul=vals[1]/vals[0];
		double rangCenter=Math.pow(factMul,xPos*1.0/lengPlot)*vals[0];
		this.rangingBoundaries=new double[] {rangCenter*(1-this.rangingFactor),rangCenter,rangCenter*(1+this.rangingFactor)};
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
	public void computeResultsAgain() {
		actualizeMeanEstimations();
		actualizeDisplayedNumbers();
		actualizeExportSentence();
		Timer tim2=new Timer();
		System.out.println("\n\nStarting estimation");
		if(multiThreaded) computeEstimationsForAllPointsMultiThreaded();
		else computeEstimationsForAllPoints();
		tim2.print("Fit T1, T2mono and T2 bi-exp on "+this.nPtsCur+" points over "+this.nTimepoints+" timepoints. Total computation time : ");
		actualizeSpectrumCurves();
		identifyRangedData();
	}

	public void actualizeMriObservationsBasedOnData() {
		//Update for all time, for all voxels selected
		if(statusRoi!=2 || userRoi==null) {
			this.dataTimelapseFull=hyperMRIT1T2.getFullMRISignalAroundThisVoxel((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick);
			this.pointsCoords=hyperMRIT1T2.getCoordinatesAroundThisVoxel((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick);
		}
		else {
			IJ.log("Actualizing data from Roi, at coordinates Z="+zCor+" , T="+tCor);
			this.pointsCoords=VitimageUtils.getRoiAsCoords(this.userRoi);
			this.dataTimelapseFull=hyperMRIT1T2.getFullMRISignalInTheseCoordinates(this.pointsCoords);
		}
		this.nPtsCur=this.dataTimelapseFull[0].length;

		//For each time, compute the mean over voxels
		double[]vals;
		for(int t=0;t<this.nTimepoints;t++) {			
			this.dataTimelapseT1[t]=new double[this.dataTimelapseFull[t][0][0].length];
			this.dataTimelapseT2[t]=new double[this.dataTimelapseFull[t][0][1].length];
			this.dataTimelapseT1Sigmas[t]=new double[this.dataTimelapseFull[t][0][1].length];
			this.dataTimelapseT2Sigmas[t]=new double[this.dataTimelapseFull[t][0][1].length];
			for(int ec=0;ec<this.dataTimelapseFull[t][0][0].length;ec++) {
				vals=new double[this.nPtsCur];
				for(int vo=0;vo<this.nPtsCur;vo++) {
					vals[vo]=this.dataTimelapseFull[t][vo][0][ec];
				}
				double []stats=VitimageUtils.statistics1D(vals);
				this.dataTimelapseT1[t][ec]=stats[0];
				this.dataTimelapseT1Sigmas[t][ec]=stats[1];						
			}
			for(int ec=0;ec<this.dataTimelapseFull[t][0][1].length;ec++) {
				vals=new double[this.nPtsCur];
				for(int vo=0;vo<this.nPtsCur;vo++) {
					vals[vo]=this.dataTimelapseFull[t][vo][1][ec];
				}
				double []stats=VitimageUtils.statistics1D(vals);
				this.dataTimelapseT2[t][ec]=stats[0];
				this.dataTimelapseT2Sigmas[t][ec]=stats[1];						
			}
		}
	}
			
	public void actualizeMeanEstimations() {
		for(int tim=0;tim< this.nTimepoints;tim++) {
			if(tim==tCor) {
				IJ.log("  --> Estimations for current time "+tim+" / "+this.nTimepoints+" at z="+zCor);
		 		IJ.log("    T1 MRI times = "+TransformUtils.stringVectorNDou(t1Times[tim], "")+" with sigma rice = "+VitimageUtils.dou(hyperMRIT1T2.tabSigmasT1Seq[tim][zCor])+"\n    T1 MRI data = "+TransformUtils.stringVectorNDou(dataTimelapseT1[tim], ""));
			}
	 		Object[] obj=fitAndEvaluate(t1Times[tim],this.timesT1Cute,dataTimelapseT1[tim],hyperMRIT1T2.tabSigmasT1Seq[tim][zCor],this.fitAlgorithm,MRUtils.T1_RECOVERY_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
	 		paramsTimelapseT1[tim]=(double[]) obj[0];
	 		//paramsTimelapseT1Sigmas[tim]=(double[]) obj[1];
	 		tabFittenT1[tim]=(double[]) obj[2];
	 		tabFittenT1Cute[tim]=(double[]) obj[3];
	 		khi2T1[tim]=(double) obj[4];
	 		pValT1[tim]=(double) obj[5];
	 		jitterT1[tim]=(double) obj[6];
	 		if(this.nPtsCur==1)this.dataTimelapseT1Sigmas[tim]=riceEstimator.estimateSigmas(hyperMRIT1T2.tabSigmasT1Seq[tim][zCor], this.tabFittenT1[tim]);
	
	 		//Estimer T2 Monocomp
	 		if(tim==tCor)IJ.log("    T2 MRI times = "+TransformUtils.stringVectorNDou(t2Times[tim], "")+" with sigma rice = "+VitimageUtils.dou(hyperMRIT1T2.tabSigmasT2Seq[tim][zCor])+"\n    T2 MRI data = "+TransformUtils.stringVectorNDou(dataTimelapseT2[tim], ""));
	 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseT2[tim],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor],this.fitAlgorithm,MRUtils.T2_RELAX_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
	 		for(int prm=0;prm<2;prm++) {
	 			paramsTimelapseT2[tim][prm]=((double[]) obj[0])[prm];
	 		//	paramsTimelapseT2Sigmas[tim][prm]=((double[]) obj[1])[prm];
	 		}
	 		tabFittenT2Mono[tim]=(double[]) obj[2];
	 		tabFittenT2MonoCute[tim]=(double[]) obj[3];
	 		khi2T2Mono[tim]=(double) obj[4];
	 		pValT2Mono[tim]=(double) obj[5];
	 		jitterT2Mono[tim]=(double) obj[6];
	 		if(this.nPtsCur==1)this.dataTimelapseT2Sigmas[tim]=riceEstimator.estimateSigmas(hyperMRIT1T2.tabSigmasT2Seq[tim][zCor], this.tabFittenT2Mono[tim]);
		
	 		//Estimer T2 Bicomp
	 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseT2[tim],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor],this.fitAlgorithm,MRUtils.MULTICOMP_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
	 		for(int prm=2;prm<6;prm++) {
	 			paramsTimelapseT2[tim][prm]=((double[]) obj[0])[prm-2];
	 		//	paramsTimelapseT2Sigmas[tim][prm]=((double[]) obj[1])[prm-2];
	 		}
	 		tabFittenT2Bicomp[tim]=(double[]) obj[2];
	 		tabFittenT2BicompCute[tim]=(double[]) obj[3];
	 		khi2T2Bicomp[tim]=(double) obj[4];
	 		pValT2Bicomp[tim]=(double) obj[5];
	 		jitterT2Bicomp[tim]=(double) obj[6];
	 		//if(this.nPtsCur==1)this.dataTimelapseT2Sigmas[tim]=riceEstimator.estimateSigmas(hyperMRIT1T2.tabSigmasT2Seq[tim][zCor], this.tabFittenT2Bicomp[tim]);
		}		
	}

	public void computeEstimationsForAllPoints() {
		pointsEstimations=new double[this.nTimepoints][this.nPtsCur][2][8];
		for(int tim=0;tim< this.nTimepoints;tim++) {
			for(int vo=0;vo<this.nPtsCur;vo++) {
				Timer tt=new Timer();
				Object[] obj=fitAndEvaluate(t1Times[tim],this.timesT1Cute,dataTimelapseFull[tim][vo][0],hyperMRIT1T2.tabSigmasT1Seq[tim][zCor],this.fitAlgorithm,MRUtils.T1_RECOVERY_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
				pointsEstimations[tim][vo][0][0]=((double[])obj[0])[0];
				pointsEstimations[tim][vo][0][1]=((double[])obj[0])[1];
				pointsEstimations[tim][vo][0][2]=((double) obj[6])>=9999 ? 0 : 1;
		
		 		//Estimer T2 Monocomp
				tt=new Timer();
		 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseFull[tim][vo][1],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor],this.fitAlgorithm,MRUtils.T2_RELAX_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
				pointsEstimations[tim][vo][1][0]=((double[])obj[0])[0];
				pointsEstimations[tim][vo][1][1]=((double[])obj[0])[1];
				pointsEstimations[tim][vo][1][2]=((double) obj[6])>=9999 ? 0 : 1;
				double khiMono=(double) obj[4];

		 		//Estimer T2 Bicomp
				tt=new Timer();
		 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseFull[tim][vo][1],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor],this.fitAlgorithm,MRUtils.MULTICOMP_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
				pointsEstimations[tim][vo][1][3]=((double[])obj[0])[0];
				pointsEstimations[tim][vo][1][4]=((double[])obj[0])[1];
				pointsEstimations[tim][vo][1][5]=((double[])obj[0])[2];
				pointsEstimations[tim][vo][1][6]=((double[])obj[0])[3];
				pointsEstimations[tim][vo][1][7]=((double) obj[6])>=9999 ? 0 : 1;
				double khiBi=(double) obj[4];

				//Select best T2
				if(khiMono<=khiBi && pointsEstimations[tim][vo][1][2]==1)pointsEstimations[tim][vo][1][7]=0;
				if(khiMono>khiBi && pointsEstimations[tim][vo][1][7]==1)pointsEstimations[tim][vo][1][2]=0;		 		
			}
		}		
	}
	
	public void computeEstimationsForAllPointsMultiThreaded() {
		//Prepare input data and output places for threads
		int nThreads=VitimageUtils.getNbCores();
		Thread[]threads=VitimageUtils.newThreadArray(nThreads);
		int[][]listVoxOnThreads=VitimageUtils.listForThreads(this.nPtsCur, nThreads);
		AtomicInteger atomNumThread=new AtomicInteger(0);
		AtomicInteger curProcessedBlock=new AtomicInteger(0);
		final double[][][][][]dataParams=new double[nThreads][][][][];
		final double[][][][][]estimatedParams=new double[nThreads][][][][];
		final int finalFitAlgorithm=this.fitAlgorithm;
		final int finalRepetMonteCarlo=this.nRepetMonteCarlo;
		final int nData=this.nPtsCur;
		final int nTimes=this.nTimepoints;
		for(int nt=0;nt<nThreads;nt++) {
			dataParams[nt]=new double[listVoxOnThreads[nt].length][][][];
			estimatedParams[nt]=new double[listVoxOnThreads[nt].length][][][];
			for(int nVo=0;nVo<listVoxOnThreads[nt].length;nVo++) {
				dataParams[nt][nVo]=new double[this.nTimepoints][][];
				estimatedParams[nt][nVo]=new double[this.nTimepoints][][];
				for(int nTime=0;nTime<this.nTimepoints;nTime++) {
					dataParams[nt][nVo][nTime]=new double[6][];
					dataParams[nt][nVo][nTime][0]=t1Times[nTime];
					dataParams[nt][nVo][nTime][1]=dataTimelapseFull[nTime][listVoxOnThreads[nt][nVo]][0];
					dataParams[nt][nVo][nTime][2]=new double[] {hyperMRIT1T2.tabSigmasT1Seq[nTime][zCor]};
					dataParams[nt][nVo][nTime][3]=t2Times[nTime];
					dataParams[nt][nVo][nTime][4]=dataTimelapseFull[nTime][listVoxOnThreads[nt][nVo]][1];
					dataParams[nt][nVo][nTime][5]=new double[] {hyperMRIT1T2.tabSigmasT2Seq[nTime][zCor]};
				}
			}
		}

	
		
		//Run threads
		for (int ithread = 0; ithread < nThreads; ithread++) {  
			threads[ithread] = new Thread() {  { setPriority(Thread.NORM_PRIORITY); }  
				public void run() {  
					try {
						int numThread=atomNumThread.getAndIncrement();
						for (int numVo=0;numVo<dataParams[numThread].length;numVo++) {
							int numBl=curProcessedBlock.getAndIncrement();
							if(nData>10 && (numBl%(nData/10)==0))IJ.log(" "+VitimageUtils.dou((numBl*100.0)/nData)+"%");
							for (int numTime=0;numTime<nTimes;numTime++) {
								double[][]tempParams=new double[2][8];
								Object[] obj=fitAndEvaluate(dataParams[numThread][numVo][numTime][0],null,dataParams[numThread][numVo][numTime][1],dataParams[numThread][numVo][numTime][2][0],finalFitAlgorithm,MRUtils.T1_RECOVERY_RICE,finalRepetMonteCarlo,1,false);
								tempParams[0][0]=((double[])obj[0])[0];
								tempParams[0][1]=((double[])obj[0])[1];
								tempParams[0][2]=((double) obj[6])>=9999 ? 0 : 1;
						
						 		//Estimer T2 Monocomp
						 		obj=fitAndEvaluate(dataParams[numThread][numVo][numTime][3],null,dataParams[numThread][numVo][numTime][4],dataParams[numThread][numVo][numTime][5][0],finalFitAlgorithm,MRUtils.T2_RELAX_RICE,finalRepetMonteCarlo,1,false);	
						 		tempParams[1][0]=((double[])obj[0])[0];
						 		tempParams[1][1]=((double[])obj[0])[1];
						 		tempParams[1][2]=((double) obj[6])>=9999 ? 0 : 1;
								double khiMono=(double) obj[4];

						 		//Estimer T2 Bicomp
						 		obj=fitAndEvaluate(dataParams[numThread][numVo][numTime][3],null,dataParams[numThread][numVo][numTime][4],dataParams[numThread][numVo][numTime][5][0],finalFitAlgorithm,MRUtils.MULTICOMP_RICE,finalRepetMonteCarlo,1,false);
						 		tempParams[1][3]=((double[])obj[0])[0];
						 		tempParams[1][4]=((double[])obj[0])[1];
						 		tempParams[1][5]=((double[])obj[0])[2];
						 		tempParams[1][6]=((double[])obj[0])[3];
						 		tempParams[1][7]=((double) obj[6])>=9999 ? 0 : 1;
								double khiBi=(double) obj[4];

								//Select best T2
								if(khiMono<=khiBi && tempParams[1][2]==1)tempParams[1][7]=0;
								if(khiMono>khiBi && tempParams[1][7]==1)tempParams[1][2]=0;		 		
								estimatedParams[numThread][numVo][numTime]=tempParams;
							}							
						}
					} catch(Exception ie) {}
				} 
			};  		
		}				
		VitimageUtils.startAndJoin(threads);

		
		
		
		//Gather results
		pointsEstimations=new double[this.nTimepoints][this.nPtsCur][2][8];
		for(int nt=0;nt<nThreads;nt++) {
			for(int nVo=0;nVo<listVoxOnThreads[nt].length;nVo++) {
				for(int nTime=0;nTime<this.nTimepoints;nTime++) {
					pointsEstimations[nTime][listVoxOnThreads[nt][nVo]]=estimatedParams[nt][nVo][nTime];
				}
			}
		}
	}

	public Object[] fitAndEvaluate(double[]tabTimes,double[]tabTimesCute,double[]tabData,double sigmaRice,int fitAlgorithm,int fitCurveType,int nbMonteCarloSimulations,int nbPts,boolean niceCurveComputation) {
		int nParams=(fitCurveType==MRUtils.T1_RECOVERY_RICE) ? 2 : (fitCurveType==MRUtils.T2_RELAX_RICE) ? 2 :(fitCurveType==MRUtils.MULTICOMP_RICE) ? 4 : 5; 
		boolean isT1=(fitCurveType==MRUtils.T1_RECOVERY_RICE); 
		double []estimatedParams=new double[nParams];
		double []estimatedSigmas=new double[nParams];
		double [][]estimation=MRUtils.makeFit(tabTimes, tabData,fitCurveType,fitAlgorithm,100,sigmaRice,nbMonteCarloSimulations,nbPts,this.riceEstimator);
		for(int i=0;i<nParams;i++) {estimatedParams[i]=estimation[0][i] ;  estimatedSigmas[i]=estimation[1][i];}
		double []tabFitten=MRUtils.fittenRelaxationCurve(tabTimes,estimatedParams,sigmaRice,fitCurveType);
		double []tabFittenCute=niceCurveComputation ? MRUtils.fittenRelaxationCurve(tabTimesCute,estimatedParams,sigmaRice,fitCurveType) : null;
		double[]accs=MRUtils.fittingAccuracies(tabData,tabTimes,sigmaRice,estimatedParams,fitCurveType,false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
		if((nParams==4) && (estimatedParams[1]>estimatedParams[3])) {
			double tmp=estimatedParams[1];estimatedParams[1]=estimatedParams[3];estimatedParams[3]=tmp;
			tmp=estimatedParams[0];estimatedParams[0]=estimatedParams[2];estimatedParams[2]=tmp;
			tmp=estimatedSigmas[1];estimatedSigmas[1]=estimatedSigmas[3];estimatedSigmas[3]=tmp;
			tmp=estimatedSigmas[0];estimatedSigmas[0]=estimatedSigmas[2];estimatedSigmas[2]=tmp;
		}		
		
		
		double jitter=0;
		double[]statsRice=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(sigmaRice);
		if((!isT1) && tabData[0]<(statsRice[0]+3*statsRice[1]) ) {jitter=9999;accs[1]=100;System.out.println("Jitter set cause 1 - "+tabData[0]+" , "+statsRice[0]+" , "+statsRice[1]);}
		if(isT1 && (estimatedParams[0]<0 || estimatedParams[1]<minAcceptableT1 || estimatedParams[1]>maxT1 || estimatedParams[0]>maxAcceptableM0T1  ) )  {jitter=9999;accs[1]=100;System.out.println("Jitter set cause 2 - "+TransformUtils.stringVectorN(estimatedParams, ""));}
		if(!isT1 && (nParams==2) && (estimatedParams[0]<0 || estimatedParams[1]<0 || estimatedParams[0]>maxAcceptableM0T2 || estimatedParams[1]>maxAcceptableT2) ) {jitter=9999;accs[1]=100;System.out.println("Jitter set cause 3 - "+TransformUtils.stringVectorN(estimatedParams, ""));}
		if(!isT1 && (nParams==4) && (estimatedParams[0]<0 || estimatedParams[1]<minAcceptableT2 || 
				estimatedParams[2]<0 || estimatedParams[3]>maxAcceptableT2 || 
				estimatedParams[0]>maxAcceptableM0T2 ||   estimatedParams[2]>maxAcceptableM0T2) ) {jitter=9999;accs[1]=100;System.out.println("Jitter set cause 4 - "+TransformUtils.stringVectorN(estimatedParams, ""));}		
		return new Object[] {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accs[0],accs[1],jitter};
	}
	
	public void identifyRangedData() {
		this.rangeRoiPoints=new ArrayList<int[]>();
		if(this.spectrumRangingModeT1) {
			//Ranging T1 values
			for(int p=0;p<pointsCoords.length;p++) {
				if ( (pointsCoords[p][2]==zCor) && /*Point is on current slice*/
						 this.pointsEstimations[tCor][p][0][2]==1 /*Point estimation is valid (jitter not too big)*/ && 
						 (this.pointsEstimations[tCor][p][0][1]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][1]>this.rangingBoundaries[0]) /*Estimated T1 is in the range*/
					 ) {
					rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
			}
			System.out.println("Ranging T1 : ["+this.rangingBoundaries[0]+"-"+this.rangingBoundaries[2]+"] . After update : "+rangeRoiPoints.size()+" points");
		}
		if(this.spectrumRangingModeT2) {
			//Ranging T2 values
			for(int p=0;p<pointsCoords.length;p++) {
				if ( (pointsCoords[p][2]==zCor) && /*Point is on current slice*/
					 (this.pointsEstimations[tCor][p][1][2]==1) /*Point estimation is valid in T2 monomodal (jitter not too big)*/ && 
					 (this.pointsEstimations[tCor][p][1][1]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][1][1]>this.rangingBoundaries[0]) /*Estimated T2 is in the range*/
				 ) {
					 rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
				if ( (pointsCoords[p][2]==zCor) && /*Point is on current slice*/
					 (this.pointsEstimations[tCor][p][1][7]==1) /*Point estimation is valid in T2 monomodal (jitter not too big)*/ && 
					 (   (this.pointsEstimations[tCor][p][1][4]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][1][4]>this.rangingBoundaries[0]) || /*Estimated short T2 is in the range*/
					     (this.pointsEstimations[tCor][p][1][6]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][1][6]>this.rangingBoundaries[0]) ) /*Estimated long T2 is in the range*/
					 ) {
					 rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
			}
			System.out.println("Ranging T2 : ["+this.rangingBoundaries[0]+"-"+this.rangingBoundaries[2]+"] . After update : "+rangeRoiPoints.size()+" points");
		}
	}
	
	public double[][] computeSpectrumCurve(int time) {
		int numberBins=150 /(1) ;//over 3 decades it makes every 1.2% or every 10 %
		final double[]binInf=new double[numberBins];
		final double[]binSup=new double[numberBins];
		final double[]binMed=new double[numberBins];
		double[]histoM0T1=new double[numberBins];
		double[]histoM0T2=new double[numberBins];
		double [][]output=new double[3][numberBins];
		double minBin=t0T2;
		double maxBin=t1T1;
		double multFromMinToMax=maxBin/minBin;
		for(int i=0;i<numberBins;i++) {
			binInf[i]=minBin*Math.pow(multFromMinToMax, i*1.0/numberBins);
			binSup[i]=minBin*Math.pow(multFromMinToMax, (i+1)*1.0/numberBins);
			binMed[i]=binSup[i]*0.5+binInf[i]*0.5;
		}
		//Bin results
		for(int p=0;p<this.nPtsCur;p++) {
			if(pointsEstimations[time][p][0][2]==1) {//Fit T1 is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][1]<binSup[bin]) && (pointsEstimations[time][p][0][1]>=binInf[bin]) ) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][0];
				}
			}

			if(pointsEstimations[time][p][1][2]>=1) {//Fit T2 mono is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][1][1]<binSup[bin]) && (pointsEstimations[time][p][1][1]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][1][0];
				}
			}

			if(pointsEstimations[time][p][1][7]>=1) {//Fit T2 bicomp is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][1][4]<binSup[bin]) && (pointsEstimations[time][p][1][4]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][1][3];
				}
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][1][6]<binSup[bin]) && (pointsEstimations[time][p][1][6]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][1][5];
				}
			}
		}
				
		//Smooth this histogram with a factor to be defined, maybe depending on the estimation error on parameters
		output[0]=smoothHisto(histoM0T1,gaussianSpectrum);
		output[1]=smoothHisto(histoM0T2,gaussianSpectrum);
		output[2]=binMed;
		return output;
	}
		





	/** Update graphs using new computed results*/		
	public void actualizeFirstPlots() {
		int incrT1=0;
		int incrT2=0;
		int deltaXt2=3;
		int bubbleSize=30;
		Color crossColor=new Color (60,60,255);
		Color t1Moy=new Color (255,50,50);
		String []tabNamesT1=new String[]{"-Fit with noise estimation"};
		if(this.autoSizingTimePlots) {
			maxPlotYT1=VitimageUtils.max(dataTimelapseT1[tCor])*1.3;
			maxPlotYT2=VitimageUtils.max(dataTimelapseT2[tCor])*1.3;
		}

		//BUBBLE SHOWING CURRENT ECHO
		if(cEchoes<this.t1Times[tCor].length) {
			plotT1.setColor(bubbleColor);
	        plotT1.setLineWidth(bubbleSize);
	        plotT1.replace(incrT1++, "line",new double[]{t1Times[tCor][cEchoes],t1Times[tCor][cEchoes]},new double[]{dataTimelapseT1[tCor][cEchoes],dataTimelapseT1[tCor][cEchoes]});//Afficher le mean+sigma
		}
		else {
			plotT2.setColor(bubbleColor);
	        plotT2.setLineWidth(bubbleSize);
	        plotT2.replace(incrT2++, "line",new double[]{t2Times[tCor][cEchoes-this.t1Times[tCor].length],t2Times[tCor][cEchoes-this.t1Times[tCor].length]},new double[]{dataTimelapseT2[tCor][cEchoes-this.t1Times[tCor].length],dataTimelapseT2[tCor][cEchoes-this.t1Times[tCor].length]});
		}
	
        //NOISE
        plotT1.setLineWidth(1);
		double[]statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMRIT1T2.tabSigmasT1Seq[tCor][zCor]);
		double nSum=statusRoi==2 ? Math.sqrt(this.nPtsCur) : (1+crossWidth);

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
		plotT1.replace(incrT1++,"x",t1Times[tCor],dataTimelapseT1[tCor]);//Afficher les points IRM
        plotT1.setLineWidth(2);
        

        //FIT T1
        if(jitterT1[tCor]>=9999)        		plotT1.setColor(Color.gray);
        else plotT1.setColor(t1Moy);
        
        plotT1.setLineWidth(1+thickCurve);
        plotT1.replace(incrT1++, "line",timesT1Cute,tabFittenT1Cute[tCor]);//Afficher la courbe 
        
        //PARAM T1 MONO
        if(jitterT1[tCor]>=9999)        		plotT2.setColor(Color.gray);
        else         plotT1.setColor(t1Moy);
        plotT1.setLineWidth(3);
		plotT1.replace(incrT1++, "line",new double[]{0,(double)(paramsTimelapseT1[tCor][1])},new double[]{maxPlotYT1*0.9,maxPlotYT1*0.9});//Afficher le T1
        plotT1.setColor(crossColor);
        plotT1.setLineWidth(1);
		for(int t=0;t<this.t1Times[tCor].length;t++)plotT1.replace(incrT1++, "line",new double[]{this.t1Times[tCor][t]-deltaXt2,this.t1Times[tCor][t]-deltaXt2},new double[]{dataTimelapseT1[tCor][t]-dataTimelapseT1Sigmas[tCor][t],dataTimelapseT1[tCor][t]+dataTimelapseT1Sigmas[tCor][t]});//Afficher le T2

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
		Color t2Moy=new Color (0,125,0);
		Color t2Bi=new Color (10,185,10);
		String[] tabNamesT2=new String[]{"One T2","Two T2"};

		//NOISE LEVEL
		plotT2.setLineWidth(1);
		statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMRIT1T2.tabSigmasT2Seq[tCor][zCor]);
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
		plotT2.replace(incrT2++,"x",t2Times[tCor],dataTimelapseT2[tCor]);//Afficher les points IRM	//	plotT2.replace(incrT2++, "line",new double[]{0,0},new double[]{0,0});
        plotT2.setLineWidth(2);
        
        //COURBES FITTED
        if(jitterT2Mono[tCor]>=9999)        		plotT2.setColor(Color.gray);
        else plotT2.setColor(t2Moy);
        plotT2.setLineWidth(1+thickCurve);
        plotT2.replace(incrT2++, "line",timesT2Cute,tabFittenT2MonoCute[tCor]);//Afficher la courbe monocomp
        if(jitterT2Bicomp[tCor]>=9999)        		plotT2.setColor(Color.gray);
        else plotT2.setColor(t2Bi);
        plotT2.setLineWidth(1+thickCurve);
        plotT2.replace(incrT2++, "line",timesT2Cute,tabFittenT2BicompCute[tCor]);//Afficher la courbe bicomp
        
        //PARAM T2 MONO ET STD MONO
        if(jitterT2Mono[tCor]>=9999)        		plotT2.setColor(Color.gray);
        else         plotT2.setColor(t2Moy);
        plotT2.setLineWidth(3);
		plotT2.replace(incrT2++, "line",new double[]{0,(double)(paramsTimelapseT2[tCor][1])},new double[]{maxPlotYT2*0.9,maxPlotYT2*0.9});//Afficher le T2
        plotT2.setLineWidth(1);
        plotT2.setColor(crossColor);
		for(int t=0;t<this.t2Times[tCor].length;t++)plotT2.replace(incrT2++, "line",new double[]{this.t2Times[tCor][t],this.t2Times[tCor][t]},new double[]{dataTimelapseT2[tCor][t]-dataTimelapseT2Sigmas[tCor][t],dataTimelapseT2[tCor][t]+dataTimelapseT2Sigmas[tCor][t]});//Afficher le T2

		
        //PARAM T2 BI
        if(jitterT2Bicomp[tCor]>=9999)        		plotT2.setColor(Color.gray);
        else         plotT2.setColor(t2Bi);
        plotT2.setLineWidth(3);
		plotT2.replace(incrT2++, "line",new double[]{0,paramsTimelapseT2[tCor][3]},new double[]{maxPlotYT2*0.83,maxPlotYT2*0.83});//Afficher le T2

        if(jitterT2Bicomp[tCor]>=9999)        		plotT2.setColor(Color.gray);
        else         plotT2.setColor(t2Bi);
        plotT2.setLineWidth(3);
		plotT2.replace(incrT2++, "line",new double[]{0,paramsTimelapseT2[tCor][5]},new double[]{maxPlotYT2*0.79,maxPlotYT2*0.79});//Afficher le T2
        plotT2.setLineWidth(1);
 
        //LEGENDE
		plotT2.setLineWidth(1);
		String strLegendT2="";
		for(int hh=0;hh<1;hh++) strLegendT2+="\n";
		strLegendT2+="\nNoise +-sigma\nMRI T2 relaxation"+"\n"+tabNamesT2[0]+"\n"+tabNamesT2[1];
		if(this.autoSizingTimePlots) {
			maxPlotYT2=VitimageUtils.max(dataTimelapseT2[tCor]);
			maxPlotYT2*=1.3;
		}
		plotT2.setLimits(0, maxT2, 0,maxPlotYT2);
		plotT2.setColor(new Color(150,150 ,150) );
		plotT2.addLegend(strLegendT2,"bottom-right");
		plotCan2.setPlot(plotT2);
	}
		
	public void actualizeSecondPlots() {
		Timer timer=new Timer();
		//washTheCurves
		for(int c=0;c<maxCurves;c++) {
		plotT21.replace(c, "line",new double[] {-200,-199},new double[] {-100,-100});
			plotT22.replace(c, "line",new double[] {-200,-199},new double[] {-100,-100});
		}
		double [][]compColMin=new double[][] {{255,0,0},{0,215,0},{0,0,255}};
		int [][]compColMinCentral=new int[][] {{90,0,0},{0,90,0},{0,0,90}};
		
		int incrT1=0;
		int incrT2=0;
		double y0T1;		double y1T1;								double y0T2;		double y1T2;	
		double x0T2;double x1T2; 
		double x0T1;double x1T1; 
		
        //HORIZONTAL GRID T21 ET T22
   		for(int tim=0;tim<this.nTimepoints;tim++) {
			plotT21.setLineWidth(4);
			plotT21.setColor(new Color(70,70,70));
			plotT21.replace(incrT1++, "line",new double[] {t0T1,t1T1},new double[] {tim,tim});
			plotT22.setLineWidth(4);
			plotT22.setColor(new Color(70,70,70));
			plotT22.replace(incrT2++, "line",new double[] {t0T2,t1T2},new double[] {tim,tim});
   		}

        //VERTICAL GRID T21
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT21.setLineWidth(1);        
			plotT21.setColor(new Color(150,150,150));
        	for(int mul=2;mul<=9;mul++)plotT21.replace(incrT1++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.nTimepoints});
	        plotT21.setLineWidth(3);        
			plotT21.setColor(new Color(250,250,250));
        	plotT21.replace(incrT1++, "line",new double[] {tt,tt},new double[] {0,this.nTimepoints});
	        plotT21.setLineWidth(1);        
        	plotT21.replace(incrT1++, "line",new double[] {tt*3,tt*3},new double[] {0,this.nTimepoints});
        }
        //VERTICAL GRID T22
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT22.setLineWidth(1);        
			plotT22.setColor(new Color(150,150,150));
			for(int mul=2;mul<=9;mul++)plotT22.replace(incrT2++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.nTimepoints});
	        plotT22.setLineWidth(3);        
			plotT22.setColor(new Color(250,250,250));
        	plotT22.replace(incrT2++, "line",new double[] {tt,tt},new double[] {0,this.nTimepoints});
	        plotT22.setLineWidth(1);        
           	plotT22.replace(incrT2++, "line",new double[] {tt*3,tt*3},new double[] {0,this.nTimepoints});
        }

         //Draw the suns
		y0T1=tCor+1-0.95;		y1T1=tCor+1-0.05;		
		x0T1=t0T1*1.15;		    x1T1=t1T1*0.92;
		x0T2=t0T2*1.15;		    x1T2=t1T2*0.92;
		y0T2=tCor+1-0.95;		y1T2=tCor+1-0.05;		
		plotT21.setLineWidth(30);
		plotT21.setColor(bubbleColor);
//		plotT21.replace(incrT1++, "line",new double[] {x0T1,x1T1},new double[] {y0T1,y0T1});
//        plotT21.replace(incrT1++, "line",new double[] {x0T1,x1T1},new double[] {y1T1,y1T1});
        plotT21.replace(incrT1++, "line",new double[] {x0T1,x0T1},new double[] {tCor+0.2,tCor+0.2});
//        plotT21.replace(incrT1++, "line",new double[] {x1T1,x1T1},new double[] {y0T1/2+y1T1/2,y0T1/2+y1T1/2});

        plotT22.setLineWidth(25);
		plotT22.setColor(new Color(220,220,0));
//        plotT22.replace(incrT2++, "line",new double[] {x0T2,x1T2},new double[] {y0T2,y0T2});
//        plotT22.replace(incrT2++, "line",new double[] {x0T2,x1T2},new double[] {y1T2,y1T2});
        plotT22.replace(incrT2++, "line",new double[] {x0T2,x0T2},new double[] {tCor+0.2,tCor+0.2});
  //      plotT22.replace(incrT2++, "line",new double[] {x1T2,x1T2},new double[] {y0T2/2+y1T2/2,y0T2/2+y1T2/2});
        
   		for(int tim=0;tim<this.nTimepoints;tim++) {
 			double maxValT1=standardCapillaryM0;
			double maxValT2Mono=standardCapillaryM0;
	 		double quota;
	 		double radius=0;
	 		boolean flagGray;
	 		int rad;
	 		double factRad=0.005;

	 		//Draw markers T1
			compColMin=new double[][] {{255,0,0},{0,215,0},{0,0,255}};
			compColMinCentral=new int[][] {{0,220,0},{220,0,0},{0,0,90}};	 		
	 		quota=0.7;flagGray=false;
	 		radius=0.05+0.5*paramsTimelapseT1[tim][0]/maxValT1;if(radius>1)radius=1;if(radius<0.05)radius=0.05;	 		rad=1+(int)(Math.round((radius-0.05)/0.015));plotT21.setLineWidth(rad);
	 		if(jitterT1[tim]>40) flagGray=true;
			if(!flagGray)plotT21.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
	 		else plotT21.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));
//	 		plotT21.setLineWidth((rad));
//			plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim][1],paramsTimelapseT1[tim][1]},new double[] {tim+0.02+0.01*rad,tim-0.01*rad});//Afficher la courbe
	 		plotT21.setLineWidth(6);
//	 		plotT21.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
			plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim][1],paramsTimelapseT1[tim][1]},new double[] {tim,tim+radius});//Afficher la courbe
	/*		if(this.nRepetMonteCarlo>0) {
				plotT21.setLineWidth(2);
		 		plotT21.setColor(new Color(255,255,255));
				plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim][1]-paramsTimelapseT1Sigmas[tim][1],paramsTimelapseT1[tim][1]+paramsTimelapseT1Sigmas[tim][1]},new double[] {tim-0.5,tim-0.5});//Afficher la courbe
			}
		*/	

			
			
			//Draw markers T2
 			compColMin=new double[][] {{0,215,0},{255,0,0},{0,0,255}};
			compColMinCentral=new int[][] {{220,0,0},{90,0,0},{0,0,90}};
	 		quota=0.7;flagGray=false;
			if(khi2T2Mono[tim]<khi2T2Bicomp[tim]) {
		 		if(jitterT2Mono[tim]>40) flagGray=true;
				if(!flagGray)plotT22.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
		 		else plotT22.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));

				radius=0.05+0.5*paramsTimelapseT2[tim][0]/maxValT2Mono;if(radius>1.0)radius=1.0;if(radius<0.05)radius=0.05;	 		rad=2+(int)(Math.round((radius-0.05)/0.015));plotT22.setLineWidth(rad);
//		 		plotT22.setLineWidth((rad));
//				plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][1],paramsTimelapseT2[tim][1]},new double[] {tim+0.1+Math.min(0.5,factRad*rad),tim-0.1-Math.min(0.5,factRad*rad)});//Afficher la courbe
		 		plotT22.setLineWidth(6);
//		 		plotT22.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
				plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][1],paramsTimelapseT2[tim][1]},new double[] {tim,tim+radius});//Afficher la courbe
	/*			if(this.nRepetMonteCarlo>0) {
			 		plotT21.setLineWidth(2);
			 		plotT21.setColor(new Color(255,255,255));
					plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim][1]-paramsTimelapseT2Sigmas[tim][1],paramsTimelapseT2[tim][1]+paramsTimelapseT2Sigmas[tim][1]},new double[] {tim-0.5,tim-0.5});//Afficher la courbe
				}*/
			}
			else {
		 		if(jitterT2Bicomp[tim]>40) flagGray=true;
				if(!flagGray)plotT22.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
		 		else plotT22.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));

				radius=0.05+0.5*paramsTimelapseT2[tim][2]/maxValT2Mono;if(radius>1.0)radius=1.0;if(radius<0.05)radius=0.05;	 		rad=2+(int)(Math.round((radius-0.05)/0.015));plotT22.setLineWidth(rad);
//		 		plotT22.setLineWidth((rad));
//				plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][3],paramsTimelapseT2[tim][3]},new double[] {tim+0.1+Math.min(0.5,factRad*rad),tim-0.1-Math.min(0.5,factRad*rad)});//Afficher la courbe
		 		plotT22.setLineWidth(6);
//		 		plotT22.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
				plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][3],paramsTimelapseT2[tim][3]},new double[] {tim,tim+radius});//Afficher la courbe
		/*		if(this.nRepetMonteCarlo>0) {
			 		plotT21.setLineWidth(2);
			 		plotT21.setColor(new Color(255,255,255));
					plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim][3]-paramsTimelapseT2Sigmas[tim][3],paramsTimelapseT2[tim][3]+paramsTimelapseT2Sigmas[tim][3]},new double[] {tim-0.5,tim-0.5});//Afficher la courbe
				}*/
		 		quota=0.7;flagGray=false;
		 		radius=0.05+0.5*paramsTimelapseT2[tim][4]/maxValT2Mono;if(radius>1.0)radius=1.0;if(radius<0.05)radius=0.05;	 		rad=2+(int)(Math.round((radius-0.05)/0.015));plotT22.setLineWidth(rad);
		 		if(jitterT2Bicomp[tim]>30) flagGray=true;
				if(!flagGray)plotT22.setColor(new Color((int)(quota*compColMin[0][0]),(int)(quota*compColMin[0][1]),(int)(quota*compColMin[0][2])));
		 		else plotT22.setColor(new Color(120+(int)(0.11*compColMin[0][0]),120+(int)(0.11*compColMin[0][1]),120+(int)(0.11*compColMin[0][2])));
//		 		plotT22.setLineWidth((rad));
//				plotT22.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim][5],paramsTimelapseT2[tim][5]},new double[] {tim+0.1+Math.min(0.5,factRad*rad),tim-0.1-Math.min(0.5,factRad*rad)});//Afficher la courbe
		 		plotT22.setLineWidth(6);
//		 		plotT22.setColor(new Color(compColMinCentral[0][0],compColMinCentral[0][1],compColMinCentral[0][2]) );
				plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][5],paramsTimelapseT2[tim][5]},new double[] {tim,tim+radius});//Afficher la courbe
			/*	if(this.nRepetMonteCarlo>0) {
			 		plotT21.setLineWidth(2);
			 		plotT21.setColor(new Color(255,255,255));
					plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT2[tim][5]-paramsTimelapseT2Sigmas[tim][5],paramsTimelapseT2[tim][5]+paramsTimelapseT2Sigmas[tim][5]},new double[] {tim-0.25,tim-0.25});//Afficher la courbe
				}*/	
			}	
	 		plotT21.setLineWidth(2);
	 		plotT21.setColor(new Color(255,0,0));
	 		plotT21.replace(incrT1++, "line", valsSpectrum[tim][2], valsSpectrum[tim][0]);
	 		plotT22.setLineWidth(2);
	 		plotT22.setColor(new Color(0,255,0));
	 		plotT22.replace(incrT2++, "line", valsSpectrum[tim][2], valsSpectrum[tim][1]);
   			//if needed, draw ranging interval

	 		if(this.spectrumRangingModeT1) {
	 			plotT21.setColor(new Color(0,0,255));	 		
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[0]}, new double[] {tim+1-0.7,tim+1-0.3});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[2],this.rangingBoundaries[2]}, new double[] {tim+1-0.7,tim+1-0.3});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[1],this.rangingBoundaries[1]}, new double[] {tim+1-0.8,tim+1-0.2});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[2]}, new double[] {tim+1-0.5,tim+1-0.5});
	 		}   			
	 		if(this.spectrumRangingModeT2) {
	 			plotT22.setColor(new Color(0,0,255));	 		
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[0]}, new double[] {tim+1-0.7,tim+1-0.3});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[2],this.rangingBoundaries[2]}, new double[] {tim+1-0.7,tim+1-0.3});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[1],this.rangingBoundaries[1]}, new double[] {tim+1-0.8,tim+1-0.2});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[2]}, new double[] {tim+1-0.5,tim+1-0.5});
	 		}   			
   		}				
		plotT21.setLimits(t0T1,t1T1, 0, nTimepoints);
		plotT22.setLimits(t0T2,t1T2, 0, nTimepoints);
	}				

	public void actualizeSpectrumCurves() {
		for(int tim=0;tim<this.nTimepoints;tim++) {
			this.valsSpectrum[tim]=computeSpectrumCurve(tim);
		}
		if(!this.separateNormalizationSpectrum) {
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double maxValT1=VitimageUtils.max(this.valsSpectrum[tim][0]);
				double maxValT2=VitimageUtils.max(this.valsSpectrum[tim][1]);
				for(int i=0;i<this.valsSpectrum[tim][0].length;i++)this.valsSpectrum[tim][0][i]=tim+this.valsSpectrum[tim][0][i]/maxValT1*normValueSpectrum;
				for(int i=0;i<this.valsSpectrum[tim][1].length;i++)this.valsSpectrum[tim][1][i]=tim+this.valsSpectrum[tim][1][i]/maxValT2*normValueSpectrum;
			}
		}
		else {
			double maxTotT1=0;
			double maxTotT2=0;
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double maxValT1=VitimageUtils.max(this.valsSpectrum[tim][0]);
				double maxValT2=VitimageUtils.max(this.valsSpectrum[tim][1]);
				if(maxTotT1<maxValT1)maxTotT1=maxValT1;
				if(maxTotT2<maxValT2)maxTotT2=maxValT2;
			}
			for(int tim=0;tim<this.nTimepoints;tim++) {
				for(int i=0;i<this.valsSpectrum[tim][0].length;i++)this.valsSpectrum[tim][0][i]=tim+this.valsSpectrum[tim][0][i]/maxTotT1*normValueSpectrum;
				for(int i=0;i<this.valsSpectrum[tim][1].length;i++)this.valsSpectrum[tim][1][i]=tim+this.valsSpectrum[tim][1][i]/maxTotT2*normValueSpectrum;
			}
		}
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
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(t+1)+",";
		rSentence+=""+nTimepoints+")";
		rSentence+=", acquisition_day=c(0,";
		for(int t=1;t<nTimepoints-1;t++)rSentence+=""+hyperMRIT1T2.actualDay[t]+",";
		rSentence+=""+hyperMRIT1T2.actualDay[nTimepoints-1]+")";
		for(int ec=0;ec<t1Times.length;ec++) {
			rSentence+=", Echo"+ec+"=c("+t1Times[ec]+",";
			for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(dataTimelapseT1[t][ec])+",";
			rSentence+=""+(dataTimelapseT1[nTimepoints-1][ec])+")";
		}
		rSentence+=")\n";
		rSentence+="t2TimesAndObservations <- data.frame(acquisition_number=c(0,";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(t+1)+",";
		rSentence+=""+nTimepoints+")";
		rSentence+=", acquisition_day=c(0,";
		for(int t=1;t<nTimepoints-1;t++)rSentence+=""+hyperMRIT1T2.actualDay[t]+",";
		rSentence+=""+hyperMRIT1T2.actualDay[nTimepoints-1]+")";
		for(int ec=0;ec<t2Times.length;ec++) {
			rSentence+=", Echo"+ec+"=c("+t2Times[ec]+",";
			for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(dataTimelapseT2[t][ec])+",";
			rSentence+=""+(dataTimelapseT2[nTimepoints-1][ec])+")";
		}
		rSentence+=")\n";
	
	
		
		rSentence+="estimatedParameters <- data.frame(acquisition_number=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(t+1)+",";
		rSentence+=""+nTimepoints+"),acquisition_day=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+hyperMRIT1T2.actualDay[t]+",";
		rSentence+=""+hyperMRIT1T2.actualDay[nTimepoints-1]+"), M0T1=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT1[t][0])+",";
		rSentence+=""+""+(paramsTimelapseT1[nTimepoints-1][0])+"), T1=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT1[t][1])+",";
		rSentence+=""+""+(paramsTimelapseT1[nTimepoints-1][1])+"), M0T2=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT2[t][0])+",";
		rSentence+=""+""+(paramsTimelapseT2[nTimepoints-1][0])+"), T2=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT2[t][1])+",";
		rSentence+=""+""+(paramsTimelapseT2[nTimepoints-1][1])+"), M0T21=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT2[t][2])+",";
		rSentence+=""+""+(paramsTimelapseT2[nTimepoints-1][2])+"), T21=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT2[t][3])+",";
		rSentence+=""+""+(paramsTimelapseT2[nTimepoints-1][3])+"), M0T22=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT2[t][4])+",";
		rSentence+=""+""+(paramsTimelapseT2[nTimepoints-1][4])+"), T22=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(paramsTimelapseT2[t][5])+",";
		rSentence+=""+""+(paramsTimelapseT2[nTimepoints-1][5])+"), JitterT1=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(jitterT1[t])+",";
		rSentence+=""+""+(jitterT1[nTimepoints-1])+"), JitterT2MonoExp=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(jitterT2Mono[t])+",";
		rSentence+=""+""+(jitterT2Mono[nTimepoints-1])+"), JitterT2BiExp=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(jitterT2Bicomp[t])+",";
		rSentence+=""+""+(jitterT2Bicomp[nTimepoints-1])+"))\n\n";
	
	
		
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
		for(int t=1;t<=nTimepoints;t++) {
			pySentence+=",\n["+t+","+hyperMRIT1T2.actualDay[t-1];
			for(int ec=0;ec<t1Times.length;ec++) pySentence+=","+dataTimelapseT1[t-1][ec];
			pySentence+="]";
		}
		pySentence+="\n])\n";
	
		pySentence+="t2TimesAndObservations=np.array([[0,0";
		for(int ec=0;ec<t2Times.length;ec++)pySentence+=","+t2Times[ec];
		pySentence+="]";					
		for(int t=1;t<=nTimepoints;t++) {
			pySentence+=",\n["+t+","+hyperMRIT1T2.actualDay[t-1];
			for(int ec=0;ec<t2Times.length;ec++) pySentence+=","+dataTimelapseT2[t-1][ec];
			pySentence+="]";
		}
		pySentence+="\n])\n";
	
		
		pySentence+="estimatedParameters=np.array([\n";
		for(int t=0;t<nTimepoints;t++) {
			pySentence+="["+(t+1)+","+hyperMRIT1T2.actualDay[t];		
			for(int val=0;val<2;val++)pySentence+=","+paramsTimelapseT1[t][val];
			for(int val=0;val<6;val++)pySentence+=","+paramsTimelapseT2[t][val];
			pySentence+=","+jitterT1[t]+","+jitterT2Mono[t]+","+jitterT2Bicomp[t];
			pySentence+="]"+(t==nTimepoints-1 ? "\n" : ",\n");
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
		for(int t=1;t<=nTimepoints;t++) {
			matSentence+=";\n["+t+","+hyperMRIT1T2.actualDay[t-1];
			for(int ec=0;ec<t1Times.length;ec++) matSentence+=","+dataTimelapseT1[t-1][ec];
			matSentence+="]";
		}
		matSentence+="\n]\n";
	
		matSentence+="t2TimesAndObservations=[\n[0,0";
		for(int ec=0;ec<t2Times.length;ec++)matSentence+=","+t2Times[ec];
		matSentence+="]";					
		for(int t=1;t<=nTimepoints;t++) {
			matSentence+=";\n["+t+","+hyperMRIT1T2.actualDay[t-1];
			for(int ec=0;ec<t2Times.length;ec++) matSentence+=","+dataTimelapseT2[t-1][ec];
			matSentence+="]";
		}
		matSentence+="\n]\n";
	
		
		matSentence+="estimatedParameters=[\n";
		for(int t=0;t<nTimepoints;t++) {
			matSentence+="["+(t+1)+","+hyperMRIT1T2.actualDay[t];		
			for(int val=0;val<2;val++)matSentence+=","+paramsTimelapseT1[t][val];
			for(int val=0;val<6;val++)matSentence+=","+paramsTimelapseT2[t][val];
			matSentence+=","+jitterT1[t]+","+jitterT2Mono[t]+","+jitterT2Bicomp[t];
			matSentence+="]"+(t==nTimepoints-1 ? "\n" : ";\n");
		}
		matSentence+="]\n";
	}
	
	public void actualizeDisplayedNumbers() {
		String infoT1 = String.format("T1 parameters estimation : M0=%5.3f ; T1=%6.1f ms ; Khi2=%5.4f",
				paramsTimelapseT1[this.tCor][0],paramsTimelapseT1[this.tCor][1],this.khi2T1[this.tCor]);
		info1.setText(infoT1);

		String infoT2 = String.format("Estimation area=%2dx%2dx%2d-%s     #Points=%d ",
				1+2*this.crossWidth,1+2*this.crossWidth,1+2*this.crossThick,
				(this.gaussianWeighting?"" : ""),this.nPtsCur);
		info2.setText(infoT2);

		String infoT3 = String.format("T2 parameters estimation : M0 (mono)=%5.3f ;  M0 (bi)=(%5.3f | %5.3f)  ; T2 (mono)=%5.1f ms ; T2 (bi)=(%5.1f ms | %5.1f ms) ; Khi2 (mono)=%4.3f ; Khi2 (bi)=%4.3f ",
				paramsTimelapseT2[this.tCor][0],paramsTimelapseT2[this.tCor][2],paramsTimelapseT2[this.tCor][4],
				paramsTimelapseT2[this.tCor][1],paramsTimelapseT2[this.tCor][3],(paramsTimelapseT2[this.tCor][5]>10000 ? 9999 : paramsTimelapseT2[this.tCor][5]),
				this.khi2T2Mono[this.tCor],this.khi2T2Bicomp[this.tCor]);
		info3.setText(infoT3);
	}

	public void displayResultsAgain() {
		Timer tim2=new Timer();
		actualizeFirstPlots();
		actualizeSecondPlots();
		actualizeDisplayedNumbers();
		actualizeCursor();
	}
	
	
	
		
		
	
	
	

	
	
	
	
	


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

	public double[]smoothHisto(double[]histo,int sig){
		if(sig<=0)return histo;
		int nbBins=histo.length;
		double[]histSmooth=new double[nbBins];

		//build gaussian kernel
		double[]kernel=new double[4*sig+1];
		double sum=0,sumWeight=0;
		for(int i=0;i<kernel.length;i++)kernel[i]=Math.abs(2*sig-i);
		for(int i=0;i<kernel.length;i++)kernel[i]=(kernel[i]*kernel[i])/(2*sig*sig);
		for(int i=0;i<kernel.length;i++) {kernel[i]=Math.exp(-kernel[i]*kernel[i])/(sig*Math.sqrt(2*Math.PI));sum+=kernel[i];}
		for(int i=0;i<kernel.length;i++)kernel[i]/=sum;
		sum=0;
		for(int i=0;i<kernel.length;i++)sum+=kernel[i];

		for(int i=0;i<nbBins;i++) {
			sum=0;
			sumWeight=0;
			for(int j=0;j<kernel.length;j++) {
				int index=i+j-2*sig;
				if(index>=0 && index<nbBins) {
					sumWeight+=kernel[j];
					sum+=kernel[j]*histo[index];
				}
			}
			histSmooth[i]=(sum/sumWeight);
		}
		return histSmooth;  
	}
	

	
	
	public void mouseClicked(MouseEvent e) {
		if(imgCan1.cursorOverImage()) currentCanvas=WIN_T1;
		if(imgCan2.cursorOverImage()) currentCanvas=WIN_T1;
		if(plotCan1.cursorOverImage()) currentCanvas=WIN_PLOT1;
		if(plotCan2.cursorOverImage()) currentCanvas=WIN_PLOT2;
		if(plotCan21.cursorOverImage()) currentCanvas=WIN_PLOT21;
		if(plotCan22.cursorOverImage()) currentCanvas=WIN_PLOT22;

		if(currentCanvas==WIN_T1 || currentCanvas==WIN_T2) {
			IJ.log("Click on MR images | Coordinates=("+ xMouse+","+yMouse+")"+"  |  zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			xMouse=e.getX();
			yMouse=e.getY();
			if(currentCanvas==WIN_T1 ) {
				xCor=imgCan1.offScreenX(xMouse);			
				yCor=imgCan1.offScreenY(yMouse);
				xD=imgCan1.offScreenXD(xMouse);
				yD=imgCan1.offScreenYD(yMouse);
			}
			if(currentCanvas==WIN_T2 ) {
				xCor=imgCan2.offScreenX(xMouse);			
				yCor=imgCan2.offScreenY(yMouse);
				xD=imgCan2.offScreenXD(xMouse);
				yD=imgCan2.offScreenYD(yMouse);
			}
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();			
			displayResultsAgain();
			repaintAll();
			return;
		}
		else if(currentCanvas==WIN_PLOT21) {
			IJ.log("Click on Plot21 | Coordinates=("+ xMouse+","+yMouse+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			spectrumRangingModeT1=true;	
			spectrumRangingModeT2=false;	
			this.xMouseRange=e.getX();
			actualizeRangingBoundaries();
			identifyRangedData();			
			displayResultsAgain();
			repaintAll();
			return;
		}

		else if(currentCanvas==WIN_PLOT22) {
			IJ.log("Click on Plot22 | Coordinates=("+ xMouse+","+yMouse+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			spectrumRangingModeT2=true;	
			spectrumRangingModeT1=false;	
			this.xMouseRange=e.getX();
			actualizeRangingBoundaries();
			identifyRangedData();			
			displayResultsAgain();
			repaintAll();
			return;
		}
		else {
			IJ.log("Click on another plot. Ranging mode stop");
			spectrumRangingModeT2=false;	
			spectrumRangingModeT1=false;	
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
			repaintAll();	
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
	

	
	

	
}

