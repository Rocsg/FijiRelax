package com.vitimage.aplimtools;


import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Event;
import java.awt.Font;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import com.jogamp.nativewindow.VisualIDHolder;
import com.vitimage.common.Timer;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

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
import ij.gui.PolygonRoi;
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


public class MRI_HyperCurvesExplorer extends PlugInFrame implements ActionListener,KeyListener, MouseListener, MouseWheelListener, MouseMotionListener {// 
	private double sigmaSmoothingBeforeEstimation=0.75;
	private int xSrcStart;
	private int ySrcStart;
	private boolean roiIsHidden=false;
    int deltaT2=0;
	JButton butNextEcho=new JButton("<Html><br>Next<br>echo<br> <Html>");
	JButton butPrevEcho=new JButton("<Html><br>Prev.<br>echo<br> <Html>");
	JButton butNextTr=new JButton("<Html><br>Next<br>Tr<br> <Html>");
	JButton butPrevTr=new JButton("<Html><br>Prev.<br>Tr<br> <Html>");
	JButton butPlay=new JButton("<Html><br>Play<br>with<br>prms<Html>");
	JButton butAllCurves=new JButton("<Html><br>All<br>data<br> <Html>");

	JButton butM0Map=new JButton("<Html> <br>M0<br>map<Html>");
	JButton butT1Map=new JButton("<Html> <br>T1<br>map<Html>");
	JButton butT2Map=new JButton("<Html> <br>T2<br>map<Html>");

	JTextField textExp=new JTextField("Image exploration");
	JTextArea textTim=new JTextArea("              Observation day");
	JButton butPrevDay=new JButton("Prev. day");
	JButton butNextDay=new JButton("Next day");
	JTextArea textSli=new JTextArea("               Current Z-slice");
	JButton butPrevSli=new JButton("Prev. slice");
	JButton butNextSli=new JButton("Next slice");
	JButton butFireLut=new JButton("Fire LUT");
	JButton butGrayLut=new JButton("Gray LUT");
	JTextArea textLuts=new JTextArea("               Colormap");

	JTextField textRoi=new JTextField("       Region of interest");
	JTextArea textSam=new JTextArea("       Sample size (5x5x1)");
	JButton butPrevRoi=new JButton("Smaller");
	JButton butNextRoi=new JButton("Larger");
	JTextArea textRan=new JTextArea("           Ranging interval");
	JButton butPrevRang=new JButton("Smaller");
	JButton butNextRang=new JButton("Larger");
	JTextArea textMisc=new JTextArea("           Miscellaneous");
	JButton butHide=new JButton("<Html>Show<br>ranging<Html>");
	JButton butRoi=new JButton("<Html>Custom<br>ROI<Html>");

	JTextArea textInfoEchoes=new JTextArea("-------T2SEQ Day 28 TR=10000 ms TE=11ms------",45,1);
	JTextArea textInfoMaps=new JTextArea("--------M0MAP Day 28 -----------",45,1);

	private boolean playMode=false;
	private int switchT1T2=1;
	private double initMag=0;
	private double xD;
	private double yD;
	private ImageJ ij;
	private HyperMRIT1T2 hyperMRIT1T2;
	private int[] dims;
	private double[] voxs;
	private int nTimepoints;
	private int cEchoes;
	private int cMaps;
	private int cT1T2TrIndex;
	private int cT1T2TeIndex;
	private int nT1T2TrAvailable;
	private int [] nT1T2TeAvailable;
	private double[][] dataTimelapseT1;
	private double[][] dataTimelapseT2;
	private double[][] dataTimelapseT1T2;
	private double arbitraryLimitBetweenT1andT2=500;
	private final double normValueSpectrum=0.9;
	private ImagePlus mapsImage;
	private ImagePlus echoesImage;
	Color bubbleColor=new Color(220,220,0);
	Color paramsUnactivated=new Color(150,150,150);
	Color paramsT1ActivatedColor=new Color(255,200,200);
	Color paramsT1TextActivatedColor=new Color(255,220,220);
	Color paramsT2ActivatedColor=new Color(200,255,200);
	Color paramsT2TextActivatedColor=new Color(220,255,220);
	Color crossColor=new Color (60,60,255);
	Color curveT1Mono=new Color (170,0,0);
	Color curveT1Bicomp=new Color (255,70,70);
	Color curveT2Mono=new Color (0,115,0);
	Color curveT2Bicomp=new Color (50,205,50);

	private static final double maxAcceptableM0=2*MRUtils.maxM0ForNormalization;
	private static final double minAcceptableT1=10;
	private static final double minAcceptableT2=8;
	private static final double maxAcceptableT1=10000;
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
	private static final int WIN_PLOT1T2=5;
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
	double[][][]sigmaDataRiceLookupTableT1T2;
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
	double meanNoiseT1T2Cur=1;
	double sigmaNoiseT1T2Cur=1;
	int nPtsCur=1;
	int thickCurve=1;

	
	
	//Gui parameters : Window constants
	public int xMouseRange=0;
	private int TARGET_WIDTH_FOR_PLOTS_T1T2;
	private int TARGET_WIDTH_FOR_PLOTS_ALONE;
	private int TARGET_HEIGHT_FOR_PLOTS_AND_IMGS;
	private int TARGET_WIDTH_FOR_PLOTS_AND_IMGS;
	private int DX_IMG;
	private int DY_IMG;
	private int DY_TEXT;
	private int DELTA_X;
	private int DELTA_Y;
	private int DX_PLOT_1;
	private int DX_PLOT_2;
	private int DY_PLOT;
	private int totalSizeX;  
	private int totalSizeY;  
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
	int crossWidth=0;//3;
	int crossThick=0;
	int xMouse=0;
	int yMouse=0;
	int xCor=1;
	int yCor=1;
	int tCor=0;
	int zCor=0;
	private Plot plotT1;
	private Plot plotT2;
	private Plot plotT1T2;
	private Plot plotT21;
	private Plot plotT22;
	private ImageCanvas imgCan1;
	private ImageCanvas imgCan2;
	private PlotCanvas plotCan1;
	private PlotCanvas plotCan2;
	private PlotCanvas plotCan21;
	private PlotCanvas plotCan22;
	private double zoomLevel=0;
	private double maxPlotYT1;
	private double maxPlotYT2;
	private int currentCanvas=1;
	JTextField []titles=new JTextField[8];
	JTextField [][]params=new JTextField[16][5];
	static PlugInFrame instance;
	boolean computeMultiComp=true;
	private static final long serialVersionUID = 1L;
	double xStepCuteT1=50;
	double xStepCuteT2=3;
	double [][]t1Times;
	double [][]t2Times;
	double [][][]t1t2Times;
	double[]timesT1Cute;
	double[]timesT2Cute;
	double[][][][][]timesT1T2Cute;

	double[][][][]dataTimelapseFull;
	int [][]pointsCoords;
	double[][][][]pointsEstimations;
	double[][]dataTimelapseT1Sigmas;
	double[][]dataTimelapseT2Sigmas;
	double[][]dataTimelapseT1T2Sigmas;

	double[][]paramsTimelapseT1;
	double[][]paramsTimelapseT2;
	double[][]paramsTimelapseT1T2;
	double[][][]valsSpectrum;

	//Updated just before plot being updated
	double [][]tabFittenT1Mono;
	double [][]tabFittenT1MonoCute;
	double[]jitterT1Mono;
	double[]khi2T1Mono;
	double[]pValT1Mono;

	double [][]tabFittenT1Bicomp;
	double [][]tabFittenT1BicompCute;
	double[]jitterT1Bicomp;
	double[]khi2T1Bicomp;
	double[]pValT1Bicomp;

	double [][]tabFittenT2Mono;
	double [][]tabFittenT2MonoCute;
	double[]jitterT2Mono;
	double[]khi2T2Mono;
	double[]pValT2Mono;

	double [][]tabFittenT2Bicomp;
	double [][]tabFittenT2BicompCute;
	double[]jitterT2Bicomp;
	double[]khi2T2Bicomp;	
	double[]pValT2Bicomp;

	double [][]tabFittenT1T2Mono;
	double [][][]tabFittenT1T2MonoCute;
	double[]jitterT1T2Mono;
	double[]khi2T1T2Mono;	
	double[]pValT1T2Mono;

	double [][]tabFittenT1T2Bicomp;
	double [][][]tabFittenT1T2BicompCute;
	double[]jitterT1T2Bicomp;
	double[]khi2T1T2Bicomp;	
	double[]pValT1T2Bicomp;

	double [][]tabFittenT1BiT2Bicomp;
	double [][][]tabFittenT1BiT2BicompCute;
	double[]jitterT1BiT2Bicomp;
	double[]khi2T1BiT2Bicomp;	
	double[]pValT1BiT2Bicomp;

	double [][]tabFittenT1T2Bionano;
	double [][][]tabFittenT1T2BionanoCute;
	double[]jitterT1T2Bionano;
	double[]khi2T1T2Bionano;	
	double[]pValT1T2Bionano;

	double statusRoi;
	int[][]correspondanceCanvas;
	private ArrayList<int[]> rangeRoiPoints;
	private int lastCount21=0;
	private int lastCount22=0;
	private int lastCount1=0;
	private int lastCount2=0;
	private int yMouseRange;
	private String[][][] mapsText;
	private String[][][] echoesText;
	private int xMouseRangeCurve;
	private int yMouseRangeCurve;
	private int xMouseDrag;
	private int yMouseDrag;
	private int xMouseStart;
	private int yMouseStart;
	private String file;
	private Roi userRoiInitial;
	private boolean T1T2MixedEstimation=false;
	private int visualFactorT2=50;
	private double maxT1T2OnPlot=0;
	private double[][][] t1t2TrTimes;
	private double[][][] t1t2TeTimes;
	private PlotCanvas plotCanT1T2;
	private double maxPlotYT1T2;
	private double[][] timesT1T2BubbleCute;
	private int lastCountT1T2;
	private boolean[] isSelectedCurve;
	private boolean[] isSelectedIndex;
	private int TARGET_HEIGHT_FOR_UP_PLOTS_AND_IMGS;
	private int TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS;
	private JPanel[][] paramsPanelRightTX;
	private boolean switchAllCurves=false;
	private int memSwitch;
	private int deltaSigma=0;

	private static boolean testT1T2=true;
	
	//TODO : mono R bi G 
	//TODO : Activate and display parameters
	
	/** Entry points and startup functions*/
	public static void main(String[]args) {
		if(!testT1T2)runExplorer("/home/fernandr/Bureau/Traitements/Sorgho/Tests/Test SSM1/Output_normalization/Normalized Hyper MRI.tif");
		else {
			boolean sorghoTime=false;
			if(sorghoTime) {
				runExplorer("/home/fernandr/Bureau/Traitements/Sorgho/RegBM1/BM1_Combined.tif");
				return;
			}
			
			
			
			boolean bouture=false;
			if(bouture) {
				runExplorer("/home/fernandr/Bureau/Recherche_diff/B099J0/Imported2/SL_RAW_hyperimage_build_time_2020-04-09_05-00.tif");
				return;
			}
			else {
				String spec="BM1"; 
				
				int nbUsed=0; //
				File f=new File("/home/fernandr/Bureau/Traitements/Sorgho/Cartes_rangees_par_specimen/"+spec);
				String []allImages=VitimageUtils.stringArraySort(f.list());
				if(nbUsed>allImages.length-1)nbUsed=allImages.length-1;
				File img=new File(f,allImages[nbUsed]);
				runExplorer(img.getAbsolutePath());
			}
		}
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
		if(hyperMRIT1T2.hasT1T2sequence)mrExplo.T1T2MixedEstimation=true;
		
		Object[]obj=hyperMRIT1T2.getMapsImage(false);
		mrExplo.mapsImage=(ImagePlus) obj[0];
		mrExplo.mapsText=(String[][][]) obj[1];
		obj=hyperMRIT1T2.getEchoesImage(false);
		mrExplo.echoesImage=(ImagePlus) obj[0];
		mrExplo.echoesText=(String[][][]) obj[1];
		mrExplo.hyperMRIT1T2=hyperMRIT1T2;
		mrExplo.dims=hyperMRIT1T2.dims;
		mrExplo.voxs=hyperMRIT1T2.voxs;
		if(mrExplo.ij==null)mrExplo.ij = new ImageJ();
		else mrExplo.ij = IJ.getInstance();
		IJ.log("Starting MRI Curve Explorer");
		mrExplo.xCor=hyperMRIT1T2.dims[0]/2;
		mrExplo.yCor=hyperMRIT1T2.dims[1]/2;
		mrExplo.zCor=0;//hyperMRIT1T2.dims[2]/2;
		mrExplo.tCor=0;
		mrExplo.nTimepoints=hyperMRIT1T2.nTimepoints;
		mrExplo.setupTimesTrTe();
		mrExplo.setupStructures();
		mrExplo.startGui();
	}
	
	public void setupTimesTrTe() {
		if(!T1T2MixedEstimation) {
			//this.t1Times=this.hyperMRIT1T2.getT1TrTimes();
			//this.t2Times=this.hyperMRIT1T2.getT2TeTimes();
		}
		else {
		}
	}
	
	public void setupStructures() {
		initializeScreenConstants();		
		riceEstimator=RiceEstimator.getDefaultRiceEstimatorForNormalizedHyperEchoesT1AndT2Images();

		//Data gathered on a voxel, on the neighbourhood of a voxel, and associated measured sigma values
		this.dataTimelapseT1=new double[this.nTimepoints][];
		this.dataTimelapseT2=new double[this.nTimepoints][];
		this.dataTimelapseT1T2=new double[this.nTimepoints][];
		this.dataTimelapseFull=new double[this.nTimepoints][][][];
		this.dataTimelapseT1Sigmas=new double[this.nTimepoints][];
		this.dataTimelapseT2Sigmas=new double[this.nTimepoints][];
		this.dataTimelapseT1T2Sigmas=new double[this.nTimepoints][];

		//Params estimated from fits
		this.paramsTimelapseT1=new double[this.nTimepoints][12];
		this.paramsTimelapseT2=new double[this.nTimepoints][12];
		this.paramsTimelapseT1T2=new double[this.nTimepoints][18];

		//Corresponding curves and nice curves (for display)
		this.tabFittenT1Mono=new double[this.nTimepoints][];
		this.tabFittenT1Bicomp=new double[this.nTimepoints][];
		this.tabFittenT1MonoCute=new double[this.nTimepoints][];
		this.tabFittenT1BicompCute=new double[this.nTimepoints][];
		this.tabFittenT2Mono=new double[this.nTimepoints][];
		this.tabFittenT2Bicomp=new double[this.nTimepoints][];
		this.tabFittenT2MonoCute=new double[this.nTimepoints][];
		this.tabFittenT2BicompCute=new double[this.nTimepoints][];
		this.tabFittenT1T2Mono=new double[this.nTimepoints][];
		this.tabFittenT1T2Bicomp=new double[this.nTimepoints][];
		this.tabFittenT1BiT2Bicomp=new double[this.nTimepoints][];
		this.tabFittenT1T2Bionano=new double[this.nTimepoints][];
		this.tabFittenT1T2MonoCute=new double[this.nTimepoints][][];
		this.tabFittenT1T2BicompCute=new double[this.nTimepoints][][];
		this.tabFittenT1BiT2BicompCute=new double[this.nTimepoints][][];
		this.tabFittenT1T2BionanoCute=new double[this.nTimepoints][][];
		this.timesT1Cute=MRUtils.getProportionalTimes(0,maxT1,xStepCuteT1);
		this.timesT2Cute=MRUtils.getProportionalTimes(0,maxT2,xStepCuteT2);
		if(T1T2MixedEstimation) {
			this.timesT1T2Cute=getCuteTimesT1T2();
		}
	
		//Estimation errors
		this.jitterT1Mono=new double[this.nTimepoints];
		this.jitterT1Bicomp=new double[this.nTimepoints];
		this.jitterT2Mono=new double[this.nTimepoints];
		this.jitterT2Bicomp=new double[this.nTimepoints];
		this.jitterT1T2Mono=new double[this.nTimepoints];
		this.jitterT1T2Bicomp=new double[this.nTimepoints];
		this.jitterT1BiT2Bicomp=new double[this.nTimepoints];
		this.jitterT1T2Bionano=new double[this.nTimepoints];

		//Khi and pvalue for estimation
		this.khi2T1Mono=new double[this.nTimepoints];
		this.khi2T1Bicomp=new double[this.nTimepoints];
		this.khi2T2Mono=new double[this.nTimepoints];
		this.khi2T2Bicomp=new double[this.nTimepoints];
		this.khi2T1T2Mono=new double[this.nTimepoints];
		this.khi2T1T2Bicomp=new double[this.nTimepoints];
		this.khi2T1BiT2Bicomp=new double[this.nTimepoints];
		this.khi2T1T2Bionano=new double[this.nTimepoints];

		this.pValT1Mono=new double[this.nTimepoints];
		this.pValT1Bicomp=new double[this.nTimepoints];
		this.pValT2Mono=new double[this.nTimepoints];
		this.pValT2Bicomp=new double[this.nTimepoints];
		this.pValT1T2Mono=new double[this.nTimepoints];
		this.pValT1T2Bicomp=new double[this.nTimepoints];
		this.pValT1BiT2Bicomp=new double[this.nTimepoints];
		this.pValT1T2Bionano=new double[this.nTimepoints];
		
		//Computed distribution of T1 and T2 values around the neighbourhood of a voxel
		this.correspondanceCanvas=new int[1000][2];
		this.valsSpectrum=new double[this.nTimepoints][][];
		this.visualFactorT2=50;
	}
	

	public double[][][][][]getCuteTimesT1T2(){
		double[][][][]t1t2trte=hyperMRIT1T2.getT1T2TrTeTimes();// [Time][Curve (T1, puis T2 succ][Choice : Tr, Te or Visual T][actual values]
		this.t1t2TrTimes=this.hyperMRIT1T2.getT1T2TrTimes();
		this.t1t2TeTimes=this.hyperMRIT1T2.getT1T2TeTimes();
		this.t1t2Times=new double[nTimepoints][dims[2]][];
		double[][][][][]tabRet=new double[nTimepoints][dims[2]][][][];
		for(int t=0;t<nTimepoints;t++) {
			for(int z=0;z<dims[2];z++) {
				//Collect all the possible transversal relaxations curves, that are successive points having the same Tr, but with varying Te
				int nt1=0;
				double memTr=-1;double memTe=-1;
				boolean isUp=false;
				this.isSelectedCurve=new boolean[t1t2trte[t][z].length];//Enough space in fact, whereas it is not exactly the good number, but it is greater...
				for(int c=0;c<t1t2trte[t][z].length;c++) {
					//System.out.println("Processing c="+c+" : T1="+t1t2trte[t][z][c][0]);
					if(t1t2trte[t][z][c][0]!=memTr) {
						//System.out.println(" -->Is different.");
						isSelectedCurve[nt1]=isUp;isUp=!isUp;memTr=t1t2trte[t][z][c][0];nt1++;
					}
					if(t1t2trte[t][z][c][0]>9000) {
						t1t2trte[t][z][c][1]+=deltaT2;
						t1t2TeTimes[t][z][c]+=deltaT2;
					}
				}
				isSelectedCurve[nt1]=true;
				double[][]valuesTrTe=new double[nt1][2];
				tabRet[t][z]=new double[nt1+1][3][];
				this.isSelectedIndex=new boolean[t1t2trte[t][0].length];
				this.t1t2Times[t][z]=new double[t1t2trte[t][z].length];
				nt1=-1;memTr=-1;
			
				for(int c=0;c<t1t2trte[t][z].length;c++) {
					t1t2Times[t][z][c]=t1t2trte[t][z][c][0] + visualFactorT2 * t1t2trte[t][z][c][1];

					if(t1t2trte[t][z][c][0]!=memTr) {memTr=t1t2trte[t][z][c][0];nt1++;valuesTrTe[nt1][0]=memTr;}
					isSelectedIndex[c]=isSelectedCurve[nt1+1];
					valuesTrTe[nt1][1]=t1t2trte[t][z][c][1];// Fatally, this will be the bigger (last) Te at the end
				}
		//		System.out.println("Bilan des curves :");
				for(int i=0;i<valuesTrTe.length;i++) {
			//		System.out.println("Curve "+i+" : Tr="+valuesTrTe[i][0]+" , Temax="+valuesTrTe[i][1]);
				}
				
				//Add the first curve : the longitudinal relaxation, constituted of succesive growing Tr from 0 to the max, that is memTr, with Te=0
				int nPts=(int)Math.ceil(valuesTrTe[valuesTrTe.length-1][0]*1.0/xStepCuteT1+1);
				tabRet[t][z][0]=new double[3][nPts];
				for(int pt=0;pt<nPts;pt++) {
					tabRet[t][z][0][0][pt]=pt*xStepCuteT1;
					tabRet[t][z][0][1][pt]=0;
				}
			
				//Add the other curves : the transversal relaxations, constituted of successive growing Tr curves, with for each Te going from 0 to the max Te at this level
				for(int cu=0;cu<valuesTrTe.length;cu++) {
					nPts=(int)Math.ceil(valuesTrTe[cu][1]*1.0/xStepCuteT2+5);
					tabRet[t][z][cu+1]=new double[3][nPts];
					for(int pt=0;pt<nPts;pt++) {
						tabRet[t][z][cu+1][1][pt]=pt*xStepCuteT2+0.0001;
						tabRet[t][z][cu+1][0][pt]=valuesTrTe[cu][0];
					}
				}
				for(int i=0;i<tabRet[t][z].length;i++) {	
					System.out.println("\nIs selected curve ? "+isSelectedCurve[i]);
					for(int j=0;j<tabRet[t][z][i][0].length;j++) {
						tabRet[t][z][i][2][j]=   tabRet[t][z][i][0][j]     +   visualFactorT2 * tabRet[t][z][i][1][j];
						//System.out.println("Converting "+tabRet[t][i][0][j]+" and "+tabRet[t][i][1][j]+" in "+tabRet[t][i][2][j]);
						if(tabRet[t][z][i][2][j]>maxT1T2OnPlot)maxT1T2OnPlot=tabRet[t][z][i][2][j];
					}
					if(z==0 && t==0) {
						System.out.println("New nice build :");
						System.out.println("Data Tr="+TransformUtils.stringVectorN(tabRet[t][z][i][0], ""));
						System.out.println("Data Te="+TransformUtils.stringVectorN(tabRet[t][z][i][1], ""));
						System.out.println("Data graphic="+TransformUtils.stringVectorN(tabRet[t][z][i][2], ""));
					}
				}
			}
		}
		//System.exit(0);
		return tabRet;
	}
	

	public double[]t1t2selectCute(double[]dat){
		int n=0;
		for(int i=0;i<dat.length;i++)if(switchAllCurves || isSelectedIndex[i])n++;
		double[]ret=new double[n];n=0;
		for(int i=0;i<dat.length;i++)if(switchAllCurves || isSelectedIndex[i])ret[n++]=dat[i];
		return ret;
	}
	
	public boolean t1t2SelectedCurve(int tran) {
		return (switchAllCurves || isSelectedCurve[tran]);
	}
			
			
	public void startGui() {
	    WindowManager.addWindow(this);
	    instance = this;
		this.imgCan1=new ImageCanvas(this.mapsImage);
		this.imgCan2=new ImageCanvas(this.echoesImage);
        startPlotsAndRoi();
 		initializeGUI();
 		cEchoes=this.echoesImage.getNChannels()/2;
		imgCan1.getImage().setPosition(1,zCor+1,tCor+1);
		imgCan2.getImage().setPosition(cEchoes,zCor+1,tCor+1);
		zoomLevel=imgCan1.getMagnification();
		xMouse=TARGET_WIDTH_FOR_PLOTS_AND_IMGS/2;
		yMouse=TARGET_HEIGHT_FOR_PLOTS_AND_IMGS/2;
		xCor=imgCan1.offScreenX(xMouse);			
		yCor=imgCan1.offScreenY(yMouse);
		xD=imgCan1.offScreenXD(xMouse);
		yD=imgCan1.offScreenYD(yMouse);
		actualizeMriObservationsBasedOnData();
		computeResultsAgain();			
		actualizeCursor();
		displayResultsAgain();
	}


		

	
	
	
	
	
	/** Gui building helpers*/	
	public void repaintAll() {
		imgCan1.repaint();
		imgCan2.repaint();
		if(T1T2MixedEstimation)plotCanT1T2.repaint();
		else{ plotCan1.repaint();plotCan2.repaint();}
		plotCan21.repaint();
		plotCan22.repaint();
	}
		

	public void initializeScreenConstants() {
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		int screenY=(int) Math.round(screenSize.getHeight());
		int screenX=(int) Math.round(screenSize.getWidth());
		if(screenX>1920)screenX=1920;
		IJ.log("Screen resolution : "+screenX+" X "+screenY);
		if(screenX<1054)IJ.showMessage("Your screen has a very low resolution : "+screenX+" X "+screenY+"\nPlease consider investing in one which have at least 1024 lines x 1480 columns.\nPlugin will run in survivor mode, thus everything can happen");
		int DY_TITLE=30;
		TARGET_HEIGHT_FOR_PLOTS_AND_IMGS=400;
		TARGET_HEIGHT_FOR_UP_PLOTS_AND_IMGS=480;
		TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS=300;
		TARGET_WIDTH_FOR_PLOTS_AND_IMGS=400;
		TARGET_WIDTH_FOR_PLOTS_T1T2=1300;
		TARGET_WIDTH_FOR_PLOTS_ALONE=650;
		DX_IMG=TARGET_WIDTH_FOR_PLOTS_AND_IMGS;
		DY_IMG=TARGET_HEIGHT_FOR_PLOTS_AND_IMGS;
		DY_TEXT=30;
		DELTA_X=10;
		DELTA_Y=10;
		totalSizeY= 2*DY_IMG+DY_TITLE;  
		totalSizeX= DX_IMG+DX_PLOT_1+DX_PLOT_2+2*DELTA_X ;  
	}
	

	
	public void initializeGUI() {
		ImageJ ij = IJ.getInstance();
		//Prepare Image T1
		imgCan1.setMagnification(1.0*DX_IMG/dims[0]);
        imgCan1.addKeyListener(this);
        imgCan1.addMouseListener(this);
        imgCan1.addMouseMotionListener(this);
        imgCan1.addMouseWheelListener(this);
        imgCan1.setSize(TARGET_WIDTH_FOR_PLOTS_AND_IMGS, TARGET_HEIGHT_FOR_PLOTS_AND_IMGS);

        //Prepare Image T2
		imgCan2.setMagnification(1.0*DX_IMG/dims[0]);
		initMag=imgCan2.getMagnification();
        imgCan2.addKeyListener(this);
        imgCan2.addMouseListener(this);
        imgCan2.addMouseMotionListener(this);
        imgCan2.addMouseWheelListener(this);
        imgCan2.setSize(TARGET_WIDTH_FOR_PLOTS_AND_IMGS, TARGET_HEIGHT_FOR_PLOTS_AND_IMGS);
 		
        //Prepare Texts
        int totalY=3*DY_TEXT+2*DELTA_Y;
        int totalX=DX_PLOT_1+DELTA_X+DX_PLOT_2;
        int eightX=totalX/8;
        int eightY=totalY/8;
        int quarterY=totalY/4;
        int []x0=new int[8]; 
        int []y0=new int[8];
        for(int i=0;i<8;i++){
        	x0[i]=DX_IMG+DELTA_X + (i * totalX)/8;
        	y0[i]= DY_PLOT+DELTA_Y + (i * totalY)/8;
       }
        String []strTitles=new String[] {" T1 mono-exponential ", " T1 mono-exponential " ," T1 bi-exponential "," T1 bionano "," T2 mono-exponential    "," T2 bi-exponential"," T2 bi-exponential"," T2 bionano"};
        String [][]strParams=new String[][] {
        	{"M0 (% cap.)"  , "T1 (ms)" , "" ,  "" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"M0 (% cap.)"  , "T1 (ms)" , "" ,  "" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"Short M0"  , "1st T1" , "Long M0" ,  "2nd T1" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"M0 (% cap.)"  , "T1" , "" ,  "" , "" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"M0 (% cap.)"  , "T2 (ms)" , "" ,  "" , "Khi2" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"Short M0"  , "Short T2" , "Long M0" ,  "Long T2" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"Short M0"  , "Short T2" , "Long M0" ,  "Long T2" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"M0 (% cap.)"  , "T2" , "" ,  "" , "" },
        	{"-"  , "-" , "" ,  "" , "-" }
        };
        for(int i=0;i<8;i++) {
        	titles[i]=new JTextField(strTitles[i],0);
        	titles[i].setEditable(false);titles[i].setFont(new Font("Helvetica", Font.BOLD, 13));titles[i].setBackground(paramsUnactivated);
    	}
        for(int i=0;i<16;i++)for(int j=0;j<5;j++) {
        	params[i][j]=new JTextField(strParams[i][j],0);
        	params[i][j].setEditable(false);params[i][j].setFont(new Font( "Helvetica" , i%2==0 ? Font.BOLD : 0, 13));params[i][j].setBackground(paramsUnactivated);
        	params[i][j].setHorizontalAlignment(JTextField.CENTER);
        }



		//Prepare plots
        if(!T1T2MixedEstimation) {
	        plotCan1.removeKeyListener(ij);
	        plotCan1.addKeyListener(this);
	        plotCan1.removeMouseListener(ij);
	        plotCan1.addMouseListener(this);
	        plotCan1.addMouseWheelListener(this);
	
	        plotCan2.removeKeyListener(ij);
	        plotCan2.addKeyListener(this);
	        plotCan2.removeMouseListener(ij);
	        plotCan2.addMouseListener(this);
	        plotCan2.addMouseWheelListener(this);
        }
        else {
	        plotCanT1T2.removeKeyListener(ij);
	        plotCanT1T2.addKeyListener(this);
	        plotCanT1T2.removeMouseListener(ij);
	        plotCanT1T2.addMouseListener(this);
	        plotCanT1T2.addMouseWheelListener(this);        	
        }
       
        plotCan21.removeKeyListener(ij);
        plotCan21.addKeyListener(this);
        plotCan21.removeMouseListener(ij);
        plotCan21.addMouseListener(this);
        plotCan21.addMouseWheelListener(this);

        plotCan22.removeKeyListener(ij);
        plotCan22.addKeyListener(this);
        plotCan22.removeMouseListener(ij);
        plotCan22.addMouseListener(this);
        plotCan22.addMouseWheelListener(this);

        
        //Image and buttons panel


		GridBagConstraints gbc;
		JTextArea text;
        
    	//Echoes
    	textInfoEchoes.setAlignmentX(CENTER_ALIGNMENT);
    	textInfoEchoes.setBackground(Color.black);
    	textInfoEchoes.setForeground(Color.white);
    	textInfoEchoes.setFont(new Font("Helvetica", Font.PLAIN, 15));
    	textInfoEchoes.setText(this.echoesText[cEchoes][zCor][tCor]);
		JPanel butEchoesPan=new JPanel();
		butEchoesPan.setLayout(new GridBagLayout());
		butEchoesPan.setBorder(BorderFactory.createEmptyBorder(0,5,0,0));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.NONE;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        butEchoesPan.add(butNextEcho,gbc);
        gbc.gridy=1;        gbc.gridheight=1;     
        butEchoesPan.add(butPrevEcho,gbc);
        if(T1T2MixedEstimation) {
    		gbc.fill = GridBagConstraints.NONE;
            gbc.gridy=2;        gbc.gridheight=1;     
            text=new JTextArea();
            text.setFont(new Font("Arial",0,20));
    		butEchoesPan.add(text,gbc);
            gbc.gridy=3;        gbc.gridheight=1;     
            text=new JTextArea();
            text.setFont(new Font("Arial",0,10));
    		butEchoesPan.add(text,gbc);
	        gbc.gridy=4;        gbc.gridheight=1;     
	        butEchoesPan.add(butNextTr,gbc);
	        gbc.gridy=5;        gbc.gridheight=1;     
	        butEchoesPan.add(butPrevTr,gbc);
            gbc.gridy=6;        gbc.gridheight=1;     
            text=new JTextArea();
            text.setFont(new Font("Arial",0,20));
    		butEchoesPan.add(text,gbc);
            gbc.gridy=7;        gbc.gridheight=1;     
	        butEchoesPan.add(butAllCurves,gbc);
        }
        
        
		JPanel butEchoesImgPan=new JPanel();
		butEchoesImgPan.setLayout(new GridBagLayout());
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        gbc.anchor=GridBagConstraints.SOUTH;
        butEchoesImgPan.add(imgCan2,gbc);
        gbc.gridy=1;        gbc.gridheight=1;     
        butEchoesImgPan.add(textInfoEchoes,gbc);
        
    	JPanel echoesPanel=new JPanel();
    	echoesPanel.setLayout(new GridBagLayout());
    	echoesPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.NONE;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        echoesPanel.add(butEchoesImgPan,gbc);
        gbc.gridx=1;        gbc.gridheight=1;     
        echoesPanel.add(butEchoesPan,gbc);

        //maps
    	textInfoMaps.setAlignmentX(CENTER_ALIGNMENT);
    	textInfoMaps.setBackground(Color.black);
    	textInfoMaps.setForeground(Color.white);
    	textInfoMaps.setFont(new Font("Helvetica", Font.PLAIN, 15));
    	textInfoMaps.setText(this.mapsText[cEchoes][zCor][tCor]);
		JPanel butMapsPan=new JPanel();
		butMapsPan.setLayout(new GridBagLayout());
		butMapsPan.setBorder(BorderFactory.createEmptyBorder(0,5,0,0));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        butMapsPan.add(butM0Map,gbc);
        gbc.gridy=1;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,10));
		gbc.fill = GridBagConstraints.NONE;
		butMapsPan.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridy=2;        gbc.gridheight=1;     
        butMapsPan.add(butT1Map,gbc);
        gbc.gridy=3;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,10));
		gbc.fill = GridBagConstraints.NONE;
		butMapsPan.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;
		gbc.gridy=4;        gbc.gridheight=1;     
        butMapsPan.add(butT2Map,gbc);
        gbc.gridy=5;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,50));
		gbc.fill = GridBagConstraints.NONE;
		butMapsPan.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;
		gbc.gridy=6;        gbc.gridheight=1;     
        butMapsPan.add(butPlay,gbc);
      
		JPanel butMapsImgPan=new JPanel();
		butMapsImgPan.setLayout(new GridBagLayout());
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        butMapsImgPan.add(textInfoMaps,gbc);
        gbc.gridy=1;        gbc.gridheight=1;     
        gbc.anchor=GridBagConstraints.NORTH;
        butMapsImgPan.add(imgCan1,gbc);
        
    	JPanel mapsPanel=new JPanel();
    	mapsPanel.setLayout(new GridBagLayout());
    	mapsPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        mapsPanel.add(butMapsImgPan,gbc);
        gbc.gridx=1;        gbc.gridheight=1;     
        mapsPanel.add(butMapsPan,gbc);

        

        
    	JPanel leftButPan=new JPanel();
    	leftButPan.setLayout(new GridBagLayout());
    	leftButPan.setBorder(BorderFactory.createEmptyBorder(0,10,0,25));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;


		
		int decY=0;
		gbc.gridwidth=2;     gbc.gridx=0;
        gbc.gridy=0+decY;        gbc.gridheight=1;     
        leftButPan.add(textTim,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=1+decY;        gbc.gridheight=1;     
        leftButPan.add(butPrevDay,gbc);
        gbc.gridwidth=1;     gbc.gridx=1;
        leftButPan.add(butNextDay,gbc);
        gbc.gridwidth=2;     gbc.gridx=0;
        gbc.gridy=2+decY;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,10));
		gbc.fill = GridBagConstraints.NONE;
        leftButPan.add(text,gbc);

		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=2;     gbc.gridx=0;gbc.gridwidth=2;
        gbc.gridy=3+decY;        gbc.gridheight=1;     
        leftButPan.add(textSli,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=4+decY;        gbc.gridheight=1;     
        leftButPan.add(butPrevSli,gbc);
        gbc.gridwidth=1;     gbc.gridx=1;
        leftButPan.add(butNextSli,gbc);
        gbc.gridwidth=2;     gbc.gridx=0;
        gbc.gridy=5+decY;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,10));
		gbc.fill = GridBagConstraints.NONE;
        leftButPan.add(text,gbc);

		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=2;     gbc.gridx=0;gbc.gridwidth=2;
        gbc.gridy=6+decY;        gbc.gridheight=1;     
        leftButPan.add(textLuts,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=7+decY;        gbc.gridheight=1;     
        leftButPan.add(butFireLut,gbc);
        gbc.gridwidth=1;     gbc.gridx=1;
        leftButPan.add(butGrayLut,gbc);


    	JPanel rightButPan=new JPanel();
    	rightButPan.setBorder(BorderFactory.createEmptyBorder(0,25,0,10));
       	rightButPan.setLayout(new GridBagLayout());
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=2;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        rightButPan.add(textSam,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=1;        gbc.gridheight=1;     
        rightButPan.add(butPrevRoi,gbc);
        gbc.gridwidth=1;     gbc.gridx=1;
        rightButPan.add(butNextRoi,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=2;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,5));
		gbc.fill = GridBagConstraints.NONE;
		rightButPan.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;

        gbc.gridwidth=2;     gbc.gridx=0;
        gbc.gridy=3;        gbc.gridheight=1;     
        rightButPan.add(textRan,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=4;        gbc.gridheight=1;     
        rightButPan.add(butPrevRang,gbc);
        gbc.gridwidth=1;     gbc.gridx=1;
        rightButPan.add(butNextRang,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=5;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,5));
		gbc.fill = GridBagConstraints.NONE;
		rightButPan.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;

        gbc.gridwidth=2;     gbc.gridx=0;
        gbc.gridy=6;        gbc.gridheight=1;     
        rightButPan.add(textMisc,gbc);
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=7;        gbc.gridheight=1;     
        rightButPan.add(butHide,gbc);
        gbc.gridwidth=1;     gbc.gridx=1;
        gbc.gridy=7;        gbc.gridheight=1;     
        rightButPan.add(butRoi,gbc);

    	
    	
    	//Combine previous in the imagesPanel
        JPanel centralButtonsPanel=new JPanel();
        centralButtonsPanel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
        centralButtonsPanel.setLayout(new GridBagLayout());
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        centralButtonsPanel.add(leftButPan,gbc);
        gbc.gridwidth=1;     gbc.gridx=1;
        centralButtonsPanel.add(new JSeparator(SwingConstants.VERTICAL),gbc);
        gbc.gridwidth=1;     gbc.gridx=2;
        centralButtonsPanel.add(rightButPan,gbc);
 
        JPanel imagesPanel=new JPanel();
		imagesPanel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		imagesPanel.setLayout(new GridBagLayout());
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        imagesPanel.add(echoesPanel,gbc);
        gbc.gridy=1;        gbc.gridheight=1;     
        imagesPanel.add(new JSeparator(),gbc);
        gbc.gridy=2;        gbc.gridheight=1;     
        imagesPanel.add(centralButtonsPanel,gbc);
        gbc.gridy=3;        gbc.gridheight=1;     
        imagesPanel.add(new JSeparator(),gbc);
        gbc.gridy=4;        gbc.gridheight=1;     
        imagesPanel.add(mapsPanel,gbc);

        

		
		
		
		
		
        //Params t1 and t2 panel
		JPanel []paramsTXPanel=new JPanel[2];
		for(int ti=0;ti<2;ti++) {
			paramsTXPanel[ti]=new JPanel();
			paramsTXPanel[ti].setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
			paramsTXPanel[ti].setLayout(new GridBagLayout());
			JPanel paramsLeftTX=new JPanel();
			paramsLeftTX.setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
			paramsLeftTX.setLayout(new GridLayout(4,1,0,0));
			JPanel paramsRightTX=new JPanel();
			paramsRightTX.setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
			paramsRightTX.setLayout(new GridLayout(4,1,0,0));
	
			paramsPanelRightTX=new JPanel [2][4];
			for(int tt=0;tt<4;tt++) {
				paramsLeftTX.add(titles[tt+ti*4]);
				titles[tt+ti*4].removeKeyListener(ij);
				titles[tt+ti*4].addMouseListener(this);
				paramsPanelRightTX[ti][tt]=new JPanel();
				paramsPanelRightTX[ti][tt].setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
				paramsPanelRightTX[ti][tt].setLayout(new GridLayout(2,5,0,0));
				for(int col=0;col<2;col++)for(int row=0;row<5;row++) {				
					paramsPanelRightTX[ti][tt].add(params[col+2*tt+8*ti][row]);
					params[col+2*tt+8*ti][row].removeKeyListener(ij);
					params[col+2*tt+8*ti][row].addMouseListener(this);
				}
				paramsRightTX.add(paramsPanelRightTX[ti][tt]);
			}
		
			
			gbc = new GridBagConstraints();
	        gbc.fill = GridBagConstraints.BOTH;
	        gbc.gridwidth=1;     gbc.gridx=0;
	
	        gbc.gridy=0;        gbc.gridheight=1;     
	        paramsTXPanel[ti].add(paramsLeftTX,gbc);
	        gbc.gridwidth=3;     gbc.gridx=1;
			paramsTXPanel[ti].add(paramsRightTX,gbc);
		}


		
		

		
		
		
		
        //T1 and T2 plots panel
		JPanel t1t2PlotsPanel=new JPanel();
		t1t2PlotsPanel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		t1t2PlotsPanel.setLayout(new GridBagLayout());
		gbc = new GridBagConstraints();
        gbc.gridy=0;        gbc.gridheight=1;     
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.fill = GridBagConstraints.BOTH;
        
		if(!T1T2MixedEstimation) {
			t1t2PlotsPanel.add(plotCan1,gbc);
	        gbc.gridx=1;        gbc.gridheight=1;     
			t1t2PlotsPanel.add(plotCan2,gbc);
		}
		else t1t2PlotsPanel.add(plotCanT1T2,gbc);
		
        //Params T1 and T2 panel
		JPanel t1t2ParamsPanel=new JPanel();
		t1t2ParamsPanel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		t1t2ParamsPanel.setLayout(new GridBagLayout());
		gbc = new GridBagConstraints();
        gbc.gridy=0;        gbc.gridheight=1;     
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.fill = GridBagConstraints.BOTH;
        
		t1t2ParamsPanel.add(paramsTXPanel[0],gbc);
        gbc.gridx=1;        gbc.gridheight=1;     
		t1t2ParamsPanel.add(paramsTXPanel[1],gbc);
		
        //Spectrum T1 and T2 panel
		JPanel t1t2SpectrumPanel=new JPanel();
		t1t2SpectrumPanel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		t1t2SpectrumPanel.setLayout(new GridBagLayout());
		gbc = new GridBagConstraints();
        gbc.gridy=0;        gbc.gridheight=1;     
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.fill = GridBagConstraints.BOTH;
        
		t1t2SpectrumPanel.add(plotCan21,gbc);
        gbc.gridx=1;        gbc.gridheight=1;     
		t1t2SpectrumPanel.add(plotCan22,gbc);
		
		
		
		int policeSize=2;
        //Right panel
		JPanel rightPanel=new JPanel();
		rightPanel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		rightPanel.setLayout(new GridBagLayout());
		gbc = new GridBagConstraints();
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.fill = GridBagConstraints.BOTH;
        
        gbc.gridy=0;        gbc.gridheight=3;     
        rightPanel.add(t1t2PlotsPanel,gbc);
        gbc.gridy=3;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,policeSize));
		gbc.fill = GridBagConstraints.NONE;
		rightPanel.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridy=4;        gbc.gridheight=1;     
        rightPanel.add(t1t2ParamsPanel,gbc);
        gbc.gridy=5;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,policeSize));
		gbc.fill = GridBagConstraints.NONE;
		rightPanel.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridy=6;        gbc.gridheight=2;     
        rightPanel.add(t1t2SpectrumPanel,gbc);
		
		
		
	

  		JPanel globalPanel=new JPanel();
		globalPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
		globalPanel.setLayout(new GridBagLayout());
		GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.fill = GridBagConstraints.NONE;
        gbc2.gridy=0;        gbc2.gridheight=1;     
        gbc2.gridwidth=1;     gbc2.gridx=0;
		globalPanel.add(imagesPanel,gbc2);
        gbc2.fill = GridBagConstraints.BOTH;
        gbc2.gridwidth=1;     gbc2.gridx=1;
		globalPanel.add(new JSeparator(SwingConstants.VERTICAL),gbc2);
        gbc2.gridwidth=2;     gbc2.gridx=2;
		globalPanel.add(rightPanel,gbc2);
        JPanel targetPanel=globalPanel;
		add(targetPanel);        

		textSli.setBackground(centralButtonsPanel.getBackground());textSli.setAlignmentX(CENTER_ALIGNMENT);
		textSam.setBackground(centralButtonsPanel.getBackground());textSam.setAlignmentX(CENTER_ALIGNMENT);
		textTim.setBackground(centralButtonsPanel.getBackground());textTim.setAlignmentX(CENTER_ALIGNMENT);
		textRan.setBackground(centralButtonsPanel.getBackground());textRan.setAlignmentX(CENTER_ALIGNMENT);
		textMisc.setBackground(centralButtonsPanel.getBackground());textMisc.setAlignmentX(CENTER_ALIGNMENT);
		textLuts.setBackground(centralButtonsPanel.getBackground());textLuts.setAlignmentX(CENTER_ALIGNMENT);
		
		butAllCurves.addActionListener(this);
		butPlay.addActionListener(this);
		butPrevDay.addActionListener(this);
		butNextDay.addActionListener(this);
		butPrevSli.addActionListener(this);
		butNextSli.addActionListener(this);
		butPrevRang.addActionListener(this);
		butNextRang.addActionListener(this);
		butPrevRoi.addActionListener(this);
		butNextRoi.addActionListener(this);
		butHide.addActionListener(this);
		butRoi.addActionListener(this);
		butNextEcho.addActionListener(this);
		butPrevEcho.addActionListener(this);
		butNextTr.addActionListener(this);
		butPrevTr.addActionListener(this);
		butM0Map.addActionListener(this);
		butT1Map.addActionListener(this);
		butT2Map.addActionListener(this);
		butFireLut.addActionListener(this);
		butGrayLut.addActionListener(this);
		if(T1T2MixedEstimation) {
			butNextEcho.setText("<Html><br>Next<br>Te<br> <Html>");
			butPrevEcho.setText("<Html><br>Prev<br>Te<br> <Html>");
		}
        GraphicsDevice device = GraphicsEnvironment.getLocalGraphicsEnvironment().getScreenDevices()[0]; 
        this.setResizable(true);
        pack();
		setVisible(true);
        setTitle("MRI Curve explorer V2");
		repaint();
	}

	public void startPlotsAndRoi(){
		IJ.log("Starting Plots");
		if(!T1T2MixedEstimation) {
			plotT1 = new Plot("T1 curve explorer","Recovery time","Spin-echo magnitude signal");
			plotT1.changeFont(new Font("Helvetica", 0, 14));
			plotCan1=new PlotCanvas(plotT1.getImagePlus());
			plotCan1.setPlot(plotT1);
	
			plotT2 = new Plot("T2 curve explorer","Echo time","Spin-echo magnitude signal");
			plotT2.changeFont(new Font("Helvetica", 0, 14));
			plotCan2=new PlotCanvas(plotT2.getImagePlus());
			plotCan2.setPlot(plotT2);
		}
		else {
			plotT1T2 = new Plot("T1-T2 curve explorer",("Repetition time + "+visualFactorT2+" x Echo time"),"Magnitude signal");
			plotT1T2.changeFont(new Font("Helvetica", 0, 14));
			plotT1T2.setSize(TARGET_WIDTH_FOR_PLOTS_T1T2, TARGET_HEIGHT_FOR_UP_PLOTS_AND_IMGS);
			plotCanT1T2=new PlotCanvas(plotT1T2.getImagePlus());
			plotCanT1T2.setPlot(plotT1T2);
//			plotCanT1T2.setSize(TARGET_WIDTH_FOR_PLOTS_T1T2, TARGET_HEIGHT_FOR_PLOTS_AND_IMGS);
		}
		
		plotT21 = new Plot("T1 timelapse tracker","T1 distribution (M0-weighted)","Observation time",Plot.DEFAULT_FLAGS-Plot.X_GRID-Plot.Y_GRID+Plot.X_LOG_TICKS+Plot.X_LOG_NUMBERS);
		plotT22 = new Plot("T2 timelapse tracker","T2 distribution (M0-weighted)","Observation time",Plot.DEFAULT_FLAGS-Plot.X_GRID-Plot.Y_GRID+Plot.X_LOG_TICKS+Plot.X_LOG_NUMBERS);
		plotT21.changeFont(new Font("Helvetica", 0, 14));
		plotT22.changeFont(new Font("Helvetica", 0, 14));
		plotT21.setBackgroundColor(new Color(0,0,0));
		plotT22.setBackgroundColor(new Color(0,0,0));
		plotT21.setSize(TARGET_WIDTH_FOR_PLOTS_ALONE, TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS);
		plotT22.setSize(TARGET_WIDTH_FOR_PLOTS_ALONE,TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS);
		plotCan21=new PlotCanvas(plotT21.getImagePlus());
		plotCan22=new PlotCanvas(plotT22.getImagePlus());
		plotCan21.setPlot(plotT21);
		plotCan22.setPlot(plotT22);
		for(int i =0;i<maxCurves;i++){
			if(!T1T2MixedEstimation) {
				plotT1.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
				plotT2.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
			}
			else plotT1T2.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
			plotT21.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
			plotT22.addPoints(new double[]{0,1},new double[]{0,0},Plot.LINE);
		}
		if(!T1T2MixedEstimation) {
			plotT1.setLimits(0, maxT1, 0, maxPlotYT1);
			plotT2.setLimits(0, maxT2, 0, maxPlotYT2);
		}
		else plotT1T2.setLimits(0, maxT1T2OnPlot+maxT2*5, 0, maxPlotYT2);
		plotT21.setLimits(t0T1,t1T1, 0, nTimepoints);
		plotT22.setLimits(t0T2,t1T2, 0, nTimepoints);
	}

	
	
	
	public void actualizeCursor() {
		Overlay overT1;
		Overlay overT2;
		if(statusRoi!=2) {
			//If no Roi is defined, draw a square with a cross in the center of it
			int xMouseCenter=(int) Math.round(xMouse+  (   Math.floor(xD)+0.5-xD  )*zoomLevel);//screen location of square center, according to xMouse	
			int yMouseCenter=(int) Math.round(yMouse+  (   Math.floor(yD)+0.5-yD  )*zoomLevel);
			double sizeOnScreen=(0.5+crossWidth)*zoomLevel;//half radius
			PointRoi prT1=new PointRoi(xMouseCenter,yMouseCenter,sizeOfCursor()+" yellow hybrid");
			overT1=new Overlay(prT1);
			overT2=new Overlay(prT1);
			Roi sq1=new Roi(xMouseCenter-sizeOnScreen-1,yMouseCenter-sizeOnScreen-1,2*sizeOnScreen+2,2*sizeOnScreen+2);
			sq1.setStrokeColor(Color.black);
			Roi sq2=new Roi(xMouseCenter-sizeOnScreen,yMouseCenter-sizeOnScreen,2*sizeOnScreen,2*sizeOnScreen);
			sq2.setStrokeColor(Color.yellow);
			overT1.add(sq1);overT2.add(sq1);overT1.add(sq2);overT2.add(sq2);

			if(rangeDisplayedFullBlocks && (this.spectrumRangingModeT1 || this.spectrumRangingModeT2) ) {
				for(int pt=0;pt<this.rangeRoiPoints.size();pt++) {
					int dx=rangeRoiPoints.get(pt)[0]-xCor;//Relative coordinates to the cross
					int dy=rangeRoiPoints.get(pt)[1]-yCor;
					Roi r=new Roi(xMouseCenter+(dx-0.5)*zoomLevel,yMouseCenter+(dy-0.5)*zoomLevel,zoomLevel,zoomLevel);
					r.setFillColor(Color.blue);
					overT1.add(r);
					overT2.add(r);
				}
			}
		}
		else {
			double[]centroid=this.userRoi.getContourCentroid();
			xCor=(int) centroid[0];
			yCor=(int) centroid[1];
			xMouse=imgCan1.screenXD(centroid[0]);
			yMouse=imgCan1.screenYD(centroid[1]);
			PointRoi prT1=new PointRoi(xMouse,yMouse,sizeOfCursor()+" yellow hybrid");
			overT1=new Overlay(prT1);
			overT2=new Overlay(prT1);

			int x0VoxOfCan=imgCan1.offScreenX(0);
			int y0VoxOfCan=imgCan1.offScreenY(0);
			int[]xVoxRoi=this.userRoi.getPolygon().xpoints;
			int[]yVoxRoi=this.userRoi.getPolygon().ypoints;
			Polygon polyScaled=new Polygon();
			for(int i=0;i<xVoxRoi.length;i++) {
				polyScaled.addPoint((int)Math.round(zoomLevel*(xVoxRoi[i]-x0VoxOfCan)),(int)Math.round(zoomLevel*(yVoxRoi[i]-y0VoxOfCan)));
			}
			Roi polyRoiScaled=new PolygonRoi(polyScaled, Roi.POLYGON);
			polyRoiScaled.setStrokeColor(Color.yellow);
			polyRoiScaled.setPosition(0);
			overT1.add(polyRoiScaled);
			overT2.add(polyRoiScaled);

			if(rangeDisplayedFullBlocks && (this.spectrumRangingModeT1 || this.spectrumRangingModeT2) ) {
				for(int pt=0;pt<this.rangeRoiPoints.size();pt++) {
					int dx=rangeRoiPoints.get(pt)[0];//Relative coordinates to the cross
					int dy=rangeRoiPoints.get(pt)[1];
					Roi r=new Roi(zoomLevel*(dx-x0VoxOfCan),zoomLevel*(dy-y0VoxOfCan),zoomLevel,zoomLevel);
					r.setFillColor(Color.blue);
					overT1.add(r);
					overT2.add(r);
				}
			}
		}

		//If spectrum ranging is active, draw squares around the pixels which estimate Ttimes falls into the selected interval
		if(roiIsHidden==false) {
			imgCan1.setOverlay(overT1);
			imgCan2.setOverlay(overT2);
		}
	}
	
	public String sizeOfCursor() {
		return "medium";
	}
		
	public int actualizeRangingBoundaries() {
		double borderLeft=76;
		double borderRight=20;
		double borderUp=18;
		double borderDown=47;
		double lengPlot=plotCan21.getWidth()-borderRight-borderLeft;
		double latPlot=(plotCan21.getHeight()-yMouseRange-borderDown)/(plotCan21.getHeight()-borderDown-borderUp);
		double xPos=this.xMouseRange-borderLeft;
		double[]vals=spectrumRangingModeT1 ? plotT21.getLimits() : plotT22.getLimits();
		double factMul=vals[1]/vals[0];
		double rangCenter=Math.pow(factMul,xPos*1.0/lengPlot)*vals[0];
		System.out.println("Range center ="+rangCenter);
		this.rangingBoundaries=new double[] {rangCenter*(1-this.rangingFactor),rangCenter,rangCenter*(1+this.rangingFactor)};
		int newTime=Math.max(0,Math.min(this.nTimepoints-1,    (int)Math.floor( this.nTimepoints*1.0*latPlot )  ));
		return newTime;
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
		IJ.log("Ending exploration");
	}

	
	
	
	public int actualizeSelectedEcho(int plotType) {
		double borderLeft=76;
		double borderRight=20;
		double borderUp=18;
		double borderBottom=44;
		double lengPlot=(plotType==WIN_PLOT1 ? plotCan1.getWidth() : (plotType==WIN_PLOT2 ? plotCan2.getWidth() : plotCanT1T2.getWidth()))-borderRight-borderLeft;
		double latPlot=(plotType==WIN_PLOT1 ? plotCan1.getHeight() : (plotType==WIN_PLOT2 ? plotCan2.getHeight() : plotCanT1T2.getHeight()))-borderUp-borderBottom;
		double xPos=(this.xMouseRangeCurve-borderLeft)/lengPlot;
		double yPos=1-(this.yMouseRangeCurve-borderUp)/latPlot;
		double[]vals=(plotType==WIN_PLOT1 ? plotT1.getLimits() : plotType==WIN_PLOT2 ? plotT2.getLimits() : plotT1T2.getLimits());
		double xPosPlot=vals[1]*xPos;
		double yPosPlot=vals[3]*yPos;
		System.out.println("Point clicked in plot : "+xPosPlot+","+yPosPlot);
		int newTime=1;double minDist=1E20;
		if(plotType==WIN_PLOT1) {
			for(int t=0;t<this.t1Times[tCor].length;t++) {
				if(Math.abs(this.t1Times[tCor][t]-xPosPlot)<minDist) {
					minDist=Math.abs(this.t1Times[tCor][t]-xPosPlot);
					newTime=t;
				}
			}
		}
		else if(plotType==WIN_PLOT2){
			for(int t=0;t<this.t2Times[tCor].length;t++) {
				if(Math.abs(this.t2Times[tCor][t]-xPosPlot)<minDist) {
					minDist=Math.abs(this.t2Times[tCor][t]-xPosPlot);
					newTime=t+this.t1Times[tCor].length;
				}
			}
		}
		else {
			for(int t=0;t<this.t1t2Times[tCor][zCor].length;t++) {
				double dist=(t1t2Times[tCor][zCor][t]-xPosPlot)*(t1t2Times[tCor][zCor][t]-xPosPlot) + (dataTimelapseT1T2[tCor][t]-yPosPlot)*(dataTimelapseT1T2[tCor][t]-yPosPlot);
				if(dist<minDist) {
					minDist=dist;
					newTime=t;
				}
			}
		}
		return newTime;
	}		
	
	

	
	/** Read new data and update estimation and fits*/	
	public void computeResultsAgain() {
		actualizeMeanEstimations();
		actualizeDisplayedNumbers();
		actualizeExportSentence();
		Timer tim2=new Timer();
		computeEstimationsForAllPointsMultiThreaded();
		if(!T1T2MixedEstimation) {
			tim2.print("Updating fit T1mono, FitT1 bi-exp, T2mono and T2 bi-exp on "+this.nPtsCur+" points over "+this.nTimepoints+" timepoints."+
					" Total fits : "+(4*this.nPtsCur*this.nTimepoints)+"  . Total computation time : ");
		}
		else 			tim2.print("Updating fit T1T21mono, and bi-exp using "+this.nPtsCur+" voxels over "+this.nTimepoints+" timepoints."+
				" Total fits : "+(2*this.nPtsCur*this.nTimepoints)+"  . Total computation time : ");
		actualizeSpectrumCurves();
		identifyRangedData();
	}

	public void actualizeMriObservationsBasedOnData() {
		//Update for all time, for all voxels selected
		if(statusRoi!=2 || userRoi==null) {
			//if(!T1T2MixedEstimation) 
				//this.dataTimelapseFull=hyperMRIT1T2.getFullMRISignalAroundThisVoxel((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick);
			//else 
				this.dataTimelapseFull=hyperMRIT1T2.getFullMRISignalAroundThisVoxelT1T2((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick,this.sigmaSmoothingBeforeEstimation);

			this.pointsCoords=hyperMRIT1T2.getCoordinatesAroundThisVoxel((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick);

		}
		else {
			IJ.log("Actualizing data from Roi, at coordinates Z="+zCor+" , T="+tCor);
			this.pointsCoords=VitimageUtils.getRoiAsCoords(this.userRoiInitial);
			//if(!T1T2MixedEstimation) 
				//this.dataTimelapseFull=hyperMRIT1T2.getFullMRISignalInTheseCoordinates((int)xCor,(int)yCor,zCor,this.pointsCoords);
			//else
				this.dataTimelapseFull=hyperMRIT1T2.getFullMRISignalInTheseCoordinatesT1T2((int)xCor,(int)yCor,zCor,this.pointsCoords,this.sigmaSmoothingBeforeEstimation);
		}
		this.nPtsCur=this.dataTimelapseFull[0].length;

		//For each time, compute the mean over voxels
		double[]vals;
		for(int t=0;t<this.nTimepoints;t++) {			
			if(!T1T2MixedEstimation) {
				this.dataTimelapseT1[t]=new double[this.dataTimelapseFull[t][0][0].length];
				this.dataTimelapseT1Sigmas[t]=new double[this.dataTimelapseFull[t][0][1].length];
				this.dataTimelapseT2[t]=new double[this.dataTimelapseFull[t][0][1].length];
				this.dataTimelapseT2Sigmas[t]=new double[this.dataTimelapseFull[t][0][1].length];
			}
			else {
				this.dataTimelapseT1T2[t]=new double[this.dataTimelapseFull[t][0][0].length];
				this.dataTimelapseT1T2Sigmas[t]=new double[this.dataTimelapseFull[t][0][0].length];				
			}
			for(int ec=0;ec<this.dataTimelapseFull[t][0][0].length;ec++) {
				vals=new double[this.nPtsCur];
				for(int vo=0;vo<this.nPtsCur;vo++) {
					vals[vo]=this.dataTimelapseFull[t][vo][0][ec];
				}
				double []stats=VitimageUtils.statistics1D(vals);
				if(!T1T2MixedEstimation) {
					this.dataTimelapseT1[t][ec]=stats[0];
					this.dataTimelapseT1Sigmas[t][ec]=stats[1];						
				}
				else {
					this.dataTimelapseT1T2[t][ec]=stats[0];
					this.dataTimelapseT1T2Sigmas[t][ec]=stats[1];											
				}
			}
			if(!T1T2MixedEstimation) {
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
	}
			
	
	
	public void actualizeMeanEstimations() {
		Object[] obj;
		for(int tim=0;tim< this.nTimepoints;tim++) {
			if(!T1T2MixedEstimation) {
		 		//Estimer T1 Monocomp
				obj=fitAndEvaluate(t1Times[tim],this.timesT1Cute,dataTimelapseT1[tim],hyperMRIT1T2.tabSigmasT1Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T1_RECOVERY_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=0;prm<2;prm++) {
		 			paramsTimelapseT1[tim][prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT1Mono[tim]=(double[]) obj[2];
		 		tabFittenT1MonoCute[tim]=(double[]) obj[3];
		 		khi2T1Mono[tim]=(double) obj[4];
		 		pValT1Mono[tim]=(double) obj[5];
		 		jitterT1Mono[tim]=(double) obj[6];
		
		 		//Estimer T1 Bicomp
				obj=fitAndEvaluate(t1Times[tim],this.timesT1Cute,dataTimelapseT1[tim],hyperMRIT1T2.tabSigmasT1Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T1_MULTICOMP_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=2;prm<6;prm++) {
		 			paramsTimelapseT1[tim][prm]=((double[]) obj[0])[prm-2];
		 		}
		 		tabFittenT1Bicomp[tim]=(double[]) obj[2];
		 		tabFittenT1BicompCute[tim]=(double[]) obj[3];
		 		khi2T1Bicomp[tim]=(double) obj[4];
		 		pValT1Bicomp[tim]=(double) obj[5];
		 		jitterT1Bicomp[tim]=(double) obj[6];
	
		 		//Estimer T2 Monocomp
		 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseT2[tim],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T2_RELAX_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=0;prm<2;prm++) {
		 			paramsTimelapseT2[tim][prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT2Mono[tim]=(double[]) obj[2];
		 		tabFittenT2MonoCute[tim]=(double[]) obj[3];
		 		khi2T2Mono[tim]=(double) obj[4];
		 		pValT2Mono[tim]=(double) obj[5];
		 		jitterT2Mono[tim]=(double) obj[6];
			
		 		//Estimer T2 Bicomp
		 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseT2[tim],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T2_MULTICOMP_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=2;prm<6;prm++) {
		 			paramsTimelapseT2[tim][prm]=((double[]) obj[0])[prm-2];
		 		}
		 		tabFittenT2Bicomp[tim]=(double[]) obj[2];
		 		tabFittenT2BicompCute[tim]=(double[]) obj[3];
		 		khi2T2Bicomp[tim]=(double) obj[4];
		 		pValT2Bicomp[tim]=(double) obj[5];
		 		jitterT2Bicomp[tim]=(double) obj[6];
			}
			else {
		 		//Estimer T1T2 Monocomp
		 		obj=fitAndEvaluateT1T2(t1t2TrTimes[tim][zCor],t1t2TeTimes[tim][zCor],this.timesT1T2Cute[tim][zCor],dataTimelapseT1T2[tim],hyperMRIT1T2.tabSigmasT1T2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T1T2_MONO_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=0;prm<3;prm++) {
		 			paramsTimelapseT1T2[tim][prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT1T2Mono[tim]=(double[]) obj[2];
		 		tabFittenT1T2MonoCute[tim]=(double[][]) obj[3];
		 		System.out.println("H1 size"+tabFittenT1T2MonoCute[tim].length);
		 		khi2T1T2Mono[tim]=(double) obj[4];
		 		pValT1T2Mono[tim]=(double) obj[5];
		 		jitterT1T2Mono[tim]=(double) obj[6];
			
		 		//Estimer T1T2 Multicomp
		 		obj=fitAndEvaluateT1T2(t1t2TrTimes[tim][zCor],t1t2TeTimes[tim][zCor],this.timesT1T2Cute[tim][zCor],dataTimelapseT1T2[tim],hyperMRIT1T2.tabSigmasT1T2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T1T2_MULTI_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=0;prm<5;prm++) {
		 			paramsTimelapseT1T2[tim][3+prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT1T2Bicomp[tim]=(double[]) obj[2];
		 		tabFittenT1T2BicompCute[tim]=(double[][]) obj[3];
		 		khi2T1T2Bicomp[tim]=(double) obj[4];
		 		pValT1T2Bicomp[tim]=(double) obj[5];
		 		jitterT1T2Bicomp[tim]=(double) obj[6];

		 		
		 		
		 		//Estimer T1T2 MultiMulticomp
		 		obj=fitAndEvaluateT1T2(t1t2TrTimes[tim][zCor],t1t2TeTimes[tim][zCor],this.timesT1T2Cute[tim][zCor],dataTimelapseT1T2[tim],hyperMRIT1T2.tabSigmasT1T2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T1T2_MULTIMULTI_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=0;prm<6;prm++) {
		 			paramsTimelapseT1T2[tim][8+prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT1BiT2Bicomp[tim]=(double[]) obj[2];
		 		tabFittenT1BiT2BicompCute[tim]=(double[][]) obj[3];
		 		khi2T1BiT2Bicomp[tim]=(double) obj[4];
		 		pValT1BiT2Bicomp[tim]=(double) obj[5];
		 		jitterT1BiT2Bicomp[tim]=(double) obj[6];
		 		
		 		
				System.out.println("\nMono : khi="+khi2T1T2Mono[tim]+" p="+pValT1T2Mono[tim]);
				System.out.println("Bi : khi="+khi2T1T2Bicomp[tim]+" p="+pValT1T2Bicomp[tim]);
				System.out.println("Bibi : khi="+khi2T1BiT2Bicomp[tim]+" p="+pValT1BiT2Bicomp[tim]);
				

		 		//Estimer Bionano
		 		obj=fitAndEvaluateT1T2(t1t2TrTimes[tim][zCor],t1t2TeTimes[tim][zCor],this.timesT1T2Cute[tim][zCor],dataTimelapseT1T2[tim],hyperMRIT1T2.tabSigmasT1T2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,MRUtils.T1T2_BIONANO,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=0;prm<4;prm++) {
		 			paramsTimelapseT1T2[tim][14+prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT1T2Bionano[tim]=(double[]) obj[2];
		 		tabFittenT1T2BionanoCute[tim]=(double[][]) obj[3];
		 		khi2T1T2Bionano[tim]=(double) obj[4];
		 		pValT1T2Bionano[tim]=(double) obj[5];
		 		jitterT1T2Bionano[tim]=(double) obj[6];
		 		
		 		
				System.out.println("\nMono : khi="+khi2T1T2Mono[tim]+" p="+pValT1T2Mono[tim]);
				System.out.println("Bi : khi="+khi2T1T2Bicomp[tim]+" p="+pValT1T2Bicomp[tim]);
				System.out.println("Bibi : khi="+khi2T1BiT2Bicomp[tim]+" p="+pValT1BiT2Bicomp[tim]);
				System.out.println("Bionano : khi="+khi2T1T2Bionano[tim]+" p="+pValT1T2Bionano[tim]);
			}
		}		
	}

/*	public void computeEstimationsForAllPoints() {
		if(!T1T2MixedEstimation)		pointsEstimations=new double[this.nTimepoints][this.nPtsCur][2][8];
		else pointsEstimations=new double[this.nTimepoints][this.nPtsCur][1][11];
		for(int tim=0;tim< this.nTimepoints;tim++) {
			for(int vo=0;vo<this.nPtsCur;vo++) {
				if(!T1T2MixedEstimation) {
			 		//Estimer T1 Monocomp
					Object[] obj=fitAndEvaluate(t1Times[tim],this.timesT1Cute,dataTimelapseFull[tim][vo][0],hyperMRIT1T2.tabSigmasT1Seq[tim][zCor],this.fitAlgorithm,MRUtils.T1_RECOVERY_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
					pointsEstimations[tim][vo][0][0]=((double[])obj[0])[0];
					pointsEstimations[tim][vo][0][1]=((double[])obj[0])[1];
					pointsEstimations[tim][vo][0][2]=((double) obj[6])>=9999 ? 0 : 1;
					double khiT1Mono=(double) obj[4];
			
			 		//Estimer T1 Bicomp
			 		obj=fitAndEvaluate(t1Times[tim],this.timesT1Cute,dataTimelapseFull[tim][vo][0],hyperMRIT1T2.tabSigmasT1Seq[tim][zCor],this.fitAlgorithm,MRUtils.T1_MULTICOMP_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
					pointsEstimations[tim][vo][0][3]=((double[])obj[0])[0];
					pointsEstimations[tim][vo][0][4]=((double[])obj[0])[1];
					pointsEstimations[tim][vo][0][5]=((double[])obj[0])[2];
					pointsEstimations[tim][vo][0][6]=((double[])obj[0])[3];
					pointsEstimations[tim][vo][0][7]=((double) obj[6])>=9999 ? 0 : 1;
					double khiT1Bi=(double) obj[4];
	
					//Estimer T2 Monocomp
			 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseFull[tim][vo][1],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor],this.fitAlgorithm,MRUtils.T2_RELAX_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
					pointsEstimations[tim][vo][1][0]=((double[])obj[0])[0];
					pointsEstimations[tim][vo][1][1]=((double[])obj[0])[1];
					pointsEstimations[tim][vo][1][2]=((double) obj[6])>=9999 ? 0 : 1;
					double khiT2Mono=(double) obj[4];
	
			 		//Estimer T2 Bicomp
			 		obj=fitAndEvaluate(t2Times[tim],this.timesT2Cute,dataTimelapseFull[tim][vo][1],hyperMRIT1T2.tabSigmasT2Seq[tim][zCor],this.fitAlgorithm,MRUtils.T2_MULTICOMP_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
					pointsEstimations[tim][vo][1][3]=((double[])obj[0])[0];
					pointsEstimations[tim][vo][1][4]=((double[])obj[0])[1];
					pointsEstimations[tim][vo][1][5]=((double[])obj[0])[2];
					pointsEstimations[tim][vo][1][6]=((double[])obj[0])[3];
					pointsEstimations[tim][vo][1][7]=((double) obj[6])>=9999 ? 0 : 1;
					double khiT2Bi=(double) obj[4];
	
					//Select best T1
					if(khiT1Mono<=khiT1Bi && pointsEstimations[tim][vo][0][2]==1)pointsEstimations[tim][vo][0][7]=0;
					if(khiT1Mono>khiT1Bi && pointsEstimations[tim][vo][0][7]==1)pointsEstimations[tim][vo][0][2]=0;		 		
					//Select best T2
					if(khiT2Mono<=khiT2Bi && pointsEstimations[tim][vo][1][2]==1)pointsEstimations[tim][vo][1][7]=0;
					if(khiT2Mono>khiT2Bi && pointsEstimations[tim][vo][1][7]==1)pointsEstimations[tim][vo][1][2]=0;		 
				}
				else {
			 		//Estimer T1T2 Monocomp
					Object[] obj=fitAndEvaluateT1T2(t1t2TrTimes[tim],t1t2TeTimes[tim],this.timesT1T2Cute[tim],dataTimelapseT1T2[tim],hyperMRIT1T2.tabSigmasT1T2Seq[tim][zCor],this.fitAlgorithm,MRUtils.T1T2_MONO_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
					pointsEstimations[tim][vo][0][0]=((double[])obj[0])[0];
					pointsEstimations[tim][vo][0][1]=((double[])obj[0])[1];
					pointsEstimations[tim][vo][0][2]=((double[])obj[0])[1];
					pointsEstimations[tim][vo][0][3]=((double) obj[6])>=9999 ? 0 : 1;
					double khiT1T2Mono=(double) obj[4];
			
			 		//Estimer T1 Bicomp
					obj=fitAndEvaluateT1T2(t1t2TrTimes[tim],t1t2TeTimes[tim],this.timesT1T2Cute[tim],dataTimelapseT1T2[tim],hyperMRIT1T2.tabSigmasT1T2Seq[tim][zCor],this.fitAlgorithm,MRUtils.T1T2_MULTI_RICE,this.nRepetMonteCarlo,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,false);
					pointsEstimations[tim][vo][0][4]=((double[])obj[0])[0];
					pointsEstimations[tim][vo][0][5]=((double[])obj[0])[1];
					pointsEstimations[tim][vo][0][6]=((double[])obj[0])[2];
					pointsEstimations[tim][vo][0][7]=((double[])obj[0])[3];
					pointsEstimations[tim][vo][0][8]=((double[])obj[0])[4];
					pointsEstimations[tim][vo][0][9]=((double[])obj[0])[5];
					pointsEstimations[tim][vo][0][10]=((double) obj[6])>=9999 ? 0 : 1;
					double khiT1T2Bi=(double) obj[4];
	
					//Select best T1
					if(khiT1T2Mono<=khiT1T2Bi && pointsEstimations[tim][vo][0][3]==1)pointsEstimations[tim][vo][0][10]=0;
					if(khiT1T2Mono>khiT1T2Bi && pointsEstimations[tim][vo][0][10]==1)pointsEstimations[tim][vo][0][3]=0;		 		
				}
			}
		}		
	}
	*/

	public void computeEstimationsForAllPointsMultiThreaded() {
		if(T1T2MixedEstimation) {computeEstimationsForAllPointsMultiThreadedT1T2();return;}
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
					dataParams[nt][nVo][nTime][2]=new double[] {hyperMRIT1T2.tabSigmasT1Seq[nTime][zCor]+deltaSigma};
					dataParams[nt][nVo][nTime][3]=t2Times[nTime];
					dataParams[nt][nVo][nTime][4]=dataTimelapseFull[nTime][listVoxOnThreads[nt][nVo]][1];
					dataParams[nt][nVo][nTime][5]=new double[] {hyperMRIT1T2.tabSigmasT2Seq[nTime][zCor]+deltaSigma};
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
						 		//Estimer T1 Monocomp
								Object[] obj=fitAndEvaluate(dataParams[numThread][numVo][numTime][0],null,dataParams[numThread][numVo][numTime][1],dataParams[numThread][numVo][numTime][2][0],finalFitAlgorithm,MRUtils.T1_RECOVERY_RICE,finalRepetMonteCarlo,1,false);
								tempParams[0][0]=((double[])obj[0])[0];
								tempParams[0][1]=((double[])obj[0])[1];
								tempParams[0][2]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1Mono=(double) obj[4];
						
						 		//Estimer T1 Bicomp
						 		obj=fitAndEvaluate(dataParams[numThread][numVo][numTime][0],null,dataParams[numThread][numVo][numTime][1],dataParams[numThread][numVo][numTime][2][0],finalFitAlgorithm,MRUtils.T1_MULTICOMP_RICE,finalRepetMonteCarlo,1,false);
						 		tempParams[0][3]=((double[])obj[0])[0];
						 		tempParams[0][4]=((double[])obj[0])[1];
						 		tempParams[0][5]=((double[])obj[0])[2];
						 		tempParams[0][6]=((double[])obj[0])[3];
						 		tempParams[0][7]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1Bi=(double) obj[4];

								//Estimer T2 Monocomp
						 		obj=fitAndEvaluate(dataParams[numThread][numVo][numTime][3],null,dataParams[numThread][numVo][numTime][4],dataParams[numThread][numVo][numTime][5][0],finalFitAlgorithm,MRUtils.T2_RELAX_RICE,finalRepetMonteCarlo,1,false);	
						 		tempParams[1][0]=((double[])obj[0])[0];
						 		tempParams[1][1]=((double[])obj[0])[1];
						 		tempParams[1][2]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT2Mono=(double) obj[4];

						 		//Estimer T2 Bicomp
						 		obj=fitAndEvaluate(dataParams[numThread][numVo][numTime][3],null,dataParams[numThread][numVo][numTime][4],dataParams[numThread][numVo][numTime][5][0],finalFitAlgorithm,MRUtils.T2_MULTICOMP_RICE,finalRepetMonteCarlo,1,false);
						 		tempParams[1][3]=((double[])obj[0])[0];
						 		tempParams[1][4]=((double[])obj[0])[1];
						 		tempParams[1][5]=((double[])obj[0])[2];
						 		tempParams[1][6]=((double[])obj[0])[3];
						 		tempParams[1][7]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT2Bi=(double) obj[4];

								//Select best T1
								if(khiT1Mono<=khiT1Bi && tempParams[0][2]==1)tempParams[0][7]=0;
								if(khiT1Mono>khiT1Bi && tempParams[0][7]==1)tempParams[0][2]=0;		 		
								//Select best T2
								if(khiT2Mono<=khiT2Bi && tempParams[1][2]==1)tempParams[1][7]=0;
								if(khiT2Mono>khiT2Bi && tempParams[1][7]==1)tempParams[1][2]=0;		 		
								estimatedParams[numThread][numVo][numTime]=tempParams;
							}							
						}
					} catch(Exception ie) {ie.printStackTrace();}
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

	
	
	
	public void computeEstimationsForAllPointsMultiThreadedT1T2() {
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
					dataParams[nt][nVo][nTime]=new double[4][];
					dataParams[nt][nVo][nTime][0]=t1t2TrTimes[nTime][zCor];
					dataParams[nt][nVo][nTime][1]=t1t2TeTimes[nTime][zCor];//TODO : danger at the slice 0 and 1 if crossThick > 0, estimation will be abusive. Cross thick have been prevented acting, thus
					dataParams[nt][nVo][nTime][2]=dataTimelapseFull[nTime][listVoxOnThreads[nt][nVo]][0];
					dataParams[nt][nVo][nTime][3]=new double[] {hyperMRIT1T2.tabSigmasT1T2Seq[nTime][zCor]};
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
								double[][]tempParams=new double[2][17];
						 		//Estimer T1T2 Monocomp
								Object[] obj=fitAndEvaluateT1T2(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,MRUtils.T1T2_MONO_RICE,finalRepetMonteCarlo,1,false);
								tempParams[0][0]=((double[])obj[0])[0];
								tempParams[0][1]=((double[])obj[0])[1];
								tempParams[0][2]=((double[])obj[0])[2];
								tempParams[0][3]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1T2Mono=(double) obj[4];
								double pT1T2Mono=(double) obj[5];
						
						 		//Estimer T1MonoT2Bi
								obj=fitAndEvaluateT1T2(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,MRUtils.T1T2_MULTI_RICE,finalRepetMonteCarlo,1,false);
								tempParams[0][4]=((double[])obj[0])[0];
								tempParams[0][5]=((double[])obj[0])[1];
								tempParams[0][6]=((double[])obj[0])[2];
								tempParams[0][7]=((double[])obj[0])[3];
								tempParams[0][8]=((double[])obj[0])[4];
								tempParams[0][9]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1T2Bi=(double) obj[4];
								double pT1T2Bi=(double) obj[5];

								//Estimer T1BiT2Bi
								obj=fitAndEvaluateT1T2(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,MRUtils.T1T2_MULTIMULTI_RICE,finalRepetMonteCarlo,1,false);
								tempParams[0][10]=((double[])obj[0])[0];
								tempParams[0][11]=((double[])obj[0])[1];
								tempParams[0][12]=((double[])obj[0])[2];
								tempParams[0][13]=((double[])obj[0])[3];
								tempParams[0][14]=((double[])obj[0])[4];
								tempParams[0][15]=((double[])obj[0])[5];
								tempParams[0][16]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1BiT2Bi=(double) obj[4];
								double pT1BiT2Bi=(double) obj[5];

								
								//Select best
								if(pT1T2Mono<=pT1T2Bi && pT1T2Mono<=pT1BiT2Bi)tempParams[0][9]=tempParams[0][16]=0;
								if(pT1T2Bi<pT1T2Mono && pT1T2Bi<=pT1BiT2Bi)tempParams[0][3]=tempParams[0][16]=0;
								if(pT1BiT2Bi<pT1T2Mono && pT1BiT2Bi<pT1T2Bi)tempParams[0][3]=tempParams[0][9]=0;
								estimatedParams[numThread][numVo][numTime]=tempParams;
							}							
						}
					} catch(Exception ie) {ie.printStackTrace();}
				} 
			};  		
		}				
		VitimageUtils.startAndJoin(threads);

		
		//Gather results
		pointsEstimations=new double[this.nTimepoints][this.nPtsCur][1][17];
		for(int nt=0;nt<nThreads;nt++) {
			for(int nVo=0;nVo<listVoxOnThreads[nt].length;nVo++) {
				for(int nTime=0;nTime<this.nTimepoints;nTime++) {
					pointsEstimations[nTime][listVoxOnThreads[nt][nVo]]=estimatedParams[nt][nVo][nTime];
				}
			}
		}
	}

	
	
	
	
	
	
	
	
	
	
	
	public Object[] fitAndEvaluate(double[]tabTimes,double[]tabTimesCute,double[]tabData,double sigmaRice,int fitAlgorithm,int fitCurveType,int nbMonteCarloSimulations,int nbPts,boolean niceCurveComputation) {
		int nParams=(fitCurveType==MRUtils.T1_RECOVERY_RICE || fitCurveType==MRUtils.T2_RELAX_RICE) ? 2 : 4; 
		boolean isT1=(fitCurveType==MRUtils.T1_RECOVERY_RICE || fitCurveType==MRUtils.T1_MULTICOMP_RICE); 
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
		if((!isT1) && tabData[0]<(statsRice[0]+3*statsRice[1]) ) {//T2 curve with the first point being very low -> no estimation possible there
			jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 computation cause first point is too low : "+tabData[0]+" , "+statsRice[0]+" , "+statsRice[1]);
		}
		if((isT1) && tabData[tabData.length-1]<(statsRice[0]+3*statsRice[1]) ) {//T1 curve with the last point being very low -> no estimation possible there
			jitter=9999;accs[1]=100;//System.out.println("Jitter set in T1 computation cause last point is too low : "+tabData[0]+" , "+statsRice[0]+" , "+statsRice[1]);
		}
		if(isT1 && nParams==2) {//T1_MONOCOMP
			if(estimatedParams[0]<0 || estimatedParams[0]>maxAcceptableM0  
			|| estimatedParams[1]<minAcceptableT1 || estimatedParams[1]>maxAcceptableT1)  {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T1 monocomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}
		}
		if(isT1 && nParams==4) {//T1_BICOMP
			if(estimatedParams[0]<0 || estimatedParams[2]<0 ||  (estimatedParams[0]+estimatedParams[2])>maxAcceptableM0  
		    || estimatedParams[1]<minAcceptableT1 || estimatedParams[1]>maxAcceptableT1   || estimatedParams[3]<minAcceptableT1 || estimatedParams[3]>maxAcceptableT1  )  {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T1 bicomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}
		}
		if(!isT1 && nParams==2){//T2 MONOCOMP
			if(estimatedParams[0]<0 || estimatedParams[0]>maxAcceptableM0 || estimatedParams[1]<minAcceptableT2 || estimatedParams[1]>maxAcceptableT2 ) {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 monocomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}
		}
		if(!isT1 && nParams==4){//T2 BICOMP
			if(estimatedParams[0]<0 || estimatedParams[2]<0 || (estimatedParams[0]+estimatedParams[2])>maxAcceptableM0 
			|| estimatedParams[1]<minAcceptableT2 || estimatedParams[1]>maxAcceptableT2   || estimatedParams[3]<minAcceptableT2 || estimatedParams[3]>maxAcceptableT2) {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 bicomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}		
		}
		return new Object[] {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accs[0],accs[1],jitter};
	}
	

	
	
	
	public Object[] fitAndEvaluateT1T2(double[]tabTimesTr,double[]tabTimesTe,double[][][]tabCuteTimesPerCurve,double[]tabData,double sigmaRice,int fitAlgorithm,int fitCurveType,int nbMonteCarloSimulations,int nbPts,boolean niceCurveComputation) {
		int type=(fitCurveType==MRUtils.T1T2_MONO_RICE ? 0 : (fitCurveType==MRUtils.T1T2_MULTI_RICE ? 1 : 2 ) );
		int nParams= (type==0) ? 3 : (type==1 ? 5 : 6);
		if(fitCurveType==MRUtils.T1T2_BIONANO) {nParams=4;type=4;};
		double []estimatedParams=new double[nParams];
		double []estimatedSigmas=new double[nParams];
		double []estimation=MRUtils.makeFitDouble(tabTimesTr, tabTimesTe, tabData, fitCurveType, MRUtils.SIMPLEX, MRUtils.N_ITER_T1T2, sigmaRice);
		for(int i=0;i<nParams;i++) {estimatedParams[i]=estimation[i] ;  estimatedSigmas[i]=estimation[i];}
		double []tabFitten=MRUtils.fittenRelaxationCurveT1T2(tabTimesTr,tabTimesTe,estimatedParams,sigmaRice,fitCurveType);
		double [][]tabFittenCute=null;
		if(niceCurveComputation) {
			tabFittenCute=new double[tabCuteTimesPerCurve.length][];
			for(int fi=0;fi<tabCuteTimesPerCurve.length;fi++) {
				tabFittenCute[fi]=MRUtils.fittenRelaxationCurveT1T2(tabCuteTimesPerCurve[fi][0],tabCuteTimesPerCurve[fi][1],estimatedParams,sigmaRice,fitCurveType);
			}
		}
		double[]accs=MRUtils.fittingAccuraciesT1T2(tabData,tabTimesTr,tabTimesTe,sigmaRice,estimatedParams,fitCurveType,false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
		if((nParams==6) && (estimatedParams[4]>estimatedParams[5])) {
			double tmp=estimatedParams[0];estimatedParams[0]=estimatedParams[1];estimatedParams[1]=tmp;
			tmp=estimatedParams[2];estimatedParams[2]=estimatedParams[3];estimatedParams[3]=tmp;
			tmp=estimatedSigmas[4];estimatedSigmas[4]=estimatedSigmas[5];estimatedSigmas[5]=tmp;
		}		
		
		if((nParams==5) && (estimatedParams[3]>estimatedParams[4])) {
			double tmp=estimatedParams[0];estimatedParams[0]=estimatedParams[1];estimatedParams[1]=tmp;
			tmp=estimatedParams[3];estimatedParams[3]=estimatedParams[4];estimatedParams[4]=tmp;
		}		
		
		double jitter=0;
		if(VitimageUtils.max(tabData)<MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*sigmaRice)  {//T2 curve with the first point being very low -> no estimation possible there
			jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 computation cause first point is too low : "+tabData[0]+" , "+statsRice[0]+" , "+statsRice[1]);
		}
		if(type==0) {
			if(estimatedParams[0]<=0 || estimatedParams[0]>maxAcceptableM0  
			|| estimatedParams[1]<minAcceptableT1 || estimatedParams[1]>maxAcceptableT1
			|| estimatedParams[2]<minAcceptableT2 || estimatedParams[2]>maxAcceptableT2)  {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T1 monocomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}
		}
		else if(type==1) {
			if(estimatedParams[0]<=0 || estimatedParams[1]<=0 || (estimatedParams[0]+estimatedParams[1])>maxAcceptableM0 
			|| estimatedParams[2]<minAcceptableT1 || estimatedParams[2]>maxAcceptableT1
			|| estimatedParams[3]<minAcceptableT2 || estimatedParams[3]>maxAcceptableT2   || estimatedParams[4]<minAcceptableT2 || estimatedParams[4]>maxAcceptableT2) {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 bicomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}		
		}
		else if(type==2){
			if(estimatedParams[0]<=0 || estimatedParams[1]<=0 || (estimatedParams[0]+estimatedParams[1])>maxAcceptableM0 
			|| estimatedParams[2]<minAcceptableT1 || estimatedParams[2]>maxAcceptableT1   || estimatedParams[3]<minAcceptableT1 || estimatedParams[3]>maxAcceptableT1
			|| estimatedParams[4]<minAcceptableT2 || estimatedParams[4]>maxAcceptableT2   || estimatedParams[5]<minAcceptableT2 || estimatedParams[5]>maxAcceptableT2) {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 bicomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}		
		}
		else {
			if(estimatedParams[0]<0 || estimatedParams[2]<0 || estimatedParams[1]<=0 || estimatedParams[0]>maxAcceptableM0 || estimatedParams[3]<=0 || estimatedParams[2]>maxAcceptableM0 
			|| estimatedParams[1]<minAcceptableT1 || estimatedParams[1]>maxAcceptableT1 
			|| estimatedParams[3]<minAcceptableT2 || estimatedParams[3]>maxAcceptableT2 ) {
				jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 bicomp. Parameters were "+TransformUtils.stringVectorN(estimatedParams, ""));
			}		
			
		}
		return new Object[] {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accs[0],accs[1],jitter};
	}
	

	
	
	
	
	
	
	public void identifyRangedData() {
		if(T1T2MixedEstimation) {
			identifyRangedDataT1T2();return;
		}
		this.rangeRoiPoints=new ArrayList<int[]>();
		if(this.spectrumRangingModeT1) {
			//Ranging T1 values
			for(int p=0;p<pointsCoords.length;p++) {
				if ( (statusRoi==2 ||  pointsCoords[p][2]==zCor) && /*Point is on current slice*/
						 this.pointsEstimations[tCor][p][0][2]==1 /*Point estimation is valid in monomodal (jitter not too big)*/ && 
						 (this.pointsEstimations[tCor][p][0][1]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][1]>this.rangingBoundaries[0]) /*Estimated T1 is in the range*/
					 ) {
					rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
				if ( (statusRoi==2 || pointsCoords[p][2]==zCor) && /*Point is on current slice*/
						 this.pointsEstimations[tCor][p][0][7]==1 /*Point estimation is valid in bimodal (jitter not too big)*/ && 
						 (this.pointsEstimations[tCor][p][0][4]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][4]>this.rangingBoundaries[0]) || /*Estimated short T1 is in the range*/
						 (this.pointsEstimations[tCor][p][0][6]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][6]>this.rangingBoundaries[0]) /*Estimated long T1 is in the range*/
					 ) {
					rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
			}
			System.out.println("Ranging T1 : ["+this.rangingBoundaries[0]+"-"+this.rangingBoundaries[2]+"] . After update : "+rangeRoiPoints.size()+" points");
		}
		if(this.spectrumRangingModeT2) {
			//Ranging T2 values
			for(int p=0;p<pointsCoords.length;p++) {
				if ( (statusRoi==2 ||  pointsCoords[p][2]==zCor) && /*Point is on current slice*/
					 (this.pointsEstimations[tCor][p][1][2]==1) /*Point estimation is valid in T2 monomodal (jitter not too big)*/ && 
					 (this.pointsEstimations[tCor][p][1][1]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][1][1]>this.rangingBoundaries[0]) /*Estimated T2 is in the range*/
				 ) {
					 rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
				if ( (statusRoi==2 || pointsCoords[p][2]==zCor) && /*Point is on current slice*/
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

	public void identifyRangedDataT1T2() {
		this.rangeRoiPoints=new ArrayList<int[]>();
		if(this.spectrumRangingModeT1) {
			//Ranging T1 values
			for(int p=0;p<pointsCoords.length;p++) {
				if ( (statusRoi==2 ||  pointsCoords[p][2]==zCor) && /*Point is on current slice*/
						 this.pointsEstimations[tCor][p][0][3]==1 /*Point estimation is valid in monomodal (jitter not too big)*/ && 
						 (this.pointsEstimations[tCor][p][0][1]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][1]>this.rangingBoundaries[0]) /*Estimated T1 is in the range*/
					 ) {
					rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
				if ( (statusRoi==2 || pointsCoords[p][2]==zCor) && /*Point is on current slice*/
						 this.pointsEstimations[tCor][p][0][9]==1 /*Point estimation is valid in bimodal (jitter not too big)*/ && 
						 (this.pointsEstimations[tCor][p][0][6]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][6]>this.rangingBoundaries[0])/*Estimated short T1 is in the range*/
					 ) {
					rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
				if ( (statusRoi==2 || pointsCoords[p][2]==zCor) && /*Point is on current slice*/
						 this.pointsEstimations[tCor][p][0][16]==1 /*Point estimation is valid in bimodal (jitter not too big)*/ && 
						 (this.pointsEstimations[tCor][p][0][12]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][12]>this.rangingBoundaries[0]) || /*Estimated short T1 is in the range*/
						 (this.pointsEstimations[tCor][p][0][13]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][13]>this.rangingBoundaries[0]) /*Estimated short T1 is in the range*/
					 ) {
					rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
			}
			System.out.println("Ranging T1 : ["+this.rangingBoundaries[0]+"-"+this.rangingBoundaries[2]+"] . After update : "+rangeRoiPoints.size()+" points");
		}
		if(this.spectrumRangingModeT2) {
			//Ranging T2 values
			for(int p=0;p<pointsCoords.length;p++) {
				if ( (statusRoi==2 ||  pointsCoords[p][2]==zCor) && /*Point is on current slice*/
					 (this.pointsEstimations[tCor][p][0][3]==1) /*Point estimation is valid in T2 monomodal (jitter not too big)*/ && 
					 (this.pointsEstimations[tCor][p][0][2]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][1][2]>this.rangingBoundaries[0]) /*Estimated T2 is in the range*/
				 ) {
					 rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
				if ( (statusRoi==2 || pointsCoords[p][2]==zCor) && /*Point is on current slice*/
					 (this.pointsEstimations[tCor][p][0][9]==1) /*Point estimation is valid in T2 monomodal (jitter not too big)*/ && 
					 (   (this.pointsEstimations[tCor][p][0][8]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][8]>this.rangingBoundaries[0]) || /*Estimated short T2 is in the range*/
					     (this.pointsEstimations[tCor][p][0][7]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][7]>this.rangingBoundaries[0]) ) /*Estimated long T2 is in the range*/
					 ) {
					 rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
				}
				if ( (statusRoi==2 || pointsCoords[p][2]==zCor) && /*Point is on current slice*/
						 (this.pointsEstimations[tCor][p][0][16]==1) /*Point estimation is valid in T2 monomodal (jitter not too big)*/ && 
						 (   (this.pointsEstimations[tCor][p][0][14]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][14]>this.rangingBoundaries[0]) || /*Estimated short T2 is in the range*/
						     (this.pointsEstimations[tCor][p][0][15]<this.rangingBoundaries[2]) && (this.pointsEstimations[tCor][p][0][15]>this.rangingBoundaries[0]) ) /*Estimated long T2 is in the range*/
						 ) {
						 rangeRoiPoints.add(new int[] {this.pointsCoords[p][0],this.pointsCoords[p][1]});
					}
			}
			System.out.println("Ranging T2 : ["+this.rangingBoundaries[0]+"-"+this.rangingBoundaries[2]+"] . After update : "+rangeRoiPoints.size()+" points");
		}
	}
	
	public double[][] computeSpectrumCurve(int time) {
		if(T1T2MixedEstimation) {return computeSpectrumCurveT1T2(time);}
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
			if(pointsEstimations[time][p][0][2]==1) {//Fit T1 mono is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][1]<binSup[bin]) && (pointsEstimations[time][p][0][1]>=binInf[bin]) ) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][0];
				}
			}
			if(pointsEstimations[time][p][0][7]>=1) {//Fit T1 bicomp is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][4]<binSup[bin]) && (pointsEstimations[time][p][0][4]>=binInf[bin]) ) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][3];
				}
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][6]<binSup[bin]) && (pointsEstimations[time][p][0][6]>=binInf[bin]) ) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][5];
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
		
	public double[][] computeSpectrumCurveT1T2(int time) {
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
			if(pointsEstimations[time][p][0][3]==1) {//Fit T1T2 mono is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][1]<binSup[bin]) && (pointsEstimations[time][p][0][1]>=binInf[bin]) ) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][0];
				}
			}
			if(pointsEstimations[time][p][0][9]>=1) {//Fit T1T2 bicomp is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][6]<binSup[bin]) && (pointsEstimations[time][p][0][6]>=binInf[bin]) ) {
					histoM0T1[bin]+=(pointsEstimations[time][p][0][4]+pointsEstimations[time][p][0][5]);
				}
			}
			if(pointsEstimations[time][p][0][16]>=1) {//Fit T1BiT2 bicomp is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][12]<binSup[bin]) && (pointsEstimations[time][p][0][12]>=binInf[bin]) ) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][10];
				}
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][13]<binSup[bin]) && (pointsEstimations[time][p][0][13]>=binInf[bin]) ) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][11];
				}
			}

			
			if(pointsEstimations[time][p][0][3]>=1) {//Fit T1T2 mono is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][2]<binSup[bin]) && (pointsEstimations[time][p][0][2]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][0];
				}
			}
			if(pointsEstimations[time][p][0][9]>=1) {//Fit T1T2 bicomp is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][7]<binSup[bin]) && (pointsEstimations[time][p][0][7]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][4];
				}
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][8]<binSup[bin]) && (pointsEstimations[time][p][0][8]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][5];
				}
			}
			if(pointsEstimations[time][p][0][16]>=1) {//Fit T1T2 bicomp is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][14]<binSup[bin]) && (pointsEstimations[time][p][0][14]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][10];
				}
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][15]<binSup[bin]) && (pointsEstimations[time][p][0][15]>=binInf[bin]) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][11];
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
		if(T1T2MixedEstimation) {actualizeFirstPlotsT1T2();return;}
		int incrT1=0;
		int incrT2=0;
		int deltaXt2=3;
		int bubbleSize=30;
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
	
		
		
		
		///////////////// T1 FITTING CURVES
        //NOISE
        plotT1.setLineWidth(1);
		double[]statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMRIT1T2.tabSigmasT1Seq[tCor][zCor]+deltaSigma);
		double nSum=statusRoi==2 ? Math.sqrt(this.nPtsCur) : (1+crossWidth);
		this.meanNoiseT1Cur=statsNoise[0];
		this.sigmaNoiseT1Cur=statsNoise[1];
		plotT1.setColor(Color.black);
        plotT1.setLineWidth(2);
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1},new double[]{statsNoise[0],statsNoise[0]});//Afficher le mean+sigma
        plotT1.setLineWidth(1);
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1},new double[]{statsNoise[0]+statsNoise[1]/nSum,statsNoise[0]+statsNoise[1]/nSum});//Afficher le mean+sigma
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1},new double[]{statsNoise[0]-statsNoise[1]/nSum,statsNoise[0]-statsNoise[1]/nSum});//Afficher le mean-sigma

		

        

		//OBSERVATIONS IRM
        plotT1.setLineWidth(2);
		plotT1.setColor(crossColor );
		plotT1.replace(incrT1++,"x",t1Times[tCor],dataTimelapseT1[tCor]);//Afficher les points IRM
        plotT1.setLineWidth(1);
		for(int t=0;t<this.t1Times[tCor].length;t++)plotT1.replace(incrT1++, "line",new double[]{this.t1Times[tCor][t]-deltaXt2,this.t1Times[tCor][t]-deltaXt2},new double[]{dataTimelapseT1[tCor][t]-dataTimelapseT1Sigmas[tCor][t],dataTimelapseT1[tCor][t]+dataTimelapseT1Sigmas[tCor][t]});//Afficher le T2

        //FIT T1 Mono
        plotT1.setColor((jitterT1Mono[tCor]>=9999) ? Color.gray : curveT1Mono);        
        plotT1.setLineWidth(1+thickCurve);
        plotT1.replace(incrT1++, "line",timesT1Cute,tabFittenT1MonoCute[tCor]);//Afficher la courbe 
        plotT1.setLineWidth(3);
        //FIT T1 Bicomp
        plotT1.setColor((jitterT1Bicomp[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
        plotT1.setLineWidth(1+thickCurve);
        plotT1.replace(incrT1++, "line",timesT1Cute,tabFittenT1BicompCute[tCor]);//Afficher la courbe 
        plotT1.setLineWidth(3);
        //Quite a hack, in order that if both are on the same trajectories, that obviously the two components will be gray, that is not handsome to be displayed
        plotT1.setColor((jitterT1Mono[tCor]>=9999) ? Color.gray : curveT1Mono);        
        plotT1.setLineWidth(1+thickCurve);
        plotT1.replace(incrT1++, "line",timesT1Cute,tabFittenT1MonoCute[tCor]);//Afficher la courbe 
        plotT1.setLineWidth(3);
       

        //PARAMS ESTIMATED
        plotT1.setLineWidth(3);
        plotT1.setColor((jitterT1Mono[tCor]>=9999) ? Color.gray : curveT1Mono);        
		plotT1.replace(incrT1++, "line",new double[]{0,(double)(paramsTimelapseT1[tCor][1])},new double[]{maxPlotYT1*0.85,maxPlotYT1*0.85});//Afficher le T1
        plotT1.setColor((jitterT1Bicomp[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
		plotT1.replace(incrT1++, "line",new double[]{0,(double)(paramsTimelapseT1[tCor][3])},new double[]{maxPlotYT1*0.90,maxPlotYT1*0.90});//Afficher le T1 court
		plotT1.replace(incrT1++, "line",new double[]{0,(double)(paramsTimelapseT1[tCor][5])},new double[]{maxPlotYT1*0.80,maxPlotYT1*0.80});//Afficher le T1

        
        //LEGENDE T1 AND GLOBAL PLOT DESIGN
        plotT1.setLineWidth(1);
		String strLegendT1="";
		for(int hh=0;hh<1;hh++) strLegendT1+="\n";
		strLegendT1+="Noise +-sigma"+"\n"+"\n"+"\n"+"Spin-echo signal"+"\n-Fit one exp.\n-Fit two exp.\nbla\nbla";
		plotT1.setLimits(0, maxT1, 0,maxPlotYT1);
		plotT1.setColor(new Color(150,150 ,150) );
		plotT1.addLegend(strLegendT1,"bottom-right");
		plotCan1.setPlot(plotT1);	

		
		
		
		
	
		///////////////// T2 FITTING CURVES
		//NOISE LEVEL
		statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMRIT1T2.tabSigmasT2Seq[tCor][zCor]+deltaSigma);
		this.meanNoiseT2Cur=statsNoise[0];
		this.sigmaNoiseT2Cur=statsNoise[1];
		nSum=Math.sqrt(this.nPtsCur);
		plotT2.setColor(Color.black);
		plotT2.setLineWidth(2);
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2},new double[]{statsNoise[0],statsNoise[0]});//Afficher le mean+sigma
		plotT2.setLineWidth(1);
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2},new double[]{statsNoise[0]+statsNoise[1]/nSum,statsNoise[0]+statsNoise[1]/nSum});//Afficher le mean+sigma
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2},new double[]{statsNoise[0]-statsNoise[1]/nSum,statsNoise[0]-statsNoise[1]/nSum});//Afficher le mean-sigma

        //DONNEES IRM
		plotT2.setLineWidth(2);
		plotT2.setColor(crossColor);
		plotT2.replace(incrT2++,"x",t2Times[tCor],dataTimelapseT2[tCor]);
        plotT2.setLineWidth(1);
		for(int t=0;t<this.t2Times[tCor].length;t++)plotT2.replace(incrT2++, "line",new double[]{this.t2Times[tCor][t],this.t2Times[tCor][t]},new double[]{dataTimelapseT2[tCor][t]-dataTimelapseT2Sigmas[tCor][t],dataTimelapseT2[tCor][t]+dataTimelapseT2Sigmas[tCor][t]});
      
        //COURBES FITTED
        plotT2.setColor((jitterT2Mono[tCor]>=9999) ? Color.gray : curveT2Mono);
        plotT2.setLineWidth(1+thickCurve);
        plotT2.replace(incrT2++, "line",timesT2Cute,tabFittenT2MonoCute[tCor]);//Afficher la courbe monocomp
        plotT2.setColor((jitterT2Bicomp[tCor]>=9999) ? Color.gray : curveT2Bicomp);
        plotT2.setLineWidth(1+thickCurve);
        plotT2.replace(incrT2++, "line",timesT2Cute,tabFittenT2BicompCute[tCor]);//Afficher la courbe bicomp
        plotT2.setColor((jitterT2Mono[tCor]>=9999) ? Color.gray : curveT2Mono);
        plotT2.setLineWidth(1+thickCurve);
        plotT2.replace(incrT2++, "line",timesT2Cute,tabFittenT2MonoCute[tCor]);//Afficher la courbe monocomp
      
        //PARAM T2 MONO ET STD MONO
        plotT2.setLineWidth(3);
        plotT2.setColor((jitterT2Mono[tCor]>=9999) ? Color.gray : curveT2Mono);
		plotT2.replace(incrT2++, "line",new double[]{0,(double)(paramsTimelapseT2[tCor][1])},new double[]{maxPlotYT2*0.85,maxPlotYT2*0.85});//Afficher le T2
        plotT2.setColor((jitterT2Bicomp[tCor]>=9999) ? Color.gray : curveT2Bicomp);
		plotT2.replace(incrT2++, "line",new double[]{0,paramsTimelapseT2[tCor][3]},new double[]{maxPlotYT2*0.90,maxPlotYT2*0.90});//Afficher le T2
		plotT2.replace(incrT2++, "line",new double[]{0,paramsTimelapseT2[tCor][5]},new double[]{maxPlotYT2*0.80,maxPlotYT2*0.80});//Afficher le T2

        //LEGENDE
		plotT2.setLineWidth(1);
		String strLegendT2="\n\nNoise +-sigma\nMRI T2 relaxation"+"\nFit one exp.\nFit two exp.";
		plotT2.setLimits(0, maxT2, 0,maxPlotYT2);
		plotT2.setColor(new Color(150,150 ,150) );
		plotT2.addLegend(strLegendT2,"bottom-right");
		plotCan2.setPlot(plotT2);

		plotT1.setLineWidth(1);
		for(int cur=incrT1;cur<lastCount1;cur++)plotT1.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		plotT2.setLineWidth(1);
		for(int cur=incrT2;cur<lastCount2;cur++)plotT2.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		lastCount1=incrT1;
		lastCount2=incrT2;

	
	}
		
	
/*	public double graphicTime(int tC,int cEc) {
		return this.t1t2TrTimes[tC][zCor][cEc] + visualFactorT2 * this.t1t2TeTimes[tC][cEc];
	}
	*/
	
	/** Update graphs using new computed results*/		
	public void actualizeFirstPlotsT1T2() {
		int incrT1T2=0;
		int deltaXt2=3;
		int bubbleSize=30;
		if(this.autoSizingTimePlots) maxPlotYT1T2=VitimageUtils.max(dataTimelapseT1T2[tCor])*1.6;

		double maxCurve=0;
		if(switchT1T2==3) maxCurve=VitimageUtils.max(tabFittenT1T2BionanoCute[tCor][0]);
		if(switchT1T2==2) maxCurve=VitimageUtils.max(tabFittenT1BiT2BicompCute[tCor][0]);
		if(switchT1T2==1) maxCurve=VitimageUtils.max(tabFittenT1T2BicompCute[tCor][0]);
		if(switchT1T2==0) maxCurve=VitimageUtils.max(tabFittenT1T2MonoCute[tCor][0]);
		if(maxCurve>0.98*maxPlotYT1T2)maxPlotYT1T2=maxCurve*1.2;
		
		//BUBBLE SHOWING CURRENT ECHO
		plotT1T2.setColor(bubbleColor);
        plotT1T2.setLineWidth(bubbleSize);
        plotT1T2.replace(incrT1T2++, "line",new double[]{t1t2Times[tCor][zCor][cEchoes],t1t2Times[tCor][zCor][cEchoes]},new double[]{dataTimelapseT1T2[tCor][cEchoes],dataTimelapseT1T2[tCor][cEchoes]});
		
		
        //NOISE
        plotT1T2.setLineWidth(1);
		double[]statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMRIT1T2.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma);
		double nSum=statusRoi==2 ? Math.sqrt(this.nPtsCur) : (1+crossWidth);
		this.meanNoiseT1T2Cur=statsNoise[0];
		this.sigmaNoiseT1T2Cur=statsNoise[1];
		plotT1T2.setColor(Color.black);
        plotT1T2.setLineWidth(2);
		plotT1T2.replace(incrT1T2++, "line",new double[]{0,maxT1T2OnPlot},new double[]{statsNoise[0],statsNoise[0]});//Afficher le mean+sigma
        plotT1T2.setLineWidth(1);
		plotT1T2.replace(incrT1T2++, "line",new double[]{0,maxT1T2OnPlot},new double[]{statsNoise[0]+statsNoise[1]/nSum,statsNoise[0]+statsNoise[1]/nSum});//Afficher le mean+sigma
		plotT1T2.replace(incrT1T2++, "line",new double[]{0,maxT1T2OnPlot},new double[]{statsNoise[0]-statsNoise[1]/nSum,statsNoise[0]-statsNoise[1]/nSum});//Afficher le mean-sigma


		//OBSERVATIONS IRM MEAN AND STDEV
        plotT1T2.setLineWidth(2);
		plotT1T2.setColor(crossColor );
		plotT1T2.replace(incrT1T2++,"x",t1t2selectCute(t1t2Times[tCor][zCor]),t1t2selectCute(dataTimelapseT1T2[tCor]));//Afficher les points IRM
        plotT1T2.setLineWidth(1);
		for(int t=0;t<t1t2selectCute(this.t1t2Times[tCor][zCor]).length;t++)plotT1T2.replace(incrT1T2++, "line",new double[]{t1t2selectCute(this.t1t2Times[tCor][zCor])[t]-deltaXt2,t1t2selectCute(this.t1t2Times[tCor][zCor])[t]-deltaXt2},new double[]{t1t2selectCute(dataTimelapseT1T2[tCor])[t]-t1t2selectCute(dataTimelapseT1T2Sigmas[tCor])[t],t1t2selectCute(dataTimelapseT1T2[tCor])[t]+t1t2selectCute(dataTimelapseT1T2Sigmas[tCor])[t]});//Afficher le T2

		
		if(switchT1T2==3) {
			//FIT T1 Bionano
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1T2Bionano[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][0][2],tabFittenT1T2BionanoCute[tCor][0]);//Afficher la courbe 
	           	
	        //FIT T2 Bionano
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1T2Bionano[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=1;tran<timesT1T2Cute[tCor][zCor].length;tran++) {
	        	 if(t1t2SelectedCurve(tran))plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][tran][2],tabFittenT1T2BionanoCute[tCor][tran]);//Afficher la courbe Bicomp
	        }
		}
		else if(switchT1T2==2) {
			//FIT T1 Bicomp
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1BiT2Bicomp[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][0][2],tabFittenT1BiT2BicompCute[tCor][0]);//Afficher la courbe 
	           	
	        //FIT T2 Bicomp
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1BiT2Bicomp[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=1;tran<timesT1T2Cute[tCor][zCor].length;tran++) {
	        	 if(t1t2SelectedCurve(tran))plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][tran][2],tabFittenT1BiT2BicompCute[tCor][tran]);//Afficher la courbe Bicomp
	        }
		}
		else if(switchT1T2==0) {
	        //FIT T1 Mono
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1T2Mono[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][0][2],tabFittenT1T2MonoCute[tCor][0]);//Afficher la courbe 
	
	        //FIT T2 Mono
	        plotT1T2.setColor((jitterT1T2Mono[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=1;tran<timesT1T2Cute[tCor][zCor].length;tran++) {
		        if(t1t2SelectedCurve(tran))plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][tran][2],tabFittenT1T2MonoCute[tCor][tran]);//Afficher la courbe Bicomp
	        }
		}        
        

		else  {
	        //FIT T1 Mono
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1T2Bicomp[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][0][2],tabFittenT1T2BicompCute[tCor][0]);//Afficher la courbe 
	
	        //FIT T2 Mono
	        plotT1T2.setColor((jitterT1T2Bicomp[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=1;tran<timesT1T2Cute[tCor][zCor].length;tran++) {
		        if(t1t2SelectedCurve(tran))plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][tran][2],tabFittenT1T2BicompCute[tCor][tran]);//Afficher la courbe Bicomp
	        }
		}        

        //LEGENDE
		plotT1T2.setLineWidth(1);
		String strLegendT2="\n\nNoise +-sigma\n\nSpin-echo signal"+"\nData std over Roi\n";
		
		plotT1T2.setLimits(0, maxT1T2OnPlot, 0,maxPlotYT1T2);
		plotT1T2.setColor(new Color(150,150 ,150) );
		plotT1T2.addLegend(strLegendT2,"bottom-up");
		plotCanT1T2.setPlot(plotT1T2);

		plotT1T2.setLineWidth(1);
		for(int cur=incrT1T2;cur<lastCountT1T2;cur++)plotT1T2.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		plotT1T2.setLineWidth(1);
		for(int cur=incrT1T2;cur<lastCount2;cur++)plotT1T2.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		lastCountT1T2=incrT1T2;
	}
		
	
	
	
	
	public void actualizeSecondPlots() {
		int incrT1=0;
		int incrT2=0;
		
        //HORIZONTAL GRID T21 ET T22
   		for(int tim=0;tim<this.nTimepoints;tim++) {
			plotT21.setLineWidth(4);
			plotT21.setColor(Color.darkGray);
			plotT21.replace(incrT1++, "line",new double[] {t0T1,t1T1},new double[] {tim,tim});
			plotT22.setLineWidth(4);
			plotT22.setColor(Color.darkGray);
			plotT22.replace(incrT2++, "line",new double[] {t0T2,t1T2},new double[] {tim,tim});
   		}

        //VERTICAL GRID T21
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT21.setLineWidth(1);        
			plotT21.setColor(Color.lightGray);
        	for(int mul=2;mul<=9;mul++)plotT21.replace(incrT1++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.nTimepoints});
			plotT21.setColor(Color.white);
        	plotT21.replace(incrT1++, "line",new double[] {tt*3,tt*3},new double[] {0,this.nTimepoints});
	        plotT21.setLineWidth(3);        
        	plotT21.replace(incrT1++, "line",new double[] {tt,tt},new double[] {0,this.nTimepoints});
        }
        //VERTICAL GRID T22
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT22.setLineWidth(1);        
			plotT22.setColor(Color.lightGray);
			for(int mul=2;mul<=9;mul++)plotT22.replace(incrT2++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.nTimepoints});
			plotT22.setColor(Color.white);
           	plotT22.replace(incrT2++, "line",new double[] {tt*3,tt*3},new double[] {0,this.nTimepoints});
	        plotT22.setLineWidth(3);        
        	plotT22.replace(incrT2++, "line",new double[] {tt,tt},new double[] {0,this.nTimepoints});
        }

         //Draw the suns		
		double x0T1=t0T1*1.15;			double x0T2=t0T2*1.15;
		plotT21.setLineWidth(30);
		plotT21.setColor(bubbleColor);
        plotT21.replace(incrT1++, "line",new double[] {x0T1,x0T1},new double[] {tCor+0.2,tCor+0.2});
		plotT22.setLineWidth(30);
		plotT22.setColor(bubbleColor);
        plotT22.replace(incrT2++, "line",new double[] {x0T2,x0T2},new double[] {tCor+0.2,tCor+0.2});
   
        
        
   		for(int tim=0;tim<this.nTimepoints;tim++) {
	 		double radius=0;
	 		plotT21.setLineWidth(6);
	 		plotT22.setLineWidth(6);

	 		//Draw markers T1
			plotT21.setColor((jitterT1Bicomp[tim]>=9999) ? Color.lightGray : curveT1Bicomp);
	 		radius=0.05+0.5*paramsTimelapseT1[tim][2]/MRUtils.maxM0ForNormalization;
	 		radius=Math.min(0.9, Math.max(0.05, radius));
			plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim][3],paramsTimelapseT1[tim][3]},new double[] {tim+0.02,tim+0.02+radius});//Afficher la courbe
	 		radius=0.05+0.5*paramsTimelapseT1[tim][4]/MRUtils.maxM0ForNormalization;
	 		radius=Math.min(0.9, Math.max(0.05, radius));
			plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim][5],paramsTimelapseT1[tim][5]},new double[] {tim+0.02,tim+0.02+radius});//Afficher la courbe

			plotT21.setColor((jitterT1Mono[tim]>=9999) ? Color.gray : curveT1Mono);
	 		radius=0.05+0.5*paramsTimelapseT1[tim][0]/MRUtils.maxM0ForNormalization;
	 		radius=Math.min(0.9, Math.max(0.05, radius));
			plotT21.replace(incrT1++, "line",new double[] {paramsTimelapseT1[tim][1],paramsTimelapseT1[tim][1]},new double[] {tim+0.02,tim+0.02+radius});//Afficher la courbe

			
	 		//Draw markers T2
			plotT22.setColor((jitterT2Mono[tim]>=9999) ? Color.gray : curveT2Mono);
	 		radius=0.05+0.5*paramsTimelapseT2[tim][0]/MRUtils.maxM0ForNormalization;
	 		radius=Math.min(0.9, Math.max(0.05, radius));
			plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][1],paramsTimelapseT2[tim][1]},new double[] {tim+0.02,tim+0.02+radius});//Afficher la courbe

			plotT22.setColor((jitterT2Bicomp[tim]>=9999) ? Color.lightGray : curveT2Bicomp);
	 		radius=0.05+0.5*paramsTimelapseT2[tim][2]/MRUtils.maxM0ForNormalization;
	 		radius=Math.min(0.9, Math.max(0.05, radius));
			plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][3],paramsTimelapseT2[tim][3]},new double[] {tim+0.02,tim+0.02+radius});//Afficher la courbe
	 		radius=0.05+0.5*paramsTimelapseT2[tim][4]/MRUtils.maxM0ForNormalization;
	 		radius=Math.min(0.9, Math.max(0.05, radius));
			plotT22.replace(incrT2++, "line",new double[] {paramsTimelapseT2[tim][5],paramsTimelapseT2[tim][5]},new double[] {tim+0.02,tim+0.02+radius});//Afficher la courbe
			
			
	 		//Draw spectrum curves
	 		plotT21.setLineWidth(2);
	 		plotT21.setColor(Color.red);
	 		plotT21.replace(incrT1++, "line", valsSpectrum[tim][2], valsSpectrum[tim][0]);
	 		plotT22.setLineWidth(2);
	 		plotT22.setColor(Color.green);
	 		plotT22.replace(incrT2++, "line", valsSpectrum[tim][2], valsSpectrum[tim][1]);
 
	 		
	 		//Draw ranging interval
	 		if(rangeDisplayedFullBlocks && this.spectrumRangingModeT1 && tim==tCor) {
	 			plotT21.setColor(Color.blue);	 		
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[0]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[2],this.rangingBoundaries[2]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[1],this.rangingBoundaries[1]}, new double[] {tim+1-0.6,tim+1-0.4});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[2]}, new double[] {tim+1-0.5,tim+1-0.5});
	 		}   			
	 		if(rangeDisplayedFullBlocks && this.spectrumRangingModeT2 && tim==tCor) {
	 			plotT22.setColor(Color.blue);	 		
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[0]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[2],this.rangingBoundaries[2]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[1],this.rangingBoundaries[1]}, new double[] {tim+1-0.6,tim+1-0.4});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[2]}, new double[] {tim+1-0.5,tim+1-0.5});
	 		}   		 		
   		}				
		plotT21.setLimits(t0T1,t1T1, 0, nTimepoints);
		plotT22.setLimits(t0T2,t1T2, 0, nTimepoints);

		plotT21.setLineWidth(1);
		for(int cur=incrT1;cur<lastCount21;cur++)plotT21.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		plotT22.setLineWidth(1);
		for(int cur=incrT2;cur<lastCount22;cur++)plotT22.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		lastCount21=incrT1;
		lastCount22=incrT2;
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
		if(T1T2MixedEstimation)return;
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
		rSentence+=""+""+(paramsTimelapseT2[nTimepoints-1][5])+"), JitterT1Mono=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(jitterT1Mono[t])+",";
		rSentence+=""+""+(jitterT1Mono[nTimepoints-1])+"), JitterT1BiExp=c(";
		for(int t=0;t<nTimepoints-1;t++)rSentence+=""+(jitterT1Bicomp[t])+",";
		rSentence+=""+""+(jitterT1Bicomp[nTimepoints-1])+"), JitterT2MonoExp=c(";
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
			pySentence+=","+jitterT1Mono[t]+","+jitterT1Bicomp[t]+","+jitterT2Mono[t]+","+jitterT2Bicomp[t];
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
			matSentence+=","+jitterT1Mono[t]+","+jitterT1Mono[t]+","+jitterT2Mono[t]+","+jitterT2Bicomp[t];
			matSentence+="]"+(t==nTimepoints-1 ? "\n" : ";\n");
		}
		matSentence+="]\n";
	}
	
	public void actualizeDisplayedNumbers() {
		//Update colors
		boolean t1MonoActive=T1T2MixedEstimation ? (switchT1T2==0) : (jitterT1Mono[tCor]<9999);
		boolean t2MonoActive=T1T2MixedEstimation ? (switchT1T2==0) : (jitterT2Mono[tCor]<9999);
		boolean t1BiActive=T1T2MixedEstimation ? (switchT1T2==1) : (jitterT1Bicomp[tCor]<9999);
		boolean t2BiActive=T1T2MixedEstimation ? (switchT1T2==1) : (jitterT2Bicomp[tCor]<9999);
		boolean t1BiBiActive=T1T2MixedEstimation ? (switchT1T2==2) : (jitterT1Bicomp[tCor]<9999);
		boolean t2BiBiActive=T1T2MixedEstimation ? (switchT1T2==2) : (jitterT2Bicomp[tCor]<9999);
		boolean t1BionanoActive=T1T2MixedEstimation ? (switchT1T2==3) : (jitterT1T2Bionano[tCor]<9999);
		boolean t2BionanoActive=T1T2MixedEstimation ? (switchT1T2==3) : (jitterT1T2Bionano[tCor]<9999);
		if(!T1T2MixedEstimation){
			if(t1BiActive && t1MonoActive) {
				if(khi2T1Mono[tCor]>khi2T1Bicomp[tCor])t1MonoActive=false;
				else t1BiActive=false;
			}
			if(t2BiActive && t2MonoActive) {
				if(khi2T2Mono[tCor]>khi2T2Bicomp[tCor])t2MonoActive=false;
				else t2BiActive=false;
			}
		}
		
		
		titles[0].setBackground(t1MonoActive ? paramsT1ActivatedColor : paramsUnactivated);
		titles[1].setBackground(t1BiActive ? paramsT1ActivatedColor : paramsUnactivated);
		titles[2].setBackground(t1BiBiActive ? paramsT1ActivatedColor : paramsUnactivated);
		titles[3].setBackground(t1BionanoActive ? paramsT1ActivatedColor : paramsUnactivated);
		titles[4].setBackground(t2MonoActive ? paramsT2ActivatedColor : paramsUnactivated);
		titles[5].setBackground(t2BiActive ? paramsT2ActivatedColor : paramsUnactivated);
		titles[6].setBackground(t2BiBiActive ? paramsT2ActivatedColor : paramsUnactivated);
		titles[7].setBackground(t2BionanoActive ? paramsT2ActivatedColor : paramsUnactivated);

		for(int j=0;j<5;j++) {
			params[0][j].setBackground(t1MonoActive ? paramsT1ActivatedColor : paramsUnactivated);
			params[1][j].setBackground(t1MonoActive ? paramsT1TextActivatedColor : paramsUnactivated);
			params[2][j].setBackground(t1BiActive ? paramsT1ActivatedColor : paramsUnactivated);
			params[3][j].setBackground(t1BiActive ? paramsT1TextActivatedColor : paramsUnactivated);
			params[4][j].setBackground(t1BiBiActive ? paramsT1ActivatedColor : paramsUnactivated);
			params[5][j].setBackground(t1BiBiActive ? paramsT1TextActivatedColor : paramsUnactivated);
			params[6][j].setBackground(t1BionanoActive ? paramsT1ActivatedColor : paramsUnactivated);
			params[7][j].setBackground(t2BionanoActive ? paramsT1TextActivatedColor : paramsUnactivated);

			params[8][j].setBackground(t2MonoActive ? paramsT2ActivatedColor : paramsUnactivated);
			params[9][j].setBackground(t2MonoActive ? paramsT2TextActivatedColor : paramsUnactivated);
			params[10][j].setBackground(t2BiActive ? paramsT2ActivatedColor : paramsUnactivated);
			params[11][j].setBackground(t2BiActive ? paramsT2TextActivatedColor : paramsUnactivated);
			params[12][j].setBackground(t2BiBiActive ? paramsT2ActivatedColor : paramsUnactivated);
			params[13][j].setBackground(t2BiBiActive ? paramsT2TextActivatedColor : paramsUnactivated);
			params[14][j].setBackground(t1BionanoActive ? paramsT2ActivatedColor : paramsUnactivated);
			params[15][j].setBackground(t2BionanoActive ? paramsT2TextActivatedColor : paramsUnactivated);
		}
		//Update values

		params[1][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][0] : paramsTimelapseT1[this.tCor][0])));
		params[1][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][1] : paramsTimelapseT1[this.tCor][1]));
		params[1][4].setText(String.format("%5.4f",T1T2MixedEstimation ? khi2T1T2Mono[this.tCor] : this.khi2T1Mono[this.tCor]));

		params[3][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? (paramsTimelapseT1T2[this.tCor][3]+paramsTimelapseT1T2[this.tCor][4]) : paramsTimelapseT1[this.tCor][0])));
		params[3][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][5] : paramsTimelapseT1[this.tCor][1]));
		params[3][4].setText(String.format("%5.4f",T1T2MixedEstimation ? khi2T1T2Bicomp[this.tCor] : this.khi2T1Mono[this.tCor]));

		params[5][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][8] : paramsTimelapseT1[this.tCor][2])));
		params[5][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][10] : paramsTimelapseT1[this.tCor][3]));
		params[5][2].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][9] : paramsTimelapseT1[this.tCor][4])));
		params[5][3].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][11] : paramsTimelapseT1[this.tCor][5]));
		params[5][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1BiT2Bicomp[this.tCor] : this.khi2T1Bicomp[this.tCor]));

		params[7][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][14] : paramsTimelapseT1[this.tCor][2])));
		params[7][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][15] : paramsTimelapseT1[this.tCor][3]));
		params[7][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1T2Bionano[this.tCor] : this.khi2T1T2Bionano[this.tCor]));

		params[9][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][0] : paramsTimelapseT2[this.tCor][0])));
		params[9][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][2] : paramsTimelapseT2[this.tCor][1]));
		params[9][4].setText(String.format("%5.4f",T1T2MixedEstimation ? khi2T1T2Mono[this.tCor] : this.khi2T2Mono[this.tCor]));

		params[11][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][3] : paramsTimelapseT2[this.tCor][2])));
		params[11][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][6] : paramsTimelapseT2[this.tCor][3]));
		params[11][2].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][4] : paramsTimelapseT2[this.tCor][4])));
		params[11][3].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][7] : paramsTimelapseT2[this.tCor][5]));
		params[11][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1T2Bicomp[this.tCor] : this.khi2T2Bicomp[this.tCor]));

		params[13][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][8] : paramsTimelapseT2[this.tCor][2])));
		params[13][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][12] : paramsTimelapseT2[this.tCor][3]));
		params[13][2].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][9] : paramsTimelapseT2[this.tCor][4])));
		params[13][3].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][13] : paramsTimelapseT2[this.tCor][5]));
		params[13][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1BiT2Bicomp[this.tCor] : this.khi2T2Bicomp[this.tCor]));
		
		params[15][0].setText(String.format("%5.3f",MRUtils.convertM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][16] : paramsTimelapseT1[this.tCor][2])));
		params[15][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][17] : paramsTimelapseT1[this.tCor][3]));
		params[15][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1T2Bionano[this.tCor] : this.khi2T1T2Bionano[this.tCor]));

		//Update texts
		textInfoEchoes.setText(this.echoesText[cEchoes][zCor][tCor]);
    	textInfoMaps.setText(this.mapsText[cMaps][zCor][tCor]);
    	textSam.setText("       Sample size ("+(1+2*crossWidth)+"x"+(1+2*crossWidth)+"x"+(1+2*crossThick)+")");

	}

	public void displayResultsAgain() {
		imgCan1.getImage().setPosition(cMaps+1,zCor+1,tCor+1);
		imgCan1.setImageUpdated();
		if(isFireLut) IJ.run(imgCan1.getImage(),"Fire","");
		else IJ.run(imgCan1.getImage(),"Grays","");
		imgCan2.getImage().setPosition(cEchoes+1,zCor+1,tCor+1);
		imgCan2.setImageUpdated();
		if(isFireLut) IJ.run(imgCan2.getImage(),"Fire","");
		else IJ.run(imgCan2.getImage(),"Grays","");
		actualizeFirstPlots();
		actualizeSecondPlots();
		actualizeDisplayedNumbers();
		actualizeCursor();
		repaintAll();
	}
	
	
	
		
		
	
	
	

	
	
	

	/** Gui updating functions and callbacks*/	
	@Override
	public void keyTyped(KeyEvent e) {
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
//			if(statusRoi==2) {}//VitiDialogs.notYet("Warning : no zoom in / out during Roi mode");return;
			imgCan1.zoomIn(xMouse,yMouse);
			imgCan2.zoomIn(xMouse,yMouse);
			zoomLevel=imgCan1.getMagnification();
			xMouse=imgCan1.screenX(xCor);
			yMouse=imgCan1.screenY(yCor);
			xD=imgCan1.offScreenXD(xMouse);
			yD=imgCan1.offScreenYD(yMouse);
			displayResultsAgain();
		}
		if (e.getKeyChar()=='-') {
//			if(statusRoi==2) {}//VitiDialogs.notYet("Warning : no zoom in / out during Roi mode");return;
			if(imgCan1.getMagnification()>(initMag*1.3333)) {
				imgCan1.setMagnification(imgCan1.getMagnification()/2.0);
				imgCan2.setMagnification(imgCan2.getMagnification()/2.0);
				xMouse=imgCan1.screenX(xCor);
				yMouse=imgCan1.screenY(yCor);
				imgCan1.zoomIn(xMouse,yMouse);
				imgCan2.zoomIn(xMouse,yMouse);
				xMouse=imgCan1.screenX(xCor);
				yMouse=imgCan1.screenY(yCor);
				zoomLevel=imgCan1.getMagnification();
			}
			else {
				imgCan1.setMagnification(initMag/1.333);
				imgCan2.setMagnification(initMag/1.333);
				xMouse=imgCan1.screenX(xCor);
				yMouse=imgCan1.screenY(yCor);
				imgCan1.zoomIn(xMouse,yMouse);
				imgCan2.zoomIn(xMouse,yMouse);
				xMouse=imgCan1.screenX(xCor);
				yMouse=imgCan1.screenY(yCor);
				zoomLevel=imgCan1.getMagnification();
			}
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
			if((e.getKeyChar()=='3') && (cMaps<2) )cMaps+=1;
			if((e.getKeyChar()=='1') && (cMaps>0) )cMaps-=1;

			//Move c of echos
			if((e.getKeyChar()=='9') && (cEchoes<hyperMRIT1T2.nChannels-1) )cEchoes+=1;
			if((e.getKeyChar()=='7') && (cEchoes>0) )cEchoes-=1;

			if(e.getKeyChar()=='4' || e.getKeyChar()=='6') {
				actualizeMriObservationsBasedOnData();
				computeResultsAgain();
			}
			else {
				identifyRangedData();
			}
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
		System.out.println("Key pressed !");
		repaintAll();
	}
	

	public void changeTr(boolean up) {
		double globalIndex=cEchoes;
		int correspondingFirstEchoGlobalIndex=cEchoes;
		int correspondingLastEchoGlobalIndex=cEchoes;
		int correspondingFirstEchoGlobalIndexPrevious=cEchoes;
		int correspondingLastEchoGlobalIndexPrevious=cEchoes;
		int correspondingFirstEchoGlobalIndexNext=cEchoes;
		int correspondingLastEchoGlobalIndexNext=cEchoes;
		double curTr=this.t1t2TrTimes[tCor][zCor][cEchoes];
		int nEchoes=this.t1t2TrTimes[tCor][zCor].length;
		
		//Reaching first echo of the current Tr
		correspondingFirstEchoGlobalIndex=cEchoes;
		double prevTr=curTr;
		int indexReach=cEchoes;		
		while(prevTr==curTr && indexReach>0 ) {
			indexReach--;
			prevTr=this.t1t2TrTimes[tCor][zCor][indexReach];
			if(prevTr==curTr)correspondingFirstEchoGlobalIndex=indexReach;
		}

		//Reaching last echo of the current Tr
		correspondingLastEchoGlobalIndex=cEchoes;
		indexReach=cEchoes;
		prevTr=curTr;
		while(prevTr==curTr && indexReach<nEchoes-1 ) {
			indexReach++;
			prevTr=this.t1t2TrTimes[tCor][zCor][indexReach];
			if(prevTr==curTr)correspondingLastEchoGlobalIndex=indexReach;
		}
		System.out.println("DEB : found curTr from "+correspondingFirstEchoGlobalIndex+" to "+correspondingLastEchoGlobalIndex);
		

		//Reaching last echo of the previous Tr
		if(correspondingFirstEchoGlobalIndex==0) {
			correspondingFirstEchoGlobalIndexPrevious=0;
			correspondingLastEchoGlobalIndexPrevious=correspondingLastEchoGlobalIndex;
		}
		else {
			//Reaching last echo of the current Tr
			correspondingLastEchoGlobalIndexPrevious=correspondingFirstEchoGlobalIndex-1;
			indexReach=correspondingLastEchoGlobalIndexPrevious;
			correspondingFirstEchoGlobalIndexPrevious=indexReach;
			curTr=this.t1t2TrTimes[tCor][zCor][indexReach];
			prevTr=curTr;
			while(prevTr==curTr && indexReach>0 ) {
				System.out.println("While with indexReach="+indexReach);
				indexReach--;
				prevTr=this.t1t2TrTimes[tCor][zCor][indexReach];
				if(prevTr==curTr)correspondingFirstEchoGlobalIndexPrevious=indexReach;
			}
		}
		System.out.println("DEB : found prevTr from "+correspondingFirstEchoGlobalIndexPrevious+" to "+correspondingLastEchoGlobalIndexPrevious);
		
		//Reaching first echo of the next Tr
		if(correspondingLastEchoGlobalIndex==(nEchoes-1)) {
			correspondingFirstEchoGlobalIndexNext=correspondingFirstEchoGlobalIndex;
			correspondingLastEchoGlobalIndexNext=correspondingLastEchoGlobalIndex;
		}
		else {
			//Reaching first echo of the current Tr
			correspondingFirstEchoGlobalIndexNext=correspondingLastEchoGlobalIndex+1;
			correspondingLastEchoGlobalIndexNext=correspondingFirstEchoGlobalIndexNext;
			indexReach=correspondingFirstEchoGlobalIndexNext;
			curTr=this.t1t2TrTimes[tCor][zCor][indexReach];
			prevTr=curTr;
			while(prevTr==curTr && indexReach<nEchoes-1 ) {
				indexReach++;
				prevTr=this.t1t2TrTimes[tCor][zCor][indexReach];
				if(prevTr==curTr)correspondingLastEchoGlobalIndexNext=indexReach;
			}
		}
		System.out.println("DEB : found nextTr from "+correspondingFirstEchoGlobalIndexNext+" to "+correspondingLastEchoGlobalIndexNext);
		
		//Game
		int indexRelativeToTr=cEchoes-correspondingFirstEchoGlobalIndex;
		int retIndex=0;
		if(up) {
			retIndex=indexRelativeToTr+correspondingFirstEchoGlobalIndexNext;
			if (retIndex>correspondingLastEchoGlobalIndexNext)retIndex=correspondingLastEchoGlobalIndexNext;
		}
		else {
			retIndex=indexRelativeToTr+correspondingFirstEchoGlobalIndexPrevious;
			if (retIndex>correspondingLastEchoGlobalIndexPrevious)retIndex=correspondingLastEchoGlobalIndexPrevious;
		}
		cEchoes=retIndex;return;//Game
	}
	



	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource()==butAllCurves) {
			switchAllCurves=!switchAllCurves;
			this.displayResultsAgain();
		}
		if(e.getSource()==butPlay) {
			if(! playMode) {
				this.memSwitch=this.switchT1T2;
				this.switchT1T2=0;
				this.computeResultsAgain();
				this.displayResultsAgain();
				playMode=true;
				for(int i=0;i<0;i++)updatePlaying(0,0);
				for(int i=0;i<0;i++)updatePlaying(1,1);
				for(int i=0;i<0;i++)updatePlaying(2,0);
				System.out.println("PlayMode : "+playMode);
				return;
			}
			else {
				playMode=false;
				System.out.println("PlayMode : "+playMode);
				this.switchT1T2=memSwitch;
				this.computeResultsAgain();
				this.displayResultsAgain();
				return;
			}			
		}
		
		if (e.getSource()==butPrevDay || e.getSource()==butNextDay || e.getSource()==butNextSli || e.getSource()==butPrevSli 
				|| e.getSource()==butM0Map || e.getSource()==butT1Map || e.getSource()==butT2Map
				|| e.getSource()==butPrevEcho || e.getSource()==butNextEcho || e.getSource()==butPrevTr || e.getSource()==butNextTr  || e.getSource()==butGrayLut || e.getSource()==butFireLut ) {
			//Move t
			if((e.getSource()==butNextDay) && (tCor<nTimepoints-1) )tCor+=1;
			if((e.getSource()==butPrevDay) && (tCor>0) )tCor-=1;
	
			//Move z
			if((e.getSource()==butNextSli) && (zCor<dims[2]-1) )zCor+=1;
			if((e.getSource()==butPrevSli) && (zCor>0) )zCor-=1;
	
			//Move c of maps
			if(e.getSource()==butM0Map) cMaps=0;
			if(e.getSource()==butT1Map) cMaps=1;
			if(e.getSource()==butT2Map) cMaps=2;
	
			//Move c of echos
			if((e.getSource()==butNextEcho) && (cEchoes<hyperMRIT1T2.nChannels-1) ) {System.out.println("Changing cEcho from "+cEchoes+" to "+cEchoes+1+" because limit="+(hyperMRIT1T2.nChannels-1));cEchoes+=1;}
			if((e.getSource()==butPrevEcho) && (cEchoes>0) )cEchoes-=1;


			//Move c of Tr
			if((e.getSource()==butNextTr))changeTr(true);
			if((e.getSource()==butPrevTr))changeTr(false);


			
			//Move lut
			if(e.getSource()==butFireLut)isFireLut=true;
			if(e.getSource()==butGrayLut)isFireLut=false;
	
			if(e.getSource()==butPrevSli || e.getSource()==butNextSli) {
				actualizeMriObservationsBasedOnData();
				computeResultsAgain();
			}
			identifyRangedData();
			displayResultsAgain();			
		}			

		if (e.getSource()==butHide) {
			if(this.spectrumRangingModeT1==false && this.spectrumRangingModeT2 ==false) {
				this.spectrumRangingModeT1=true;
				xMouseRange=177;yMouseRange=(int) (18+(340.0*(this.nTimepoints-1-this.tCor)/(1.0*this.nTimepoints)) +20);
				System.out.println("Brute starting ranging coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
				this.rangeDisplayedFullBlocks=true;
			}
			else this.rangeDisplayedFullBlocks=!this.rangeDisplayedFullBlocks;
			butHide.setText(this.rangeDisplayedFullBlocks ? "<Html>Hide<br>ranging<Html>" : "<Html>Show<br>ranging<Html>");

			actualizeRangingBoundaries();
			this.identifyRangedData();
			displayResultsAgain();
		}
		
		
		
		if (e.getSource()==butRoi) {
			if(statusRoi==0) {				
				zoomLevel=imgCan1.getMagnification();
				IJ.log("Please select a ROI now, add it to the Roi manager, then hit the 'r' strike again");
				file=VitiDialogs.chooseOneRoiPathUI("Choose a .roi file", "");
				if(file==null)return;
				this.userRoi=new Opener().openRoi(file);
				this.userRoiInitial=new Opener().openRoi(file);
				IJ.selectWindow("MRI Curve explorer V2");
				statusRoi=2;
				this.userRoi.setPosition(0);				
				actualizeMriObservationsBasedOnData();
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
		
		if (e.getSource()==butNextRoi) {
			if( statusRoi!=2) {
				if(this.crossWidth<dims[0])this.crossWidth++;
				actualizeMriObservationsBasedOnData();
				computeResultsAgain();
				displayResultsAgain();
			}
		}
		if (e.getSource()==butPrevRoi) {
			if( statusRoi!=2) {
				if(this.crossWidth>0)this.crossWidth--;
				actualizeMriObservationsBasedOnData();
				computeResultsAgain();
				displayResultsAgain();
			}
		}
		if (e.getSource()==butNextRang) {
			if(this.spectrumRangingModeT1==false && this.spectrumRangingModeT2 ==false) {
				this.spectrumRangingModeT1=true;
				xMouseRange=177;yMouseRange=(int) (18+(340.0*(this.nTimepoints-1-this.tCor)/(1.0*this.nTimepoints)) +20);
				System.out.println("Brute starting ranging coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			}
			else {this.rangingFactor=this.rangingFactor*1.2;}

			actualizeRangingBoundaries();
			identifyRangedData();
			displayResultsAgain();
		}
		if (e.getSource()==butPrevRang) {
			if(this.spectrumRangingModeT1==false && this.spectrumRangingModeT2 ==false) {
				this.spectrumRangingModeT1=true;
				xMouseRange=177;yMouseRange=(int) (18+(340.0*(this.nTimepoints-1-this.tCor)/(1.0*this.nTimepoints)) +20);
				System.out.println("Brute starting ranging coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			}
			else {this.rangingFactor=this.rangingFactor*0.8;}
			this.actualizeRangingBoundaries();
			this.identifyRangedData();
			displayResultsAgain();
		}
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
	

	
	@Override
	public void mousePressed(MouseEvent e) {
	}
	
	@Override
	public void mouseReleased(MouseEvent e) {
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

    public void setupScroll(int ox, int oy) {
    	Rectangle srcRect = new Rectangle(0, 0, dims[0], dims[1]);
        xMouseStart = ox;
        yMouseStart = oy;
        xSrcStart = srcRect.x;
        ySrcStart = srcRect.y;
    }

    
    public void updatePlaying(int param,int action) {
		System.out.println("Update playing "+param+" , "+action);
    	int type=0;		int nParams=  3; 
		paramsTimelapseT1T2[tCor][0]=paramsTimelapseT1T2[tCor][0]*(param!=0 ? 1 : (action==0 ? 1.07 : (1/1.07)));
		paramsTimelapseT1T2[tCor][1]=paramsTimelapseT1T2[tCor][1]*(param!=1 ? 1 : (action==0 ? 1.07 : (1/1.07)));
		paramsTimelapseT1T2[tCor][2]=paramsTimelapseT1T2[tCor][2]*(param!=2 ? 1 : (action==0 ? 1.07 : (1/1.07)));

		if(param==3) {
			deltaT2+=(action==0 ? -1 : 1);
			params[7][3].setText(""+(((this.t1t2TeTimes[tCor][dims[2]-1][this.t1t2TeTimes[tCor][dims[2]-1].length-1]>12.1*16  || deltaT2<-10)? (zCor==0 ? 0 : 30) : 0)+deltaT2));
			params[6][3].setText("DeltaT2");
			this.timesT1T2Cute=getCuteTimesT1T2();
		}

		if(param==4) {
			deltaSigma+=(action==0 ? -10 : 10);
			params[7][2].setText(""+deltaSigma);
			params[6][2].setText("Delta sig");
			this.timesT1T2Cute=getCuteTimesT1T2();
		}

		
		
		System.out.println("Give : "+paramsTimelapseT1T2[tCor][0]+" , "+paramsTimelapseT1T2[tCor][1]+" , "+paramsTimelapseT1T2[tCor][2]+" , " +deltaT2);
		
		double []tabFitten=MRUtils.fittenRelaxationCurveT1T2(t1t2TrTimes[tCor][zCor],t1t2TeTimes[tCor][zCor],new double[] {paramsTimelapseT1T2[tCor][0],paramsTimelapseT1T2[tCor][1],paramsTimelapseT1T2[tCor][2]},hyperMRIT1T2.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma,MRUtils.T1T2_MONO_RICE);
		double[]accs=MRUtils.fittingAccuraciesT1T2(dataTimelapseT1T2[tCor],t1t2TrTimes[tCor][zCor],t1t2TeTimes[tCor][zCor],hyperMRIT1T2.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma,new double[] {paramsTimelapseT1T2[tCor][0],paramsTimelapseT1T2[tCor][1],paramsTimelapseT1T2[tCor][2]},MRUtils.T1T2_MONO_RICE,false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);		
		double [][]tabFittenCute=null;
		
		tabFittenCute=new double[this.timesT1T2Cute[tCor][zCor].length][];
		for(int fi=0;fi<this.timesT1T2Cute[tCor][zCor].length;fi++) {
			tabFittenCute[fi]=MRUtils.fittenRelaxationCurveT1T2(this.timesT1T2Cute[tCor][zCor][fi][0],this.timesT1T2Cute[tCor][zCor][fi][1],new double[] {paramsTimelapseT1T2[tCor][0],paramsTimelapseT1T2[tCor][1],paramsTimelapseT1T2[tCor][2]},hyperMRIT1T2.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma,MRUtils.T1T2_MONO_RICE);
		}
    	
   		tabFittenT1T2Mono[tCor]=tabFitten;
 		tabFittenT1T2MonoCute[tCor]=tabFittenCute;
 		khi2T1T2Mono[tCor]=accs[0];
 		pValT1T2Mono[tCor]=accs[1];
 		
 		displayResultsAgain();
    }
    
	public void mouseClicked(MouseEvent e) {
		for(int i=0;i<titles.length;i++) {
			if(e.getSource()==titles[i]) {
				switchT1T2=i%4;System.out.println("Changing to "+switchT1T2);displayResultsAgain();return;
			}
		}
		for(int i=0;i<params.length;i++)for(int j=0;j<params[i].length;j++) {
			if(e.getSource()==params[i][j]) {
				if(!playMode) {
					switchT1T2=(i/2)%4;System.out.println("Changing to "+switchT1T2);displayResultsAgain();return;
				}
				else {
					if(j==0 && (i%2)==1)updatePlaying(0,1);
					if(j==0 && (i%2)==0)updatePlaying(0,0);
					if(j==1 && (i%2)==1)updatePlaying(1,1);
					if(j==1 && (i%2)==0)updatePlaying(1,0);
					if(j==1 && (i>=6) && (i%2)==1)updatePlaying(2,1);
					if(j==1 && (i>=6) && (i%2)==0)updatePlaying(2,0);					
					if(j==3 && (i==6)) updatePlaying(3,1);
					if(j==3 && (i==7)) updatePlaying(3,0);
					if(j==2 && (i==6)) updatePlaying(4,1);
					if(j==2 && (i==7)) updatePlaying(4,0);
					return;
				}
			}
		}
		if(imgCan1.cursorOverImage()) currentCanvas=WIN_T1;
		if(imgCan2.cursorOverImage()) currentCanvas=WIN_T1;
		if(T1T2MixedEstimation) {
			if(plotCanT1T2.cursorOverImage()) {
				currentCanvas=WIN_PLOT1T2;
			}
		}
		else{
			if(plotCan1.cursorOverImage()) currentCanvas=WIN_PLOT1;
			if(plotCan2.cursorOverImage()) currentCanvas=WIN_PLOT2;
		}

		if(plotCan21.cursorOverImage()) currentCanvas=WIN_PLOT21;
		if(plotCan22.cursorOverImage()) currentCanvas=WIN_PLOT22;

		if(currentCanvas==WIN_T1 || currentCanvas==WIN_T2) {
			statusRoi=0;
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
			System.out.println("Click on MR images | Coordinates=("+ xMouse+","+yMouse+")"+"  |  zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();			
			displayResultsAgain();
			return;
		}
		else if(currentCanvas==WIN_PLOT1 || currentCanvas==WIN_PLOT2 || currentCanvas==WIN_PLOT1T2) {
			this.xMouseRangeCurve=e.getX();
			this.yMouseRangeCurve=e.getY();
			System.out.println("Click on Plot1 | Coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			int newEcho=actualizeSelectedEcho(currentCanvas);
			if(cEchoes!=newEcho) {
				cEchoes=newEcho;
				displayResultsAgain();
			}
			return;
		}
		else if(currentCanvas==WIN_PLOT21) {
			this.xMouseRange=e.getX();
			this.yMouseRange=e.getY();
			System.out.println("Click on Plot21 | Coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			spectrumRangingModeT1=true;	
			spectrumRangingModeT2=false;	
			int newTime=actualizeRangingBoundaries();
			if(tCor!=newTime)tCor=newTime;
			identifyRangedData();			
			displayResultsAgain();
			return;
		}

		else if(currentCanvas==WIN_PLOT22) {
			this.xMouseRange=e.getX();
			this.yMouseRange=e.getY();
			System.out.println("Click on Plot22 | Coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			spectrumRangingModeT2=true;	
			spectrumRangingModeT1=false;	
			int newTime=actualizeRangingBoundaries();
			if(tCor!=newTime)tCor=newTime;
			identifyRangedData();			
			displayResultsAgain();
			return;
		}

		else {
			displayResultsAgain();
		}
	}
	

	
	
	
	@Override
	public void keyReleased(KeyEvent e) {
		if(imgCan1.cursorOverImage()) currentCanvas=WIN_T1;
		if(imgCan2.cursorOverImage()) currentCanvas=WIN_T1;
		if(e.getKeyChar()==' ') {
			if(currentCanvas==WIN_T1)imgCan1.setSourceRect(imgCan2.getSrcRect());
			else imgCan2.setSourceRect(imgCan1.getSrcRect());
			displayResultsAgain();
		}
	}

	public void windowClosing(WindowEvent paramWindowEvent)
	  {
	    super.windowClosing(paramWindowEvent);
	    instance = null;
	  }




	@Override
	public void mouseWheelMoved(MouseWheelEvent e) {
		if(imgCan1.cursorOverImage()) currentCanvas=WIN_T1;
		if(imgCan2.cursorOverImage()) currentCanvas=WIN_T1;
		if(!T1T2MixedEstimation) {
			if(plotCan1.cursorOverImage()) currentCanvas=WIN_PLOT1;
			if(plotCan2.cursorOverImage()) currentCanvas=WIN_PLOT2;
		}
		else {
			if(plotCanT1T2.cursorOverImage()) currentCanvas=WIN_PLOT1T2;
		}
		if(plotCan21.cursorOverImage()) currentCanvas=WIN_PLOT21;
		if(plotCan22.cursorOverImage()) currentCanvas=WIN_PLOT22;
		//System.out.println(e);
		double init=imgCan1.getMagnification();		
		if(currentCanvas==WIN_T1 || currentCanvas==WIN_T2) {		
			//xMouse=e.getX();
			//yMouse=e.getY();
			if (e.getUnitsToScroll()<0) {
	//			if(statusRoi==2) {}//VitiDialogs.notYet("Warning : no zoom in / out during Roi mode");return;
				imgCan1.zoomIn(xMouse,yMouse);
				imgCan2.zoomIn(xMouse,yMouse);
				zoomLevel=imgCan1.getMagnification();
				xMouse=imgCan1.screenX(xCor);
				yMouse=imgCan1.screenY(yCor);
				xD=imgCan1.offScreenXD(xMouse);
				yD=imgCan1.offScreenYD(yMouse);
				displayResultsAgain();
			}
			else {
	//			if(statusRoi==2) {}//VitiDialogs.notYet("Warning : no zoom in / out during Roi mode");return;
				if(imgCan1.getMagnification()>(initMag*1.3333)) {
					imgCan1.setMagnification(imgCan1.getMagnification()/2.0);
					imgCan2.setMagnification(imgCan2.getMagnification()/2.0);
					xMouse=imgCan1.screenX(xCor);
					yMouse=imgCan1.screenY(yCor);
					imgCan1.zoomIn(xMouse,yMouse);
					imgCan2.zoomIn(xMouse,yMouse);
					xMouse=imgCan1.screenX(xCor);
					yMouse=imgCan1.screenY(yCor);
					zoomLevel=imgCan1.getMagnification();
				}
				else {
					imgCan1.setMagnification(initMag/1.333);
					imgCan2.setMagnification(initMag/1.333);
					xMouse=imgCan1.screenX(xCor);
					yMouse=imgCan1.screenY(yCor);
					imgCan1.zoomIn(xMouse,yMouse);
					imgCan2.zoomIn(xMouse,yMouse);
					xMouse=imgCan1.screenX(xCor);
					yMouse=imgCan1.screenY(yCor);
					zoomLevel=imgCan1.getMagnification();
				}
			}
		}
		if(currentCanvas==WIN_PLOT1 || currentCanvas==WIN_PLOT2 || currentCanvas==WIN_PLOT1T2) {		
			if (e.getUnitsToScroll()<0) { if (cEchoes<hyperMRIT1T2.nChannels-1)cEchoes+=1;}
			else  { if (cEchoes>0) cEchoes-=1;}
			identifyRangedData();
		}
		if(currentCanvas==WIN_PLOT21 || currentCanvas==WIN_PLOT22) {		
			if (e.getUnitsToScroll()<0) { if (tCor<nTimepoints-1)tCor+=1;}
			else  { if (tCor>0) tCor-=1;}
			identifyRangedData();
		}
		displayResultsAgain();
	}




	@Override
	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub

	}




	@Override
	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub 
		
	}

	

	
}

