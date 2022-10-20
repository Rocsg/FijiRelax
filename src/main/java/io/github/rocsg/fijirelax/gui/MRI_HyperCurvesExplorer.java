/*
 * 
 */
package io.github.rocsg.fijirelax.gui;


import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Polygon;
import java.awt.Rectangle;
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

//import com.jogamp.nativewindow.VisualIDHolder;
import io.github.rocsg.fijiyama.common.Timer;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import io.github.rocsg.fijiyama.registration.TransformUtils;
import io.github.rocsg.fijiyama.common.VitiDialogs;
import io.github.rocsg.fijiyama.common.VitimageUtils;
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
import io.github.rocsg.fijirelax.mrialgo.HyperMap;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;
import io.github.rocsg.fijirelax.mrialgo.RiceEstimator;

/**
 * MRI_HyperCurvesExplorer is the PlugInFrame describing the curve explorer of FijiRelax, an user-friendly tool for exploration of T1 T2 relaxation curves coming from T1 and T2 sequence
 */


public class MRI_HyperCurvesExplorer extends PlugInFrame implements ActionListener,KeyListener, MouseListener, MouseWheelListener, MouseMotionListener {/** The init timer. */
// 
	Timer initTimer;
	
	/** The has T 1. */
	boolean hasT1=false;
	
	/** The has T 2. */
	boolean hasT2=false;
	
	/** The sexy for bio info. */
	boolean sexyForBioInfo=true;
	
	/** The number bins. */
	int numberBins=150 /(1) ;//over 3 decades it makes one bin every 1.2%
	
	/** The sigma smoothing before estimation. */
	private double sigmaSmoothingBeforeEstimation=0.75;
	
	/** The timer. */
	private Timer timer;
	
	/** The x src start. */
	private int xSrcStart;
	
	/** The y src start. */
	private int ySrcStart;
	
	/** The roi is hidden. */
	private boolean roiIsHidden=false;
	
	/** The is bionano display. */
	private boolean isBionanoDisplay=false;//false;
	
	/** The but next echo. */
	JButton butNextEcho=new JButton("<Html>Next<br>echo<br> <Html>");
	
	/** The button switch mono. */
	JButton buttonSwitchMono=new JButton("<Html>Single T2 component<br>(Mono-exp)<br><Html>");
	
	/** The button switch bi. */
	JButton buttonSwitchBi=new JButton("<Html>Two T2 components<br>(Bi-exp)<br><Html>");
	
	/** The but prev echo. */
	JButton butPrevEcho=new JButton("<Html>Prev.<br>echo<br> <Html>");
	
	/** The but next tr. */
	JButton butNextTr=new JButton("<Html><br>Next<br>Tr<br> <Html>");
	
	/** The but prev tr. */
	JButton butPrevTr=new JButton("<Html><br>Prev.<br>Tr<br> <Html>");
	
	/** The but play. */
	JButton butPlay=new JButton("<Html><br>Play<br>with<br>prms<Html>");
	
	/** The but all curves. */
	JButton butAllCurves=new JButton("<Html><br>All<br>data<br> <Html>");

	/** The but M 0 map. */
	JButton butM0Map=new JButton("<Html> PD<br>map<Html>");
	
	/** The but T 1 map. */
	JButton butT1Map=new JButton("<Html> T1<br>map<Html>");
	
	/** The but T 2 map. */
	JButton butT2Map=new JButton("<Html> T2<br>map<Html>");

	/** The text exp. */
	JTextField textExp=new JTextField("Image exploration");
	
	/** The text tim. */
	JTextArea textTim=new JTextArea("              Observation day");
	
	/** The but prev day. */
	JButton butPrevDay=new JButton("Prev. day");
	
	/** The but next day. */
	JButton butNextDay=new JButton("Next day");
	
	/** The button exp. */
	JButton buttonExp=new JButton("Exp. fit");
	
	/** The button export. */
	JButton buttonExport=new JButton("Export results");
	
	/** The button off. */
	JButton buttonOff=new JButton("Off. fit");
	
	/** The button rice. */
	JButton buttonRice=new JButton("Rice. fit");
	
	/** The text sli. */
	JTextArea textSli=new JTextArea("               Current Z-slice");
	
	/** The but prev sli. */
	JButton butPrevSli=new JButton("Prev. slice");
	
	/** The but next sli. */
	JButton butNextSli=new JButton("Next slice");
	
	/** The but fire lut. */
	JButton butFireLut=new JButton("Fire LUT");
	
	/** The but gray lut. */
	JButton butGrayLut=new JButton("Gray LUT");
	
	/** The text luts. */
	JTextArea textLuts=new JTextArea("               Colormap");

	/** The text roi. */
	JTextField textRoi=new JTextField("       Region of interest");
	
	/** The text sam. */
	JTextArea textSam=new JTextArea("       Sample size (5x5x1)");
	
	/** The but prev roi. */
	JButton butPrevRoi=new JButton("Smaller");
	
	/** The but next roi. */
	JButton butNextRoi=new JButton("Larger");
	
	/** The text ran. */
	JTextArea textRan=new JTextArea("           Ranging interval");
	
	/** The but prev rang. */
	JButton butPrevRang=new JButton("Smaller");
	
	/** The but next rang. */
	JButton butNextRang=new JButton("Larger");
	
	/** The text misc. */
	JTextArea textMisc=new JTextArea("           Miscellaneous");
	
	/** The but hide. */
	JButton butHide=new JButton("<Html>Show<br>ranging<Html>");
	
	/** The but roi. */
	JButton butRoi=new JButton("<Html>Custom<br>ROI<Html>");

	/** The text info echoes. */
	JTextArea textInfoEchoes=new JTextArea("-------T2SEQ Day 28 TR=10000 ms TE=11ms------",45,1);
	
	/** The text info maps. */
	JTextArea textInfoMaps=new JTextArea("--------PDMAP Day 28 -----------",45,1);

	/** The play mode. */
	private boolean playMode=false;
	
	/** The switch T 1 T 2. */
	private int switchT1T2=0;
	
	/** The init mag. */
	private double initMag=0;
	
	/** The x D. */
	private double xD;
	
	/** The y D. */
	private double yD;
	
	/** The ij. */
	private ImageJ ij;
	
	/** The hyper map. */
	private HyperMap hyperMap;
	
	/** The dims. */
	private int[] dims;
	
	/** The voxs. */
	private double[] voxs;
	
	/** The n timepoints. */
	private int nTimepoints;
	
	/** The c echoes. */
	private int cEchoes;
	
	/** The c maps. */
	private int cMaps;
	
	/** The c T 1 T 2 tr index. */
	private int cT1T2TrIndex;
	
	/** The c T 1 T 2 te index. */
	private int cT1T2TeIndex;
	
	/** The n T 1 T 2 tr available. */
	private int nT1T2TrAvailable;
	
	/** The n T 1 T 2 te available. */
	private int [] nT1T2TeAvailable;
	
	/** The data timelapse T 1. */
	private double[][] dataTimelapseT1;
	
	/** The data timelapse T 2. */
	private double[][] dataTimelapseT2;
	
	/** The data timelapse T 1 T 2. */
	private double[][] dataTimelapseT1T2;
	
	/** The arbitrary limit between T 1 and T 2. */
	private double arbitraryLimitBetweenT1andT2=500;
	
	/** The norm value spectrum. */
	private final double normValueSpectrum=0.9;
	
	/** The maps image. */
	private ImagePlus mapsImage;
	
	/** The echoes image. */
	private ImagePlus echoesImage;
	
	/** The bubble color. */
	Color bubbleColor=new Color(220,220,0);
	
	/** The params unactivated. */
	Color paramsUnactivated=new Color(150,150,150);
	
	/** The params T 1 activated color. */
	Color paramsT1ActivatedColor=new Color(255,200,200);
	
	/** The params T 1 text activated color. */
	Color paramsT1TextActivatedColor=new Color(255,220,220);
	
	/** The params T 2 activated color. */
	Color paramsT2ActivatedColor=new Color(200,255,200);
	
	/** The params T 2 text activated color. */
	Color paramsT2TextActivatedColor=new Color(220,255,220);
	
	/** The cross color. */
	Color crossColor=new Color (60,60,255);
	
	/** The curve T 1 mono. */
	Color curveT1Mono=new Color (170,0,0);
	
	/** The curve T 1 bicomp. */
	Color curveT1Bicomp=new Color (255,70,70);
	
	/** The curve T 2 mono. */
	Color curveT2Mono=new Color (0,115,0);
	
	/** The curve T 2 bicomp. */
	Color curveT2Bicomp=new Color (50,205,50);

	/** The Constant maxT1. */
	private static final double maxT1=10500;
	
	/** The Constant maxT2. */
	private static final double maxT2=350;
	
	/** The Constant maxDisplayedT1. */
	private static final double maxDisplayedT1=4000;
	
	/** The Constant maxDisplayedT2. */
	private static final double maxDisplayedT2=120;
	
	/** The Constant t0T1. */
	private static final int t0T1=100;
	
	/** The Constant t1T1. */
	private static final int t1T1=10000;
	
	/** The Constant t0T2. */
	private static final int t0T2=10;
	
	/** The Constant t1T2. */
	private static final int t1T2=1000;
	
	/** The Constant WIN_T1. */
	private static final int WIN_T1=1;
	
	/** The Constant WIN_T2. */
	private static final int WIN_T2=2;
	
	/** The Constant WIN_PLOT1. */
	private static final int WIN_PLOT1=3;
	
	/** The Constant WIN_PLOT2. */
	private static final int WIN_PLOT2=4;
	
	/** The Constant WIN_PLOT1T2. */
	private static final int WIN_PLOT1T2=5;
	
	/** The Constant WIN_PLOT21. */
	private static final int WIN_PLOT21=6;
	
	/** The Constant WIN_PLOT22. */
	private static final int WIN_PLOT22=7;

	

	/** The range displayed full blocks. */
	//Parameters of graph and distribution estimation
	boolean rangeDisplayedFullBlocks=true;
	
	/** The ranging boundaries. */
	public double[]rangingBoundaries=new double[3];
	
	/** The ranging factor. */
	public double rangingFactor=0.2;
	
	/** The spectrum ranging mode T 1. */
	public boolean spectrumRangingModeT1=false;
	
	/** The spectrum ranging mode T 2. */
	public boolean spectrumRangingModeT2=false;
	
	/** The gaussian spectrum. */
	public int gaussianSpectrum=2;
	
	/** The separate normalization spectrum mode. */
	public int separateNormalizationSpectrumMode=1;//0 : all separated, 1= normalize along time, 2 : global normalization
	
	/** The multi threaded. */
	boolean multiThreaded=true;
	
	/** The sigma khi 2 weighted by number points. */
	boolean sigmaKhi2WeightedByNumberPoints=false;
	
	/** The rice estimator. */
	RiceEstimator riceEstimator;
	
	/** The sigma data rice lookup table T 1. */
	double[][][]sigmaDataRiceLookupTableT1;
	
	/** The sigma data rice lookup table T 2. */
	double[][][]sigmaDataRiceLookupTableT2;
	
	/** The sigma data rice lookup table T 1 T 2. */
	double[][][]sigmaDataRiceLookupTableT1T2;
	
	/** The fit algorithm. */
	private int fitAlgorithm=MRUtils.SIMPLEX;
	
	/** The computation rate. */
	int computationRate=100;
	
	/** The cur ec 1. */
	int curEc1;
	
	/** The cur ec 2. */
	int curEc2;

/** The gaussian weighting. */
	boolean gaussianWeighting=false;
	
	/** The max curves. */
	int maxCurves=170;
	
	/** The mean noise T 1 cur. */
	double meanNoiseT1Cur=1;
	
	/** The sigma noise T 1 cur. */
	double sigmaNoiseT1Cur=1;
	
	/** The mean noise T 2 cur. */
	double meanNoiseT2Cur=1;
	
	/** The sigma noise T 2 cur. */
	double sigmaNoiseT2Cur=1;
	
	/** The mean noise T 1 T 2 cur. */
	double meanNoiseT1T2Cur=1;
	
	/** The sigma noise T 1 T 2 cur. */
	double sigmaNoiseT1T2Cur=1;
	
	/** The n pts cur. */
	int nPtsCur=1;
	
	/** The thick curve. */
	int thickCurve=1;

	
	
	/** The x mouse range. */
	//Gui parameters : Window constants
	public int xMouseRange=0;
	
	/** The target width for plots t1t2. */
	private int TARGET_WIDTH_FOR_PLOTS_T1T2;
	
	/** The target width for plots alone. */
	private int TARGET_WIDTH_FOR_PLOTS_ALONE;
	
	/** The target height for plots and imgs. */
	private int TARGET_HEIGHT_FOR_PLOTS_AND_IMGS;
	
	/** The target width for plots and imgs. */
	private int TARGET_WIDTH_FOR_PLOTS_AND_IMGS;
	
	/** The dx img. */
	private int DX_IMG;
	
	/** The dy img. */
	private int DY_IMG;
	
	/** The dy text. */
	private int DY_TEXT;
	
	/** The delta x. */
	private int DELTA_X;
	
	/** The delta y. */
	private int DELTA_Y;
	
	/** The dx plot 1. */
	private int DX_PLOT_1;
	
	/** The dx plot 2. */
	private int DX_PLOT_2;
	
	/** The dy plot. */
	private int DY_PLOT;
	
	/** The total size X. */
	private int totalSizeX;  
	
	/** The total size Y. */
	private int totalSizeY;  
	
	/** The auto sizing time plots. */
	boolean autoSizingTimePlots=true;
	
	/** The is fire lut. */
	boolean isFireLut=true;
	

	/** The target number chars in texts. */
	//Params for Gui : Objects and texts to display
	int targetNumberCharsInTexts=50;
	
	/** The size based on image. */
	boolean sizeBasedOnImage=false;
	
	/** The debug display. */
	boolean debugDisplay=true;
	
	/** The user roi. */
	Roi userRoi;
	
	/** The r sentence. */
	String rSentence="";
	
	/** The py sentence. */
	String pySentence="";
	
	/** The mat sentence. */
	String matSentence="";
	
	/** The cross width. */
	int crossWidth=3;//3;//3;//3;
	
	/** The cross thick. */
	int crossThick=0;
	
	/** The x mouse. */
	int xMouse=0;
	
	/** The y mouse. */
	int yMouse=0;
	
	/** The x cor. */
	int xCor=1;
	
	/** The y cor. */
	int yCor=1;
	
	/** The t cor. */
	int tCor=0;
	
	/** The z cor. */
	int zCor=0;
	
	/** The plot T 1. */
	private Plot plotT1;
	
	/** The plot T 2. */
	private Plot plotT2;
	
	/** The plot T 1 T 2. */
	private Plot plotT1T2;
	
	/** The plot T 21. */
	private Plot plotT21;
	
	/** The plot T 22. */
	private Plot plotT22;
	
	/** The img can 1. */
	private ImageCanvas imgCan1;
	
	/** The img can 2. */
	private ImageCanvas imgCan2;
	
	/** The plot can 1. */
	private PlotCanvas plotCan1;
	
	/** The plot can 2. */
	private PlotCanvas plotCan2;
	
	/** The plot can 21. */
	private PlotCanvas plotCan21;
	
	/** The plot can 22. */
	private PlotCanvas plotCan22;
	
	/** The zoom level. */
	private double zoomLevel=0;
	
	/** The max plot YT 1. */
	private double maxPlotYT1;
	
	/** The max plot YT 2. */
	private double maxPlotYT2;
	
	/** The current canvas. */
	private int currentCanvas=1;
	
	/** The titles. */
	JTextField []titles=new JTextField[8];
	
	/** The params. */
	JTextField [][]params=new JTextField[16][5];
	
	/** The instance. */
	static PlugInFrame instance;
	
	/** The compute multi comp. */
	boolean computeMultiComp=true;
	
	/** The Constant serialVersionUID. */
	private static final long serialVersionUID = 1L;
	
	/** The x step cute T 1. */
	double xStepCuteT1=50;
	
	/** The x step cute T 2. */
	double xStepCuteT2=3;
	
	/** The t 1 times. */
	double [][][]t1Times;
	
	/** The t 2 times. */
	double [][][]t2Times;
	
	/** The t 1 t 2 times. */
	double [][][]t1t2Times;
	
	/** The times T 1 cute. */
	double[][][][][]timesT1Cute;
	
	/** The times T 2 cute. */
	double[][][][][]timesT2Cute;
	
	/** The times T 1 T 2 cute. */
	double[][][][][]timesT1T2Cute;

	/** The data timelapse full. */
	double[][][][]dataTimelapseFull;
	
	/** The points coords. */
	int [][]pointsCoords;
	
	/** The points estimations. */
	double[][][][]pointsEstimations;
	
	/** The data timelapse T 1 sigmas. */
	double[][]dataTimelapseT1Sigmas;
	
	/** The data timelapse T 2 sigmas. */
	double[][]dataTimelapseT2Sigmas;
	
	/** The data timelapse T 1 T 2 sigmas. */
	double[][]dataTimelapseT1T2Sigmas;

	/** The params timelapse T 1. */
	double[][]paramsTimelapseT1;
	
	/** The params timelapse T 2. */
	double[][]paramsTimelapseT2;
	
	/** The params timelapse T 1 T 2. */
	double[][]paramsTimelapseT1T2;
	
	/** The vals spectrum. */
	double[][][]valsSpectrum;
	
	/** The vals spectrum T 1. */
	double[][][]valsSpectrumT1;
	
	/** The vals spectrum T 2. */
	double[][][]valsSpectrumT2;

	/** The tab fitten T 1 mono. */
	//Updated just before plot being updated
	double [][]tabFittenT1Mono;
	
	/** The tab fitten T 1 mono cute. */
	double [][][]tabFittenT1MonoCute;
	
	/** The jitter T 1 mono. */
	double[]jitterT1Mono;
	
	/** The khi 2 T 1 mono. */
	double[]khi2T1Mono;
	
	/** The p val T 1 mono. */
	double[]pValT1Mono;

	/** The tab fitten T 1 bicomp. */
	double [][]tabFittenT1Bicomp;
	
	/** The tab fitten T 1 bicomp cute. */
	double [][][]tabFittenT1BicompCute;
	
	/** The jitter T 1 bicomp. */
	double[]jitterT1Bicomp;
	
	/** The khi 2 T 1 bicomp. */
	double[]khi2T1Bicomp;
	
	/** The p val T 1 bicomp. */
	double[]pValT1Bicomp;

	/** The tab fitten T 2 mono. */
	double [][]tabFittenT2Mono;
	
	/** The tab fitten T 2 mono cute. */
	double [][][]tabFittenT2MonoCute;
	
	/** The jitter T 2 mono. */
	double[]jitterT2Mono;
	
	/** The khi 2 T 2 mono. */
	double[]khi2T2Mono;
	
	/** The p val T 2 mono. */
	double[]pValT2Mono;

	/** The tab fitten T 2 bicomp. */
	double [][]tabFittenT2Bicomp;
	
	/** The tab fitten T 2 bicomp cute. */
	double [][][]tabFittenT2BicompCute;
	
	/** The jitter T 2 bicomp. */
	double[]jitterT2Bicomp;
	
	/** The khi 2 T 2 bicomp. */
	double[]khi2T2Bicomp;	
	
	/** The p val T 2 bicomp. */
	double[]pValT2Bicomp;

	/** The tab fitten T 1 T 2 mono. */
	double [][]tabFittenT1T2Mono;
	
	/** The tab fitten T 1 T 2 mono cute. */
	double [][][]tabFittenT1T2MonoCute;
	
	/** The jitter T 1 T 2 mono. */
	double[]jitterT1T2Mono;
	
	/** The khi 2 T 1 T 2 mono. */
	double[]khi2T1T2Mono;	
	
	/** The p val T 1 T 2 mono. */
	double[]pValT1T2Mono;

	/** The tab fitten T 1 T 2 bicomp. */
	double [][]tabFittenT1T2Bicomp;
	
	/** The tab fitten T 1 T 2 bicomp cute. */
	double [][][]tabFittenT1T2BicompCute;
	
	/** The jitter T 1 T 2 bicomp. */
	double[]jitterT1T2Bicomp;
	
	/** The khi 2 T 1 T 2 bicomp. */
	double[]khi2T1T2Bicomp;	
	
	/** The p val T 1 T 2 bicomp. */
	double[]pValT1T2Bicomp;

	/** The tab fitten T 1 bi T 2 bicomp. */
	double [][]tabFittenT1BiT2Bicomp;
	
	/** The tab fitten T 1 bi T 2 bicomp cute. */
	double [][][]tabFittenT1BiT2BicompCute;
	
	/** The jitter T 1 bi T 2 bicomp. */
	double[]jitterT1BiT2Bicomp;
	
	/** The khi 2 T 1 bi T 2 bicomp. */
	double[]khi2T1BiT2Bicomp;	
	
	/** The p val T 1 bi T 2 bicomp. */
	double[]pValT1BiT2Bicomp;

	/** The tab fitten T 1 T 2 bionano. */
	double [][]tabFittenT1T2Bionano;
	
	/** The tab fitten T 1 T 2 bionano cute. */
	double [][][]tabFittenT1T2BionanoCute;
	
	/** The jitter T 1 T 2 bionano. */
	double[]jitterT1T2Bionano;
	
	/** The khi 2 T 1 T 2 bionano. */
	double[]khi2T1T2Bionano;	
	
	/** The p val T 1 T 2 bionano. */
	double[]pValT1T2Bionano;

	/** The status roi. */
	int statusRoi;
	
	/** The correspondance canvas. */
	int[][]correspondanceCanvas;
	
	/** The range roi points. */
	private ArrayList<int[]> rangeRoiPoints;
	
	/** The last count 21. */
	private int lastCount21=0;
	
	/** The last count 22. */
	private int lastCount22=0;
	
	/** The last count 1. */
	private int lastCount1=0;
	
	/** The last count 2. */
	private int lastCount2=0;
	
	/** The y mouse range. */
	private int yMouseRange;
	
	/** The maps text. */
	private String[][][] mapsText;
	
	/** The echoes text. */
	private String[][][] echoesText;
	
	/** The x mouse range curve. */
	private int xMouseRangeCurve;
	
	/** The y mouse range curve. */
	private int yMouseRangeCurve;
	
	/** The x mouse drag. */
	private int xMouseDrag;
	
	/** The y mouse drag. */
	private int yMouseDrag;
	
	/** The x mouse start. */
	private int xMouseStart;
	
	/** The y mouse start. */
	private int yMouseStart;
	
	/** The file. */
	private String file;
	
	/** The user roi initial. */
	private Roi userRoiInitial;
	
	/** The T 1 T 2 mixed estimation. */
	private boolean T1T2MixedEstimation=false;
	
	/** The visual factor T 2. */
	private int visualFactorT2=50;
	
	/** The max T 1 T 2 on plot. */
	private double maxT1T2OnPlot=0;
	
	/** The max T 2 on plot. */
	private double maxT2OnPlot=0;
	
	/** The max T 1 on plot. */
	private double maxT1OnPlot=0;
	
	/** The t 1 t 2 tr times. */
	private double[][][] t1t2TrTimes;
	
	/** The t 1 t 2 te times. */
	private double[][][] t1t2TeTimes;
	
	/** The t 1 tr times. */
	private double[][][] t1TrTimes;
	
	/** The t 1 te times. */
	private double[][][] t1TeTimes;
	
	/** The t 2 tr times. */
	private double[][][] t2TrTimes;
	
	/** The t 2 te times. */
	private double[][][] t2TeTimes;
	
	/** The plot can T 1 T 2. */
	private PlotCanvas plotCanT1T2;
	
	/** The max plot YT 1 T 2. */
	private double maxPlotYT1T2;
	
	/** The times T 1 T 2 bubble cute. */
	private double[][] timesT1T2BubbleCute;
	
	/** The last count T 1 T 2. */
	private int lastCountT1T2;
	
	/** The last count T 1. */
	private int lastCountT1;
	
	/** The last count T 2. */
	private int lastCountT2;
	
	/** The is selected curve. */
	private boolean[][] isSelectedCurve;
	
	/** The is selected index. */
	private boolean[][] isSelectedIndex;
	
	/** The target height for up plots and imgs. */
	private int TARGET_HEIGHT_FOR_UP_PLOTS_AND_IMGS;
	
	/** The target height for down plot and imgs. */
	private int TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS;
	
	/** The params panel right TX. */
	private JPanel[][] paramsPanelRightTX;
	
	/** The switch all curves. */
	private boolean switchAllCurves=false;
	
	/** The mem switch. */
	private int memSwitch;
	
	/** The delta sigma. */
	private int deltaSigma=0;
	
	/** The memory bubble text X. */
	private double memoryBubbleTextX;
	
	/** The memory bubble text Y. */
	private double memoryBubbleTextY;
	
	/** The memory bubble text. */
	private String memoryBubbleText;
	
	/** The memory time click plot. */
	private double memoryTimeClickPlot=0;
	
	/** The double click delay. */
	private double doubleClickDelay=0.4;
	
	/** The is focus activated on spectrum. */
	private boolean isFocusActivatedOnSpectrum;
	
	/** The first draw bubbles. */
	private boolean firstDrawBubbles=true;
	
	/** The img roi. */
	private ImagePlus imgRoi;
	
	/** The sig smo bionano. */
	private String sigSmoBionano;
	
	/** The delta T 2. */
	private int deltaT2=0;
	
	/** The noise handling. */
	private int noiseHandling=2;
	
	/** The bin med. */
	private double[] binMed;
	
	/** The bin med T 1. */
	private double[] binMedT1;
	
	/** The normalize PD in middle tab. */
	private boolean normalizePDInMiddleTab=false;

		
	/**
	 *  Main function 
	 *  
	 * @param args the arguments
	 */
	public static void main(String[]args) {
		ImageJ ij=new ImageJ();
		//				new MRI_HyperCurvesExplorer().runExplorerFromHyperMap(new HyperMap(IJ.openImage("/home/fernandr/Bureau/testWithMapsT1.tif")));
	//new MRI_HyperCurvesExplorer().runExplorerFromHyperMap(new HyperMap(IJ.openImage("/home/fernandr/Bureau/Traitements/Sorgho/Series_temporelles/All_timeseries/SSM1_Timeseries.tif")));
//		new MRI_HyperCurvesExplorer().runExplorerFromHyperMap(new HyperMap(IJ.openImage("/home/fernandr/Bureau/test.tif")));
//		new MRI_HyperCurvesExplorer().runExplorerFromHyperMap(new HyperMap(IJ.openImage("/home/fernandr/Bureau/test.tif")));
		ImagePlus img=IJ.openImage("/home/fernandr/Bureau/FijiRelax_DOI/Tutorial_01_First_steps/Input_data/Images/HyperMap_simulated_fantoms.tif");
//		ImagePlus img=IJ.openImage("/home/fernandr/Bureau/test.tif");
		HyperMap hyp=new HyperMap(img);
		hyp.computeMaps();
		new MRI_HyperCurvesExplorer().runExplorerFromHyperMap(hyp);
		if(true)return;
	}
	
	/**
	 * Instantiates a new MRI hyper curves explorer.
	 */
	public MRI_HyperCurvesExplorer() {
		super("Vitimage MRI Water tracker ");
		initTimer=new Timer();
	}
	
	/**
	 * Run the explorer
	 *
	 * @param arg the arg
	 */
	public void run(String arg) {
		System.out.println("RUN CALLED !");
		runExplorer(null);
	}
	
	/**
	 * Run explorer by specifying an Hypermap to explore
	 *
	 * @param hyp the hyp
	 */
	public void runExplorerFromHyperMap(HyperMap hyp) {
		runExplorerFromImage(hyp.getAsImagePlus());
		VitiDialogs.getYesNoUI("Curve explorer help",getHelpString());	
	}
	
	/**
	 * Run explorer by specifying an HyperMap to explore
	 *
	 * @param imgPath the path to the HyperMap
	 */
	public void runExplorer(String imgPath) {
		ImagePlus fullHyp=null;
		if (imgPath!=null)fullHyp=IJ.openImage(imgPath);
		else fullHyp=VitiDialogs.chooseOneImageUI("Open hyperimage", "Open a hyperimage built using T1T2MapImporter");		
		runExplorerFromImage(fullHyp);
	}
		
	/**
	 * Run explorer from a ImagePlus containing an HyperMap
	 *
	 * @param fullHyp the ImagePlus
	 */
	public void runExplorerFromImage(ImagePlus fullHyp) {
		if(!(fullHyp.getType()==ImagePlus.GRAY32))IJ.run(fullHyp,"32-bit","");
		if(VitimageUtils.isBouture(fullHyp))this.deltaT2=25;
		IJ.log("Opening hyperimage "+VitimageUtils.imageResume(fullHyp));
		hyperMap=new HyperMap(fullHyp);
		if(!hyperMap.hasMaps)hyperMap.computeMaps();
		this.hasT1=hyperMap.hasT1sequence;
		this.hasT2=hyperMap.hasT2sequence;
		
		if(new File("/users/bionanonmri/").exists()) {
			this.isBionanoDisplay=true;
			IJ.showMessage("Detected Bionano server. \nSurvival display, but numerous cores");
		}
		if(hyperMap.hasT1T2sequence) {T1T2MixedEstimation=true;separateNormalizationSpectrumMode=2;}
		
		Object[]obj=null;
		mapsImage=hyperMap.getMapsImage();
		mapsText=null;
		echoesImage=hyperMap.getEchoesImage();
		echoesText=hyperMap.getEchoesImageText();
		dims=hyperMap.dims;
		voxs=hyperMap.voxs;
		/*if(ij==null)ij = new ImageJ();https://imagej.github.io/
		else ij = IJ.getInstance();*/
		IJ.log("Starting MRI Curve Explorer");
		xCor=hyperMap.dims[0]/2;
		yCor=hyperMap.dims[1]/2;
		zCor=0;
		tCor=0;
		nTimepoints=hyperMap.T;
		//setupTimesTrTe();
		setupStructures();
		startGui();
	}
	

	/**
	 * Setup the structures of the GUI.
	 */
	public void setupStructures() {
		timer=new Timer();
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
		this.tabFittenT1MonoCute=new double[this.nTimepoints][][];
		this.tabFittenT1BicompCute=new double[this.nTimepoints][][];
		this.tabFittenT2Mono=new double[this.nTimepoints][];
		this.tabFittenT2Bicomp=new double[this.nTimepoints][];
		this.tabFittenT2MonoCute=new double[this.nTimepoints][][];
		this.tabFittenT2BicompCute=new double[this.nTimepoints][][];
		this.tabFittenT1T2Mono=new double[this.nTimepoints][];
		this.tabFittenT1T2Bicomp=new double[this.nTimepoints][];
		this.tabFittenT1BiT2Bicomp=new double[this.nTimepoints][];
		this.tabFittenT1T2Bionano=new double[this.nTimepoints][];
		this.tabFittenT1T2MonoCute=new double[this.nTimepoints][][];
		this.tabFittenT1T2BicompCute=new double[this.nTimepoints][][];
		this.tabFittenT1BiT2BicompCute=new double[this.nTimepoints][][];
		this.tabFittenT1T2BionanoCute=new double[this.nTimepoints][][];
		if(T1T2MixedEstimation) {
			this.timesT1T2Cute=getCuteTimesT1T2();
		}
		else {
			if(hasT1) {
				this.timesT1Cute=getCuteTimesT1();				
			}
			if(hasT2) {
				this.timesT2Cute=getCuteTimesT2();				
			}
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
	
	/**
	 * Gets the cute times T1 T2 which are an interpolation of the actual MRI data, in order to generate smooth curves corresponding to estimate exponential parameters
	 *
	 * @return the cute times T1 T2
	 */
	public double[][][][][]getCuteTimesT1T2(){
		double[][][][]t1t2trte=hyperMap.getT1T2TrTeTimes();// [Time][Curve (T1, puis T2 succ][Choice : Tr, Te or Visual T][actual values]
		this.t1t2TrTimes=this.hyperMap.getT1T2TrTimes();
		
		this.t1t2TeTimes=this.hyperMap.getT1T2TeTimes();
		this.t1t2Times=new double[nTimepoints][dims[2]][];
		if(!T1T2MixedEstimation) {
			//this.t1Times=this.hyperMRIT1T2.getT1TrTimes();
			//this.t2Times=this.hyperMRIT1T2.getT2TeTimes();
		}
		double[][][][][]tabRet=new double[nTimepoints][dims[2]][][][];
		isSelectedCurve=new boolean[nTimepoints][];
		isSelectedIndex=new boolean[nTimepoints][];
		for(int t=0;t<nTimepoints;t++) {
			for(int z=0;z<dims[2];z++) {
				//Collect all the possible transversal relaxations curves, that are successive points having the same Tr, but with varying Te
				int nt1=0;
				double memTr=-1;double memTe=-1;
				boolean isUp=false;
				this.isSelectedCurve[t]=new boolean[t1t2trte[t][z].length];//Enough space in fact, whereas it is not exactly the good number, but it is greater...
				for(int c=0;c<t1t2trte[t][z].length;c++) {
					if(t1t2trte[t][z][c][0]!=memTr) {
						isSelectedCurve[t][nt1]=isUp;isUp=!isUp;memTr=t1t2trte[t][z][c][0];nt1++;
					}
					if(t1t2trte[t][z][c][0]>9000) {
						t1t2trte[t][z][c][1]+=this.deltaT2;
						t1t2TeTimes[t][z][c]+=this.deltaT2;
					}
				}
				isSelectedCurve[t][nt1]=true;
				double[][]valuesTrTe=new double[nt1][2];
				tabRet[t][z]=new double[nt1+1][3][];
				this.isSelectedIndex[t]=new boolean[t1t2trte[t][0].length];
				this.t1t2Times[t][z]=new double[t1t2trte[t][z].length];
				nt1=-1;memTr=-1;
			
				for(int c=0;c<t1t2trte[t][z].length;c++) {
					t1t2Times[t][z][c]=t1t2trte[t][z][c][0] + visualFactorT2 * t1t2trte[t][z][c][1];

					if(t1t2trte[t][z][c][0]!=memTr) {memTr=t1t2trte[t][z][c][0];nt1++;valuesTrTe[nt1][0]=memTr;}
					isSelectedIndex[t][c]=isSelectedCurve[t][nt1+1];
					valuesTrTe[nt1][1]=t1t2trte[t][z][c][1];// Fatally, this will be the bigger (last) Te at the end
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
					for(int j=0;j<tabRet[t][z][i][0].length;j++) {
						tabRet[t][z][i][2][j]=   tabRet[t][z][i][0][j]     +   visualFactorT2 * tabRet[t][z][i][1][j];
						if(tabRet[t][z][i][2][j]>maxT1T2OnPlot)maxT1T2OnPlot=tabRet[t][z][i][2][j];
					}
				}
			}
		}
		//System.exit(0);
		return tabRet;
	}

	/**
	 * Gets the cute times T2 which are an interpolation of the actual MRI data, in order to generate smooth curves corresponding to estimate exponential parameters
	 *
	 * @return the cute times T2
	 */
	public double[][][][][]getCuteTimesT2(){
		double[][][][]t2trte=hyperMap.getT2TrTeTimes();// [Time][Curve (T1, puis T2 succ][Choice : Tr, Te or Visual T][actual values]
		this.t2TrTimes=this.hyperMap.getT2TrTimes();
		this.t2TeTimes=this.hyperMap.getT2TeTimes();
//		this.t2Times=new double[nTimepoints][dims[2]][];
//		this.t1Times=this.hyperMap.getT1TrTimes();
		this.t2Times=this.hyperMap.getT2TeTimes();
		double[][][][][]tabRet=new double[nTimepoints][dims[2]][][][];
		isSelectedCurve=new boolean[nTimepoints][];
		isSelectedIndex=new boolean[nTimepoints][];
		for(int t=0;t<nTimepoints;t++) {
			for(int z=0;z<dims[2];z++) {
				//Collect all the possible transversal relaxations curves, that are successive points having the same Tr, but with varying Te
				int nt1=0;
				double memTr=-1;double memTe=-1;
				boolean isUp=false;
				this.isSelectedCurve[t]=new boolean[t2trte[t][z].length+1];//Enough space in fact, whereas it is not exactly the good number, but it is greater...
				for(int c=0;c<t2trte[t][z].length;c++) {
					if(t2trte[t][z][c][0]!=memTr) {
						isSelectedCurve[t][nt1]=isUp;isUp=!isUp;memTr=t2trte[t][z][c][0];nt1++;
					}
					if(t2trte[t][z][c][0]>9000) {
						t2trte[t][z][c][1]+=this.deltaT2;
						t2TeTimes[t][z][c]+=this.deltaT2;
					}
				}
				isSelectedCurve[t][nt1]=true;
				double[][]valuesTrTe=new double[nt1][2];
				tabRet[t][z]=new double[nt1][3][];
				this.isSelectedIndex[t]=new boolean[t2trte[t][0].length];
				this.t2Times[t][z]=new double[t2trte[t][z].length];
				nt1=-1;memTr=-1;
			
				for(int c=0;c<t2trte[t][z].length;c++) {
					t2Times[t][z][c]=t2trte[t][z][c][1];//DEBUG t2trte[t][z][c][0] + visualFactorT2 * t2trte[t][z][c][1];

					if(t2trte[t][z][c][0]!=memTr) {memTr=t2trte[t][z][c][0];nt1++;valuesTrTe[nt1][0]=memTr;}
					isSelectedIndex[t][c]=isSelectedCurve[t][nt1+1];
					valuesTrTe[nt1][1]=t2trte[t][z][c][1];// Fatally, this will be the bigger (last) Te at the end
				}
				
			
				//Add the other curves : the transversal relaxations, constituted of successive growing Tr curves, with for each Te going from 0 to the max Te at this level
				for(int cu=0;cu<valuesTrTe.length;cu++) {
					int nPts=(int)Math.ceil(valuesTrTe[cu][1]*1.0/xStepCuteT2+5);
					tabRet[t][z][cu]=new double[3][nPts];
					for(int pt=0;pt<nPts;pt++) {
						tabRet[t][z][cu][1][pt]=pt*xStepCuteT2+0.0001;
						tabRet[t][z][cu][0][pt]=valuesTrTe[cu][0];
					}
				}
				for(int i=0;i<tabRet[t][z].length;i++) {	
					for(int j=0;j<tabRet[t][z][i][0].length;j++) {
						tabRet[t][z][i][2][j]= tabRet[t][z][i][1][j];//DEBUG  tabRet[t][z][i][0][j]     +   visualFactorT2 * tabRet[t][z][i][1][j];
						if(tabRet[t][z][i][2][j]>maxT2OnPlot)maxT2OnPlot=tabRet[t][z][i][2][j];
					}
				}
			}
		}
		return tabRet;
	}

	
	
	/**
	 * Gets the cute times T1 which are an interpolation of the actual MRI data, in order to generate smooth curves corresponding to estimate exponential parameters
	 *
	 * @return the cute times T 1
	 */
	public double[][][][][]getCuteTimesT1(){
		double[][][][]t1trte=hyperMap.getT1TrTeTimes();// [Time][Curve (T1, puis T2 succ][Choice : Tr, Te or Visual T][actual values]
		this.t1TrTimes=this.hyperMap.getT1TrTimes();	
		this.t1TeTimes=this.hyperMap.getT1TeTimes();
//		this.t1Times=new double[nTimepoints][dims[2]][];
		this.t1Times=this.hyperMap.getT1TrTimes();
		
		double[][][][][]tabRet=new double[nTimepoints][dims[2]][][][];
		for(int t=0;t<nTimepoints;t++) {
			for(int z=0;z<dims[2];z++) {
				/*collect all the possible transversal relaxations curves, that are successive points having the same Tr, but with varying Te
				int nt1=0;
				double memTr=-1;double memTe=-1;
				boolean isUp=false;
				this.isSelectedCurve[t]=new boolean[t1trte[t][z].length];//Enough space in fact, whereas it is not exactly the good number, but it is greater...
				for(int c=0;c<t1trte[t][z].length;c++) {
					if(t1trte[t][z][c][0]!=memTr) {
						isSelectedCurve[t][nt1]=isUp;isUp=!isUp;memTr=t1trte[t][z][c][0];nt1++;
					}
				}*/
				double[][]valuesTrTe=new double[1][2];
				tabRet[t][z]=new double[1][3][];
				this.t1Times[t][z]=new double[t1trte[t][z].length];
				
				for(int c=0;c<t1trte[t][z].length;c++) {
					t1Times[t][z][c]=t1trte[t][z][c][0];//DEBUG + visualFactorT2 * t1trte[t][z][c][1];
				}
				
				int nPts=(int)Math.ceil(t1Times[t][z][t1trte[t][z].length-1]*1.0/xStepCuteT1+1);
				tabRet[t][z][0]=new double[3][nPts];
				for(int pt=0;pt<nPts;pt++) {
					tabRet[t][z][0][0][pt]=pt*xStepCuteT1;
					tabRet[t][z][0][1][pt]=0;
				} 
			
				//Calculate point on graph
				for(int i=0;i<tabRet[t][z].length;i++) {	
					for(int j=0;j<tabRet[t][z][i][0].length;j++) {
						tabRet[t][z][i][2][j]=   tabRet[t][z][i][0][j] ;//DEBUG   +   visualFactorT2 * tabRet[t][z][i][1][j];
						if(tabRet[t][z][i][2][j]>maxT1OnPlot)maxT1OnPlot=tabRet[t][z][i][2][j];
					}
				}
			}
		}
		//System.exit(0);
		return tabRet;
	}

	
	
	
	/**
	 * Selecting which T1/T2 data should be displayed
	 *
	 * @param dat the dat
	 * @param time the time
	 * @return the double[]
	 */
	public double[]t1t2selectCute(double[]dat,int time){
		int n=0;
		for(int i=0;i<dat.length && i<isSelectedIndex[time].length;i++)if(switchAllCurves || isSelectedIndex[time][i])n++;
		double[]ret=new double[n];n=0;
		for(int i=0;i<dat.length && i<isSelectedIndex[time].length;i++)if(switchAllCurves || isSelectedIndex[time][i])ret[n++]=dat[i];
		return ret;
	}
	
	/**
	 * Getting if T1/T2 data should be displayed
	 *
	 * @param tran the tran
	 * @param time the time
	 * @return true, if successful
	 */
	public boolean t1t2SelectedCurve(int tran,int time) {
		return (switchAllCurves || isSelectedCurve[time][tran]);
	}
	
	/**
	 * Getting if T1/T2 data should be displayed
	 *
	 * @param tran the tran
	 * @param time the time
	 * @return true, if successful
	 */
	public boolean t1SelectedCurve(int tran,int time) {
		return (true);
	}

	/**
	 * Getting if T1/T2 data should be displayed
	 *
	 * @param dat the dat
	 * @param time the time
	 * @return the double[]
	 */
	public double[]t1selectCute(double[]dat,int time){
		int n=0;
		for(int i=0;i<dat.length && i<isSelectedIndex[time].length;i++)if(switchAllCurves || isSelectedIndex[time][i])n++;
		double[]ret=new double[n];n=0;
		for(int i=0;i<dat.length && i<isSelectedIndex[time].length;i++)if(switchAllCurves || isSelectedIndex[time][i])ret[n++]=dat[i];
		return ret;
	}
	
	/**
	 * Getting if T1/T2 data should be displayed
	 *
	 * @param tran the tran
	 * @param time the time
	 * @return true, if successful
	 */
	public boolean t2SelectedCurve(int tran,int time) {
		return (switchAllCurves || isSelectedCurve[time][tran]);
	}
	
	/**
	 * Getting if T1/T2 data should be displayed
	 *
	 * @param dat the dat
	 * @param time the time
	 * @return the double[]
	 */
	public double[]t2selectCute(double[]dat,int time){
		int n=0;
		for(int i=0;i<dat.length && i<isSelectedIndex[time].length;i++)if(switchAllCurves || isSelectedIndex[time][i])n++;
		double[]ret=new double[n];n=0;
		for(int i=0;i<dat.length && i<isSelectedIndex[time].length;i++)if(switchAllCurves || isSelectedIndex[time][i])ret[n++]=dat[i];
		return ret;
	}
	

	
	
	/**
	 * Start the GUI.
	 */
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
		updateSwitchButtons();
	}


		

	
	
	
	
	
	
	
	
	
	/**
	 *  Calling the sequence of GUI building helpers--------------------------------------------------------------------------------.
	 */	
	public void repaintAll() {
		imgCan1.repaint();
		imgCan2.repaint();
		if(T1T2MixedEstimation)plotCanT1T2.repaint();
		else{ plotCan1.repaint();plotCan2.repaint();}
		plotCan21.repaint();
		plotCan22.repaint();
	}
		
	/**
	 * Initialize screen constants.
	 */
	public void initializeScreenConstants() {
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		int screenY=(int) Math.round(screenSize.getHeight());
		int screenX=(int) Math.round(screenSize.getWidth());
		if(screenX>1920)screenX=1920;
		if(screenX<=1600 && false) {
			IJ.showMessage("Your screen has a very low resolution : "+screenX+" X "+screenY+"\nPlease consider investing in one which have at least 1024 lines x 1480 columns.\nPlugin will run in survivor mode, thus everything can happen");
			isBionanoDisplay=true;
		}
		if(isBionanoDisplay) {screenX=1260; screenY=768;}
		IJ.log("Screen resolution : "+screenX+" X "+screenY);
		int DY_TITLE=isBionanoDisplay ? 5 : 10;
		TARGET_HEIGHT_FOR_PLOTS_AND_IMGS=isBionanoDisplay ? 233 : 400;
		TARGET_HEIGHT_FOR_UP_PLOTS_AND_IMGS=isBionanoDisplay ? 280 : 480;
		TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS=isBionanoDisplay ? 230 : 300;
		TARGET_WIDTH_FOR_PLOTS_AND_IMGS=isBionanoDisplay ? 235 : 400;
		TARGET_WIDTH_FOR_PLOTS_T1T2=isBionanoDisplay ? 800 : 1300;
		TARGET_WIDTH_FOR_PLOTS_ALONE=isBionanoDisplay ? 400 : 650;
		DX_IMG=TARGET_WIDTH_FOR_PLOTS_AND_IMGS;
		DY_IMG=TARGET_HEIGHT_FOR_PLOTS_AND_IMGS;
		DY_TEXT=isBionanoDisplay ? 18 : 30;
		DELTA_X=isBionanoDisplay ? 5 : 10;
		DELTA_Y=isBionanoDisplay ? 5 : 10;
		totalSizeY= 2*DY_IMG+DY_TITLE;  
		totalSizeX= DX_IMG+DX_PLOT_1+DX_PLOT_2+2*DELTA_X ;  
	}

	/**
	 * Initialize GUI.
	 */
	public void initializeGUI() {
		String name=butNextEcho.getFont().getFamily();
		if(isBionanoDisplay) {
			Font fontGui=new Font(name,0,9);
			butNextEcho.setFont(fontGui);
			butPrevEcho.setFont(fontGui);
			butNextTr.setFont(fontGui);
			butPrevTr.setFont(fontGui);
			butPlay.setFont(fontGui);
			butAllCurves.setFont(fontGui);
	
			butM0Map.setFont(fontGui);
			butT1Map.setFont(fontGui);
			butT2Map.setFont(fontGui);
	
			textExp.setFont(fontGui);
			textTim.setFont(fontGui);
			butPrevDay.setFont(fontGui);
			butNextDay.setFont(fontGui);
			textSli.setFont(fontGui);
			butPrevSli.setFont(fontGui);
			butNextSli.setFont(fontGui);
			butFireLut.setFont(fontGui);
			butGrayLut.setFont(fontGui);
			textLuts.setFont(fontGui);
	
			buttonSwitchBi.setFont(fontGui);
			buttonSwitchMono.setFont(fontGui);
			buttonExp.setFont(fontGui);
			buttonExport.setFont(fontGui);
			buttonOff.setFont(fontGui);
			buttonRice.setFont(fontGui);
			textRoi.setFont(fontGui);
			textSam.setFont(fontGui);
			butPrevRoi.setFont(fontGui);
			butNextRoi.setFont(fontGui);
			textRan.setFont(fontGui);
			butPrevRang.setFont(fontGui);
			butNextRang.setFont(fontGui);
			textMisc.setFont(fontGui);
			butHide.setFont(fontGui);
			butRoi.setFont(fontGui);
	
			textInfoEchoes.setFont(fontGui);
			textInfoMaps.setFont(fontGui);

		}		
		ImageJ ij = IJ.getInstance();
		//Prepare Image T1
		imgCan1.setMagnification(Math.min(1.0*DX_IMG/dims[0],1.0*DY_IMG/dims[1]));
        imgCan1.addKeyListener(this);
        imgCan1.addMouseListener(this);
        imgCan1.addMouseMotionListener(this);
        imgCan1.addMouseWheelListener(this);
        imgCan1.setSize(TARGET_WIDTH_FOR_PLOTS_AND_IMGS, TARGET_HEIGHT_FOR_PLOTS_AND_IMGS);

        //Prepare Image T2
//		imgCan2.setMagnification(1.0*DX_IMG/dims[0]);
		imgCan2.setMagnification(Math.min(1.0*DX_IMG/dims[0],1.0*DY_IMG/dims[1]));
		initMag=imgCan2.getMagnification();
        imgCan2.addKeyListener(this);
        imgCan2.addMouseListener(this);
        imgCan2.addMouseMotionListener(this);
        imgCan2.addMouseWheelListener(this);
        imgCan2.setSize(TARGET_WIDTH_FOR_PLOTS_AND_IMGS, TARGET_HEIGHT_FOR_PLOTS_AND_IMGS);
 		
        //Prepare Texts
        int totalY=3*DY_TEXT+2*DELTA_Y;
        int totalX=DX_PLOT_1+DELTA_X+DX_PLOT_2;
        int []x0=new int[8]; 
        int []y0=new int[8];
        for(int i=0;i<8;i++){
        	x0[i]=DX_IMG+DELTA_X + (i * totalX)/8;
        	y0[i]= DY_PLOT+DELTA_Y + (i * totalY)/8;
       }
        String []strTitles=new String[] {" T1 mono", " T1 mono" ," T1 bi"," Bionano ","   T2 mono","   T2 bi","   T2 bi","   Bionano"};
        String [][]strParams=new String[][] {
        	{"PD "  , "T1 (ms)" , "" ,  "" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"PD "  , "T1 (ms)" , "" ,  "" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"Short PD"  , "1st T1" , "Long PD" ,  "2nd T1" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"PD "  , "T1" , "" ,  "" , "" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"PD "  , "T2 (ms)" , "" ,  "" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"Short PD"  , "Short T2" , "Long PD" ,  "Long T2" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"Short PD"  , "Short T2" , "Long PD" ,  "Long T2" , "Khi2/N" },
        	{"-"  , "-" , "" ,  "" , "-" },
        	{"PD "  , "T2" , "" ,  "" , "" },
        	{"-"  , "-" , "" ,  "" , "-" }
        };
        for(int i=0;i<8;i++) {
        	titles[i]=new JTextField(strTitles[i],0);
        	titles[i].setEditable(false);titles[i].setFont(new Font("Helvetica", Font.BOLD, isBionanoDisplay ? 10 : 13));titles[i].setBackground(paramsUnactivated);
    	}
        for(int i=0;i<16;i++)for(int j=0;j<5;j++) {
        	params[i][j]=new JTextField(strParams[i][j],0);
        	params[i][j].setEditable(false);params[i][j].setFont(new Font( "Helvetica" , i%2==0 ? Font.BOLD : 0, isBionanoDisplay ? 10 : 13));params[i][j].setBackground(paramsUnactivated);
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
    	textInfoEchoes.setFont(new Font("Helvetica", Font.PLAIN, isBionanoDisplay ? 10 : 15));
    	textInfoEchoes.setText(this.echoesText[cEchoes][zCor][tCor]);
		JPanel butEchoesPan=new JPanel();
		butEchoesPan.setLayout(new GridBagLayout());
		butEchoesPan.setBorder(BorderFactory.createEmptyBorder(0,isBionanoDisplay ? 0 :5,0,0));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.NONE;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        butEchoesPan.add(butNextEcho,gbc);
        gbc.gridy=1;        gbc.gridheight=1;     
        butEchoesPan.add(new JButton(""),gbc);
        gbc.gridy=2;        gbc.gridheight=1;     
        butEchoesPan.add(butPrevEcho,gbc);
        if(T1T2MixedEstimation) {
    		gbc.fill = GridBagConstraints.NONE;
            gbc.gridy=2;        gbc.gridheight=1;     
            text=new JTextArea();
            text.setFont(new Font("Arial",0,isBionanoDisplay ? 12 : 20));
    		butEchoesPan.add(text,gbc);
            gbc.gridy=3;        gbc.gridheight=1;     
            text=new JTextArea();
            text.setFont(new Font("Arial",0,isBionanoDisplay ? 7 : 10));
    		butEchoesPan.add(text,gbc);
	        gbc.gridy=4;        gbc.gridheight=1;     
	        butEchoesPan.add(butNextTr,gbc);
	        gbc.gridy=5;        gbc.gridheight=1;     
	        butEchoesPan.add(butPrevTr,gbc);
            gbc.gridy=6;        gbc.gridheight=1;     
            text=new JTextArea();
            text.setFont(new Font("Arial",0,isBionanoDisplay ? 12 : 20));
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
    	echoesPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :5,isBionanoDisplay ? 0 :5,isBionanoDisplay ? 5 :5,isBionanoDisplay ? 0 :5));
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
    	textInfoMaps.setFont(new Font("Helvetica", Font.PLAIN, isBionanoDisplay ? 10 : 15));
    	textInfoMaps.setText("");//this.mapsText[cEchoes][zCor][tCor]
		JPanel butMapsPan=new JPanel();
		butMapsPan.setLayout(new GridBagLayout());
		butMapsPan.setBorder(BorderFactory.createEmptyBorder(0,isBionanoDisplay ? 0 :5,0,0));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        butMapsPan.add(butM0Map,gbc);
        gbc.gridy=1;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,isBionanoDisplay ? 7 : 10));
		gbc.fill = GridBagConstraints.NONE;
		butMapsPan.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridy=2;        gbc.gridheight=1;     
        butMapsPan.add(butT1Map,gbc);
        boolean survivorT2=false;
        if((!hasT1) && (!T1T2MixedEstimation))survivorT2=true;
        gbc.gridy=3;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,isBionanoDisplay ? 7 : 10));
		gbc.fill = GridBagConstraints.NONE;
		butMapsPan.add(text,gbc);
		gbc.fill = GridBagConstraints.BOTH;
		gbc.gridy=4;        gbc.gridheight=1;     
        if(survivorT2)butT1Map.setText("<Html>T2<br>map<br> <Html>");
        else if(hasT2)butMapsPan.add(butT2Map,gbc);
        gbc.gridy=5;        gbc.gridheight=1;     
        text=new JTextArea();
        text.setFont(new Font("Arial",0,isBionanoDisplay ? 30 : 50));
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
    	mapsPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 5 :5,isBionanoDisplay ? 0 :5,isBionanoDisplay ? 0 :5,isBionanoDisplay ? 0 :5));
        gbc=new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.gridy=0;        gbc.gridheight=1;     
        mapsPanel.add(butMapsImgPan,gbc);
        gbc.gridx=1;        gbc.gridheight=1;     
        mapsPanel.add(butMapsPan,gbc);

        

        
    	JPanel leftButPan=new JPanel();
    	leftButPan.setLayout(new GridBagLayout());
    	leftButPan.setBorder(BorderFactory.createEmptyBorder(0,isBionanoDisplay ? 2 :10,0,isBionanoDisplay ? 5 :25));
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
        text.setFont(new Font("Arial",0,isBionanoDisplay ? 7 : 10));
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
        text.setFont(new Font("Arial",0,isBionanoDisplay ? 7 : 10));
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
    	rightButPan.setBorder(BorderFactory.createEmptyBorder(0,isBionanoDisplay ? 5 :25,0,isBionanoDisplay ? 2 :10));
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
        text.setFont(new Font("Arial",0,isBionanoDisplay ? 3 : 5));
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
        text.setFont(new Font("Arial",0,isBionanoDisplay ? 3 : 5));
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
        centralButtonsPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
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
		imagesPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
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
			paramsLeftTX.setLayout(new GridLayout(2,1,0,0));
			JPanel paramsRightTX=new JPanel();
			paramsRightTX.setBorder(BorderFactory.createEmptyBorder(0,0,0,0));
			paramsRightTX.setLayout(new GridLayout(2,1,0,0));
	
			paramsPanelRightTX=new JPanel [2][4];
			for(int tt=0;tt<2;tt++) {
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
		t1t2PlotsPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
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
		t1t2ParamsPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
		t1t2ParamsPanel.setLayout(new GridBagLayout());
		gbc = new GridBagConstraints();
        gbc.gridy=0;        gbc.gridheight=1;     
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.fill = GridBagConstraints.BOTH;
        
        
        JPanel switchModelPanel=new JPanel();
        switchModelPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
        switchModelPanel.setLayout(new GridLayout(5,1,10,0));
        switchModelPanel.add(buttonExp);
        switchModelPanel.add(new JLabel());
        switchModelPanel.add(buttonOff);
        switchModelPanel.add(new JLabel());
        switchModelPanel.add(buttonRice);
        JPanel exportPanel=new JPanel();
        exportPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
        exportPanel.setLayout(new GridLayout(5,1,10,0));
        exportPanel.add(new JLabel());
        exportPanel.add(new JLabel());
        exportPanel.add(buttonExport);
        exportPanel.add(new JLabel());
        exportPanel.add(new JLabel());

        JPanel switchMonoBiPanel=new JPanel();
        switchMonoBiPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
        switchMonoBiPanel.setLayout(new GridLayout(2,1,10,0));
        if(hasT2) {
	        switchMonoBiPanel.add(buttonSwitchMono);
	        switchMonoBiPanel.add(buttonSwitchBi);
        }
        
        gbc.gridx=0;        gbc.gridheight=1;     
		t1t2ParamsPanel.add(switchModelPanel);
        gbc.gridx=1;        gbc.gridheight=1;     
		t1t2ParamsPanel.add(switchMonoBiPanel);
        gbc.gridx=2;        gbc.gridheight=1;     
		t1t2ParamsPanel.add(paramsTXPanel[0],gbc);
        gbc.gridx=3;        gbc.gridheight=1;     
		t1t2ParamsPanel.add(paramsTXPanel[1],gbc);
        gbc.gridx=4;        gbc.gridheight=1;     
		t1t2ParamsPanel.add(exportPanel,gbc);
		
        //Spectrum T1 and T2 panel
		JPanel t1t2SpectrumPanel=new JPanel();
		t1t2SpectrumPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
		t1t2SpectrumPanel.setLayout(new GridBagLayout());
		gbc = new GridBagConstraints();
        gbc.gridy=0;        gbc.gridheight=1;     
        gbc.gridwidth=1;     gbc.gridx=0;
        gbc.fill = GridBagConstraints.BOTH;
        
		t1t2SpectrumPanel.add(plotCan21,gbc);
        gbc.gridx=1;        gbc.gridheight=1;     
		t1t2SpectrumPanel.add(plotCan22,gbc);
		
		
		
		int policeSize=isBionanoDisplay ? 1 : 2;
        //Right panel
		JPanel rightPanel=new JPanel();
		rightPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10,isBionanoDisplay ? 0 :10));
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
		globalPanel.setBorder(BorderFactory.createEmptyBorder(isBionanoDisplay ? 0 :5,isBionanoDisplay ? 0 :5,isBionanoDisplay ? 0 :5,isBionanoDisplay ? 0 :5));
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
		
		buttonExp.addActionListener(this);
		buttonExport.addActionListener(this);
		buttonOff.addActionListener(this);
		buttonRice.addActionListener(this);
		buttonSwitchMono.addActionListener(this);
		buttonSwitchBi.addActionListener(this);
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
        this.setResizable(true);
        pack();
        if(isBionanoDisplay) {
        	this.setSize(1180, 680);
        }
        setVisible(true);
        setTitle("MRI Curve explorer V2");
		repaint();
		System.out.println(this.getSize());
	}

	/**
	 * Start plotting the spectral density and the ROI over HyperMap.
	 */
	public void startPlotsAndRoi(){
		IJ.log("Starting Plots");
		if(!T1T2MixedEstimation) {
			plotT1 = new Plot("T1 curve explorer","Recovery time","Spin-echo magnitude signal");
			plotT1.changeFont(new Font("Helvetica", 0, isBionanoDisplay ? 10 : 14));
			plotCan1=new PlotCanvas(plotT1.getImagePlus());
			plotCan1.setPlot(plotT1);
	
			plotT2 = new Plot("T2 curve explorer","Echo time","Spin-echo magnitude signal");
			plotT2.changeFont(new Font("Helvetica", 0, isBionanoDisplay ? 10 : 14));
			plotCan2=new PlotCanvas(plotT2.getImagePlus());
			plotCan2.setPlot(plotT2);
		}
		else {
			plotT1T2 = new Plot("T1-T2 curve explorer",("Repetition time + "+visualFactorT2+" x Echo time"),"Magnitude signal");
			plotT1T2.changeFont(new Font("Helvetica", 0, isBionanoDisplay ? 10 : 14));
			plotT1T2.setSize(TARGET_WIDTH_FOR_PLOTS_T1T2, TARGET_HEIGHT_FOR_UP_PLOTS_AND_IMGS);
			plotCanT1T2=new PlotCanvas(plotT1T2.getImagePlus());
			plotCanT1T2.setPlot(plotT1T2);
		}
		
		plotT21 = new Plot("T1 timelapse tracker","T1 distribution (PD-weighted)","PD-weighted distribution",Plot.DEFAULT_FLAGS-Plot.X_GRID-Plot.Y_GRID+Plot.X_LOG_TICKS+Plot.X_LOG_NUMBERS);
		plotT22 = new Plot("T2 timelapse tracker","T2 distribution (PD-weighted)","PD-weighted distribution",Plot.DEFAULT_FLAGS-Plot.X_GRID-Plot.Y_GRID+Plot.X_LOG_TICKS+Plot.X_LOG_NUMBERS);
		plotT21.changeFont(new Font("Helvetica", 0, isBionanoDisplay ? 10 : 14));
		plotT22.changeFont(new Font("Helvetica", 0, isBionanoDisplay ? 10 : 14));
		plotT21.setBackgroundColor(sexyForBioInfo ? new Color(255,255,255) : new Color(0,0,0));
		plotT22.setBackgroundColor(sexyForBioInfo ? new Color(255,255,255) : new Color(0,0,0));
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
	
	/**
	 * Actualize the ROI Cursor
	 */
	public void actualizeCursor() {
		Overlay overT1;
		Overlay overT2;
		if(statusRoi!=2 || userRoi==null) {
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
	
	/**
	 * Get cursor size.
	 *
	 * @return the string
	 */
	public String sizeOfCursor() {
		return "medium";
	}
		
	/**
	 * Actualize ranging boundaries for display selection in the spectrum
	 *
	 * @return the selected time (corresponding to click coordinates)
	 */
	public int actualizeRangingBoundaries() {
		double borderLeft=isBionanoDisplay ? 56 : 76;
		double borderRight=isBionanoDisplay ? 14 : 20;
		double borderUp=isBionanoDisplay ? 14 : 18;
		double borderDown=isBionanoDisplay ? 32 : 47;
		double lengPlot=plotCan21.getWidth()-borderRight-borderLeft;
		double latPlot=(plotCan21.getHeight()-yMouseRange-borderDown)/(plotCan21.getHeight()-borderDown-borderUp);
		double xPos=this.xMouseRange-borderLeft;
		double[]vals=spectrumRangingModeT1 ? plotT21.getLimits() : plotT22.getLimits();
		double factMul=vals[1]/vals[0];
		double rangCenter=Math.pow(factMul,xPos*1.0/lengPlot)*vals[0];
		IJ.log("Event on histogram. Cursor position = "+rangCenter+" ms");
		this.rangingBoundaries=new double[] {rangCenter*(1-this.rangingFactor),rangCenter,rangCenter*(1+this.rangingFactor)};
		int newTime=Math.max(0,Math.min(this.nTimepoints-1,    (int)Math.floor( this.nTimepoints*1.0*latPlot )  ));
		System.out.println("Click on "+xMouseRange+" , "+yMouseRange+" -> Time="+newTime+" Range center ="+rangCenter);
		return newTime;
	}
	
	/**
	 * Close plots and roi.
	 */
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
	
	
	
	
	
	
	
	

	
	/**
	 *  Read new data and update estimation and fits--------------------------------------------------------------------------------.
	 */	
	public void computeResultsAgain() {
		actualizeMeanEstimations();
		actualizeDisplayedNumbers();
		Timer tim2=new Timer();
		computeEstimationsForAllPointsMultiThreaded();
		if(!T1T2MixedEstimation) {
			tim2.print("Fit updated on "+this.nPtsCur+" voxels over "+this.nTimepoints+" timepoints."+
					" Total fits : "+(((hasT1 ? 1 : 0 )+(hasT2 ? 2 : 0))*this.nPtsCur*this.nTimepoints)+"  . Total computation time : ");
		}
		else 			tim2.print("Fit updated on "+this.nPtsCur+" voxels over "+this.nTimepoints+" timepoints."+
				" Total fits : "+(2*this.nPtsCur*this.nTimepoints)+"  . Total computation time : ");
		actualizeSpectrumCurves();
		identifyRangedData();
	}
	
	
	/**
	 * Actualize spectrum curves. Chain-called when clicking
	 */
	public void actualizeSpectrumCurves() {
		if(T1T2MixedEstimation)actualizeSpectrumCurvesT1T2();
		else {
			if(hasT1)actualizeSpectrumCurvesT1();
			if(hasT2)actualizeSpectrumCurvesT2();
		}
	}

	/**
	 * Actualize mri observations based on data of the ROI or the clicked area
	 */
	public void actualizeMriObservationsBasedOnData() {
		System.out.println("ACTUALISATION");
		System.out.println("Status="+statusRoi);
		//Update for all time, for all voxels selected
		if(statusRoi!=2 || (userRoi==null && imgRoi==null)) {
			//if(!T1T2MixedEstimation) 
			//this.dataTimelapseFull=hyperMap.getFullMRISignalAroundThisVoxel((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick);
			//else 
			this.dataTimelapseFull=hyperMap.getFullMRISignalAroundThisVoxelT1T2((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick,this.sigmaSmoothingBeforeEstimation);
			this.pointsCoords=hyperMap.getCoordinatesAroundThisVoxel((int)xCor,(int)yCor,zCor,this.crossWidth,this.crossThick);
		}
		else {
			//if(!T1T2MixedEstimation) 
			//this.dataTimelapseFull=hyperMap.getFullMRISignalInTheseCoordinates((int)xCor,(int)yCor,zCor,this.pointsCoords);
			//else
			System.out.println("DOING THIS BECAUSE ROI");
			IJ.log("Actualizing data from Roi, at coordinates Z="+zCor+" , T="+tCor);
			if(imgRoi!=null) {this.pointsCoords=VitimageUtils.getMaskAsCoords(this.imgRoi);System.out.println("DOING THIS BECAUSE IMG ROI");}
			else this.pointsCoords=VitimageUtils.getRoiAsCoords(this.userRoi);
			this.dataTimelapseFull=hyperMap.getFullMRISignalInTheseCoordinatesT1T2((int)xCor,(int)yCor,zCor,this.pointsCoords,this.sigmaSmoothingBeforeEstimation);
			System.out.println("Reading "+this.pointsCoords.length+" points from Roi or mask");
		}
		this.nPtsCur=this.dataTimelapseFull[0].length;

		//For each time, compute the mean over voxels
		double[]vals;
		for(int t=0;t<this.nTimepoints;t++) {			
			if(!T1T2MixedEstimation) {
				if(hasT1)this.dataTimelapseT1[t]=new double[this.dataTimelapseFull[t][0][0].length];
				if(hasT1)this.dataTimelapseT1Sigmas[t]=new double[this.dataTimelapseFull[t][0][0].length];
				if(hasT2)this.dataTimelapseT2[t]=new double[this.dataTimelapseFull[t][0][0].length];
				if(hasT2)this.dataTimelapseT2Sigmas[t]=new double[this.dataTimelapseFull[t][0][0].length];
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
					if(hasT1)this.dataTimelapseT1[t][ec]=stats[0];
					if(hasT1)this.dataTimelapseT1Sigmas[t][ec]=stats[1];						
				}
				else {
					this.dataTimelapseT1T2[t][ec]=stats[0];
					this.dataTimelapseT1T2Sigmas[t][ec]=stats[1];											
				}
			}
			if(!T1T2MixedEstimation) {
				for(int ec=0;ec<this.dataTimelapseFull[t][0][0].length;ec++) {
					vals=new double[this.nPtsCur];
					for(int vo=0;vo<this.nPtsCur;vo++) {
						vals[vo]=this.dataTimelapseFull[t][vo][0][ec];
					}
					double []stats=VitimageUtils.statistics1D(vals);
					if(hasT2)this.dataTimelapseT2[t][ec]=stats[0];
					if(hasT2)this.dataTimelapseT2Sigmas[t][ec]=stats[1];						
				}
			}
		}
	}
			


	/**
	 * Actualize selected echo, when scrolling the mouse
	 *
	 * @param plotType the plot type
	 * @return the int
	 */
	public int actualizeSelectedEcho(int plotType) {
		double borderLeft=isBionanoDisplay ? 56 : 76;
		double borderRight=isBionanoDisplay ? 14 : 20;
		double borderUp=isBionanoDisplay ? 15 : 18;
		double borderBottom=isBionanoDisplay ? 35 : 44;
		double lengPlot=(plotType==WIN_PLOT1 ? plotCan1.getWidth() : (plotType==WIN_PLOT2 ? plotCan2.getWidth() : plotCanT1T2.getWidth()))-borderRight-borderLeft;
		double latPlot=(plotType==WIN_PLOT1 ? plotCan1.getHeight() : (plotType==WIN_PLOT2 ? plotCan2.getHeight() : plotCanT1T2.getHeight()))-borderUp-borderBottom;
		double xPos=(this.xMouseRangeCurve-borderLeft)/lengPlot;
		double yPos=1-(this.yMouseRangeCurve-borderUp)/latPlot;
		double[]vals=(plotType==WIN_PLOT1 ? plotT1.getLimits() : plotType==WIN_PLOT2 ? plotT2.getLimits() : plotT1T2.getLimits());
		double xPosPlot=vals[1]*xPos;
		double yPosPlot=vals[3]*yPos;
		System.out.println("Point clicked in plot : "+this.xMouseRangeCurve+" , "+this.yMouseRangeCurve+" -> "+xPosPlot+","+yPosPlot);
		int newTime=1;double minDist=1E20;
		if(plotType==WIN_PLOT1 && (hasT1 || T1T2MixedEstimation)) {
			for(int t=0;t<this.t1Times[tCor][zCor].length;t++) {
				if(Math.abs(this.t1Times[tCor][zCor][t]-xPosPlot)<minDist) {
					minDist=Math.abs(this.t1Times[tCor][zCor][t]-xPosPlot);
					newTime=t;
				}
			}
		}
		else if(plotType==WIN_PLOT2 && (hasT2 || T1T2MixedEstimation)){
			for(int t=0;t<this.t2Times[tCor][zCor].length;t++) {
				if(Math.abs(this.t2Times[tCor][zCor][t]-xPosPlot)<minDist) {
					minDist=Math.abs(this.t2Times[tCor][zCor][t]-xPosPlot);
					newTime=t;
				}
			}
		}
		else if(plotType==WIN_PLOT1T2){
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
		
	/**
	 * Actualize mean when data changed (ROI or cursor clicked)
	 */
	public void actualizeMeanEstimations() {
		Object[] obj;
		Timer tt=new Timer();
		for(int tim=0;tim< this.nTimepoints;tim++) {
			if(!T1T2MixedEstimation) {
				paramsTimelapseT1[tim]=new double[20];
				paramsTimelapseT2[tim]=new double[20];
				//Estimer T1 Monocomp
				if(hasT1) {
					obj=fitAndEvaluateT1(t1TrTimes[tim][zCor],t1TeTimes[tim][zCor],this.timesT1Cute[tim][zCor],dataTimelapseT1[tim],hyperMap.tabSigmasT1Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,getFitType(MRUtils.T1_MONO_RICE),this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
					for(int prm=0;prm<((double[]) obj[0]).length;prm++) {
			 			paramsTimelapseT1[tim][prm]=((double[]) obj[0])[prm];
			 		}
			 		tabFittenT1Mono[tim]=(double[]) obj[2];
			 		tabFittenT1MonoCute[tim]=(double[][]) obj[3];
			 		khi2T1Mono[tim]=(double) obj[4];
			 		pValT1Mono[tim]=(double) obj[5];
			 		jitterT1Mono[tim]=(double) obj[6];
				}
		
		 		//Estimer T2 Monocomp
				if(hasT2) {
					System.out.println(TransformUtils.stringVectorN(t2TrTimes[tim][zCor],"1"));
					System.out.println(TransformUtils.stringVectorN(t2TeTimes[tim][zCor],"2"));
					obj=fitAndEvaluateT2(t2TrTimes[tim][zCor],t2TeTimes[tim][zCor],this.timesT2Cute[tim][zCor],dataTimelapseT2[tim],hyperMap.tabSigmasT2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,getFitType(MRUtils.T2_MONO_RICE),this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
			 		for(int prm=0;prm<((double[]) obj[0]).length;prm++) {
			 			paramsTimelapseT2[tim][prm]=((double[]) obj[0])[prm];
			 		}
			 		tabFittenT2Mono[tim]=(double[]) obj[2];
			 		tabFittenT2MonoCute[tim]=(double[][]) obj[3];
			 		khi2T2Mono[tim]=(double) obj[4];
			 		pValT2Mono[tim]=(double) obj[5];
			 		jitterT2Mono[tim]=(double) obj[6];
				
			 		//Estimer T2 Bicomp
					obj=fitAndEvaluateT2(t2TrTimes[tim][zCor],t2TeTimes[tim][zCor],this.timesT2Cute[tim][zCor],dataTimelapseT2[tim],hyperMap.tabSigmasT2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,getFitType(MRUtils.T2_MULTI_RICE),this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
					for(int prm=0;prm<((double[]) obj[0]).length;prm++) {
				 			paramsTimelapseT2[tim][6+prm]=((double[]) obj[0])[prm];
			 		}
			 		tabFittenT2Bicomp[tim]=(double[]) obj[2];
			 		tabFittenT2BicompCute[tim]=(double[][]) obj[3];
			 		khi2T2Bicomp[tim]=(double) obj[4];
			 		pValT2Bicomp[tim]=(double) obj[5];
			 		jitterT2Bicomp[tim]=(double) obj[6];
				}
			}
			else {
		 		//Estimer T1T2 Monocomp
				paramsTimelapseT1T2[tim]=new double[20];
		 		
				obj=fitAndEvaluateT1T2(t1t2TrTimes[tim][zCor],t1t2TeTimes[tim][zCor],this.timesT1T2Cute[tim][zCor],dataTimelapseT1T2[tim],hyperMap.tabSigmasT1T2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,getFitType(MRUtils.T1T2_MONO_RICE),this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
				for(int prm=0;prm<((double[]) obj[0]).length;prm++) {
			 		paramsTimelapseT1T2[tim][prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT1T2Mono[tim]=(double[]) obj[2];
		 		tabFittenT1T2MonoCute[tim]=(double[][]) obj[3];
		 		khi2T1T2Mono[tim]=(double) obj[4];
		 		pValT1T2Mono[tim]=(double) obj[5];
		 		jitterT1T2Mono[tim]=(double) obj[6];
			
		 		//Estimer T1T2 Multicomp
		 		obj=fitAndEvaluateT1T2(t1t2TrTimes[tim][zCor],t1t2TeTimes[tim][zCor],this.timesT1T2Cute[tim][zCor],dataTimelapseT1T2[tim],hyperMap.tabSigmasT1T2Seq[tim][zCor]+deltaSigma,this.fitAlgorithm,getFitType(MRUtils.T1T2_MULTI_RICE),this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1,true);
		 		for(int prm=0;prm<((double[]) obj[0]).length;prm++) {
				 		paramsTimelapseT1T2[tim][6+prm]=((double[]) obj[0])[prm];
		 		}
		 		tabFittenT1T2Bicomp[tim]=(double[]) obj[2];
		 		tabFittenT1T2BicompCute[tim]=(double[][]) obj[3];
		 		khi2T1T2Bicomp[tim]=(double) obj[4];
		 		pValT1T2Bicomp[tim]=(double) obj[5];
		 		jitterT1T2Bicomp[tim]=(double) obj[6];
			}
		}		
		System.out.println("||");
		System.out.println("||");
	}
	
	

	/**
	 * Compute estimations for all points multi threaded.
	 */
	public void computeEstimationsForAllPointsMultiThreaded() {
		if(T1T2MixedEstimation) {computeEstimationsForAllPointsMultiThreadedT1T2();return;}
		else {
			if(hasT1) {
				System.out.println("Compute T1 multi thread");
				computeEstimationsForAllPointsMultiThreadedT1();
			}
			
			if(hasT2) {
				System.out.println("Compute T2 multi thread");
				computeEstimationsForAllPointsMultiThreadedT2();
			}
			return;
		}
	}

	/**
	 * Multi-threaded estimations of curve parameters using all data points, estimating a T1T2 joint fit
	 */
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
					dataParams[nt][nVo][nTime][3]=new double[] {hyperMap.tabSigmasT1T2Seq[nTime][zCor]};
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
							if(nData>10 && (numBl%(nData/10)==0))System.out.println(" "+VitimageUtils.dou((numBl*100.0)/nData)+"%");
							for (int numTime=0;numTime<nTimes;numTime++) {
								double[][]tempParams=new double[1][17];
						 		//Estimer T1T2 Monocomp
								Object[] obj=fitAndEvaluateT1T2(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,getFitType(MRUtils.T1T2_MONO_RICE),1,false);
								for(int n=0;n<((double[])obj[0]).length;n++)tempParams[0][n]=((double[])obj[0])[n];
								tempParams[0][15]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1T2Mono=(double) obj[4];
								double pT1T2Mono=(double) obj[5];
						
						 		//Estimer T1MonoT2Bi
								obj=fitAndEvaluateT1T2(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,getFitType(MRUtils.T1T2_MULTI_RICE),1,false);
								for(int n=0;n<((double[])obj[0]).length;n++)tempParams[0][n+6]=((double[])obj[0])[n];
								tempParams[0][16]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1T2Bi=(double) obj[4];
								double pT1T2Bi=(double) obj[5];

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
					//if(nTime==tCor)System.out.println("Effectivement à vo="+nVo+" : "+pointsEstimations[nTime][listVoxOnThreads[nt][nVo]][0][8]+
					//		"\n  Estimated using "+TransformUtils.stringVectorN(dataParams[nt][nVo][nTime][2],""));
					
				}
			}
		}
	}
		
	/**
	 * Multi-threaded estimations of curve parameters using all data points, estimating a T1 recovery fit
	 */
	public void computeEstimationsForAllPointsMultiThreadedT1() {
		//Prepare input data and output places for threads
		int nThreads=VitimageUtils.getNbCores();
		Thread[]threads=VitimageUtils.newThreadArray(nThreads);
		int[][]listVoxOnThreads=VitimageUtils.listForThreads(this.nPtsCur, nThreads);
		AtomicInteger atomNumThread=new AtomicInteger(0);
		AtomicInteger curProcessedBlock=new AtomicInteger(0);
		final double[][][][][]dataParams=new double[nThreads][][][][];
		final double[][][][][]estimatedParams=new double[nThreads][][][][];
		final int finalFitAlgorithm=this.fitAlgorithm;
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
					dataParams[nt][nVo][nTime][0]=t1TrTimes[nTime][zCor];
					dataParams[nt][nVo][nTime][1]=t1TeTimes[nTime][zCor];//TODO : danger at the slice 0 and 1 if crossThick > 0, estimation will be abusive. Cross thick have been prevented acting, thus
					dataParams[nt][nVo][nTime][2]=dataTimelapseFull[nTime][listVoxOnThreads[nt][nVo]][0];
					dataParams[nt][nVo][nTime][3]=new double[] {hyperMap.tabSigmasT1Seq[nTime][zCor]};
				}
			}
		}

	
		
		for (int ithread = 0; ithread < nThreads; ithread++) {  
			threads[ithread] = new Thread() {  { setPriority(Thread.NORM_PRIORITY); }  
				public void run() {  
					try {
						int numThread=atomNumThread.getAndIncrement();
						for (int numVo=0;numVo<dataParams[numThread].length;numVo++) {
							int numBl=curProcessedBlock.getAndIncrement();
							if(nData>10 && (numBl%(nData/10)==0))IJ.log(" "+VitimageUtils.dou((numBl*100.0)/nData)+"%");
							if(nData>10 && (numBl%(nData/10)==0))System.out.println(" "+VitimageUtils.dou((numBl*100.0)/nData)+"%");
							for (int numTime=0;numTime<nTimes;numTime++) {
								double[][]tempParams=new double[1][17];
						 		//Estimer T1 Monocomp
								Object[] obj=fitAndEvaluateT1(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,getFitType(MRUtils.T1_MONO_RICE),1,false);
								for(int n=0;n<((double[])obj[0]).length;n++)tempParams[0][n]=((double[])obj[0])[n];
								tempParams[0][15]=((double) obj[6])>=9999 ? 0 : 1;
								double khiT1T2Mono=(double) obj[4];
								double pT1T2Mono=(double) obj[5];
						
								estimatedParams[numThread][numVo][numTime]=tempParams;
							}							
						}
					} catch(Exception ie) {ie.printStackTrace();}
				} 
			};  		
		}				
		VitimageUtils.startAndJoin(threads);
		
		
		//Gather results
		pointsEstimations=new double[this.nTimepoints][this.nPtsCur][2][17];
		for(int nt=0;nt<nThreads;nt++) {
			for(int nVo=0;nVo<listVoxOnThreads[nt].length;nVo++) {
				for(int nTime=0;nTime<this.nTimepoints;nTime++) {
					pointsEstimations[nTime][listVoxOnThreads[nt][nVo]][0]=estimatedParams[nt][nVo][nTime][0];
					//if(nTime==tCor)System.out.println("Effectivement à vo="+nVo+" : "+pointsEstimations[nTime][listVoxOnThreads[nt][nVo]][0][8]+
					//		"\n  Estimated using "+TransformUtils.stringVectorN(dataParams[nt][nVo][nTime][2],""));
					
				}
			}
		}
	}
	
	/**
	 * Multi-threaded estimations of curve parameters using all data points, estimating a T2 recovery fit
	 */
	public void computeEstimationsForAllPointsMultiThreadedT2() {
		Timer t=new Timer();
		//Prepare input data and output places for threads
		int nThreads=VitimageUtils.getNbCores();
		Thread[]threads=VitimageUtils.newThreadArray(nThreads);
		int[][]listVoxOnThreads=VitimageUtils.listForThreads(this.nPtsCur, nThreads);
		AtomicInteger atomNumThread=new AtomicInteger(0);
		AtomicInteger curProcessedBlock=new AtomicInteger(0);
		final double[][][][][]dataParams=new double[nThreads][][][][];
		final double[][][][][]estimatedParams=new double[nThreads][][][][];
		final int finalFitAlgorithm=this.fitAlgorithm;
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
					dataParams[nt][nVo][nTime][0]=t2TrTimes[nTime][zCor];
					dataParams[nt][nVo][nTime][1]=t2TeTimes[nTime][zCor];//TODO : danger at the slice 0 and 1 if crossThick > 0, estimation will be abusive. Cross thick have been prevented acting, thus
					dataParams[nt][nVo][nTime][2]=dataTimelapseFull[nTime][listVoxOnThreads[nt][nVo]][0];
					dataParams[nt][nVo][nTime][3]=new double[] {hyperMap.tabSigmasT2Seq[nTime][zCor]};
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
							if(nData>10 && (numBl%(nData/10)==0))System.out.println(" "+VitimageUtils.dou((numBl*100.0)/nData)+"%");
							for (int numTime=0;numTime<nTimes;numTime++) {
								double[][]tempParams=new double[2][17];
								Object[] obj=null;

								if(true) {
									//Estimer T2 Monocomp
						 			obj=fitAndEvaluateT2(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,getFitType(MRUtils.T2_MONO_RICE),1,false);
									for(int n=0;n<((double[])obj[0]).length;n++)tempParams[0][n]=((double[])obj[0])[n];
									tempParams[0][15]=((double) obj[6])>=9999 ? 0 : 1;
									double khiT2Mono=(double) obj[4];
									double pT2Mono=(double) obj[5];
								}
	
								if(true) {
							 		//Estimer T1MonoT2Bi
									obj=fitAndEvaluateT2(dataParams[numThread][numVo][numTime][0],dataParams[numThread][numVo][numTime][1],null,dataParams[numThread][numVo][numTime][2],dataParams[numThread][numVo][numTime][3][0],finalFitAlgorithm,getFitType(MRUtils.T2_MULTI_RICE),1,false);
									for(int n=0;n<((double[])obj[0]).length;n++)tempParams[0][n+6]=((double[])obj[0])[n];
									tempParams[0][16]=((double) obj[6])>=9999 ? 0 : 1;
									double khiT1T2Bi=(double) obj[4];
									double pT1T2Bi=(double) obj[5];
								}
								//Select best
								estimatedParams[numThread][numVo][numTime]=tempParams;
							}							
						}
					} catch(Exception ie) {ie.printStackTrace();}
				} 
			};  		
		}				
		VitimageUtils.startAndJoin(threads);

		//Gather results
		if(!hasT1)pointsEstimations=new double[this.nTimepoints][this.nPtsCur][2][17];
		for(int nt=0;nt<nThreads;nt++) {
			for(int nVo=0;nVo<listVoxOnThreads[nt].length;nVo++) {
				for(int nTime=0;nTime<this.nTimepoints;nTime++) {
					pointsEstimations[nTime][listVoxOnThreads[nt][nVo]][1]=estimatedParams[nt][nVo][nTime][0];
					//if(nTime==tCor)System.out.println("Effectivement à vo="+nVo+" : "+pointsEstimations[nTime][listVoxOnThreads[nt][nVo]][0][8]+
					//		"\n  Estimated using "+TransformUtils.stringVectorN(dataParams[nt][nVo][nTime][2],""));
					
				}
			}
		}
	}

	
	
	
	
	
	
	/**
	 * Fit selected equations and evaluate their accuracy by updating displayed Khi2 values.
	 *
	 * @param tabTimesTr the recovery times
	 * @param tabTimesTe the echo times
	 * @param tabCuteTimesPerCurve tab of interpolated time points for esthetic curve display
	 * @param tabData the actual magnitude data
	 * @param sigmaRice the sigma rice
	 * @param fitAlgorithm the fit algorithm
	 * @param fitCurveType the fit curve type (among fijirelax.mrialgo.MRDataType)
	 * @param nbPts the number of points involved in the computation
	 * @param niceCurveComputation boolean flag to set for getting smooth curve computation
	 * @return an Object[] containing {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accuracy1, accuracy2,jitter};
	 */
	public Object[] fitAndEvaluateT1T2(double[]tabTimesTr,double[]tabTimesTe,double[][][]tabCuteTimesPerCurve,double[]tabData,double sigmaRice,int fitAlgorithm,int fitCurveType,int nbPts,boolean niceCurveComputation) {
		int nParams= MRUtils.getNparams(fitCurveType);//(type==0) ? 3 : (type==1 ? 5 : 6);
		double []estimatedParams=new double[nParams];
		double []estimatedSigmas=new double[nParams];	

		double []estimation=(double[]) MRUtils.makeFit(tabTimesTr, tabTimesTe, tabData, fitCurveType, MRUtils.SIMPLEX, MRUtils.N_ITER_T1T2, sigmaRice,false)[0];
		for(int i=0;i<nParams;i++) {estimatedParams[i]=estimation[i] ;  estimatedSigmas[i]=estimation[i];}
		double []tabFitten=MRUtils.fittenRelaxationCurve(tabTimesTr,tabTimesTe,estimatedParams,sigmaRice,fitCurveType);
		double [][]tabFittenCute=null;
		if(niceCurveComputation) {
			tabFittenCute=new double[tabCuteTimesPerCurve.length][];
			for(int fi=0;fi<tabCuteTimesPerCurve.length;fi++) {
				tabFittenCute[fi]=MRUtils.fittenRelaxationCurve(tabCuteTimesPerCurve[fi][0],tabCuteTimesPerCurve[fi][1],estimatedParams,sigmaRice,fitCurveType);
			}
		}
		double[]accs=MRUtils.fittingAccuracies(tabData,tabTimesTr,tabTimesTe,sigmaRice,estimatedParams,fitCurveType,false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
		if(((fitCurveType%10)==MRUtils.T2_MULTI) && (estimatedParams[2]>estimatedParams[3])) {
			double tmp=estimatedParams[0];estimatedParams[0]=estimatedParams[1];estimatedParams[1]=tmp;
			tmp=estimatedParams[2];estimatedParams[2]=estimatedParams[3];estimatedParams[3]=tmp;
		}		
		
		if(((fitCurveType%10)==MRUtils.T1T2_MULTI) && (estimatedParams[3]>estimatedParams[4])) {
			double tmp=estimatedParams[0];estimatedParams[0]=estimatedParams[1];estimatedParams[1]=tmp;
			tmp=estimatedParams[3];estimatedParams[3]=estimatedParams[4];estimatedParams[4]=tmp;
		}		
		
		double jitter=0;
		for(int i=0;i<estimatedParams.length;i++) {if(estimatedParams[i]<=0) {jitter=9999;accs[1]=100;}}
		if(VitimageUtils.max(tabData)<MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*sigmaRice)  {//T2 curve with the first point being very low -> no estimation possible there
			jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 computation cause first point is too low : "+tabData[0]+" , "+statsRice[0]+" , "+statsRice[1]);
		}
		
		return new Object[] {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accs[0],accs[1],jitter};
	}
		
	/**
	 * Fit selected equations and evaluate their accuracy by updating displayed Khi2 values.
	 *
	 * @param tabTimesTr the recovery times
	 * @param tabTimesTe the echo times
	 * @param tabCuteTimesPerCurve tab of interpolated time points for esthetic curve display
	 * @param tabData the actual magnitude data
	 * @param sigmaRice the sigma rice
	 * @param fitAlgorithm the fit algorithm
	 * @param fitCurveType the fit curve type (among fijirelax.mrialgo.MRDataType)
	 * @param nbPts the number of points involved in the computation
	 * @param niceCurveComputation boolean flag to set for getting smooth curve computation
	 * @return an Object[] containing {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accuracy1, accuracy2,jitter};
	 */
	public Object[] fitAndEvaluateT1(double[]tabTimesTr,double[]tabTimesTe,double[][][]tabCuteTimesPerCurve,double[]tabData,double sigmaRice,int fitAlgorithm,int fitCurveType,int nbPts,boolean niceCurveComputation) {
		int nParams= MRUtils.getNparams(fitCurveType);
		double []estimatedParams=new double[nParams];
		double []estimatedSigmas=new double[nParams];	

		double []estimation=(double[]) MRUtils.makeFit(tabTimesTr, tabTimesTe, tabData, fitCurveType, MRUtils.SIMPLEX, MRUtils.N_ITER_T1T2, sigmaRice,false)[0];
		for(int i=0;i<nParams;i++) {estimatedParams[i]=estimation[i] ;  estimatedSigmas[i]=estimation[i];}
		double []tabFitten=MRUtils.fittenRelaxationCurve(tabTimesTr,tabTimesTe,estimatedParams,sigmaRice,fitCurveType);
		double [][]tabFittenCute=null;
		if(niceCurveComputation) {
			tabFittenCute=new double[tabCuteTimesPerCurve.length][];
			for(int fi=0;fi<tabCuteTimesPerCurve.length;fi++) {
				tabFittenCute[fi]=MRUtils.fittenRelaxationCurve(tabCuteTimesPerCurve[fi][0],tabCuteTimesPerCurve[fi][1],estimatedParams,sigmaRice,fitCurveType);
			}
		}
		double[]accs=MRUtils.fittingAccuracies(tabData,tabTimesTr,tabTimesTe,sigmaRice,estimatedParams,fitCurveType,false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
		
		double jitter=0;
		for(int i=0;i<estimatedParams.length;i++) {if(estimatedParams[i]<=0) {jitter=9999;accs[1]=100;}}
		if(VitimageUtils.max(tabData)<MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*sigmaRice)  {//T2 curve with the first point being very low -> no estimation possible there
			jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 computation cause first point is too low : "+tabData[0]+" , "+statsRice[0]+" , "+statsRice[1]);
		}
		return new Object[] {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accs[0],accs[1],jitter};
	}
	
	/**
	 * Fit selected equations and evaluate their accuracy by updating displayed Khi2 values.
	 *
	 * @param tabTimesTr the recovery times
	 * @param tabTimesTe the echo times
	 * @param tabCuteTimesPerCurve tab of interpolated time points for esthetic curve display
	 * @param tabData the actual magnitude data
	 * @param sigmaRice the sigma rice
	 * @param fitAlgorithm the fit algorithm
	 * @param fitCurveType the fit curve type (among fijirelax.mrialgo.MRDataType)
	 * @param nbPts the number of points involved in the computation
	 * @param niceCurveComputation boolean flag to set for getting smooth curve computation
	 * @return an Object[] containing {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accuracy1, accuracy2,jitter};
	 */
	public Object[] fitAndEvaluateT2(double[]tabTimesTr,double[]tabTimesTe,double[][][]tabCuteTimesPerCurve,double[]tabData,double sigmaRice,int fitAlgorithm,int fitCurveType,int nbPts,boolean niceCurveComputation) {
		int nParams= MRUtils.getNparams(fitCurveType);
		double []estimatedParams=new double[nParams];
		double []estimatedSigmas=new double[nParams];	

		double []estimation=(double[]) MRUtils.makeFit(tabTimesTr, tabTimesTe, tabData, fitCurveType, MRUtils.SIMPLEX, MRUtils.N_ITER_T1T2, sigmaRice,false)[0];
		for(int i=0;i<nParams;i++) {estimatedParams[i]=estimation[i] ;  estimatedSigmas[i]=estimation[i];}
		double []tabFitten=MRUtils.fittenRelaxationCurve(tabTimesTr,tabTimesTe,estimatedParams,sigmaRice,fitCurveType);
		double [][]tabFittenCute=null;
		if(niceCurveComputation) {
			tabFittenCute=new double[tabCuteTimesPerCurve.length][];
			for(int fi=0;fi<tabCuteTimesPerCurve.length;fi++) {
				tabFittenCute[fi]=MRUtils.fittenRelaxationCurve(tabCuteTimesPerCurve[fi][0],tabCuteTimesPerCurve[fi][1],estimatedParams,sigmaRice,fitCurveType);
			}
		}
		double[]accs=MRUtils.fittingAccuracies(tabData,tabTimesTr,tabTimesTe,sigmaRice,estimatedParams,fitCurveType,false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);
		if(((fitCurveType%10)==MRUtils.T2_MULTI) && (estimatedParams[2]>estimatedParams[3])) {
			double tmp=estimatedParams[0];estimatedParams[0]=estimatedParams[1];estimatedParams[1]=tmp;
			tmp=estimatedParams[2];estimatedParams[2]=estimatedParams[3];estimatedParams[3]=tmp;
		}		
				
		double jitter=0;
		for(int i=0;i<estimatedParams.length;i++) {if(estimatedParams[i]<=0) {jitter=9999;accs[1]=100;}}
		if(VitimageUtils.max(tabData)<MRUtils.THRESHOLD_RATIO_BETWEEN_MAX_ECHO_AND_SIGMA_RICE*sigmaRice)  {//T2 curve with the first point being very low -> no estimation possible there
			jitter=9999;accs[1]=100;//System.out.println("Jitter set in T2 computation cause first point is too low : "+tabData[0]+" , "+statsRice[0]+" , "+statsRice[1]);
		}
		return new Object[] {estimatedParams,estimatedSigmas, tabFitten,tabFittenCute,accs[0],accs[1],jitter};
	}

	
	
	
	

	
	/**
	 * Gets the fit type as an int:(fitType %10)+10*this.noiseHandling)
	 *
	 * @param fitType the fit type
	 * @return the fit type
	 */
	public int getFitType(int fitType) {
		return ((fitType %10)+10*this.noiseHandling);
	}
	
	
	
	/**
	 * Identify the ranged data in the spectrum, depending on the actual defined interval (in blue in the plot)
	 */
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

	/**
	 * Identify the ranged data in the spectrum, depending on the actual defined interval (in blue in the plot), when using cross-fitting T1T2
	 */
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
	

		
/**
 * Compute spectrum curve T1 T2, used when cross-fitting
 *
 * @param time the selected timepoint
 * @return the double[][] containing spectral data of T1 and T2 relaxations
 */
public double[][] computeSpectrumCurveT1T2(int time) {
		
		final double[]binInf=new double[numberBins];
		final double[]binSup=new double[numberBins];
		binMedT1=new double[numberBins];
		double[]histoM0T1=new double[numberBins];
		double[]histoM0T2=new double[numberBins];
		double [][]output=new double[3][numberBins];
		double minBin=t0T2;
		double maxBin=t1T1;
		double multFromMinToMax=maxBin/minBin;
		for(int i=0;i<numberBins;i++) {
			binInf[i]=minBin*Math.pow(multFromMinToMax, i*1.0/numberBins);
			binSup[i]=minBin*Math.pow(multFromMinToMax, (i+1)*1.0/numberBins);
			binMedT1[i]=binSup[i]*0.5+binInf[i]*0.5;
		}
		for(int p=0;p<this.nPtsCur;p++) {
		}
		//Bin results
		for(int p=0;p<this.nPtsCur;p++) {

			if(switchT1T2==0) {//Fit T1T2 mono is ok
				//T1 binning
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][1]<binSup[bin]) && (pointsEstimations[time][p][0][1]>=binInf[bin])  && (pointsEstimations[time][p][0][15]>=1)) {
					histoM0T1[bin]+=pointsEstimations[time][p][0][0];
				}
				//T2 binning
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][2]<binSup[bin]) && (pointsEstimations[time][p][0][2]>=binInf[bin])  && (pointsEstimations[time][p][0][15]>=1)) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][0];
				}
			}
			else {//Fit T1T2 bicomp is ok
				//T1 binning
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][8]<binSup[bin]) && (pointsEstimations[time][p][0][8]>=binInf[bin])   && (pointsEstimations[time][p][0][16]>=1)) {
					histoM0T1[bin]+=(pointsEstimations[time][p][0][6]+pointsEstimations[time][p][0][7]);
				}
				//T2 binning
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][9]<binSup[bin]) && (pointsEstimations[time][p][0][9]>=binInf[bin])  && (pointsEstimations[time][p][0][16]>=1)) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][6];
				}
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][10]<binSup[bin]) && (pointsEstimations[time][p][0][10]>=binInf[bin])  && (pointsEstimations[time][p][0][16]>=1)) {
					histoM0T2[bin]+=pointsEstimations[time][p][0][7];
				}
			}
		}
				
		//Smooth this histogram with a factor to be defined, maybe depending on the estimation error on parameters
		output[0]=smoothHisto(histoM0T1,gaussianSpectrum);
		output[1]=smoothHisto(histoM0T2,gaussianSpectrum);
		output[2]=binMedT1;
		return output;
	}


/**
 * Compute spectrum curve T1 .
 *
 * @param time the selected timepoint
 * @return the double[][] containing spectral data of T1 and T2 relaxations
 */
	public double[][] computeSpectrumCurveT1(int time) {
		final double[]binInf=new double[numberBins];
		final double[]binSup=new double[numberBins];
		final double[]binMed=new double[numberBins];
		double[]histoM0T1=new double[numberBins];
		double[]histoM0T2=new double[numberBins];
		double [][]output=new double[3][numberBins];
		double minBin=t0T1;
		double maxBin=t1T1;
		double multFromMinToMax=maxBin/minBin;
		for(int i=0;i<numberBins;i++) {
			binInf[i]=minBin*Math.pow(multFromMinToMax, i*1.0/numberBins);
			binSup[i]=minBin*Math.pow(multFromMinToMax, (i+1)*1.0/numberBins);
			binMed[i]=binSup[i]*0.5+binInf[i]*0.5;
		}
		for(int p=0;p<this.nPtsCur;p++) {
		}
		//Bin results
		for(int p=0;p<this.nPtsCur;p++) {
			for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][0][1]<binSup[bin]) && (pointsEstimations[time][p][0][1]>=binInf[bin]) && (pointsEstimations[time][p][0][15]>=1)) {
				histoM0T1[bin]+=pointsEstimations[time][p][0][0];
			}
		}
	
		//Smooth this histogram with a factor to be defined, maybe depending on the estimation error on parameters
		output[0]=smoothHisto(histoM0T1,gaussianSpectrum);
		output[1]=binMed;
		return output;
	}
	

	/**
	 * Compute spectrum curve T1 .
	 *
	 * @param time the selected timepoint
	 * @return the double[][] containing spectral data of T1 and T2 relaxations
	 */
	public double[][] computeSpectrumCurveT2(int time) {		
		final double[]binInf=new double[numberBins];
		final double[]binSup=new double[numberBins];
		final double[]binMed=new double[numberBins];
		double[]histoM0T1=new double[numberBins];
		double[]histoM0T2=new double[numberBins];
		double [][]output=new double[3][numberBins];
		double minBin=t0T2;
		double maxBin=t1T2;
		double multFromMinToMax=maxBin/minBin;
		for(int i=0;i<numberBins;i++) {
			binInf[i]=minBin*Math.pow(multFromMinToMax, i*1.0/numberBins);
			binSup[i]=minBin*Math.pow(multFromMinToMax, (i+1)*1.0/numberBins);
			binMed[i]=binSup[i]*0.5+binInf[i]*0.5;
		}
		for(int p=0;p<this.nPtsCur;p++) {
		}
		//Bin results
		for(int p=0;p<this.nPtsCur;p++) {
			if(switchT1T2==0) {//Fit T2 mono is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][1][1]<binSup[bin]) && (pointsEstimations[time][p][1][1]>=binInf[bin]) && (pointsEstimations[time][p][1][15]>=1) ) {
					histoM0T2[bin]+=pointsEstimations[time][p][1][0];
				}
			}
			else {//Fit T1T2 bicomp is ok
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][1][8]<binSup[bin]) && (pointsEstimations[time][p][1][8]>=binInf[bin])  && (pointsEstimations[time][p][1][16]>=1)) {
					histoM0T2[bin]+=pointsEstimations[time][p][1][6];
				}
				for(int bin=0;bin<numberBins;bin++)if( (pointsEstimations[time][p][1][9]<binSup[bin]) && (pointsEstimations[time][p][1][9]>=binInf[bin])  && (pointsEstimations[time][p][1][16]>=1)) {
					histoM0T2[bin]+=pointsEstimations[time][p][1][7];
				}
			}
		}
		//Smooth this histogram with a factor to be defined, maybe depending on the estimation error on parameters
		output[0]=smoothHisto(histoM0T2,gaussianSpectrum);
		output[1]=binMed;
		return output;
	}


	
	
	
	
	
	
	


	/**
	 *  Update graphs using new computed results--------------------------------------------------------------------------------.
	 */		
	public void actualizeFirstPlots() {
		if(T1T2MixedEstimation) {actualizeFirstPlotsT1T2();return;}
		else {
			if(hasT1)actualizeFirstPlotsT1();
			if(hasT2)actualizeFirstPlotsT2();
			return;
		}
	}
			
	
	/**
	 * Actualize plots of T2 spectrum.
	 */
	public void actualizeFirstPlotsT2() {
		int incrT2=0;
		int deltaXt2=0;
		int bubbleSize=isBionanoDisplay ? 30 : 40;

		if(this.autoSizingTimePlots) maxPlotYT2=VitimageUtils.max(dataTimelapseT2[tCor])*1.6;

		double maxCurve=0;
		if(switchT1T2==1) maxCurve=VitimageUtils.max(tabFittenT2BicompCute[tCor][0]);
		if(switchT1T2==0) maxCurve=VitimageUtils.max(tabFittenT2MonoCute[tCor][0]);
		if(maxCurve>0.98*maxPlotYT1T2)maxPlotYT1T2=maxCurve*1.2;
		
		//BUBBLE SHOWING CURRENT ECHO
		plotT2.setColor(bubbleColor);
        plotT2.setLineWidth(bubbleSize);
        plotT2.replace(incrT2++, "line",new double[]{t2Times[tCor][zCor][cEchoes],t2Times[tCor][zCor][cEchoes]},new double[]{dataTimelapseT2[tCor][cEchoes],dataTimelapseT2[tCor][cEchoes]});

		
        //NOISE
        plotT2.setLineWidth(1);
		double[]statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMap.tabSigmasT2Seq[tCor][zCor]+deltaSigma);
		double nSum=statusRoi==2 ? Math.sqrt(this.nPtsCur) : (1+crossWidth);
		this.meanNoiseT2Cur=statsNoise[0];
		this.sigmaNoiseT2Cur=statsNoise[1];
		plotT2.setColor(Color.black);
        plotT2.setLineWidth(2);
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2OnPlot},new double[]{statsNoise[0],statsNoise[0]});//Afficher le mean+sigma
        plotT2.setLineWidth(1);
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2OnPlot},new double[]{statsNoise[0]+statsNoise[1]/nSum,statsNoise[0]+statsNoise[1]/nSum});//Afficher le mean+sigma
		plotT2.replace(incrT2++, "line",new double[]{0,maxT2OnPlot},new double[]{statsNoise[0]-statsNoise[1]/nSum,statsNoise[0]-statsNoise[1]/nSum});//Afficher le mean-sigma


		//OBSERVATIONS IRM MEAN AND STDEV
        plotT2.setLineWidth(2);
		plotT2.setColor(crossColor );
		plotT2.replace(incrT2++,"x",t2selectCute(t2Times[tCor][zCor],tCor),t2selectCute(dataTimelapseT2[tCor],tCor));//Afficher les points IRM
        plotT2.setLineWidth(1);
		for(int t=0;t<t2selectCute(this.t2Times[tCor][zCor],tCor).length;t++)plotT2.replace(incrT2++, "line",new double[]{t2selectCute(this.t2Times[tCor][zCor],tCor)[t]-deltaXt2,t2selectCute(this.t2Times[tCor][zCor],tCor)[t]-deltaXt2},new double[]{t2selectCute(dataTimelapseT2[tCor],tCor)[t]-t2selectCute(dataTimelapseT2Sigmas[tCor],tCor)[t],t2selectCute(dataTimelapseT2[tCor],tCor)[t]+t2selectCute(dataTimelapseT2Sigmas[tCor],tCor)[t]});//Afficher le T2

		
		if(switchT1T2==0) {	
	        //FIT T2 Mono
	        plotT2.setColor(curveT2Bicomp);
	        plotT2.setLineWidth(1+thickCurve);
//	        plotT2.setColor((jitterT2Mono[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=0;tran<timesT2Cute[tCor][zCor].length;tran++) {
		        if(true || t2SelectedCurve(tran,tCor))plotT2.replace(incrT2++, "line",timesT2Cute[tCor][zCor][tran][2],tabFittenT2MonoCute[tCor][tran]);//Afficher la courbe Bicomp
	        }
		}        
        

		else  {	
	        //FIT T2 Mono
	        plotT2.setColor(curveT2Bicomp);
	        plotT2.setLineWidth(1+thickCurve);
//	        plotT2.setColor((jitterT2Bicomp[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=0;tran<timesT2Cute[tCor][zCor].length;tran++) {
		        if(true || t2SelectedCurve(tran,tCor)) {
		        	System.out.println("IT IS ");
		        	plotT2.replace(incrT2++, "line",timesT2Cute[tCor][zCor][tran][2],tabFittenT2BicompCute[tCor][tran]);//Afficher la courbe Bicomp
		        }
	        }
		}        

        //LEGENDE
		plotT2.setLineWidth(1);
		String strLegendT2="\n\nNoise +-sigma\n\nSpin-echo signal"+"\nData std over Roi\n";
		
		plotT2.setLimits(0, maxT2OnPlot, 0,maxPlotYT2);
		plotT2.setColor(new Color(150,150 ,150) );
		plotT2.addLegend(strLegendT2,"bottom-up");
		plotCan2.setPlot(plotT2);

		plotT2.setLineWidth(1);
		for(int cur=incrT2;cur<lastCountT2;cur++)plotT2.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		plotT2.setLineWidth(1);
		for(int cur=incrT2;cur<lastCount2;cur++)plotT2.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		lastCountT2=incrT2;
	}

	/**
	 * Actualize first plots T 1 T 2.
	 */
	public void actualizeFirstPlotsT1T2() {
		int incrT1T2=0;
		int deltaXt2=3;
		int bubbleSize=isBionanoDisplay ? 30 : 40;

		if(this.autoSizingTimePlots) maxPlotYT1T2=VitimageUtils.max(dataTimelapseT1T2[tCor])*1.6;

		double maxCurve=0;
		if(switchT1T2==1) maxCurve=VitimageUtils.max(tabFittenT1T2BicompCute[tCor][0]);
		if(switchT1T2==0) maxCurve=VitimageUtils.max(tabFittenT1T2MonoCute[tCor][0]);
		if(maxCurve>0.98*maxPlotYT1T2)maxPlotYT1T2=maxCurve*1.2;
		
		//BUBBLE SHOWING CURRENT ECHO
		plotT1T2.setColor(bubbleColor);
        plotT1T2.setLineWidth(bubbleSize);
        plotT1T2.replace(incrT1T2++, "line",new double[]{t1t2Times[tCor][zCor][cEchoes],t1t2Times[tCor][zCor][cEchoes]},new double[]{dataTimelapseT1T2[tCor][cEchoes],dataTimelapseT1T2[tCor][cEchoes]});

		
        //NOISE
        plotT1T2.setLineWidth(1);
		double[]statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMap.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma);
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
		plotT1T2.replace(incrT1T2++,"x",t1t2selectCute(t1t2Times[tCor][zCor],tCor),t1t2selectCute(dataTimelapseT1T2[tCor],tCor));//Afficher les points IRM
        plotT1T2.setLineWidth(1);
		for(int t=0;t<t1t2selectCute(this.t1t2Times[tCor][zCor],tCor).length;t++)plotT1T2.replace(incrT1T2++, "line",new double[]{t1t2selectCute(this.t1t2Times[tCor][zCor],tCor)[t]-deltaXt2,t1t2selectCute(this.t1t2Times[tCor][zCor],tCor)[t]-deltaXt2},new double[]{t1t2selectCute(dataTimelapseT1T2[tCor],tCor)[t]-t1t2selectCute(dataTimelapseT1T2Sigmas[tCor],tCor)[t],t1t2selectCute(dataTimelapseT1T2[tCor],tCor)[t]+t1t2selectCute(dataTimelapseT1T2Sigmas[tCor],tCor)[t]});//Afficher le T1T2

		
		if(switchT1T2==0) {
	        //FIT T1 Mono
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1T2Mono[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][0][2],tabFittenT1T2MonoCute[tCor][0]);//Afficher la courbe 
	
	        //FIT T1T2 Mono
	        plotT1T2.setColor((jitterT1T2Mono[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=1;tran<timesT1T2Cute[tCor][zCor].length;tran++) {
		        if(t1t2SelectedCurve(tran,tCor))plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][tran][2],tabFittenT1T2MonoCute[tCor][tran]);//Afficher la courbe Bicomp
	        }
		}        
        

		else  {
	        //FIT T1 Mono
	        plotT1T2.setLineWidth(1+thickCurve);
	        plotT1T2.setColor((jitterT1T2Bicomp[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][0][2],tabFittenT1T2BicompCute[tCor][0]);//Afficher la courbe 
	
	        //FIT T1T2 Mono
	        plotT1T2.setColor((jitterT1T2Bicomp[tCor]>=9999) ? Color.gray : curveT2Bicomp);
	        for(int tran=1;tran<timesT1T2Cute[tCor][zCor].length;tran++) {
		        if(t1t2SelectedCurve(tran,tCor))plotT1T2.replace(incrT1T2++, "line",timesT1T2Cute[tCor][zCor][tran][2],tabFittenT1T2BicompCute[tCor][tran]);//Afficher la courbe Bicomp
	        }
		}        

        //LEGENDE
		plotT1T2.setLineWidth(1);
		String strLegendT1T2="\n\nNoise +-sigma\n\nSpin-echo signal"+"\nData std over Roi\n";
		
		plotT1T2.setLimits(0, maxT1T2OnPlot, 0,maxPlotYT1T2);
		plotT1T2.setColor(new Color(150,150 ,150) );
		plotT1T2.addLegend(strLegendT1T2,"bottom-up");
		plotCanT1T2.setPlot(plotT1T2);

		plotT1T2.setLineWidth(1);
		for(int cur=incrT1T2;cur<lastCountT1T2;cur++)plotT1T2.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		plotT1T2.setLineWidth(1);
		for(int cur=incrT1T2;cur<lastCount2;cur++)plotT1T2.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		lastCountT1T2=incrT1T2;
	}

	/**
	 * Actualize plots of T1 spectrum.
	 */
	public void actualizeFirstPlotsT1() {
		int incrT1=0;
		int deltaXt2=3;
		int bubbleSize=isBionanoDisplay ? 30 : 40;

		if(this.autoSizingTimePlots) maxPlotYT1=VitimageUtils.max(dataTimelapseT1[tCor])*1.6;

		double maxCurve=0;
		maxCurve=VitimageUtils.max(tabFittenT1MonoCute[tCor][0]);
		if(maxCurve>0.98*maxPlotYT1)maxPlotYT1=maxCurve*1.2;
		
		//BUBBLE SHOWING CURRENT ECHO
		plotT1.setColor(bubbleColor);
        plotT1.setLineWidth(bubbleSize);
        plotT1.replace(incrT1++, "line",new double[]{t1Times[tCor][zCor][cEchoes],t1Times[tCor][zCor][cEchoes]},new double[]{dataTimelapseT1[tCor][cEchoes],dataTimelapseT1[tCor][cEchoes]});

		
        //NOISE
        plotT1.setLineWidth(1);
		double[]statsNoise=RiceEstimator.computeSigmaAndMeanBgFromRiceSigmaStatic(hyperMap.tabSigmasT1Seq[tCor][zCor]+deltaSigma);
		double nSum=statusRoi==2 ? Math.sqrt(this.nPtsCur) : (1+crossWidth);
		this.meanNoiseT1Cur=statsNoise[0];
		this.sigmaNoiseT1Cur=statsNoise[1];
		plotT1.setColor(Color.black);
        plotT1.setLineWidth(2);
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1OnPlot},new double[]{statsNoise[0],statsNoise[0]});//Afficher le mean+sigma
        plotT1.setLineWidth(1);
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1OnPlot},new double[]{statsNoise[0]+statsNoise[1]/nSum,statsNoise[0]+statsNoise[1]/nSum});//Afficher le mean+sigma
		plotT1.replace(incrT1++, "line",new double[]{0,maxT1OnPlot},new double[]{statsNoise[0]-statsNoise[1]/nSum,statsNoise[0]-statsNoise[1]/nSum});//Afficher le mean-sigma


		//OBSERVATIONS IRM MEAN AND STDEV
        plotT1.setLineWidth(2);
		plotT1.setColor(crossColor );
		//plotT1.replace(incrT1++,"x",t1selectCute(t1Times[tCor][zCor],tCor),t1selectCute(dataTimelapseT1[tCor],tCor));//Afficher les points IRM
 		plotT1.replace(incrT1++,"x",t1Times[tCor][zCor],dataTimelapseT1[tCor]);//Afficher les points IRM
        plotT1.setLineWidth(1);
	//	for(int t=0;t<t1selectCute(this.t1Times[tCor][zCor],tCor).length;t++)plotT1.replace(incrT1++, "line",new double[]{t1selectCute(this.t1Times[tCor][zCor],tCor)[t]-deltaXt2,t1selectCute(this.t1Times[tCor][zCor],tCor)[t]-deltaXt2},new double[]{t1selectCute(dataTimelapseT1[tCor],tCor)[t]-t1selectCute(dataTimelapseT1Sigmas[tCor],tCor)[t],t1selectCute(dataTimelapseT1[tCor],tCor)[t]+t1selectCute(dataTimelapseT1Sigmas[tCor],tCor)[t]});//Afficher le T2
		for(int t=0;t<this.t1Times[tCor][zCor].length;t++)plotT1.replace(incrT1++, "line",new double[]{this.t1Times[tCor][zCor][t]-deltaXt2,this.t1Times[tCor][zCor][t]-deltaXt2},new double[]{dataTimelapseT1[tCor][t]-dataTimelapseT1Sigmas[tCor][t],dataTimelapseT1[tCor][t]+dataTimelapseT1Sigmas[tCor][t]});//Afficher le T2

		
		if(switchT1T2==0) {
	        //FIT T1 Mono
	        plotT1.setLineWidth(1+thickCurve);
	        plotT1.setColor(curveT1Bicomp);        
//	        plotT1.setColor((jitterT1Mono[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1.replace(incrT1++, "line",timesT1Cute[tCor][zCor][0][2],tabFittenT1MonoCute[tCor][0]);//Afficher la courbe 
	
		}        
        

		else  {
	        //FIT T1 Mono
	        plotT1.setLineWidth(1+thickCurve);
	        plotT1.setColor(curveT1Bicomp);        
//	        plotT1.setColor((jitterT1Bicomp[tCor]>=9999) ? Color.gray : curveT1Bicomp);        
	        plotT1.replace(incrT1++, "line",timesT1Cute[tCor][zCor][0][2],tabFittenT1BicompCute[tCor][0]);//Afficher la courbe 
	
		}        

        //LEGENDE
		plotT1.setLineWidth(1);
		String strLegendT2="\n\nNoise +-sigma\n\nSpin-echo signal"+"\nData std over Roi\n";
		
		plotT1.setLimits(0, maxT1OnPlot, 0,maxPlotYT1);
		plotT1.setColor(new Color(150,150 ,150) );
		plotT1.addLegend(strLegendT2,"bottom-up");
		plotCan1.setPlot(plotT1);

		plotT1.setLineWidth(1);
		for(int cur=incrT1;cur<lastCountT1;cur++)plotT1.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		plotT1.setLineWidth(1);
		for(int cur=incrT1;cur<lastCount2;cur++)plotT1.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		lastCountT1=incrT1;
	}

	
	
	/**
	 * Actualize Secondary plots.
	 */
	public void actualizeSecondPlots() {
		int incrT1=0;
		int incrT2=0;
		int bubbleSize=isBionanoDisplay ? 30 : 50;

        //HORIZONTAL GRID T21 ET T22
   		for(int tim=0;tim<this.nTimepoints;tim++) {
			plotT21.setLineWidth(sexyForBioInfo ? 2: 4);
			plotT21.setColor(Color.darkGray);
			plotT21.replace(incrT1++, "line",new double[] {t0T1,t1T1},new double[] {tim,tim});
			plotT22.setLineWidth(sexyForBioInfo ? 2: 4);
			plotT22.setColor(Color.darkGray);
			plotT22.replace(incrT2++, "line",new double[] {t0T2,t1T2},new double[] {tim,tim});
   		}

        //VERTICAL GRID T21
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT21.setLineWidth(1);        
			plotT21.setColor(Color.lightGray);
        	for(int mul=2;mul<=9;mul++) {
    	        plotT21.setLineWidth(1+(mul-1)%2);        
        		plotT21.replace(incrT1++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.nTimepoints});
        	}
/*			plotT21.setColor(Color.white);
        	plotT21.replace(incrT1++, "line",new double[] {tt*3,tt*3},new double[] {0,this.nTimepoints});
	*/        plotT21.setLineWidth(3);        
			plotT21.setColor(Color.darkGray);
        	plotT21.replace(incrT1++, "line",new double[] {tt,tt},new double[] {0,this.nTimepoints});
        }
        //VERTICAL GRID T22
        for(int tt=t0T2;tt<t1T1;tt*=10) {
	        plotT22.setLineWidth(1);        
			plotT22.setColor(Color.lightGray);
        	for(int mul=2;mul<=9;mul++) {
    	        plotT22.setLineWidth(1+(mul-1)%2);        
        		plotT22.replace(incrT2++, "line",new double[] {tt*mul,tt*mul},new double[] {0,this.nTimepoints});
        	}
			/*			plotT22.setColor(Color.white);
           	plotT22.replace(incrT2++, "line",new double[] {tt*3,tt*3},new double[] {0,this.nTimepoints});
*/	        plotT22.setLineWidth(3);        
			plotT22.setColor(Color.darkGray);
        	plotT22.replace(incrT2++, "line",new double[] {tt,tt},new double[] {0,this.nTimepoints});
        }
  
		double x0T1=t0T1*1.15;			double x0T2=t0T2*1.15;
		plotT21.setLineWidth(bubbleSize);
		plotT21.setColor(bubbleColor);
		double yBub=(!isFocusActivatedOnSpectrum? 0 : tCor*1.0/this.nTimepoints)+(tCor)+(0.2)/(isFocusActivatedOnSpectrum? this.nTimepoints: 1);
        plotT21.replace(incrT1++, "line",new double[] {x0T1,x0T1},new double[] {yBub,yBub});
        if(!sexyForBioInfo)        plotT22.setLineWidth(bubbleSize);
        if(!sexyForBioInfo)plotT22.setColor(bubbleColor);
        if(!sexyForBioInfo)plotT22.replace(incrT2++, "line",new double[] {x0T2,x0T2},new double[] {yBub,yBub});



		if(firstDrawBubbles) {
			firstDrawBubbles=false;
			for(int tt=0;tt<this.nTimepoints;tt++) {
				double bubblePosFrameX=isBionanoDisplay ? ((2)*1.0/TARGET_WIDTH_FOR_PLOTS_ALONE) : ((12)*1.0/TARGET_WIDTH_FOR_PLOTS_ALONE);
				double[]vals=spectrumRangingModeT1 ? plotT21.getLimits() : plotT22.getLimits();
				double plotValDown=vals[2];
				double plotValUp=vals[3];
				double borderUp=isBionanoDisplay ? 14 : 18;
				double borderDown=isBionanoDisplay ? 32 : 47;
				double bubblePlotY=(((tt)+(0.2)/(isFocusActivatedOnSpectrum? this.nTimepoints: 1))-plotValDown)/(plotValUp-plotValDown);
				double bubblePosFrameY=((1-bubblePlotY));//*(TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS-borderUp-borderDown)+borderUp)/TARGET_HEIGHT_FOR_DOWN_PLOT_AND_IMGS;

				if(!sexyForBioInfo)plotT22.setColor(Color.black);
				if(!sexyForBioInfo)plotT22.addLabel(bubblePosFrameX, bubblePosFrameY, "D "+tt);incrT2++;
				if(!sexyForBioInfo)plotT21.setColor(Color.black);
				if(!sexyForBioInfo)plotT21.addLabel(bubblePosFrameX,bubblePosFrameY, "D "+tt);incrT1++;
			}
		}
		
		
   		for(int tim=0;tim<this.nTimepoints;tim++) {
	 		double radius=0;
	 		plotT21.setLineWidth(6);
	 		plotT22.setLineWidth(6);

	 		//Draw spectrum curves
			if(T1T2MixedEstimation) {
				plotT21.setLineWidth(2);
		 		plotT21.setColor(Color.red);
		 		plotT21.replace(incrT1++, "line", valsSpectrum[tim][2], valsSpectrum[tim][0]);
		 		plotT22.setLineWidth(2);
		 		if(sexyForBioInfo)plotT22.setColor(curveT2Bicomp);
		 		else plotT22.setColor(Color.green);
		 		plotT22.replace(incrT2++, "line", valsSpectrum[tim][2], valsSpectrum[tim][1]);
			}
			else {
				if(hasT1) {
					plotT21.setLineWidth(2);
			 		plotT21.setColor(Color.red);
			 		plotT21.replace(incrT1++, "line", valsSpectrumT1[tim][1], valsSpectrumT1[tim][0]);
				}
				if(hasT2) {
			 		plotT22.setLineWidth(2);
			 		if(sexyForBioInfo)plotT22.setColor(curveT2Bicomp);
			 		else plotT22.setColor(Color.green);
			 		plotT22.replace(incrT2++, "line", valsSpectrumT2[tim][1], valsSpectrumT2[tim][0]);
				}
			}
 
	 		
	 		//Draw ranging interval
	 		if(rangeDisplayedFullBlocks && this.spectrumRangingModeT1 && tim==tCor && (!sexyForBioInfo)) {
	 			plotT21.setColor(Color.blue);	 		
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[0]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[2],this.rangingBoundaries[2]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[1],this.rangingBoundaries[1]}, new double[] {tim+1-0.6,tim+1-0.4});
	   			plotT21.replace(incrT1++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[2]}, new double[] {tim+1-0.5,tim+1-0.5});
	 		}   			
	 		if(rangeDisplayedFullBlocks && this.spectrumRangingModeT2 && tim==tCor && (!sexyForBioInfo)) {
	 			plotT22.setColor(Color.blue);	 		
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[0]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[2],this.rangingBoundaries[2]}, new double[] {tim+1-0.9,tim+1-0.1});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[1],this.rangingBoundaries[1]}, new double[] {tim+1-0.6,tim+1-0.4});
	   			plotT22.replace(incrT2++, "line", new double[] {this.rangingBoundaries[0],this.rangingBoundaries[2]}, new double[] {tim+1-0.5,tim+1-0.5});
	 		}   		 		
   		}				
		if(! isFocusActivatedOnSpectrum) {
			System.out.println("Setting limits in multi mode. Update ok at "+initTimer.getTime()+" s\n");
			plotT21.setLimits(t0T1,t1T1, 0, nTimepoints);
			plotT22.setLimits(t0T2,t1T2, 0, nTimepoints);
		}
		else {
			System.out.println("Setting limits in focus mode");
			plotT21.setLimits(t0T1,t1T1, tCor, tCor+1);
			plotT22.setLimits(t0T2,t1T2, tCor, tCor+1);
		}

		plotT21.setLineWidth(1);
		for(int cur=incrT1;cur<lastCount21;cur++)plotT21.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		plotT22.setLineWidth(1);
		for(int cur=incrT2;cur<lastCount22;cur++)plotT22.replace(cur, "line", new double[] {-1,-1}, new double[] {-1,-1});
		lastCount21=incrT1;
		lastCount22=incrT2;
	}				

	/**
	 * Actualize spectrum curves when computing cross-fitting T1T2
	 */
	public void actualizeSpectrumCurvesT1T2() {
		for(int tim=0;tim<this.nTimepoints;tim++) {
			this.valsSpectrum[tim]=computeSpectrumCurveT1T2(tim);
		}
		double d1=0;
		for(int i=0;i<this.valsSpectrum[0][0].length;i++) {d1+=valsSpectrum[0][0][i];}
		if(this.separateNormalizationSpectrumMode==0) {
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double maxValT1=VitimageUtils.max(this.valsSpectrum[tim][0]);
				for(int i=0;i<this.valsSpectrum[tim][0].length;i++)this.valsSpectrum[tim][0][i]=tim+this.valsSpectrum[tim][0][i]/maxValT1*normValueSpectrum;
				double maxValT2=VitimageUtils.max(this.valsSpectrum[tim][1]);
				for(int i=0;i<this.valsSpectrum[tim][1].length;i++)this.valsSpectrum[tim][1][i]=tim+this.valsSpectrum[tim][1][i]/maxValT2*normValueSpectrum;
			}
		}
		else if(this.separateNormalizationSpectrumMode==1) {
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
		else if(this.separateNormalizationSpectrumMode==2) {
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double factorViewHistoT1=5;
				double factorViewHistoT2=10;
				for(int i=0;i<this.valsSpectrum[tim][0].length;i++)this.valsSpectrum[tim][0][i]=tim+this.valsSpectrum[tim][0][i]*factorViewHistoT1/(this.nPtsCur*MRUtils.maxM0BionanoForNormalization);
				for(int i=0;i<this.valsSpectrum[tim][1].length;i++)this.valsSpectrum[tim][1][i]=tim+this.valsSpectrum[tim][1][i]*factorViewHistoT2/(this.nPtsCur*MRUtils.maxM0BionanoForNormalization);
			}
		}
		double d=0;
		for(int i=0;i<this.valsSpectrum[0][0].length;i++) {d+=valsSpectrum[0][0][i];}
	}

	/**
	 * Actualize spectrum curves 
	 */
	public void actualizeSpectrumCurvesT1() {
		for(int tim=0;tim<this.nTimepoints;tim++) {
			if(this.valsSpectrumT1==null)this.valsSpectrumT1=new double[this.nTimepoints][][];
			if(this.valsSpectrumT1[tim]==null)this.valsSpectrumT1[tim]=new double[2][];
			this.valsSpectrumT1[tim]=computeSpectrumCurveT1(tim);
		}
		double d1=0;
		for(int i=0;i<this.valsSpectrumT1[0][0].length;i++) {d1+=valsSpectrumT1[0][0][i];}
		if(this.separateNormalizationSpectrumMode==0) {
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double maxValT1=VitimageUtils.max(this.valsSpectrumT1[tim][0]);
				for(int i=0;i<this.valsSpectrumT1[tim][0].length;i++)this.valsSpectrumT1[tim][0][i]=tim+this.valsSpectrumT1[tim][0][i]/maxValT1*normValueSpectrum;
			}
		}
		else if(this.separateNormalizationSpectrumMode==1) {
			double maxTotT1=0;
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double maxValT1=VitimageUtils.max(this.valsSpectrumT1[tim][0]);
				if(maxTotT1<maxValT1)maxTotT1=maxValT1;
			}
			for(int tim=0;tim<this.nTimepoints;tim++) {
				for(int i=0;i<this.valsSpectrumT1[tim][0].length;i++)this.valsSpectrumT1[tim][0][i]=tim+this.valsSpectrumT1[tim][0][i]/maxTotT1*normValueSpectrum;
			}
		}
		else if(this.separateNormalizationSpectrumMode==2) {
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double factorViewHistoT1=5;
				for(int i=0;i<this.valsSpectrumT1[tim][0].length;i++)this.valsSpectrumT1[tim][0][i]=tim+this.valsSpectrumT1[tim][0][i]*factorViewHistoT1/(this.nPtsCur*MRUtils.maxM0BionanoForNormalization);
			}
		}
	}

	/**
	 * Actualize spectrum curves T 2.
	 */
	public void actualizeSpectrumCurvesT2() {
		for(int tim=0;tim<this.nTimepoints;tim++) {
			if(this.valsSpectrumT2==null)this.valsSpectrumT2=new double[this.nTimepoints][][];
			if(this.valsSpectrumT2[tim]==null)this.valsSpectrumT2[tim]=new double[2][];
			this.valsSpectrumT2[tim]=computeSpectrumCurveT2(tim);
		}
		double d1=0;
		for(int i=0;i<this.valsSpectrumT2[0][0].length;i++) {d1+=valsSpectrumT2[0][0][i];}
		if(this.separateNormalizationSpectrumMode==0) {
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double maxValT2=VitimageUtils.max(this.valsSpectrumT2[tim][0]);
				for(int i=0;i<this.valsSpectrumT2[tim][0].length;i++)this.valsSpectrumT2[tim][0][i]=tim+this.valsSpectrumT2[tim][0][i]/maxValT2*normValueSpectrum;
			}
		}
		else if(this.separateNormalizationSpectrumMode==1) {
			double maxTotT2=0;
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double maxValT2=VitimageUtils.max(this.valsSpectrumT2[tim][0]);
				if(maxTotT2<maxValT2)maxTotT2=maxValT2;
			}
			for(int tim=0;tim<this.nTimepoints;tim++) {
				for(int i=0;i<this.valsSpectrumT2[tim][0].length;i++)this.valsSpectrumT2[tim][0][i]=tim+this.valsSpectrumT2[tim][0][i]/maxTotT2*normValueSpectrum;
			}
		}
		else if(this.separateNormalizationSpectrumMode==2) {
			for(int tim=0;tim<this.nTimepoints;tim++) {
				double factorViewHistoT2=10;
				for(int i=0;i<this.valsSpectrumT2[tim][0].length;i++)this.valsSpectrumT2[tim][0][i]=tim+this.valsSpectrumT2[tim][0][i]*factorViewHistoT2/(this.nPtsCur*MRUtils.maxM0BionanoForNormalization);
			}
		}
	}

	 

	/**
 * Actualize displayed numbers in the central tab
 */
public void actualizeDisplayedNumbers() {
		//Update colors
		boolean t1MonoActive=T1T2MixedEstimation ? (switchT1T2==0) : (hasT1 && (switchT1T2==0));
		boolean t2MonoActive=T1T2MixedEstimation ? (switchT1T2==0) : (hasT2 && (switchT1T2==0));
		boolean t1BiActive=T1T2MixedEstimation ? (switchT1T2==1) : (hasT1 && (switchT1T2==1));
		boolean t2BiActive=T1T2MixedEstimation ? (switchT1T2==1) : (hasT2 && (switchT1T2==1));
		
		
		titles[0].setBackground(t1MonoActive ? paramsT1ActivatedColor : paramsUnactivated);
		titles[1].setBackground(t1BiActive ? paramsT1ActivatedColor : paramsUnactivated);
		titles[4].setBackground(t2MonoActive ? paramsT2ActivatedColor : paramsUnactivated);
		titles[5].setBackground(t2BiActive ? paramsT2ActivatedColor : paramsUnactivated);

		for(int j=0;j<5;j++) {
			params[0][j].setBackground(t1MonoActive ? paramsT1ActivatedColor : paramsUnactivated);
			params[1][j].setBackground(t1MonoActive ? paramsT1TextActivatedColor : paramsUnactivated);
			params[2][j].setBackground(t1BiActive ? paramsT1ActivatedColor : paramsUnactivated);
			params[3][j].setBackground(t1BiActive ? paramsT1TextActivatedColor : paramsUnactivated);

			params[8][j].setBackground(t2MonoActive ? paramsT2ActivatedColor : paramsUnactivated);
			params[9][j].setBackground(t2MonoActive ? paramsT2TextActivatedColor : paramsUnactivated);
			params[10][j].setBackground(t2BiActive ? paramsT2ActivatedColor : paramsUnactivated);
			params[11][j].setBackground(t2BiActive ? paramsT2TextActivatedColor : paramsUnactivated);
		}
		//Update values

		params[1][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][0] : paramsTimelapseT1[this.tCor][0],this.normalizePDInMiddleTab ) ) );
		params[1][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][1] : paramsTimelapseT1[this.tCor][1]));
		params[1][4].setText(String.format("%5.4f",T1T2MixedEstimation ? khi2T1T2Mono[this.tCor] : this.khi2T1Mono[this.tCor]));

		params[3][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? (paramsTimelapseT1T2[this.tCor][6]+paramsTimelapseT1T2[this.tCor][7]) : paramsTimelapseT1[this.tCor][0],this.normalizePDInMiddleTab)));
		params[3][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][8] : paramsTimelapseT1[this.tCor][1]));
		params[3][4].setText(String.format("%5.4f",T1T2MixedEstimation ? khi2T1T2Bicomp[this.tCor] : this.khi2T1Mono[this.tCor]));
/*
		params[5][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][8] : paramsTimelapseT1[this.tCor][2])));
		params[5][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][10] : paramsTimelapseT1[this.tCor][3]));
		params[5][2].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][9] : paramsTimelapseT1[this.tCor][4])));
		params[5][3].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][11] : paramsTimelapseT1[this.tCor][5]));
		params[5][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1BiT2Bicomp[this.tCor] : this.khi2T1Bicomp[this.tCor]));

		params[7][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][14] : paramsTimelapseT1[this.tCor][2])));
		params[7][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][15] : paramsTimelapseT1[this.tCor][3]));
		params[7][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1T2Bionano[this.tCor] : this.khi2T1T2Bionano[this.tCor]));
*/
		params[9][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][0] : paramsTimelapseT2[this.tCor][0],this.normalizePDInMiddleTab)));
		params[9][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][2] : paramsTimelapseT2[this.tCor][1]));
		params[9][4].setText(String.format("%5.4f",T1T2MixedEstimation ? khi2T1T2Mono[this.tCor] : this.khi2T2Mono[this.tCor]));

		params[11][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][6] : paramsTimelapseT2[this.tCor][6],this.normalizePDInMiddleTab)));
		params[11][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][9] : paramsTimelapseT2[this.tCor][8]));
		params[11][2].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][7] : paramsTimelapseT2[this.tCor][7],this.normalizePDInMiddleTab)));
		params[11][3].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][10] : paramsTimelapseT2[this.tCor][9]));
		params[11][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1T2Bicomp[this.tCor] : this.khi2T2Bicomp[this.tCor]));
		/*
		params[13][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][8] : paramsTimelapseT2[this.tCor][2])));
		params[13][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][12] : paramsTimelapseT2[this.tCor][3]));
		params[13][2].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][9] : paramsTimelapseT2[this.tCor][4])));
		params[13][3].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][13] : paramsTimelapseT2[this.tCor][5]));
		params[13][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1BiT2Bicomp[this.tCor] : this.khi2T2Bicomp[this.tCor]));
*/		
		params[15][0].setText(String.format("%5.3f",MRUtils.convertBionanoM0InPercentage(T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][16] : paramsTimelapseT1[this.tCor][2],normalizePDInMiddleTab)));
		params[15][1].setText(String.format("%5.0f",T1T2MixedEstimation ? paramsTimelapseT1T2[this.tCor][17] : paramsTimelapseT1[this.tCor][3]));
		params[15][4].setText(String.format("%5.4f",T1T2MixedEstimation ? this.khi2T1T2Bionano[this.tCor] : this.khi2T1T2Bionano[this.tCor]));

		//Update texts
		textInfoEchoes.setText(this.echoesText[cEchoes][zCor][tCor]);
    	//textInfoMaps.setText(this.mapsText[cMaps][zCor][tCor]+" X="+xCor+" Y="+yCor);
    	textSam.setText("       Sample size ("+(1+2*crossWidth)+"x"+(1+2*crossWidth)+"x"+(1+2*crossThick)+")");

	}

	
	/**
	 * Export fitting data and spectrum to csv.
	 *
	 * @param dirPath export directory
	 * @param basename the basename of the csv file to export
	 */
	//TODO
	public void exportDataToCsv(String dirPath,String basename) {
		String[][]tab;
		boolean debug=false;
		if(T1T2MixedEstimation) {
			//Export T1 T2 data			
			int maxPts=0;
			int T=nTimepoints;
			for(int t=0;t<T;t++) {
				if(t1t2TrTimes[t][0].length>maxPts)maxPts=t1t2TrTimes[t][0].length;
			}

			tab=new String[1+3*T][maxPts+1];
			tab[0][0]="Label ";for(int i=0;i<maxPts;i++)tab[0][1+i]="Value_"+i;
			for(int t=0;t<T;t++) {
				int n=t1t2TrTimes[t][0].length;
				tab[1+3*t][0]="Day "+t+" TR value";for(int i=0;i<n;i++)tab[1+3*t][1+i]=""+t1t2TrTimes[t][0][i];
				tab[2+3*t][0]="Day "+t+" TE value";for(int i=0;i<n;i++)tab[2+3*t][1+i]=""+t1t2TeTimes[t][0][i];
				tab[3+3*t][0]="Day "+t+" Magnitude value";for(int i=0;i<n;i++)tab[3+3*t][1+i]=""+dataTimelapseT1T2[t][i];
			}
			if(debug) {
				System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
			}
			else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T1T2Seq_magnitude_data.csv");
			//Line 1 : Day 0 "Nameday" - TR value ; TR1 ;  TR2 ;  TR3 ....
			//Line 2 : Day 0 "Nameday" - TE value ; TE1 ;  TE2 ;  TE3 ....
			//Line 3 : Day 0 "Nameday" - Magnitude value ; TE1 ;  TE2 ;  TE3 ....
			//Line 4 : Day 1 "Nameday" - TR value ; TR1 ;  TR2 ;  TR3 ....
			//...

			
			
			tab=new String[1+T][11];
			String []labels=new String[] {"Label","PD_tot", "T1", "T2", "Khi2_mono","PD_short", "PD_long","T1", "T2_short", "T2_long", "Khi2_biexp"};
			for(int i=0;i<labels.length;i++)tab[0][i]=labels[i];
			for(int t=0;t<T;t++) {
				labels=new String[] {"Day "+t,""+paramsTimelapseT1T2[t][0],""+paramsTimelapseT1T2[t][1],""+paramsTimelapseT1T2[t][2],""+this.khi2T1T2Mono[t],
						""+paramsTimelapseT1T2[t][6],""+paramsTimelapseT1T2[t][7],""+paramsTimelapseT1T2[t][8],
						""+paramsTimelapseT1T2[t][9],""+paramsTimelapseT1T2[t][10],""+this.khi2T1T2Bicomp[t]};
				for(int l=0;l<labels.length;l++)tab[1+t][l]=labels[l];	
			}

			if(debug) {
				System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
			}
			else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T1T2Seq_estimated_params.csv");
			//Export T1 T2 estimated params			
			//Line 0 : Label  ; PD_tot ; T1 ; T2 ; Khi2_mono ; PD_short ; PD_long ; T1 ; T2_short ; T2_long ; Khi2_biexp 
			//Line 1 : Day 0 ; val ; val ; val ; val ; val
			//Line 1 : Day 1 ; val ; val ; val ; val ; val

			
			
			tab=new String[1+4*T][this.numberBins+1];
			labels=new String[] {"Label","Val_i"};
			for(int i=0;i<labels.length;i++)tab[0][i]=labels[i];
			for(int t=0;t<T;t++) {
				tab[1+4*t][0]="Day "+t+" T1_value(x_axis)";for(int i=0;i<this.numberBins;i++)tab[1+4*t][1+i]=""+valsSpectrum[t][2][i];
				tab[2+4*t][0]="Day "+t+" PD_weighted_T1_distribution(y_axis)";for(int i=0;i<this.numberBins;i++)tab[2+4*t][1+i]=""+(valsSpectrum[t][0][i]-t);
				tab[3+4*t][0]="Day "+t+" T2_value(x_axis)";for(int i=0;i<this.numberBins;i++)tab[3+4*t][1+i]=""+valsSpectrum[t][2][i];
				tab[4+4*t][0]="Day "+t+" PD_weighted_T2_distribution(y_axis)";for(int i=0;i<this.numberBins;i++)tab[4+4*t][1+i]=""+(valsSpectrum[t][1][i]-t);
			}
			if(debug) {
				System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
			}
			else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T1T2Seq_estimated_T1T2distributions.csv");
			//Export spectrum		
			//Line 0 : Label  ; val ; val ; val 
			//Day 0 T1_value(x_axis) ; 10 ; 11 ; val ; val ; val
			//Day 0 PD_weighted_T1_distribution(y_axis) ; 0.8 ; 0.9 ; 0.7
			//Day 0 T2_value(x_axis) ; 10 ; 11 ; val ; val ; val
			//Day 0 PD_weighted_T2_distribution(y_axis) ; 0.8 ; 0.9 ; 0.7
		}
		
		else {
			if(hasT1) {
				//Export T1 data			
				int maxPts=0;
				int T=nTimepoints;
				for(int t=0;t<T;t++) {
					if(t1TrTimes[t][0].length>maxPts)maxPts=t1TrTimes[t][0].length;
				}

				tab=new String[1+2*T][maxPts+1];
				tab[0][0]="Label ";for(int i=0;i<maxPts;i++)tab[0][1+i]="Value_"+i;
				for(int t=0;t<T;t++) {
					int n=t1TrTimes[t][0].length;
					tab[1+2*t][0]="Day "+t+" TR value";for(int i=0;i<n;i++)tab[1+2*t][1+i]=""+t1TrTimes[t][0][i];
					tab[2+2*t][0]="Day "+t+" Magnitude value";for(int i=0;i<n;i++)tab[2+2*t][1+i]=""+dataTimelapseT1[t][i];
				}
				if(debug) {
					System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
				}
				else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T1Seq_magnitude_data.csv");
				//Line 1 : Day 0 "Nameday" - TR value ; TR1 ;  TR2 ;  TR3 ....
				//Line 3 : Day 0 "Nameday" - Magnitude value ; TE1 ;  TE2 ;  TE3 ....
				//...

				
				
				tab=new String[1+T][4];
				String []labels=new String[] {"Label","PD_tot", "T1",  "Khi2_mono"};
				for(int i=0;i<labels.length;i++)tab[0][i]=labels[i];
				for(int t=0;t<T;t++) {
					labels=new String[] {"Day "+t,""+paramsTimelapseT1[t][0],""+paramsTimelapseT1[t][1],""+this.khi2T1Mono[t]};
					for(int l=0;l<labels.length;l++)tab[1+t][l]=labels[l];	
				}
				if(debug) {
					System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
				}
				else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T1Seq_estimated_params.csv");
				//Export T1 T2 estimated params			
				//Line 0 : Label  ; PD_tot ; T1 ; T2 ; Khi2_mono ; PD_short ; PD_long ; T1 ; T2_short ; T2_long ; Khi2_biexp 
				//Line 1 : Day 0 ; val ; val ; val ; val ; val
				//Line 1 : Day 1 ; val ; val ; val ; val ; val

				
				
				tab=new String[1+2*T][this.numberBins+1];
				labels=new String[] {"Label","Val_i"};
				for(int i=0;i<labels.length;i++)tab[0][i]=labels[i];
				for(int t=0;t<T;t++) {
					tab[1+2*t][0]="Day "+t+" T1_value(x_axis)";for(int i=0;i<this.numberBins;i++)tab[1+2*t][1+i]=""+valsSpectrumT1[t][1][i];
					tab[2+2*t][0]="Day "+t+" PD_weighted_T1_distribution(y_axis)";for(int i=0;i<this.numberBins;i++)tab[2+4*t][1+i]=""+(valsSpectrumT1[t][0][i]-t);
				}
				if(debug) {
					System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
				}
				else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T1Seq_estimated_T1distributions.csv");
				//Export spectrum		
				//Line 0 : Label  ; val ; val ; val 
				//Day 0 T1_value(x_axis) ; 10 ; 11 ; val ; val ; val
				//Day 0 PD_weighted_T1_distribution(y_axis) ; 0.8 ; 0.9 ; 0.7
			}

			
			
			if(hasT2) {
				//Export T1 data			
				int maxPts=0;
				int T=nTimepoints;
				for(int t=0;t<T;t++) {
					if(t2TeTimes[t][0].length>maxPts)maxPts=t2TeTimes[t][0].length;
				}
				System.out.println("Maxpts="+maxPts);
				tab=new String[1+2*T][maxPts+1];
				tab[0][0]="Label ";for(int i=0;i<maxPts;i++)tab[0][1+i]="Value_"+i;
				for(int t=0;t<T;t++) {
					int n=t2TeTimes[t][0].length;
					System.out.println("n="+n);
					tab[1+2*t][0]="Day "+t+" TE value";for(int i=0;i<n;i++)tab[1+2*t][1+i]=""+t2TeTimes[t][0][i];
					tab[2+2*t][0]="Day "+t+" Magnitude value";for(int i=0;i<n;i++)tab[2+2*t][1+i]=""+dataTimelapseT2[t][i];
				}
				if(debug) {
					System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
				}
				else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T2Seq_magnitude_data.csv");
				//Line 1 : Day 0 "Nameday" - TR value ; TR1 ;  TR2 ;  TR3 ....
				//Line 3 : Day 0 "Nameday" - Magnitude value ; TE1 ;  TE2 ;  TE3 ....
				//...

				
				
				tab=new String[1+T][9];
				String []labels=new String[] {"Label","PD_tot", "T2",  "Khi2_mono","PD_short","PD_long", "T2_short", "T2_long",  "Khi2_biexp"};
				for(int i=0;i<labels.length;i++)tab[0][i]=labels[i];
				for(int t=0;t<T;t++) {
					labels=new String[] {"Day "+t,""+paramsTimelapseT2[t][0],""+paramsTimelapseT2[t][1],""+this.khi2T2Mono[t],
							""+paramsTimelapseT2[t][6],""+paramsTimelapseT2[t][7],""+paramsTimelapseT2[t][8],""+paramsTimelapseT2[t][9],""+this.khi2T2Bicomp[t]};
					for(int l=0;l<labels.length;l++)tab[1+t][l]=labels[l];	
				}
				if(debug) {
					System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
				}
				else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T2Seq_estimated_params.csv");
				//Export T1 T2 estimated params			
				//Line 0 : Label  ; PD_tot ; T1 ; T2 ; Khi2_mono ; PD_short ; PD_long ; T1 ; T2_short ; T2_long ; Khi2_biexp 
				//Line 1 : Day 0 ; val ; val ; val ; val ; val
				//Line 1 : Day 1 ; val ; val ; val ; val ; val

				
				
				tab=new String[1+2*T][this.numberBins+1];
				labels=new String[] {"Label","Val_i"};
				for(int i=0;i<labels.length;i++)tab[0][i]=labels[i];
				for(int t=0;t<T;t++) {
					tab[1+2*t][0]="Day "+t+" T2_value(x_axis)";for(int i=0;i<this.numberBins;i++)tab[1+2*t][1+i]=""+valsSpectrumT2[t][1][i];
					tab[2+2*t][0]="Day "+t+" PD_weighted_T2_distribution(y_axis)";for(int i=0;i<this.numberBins;i++)tab[2+4*t][1+i]=""+(valsSpectrumT2[t][0][i]-t);
				}
				if(debug) {
					System.out.println();for(int i=0;i<tab.length;i++) { System.out.println();for(int j=0;j<tab[i].length;j++) System.out.print("["+tab[i][j]+"]");}
				}
				else VitimageUtils.writeStringTabInCsv2(tab, dirPath+"_T2Seq_estimated_T2distributions.csv");
				//Export spectrum		
				//Line 0 : Label  ; val ; val ; val 
				//Day 0 T1_value(x_axis) ; 10 ; 11 ; val ; val ; val
				//Day 0 PD_weighted_T1_distribution(y_axis) ; 0.8 ; 0.9 ; 0.7
			}
			
			
		}
		IJ.showMessage("CSV data is now exported in "+new File(dirPath,basename).getAbsolutePath());
	}


	/**
	 * Update display of both results
	 */
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
		//actualizeExportSentence();
		repaintAll();
	}
	
	
	
		
		
	
	
	

	
	
	

	/**
	 *  Gui updating functions and callbacks--------------------------------------------------------------------------------.
	 *
	 * @param e the event
	 */	
	@Override
	public void keyTyped(KeyEvent e) {
		///////ACTIONS TO CHANGE SIZE OF AREA
		if (e.getKeyChar()=='p' && statusRoi!=2) {
			this.autoSizingTimePlots=!this.autoSizingTimePlots;
			actualizeMriObservationsBasedOnData();
			computeResultsAgain();
			displayResultsAgain();
		}
		if (e.getKeyChar()=='y') {
			this.imgCan1.getImage().setDisplayRange(0, 500);
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
			if((e.getKeyChar()=='9') && (cEchoes<hyperMap.C-1) )cEchoes+=1;
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
		if (e.getKeyChar()=='n') {
			normalizePDInMiddleTab=!normalizePDInMiddleTab;
			displayResultsAgain();						
		}
		if (e.getKeyChar()=='e') {
			IJ.log("\n\nExporting data in CSV files for importation in R, Python, or Matlab/Octave");
			String dirPath=null;
			String basename=null;
			if(true) {
				IJ.showMessage("Export data in CSV files. You will select a directory and a basename.");
				dirPath=VitiDialogs.chooseDirectoryUI("Choose directory for results export", "Directory chosen");
				basename=VitiDialogs.getStringUI("Choose basename (without extension)", "Basename :","my_fijirelax_csv" , false);
			}
			if((basename==null) || (dirPath==null) || (!new File(dirPath).exists())) {int a=1;} 
			else exportDataToCsv(dirPath,basename);
		}		
		System.out.println("Key pressed !");
		repaintAll();
	}
	
	
	/**
	 * Update switch buttons.
	 */
	public void updateSwitchButtons() {
		
		if(this.noiseHandling==0) {
			buttonExp.setForeground(Color.red);
			buttonOff.setForeground(Color.black);
			buttonRice.setForeground(Color.black);		
		}
		if(this.noiseHandling==1) {
			buttonExp.setForeground(Color.black);
			buttonOff.setForeground(Color.red);
			buttonRice.setForeground(Color.black);		
		}
		if(this.noiseHandling==2) {
			buttonExp.setForeground(Color.black);
			buttonOff.setForeground(Color.black);
			buttonRice.setForeground(Color.red);		
		}
		if(this.switchT1T2==0) {
			buttonSwitchMono.setForeground(Color.red);
			buttonSwitchBi.setForeground(Color.black);
		}
		if(this.switchT1T2==1) {
			buttonSwitchMono.setForeground(Color.black);
			buttonSwitchBi.setForeground(Color.red);
		}
	}
	
	
	
	/**
	 * Change tr. Button to move Tr
	 *
	 * @param up the up
	 */
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
	
	/**
	 * Action performed.
	 *
	 * @param e the event
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource()==buttonSwitchMono) {
			switchT1T2=0;System.out.println("Changing to "+switchT1T2);updateSwitchButtons();actualizeSpectrumCurves();		actualizeDisplayedNumbers();
			displayResultsAgain();return;
		}
		if(e.getSource()==buttonSwitchBi) {
			switchT1T2=1;System.out.println("Changing to "+switchT1T2);updateSwitchButtons();actualizeSpectrumCurves();		actualizeDisplayedNumbers();
			displayResultsAgain();return;
		}
		if(e.getSource()==buttonExp) {
			this.noiseHandling=0;
			updateSwitchButtons();this.computeResultsAgain();actualizeSpectrumCurves();		actualizeDisplayedNumbers();
			displayResultsAgain();return;
		}
		if(e.getSource()==buttonExport) {
			IJ.log("\n\nExporting data in CSV files for importation in R, Python, or Matlab/Octave");
			String dirPath=null;
			String basename=null;
			if(true) {
				IJ.showMessage("Export data in CSV files. You will select a directory and a basename.");
				dirPath=VitiDialogs.chooseDirectoryUI("Choose directory for results export", "Directory chosen");
				basename=VitiDialogs.getStringUI("Choose basename (without extension)", "Basename :","my_fijirelax_csv" , false);
			}
			if((basename==null) || (dirPath==null) || (!new File(dirPath).exists())) {int a=1;} 
			else exportDataToCsv(dirPath,basename);
		}
		if(e.getSource()==buttonOff) {
			this.noiseHandling=1;
			updateSwitchButtons();this.computeResultsAgain();actualizeSpectrumCurves();		actualizeDisplayedNumbers();
			displayResultsAgain();return;
		}
		if(e.getSource()==buttonRice) {
			this.noiseHandling=2;
			updateSwitchButtons();this.computeResultsAgain();actualizeSpectrumCurves();		actualizeDisplayedNumbers();
			displayResultsAgain();return;
		}
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
			if((e.getSource()==butNextEcho) && (cEchoes<hyperMap.C-1) ) {System.out.println("Changing cEcho from "+cEchoes+" to "+cEchoes+1+" because limit="+(hyperMap.C-1));cEchoes+=1;}
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
				if(file.contains(".roi")) {
					this.userRoi=new Opener().openRoi(file);
					this.userRoiInitial=new Opener().openRoi(file);
					IJ.selectWindow("MRI Curve explorer V2");
					this.userRoi.setPosition(0);				
				}
				else this.imgRoi=IJ.openImage(file);
				statusRoi=2;
				this.imgRoi.show();
				actualizeMriObservationsBasedOnData();
				computeResultsAgain();
				displayResultsAgain();
				IJ.log("Please hit the 'r' strike again to quit the Roi mode");
			}
			else {
				statusRoi=0;
				this.userRoi=null;
				this.imgRoi=null;
				this.nPtsCur=1;
				actualizeMriObservationsBasedOnData();
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

	/**
	 * Gets the help string.
	 *
	 * @return the help string
	 */
	public String getHelpString() {
		return (
		"Curve explorer is a java interface dedicated to T1 and T2 characterization in tissues from packed 4D or 5D HyperMaps (X-Y-Z-Echo-Times).\nType the 'h' key at any time to display this help"+
		" \n"+
		" \n\n"+
		"* Click on T1 / T2 image to set the center of the area of interest. \n"+
		"     - The area center is showed as a cross.\n"+
		"     - This area is a 2D square, centered around the cross. Its size is ( (1+2*radiusXY) x (1+2*radiusXY))\n"+
		"     - Data from the whole cube is used for T1 / T2 estimation : the MRI spin-echo values are averaged over this area. \n"+
		"     - Use '+' or '-' or the mouse scroll to zoom in / out\n"+
		" \n\n"+
		"* The data can also be gathered from a user-defined Roi.\n"+
		"     - Use the button 'Custom ROI' to load a ROI.\n"+
		"     - Click elsewhere to go back to square ROI mode_\n"+
		" \n\n"+					
		"* The plots display the current data from the area of interest : \nand the associated estimated parameters, in mono-T2 or bi-T2 (tabs in the middle)\n"+
		"     - The parameters are computed using an exponential fit with correction of the rice noise (sigma rice estimation is handled automatically)\n"+
		"     - Default estimation is bi-exponential (short T2 and long T2). Clik on the tabs or the switch button to switch the model.\n"+
		"     - Upper :\n"+
		"         .T1/T2 recovery data from the successive recovery time end echo time\n"+
		"         .Use 't' to switch between thick-curves-mode and thin-curves-mode to check visually the estimation accuracy\n"+
		"         .Use the mouse scroll or click on a cross to navigate echoes (see the yellow bubble running onto the curves)\n"+
		"     - Lower plots :\n"+
		"         .Left : estimated PD-weighted distrinbution of the T1 times in the region of interest from the successive timepoints (if any).\n"+
		"         .Right : estimated PD-weighted distrinbution of the T2 times in the region of interest from the successive timepoints (if any).\n"+
		"         .Use the mouse scroll or click on a curve to navigate observation times. Use the button 'show rangin' to range data and explore the pixels associated to values of T1 or T2\n"+
		"  \n\n"+
		" \n"+
		"* Use the right button panel to move into the data and to change the shape of the ROI\n"+
		" \n"+
		"* Export data and computed parameters : hit 'e' to export current data in matlab/octave/python/r format in the ImageJ log window\n"+
		" \n"+
		"* The tabs in the middle display continuously precious data associated to each model, computed on the mean magnitude values over the ROI:\n"+
		" \n");
	}

	/**
	 * Smooth histo. Compute smoothing of spectrum curves given as a parameter
	 *
	 * @param histo the histo
	 * @param sig the sig
	 * @return the double[]
	 */
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
	
	/**
	 * Mouse pressed.
	 *
	 * @param e the e
	 */
	@Override
	public void mousePressed(MouseEvent e) {
	}
	
	/**
	 * Mouse released.
	 *
	 * @param e the e
	 */
	@Override
	public void mouseReleased(MouseEvent e) {
	}
	
	/**
	 * Mouse entered.
	 *
	 * @param e the e
	 */
	@Override
	public void mouseEntered(MouseEvent e) {
	}

	/**
	 * Mouse exited.
	 *
	 * @param e the e
	 */
	@Override
	public void mouseExited(MouseEvent e) {
	}
		
	/**
	 * Key pressed.
	 *
	 * @param e the e
	 */
	@Override
	public void keyPressed(KeyEvent e) {
	}

    /**
     * Setup scroll.
     *
     * @param ox the ox
     * @param oy the oy
     */
    public void setupScroll(int ox, int oy) {
    	Rectangle srcRect = new Rectangle(0, 0, dims[0], dims[1]);
        xMouseStart = ox;
        yMouseStart = oy;
        xSrcStart = srcRect.x;
        ySrcStart = srcRect.y;
    }

    /**
     * Update the parameters displayed in the tab, depending on the current action
     *
     * @param param the param
     * @param action the action
     */
    public void updatePlaying(int param,int action) {
		System.out.println("Update playing "+param+" , "+action);
    	int type=0;		int nParams=  3; 
		paramsTimelapseT1T2[tCor][0]=paramsTimelapseT1T2[tCor][0]*(param!=0 ? 1 : (action==0 ? 1.07 : (1/1.07)));
		paramsTimelapseT1T2[tCor][1]=paramsTimelapseT1T2[tCor][1]*(param!=1 ? 1 : (action==0 ? 1.07 : (1/1.07)));
		paramsTimelapseT1T2[tCor][2]=paramsTimelapseT1T2[tCor][2]*(param!=2 ? 1 : (action==0 ? 1.07 : (1/1.07)));

		if(param==3) {
			this.deltaT2+=(action==0 ? -1 : 1);
			params[9][3].setText(""+(((this.t1t2TeTimes[tCor][dims[2]-1][this.t1t2TeTimes[tCor][dims[2]-1].length-1]>12.1*16  || this.deltaT2<-10)? (zCor==0 ? 0 : 0) : 0)+this.deltaT2));
			params[8][3].setText("DeltaT2");
			this.timesT1T2Cute=getCuteTimesT1T2();
		}

		if(param==4) {
			deltaSigma+=(action==0 ? -10 : 10);
			params[9][2].setText(""+deltaSigma);
			params[8][2].setText("Delta sig");
			params[7][2].setText(""+sigSmoBionano);
			params[6][2].setText("SigSmo");
			this.timesT1T2Cute=getCuteTimesT1T2();
		}

		
		
		System.out.println("Give : "+paramsTimelapseT1T2[tCor][0]+" , "+paramsTimelapseT1T2[tCor][1]+" , "+paramsTimelapseT1T2[tCor][2]+" , " +this.deltaT2);
		
		double []tabFitten=MRUtils.fittenRelaxationCurve(t1t2TrTimes[tCor][zCor],t1t2TeTimes[tCor][zCor],new double[] {paramsTimelapseT1T2[tCor][0],paramsTimelapseT1T2[tCor][1],paramsTimelapseT1T2[tCor][2]},hyperMap.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma,getFitType(MRUtils.T1T2_MONO_RICE));
		double[]accs=MRUtils.fittingAccuracies(dataTimelapseT1T2[tCor],t1t2TrTimes[tCor][zCor],t1t2TeTimes[tCor][zCor],hyperMap.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma,new double[] {paramsTimelapseT1T2[tCor][0],paramsTimelapseT1T2[tCor][1],paramsTimelapseT1T2[tCor][2]},getFitType(MRUtils.T1T2_MONO_RICE),false,riceEstimator,this.sigmaKhi2WeightedByNumberPoints ? this.nPtsCur : 1);		
		double [][]tabFittenCute=null;
		
		tabFittenCute=new double[this.timesT1T2Cute[tCor][zCor].length][];
		for(int fi=0;fi<this.timesT1T2Cute[tCor][zCor].length;fi++) {
			tabFittenCute[fi]=MRUtils.fittenRelaxationCurve(this.timesT1T2Cute[tCor][zCor][fi][0],this.timesT1T2Cute[tCor][zCor][fi][1],new double[] {paramsTimelapseT1T2[tCor][0],paramsTimelapseT1T2[tCor][1],paramsTimelapseT1T2[tCor][2]},hyperMap.tabSigmasT1T2Seq[tCor][zCor]+deltaSigma,getFitType(MRUtils.T1T2_MONO_RICE));
		}
    	
   		tabFittenT1T2Mono[tCor]=tabFitten;
 		tabFittenT1T2MonoCute[tCor]=tabFittenCute;
 		khi2T1T2Mono[tCor]=accs[0];
 		pValT1T2Mono[tCor]=accs[1];
 		
 		displayResultsAgain();
    }
    
	/**
	 * Mouse clicked.
	 *
	 * @param e the e
	 */
	public void mouseClicked(MouseEvent e) {
		for(int i=0;i<titles.length;i++) {
			if(e.getSource()==titles[i]) {
				switchT1T2=i%4;System.out.println("Changing to "+switchT1T2);actualizeSpectrumCurves();displayResultsAgain();return;
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
					if(j==3 && (i==8)) updatePlaying(3,1);
					if(j==3 && (i==9)) updatePlaying(3,0);
					if(j==2 && (i==8)) updatePlaying(4,1);
					if(j==2 && (i==9)) updatePlaying(4,0);
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
			imgRoi=null;
			userRoi=null;
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
			double newTimer=timer.getTime();
			boolean doubleClick=((newTimer-memoryTimeClickPlot)<doubleClickDelay);
			memoryTimeClickPlot=newTimer;
			this.xMouseRange=e.getX();
			this.yMouseRange=e.getY();
			System.out.println("Click on Plot21 | Coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			spectrumRangingModeT1=true;	
			spectrumRangingModeT2=false;	
			int newTime=actualizeRangingBoundaries();
			if(tCor!=newTime && !isFocusActivatedOnSpectrum) {
				tCor=newTime;
				if(cEchoes>=this.t1t2Times[tCor][zCor].length)cEchoes=this.t1t2Times[tCor][zCor].length-1;
			}
			if(doubleClick)isFocusActivatedOnSpectrum=!isFocusActivatedOnSpectrum;
			identifyRangedData();			
			displayResultsAgain();
			return;
		}

		else if(currentCanvas==WIN_PLOT22) {
			double newTimer=timer.getTime();
			boolean doubleClick=((newTimer-memoryTimeClickPlot)<doubleClickDelay);
			memoryTimeClickPlot=newTimer;
			this.xMouseRange=e.getX();
			this.yMouseRange=e.getY();
			System.out.println("Click on Plot22 | Coordinates=("+ xMouseRange+","+yMouseRange+")"+"  |  Zoom level="+zoomLevel+" | x="+this.xCor+" y="+this.yCor+" z="+zCor+"t="+tCor);
			spectrumRangingModeT2=true;	
			spectrumRangingModeT1=false;	
			int newTime=actualizeRangingBoundaries();
			if(tCor!=newTime && !isFocusActivatedOnSpectrum) {
				tCor=newTime;
				if(cEchoes>=this.t1t2Times[tCor][zCor].length)cEchoes=this.t1t2Times[tCor][zCor].length-1;
			}
			if(doubleClick)isFocusActivatedOnSpectrum=!isFocusActivatedOnSpectrum;
			identifyRangedData();			
			displayResultsAgain();
			return;
		}

		else {
			displayResultsAgain();
		}
	}
		
	/**
	 * Key released.
	 *
	 * @param e the e
	 */
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

	/**
	 * Window closing.
	 *
	 * @param paramWindowEvent the param window event
	 */
	public void windowClosing(WindowEvent paramWindowEvent)
	  {
	    super.windowClosing(paramWindowEvent);
	    instance = null;
	  }

	/**
	 * Mouse wheel moved.
	 *
	 * @param e the e
	 */
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
			System.out.println("Scroll "+e.getUnitsToScroll());
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
			if (e.getUnitsToScroll()<0) { if (cEchoes<hyperMap.C-1)cEchoes+=1;}
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

	/**
	 * Mouse dragged.
	 *
	 * @param e the e
	 */
	@Override
	public void mouseDragged(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	/**
	 * Mouse moved.
	 *
	 * @param e the e
	 */
	@Override
	public void mouseMoved(MouseEvent e) {
		// TODO Auto-generated method stub 
		
	}

	

	
}
