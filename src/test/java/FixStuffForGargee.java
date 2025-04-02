import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import io.github.rocsg.fijirelax.mrialgo.MRUtils;
import io.github.rocsg.fijirelax.mrialgo.RiceEstimator;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class FixStuffForGargee {
    public static void main(String[] args) {
        String imgPath="/home/rfernandez/Bureau/A_Test/FijiRelax_track_issues/gargee/data/206_hypermap_crop.tif";
        ImagePlus img=IJ.openImage(imgPath);
        int C=img.getNChannels();
        int Z=img.getNSlices();
        int T=img.getNFrames();


        double[]normFactors=new double[]{65956,61754,63562,57040};
//        normFactors=new double[]{30000};

        for(int c=0; c<C;c++){
            for(int z=0; z<Z;z++){
                for(int t=0; t<T;t++){
                    int sli=VitimageUtils.getCorrespondingSliceInHyperImage(img,c, z, t);
                    String chain=img.getStack().getSliceLabel(sli);
                    System.out.println(chain);
                    String[]strs=chain.split("_");
                    String totSig=strs[strs.length-1];
                    String valueString=totSig.split("=")[1];
                    double d=Double.parseDouble(valueString);
                    System.out.println(d);
                    d/=normFactors[t];
                    System.out.println("Toto1");
                    String totFutureLabel="";
                    for(int s=0;s<strs.length-1;s++)totFutureLabel=totFutureLabel+=strs[s]+"_";
                    totFutureLabel+="SIGMARICE="+VitimageUtils.dou(d);
                    img.getStack().setSliceLabel(totFutureLabel, sli);
                    System.out.println("Toto3");
                }
            }
        }
        img.show();
        IJ.saveAsTiff(img,imgPath);
    }




}
