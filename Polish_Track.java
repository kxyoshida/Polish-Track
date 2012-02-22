import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;

import ij.measure.*;
import ij.text.*;
import ij.io.*;
import java.text.DecimalFormat;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;


public class Polish_Track implements PlugInFilter {
    ImagePlus imp;
    static int hw = 3; // half of window size in pixels 
    double baseIntensity = 1.0;  // as percentile
    double peakIntensity = 95.0; // as percentile
    double baseRange = 25.0;   // Permitted half range of background level
    //    double minPeak = 67.0; // Minimal height of the peak
    double sigmaRatio = 2.0; // Ration of sigma to half window size

    public int setup(String arg, ImagePlus imp) {
	this.imp = imp;
	return DOES_8G + DOES_16 + DOES_32;
    }

    static {
	System.loadLibrary("MP2DGaussForTrack");
    }

    public native int calc(double[] params, double[] data);

    private double toD(boolean flag) {
	return flag ? 1.0 : 0.0;
    }

    private class params {
	public double a[];
	public double amin[];
	public double amax[];
	public boolean afix[];
	public boolean llim[];
	public boolean ulim[];

	public double exparams[];

	public double error[];
	public double chisq;
	public double dof;

	public params() {
	    a = new double[7];
	    amin = new double[7];
	    amax = new double[7];
	    afix = new boolean[7];
	    llim = new boolean[7];
	    ulim = new boolean[7];
	    exparams = new double[42];
	    error = new double[7];
	    chisq = 0;
	    dof = 0;
	}
	public double[] flatten() {
	    for (int i=0; i<7; i++) {
		exparams[i] = a[i];
		exparams[i+7] = toD(afix[i]);
		exparams[i+14] = amin[i];
		exparams[i+21] = toD(llim[i]);
		exparams[i+28] = amax[i];
		exparams[i+35] = toD(ulim[i]);
	    }
	    return exparams;
	}

	public void recallParams() {
	    for (int i=0; i<7; i++) {
		a[i] = exparams[i];
		error[i] = exparams[i+7];
	    }
	    chisq = exparams[14];
	    dof = exparams[15];
	    return;
	}
    }

    boolean withinCircle(double x, double y) {
	if ( (x-hw+1)*(x-hw+1)+(y-hw+1)*(y-hw+1)<=hw*hw ) 
	    return true; 
	else 
	    return false;
    }

    private void guessParams(double[] x, double[] y, double[] g, params par) 
    {
	DescriptiveStatistics ds = new DescriptiveStatistics(g);
	//	par.a[0] = ds.getPercentile(baseIntensity);
	par.a[1] = ds.getPercentile(peakIntensity) - par.a[0];
	double sumx = 0, sumy = 0;

	par.a[1] = (par.a[1]>1.0)? par.a[1] : 1.0;

	for (int i=0; i<x.length; i++) {
	    if (withinCircle(x[i], y[i])) {  // restriction added on 14 Feb 2012
		sumx += x[i]*g[i];
		sumy += y[i]*g[i];
	    }
	}
	par.a[2] = sumx/ds.getSum();
	par.a[3] = sumy/ds.getSum();
	double sumsqx = 0, sumsqy =0;
	for (int i=0; i<x.length; i++) {
	    if (withinCircle(x[i], y[i])) {    // restriction added on 14 Feb 2012
		sumsqx += (x[i]-par.a[2])*(x[i]-par.a[3])*g[i];
		sumsqy += (y[i]-par.a[3])*(y[i]-par.a[3])*g[i];
	    }
	}
	par.a[4] = par.a[5] = Math.sqrt(0.5*(sumsqx+sumsqy)/ds.getSum());
	return;
    }

    private void fixParams(params par) {

	/* Absolute Requirements: */
	par.afix[0] = false;
	par.afix[1] = false;
	par.afix[2] = false;
	par.afix[3] = false;
	par.afix[4] = par.afix[5] = false;
	par.afix[6] = false;

	/* For fine tunings */
	par.llim[0] = par.ulim[0] = true;
	par.amin[0] = par.a[0] - Math.abs(par.a[0]) * baseRange / 100.0;
	par.amax[0] = par.a[0] + Math.abs(par.a[0]) * baseRange / 100.0;

	par.llim[1] = true;
	par.amin[1] = 0;     //par.a[1] * minPeak / 100.0;
	par.ulim[1] = false;

	/* Absolute Requirements: */
	par.llim[2] = par.ulim[2] = true;
	par.llim[3] = par.ulim[3] = true;
	par.amin[2] = par.amin[3] = 0;
	par.amax[2] = 2.0 * hw;
	par.amax[3] = 2.0 * hw;

	par.llim[4] = par.llim[5] = true;
	par.amin[4] = par.amin[5] = 0;

	/* For fine tunings */
	par.ulim[4] = par.ulim[5] = true;
	par.amax[4] = sigmaRatio * hw;
	par.amax[5] = sigmaRatio * hw;

	/* Absolute Requirements: */
	par.llim[6] = par.ulim[6] = true;
	par.amin[6] = 0;
	par.amax[6] = 360;

	return;
    }


    public int callMP(ImageProcessor ip, params par) 
    {
	int width = ip.getWidth();
	int height = ip.getHeight();
	int m = width*height;
	double x[] = new double[m];
	double y[] = new double[m];
	double g[] = new double[m];

	float[] pixels = (float[])ip.getPixels();

	for (int i = 0; i<height; i++) {
	    for (int j = 0; j<width; j++) {
		int l=i*width+j;
		x[l] = j;
		y[l] = i;
		g[l] = pixels[l];
	    }
	}

	double points[] = new double[3*m];

	// concatenate y to x in points
	System.arraycopy(x, 0, points, 0, m);
	System.arraycopy(y, 0, points, m, m);
	System.arraycopy(g, 0, points, 2*m, m);

	guessParams(x, y, g, par);
	fixParams(par);

	return calc(par.flatten(), points);
    }

    private void fixResultsTable(ResultsTable rt, int slice, int tr, params par, int status) {
	if (tr >= 0 && tr < rt.getCounter()) {
	    //	    rt.setValue("slice", tr, slice);
	    rt.setValue("cx", tr, rt.getValueAsDouble(2, tr) - hw + par.a[2]);
	    rt.setValue("cy", tr, rt.getValueAsDouble(3, tr) - hw + par.a[3]);
	    //	    rt.setValue("sigmax", tr, par.a[4]);
	    //	    rt.setValue("sigmay", tr, par.a[5]);
	    //	    rt.setValue("basal", tr, par.a[0]);
	    //	    rt.setValue("peak", tr, par.a[1]);
	    //	    rt.setValue("chisq", tr, par.chisq);
	    //	    rt.setValue("dof", tr, par.dof);
	    rt.setValue("status", tr, status);
	} else {
	    IJ.log("No corresponding rows in Results Tabe.");
	}
	return;
    }

    private void fillResultsTable(ResultsTable rt, int slice, int tr) {
	if (tr >= 0 && tr < rt.getCounter()) {
	    //	    rt.setValue("slice", tr, slice);
	    rt.setValue("cx", tr, rt.getValueAsDouble(2, tr));
	    rt.setValue("cy", tr, rt.getValueAsDouble(3, tr));
	    //	    rt.setValue("sigmax", tr, 0);
	    //	    rt.setValue("sigmay", tr, 0);
	    //	    rt.setValue("basal", tr, 0);
	    //	    rt.setValue("peak", tr, 0);
	    //	    rt.setValue("chisq", tr, 0);
	    //	    rt.setValue("dof", tr, 0);
	    rt.setValue("status", tr, -1);
	} else {
	    IJ.log("No corresponding rows in Results Tabe.");
	}
	return;
    }

    public void run(ImageProcessor ip) {


	int width = ip.getWidth();
	int height = ip.getHeight();

	params par = new params();
	ImageStack stack = imp.getStack();
	int nSlices = stack.getSize();

	GenericDialog gd = new GenericDialog("MP_2D-Gaussian_Fit");
	gd.addNumericField("Half window size", hw, 0);
	gd.addNumericField("Background intensity in percentile", baseIntensity, 1);
	gd.addNumericField("+- (%)", baseRange, 1);
	gd.addNumericField("Peak intensity in percentile", peakIntensity, 1);
	//	gd.addNumericField("Minimal height of peaks (%)", minPeak, 1);
	gd.addNumericField("Ratio of sigma to half window", sigmaRatio, 1);
	gd.addCheckbox("Logging",true);
	gd.showDialog();
	if (gd.wasCanceled()) return;
	hw = (int) gd.getNextNumber();
	baseIntensity = (double) gd.getNextNumber();
	baseRange = (double) gd.getNextNumber();
	peakIntensity = (double) gd.getNextNumber();
	//	minPeak = (double) gd.getNextNumber();
	sigmaRatio = (double) gd.getNextNumber();
	boolean logging = gd.getNextBoolean();

	ResultsTable rt = ResultsTable.getResultsTable();
	if (rt == null) {
	    IJ.error("Can't find Results Table.");
	    return;
	}

	int nResults = rt.getCounter();
	float[] fr = rt.getColumn(1);

	for (int n = 1 ; n <= nSlices; n++) {
	    ip = stack.getProcessor(n).convertToFloat();

	    float[] pixels = (float[]) ip.getPixels();
	    DescriptiveStatistics ds = new DescriptiveStatistics();
	    for (int i = 0; i < pixels.length; i++)
		ds.addValue((double)pixels[i]);
	    par.a[0] = ds.getPercentile(baseIntensity);

	    for (int i = 0; i < nResults; i++) {
		if (n == (int)fr[i]) {
		    int X = (int)rt.getValueAsDouble(2, i);
		    int Y = (int)rt.getValueAsDouble(3, i);
		    if ( X>=hw && Y>=hw && X+hw<width && Y+hw<height ) {
			ip.setRoi(X-hw, Y-hw, 2*hw, 2*hw);
			int status = callMP(ip.crop(), par);
			if (status <1 || status > 8) {
			    for (int j = 0; j < 7; j++) {
				IJ.log("a_init["+j+"]= "+par.a[j]);
				IJ.log("amin["+j+"]="+par.amin[j]);
				IJ.log("amax["+j+"]="+par.amax[j]);
			    }
			    IJ.error("Native Method Failed =" + status);
			    return;
			}
			par.recallParams();
			if (par.a[2] == 0 || par.a[2] == hw || par.a[3] == 0 || par.a[3] == hw) {
			    fillResultsTable(rt, n, i);
			    if (logging) {
				for (int j = 0; j < 7; j++) 
				    IJ.log("a["+j+"]= "+par.a[j]+"+-"+par.error[j]);
				IJ.log("Native Method Apparently Diverged =" + status);
			    }
			} 
			else if (status>=1 && status<=8) {
			    fixResultsTable(rt, n, i, par, status);
			    if (logging) {
				for (int j = 0; j < 7; j++) 
				    IJ.log("a["+j+"]= "+par.a[j]+"+-"+par.error[j]);
				IJ.log("Native Method was successful =" + status);
			    }
			} 
		    } else {
			fillResultsTable(rt, n, i);
		    }
		}
	    }
	    rt.show("Results");
	}

    }
}
