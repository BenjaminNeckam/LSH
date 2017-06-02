package at.ac.univie.main;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RefineryUtilities;

import at.ac.unive.utils.Bucket;
import at.ac.unive.utils.CSVParser;
import at.ac.unive.utils.KMeans;
import at.ac.unive.utils.Plot;
import at.ac.unive.utils.Point;

public class Main {

	public static void main(String[] args) {
		ArrayList<Point> points = CSVParser.CSVToPoint("C:/Users/Benni/git/LSH/src/test/resources/LSH-nmi-small.csv");
		try {
			KMeans.lsh(points, null);

			// KMeans k = new KMeans();
			// System.out.println(">>>>>>>>>>>>>>>>Centroid-Method
			// 2<<<<<<<<<<<<<<<<<<");
			// long startTime = System.nanoTime();
			// ArrayList a = k.initCentroids(points,15);
			// double estimatedTime = (System.nanoTime() - startTime)/
			// 1000000000.0;
			// System.out.println("\nElapsed Time: " + estimatedTime + "
			// seconds");
			//
			// for(int i=0;i<a.size();++i){
			// System.out.println(a.get(i).toString());
			// }
			//
			// Plot scatterplotdemo4 = new Plot("K-Means Start",a, 6);
			// scatterplotdemo4.pack();
			// RefineryUtilities.centerFrameOnScreen(scatterplotdemo4);
			// scatterplotdemo4.setVisible(true);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();

		}
		// ArrayList<Float> point = new ArrayList<>();
		// Random random = new Random();
		// for(int i=0;i<50;i++){
		// point.add((-100 + (0+100)*random.nextFloat()));
		// }
		// Collections.sort(point, new Comparator<Float>() {
		// public int compare(Float o1, Float o2) {
		// if (o1 > o2) {
		// return -1;
		// }
		// if (o1 < o2) {
		// return 1;
		// }
		// return 0;
		// }
		// });
		//
		// for(Float value:point){
		// System.out.println(value);
		// }
		// for(Float value:point){
		// System.out.println(Float.hashCode(value));
		// }

	}
}
