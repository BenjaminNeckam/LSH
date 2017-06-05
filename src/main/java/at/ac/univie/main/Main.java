package at.ac.univie.main;

import java.util.ArrayList;

import org.jfree.ui.RefineryUtilities;

import at.ac.unive.utils.CSVParser;
import at.ac.unive.utils.KMeans;
import at.ac.unive.utils.Plot;
import at.ac.unive.utils.Point;

public class Main {

	public static void main(String[] args) {
		ArrayList<Point> points = CSVParser.CSVToPoint("LSH-nmi-corrected.csv");
		ArrayList<Integer> compare = CSVParser.CSVToPointCompare("LSH-nmi-compare.csv");
		
		if(args[0].equals("1")){
			try{
				KMeans.lloyd(points,15);
				ArrayList<Integer> result = KMeans.PointsToIntegerList(points);
				KMeans.NMI(compare, result);
				
				Plot scatterplotdemo4 = new Plot("K-Means Start",points, 15);
				scatterplotdemo4.pack();
				RefineryUtilities.centerFrameOnScreen(scatterplotdemo4);
				scatterplotdemo4.setVisible(true);
			}catch (Exception e) {
				 System.err.println("An exception was thrown (most of the time for an unknow reason which we couldn't fix)");

			}
		}else if(args[0].equals("2")){
			try{
				Float bucketSize;
				if(args[1].isEmpty()){
					bucketSize = null;
				}else{
					bucketSize = Float.valueOf(args[1]);
				}
				KMeans.lsh(points, bucketSize, (float)7.0);
				ArrayList<Integer> result2 = KMeans.PointsToIntegerList(points);
				KMeans.NMI(compare, result2);
			}catch (Exception e) {
				 System.err.println("An exception was thrown (most of the time for an unknow reason which we couldn't fix)");

			}
			
		}else if(args[0].equals("3")){
			try{
				Float bucketSize;
				if(args[1].isEmpty()){
					bucketSize = null;
				}else{
					bucketSize = Float.valueOf(args[1]);
				}
				KMeans.lloyd(points,15);
				ArrayList<Integer> result = KMeans.PointsToIntegerList(points);
				KMeans.NMI(compare, result);
				
				Plot scatterplotdemo4 = new Plot("K-Means Start",points, 15);
				scatterplotdemo4.pack();
				RefineryUtilities.centerFrameOnScreen(scatterplotdemo4);
				scatterplotdemo4.setVisible(true);
				
				KMeans.lsh(points, bucketSize, (float)7.0);
				ArrayList<Integer> result2 = KMeans.PointsToIntegerList(points);
				KMeans.NMI(compare, result2);
			}catch (Exception e) {
				 System.err.println("An exception was thrown (most of the time for an unknow reason which we couldn't fix)");

			}
		}
		
		
	}
}
