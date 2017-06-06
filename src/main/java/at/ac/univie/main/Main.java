package at.ac.univie.main;

import java.net.URL;
import java.util.ArrayList;
import java.util.Properties;

import org.jfree.ui.RefineryUtilities;

import at.ac.unive.utils.CSVParser;
import at.ac.unive.utils.KMeans;
import at.ac.unive.utils.Plot;
import at.ac.unive.utils.Point;

public class Main {

	public static void main(String[] args) {
		

		ArrayList<Point> points;
		ArrayList<Integer> compare;
		
		
		if(args[0].equals("1")){
			try{
				points = CSVParser.CSVToPoint(args[1]);
				compare = CSVParser.CSVToPointCompare(args[2]);
				
				KMeans.lloyd(points,15);
				ArrayList<Integer> result = KMeans.PointsToIntegerList(points);
				KMeans.NMI(compare, result);
				
				Plot scatterplotdemo4 = new Plot("K-Means Start",points, 15);
				scatterplotdemo4.pack();
				RefineryUtilities.centerFrameOnScreen(scatterplotdemo4);
				scatterplotdemo4.setVisible(true);
			}catch (Exception e) {
				e.printStackTrace();
				 System.err.println("An exception was thrown (most of the time for an unknow reason which we couldn't fix)");

			}
		}else if(args[0].equals("2")){
			try{
				Float bucketSize = null;
				if(args.length == 3){
					bucketSize = null;
					points = CSVParser.CSVToPoint(args[1]);
					compare = CSVParser.CSVToPointCompare(args[2]);
				}else{
					bucketSize = Float.valueOf(args[3]);
					points = CSVParser.CSVToPoint(args[1]);
					compare = CSVParser.CSVToPointCompare(args[2]);
				}
				KMeans.lsh(points, bucketSize, (float)7.0);
				ArrayList<Integer> result2 = KMeans.PointsToIntegerList(points);
				KMeans.NMI(compare, result2);
			}catch (Exception e) {
				e.printStackTrace();
				 System.err.println("An exception was thrown (most of the time for an unknow reason which we couldn't fix)");

			}
			
		}else if(args[0].equals("3")){
			try{
				Float bucketSize = null;
				if(args.length == 3){
					bucketSize = null;
					points = CSVParser.CSVToPoint(args[1]);
					compare = CSVParser.CSVToPointCompare(args[2]);
				}else{
					bucketSize = Float.valueOf(args[3]);
					points = CSVParser.CSVToPoint(args[1]);
					compare = CSVParser.CSVToPointCompare(args[2]);
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
