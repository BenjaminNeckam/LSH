package at.ac.univie.main;

import java.util.ArrayList;

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

	public static void main(String[] args){
		ArrayList<Point> points = CSVParser.CSVToPoint("C:/Users/sipic/Desktop/SDM_csv/LSH-nmi-corrected.csv");
		
		KMeans k = new KMeans();
		
		ArrayList a = k.initCentroids(points,6);
		
		for(int i=0;i<a.size();++i){
			System.out.println(a.get(i));
		}
		
		
		Plot scatterplotdemo4 = new Plot("K-Means Start",a, 6);
		scatterplotdemo4.pack();
		RefineryUtilities.centerFrameOnScreen(scatterplotdemo4);
		scatterplotdemo4.setVisible(true);

	
		

	}
}
