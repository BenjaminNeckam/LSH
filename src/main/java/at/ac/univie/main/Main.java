package at.ac.univie.main;

import java.util.ArrayList;

import at.ac.unive.utils.CSVParser;
import at.ac.unive.utils.KMeans;
import at.ac.unive.utils.Point;

public class Main {

	public static void main(String[] args){
		ArrayList<Point> points = CSVParser.CSVToPoint("C:/Users/Benni/git/LSH/src/test/resources/LSH-nmi-corrected.csv");
//		KMeans.hashPoints(points);

	}
}
