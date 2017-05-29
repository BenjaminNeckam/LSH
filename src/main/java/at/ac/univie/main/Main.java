package at.ac.univie.main;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import at.ac.unive.utils.CSVParser;
import at.ac.unive.utils.KMeans;
import at.ac.unive.utils.Point;

public class Main {

	public static void main(String[] args){
		ArrayList<Point> points = CSVParser.CSVToPoint("C:/Users/Benni/git/LSH/src/test/resources/LSH-nmi-corrected.csv");
		try {
			KMeans.hashPoints(points);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
//		ArrayList<Float> point = new ArrayList<>();
//		Random random = new Random();
//		for(int i=0;i<50;i++){
//			point.add((-100 + (0+100)*random.nextFloat()));
//		}
//		Collections.sort(point, new Comparator<Float>() {
//            public int compare(Float o1, Float o2) {
//            	if (o1 > o2) {
//                    return -1;
//                  }
//                  if (o1 < o2) {
//                    return 1;
//                  }
//                  return 0;
//            }
//        });
//		
//		for(Float value:point){
//			System.out.println(value);
//		}
//		for(Float value:point){
//			System.out.println(Float.hashCode(value));
//		}
	}
}
