package at.ac.unive.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class CSVParser {

	/**
	 * Parses CSV-File in ArrayList of Points
	 * 
	 * @param file
	 * @return
	 */
	public static ArrayList<Point> CSVToPoint(String file) {
		ArrayList<Point> points = new ArrayList<Point>();
		ArrayList<Float> coordinates;
		String line = "";
		String splitter = ",";
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			while ((line = br.readLine()) != null) {
				coordinates = new ArrayList<>();
				String[] StringCoordinates = line.split(splitter);
				for (int i = 0; i < StringCoordinates.length; i++) {
					coordinates.add(Float.valueOf(StringCoordinates[i]));
				}
				points.add(new Point(coordinates));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return points;
	}
	
	public static ArrayList<Integer> CSVToPointCompare(String file) {
		ArrayList<Integer> points = new ArrayList<>();
		float tmp;
		Integer compareClass;
		String line = "";
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			while ((line = br.readLine()) != null) {
				tmp = Float.valueOf(line);
				compareClass = Math.round(tmp);
				points.add(compareClass);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		return points;
	}
}
