package at.ac.unive.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class CSVParser {
	
	public static ArrayList<Point> CSVToPoint(String file){
		ArrayList<Point> points = new ArrayList<>();
		ArrayList<Float> coordinates;
		String line="";
		String splitter = ",";
		try(BufferedReader br = new BufferedReader(new FileReader(file))){
			while((line = br.readLine()) != null){
				coordinates = new ArrayList<>();
				String[] StringCoordinates = line.split(splitter);
				for(int i=0;i<StringCoordinates.length;i++){
					coordinates.add(Float.valueOf(StringCoordinates[i]));
				}
				points.add(new Point(coordinates));
			}
		} catch (IOException e){
			e.printStackTrace();
		}
		return points;
	}

}
