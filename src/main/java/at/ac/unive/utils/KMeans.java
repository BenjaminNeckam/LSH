package at.ac.unive.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class KMeans {
	
	public static void hashPoints(ArrayList<Point> points){
		Point v = generateHashVector();
		//Aktuelle Bucketsize
		int w = 1;
		ArrayList<Bucket> buckets = new ArrayList<>();
		//First hash
		for(Point point: points){
			hash(point, v);
		}
		
	}
	
	/**
	 * Function to compare results from Martin Perdacher
	 * @param one
	 * @param two
	 * @return
	 */
	public static double NMI(ArrayList<Integer> one, ArrayList<Integer> two){
		if(one.size()!=two.size()){
			throw new IllegalArgumentException("Sizes don't match!");
		}
		int maxone = Collections.max(one);
		int maxtwo = Collections.max(two);

		double[][] count = new double[maxone+1][maxtwo+1];
		//System.out.println(count[1][2]);
		for(int i=0;i<one.size();i++){
			count[one.get(i)][two.get(i)]++;
		}
		//i<maxone=R
		//j<maxtwo=C
		double[] bj = new double[maxtwo+1];
		double[] ai = new double[maxone+1];

		for(int m=0;m<(maxtwo+1);m++){
			for(int l=0;l<(maxone+1);l++){
				bj[m]=bj[m]+count[l][m];
			}
		}
		for(int m=0;m<(maxone+1);m++){
			for(int l=0;l<(maxtwo+1);l++){
				ai[m]=ai[m]+count[m][l];
			}
		}

		double N=0;
		for(int i=0;i<ai.length;i++){
			N=N+ai[i];
		}
		double HU = 0;
		for(int l=0;l<ai.length;l++){
			double c=0;
			c=(ai[l]/N);
			if(c>0){
				HU=HU-c*Math.log(c);
			}
		}

		double HV = 0;
		for(int l=0;l<bj.length;l++){
			double c=0;
			c=(bj[l]/N);
			if(c>0){
				HV=HV-c*Math.log(c);
			}
		}
		double HUstrichV=0;
		for(int i=0;i<(maxone+1);i++){
			for(int j=0;j<(maxtwo+1);j++){
				if(count[i][j]>0){
					HUstrichV=HUstrichV-count[i][j]/N*Math.log(((count[i][j])/(bj[j])));
				}
			}
		}
		double IUV = HU-HUstrichV;
		double reto = IUV/(Math.max(HU, HV));

		System.out.println("NMI: "+reto);
		return reto;
	}
	
	/**
	 * Normaldistributed vector for hashfunction
	 * @return
	 */
	public static Point generateHashVector(){
		Point v;
		ArrayList<Float> coordinates = new ArrayList<>();
		Random random = new Random();
		float temp;
		//10 weil Dim 10 im Moment
		for(int i=0; i<10;i++) {
			temp = (float) (0+random.nextGaussian());
			coordinates.add(temp);
		}
		v = new Point(coordinates);
		return v;
	}
	
	public static void hash(Point point, Point hashVector){
		float sum = 0;
		int size = point.getCoordinates().size();
		for(int i=0; i<size; i++){
			sum += point.getCoordinates().get(i)*hashVector.getCoordinates().get(i);
		}
		point.setHashValue(sum);
	}

}
