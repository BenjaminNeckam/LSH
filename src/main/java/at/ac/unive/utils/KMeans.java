package at.ac.unive.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableSet;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;



public class KMeans {
	
	public static void hashPoints(ArrayList<Point> points) throws Exception{
		long startTime = System.nanoTime();
		Point v = generateHashVector();
		//Current Bucketsize 
		float w;
		TreeMap<Object, Bucket> buckets = new TreeMap<>();
		float min = Float.MAX_VALUE;
		float max = Float.MIN_VALUE;
		//First hash and computing Bucketsize


		for(Point point: points){
			hash(point, v);
			if(point.getHashValue() >= max){
				max = point.getHashValue();
			}
			if(point.getHashValue() <= min){
				min = point.getHashValue();
			}
		}
		//Set intervalls of buckets
		//TODO Optimize?
		w = (float)((Math.abs(max-min))/(Math.sqrt(points.size()/Math.log(points.size()))));
		System.out.println("MIN: " + min);
		System.out.println("MAX: " + max);
		int bucketNumbs = (int) Math.ceil((Math.sqrt(points.size()/Math.log(points.size()))));
		System.out.println("Bucketnumbers: " + bucketNumbs);
		float intervall = min;
		for(int i=0;i<bucketNumbs;i++){
			buckets.put(Float.hashCode(intervall),new Bucket(new float[]{intervall,(float) (intervall+(w-0.0001))}));
			intervall+=w;
		}
		
		//Set points to buckets
		float hashValue = 0;
		
		for(Point point:points){
			hashValue = point.getHashValue();
			if(hashValue < min || hashValue > max){
				throw new Exception("Hashcode is not in range");
			}else{
//				buckets.get(hashCode).addPoint(point);
				if(Float.hashCode(hashValue)>0){
					buckets.floorEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
				}else{
					buckets.ceilingEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
				}
				
			}
		}
		
		//Swap Key Value Pairs to get k-Buckets with most elements in it (Our case 15)
		TreeMap<Integer,Object> rev = new TreeMap<>();
	    for(Map.Entry<Object,Bucket> entry : buckets.entrySet()){
	    	rev.put(entry.getValue().getPoints().size(), entry.getKey());
	    }
	    NavigableSet<Integer> orderedValues = rev.descendingKeySet();
	    double estimatedTime = (System.nanoTime() - startTime)/ 1000000000.0;
		System.out.println("\nElapsed Time: " + estimatedTime + " seconds");
		
	    for(int i=0;i<15;i++){
	    	//TODO Get 15 biggest buckets, calculate Centroid of it
	    }
	
	}
	
	public ArrayList<Point> initCentroids(ArrayList<Point> points,int numOfCtr){
		ArrayList<Point> pointArray = new ArrayList(points);
		ArrayList<Point> centroids = new ArrayList();
		Random r = new Random();
		
		Point firstCentroid = pointArray.get(r.nextInt(pointArray.size()));
		centroids.add(firstCentroid);
		pointArray.remove(firstCentroid);
		
		
		double[] d2 = new double[pointArray.size()];
		
		while(centroids.size()<numOfCtr){
			double sum = 0;
			for(int i=0;i<centroids.size();++i){
				Point p = pointArray.get(i);
				int index = nearestCentroid(centroids, p);
				double d = dist(p, centroids.get(index));
				sum+=d*d;
				d2[i]=sum;
			}
			
			double d = r.nextDouble() * sum;
			for(int i=0;i<d2.length;++i){
				if(d2[i]>=d){
					Point p = pointArray.get(i);
					centroids.add(p);
					pointArray.remove(p);
					break;
				}
			}
		}
		
		return centroids;
		
	}
	
	public int nearestCentroid(ArrayList<Point> centroids,Point p){
		double minDistance = Double.MAX_VALUE;
        int indexOfCentroid = 0;
        
        for (int i=0;i<centroids.size();++i){
        	float distance = dist(p, centroids.get(i));
        	if (distance < minDistance) {
                minDistance = distance;
                indexOfCentroid = i;
            }
        }
        
        return indexOfCentroid;
	}
	
	
	public static float dist(Point point1, Point point2){
		float sum=0;
		for(int i=0;i<point1.getCoordinates().size();i++){
			sum+=Math.pow(point1.getCoordinates().get(i)-point2.getCoordinates().get(i), 2);
		}
		float norm=(float)Math.sqrt(sum);
		return norm;

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
		ArrayList<Float> coordinates = new ArrayList<Float>();
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
