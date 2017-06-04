package at.ac.unive.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;



import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

public class KMeans {

	public static void lsh(ArrayList<Point> points,Float bucketSize) throws Exception{
		long startTime = System.nanoTime();
		long startTimeV = System.nanoTime();
		Point v = generateHashVector();
		double estimatedTimeV = (System.nanoTime() - startTimeV)/ 1000000000.0;
		System.out.println("\nElapsed Time Generating HashVector: " + estimatedTimeV + " seconds");
		
		long startTimeB = System.nanoTime();
		//TODO more than one hash
		TreeMap<Object,Bucket> buckets = initBucketsAndHashPoints(points, v, bucketSize);
		double estimatedTimeB = (System.nanoTime() - startTimeB)/ 1000000000.0;
		System.out.println("\nElapsed Time Init Buckets and Hash Points: " + estimatedTimeB + " seconds");
		
		long startTimeS = System.nanoTime();
		//Returns a sorted HashMap by Values (size of pointslist)
		HashMap<Object,Bucket> sortedBuckets = sortByValues(buckets);
		Set<Object> sortedKeySet= sortedBuckets.keySet();
		Object[] keys = sortedKeySet.toArray();
		double estimatedTimeS = (System.nanoTime() - startTimeS)/ 1000000000.0;
		System.out.println("\nElapsed Time Sorting: " + estimatedTimeS + " seconds");

		long startTimeC = System.nanoTime();
		//Init Centroids (calcultateCentroidOfBucket just for first centroid calculation)
		ArrayList<Centroid> centroid = new ArrayList<>();
		Centroid tmpCentroid;
		for(int i=0;i<15;i++){
			tmpCentroid =  new Centroid();
			buckets.get(keys[i]).setClusterNumb(i);
			tmpCentroid.addBucket(buckets.get(keys[i]));
			tmpCentroid.setBucketHashCode(buckets.get(keys[i]).getBucketHashCode());
			tmpCentroid.setClusterNumb(i);
			tmpCentroid.updateCentroid();
			centroid.add(tmpCentroid);
			tmpCentroid=null;
	    }

//		System.out.println("Centroids before update: ");
//		for(Centroid c:centroid) {
//			System.out.println(c.toString());
//			System.out.println(c.getBuckets().size());
//		}

		//Calculate remaining buckets
		for(int i=15;i<keys.length;i++) {
			if(buckets.get(keys[i]).getPoints() != null || buckets.get(keys[i]).getPoints().size() > 0) {
				if(buckets.get(keys[i]).getPoints().size()!=0){
					computeDistAndAssignCenter(buckets.get(keys[i]), centroid);
				}
				
			}
		}



//		System.out.println("Centroids after update: ");
		for(Centroid c:centroid) {
			c.updateCentroid();
//			System.out.println(c.toString());
//			System.out.println(c.getBuckets().size());
		}
		double estimatedTimeC = (System.nanoTime() - startTimeC)/ 1000000000.0;
		System.out.println("\nElapsed Time Init Centroids: " + estimatedTimeC + " seconds");
		
		
//		//TEST
//		ArrayList<Point> centroid2 = initCentroids(points, 15);
//		float hashValue = 0;
//		for(Point point:centroid2){
//			hashValue = point.getHashValue();
//		
////				buckets.get(hashCode).addPoint(point);
//				if(Float.hashCode(hashValue)>0){
//					buckets.floorEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
//				}else{
//					buckets.ceilingEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
//				}
//				System.out.println("Centroids2 hashBucketCode: " + point.getBucketHashCode());
//		}
//		
//		for(Point point:centroid2){
//			for(Point point2:centroid2){
//				if(!point.equals(point2)){
//					if(point.getBucketHashCode() == point2.getBucketHashCode()){
//						System.out.println("SAME BUCKET");
//					}
//				}
//					
//			}
//		}
//		
		//TEST ENDE

	    double estimatedTime = (System.nanoTime() - startTime)/ 1000000000.0;
		System.out.println("\nElapsed Time: " + estimatedTime + " seconds");

	}
	
	public static TreeMap<Object,Bucket> initBucketsAndHashPoints(ArrayList<Point> points,Point v, Float w) throws Exception{
		//Current Bucketsize
				TreeMap<Object, Bucket> buckets = new TreeMap<>();
				float min = Float.MAX_VALUE;
				float max = Float.MIN_VALUE;
				float hashValue = 0;
				//First hash and computing Bucketsize

				for(Point point: points){
					hash(point, v, 4);
					if(point.getHashValue() >= max){
						max = point.getHashValue();
					}
					if(point.getHashValue() <= min){
						min = point.getHashValue();
					}
				}

				//Set intervals of buckets
				//TODO Optimize?
				if(w==null){
					w = (float)(Math.abs(max-min)/(Math.sqrt(points.size())/Math.log(points.size())));
				}
				System.out.println("Bucketsize: " + w);
				System.out.println("MIN: " + min);
				System.out.println("MAX: " + max);
				int bucketNumbs = (int) Math.ceil(Math.sqrt(points.size())/Math.log(points.size()));
				System.out.println("Bucketnumbers: " + bucketNumbs);
				float interval = min;
				for(int i=0;i<bucketNumbs;i++){
					buckets.put(Float.hashCode(interval),new Bucket(new float[]{interval,(float) (interval+(w-0.0001))},Float.hashCode(interval)));
//					System.out.println("IntevalHash: " + Float.hashCode(intervall));
					interval+=w;
				}
				
				//Set points to buckets
				for(Point point:points){
					hashValue = point.getHashValue();
					if(hashValue < min || hashValue > max){
						throw new Exception("Hashcode is not in range");
					}else{
//						buckets.get(hashCode).addPoint(point);
						if(Float.hashCode(hashValue)>0){
							buckets.floorEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
						}else{
							buckets.ceilingEntry(Float.hashCode(hashValue)).getValue().addPoint(point);
						}

					}
				}
				return buckets;

	}
	
	public static ArrayList<Point> initCentroids(ArrayList<Point> points,int numOfCtr){
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
	
	public static int nearestCentroid(ArrayList<Point> centroids,Point p){
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

		for(int m=0;m<maxtwo+1;m++){
			for(int l=0;l<maxone+1;l++){
				bj[m]=bj[m]+count[l][m];
			}
		}
		for(int m=0;m<maxone+1;m++){
			for(int l=0;l<maxtwo+1;l++){
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
			c=ai[l]/N;
			if(c>0){
				HU=HU-c*Math.log(c);
			}
		}

		double HV = 0;
		for(int l=0;l<bj.length;l++){
			double c=0;
			c=bj[l]/N;
			if(c>0){
				HV=HV-c*Math.log(c);
			}
		}
		double HUstrichV=0;
		for(int i=0;i<maxone+1;i++){
			for(int j=0;j<maxtwo+1;j++){
				if(count[i][j]>0){
					HUstrichV=HUstrichV-count[i][j]/N*Math.log(count[i][j]/bj[j]);
				}
			}
		}
		double IUV = HU-HUstrichV;
		double reto = IUV/Math.max(HU, HV);

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

	public static void hash(Point point, Point hashVector, int w){
		float sum = 0;
		int size = point.getCoordinates().size();
		for(int i=0; i<size; i++){
			sum += point.getCoordinates().get(i)*hashVector.getCoordinates().get(i);
		}
		Random r= new Random(2);
		double bias = r.nextDouble();
		point.setHashValue((float)(sum+bias)/w);
	
	
	}

	@SuppressWarnings("unchecked")
	private static HashMap<Object,Bucket> sortByValues(TreeMap<Object,Bucket> map) {
		List list = new LinkedList(map.entrySet());
		Collections.sort(list, new Comparator<Object>() {

			@Override
			public int compare(Object o1, Object o2) {
				Map.Entry<Object, Bucket> entry1 = (Entry<Object, Bucket>) o1;
				Map.Entry<Object, Bucket> entry2 = (Entry<Object, Bucket>) o2;
				if (entry1.getValue().getPoints().size() > entry2.getValue().getPoints().size()) {
                  return -1;
                }
                if (entry1.getValue().getPoints().size() < entry2.getValue().getPoints().size()) {
                  return 1;
                }
                return 0;
				}
			});


		HashMap<Object,Bucket> sortedTree = new LinkedHashMap();
		for(Object obj: list) {
			Map.Entry<Object, Bucket> entry = (Entry<Object, Bucket>) obj;
			sortedTree.put(entry.getKey(), entry.getValue());
		}
		return sortedTree;
	}

	public static float dist(Point point1, Point point2){
		float sum=0;
		for(int i=0;i<point1.getCoordinates().size();i++){
			sum+=Math.pow(point1.getCoordinates().get(i)-point2.getCoordinates().get(i), 2);
		}
		float norm=(float)Math.sqrt(sum);
		return norm;
	}

	//At the moment radical usage: compute similarity of a single or few random point to all center and assign all points of the bucket to the center
	public static void computeDistAndAssignCenter(Bucket bucket,ArrayList<Centroid> centroid) {
		ArrayList<Point> bucketPoints = bucket.getPoints();
			float tmp = dist(bucketPoints.get(0),centroid.get(0));
			float tmpDist=0;
			//Just 10% of the points will get checked
			int size = (int) Math.ceil(bucketPoints.size()*0.1);
			for(int i=0;i<size;i++) {
				for(int j=0;j<centroid.size();j++) {
					tmpDist= dist(bucketPoints.get(i),centroid.get(j));
					if(tmpDist<tmp) {
						tmp=tmpDist;
						bucket.setClusterNumb(centroid.get(j).getClusterNumb());
					}
				}
			}
			centroid.get(bucket.getClusterNumb()).addBucket(bucket);

	}
	
	public static ArrayList<Bucket> initBuckets(ArrayList<Point> listOfPoints,int numOfBuckets, int w){
		
		Point v = generateHashVector();
		ArrayList<Float> hashedP = new ArrayList<>();
		ArrayList<Bucket> buckets = new ArrayList<>();
		
		
		
		for(Point p:listOfPoints){
			hash(p, v, w);
			hashedP.add(p.getHashValue());
		}

		float hMin = hashedP.get(0);
		float hMax = hashedP.get(0);
		
		
		for(int i=0;i<hashedP.size();++i){
			if(hashedP.get(i)<hMin) hMin = hashedP.get(i);
			if(hashedP.get(i)>hMax) hMax = hashedP.get(i);
		}
		
		for(float i=hMin;i<hMax;){
			float[] f = {i,i+((hMax-hMin)/numOfBuckets)};
			int hk = (int)(f[1] - f[0]);
			Bucket b = new Bucket(f,hk);
			buckets.add(b);
			i = (float) (((i+((hMax-hMin)/numOfBuckets))<hMax) ? i+((hMax-hMin)/numOfBuckets):hMax);
		}
		
		
		for(Point p:listOfPoints){
			for(int i=0;i<numOfBuckets;++i){
				float[] interval = buckets.get(i).getBucketInterval();
				if(p.getHashValue() > interval[0] && p.getHashValue() < interval[1]){
					buckets.get(i).addPoint(p);
					p.setClusterNumb(i);
				}
			}
		}
		
		return buckets;
	}
	
	public static void lloyd(ArrayList points,int numbOfBuckets){
		long startTime = System.nanoTime();
		long startTimeV = System.nanoTime();
		ArrayList<Bucket> b = new ArrayList();
		b = initBuckets(points, numbOfBuckets, 11);
		for(int i=0;i<b.size();++i){
			System.out.println(b.get(i).getPoints().size());
		}
		
		double estimatedTimeV = (System.nanoTime() - startTimeV)/ 1000000000.0;
		System.out.println("\nElapsed Time Bucket creation: " + estimatedTimeV + " seconds");
		
		
		Point centroid = getCentroid(b.get(0).getPoints());
		Point centoridtmp = null;
		int counter = 0;
		while(counter!=b.size()){
			if((counter+1) < b.size()){
				float maxDist = dist(b.get(counter).getPoints().get(0), centroid);
				for(int i=0;i<b.get(counter).getPoints().size();++i){
					float distance = dist(b.get(counter).getPoints().get(i), centroid);
					if(distance > maxDist) maxDist = distance;
				}
				
				estimatedTimeV = (System.nanoTime() - startTimeV)/ 1000000000.0;
				System.out.println("\nElapsed Time maxDist calculation " + estimatedTimeV + " seconds");
				
				
				do{
			
					for(int i=0;i<b.get(counter+1).getPoints().size();++i){
						float distance = dist(b.get(counter+1).getPoints().get(0), centroid);
						if(distance < maxDist){
							b.get(counter+1).getPoints().get(i).setClusterNumb(counter);
							b.get(counter).getPoints().add(b.get(counter+1).getPoints().get(i));
							b.get(counter+1).getPoints().remove(b.get(counter+1).getPoints().get(i));
						}
					}
					centoridtmp = centroid;
					centroid = getCentroid(b.get(counter).getPoints());
					
					estimatedTimeV = (System.nanoTime() - startTimeV)/ 1000000000.0;
					System.out.println("\nElapsed Time do while " + estimatedTimeV + " seconds");
			
				}while(centroid.equals(centoridtmp));
				
				counter++;
				System.out.println(counter);
			
			}else{
				do{
					
					centoridtmp = centroid;
					centroid = getCentroid(b.get(counter).getPoints());
					
			
				}while(centroid.equals(centoridtmp));
				counter++;
			}
		}
		for(int i=0;i<b.size();++i){
			System.out.println(b.get(i).getPoints().size());
		}
		
	}
	
	
	public static Point getCentroid(ArrayList<Point> pointsInBucket){
		ArrayList<Float> coor = new ArrayList();
		Point p = null;
		int counter = 0;
		float sum = 0;
		
		while(counter!=pointsInBucket.get(0).getCoordinates().size()){
			for(Point point:pointsInBucket){
				sum+=point.getCoordinates().get(counter);
			}
			coor.add(sum);
		counter++;
		}
		p = new Point(coor);
		
		float minDist = dist(p, pointsInBucket.get(0));
		int index = 0;
		
		for(Point point:pointsInBucket){
			float distance = dist(p, point);
			if(distance < minDist) index = pointsInBucket.indexOf(point);
		}
		
		return pointsInBucket.get(index);
	}

}
