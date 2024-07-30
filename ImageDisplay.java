import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import java.util.Arrays;
import javax.swing.*;
import java.util.HashMap;
import java.util.Map;

import java.util.ArrayList;
import java.util.List;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map.Entry;
import java.util.LinkedList;
import java.io.File;

class Pair {
    double x, y;

    public Pair(double x, double y) {
        this.x = x;
        this.y = y;
    }
}

class Cluster {
    int id;
    List<Pair> pixels = new ArrayList<>();
    public void clear() {
        this.pixels.clear();
        // You can add more clearing logic if needed
    }
}

public class ImageDisplay {

    JFrame frame;
    JLabel lbIm1;
    BufferedImage imgOne;
    BufferedImage imgTwo;

    int width = 640; // default image width and height
    int height = 480;

    Map<Pair, Double> fh = new HashMap<>();
    Map<Pair, Double> fs = new HashMap<>();
    Map<Pair, Double> fv = new HashMap<>();
    Map<Pair, Double> inputfh = new HashMap<>();
    Map<Pair, Double> inputfs = new HashMap<>();
    Map<Pair, Double> inputfv = new HashMap<>();
    Map<Double, Double> hueFrequency = new HashMap<>();
    Map<Double, Double> inputHueFrequency = new HashMap<>();
    // Map<Double, Double> filteredMap = new HashMap<>();
    double[] histogram;
    double[] saturationHistogram;
    double[] valueHistogram;

    // double thresholdFrequency = 0.0001;
    
    
        
    private static final int MIN_CLUSTER_SIZE = 30;



    /** Read Image RGB
     *  Reads the image of given width and height at the given imgPath into the provided BufferedImage.
     */
    private void readImageRGB(int width, int height, String imgPath, BufferedImage img)
    {
        try
        {
            int frameLength = width*height*3;

            File file = new File(imgPath);
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            raf.seek(0);

            long len = frameLength;
            byte[] bytes = new byte[(int) len];

            raf.read(bytes);

            int ind = 0;
            for(int y = 0; y < height; y++)
            {
                for(int x = 0; x < width; x++)
                {
                    byte a = 0;
                    byte r = bytes[ind];
                    byte g = bytes[ind+height*width];
                    byte b = bytes[ind+height*width*2]; 

                    int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                    img.setRGB(x,y,pix);
                    ind++;
                }
            }
        }
        catch (FileNotFoundException e) 
        {
            e.printStackTrace();
        } 
        catch (IOException e) 
        {
            e.printStackTrace();
        }
    }

    public static void generateHueHistogramToCSV(double[] hueValues, int numBins, String outputPath) {
        double[] histogram = new double[numBins];
        double binSize = 360.0 / numBins;
    
        for (double hue : hueValues) {
            // Normalize hue to [0, 360)
            hue = (hue + 360) % 360;
    
            // Calculate the bin index for the hue
            int binIndex = (int) (hue / binSize);
    
            // Increment the corresponding bin in the histogram
            histogram[binIndex]++;
        }
    
        // Normalize the histogram by dividing each bin count by the total number of hue values
        int totalHueValues = hueValues.length;
        for (int i = 0; i < numBins; i++) {
            histogram[i] /= totalHueValues;
        }
    
        // Write the histogram data to a CSV file
        try (FileWriter writer = new FileWriter(outputPath)) {
            // Write CSV header
            writer.append("Bin,Value\n");
    
            // Write histogram data
            for (int i = 0; i < numBins; i++) {
                writer.append(i + "," + histogram[i] + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void saveHistogramToCSV(Map<Pair, Double> histogram, String filePath) {
        try (FileWriter writer = new FileWriter(filePath)) {
            // Write CSV header
            writer.append("x,y,value\n");

            // Calculate the total number of pixels
            int totalPixels = width * height;

            // Write normalized histogram values
            for (Entry<Pair, Double> entry : histogram.entrySet()) {
                Pair p = entry.getKey();
                double value = entry.getValue() / totalPixels; // Normalize by dividing by the total pixels
                writer.append(p.x + "," + p.y + "," + value + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static double[][] extractTopBins(double[] histogram, int topCount) {
        // Clone the original histogram to avoid modifying it
        double[] sortedHistogram = Arrays.copyOf(histogram, histogram.length);
    
        // Sort the cloned histogram in descending order
        Arrays.sort(sortedHistogram);
        for (int i = 0; i < sortedHistogram.length / 2; i++) {
            double temp = sortedHistogram[i];
            sortedHistogram[i] = sortedHistogram[sortedHistogram.length - 1 - i];
            sortedHistogram[sortedHistogram.length - 1 - i] = temp;
        }
    
        // Extract the top 'topCount' values along with their bin indices
        double[][] topValuesWithIndices = new double[topCount][2];
        for (int i = 0; i < topCount; i++) {
            double topValue = sortedHistogram[i];
            int binIndex = -1;
    
            // Find the bin index corresponding to the top value
            for (int j = 0; j < histogram.length; j++) {
                if (histogram[j] == topValue) {
                    binIndex = j;
                    break;
                }
            }
    
            topValuesWithIndices[i][0] = topValue;
            topValuesWithIndices[i][1] = binIndex;
        }
    
        return topValuesWithIndices;
    }
    

    public static double[] generateSaturationHistogram(double[] saturationValues, int numBins) {
        double[] saturhistogram = new double[numBins];
        double binSize = 1.0 / numBins; // Saturation values are typically in the range [0, 1]
    
        for (double saturation : saturationValues) {
            // Calculate the bin index for the saturation
            int binIndex = (int) (saturation / binSize);
    
            // Ensure that binIndex stays within the valid range
            binIndex = Math.min(numBins - 1, binIndex);
    
            // Increment the corresponding bin in the histogram
            saturhistogram[binIndex]++;
        }
    
        // Normalize the histogram by dividing each bin count by the total number of saturation values
        int totalSaturationValues = saturationValues.length;
        for (int i = 0; i < numBins; i++) {
            saturhistogram[i] /= totalSaturationValues;
        }
    
        return saturhistogram;
    }
    
    public static double[] generateValueHistogram(double[] Values, int numBins) {
        double[] valuehistogram = new double[numBins];
        double binSize = 1.0 / numBins; // Saturation values are typically in the range [0, 1]
    
        for (double val : Values) {
            // Calculate the bin index for the saturation
            int binIndex = (int) (val / binSize);
    
            // Ensure that binIndex stays within the valid range
            binIndex = Math.min(numBins - 1, binIndex);
    
            // Increment the corresponding bin in the histogram
            valuehistogram[binIndex]++;
        }
    
        // Normalize the histogram by dividing each bin count by the total number of saturation values
        int totalValues = Values.length;
        for (int i = 0; i < numBins; i++) {
            valuehistogram[i] /= totalValues;
        }
    
        return valuehistogram;
    }
    
    public static double[] generateHueHistogram(double[] hueValues, int numBins) {
        double[] histogram = new double[numBins];
        double binSize = 360.0 / numBins;

        for (double hue : hueValues) {
            // Normalize hue to [0, 360)
            hue = (hue + 360) % 360;

            // Calculate the bin index for the hue
            int binIndex = (int) (hue / binSize);

            // Increment the corresponding bin in the histogram
            histogram[binIndex]++;
        }

        // Normalize the histogram by dividing each bin count by the total number of hue values
        int totalHueValues = hueValues.length;
        for (int i = 0; i < numBins; i++) {
            histogram[i] /= totalHueValues;
        }

        return histogram;
    }  // System.out.println(Arrays.toString(histogram));

    private void saveHueFrequencyToCSV(Map<Double, Double> hueFrequency, String filePath) {
        try (FileWriter writer = new FileWriter(filePath)) {
            // Write CSV header
            writer.append("hue,frequency\n");
    
            // Write hue frequencies
            for (Entry<Double, Double> entry : hueFrequency.entrySet()) {
                writer.append(entry.getKey() + "," + entry.getValue() + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



    public void showIms(String[] args){

        if (args.length < 2) {
            System.out.println("Please provide two image file paths as arguments.");
            return;
        }
        
        inputImage(args[0]);

        // Use label to display the image
        this.frame = new JFrame();
        GridBagLayout gLayout = new GridBagLayout();
        this.frame.getContentPane().setLayout(gLayout);
        // lbIm1 = new JLabel(new ImageIcon(imgOne));
        lbIm1 = new JLabel(new ImageIcon(this.imgTwo));
        
        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 0.5;
        c.gridx = 0;
        c.gridy = 0;

        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridx = 0;
        c.gridy = 1;
        this.frame.getContentPane().add(lbIm1, c);

        this.frame.pack();
        
        this.frame.setVisible(true);
        this.frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
       
        // objectImage(args[1]);
        for (int i = 1; i < args.length; i++) {
            fh.clear();
            fs.clear();
            fv.clear();
            // histogram.clear();
            // saturationHistogram.clear();
            // valueHistogram.clear();
            objectImage(args[i]);
        }
        // frame.setResizable(false);

    }


    private void objectImage(String imagePath) {
        
        double h = 0;
        double s = 0;
        double v = 0;
        double cmax, cmin, diff;

        imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        readImageRGB(width, height, imagePath, imgOne);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {

                int rgbb = imgOne.getRGB(x, y);
                double r = (rgbb >> 16) & 0xFF;
                double g = (rgbb >> 8) & 0xFF;
                double b = rgbb & 0xFF;

                r = r / 255.0;
                g = g / 255.0;
                b = b / 255.0;

                // if (r != 0.0 && g != 1.0 && b != 0.0 && (r != g && g != b) ) { //not green 
                if (r != 0.0 && g != 1.0 && b != 0.0) { //not green 

                    cmax = Math.max(r, Math.max(g, b));
                    cmin = Math.min(r, Math.min(g, b));
                    diff = cmax - cmin;

                    if (diff > 0) {
                        if (cmax == cmin)
                            h = 0;
                        else if (cmax == r)
                            h = (60 * ((g - b) / diff) + 360) % 360;
                        else if (cmax == g)
                            h = (60 * ((b - r) / diff) + 120) % 360;
                        else if (cmax == b)
                            h = (60 * ((r - g) / diff) + 240) % 360;

                        if (cmax == 0)
                            s = 0;
                        else
                            s = (diff / cmax);
                    } else {
                        s = 0;
                    }
                    v = cmax;

                    if (h < 0) {
                        h = 360 + h;
                    }

                    Pair p = new Pair(x, y);

                    fh.put(p, h);
                    fs.put(p, s);
                    fv.put(p, v);
                    // hueFrequency.put(h, hueFrequency.getOrDefault(h, 0.0) + 1.0);
                    
                }
            }
        }


        System.out.println("Fh size " + fh.size());

        //Hue values array
        double[] hueValues = new double[fh.size()]; 
        int index = 0;
        for (double hue : fh.values()) {
            // Fill the array with hue values from the fh map
            hueValues[index++] = hue; 
        }
        int numBins = 36; // Number of bins for the histogram
        histogram = generateHueHistogram(hueValues, numBins);
        // String outputPath = "histogramobjectnew.csv";
        // generateHueHistogramToCSV(hueValues, numBins, outputPath);

        // Extract the saturation values from the inputfs map
        double[] saturationValues = new double[fs.size()]; 
        int satuindex=0;
        for (double saturation : fs.values()) {
            // Fill the array with hue values from the fh map
            saturationValues[satuindex++] = saturation; 
        }
        int saturationnumBins = 10;
        saturationHistogram = generateSaturationHistogram(saturationValues, saturationnumBins);

        double[] Vals = new double[fv.size()]; 

        int val=0;
        for (double valu : fv.values()) {
            // Fill the array with hue values from the fh map
            Vals[val++] = valu; 
        }
        int valnumBins = 10;
        valueHistogram = generateValueHistogram(Vals, valnumBins);
        // for (int i = 0; i < valueHistogram.length; i++) {
        //     System.out.printf("val Bin %d: %.4f%n", i, valueHistogram[i]);
        // }

        
        //frequency map creation starts here
        // for (double hue : fh.values()) {
        //     hueFrequency.put(hue, hueFrequency.getOrDefault(hue, 0.0) + 1.0);
        // }

        // Normalize the hue frequencies
        // double fhSize = fh.size();
        // for (Map.Entry<Double, Double> entry : hueFrequency.entrySet()) {
        //     double normalizedFrequency = (double) entry.getValue() / fhSize;
        //     hueFrequency.put(entry.getKey(), normalizedFrequency);
        // }

        // double maxHueValue = -1.0;
        // for (Double hueValue : hueFrequency.keySet()) {
        //     if (hueValue > maxHueValue) {
        //         maxHueValue = hueValue;
        //     }
        // }
        // System.out.println("in obj  "+maxHueValue);
        
        // double thresholdFrequency =  0.1 * maxHueValue ;

        // for (Map.Entry<Double, Double> entry : hueFrequency.entrySet()) {
        //     if (entry.getValue() > thresholdFrequency) {
        //         filteredMap.put(entry.getKey(), entry.getValue());
        //     }
        // }
 
        

        // saveHistogramToCSV(fh, "hue_histogram.csv");
        // saveHistogramToCSV(fs, "saturation_histogram.csv");
        // saveHistogramToCSV(fv, "value_histogram.csv");
        // saveHueFrequencyToCSV(hueFrequency, "huevalueAND_frequency.csv");
        // saveHueFrequencyToCSV(filteredMap, "filtered_hue_frequency.csv");
        // System.out.println(imagePath);
        File file = new File(imagePath);
        String fileName = file.getName();
        System.out.println("File Name: " + fileName);
        matrixbuild(fileName);

    }


    private void inputImage(String imagePath) {
        
        double h = 0;
        double s = 0;
        double v = 0;
        double cmax, cmin, diff;

        imgTwo = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        readImageRGB(width, height, imagePath, imgTwo);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                

                int rgbb = imgTwo.getRGB(x, y);
                double r = (rgbb >> 16) & 0xFF;
                double g = (rgbb >> 8) & 0xFF;
                double b = rgbb & 0xFF;
                if (r != g && g != b) {

                r = r / 255.0;
                g = g / 255.0;
                b = b / 255.0;
                cmax = Math.max(r, Math.max(g, b));
                cmin = Math.min(r, Math.min(g, b));
                diff = cmax - cmin;

                if (diff > 0) {
                    if (cmax == cmin)
                        h = 0;
                    else if (cmax == r)
                        h = (60 * ((g - b) / diff) + 360) % 360;
                    else if (cmax == g)
                        h = (60 * ((b - r) / diff) + 120) % 360;
                    else if (cmax == b)
                        h = (60 * ((r - g) / diff) + 240) % 360;

                    if (cmax == 0)
                        s = 0;
                    else
                        s = (diff / cmax);
                } else {
                    s = 0;
                }
                v = cmax;

                if (h < 0) {
                    h = 360 + h;
                }

                Pair p = new Pair(x, y);
                
                // if (!((h == 0 && s == 0 && v == 100) || (h == 0 && s == 100 && v == 0))) {
                double yourThreshold = 0.25;
                // if (!((h == 0 && s == 0 && v == 1) || (h == 0 && s == 1 && v == 0)) && s > yourThreshold) {
                
                inputfh.put(p, h);
                inputfs.put(p, s);
                inputfv.put(p, v);
                // }
                }
            }
        }
        
        System.out.print("image size "+ inputfh.size());

        /////////INPUT IMAGE HISTOGRAM NOT NEEDED RIGHT NOW
        // double[] inputhueValues = new double[inputfh.size()]; // Create an array to hold hue values\
        // int index = 0;
        // for (double hues : inputfh.values()) {
        //     inputhueValues[index++] = hues; // Fill the array with hue values from the fh map
        // }

        // int numBins = 36; // Number of bins for the histogram
        // double[] inputhistogram = generateHueHistogram(inputhueValues, numBins);
        // String outputPath = "histograminputnew.csv";
        // generateHueHistogramToCSV(inputhueValues, numBins, outputPath);

        // System.out.println("Histogram values are " + histogram);
        // for (int i = 0; i < inputhistogram.length; i++) {
        //     System.out.printf("Bins of input img %d: %.4f%n", i, inputhistogram[i]);
        // }


        // for (double hue : inputfh.values()) {
        //     inputHueFrequency.put(hue, inputHueFrequency.getOrDefault(hue, 0.0) + 1.0);
        // }
    
        // // Normalize the hue frequencies
        // double inputfhSize = inputfh.size();
        // for (Map.Entry<Double, Double> entry : inputHueFrequency.entrySet()) {
        //     double normalizedFrequency = (double) entry.getValue() / inputfhSize;
        //     inputHueFrequency.put(entry.getKey(), normalizedFrequency);
        // }
        
        // saveHistogramToCSV(inputfh, "input_hue_histogram.csv");
        // saveHistogramToCSV(inputfs, "input_saturation_histogram.csv");
        // saveHistogramToCSV(inputfv, "input_value_histogram.csv");

        // saveHueFrequencyToCSV(inputHueFrequency, "InputHueValueAND_Frequency.csv");
  
        // return inputfh;
    }


    private void matrixbuild(String objectName) {
        int[][] inputMatrix = new int[width][height];

        for (int x = 0; x < width; x++) { //initially the matrix
            for (int y = 0; y < height; y++) {
                inputMatrix[x][y] = 0;
            }
        }

        double[][] topBinsWithIndices = extractTopBins(histogram, 4);
        double[][] topSaturationBins = extractTopBins(saturationHistogram, 5);
        double[][] topValueBins = extractTopBins(valueHistogram, 8);

        
        
        for (Map.Entry<Pair, Double> entry : inputfh.entrySet()) {
            Pair inputp = entry.getKey();
            double hueval = entry.getValue();
       

            for (int i = 0; i < topBinsWithIndices.length; i++) {
                for (int j = 0; j < topSaturationBins.length; j++) {
                    for (int k = 0; k < topValueBins.length; k++){

                    int topBinIndex = (int) topBinsWithIndices[i][1];
                    int topSaturBinIndex = (int) topSaturationBins[j][1];
                    int topValueBinIndex = (int) topValueBins[k][1];


        
                    double binStart = topBinIndex * (360.0 / 36.0); // Calculate the start of the bin
                    double binEnd = (topBinIndex + 1) * (360.0 / 36.0); // Calculate the end of the bin
                    double saturbinStart = topSaturBinIndex * (1.0 / 10.0); // Calculate the start of the bin
                    double saturbinEnd = (topSaturBinIndex + 1) * (1.0 / 10.0); // Calculate the end of the bin
                    double valbinStart = topValueBinIndex * (1.0 / 10.0); // Calculate the start of the bin
                    double valbinEnd = (topValueBinIndex + 1) * (1.0 / 10.0); // Calculate the end of the bin
                    // Check if the hue value is within the bin
                    if (hueval >= binStart && hueval < binEnd) {
                        double saturationValue = inputfs.get(inputp); // Get the saturation value for the current pixel
                        double valValue = inputfv.get(inputp);
                        
                            if ((saturationValue >= saturbinStart && saturationValue < saturbinEnd) && (valValue >= valbinStart && valValue < valbinEnd)) {
                                int x = (int) inputp.x;
                                int y = (int) inputp.y;
                                inputMatrix[x][y] = 1;
                                break; // No need to check other bins
                        }
                    }
                }
                }
            }

            
        }

        //////////////////////////////////////////////////////////////////////////////////clustering LOGIC here please
        int[][] labelMatrix = new int[width][height];
        int currentLabel = 1;

        List<Cluster> clusters = new ArrayList<>();
        clusters.clear();

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (inputMatrix[x][y] == 1 && labelMatrix[x][y] == 0) {
                    Cluster cluster = new Cluster();
                    cluster.id = currentLabel++;

                    bfsLabel(inputMatrix, labelMatrix, x, y, cluster);

                    clusters.add(cluster);

                }
            }
        }
        // System.out.println("Cluster thingy" + clusters.size());
        
        // Create a map to store the count of pixels per cluster ID
        Map<Integer, Integer> clusterPixelCountMap = new HashMap<>();

        // Count pixels for each cluster
        for (Cluster cluster : clusters) {
            int clusterId = cluster.id;
            int pixelCount = cluster.pixels.size();
            
            // Update the count in the map
            clusterPixelCountMap.put(clusterId, pixelCount);
        }

        // Print the count of pixels per cluster ID
        for (Map.Entry<Integer, Integer> entry : clusterPixelCountMap.entrySet()) {
            int clusterId = entry.getKey();
            int pixelCount = entry.getValue();
            // System.out.println("Cluster " + clusterId + ": Number of Pixels = " + pixelCount);
        }


        for (Cluster cluster : clusters) {
            if (cluster.pixels.size() >= MIN_CLUSTER_SIZE) {
                // System.out.println("Number of Pixels: " + cluster.pixels.size());
                for (Pair pixel : cluster.pixels) {
                    int x = (int) pixel.x;
                    int y = (int) pixel.y;
                    // imgTwo.setRGB(x, y, 0); /////////////////////////////////////////////////////////UNCOMMENT HERE TO SEE THE CLUSTERS GETTING FORMED
                }
            }
        }

        Cluster largestCluster = null;
        int largestClusterSize = 0;
        // System.out.println("Number of Pixels in Largest Cluster: " + largestClusterSize);
        
        // for (Cluster cluster : clusters) {
        //     if (cluster.pixels.size() > largestClusterSize) {
        //         largestCluster = cluster;
        //         largestClusterSize = cluster.pixels.size();
        //     }
        // }
        for (Cluster cluster : clusters) {
            if (cluster.pixels.size() > largestClusterSize) {
                if (largestCluster != null) {
                    largestCluster.clear(); // Clear the previous largest cluster
                }
                largestCluster = cluster;
                largestClusterSize = cluster.pixels.size();
            }
        }
        // for (Cluster cluster : clusters) {
        //     System.out.println("Cluster ID clust: " + cluster.id);
        //     System.out.println("Number of Pixels in Cluster clust: " + cluster.pixels.size());
        //     // You can print more information about the cluster if needed
        // }


        List<Cluster> selectedClusters = new ArrayList<>();

        // Define the size threshold for inclusion (80% of the largest cluster)
        double sizeThreshold = 0.8 * largestClusterSize;
        for (Cluster cluster : clusters) {
            int clusterSize = cluster.pixels.size();
            if (clusterSize >= sizeThreshold) {
                selectedClusters.add(cluster);
            }
        }
        //print selected clusters
        // for (Cluster cluster : selectedClusters) {
        //     System.out.println("Cluster ID select: " + cluster.id);
        //     System.out.println("Number of Pixels in Cluster select: " + cluster.pixels.size());
        //     // You can print more information about the cluster if needed
        // }
        
        clusters.clear();
        labelMatrix = new int[width][height];

        // System.out.println("Number of Pixels in Largest Cluster2: " + largestClusterSize);

        // if (largestCluster != null) {
        for (Cluster clust : selectedClusters) {
            // System.out.println("Largest Cluster ID: " + clust.id);
            // System.out.println("Number of Pixels in Largest Cluster3: " + largestClusterSize);

            int minX = Integer.MAX_VALUE;
            int minY = Integer.MAX_VALUE;
            int maxX = Integer.MIN_VALUE;
            int maxY = Integer.MIN_VALUE;

            for (Pair pixel : clust.pixels) {
                int x = (int) pixel.x;
                int y = (int) pixel.y;

                minX = Math.min(minX, x);
                minY = Math.min(minY, y);
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
            
        }
        // largestCluster.clear();
        // if (largestCluster != null) {
        //     largestCluster.clear();
        // }
            
            // Draw bounding box around the largest cluster in red (0xFFFF0000)
            for (int x = minX; x <= maxX; x++) {
                imgTwo.setRGB(x, minY, 0xFFFF0000); // Red, top border
                imgTwo.setRGB(x, maxY, 0xFFFF0000); // Red, bottom border
            }
            for (int y = minY; y <= maxY; y++) {
                imgTwo.setRGB(minX, y, 0xFFFF0000); // Red, left border
                imgTwo.setRGB(maxX, y, 0xFFFF0000); // Red, right border
            }

            Graphics g = this.imgTwo.getGraphics();
            Font font = new Font("Arial", Font.BOLD, 15);
            g.setColor(Color.DARK_GRAY);
            g.setFont(font);
            g.drawString(objectName, minX, maxY);
            this.frame.repaint();
        }

    }


    private void bfsLabel(int[][] inputMatrix, int[][] labelMatrix, int x, int y, Cluster cluster) {
        int[] dx = {-1, 1, 0, 0}; // Possible x-direction movements
        int[] dy = {0, 0, -1, 1}; // Possible y-direction movements

        

        LinkedList<Pair> queue = new LinkedList<>();
        queue.clear();
        queue.add(new Pair(x, y)); // Start from the initial pixel

        // FileWriter csvWriter = null;

        // try {
            // csvWriter = new FileWriter("input_matrix.csv");

            // for (int i = 0; i < width; i++) {
            //     for (int j = 0; j < height; j++) {
            //         csvWriter.append(inputMatrix[i][j] + ",");
            //     }
            //     csvWriter.append("\n");
            // }

            while (!queue.isEmpty()) {
                Pair currentPixel = queue.poll();
                x = (int) currentPixel.x;
                y = (int) currentPixel.y;

                if (x < 0 || x >= width || y < 0 || y >= height || inputMatrix[x][y] != 1 || labelMatrix[x][y] != 0) {
                    continue;
                }

                // Mark the pixel as visited and assign the cluster's label
                labelMatrix[x][y] = cluster.id;

                // Add the pixel to the current cluster
                cluster.pixels.add(new Pair(x, y));

                // Enqueue neighboring pixels
                for (int i = 0; i < 4; i++) {
                    int newX = x + dx[i];
                    int newY = y + dy[i];
                    if (newX >= 0 && newX < width && newY >= 0 && newY < height) {
                        queue.add(new Pair(newX, newY));
                    }
                }
            }
        // } catch (IOException e) {
        //     e.printStackTrace();
        // } 
        // finally {
            // try {
            //     if (csvWriter != null) {
            //         csvWriter.flush();
            //         csvWriter.close();
            //     }
            // } catch (IOException e) {
            //     e.printStackTrace();
            // }
        // }
        // System.out.println("Cluster ID: " + cluster.id);
        // System.out.println("Number of Pixels in Cluster: " + cluster.pixels.size());
    }



    public static void main(String[] args) {
        ImageDisplay ren = new ImageDisplay();
        ren.showIms(args);
    }
}
