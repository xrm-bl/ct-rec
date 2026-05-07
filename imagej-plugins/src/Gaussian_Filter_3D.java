import ij.*;
import ij.plugin.PlugIn;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.*;

public class Gaussian_Filter_3D implements PlugIn {
    
    private JTextField inputDirField, outputDirField;
    private JTextField sigmaField;
    private PrintWriter logWriter;
    private long startTime;
    private String logFilePath;
    
    public void run(String arg) {
        showDialog();
    }
    
    private void showDialog() {
        // Create main dialog
        JDialog dialog = new JDialog();
        dialog.setTitle("3D Gaussian Filter Plugin");
        dialog.setModal(true);
        dialog.setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.insets = new Insets(5, 5, 5, 5);
        
        // Input directory
        gbc.gridx = 0; gbc.gridy = 0;
        dialog.add(new JLabel("Input Directory:"), gbc);
        gbc.gridx = 1; gbc.gridwidth = 2;
        inputDirField = new JTextField(30);
        dialog.add(inputDirField, gbc);
        gbc.gridx = 3; gbc.gridwidth = 1;
        JButton inputBrowseBtn = new JButton("Browse...");
        inputBrowseBtn.addActionListener(e -> browseDirectory(inputDirField));
        dialog.add(inputBrowseBtn, gbc);
        
        // Output directory
        gbc.gridx = 0; gbc.gridy = 1;
        dialog.add(new JLabel("Output Directory:"), gbc);
        gbc.gridx = 1; gbc.gridwidth = 2;
        outputDirField = new JTextField(30);
        dialog.add(outputDirField, gbc);
        gbc.gridx = 3; gbc.gridwidth = 1;
        JButton outputBrowseBtn = new JButton("Browse...");
        outputBrowseBtn.addActionListener(e -> browseDirectory(outputDirField));
        dialog.add(outputBrowseBtn, gbc);
        
        // Sigma (kernel size)
        gbc.gridx = 0; gbc.gridy = 2;
        dialog.add(new JLabel("Sigma (s):"), gbc);
        gbc.gridx = 1;
        sigmaField = new JTextField("2.0", 10);
        dialog.add(sigmaField, gbc);
        gbc.gridx = 2;
        dialog.add(new JLabel("(Standard deviation for Gaussian kernel)"), gbc);
        
        // Buttons
        JPanel buttonPanel = new JPanel();
        JButton runBtn = new JButton("Run");
        JButton cancelBtn = new JButton("Cancel");
        
        runBtn.addActionListener(e -> {
            if (validateInputs()) {
                dialog.dispose();
                processFilter();
            }
        });
        
        cancelBtn.addActionListener(e -> dialog.dispose());
        
        buttonPanel.add(runBtn);
        buttonPanel.add(cancelBtn);
        
        gbc.gridx = 0; gbc.gridy = 3;
        gbc.gridwidth = 4;
        dialog.add(buttonPanel, gbc);
        
        // Set dialog properties
        dialog.pack();
        dialog.setLocationRelativeTo(null);
        dialog.setVisible(true);
    }
    
    private void browseDirectory(JTextField field) {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
            field.setText(chooser.getSelectedFile().getAbsolutePath());
        }
    }
    
    private boolean validateInputs() {
        try {
            // Check directories
            File inputDir = new File(inputDirField.getText());
            File outputDir = new File(outputDirField.getText());
            
            if (!inputDir.exists() || !inputDir.isDirectory()) {
                IJ.error("Invalid input directory");
                return false;
            }
            
            if (!outputDir.exists()) {
                outputDir.mkdirs();
            }
            
            // Parse sigma
            double sigma = Double.parseDouble(sigmaField.getText());
            
            // Validate sigma
            if (sigma <= 0) {
                IJ.error("Sigma must be positive");
                return false;
            }
            
            return true;
        } catch (NumberFormatException e) {
            IJ.error("Invalid sigma value. Please enter a valid number.");
            return false;
        }
    }
    
    private void processFilter() {
        try {
            // Initialize logging
            initializeLogging();
            
            File inputDir = new File(inputDirField.getText());
            File outputDir = new File(outputDirField.getText());
            
            double sigma = Double.parseDouble(sigmaField.getText());
            
            logCommand("Starting 3D Gaussian Filter process");
            logCommand("Input directory: " + inputDir.getAbsolutePath());
            logCommand("Output directory: " + outputDir.getAbsolutePath());
            logCommand("Sigma (s): " + sigma);
            
            // Get TIFF files sorted by name
            File[] tiffFiles = inputDir.listFiles((dir, name) -> 
                name.toLowerCase().endsWith(".tif") || name.toLowerCase().endsWith(".tiff"));
            
            if (tiffFiles == null || tiffFiles.length == 0) {
                IJ.error("No TIFF files found in input directory");
                logCommand("ERROR: No TIFF files found in input directory");
                closeLogging();
                return;
            }
            
            // Sort files by name to ensure proper Z order
            Arrays.sort(tiffFiles, Comparator.comparing(File::getName));
            logCommand("Found " + tiffFiles.length + " TIFF files");
            
            // Load all images into memory for 3D processing
            logCommand("Loading images into memory...");
            ImagePlus[] images = new ImagePlus[tiffFiles.length];
            int width = 0, height = 0;
            
            for (int i = 0; i < tiffFiles.length; i++) {
                IJ.showStatus("Loading image " + (i + 1) + " of " + tiffFiles.length);
                IJ.showProgress(i, tiffFiles.length);
                
                images[i] = new Opener().openImage(tiffFiles[i].getAbsolutePath());
                if (images[i] == null) {
                    IJ.error("Failed to open: " + tiffFiles[i].getName());
                    logCommand("ERROR: Failed to open " + tiffFiles[i].getName());
                    closeLogging();
                    return;
                }
                
                if (i == 0) {
                    width = images[i].getWidth();
                    height = images[i].getHeight();
                    logCommand("Image dimensions: " + width + " x " + height + " x " + tiffFiles.length);
                    logCommand("Bit depth: " + images[i].getBitDepth() + " bits");
                }
            }
            
            // Create 3D array for processing
            logCommand("Creating 3D data structure...");
            float[][][] volume = new float[tiffFiles.length][height][width];
            
            // Convert to float array
            for (int z = 0; z < tiffFiles.length; z++) {
                FloatProcessor fp = images[z].getProcessor().toFloat(0, null);
                float[] pixels = (float[]) fp.getPixels();
                for (int y = 0; y < height; y++) {
                    System.arraycopy(pixels, y * width, volume[z][y], 0, width);
                }
            }
            
            // Apply 3D Gaussian filter
            logCommand("Applying 3D Gaussian filter...");
            float[][][] filteredVolume = apply3DGaussianFilter(
                volume, width, height, tiffFiles.length, sigma
            );
            
            // Save filtered images
            logCommand("Saving filtered images...");
            for (int z = 0; z < tiffFiles.length; z++) {
                IJ.showStatus("Saving image " + (z + 1) + " of " + tiffFiles.length);
                IJ.showProgress(z, tiffFiles.length);
                
                // Convert back to ImagePlus
                float[] pixels = new float[width * height];
                for (int y = 0; y < height; y++) {
                    System.arraycopy(filteredVolume[z][y], 0, pixels, y * width, width);
                }
                
                FloatProcessor fp = new FloatProcessor(width, height, pixels);
                ImagePlus filteredImage = new ImagePlus("Filtered", fp);
                
                // Preserve original bit depth
                if (images[z].getBitDepth() == 8) {
                    filteredImage = new ImagePlus("Filtered", fp.convertToByte(true));
                } else if (images[z].getBitDepth() == 16) {
                    filteredImage = new ImagePlus("Filtered", fp.convertToShort(true));
                } else if (images[z].getBitDepth() == 32) {
                    // Keep as 32-bit float
                    filteredImage = new ImagePlus("Filtered", fp);
                }
                
                // Save with original filename
                File outputFile = new File(outputDir, tiffFiles[z].getName());
                FileSaver fs = new FileSaver(filteredImage);
                
                if (fs.saveAsTiff(outputFile.getAbsolutePath())) {
                    logCommand("Saved: " + tiffFiles[z].getName());
                } else {
                    logCommand("ERROR: Failed to save " + tiffFiles[z].getName());
                }
                
                // Clean up
                images[z].close();
                filteredImage.close();
            }
            
            long totalTime = System.currentTimeMillis() - startTime;
            logCommand("Process completed successfully");
            logCommand("Total processing time: " + totalTime + " ms");
            
            // Clear progress
            IJ.showProgress(1.0);
            IJ.showStatus("Gaussian filter completed");
            
            IJ.showMessage("3D Gaussian Filter Complete", 
                          "Processed " + tiffFiles.length + " images\n" +
                          "Output directory: " + outputDir.getAbsolutePath() + "\n" +
                          "Log file: " + logFilePath);
            
        } catch (Exception e) {
            IJ.error("Error during processing: " + e.getMessage());
            logCommand("ERROR: " + e.getMessage());
            e.printStackTrace();
        } finally {
            closeLogging();
        }
    }
    
    private float[][][] apply3DGaussianFilter(float[][][] volume, int width, int height, int depth,
                                              double sigma) {
        
        float[][][] output = new float[depth][height][width];
        
        // Calculate kernel radius (3 sigma coverage)
        int radius = (int)Math.ceil(3 * sigma);
        int kernelSize = 2 * radius + 1;
        
        logCommand("Kernel radius: " + radius + " (kernel size: " + kernelSize + "x" + kernelSize + "x" + kernelSize + ")");
        
        // Create 3D Gaussian kernel
        double[][][] kernel = create3DGaussianKernel(radius, sigma);
        
        // Apply convolution using parallel processing
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        logCommand("Using " + numThreads + " threads for parallel processing");
        
        // Process each slice in parallel
        CountDownLatch latch = new CountDownLatch(depth);
        
        for (int z = 0; z < depth; z++) {
            final int currentZ = z;
            executor.submit(() -> {
                try {
                    IJ.showStatus("Processing slice " + (currentZ + 1) + " of " + depth);
                    IJ.showProgress(currentZ, depth);
                    
                    for (int y = 0; y < height; y++) {
                        for (int x = 0; x < width; x++) {
                            output[currentZ][y][x] = applyGaussianAtVoxel(
                                volume, x, y, currentZ, width, height, depth,
                                radius, kernel
                            );
                        }
                    }
                    
                    logCommand("Completed slice " + (currentZ + 1));
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await();
        } catch (InterruptedException e) {
            logCommand("ERROR: Processing interrupted");
        }
        
        executor.shutdown();
        
        return output;
    }
    
    private double[][][] create3DGaussianKernel(int radius, double sigma) {
        int size = 2 * radius + 1;
        double[][][] kernel = new double[size][size][size];
        double sum = 0.0;
        double sigma2 = 2.0 * sigma * sigma;
        
        // Calculate kernel values
        for (int z = -radius; z <= radius; z++) {
            for (int y = -radius; y <= radius; y++) {
                for (int x = -radius; x <= radius; x++) {
                    double distance2 = x*x + y*y + z*z;
                    double value = Math.exp(-distance2 / sigma2);
                    kernel[z+radius][y+radius][x+radius] = value;
                    sum += value;
                }
            }
        }
        
        // Normalize kernel
        for (int z = 0; z < size; z++) {
            for (int y = 0; y < size; y++) {
                for (int x = 0; x < size; x++) {
                    kernel[z][y][x] /= sum;
                }
            }
        }
        
        return kernel;
    }
    
    private float applyGaussianAtVoxel(float[][][] volume, int x, int y, int z,
                                       int width, int height, int depth, int radius,
                                       double[][][] kernel) {
        
        double sum = 0.0;
        double weightSum = 0.0;
        
        // Apply convolution
        for (int dz = -radius; dz <= radius; dz++) {
            int nz = z + dz;
            if (nz < 0 || nz >= depth) continue;
            
            for (int dy = -radius; dy <= radius; dy++) {
                int ny = y + dy;
                if (ny < 0 || ny >= height) continue;
                
                for (int dx = -radius; dx <= radius; dx++) {
                    int nx = x + dx;
                    if (nx < 0 || nx >= width) continue;
                    
                    double weight = kernel[dz+radius][dy+radius][dx+radius];
                    sum += weight * volume[nz][ny][nx];
                    weightSum += weight;
                }
            }
        }
        
        // Handle edge cases where not all kernel weights were used
        return (float)(sum / weightSum);
    }
    
    private void initializeLogging() {
        try {
            startTime = System.currentTimeMillis();
            
            // Create log directory
            String userHome = System.getProperty("user.home");
            File logDir = new File(userHome, "com-log");
            if (!logDir.exists()) {
                logDir.mkdirs();
            }
            
            // Create log file with timestamp
            SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMdd_HHmmss");
            String timestamp = sdf.format(new Date());
            File logFile = new File(logDir, timestamp + ".log");
            logFilePath = logFile.getAbsolutePath();
            
            logWriter = new PrintWriter(new FileWriter(logFile));
            logCommand("Gaussian Filter 3D Plugin - Log started");
            logCommand("Timestamp: " + new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
        } catch (IOException e) {
            IJ.error("Failed to create log file: " + e.getMessage());
        }
    }
    
    private void logCommand(String message) {
        if (logWriter != null) {
            long elapsed = System.currentTimeMillis() - startTime;
            String logEntry = String.format("[%8d ms] %s", elapsed, message);
            logWriter.println(logEntry);
            logWriter.flush();
            IJ.log(logEntry);
        }
    }
    
    private void closeLogging() {
        if (logWriter != null) {
            logCommand("Log closed");
            logWriter.close();
        }
    }
}