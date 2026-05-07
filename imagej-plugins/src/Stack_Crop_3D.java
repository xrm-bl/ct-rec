import ij.*;
import ij.plugin.PlugIn;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.ImageProcessor;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Arrays;
import java.util.Comparator;

public class Stack_Crop_3D implements PlugIn {
    
    private JTextField inputDirField, outputDirField;
    private JTextField x1Field, y1Field, z1Field;
    private JTextField x2Field, y2Field, z2Field;
    private PrintWriter logWriter;
    private long startTime;
    private String logFilePath;
    
    public void run(String arg) {
        showDialog();
    }
    
    private void showDialog() {
        // Create main dialog
        JDialog dialog = new JDialog();
        dialog.setTitle("3D Stack Crop Plugin");
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
        
        // Start point coordinates
        gbc.gridx = 0; gbc.gridy = 2;
        dialog.add(new JLabel("Start Point (X1, Y1, Z1):"), gbc);
        gbc.gridx = 1;
        x1Field = new JTextField(10);
        dialog.add(x1Field, gbc);
        gbc.gridx = 2;
        y1Field = new JTextField(10);
        dialog.add(y1Field, gbc);
        gbc.gridx = 3;
        z1Field = new JTextField(10);
        dialog.add(z1Field, gbc);
        
        // End point coordinates
        gbc.gridx = 0; gbc.gridy = 3;
        dialog.add(new JLabel("End Point (X2, Y2, Z2):"), gbc);
        gbc.gridx = 1;
        x2Field = new JTextField(10);
        dialog.add(x2Field, gbc);
        gbc.gridx = 2;
        y2Field = new JTextField(10);
        dialog.add(y2Field, gbc);
        gbc.gridx = 3;
        z2Field = new JTextField(10);
        dialog.add(z2Field, gbc);
        
        // Buttons
        JPanel buttonPanel = new JPanel();
        JButton runBtn = new JButton("Run");
        JButton cancelBtn = new JButton("Cancel");
        
        runBtn.addActionListener(e -> {
            if (validateInputs()) {
                dialog.dispose();
                processCrop();
            }
        });
        
        cancelBtn.addActionListener(e -> dialog.dispose());
        
        buttonPanel.add(runBtn);
        buttonPanel.add(cancelBtn);
        
        gbc.gridx = 0; gbc.gridy = 4;
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
            
            // Parse coordinates
            int x1 = Integer.parseInt(x1Field.getText());
            int y1 = Integer.parseInt(y1Field.getText());
            int z1 = Integer.parseInt(z1Field.getText());
            int x2 = Integer.parseInt(x2Field.getText());
            int y2 = Integer.parseInt(y2Field.getText());
            int z2 = Integer.parseInt(z2Field.getText());
            
            // Validate coordinate order
            if (x1 > x2 || y1 > y2 || z1 > z2) {
                IJ.error("Start coordinates must be less than or equal to end coordinates");
                return false;
            }
            
            return true;
        } catch (NumberFormatException e) {
            IJ.error("Invalid coordinate values. Please enter integers only.");
            return false;
        }
    }
    
    private void processCrop() {
        try {
            // Initialize logging
            initializeLogging();
            
            File inputDir = new File(inputDirField.getText());
            File outputDir = new File(outputDirField.getText());
            
            int x1 = Integer.parseInt(x1Field.getText());
            int y1 = Integer.parseInt(y1Field.getText());
            int z1 = Integer.parseInt(z1Field.getText());
            int x2 = Integer.parseInt(x2Field.getText());
            int y2 = Integer.parseInt(y2Field.getText());
            int z2 = Integer.parseInt(z2Field.getText());
            
            logCommand("Starting 3D Stack Crop process");
            logCommand("Input directory: " + inputDir.getAbsolutePath());
            logCommand("Output directory: " + outputDir.getAbsolutePath());
            logCommand("Crop region: (" + x1 + "," + y1 + "," + z1 + ") to (" + x2 + "," + y2 + "," + z2 + ")");
            
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
            
            // Check Z range
            if (z1 < 1 || z2 > tiffFiles.length) {
                IJ.error("Z coordinates out of range. Valid range: 1-" + tiffFiles.length);
                logCommand("ERROR: Z coordinates out of range");
                closeLogging();
                return;
            }
            
            // Read first image to get dimensions
            ImagePlus firstImage = new Opener().openImage(tiffFiles[0].getAbsolutePath());
            if (firstImage == null) {
                IJ.error("Failed to open first image");
                logCommand("ERROR: Failed to open first image");
                closeLogging();
                return;
            }
            
            int width = firstImage.getWidth();
            int height = firstImage.getHeight();
            logCommand("Image dimensions: " + width + " x " + height);
            
            // Validate XY coordinates
            if (x1 < 1 || x2 > width || y1 < 1 || y2 > height) {
                IJ.error("XY coordinates out of range. Valid range: X(1-" + width + "), Y(1-" + height + ")");
                logCommand("ERROR: XY coordinates out of range");
                closeLogging();
                return;
            }
            
            // Process crop
            int cropWidth = x2 - x1 + 1;
            int cropHeight = y2 - y1 + 1;
            int cropDepth = z2 - z1 + 1;
            
            logCommand("Crop dimensions: " + cropWidth + " x " + cropHeight + " x " + cropDepth);
            logCommand("Starting image processing...");
            
            int outputIndex = 0;
            
            // Process each Z slice in the specified range
            for (int z = z1 - 1; z < z2; z++) {
                long sliceStart = System.currentTimeMillis();
                
                // Show progress
                IJ.showStatus("Processing slice " + (z - z1 + 2) + " of " + cropDepth);
                IJ.showProgress(z - z1 + 1, cropDepth);
                
                File tiffFile = tiffFiles[z];
                logCommand("Processing: " + tiffFile.getName() + " (slice " + (z + 1) + ")");
                
                ImagePlus imp = new Opener().openImage(tiffFile.getAbsolutePath());
                if (imp == null) {
                    logCommand("WARNING: Failed to open " + tiffFile.getName());
                    continue;
                }
                
                // Crop the image
                ImageProcessor ip = imp.getProcessor();
                ip.setRoi(x1 - 1, y1 - 1, cropWidth, cropHeight);
                ImageProcessor croppedIp = ip.crop();
                
                ImagePlus croppedImp = new ImagePlus("Cropped", croppedIp);
                
                // Save with zero-padded filename
                String outputFileName = String.format("%05d.tif", outputIndex);
                File outputFile = new File(outputDir, outputFileName);
                
                FileSaver fs = new FileSaver(croppedImp);
                if (fs.saveAsTiff(outputFile.getAbsolutePath())) {
                    logCommand("Saved: " + outputFileName + " (" + 
                              (System.currentTimeMillis() - sliceStart) + " ms)");
                    outputIndex++;
                } else {
                    logCommand("ERROR: Failed to save " + outputFileName);
                }
                
                imp.close();
                croppedImp.close();
            }
            
            long totalTime = System.currentTimeMillis() - startTime;
            logCommand("Process completed successfully");
            logCommand("Total images processed: " + outputIndex);
            logCommand("Total processing time: " + totalTime + " ms");
            
            // Clear progress
            IJ.showProgress(1.0);
            IJ.showStatus("Crop completed: " + outputIndex + " images processed");
            
            IJ.showMessage("3D Stack Crop Complete", 
                          "Processed " + outputIndex + " images\n" +
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
            logCommand("Stack Crop 3D Plugin - Log started");
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
