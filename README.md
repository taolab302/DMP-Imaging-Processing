# DMP Imaging processing

Scripts for behavior & calcium imaging analysis were written using MATLAB R2018a (MathWorks). You can run the script `main.m` by code sections to go through the following processing pipeline. All output files will be stored in `waveFolder` that you set in `main.m`. Some of results generated in one section might need proofreading before the next section is executed.

1. Worm segmentation   
   Output: 
   - Generate binary images of worm body (`worm_region\`). 
   - Merged images of the original images and segmented worm region (`worm_region_check\`). 
   - Starting point of the centerline in each frame (`centerline_start_xy.txt`). 

2. Centerline calculation   
   Output:
   - Centerline of the worm in each frame (`centerline\*.mat`).
   - Images where centerlines are plotted  (`fig\`).
   - Lengths of centerlines & Lengths of the anterior and posterior halves. (`wormLength.mat`).

3. AVL & DVB tracking
   Before tracking, use `neuron_seg\ExtractNeurons.m` to segment candidate points for neuron positions in each frame.  
   Output: 
   - AVL & DVB positions in GCaMP and RFP images (`neuron_pos\red\`, `neuron_pos\green\`).

4. AVL & DVB fluorescence intensities  
   Output:  
   - AVL & DVB fluorescence intensities in GCaMP and RFP images (stored in `neuron_pos\` as `.mat` files).

5. Intestinal fluorescence intensities  
   Output:
   - Intestinal fluorescence intensities in GCaMP and RFP images (`waveIntensity.mat`).

# DMP-Imaging-Processing
