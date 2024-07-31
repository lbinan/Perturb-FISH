# Perturb-FISH
Contains all scripts used in the Perturb-FISH publication.<br/>
<br/>
Scripts are grouped by dataset: THP1 (1st dataset, shows analysis of spatial effects: density and neighborhood), astrocyte (2nd dataset, shows matching perturbation, transcriptome and time lapse live cell phenotyping), and finally tumor xenograft (shows in vivo data, with spatial effects analysis).<br/>
<br/>
For each dataset, scripts are grouped into preprocessing and analysis. Preprocessing takes the user from images to count tables. It filers the images, decodes the identity of the perturbations,  segments cells, and quantifies signal in each cell. It also stiches images into a mosaic. Analysis countains the scripts that take the count tabels as input, align data between modalities for astrocytes and tumor data, and analyse perturbation effects, check the consistency of the effects between experimental replicates, and between randomly designed subsets of the data, compute spatial and phenotypic effects, clusters cell, and prepare plots.<br/>
<br/>
When working from images, in preprocessing, users should:<br/>
-run the filter script on the images of the perturbation. (alligns images between rounds of imaging, evens the shape of the illumination profile, and applies gaussian filters to help spot detection)<br/>
-set threshold values to identify gRNA spots on images, either manually or using the getvaluesforthresholdfilteredimagesclearedthp1.m script from the dataexploration folder, or a combination of both.<br/>
-run the decodeHW4 script (finds the binary barcode each detected spot corresponds to)<br/>
-run the get geneID (matches the binary barcode to a perturbation from the codebook)<br/>
-run the spotsQC (this allows checking the quality of decoding, and can help iteratively improving the thresholds set at step 2. it also collects average intensity values and standard deviation in the intensity of each spot through the rounds of imaging. this info is used to filter spots based on decoding confidence.<br/>
-run the makeseedsimageAGAIN script, which stiches images into a mosaic, and prepares a mosaic of cell centers (from nuclei) and contours (from transcript density) to run watershed and segment all cells (the first time, uncomment the first part that determines how big the mosaic will be)<br/>
-run the finishglobalcellmosaic, which finishes the segmentation, and converts the spatial coordinates of devcoded gRNAs and RNA transcripts into global coordinates, then build 2 count tables, one for perturbations, one for RNA transcripts<br/>
-for the astrocyte and the tumor data, data modalities are acquired on different micropscopes and therefore reuquire registration. all images are registered to the MERFISH (cells are segmented from stiched MERFISH images, the other modalites are alligned to this mask, because the data most sensitive to alignment is the MERFISH data). the scripts getTranslationMtocalcium align the mask from merfish onto the calcium images, to quaantify calcium traces through time. the scripts getTranslationMtoZ align the MERFISH mask onto the stiched images of the perutrbations<br/>
<br/>
the analysis folder contains the scripts that split the data randomly to check correlations and verify that the detected signal remains consitant. it also contains the scripts used to infer perturbation effects (in python, called using the scripts in the txt files from a terminal)


