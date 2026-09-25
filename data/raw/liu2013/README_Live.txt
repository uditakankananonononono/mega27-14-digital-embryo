LiveImaging_data_structure.mat contains a data structure named "FlyLines" with Bicoid-GFP intensity gradient and cephalic furrow (CF) position measurements of fly embryos from live two-photon imaging for 29 fly lines.

Each record in the structure contains the following fields:
1) FlyLineNumber: fly line number as listed in column 1 in Table S1 of Liu et al. PNAS 2013.
2) FlyLineName: fly line name as listed in column 2 of Table S1 of Liu et al. PNAS 2013.
3) Embryos: structure that contains the source data for one representative live imaging session. Size of this structure equals the number of Embryos N imaged in this session. 
---Embryos.Gradient.left and Em.Gradient.right contain data for N gradients on left and right sides of the center line. Each has two columns: column1 contains the position of each detected nucleus in unit of egg length; column2 contains the apparent nuclear intensity of each nucleus after background subtraction. Bcd gradients were measured in the mid-coronal plane under dorsal view at about 16 min into nuclear cycle 14.  
---Embryos.CF contains N cephalic furrow positions in unit of egg length measured ~50 min after Bcd gradient measurement.NaN means no data is available (typically because only the gradients was measured for this particular embryo). 
---Embryos.EggLength: Egg lengths for N embryos in micrometer.
---Embryos.Orientation.LR: Rotational asymmetry around the left-right axis. Three classes of embryos identified based on the range of faint membrane segments: LR=0 (<1%EL), LR=1 or -1(<3%EL), and LR=2 or -2 (>3%EL). LR<0 and LR>0 indicates that the faint membrane segments are in the anterior pole and posterior pole, respectively. Only the embryos with LR=0 were selected for further data analysis. (See Liu et al. PNAS 2013.)
---Embryos.Orientation.AP: Rotational asymmetry around the AP axis. Two groups of embryos identified based on shape symmetry. AP=1 if the image plane is closer to midsagittal plane, and  AP=0 if it is closer to the coronal plane. Only embryos with AP=0 were selected for further data analysis. (See Liu et al. PNAS 2013.)

