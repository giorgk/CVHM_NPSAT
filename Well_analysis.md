Well Data Analysis
===============

The following is a statistical analysis on the available Central Valley well data. The data will be analyzed with respect to the parameters that are usefull for the Random well generation algorithm. 

- Well density
- Spatial variability of pumping capacity
- Relation between well depth and pumping
- Relation between well screen length depth and pumping

####Overview of data


The total number of well records in Central Valley that are characterized as public or agricultural is **50,627**. 
<img src="CVallWells1.png" alt="Wells in Central Valley" width="150"/>

From these records **9,704** have also well yield. 
<img src="wells_with_pumping1.png" alt="Wells in Central Valley" width="700"/>

In the well generation algorithm we are interested in developing a relation between Depth and Pumping rate. Out of the 9,704 records **9,348** have both pumping rates and depth. The plots below show the spatial distribution (left), the historgam of depth across Central Valley (right top) and the relation between pumping rate and depth. The red solid line corresponds to a fitted linear model, while the dashed lines corresponds to 95 prediction intervals.
<img src="Wells_with_D_Q1.png" alt="Wells in Central Valley" width="600"/>

Last for the well generation we need a relation between pumping, depth and screen length. In total there are **6,676** records that have all the information available. 
<img src="Wells_with_D_Q_SL1.png" alt="Wells in Central Valley" width="600"/>

####Well density
To compute the density map for Central Valley we used a Kernel smooth function. For this data set we used a gaussian kernel function with 3 km
bandwith. The bandwidth was found empirically by trial & error. 


 