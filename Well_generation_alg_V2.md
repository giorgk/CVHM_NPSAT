## Well Generation Algorithm V2
#### Data description
In this version the pumping will be distributed according to the CVHM pumping that is showing in the figure. In CVHM there are two types of wells. Multi node and farm wells. The amount of water pumped from each cell at every stress is written in the cell by cell flow output file. In our analysis we averaged the pumping for the period between October 1977 and February 2003. To make a distinction between pumped water used for urban ang agriculture we overlaid the pumping distribution map to the 2000 land use map. 84% of the pumped water was estimated that is used for agricultural  and 16% for urban uses. The ratio between Ag and urban water use is shown for each farm as piechart below.

<img src="CVHM_pumping.png" alt="CVHM pumping" width="700"/>


The well data base contains **50,627** records. When we exclude the wells with unknown or erroneous (190, 2028 etc) construction year the number of wells is **47,510**. The majority of the wells have been constructed after 1950. However in the algorithm we use only the relatively new wells (e.g less than 20-30 years old). The second plot shows the number of wells with respect to the well age. For example we see that there are approximately 5,000 wells that are 10 years old. Here we set the age threshold for public and agricultural wells equal to 30 years. We can see that there are about 18,000 wells that are maximum 30 years old. In addition the spatial distribution of the 18,000 wells covers the entire Central Valley.

<img src="HistyearConstruct.png" alt="Histogram of well age" width="700"/>
<img src="Nwells_vs_WellAge.png" alt="Number of wells vs well age" width="700"/>
<img src="SpatialCovNewwells.png" alt="Spatial Coverage of new wells" width="400"/>

#### Well generation algorithm

The algorithm loops through the 21 CVHM Farms.
For each farm:
1. Find the wells that correspond to the farm.
2. Identify which wells are characterized as public and agricultural. For example for farm 21 there are 466 wells with 413 aggricultiral and 53 public.
<img src="farm_21_wells.png" alt="Well locations for farm 21" width="500"/>
3. For each of the two sets of agricultural and public wells populate the missing records. The following statistics are executed once for each set.

- Pumping rates.
For the pumping rates compute the empirical cumulative distribution function based on the data that have known rates. Using the ECDF generate random pumping for the records that dont have pumping.
Generate a random probability between [0,1]. Assign the corresponding pumping to that probability. 
 <img src="farm21_Q_ecdf.png" alt="ECDF of pumping" width="500"/>
 - Depth
 For the depth we identify the records that have both pumping and depth assigned. Note that we do not use the random pumping in that calculation. We fit a linear model and calculate the prediction interval that corresponds to standard deviation. Then for the missing depths we calculate a mean and standard deviation depth as a function of pumping. For this calculation we use also the gerenated pumpings. A normally distributed random depth is generated based on the mean and standard that were calculated as function of assigned pumping. 
 <img src="farm21_QD_fit.png" alt="Q vs depth" width="500"/>
 - Screen length
 A similar approach was used do the screen length with the main difference beign that the mean and standard deviation of screen length is a function of a pumping and depth.
 
 The above algorithm populates the missing records using independent statistics for each farm and each set of public and agricultural wells.
 However there are some cases were that was not possible. For the public wells, for example, there were farms with very few records that was not possible to compute the required statistics. In those cases we expanded the data set and use the all subbasins records (e.g. Sacramento Valley, San Joaquin Valley Tularey lake basin). In addition, for farm 19 there are not available records with all three of the well properties available and we included for the statistic calclulation the farms 20 and 21. 


#### Correct Pumping
The above algorithm populated all the needed missing records. However the total amount of pumping does not match the amount of pumping that was estimated from the CVHM model. To correct the pumping we scaled the pumping rates for each farm and each ag or public well set to match the total pumping per farm.

#### Correct location
The well sets contains six methods of coordinate determination. For two of them we use the provided coordinates, while for the remaining the coordinates were translated up to 800 m on each direction:

| Method of determination | In Algorithm |
| -----------------------------     | ---------------- |
| Derived from Address      | Use 
| Derived from TRS       	   | Translate
| GPS				   | USE
| NA					   | Translate
| Other				| Translate
|Unknown			| Translate

**Clarify what the other methods mean**
<img src="translatedWellPos.png" alt="Q vs depth" width="800"/>
  