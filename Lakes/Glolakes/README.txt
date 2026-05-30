Overview:
The GloLakes dataset provides the first near-global, multi-decadal record of lake and reservoir water storage dynamics by combining satellite altimetry and optical remote sensing to monitor changes in water volume for more than 27 000 lakes worldwide from 1984 to the present. The dataset integrates multi-source satellite observations, including optical imagery from Landsat and Sentinel-2 and satellite altimetry from TOPEX/Poseidon, Jason-1/2/3, Sentinel-3, Sentinel-6, and ICESat-2, to estimate changes in lake surface area, water level, and total water storage. GloLakes supports ongoing monitoring of lakes and reservoirs by incorporating newly available remote sensing observations as they become available.

By providing consistent, spatially extensive, and openly accessible lake-storage information, GloLakes addresses the scarcity of in-situ observations in many regions and supports applications in hydrology, climate change analysis, water resources management, ecosystem studies, and drought and flood assessment. The dataset is distributed in standard formats and is accompanied by comprehensive metadata and documentation to facilitate direct use in scientific research and operational workflows.

The latest version was released in December 2025. Additional background information and frequently asked questions are available in the NCAR Climate Data Guide: https://climatedataguide.ucar.edu/climate-data/lake-water-storage-glolakes. 


Data Availability:
The GloLakes dataset is publicly available through https://doi.org/10.25914/K8ZF-6G46 and a web-based data explorer
(Global Water Monitor (GWM): http://www.globalwater.online). By using GWM, you can check water storage time series and statistics either for indivdial lake, basin/catchment or country/state(province).  


Usage Instructions:
It has six products in NetCDF format and each product provides lake storage time series, lake ID (linked to the HydroLAKES database) and name, data quality, latitidue, longitude, and ID (linked to the HydroBASINS database) and name of basin, catchment, country, and state(province). The data quality labels (Q1-Q4) represent the methods used to calculate lake storage for the specific period:

Q1: absolute volume estimated using geostatistical model and satellite-derived lake extents; 
Q2: absolute volume estimated based on the volume-heihgt (V-H) relationship; 
Q3: relative volume estimated based on both satellite-derived heights and extents; 
Q4: relative volume estimated based on heights obtained from satellite measurements, combined with the extents of the lake derived from the area-height (A-H) relationship


Updates:
Dec-2025: Based on feedback from end users, changes in data accessibility for certain sources, and improvements in other data products, the GloLakes dataset was refined to focus on maintaining two high-quality products: Landsat + ICESat-2 and Landsat + G-REALM. Given resource constraints, development efforts were concentrated on these two products to ensure regular updates, as well as consistency, reliability, and high data quality. 


If you use this dataset, please include the following citation in your publication:

Hou, J., Van Dijk, A. I. J. M., Renzullo, L. J., and Larraondo, P. R.: GloLakes: a database of global lake water storage dynamics from 1984 to present derived using laser and radar altimetry and optical remote sensing, Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2022-266, in review, 2022.


For further inquiries or assistance, please contact:

Dr. Jiawei Hou
Fenner School of Environment and Society
College of Science
Australian National University
Canberra ACT 2601, AUSTRALIA
Email: jiawei.hou@anu.edu.au
