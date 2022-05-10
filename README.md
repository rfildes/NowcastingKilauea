# NowcastingKilauea
This repository contains the python code and data files for the submitted manuscript Fildes et al. (2022) on nowcasting caldera collapse events during the 2018 Kīlauea volcano eruption.

KilaueaEQs.py is a python file written in Spyder (Spyder-IDE.org ) 4.0.1 with Python 3.7. This file contains the functions used to carry out the data analysis in the paper. Within KilaeuaEQs.py, choose which txt file of data to read in based on the subcase you want to run. 

Kilauea4.txt is the data file of earthquakes for the summit focused data (subcases S) and HawaiiIsland.txt is the data file of earthquakes for the whole island (subcases I). Earthquake data acquired on 8/21/2019 (summit focused data) and 11/22/2019 (whole island data) from the United States Geological Survey Earthquake Catalog: https://earthquake.usgs.gov/earthquakes/search/.

Kilauea4.txt:  Earthquakes of Magnitude >= 2, starting 2018-05-04, 00:00:00.00 UTC, within 5 km radius of 19.406705 deg latitude and -155.283386 deg longitude

HawaiiIsland.txt: Earthquakes of Magnitude >=2, starting 2018-04-15, 00:00:00.00 UTC, within the full island of Hawai'i

For both of the .txt files, the 1st column is the origin date and the 2nd column is the origin time (UTC). The 3rd column is the earthquake magnitude. 
