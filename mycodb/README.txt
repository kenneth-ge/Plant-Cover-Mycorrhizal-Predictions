Title: MycoDB, a global database of plant response to mycorrhizal fungi
Contact: bala.chaudhary@depaul.edu
=====================================================================================
MycoDB Version 1 – published April 22, 2016

See original publication, specifically Table 1, for full meta-data and description of variables:
Chaudhary VB, Rúa MA, Antoninka A, Bever JD, Cannon J, Craig A, Duchicela J, Frame A, Gardes M, Gehring C, Ha M, Hart M, Hopkins J, Ji B, Johnson NC, Kaonongbua W, Karst J, Koide RT, Lamit LJ, Meadow J, Milligan BG, Moore JC, Pendergast IV TH, Piculell B, Ramsby B, Simard S, Shrestha S, Umbanhowar J, Viechtbauer W, Walters L, Wilson GWT, Zee PC, Hoeksema JD (2016) MycoDB, a global database of plant response to mycorrhizal fungi. Scientific Data 3: 160028. https://doi.org/10.1038/sdata.2016.28

Associated files FungalTree_version1 and PlantTree_version1 can be used to conduct phylogenetic meta-analyses with MycoDB_version1.

=====================================================================================
MycoDB Version 2 – published April 28, 2017

* Version 2 has two additional variables (columns):
o “Rua2017” – indicates whether a study was used in a 2017 meta-analysis conducted by Rua et al. 20XX.
o “PlantSpeciesMar2017” – Species names for ectomycorrhizal plant species 
* Additional ectomycorrhizal studies (573 rows) were added to increase sample size for studies published in Rua et al. 2017 (in review) and Karst et al. 2018 (in review). Full details for the data amendments made in MycoDB_version2 can be found in these publications.

=====================================================================================
MycoDB Version 3 – published November 3, 2017

* Fixed an issue with the NONCTLTRTSETID and CTLTRTSETID columns being swapped for entries associated with EXPERIMENTIDs: 117, 144, 146, 206, 218, 221, 226, 232, 235, 302, 322, 327, 387, 388, 389, 391, 392, 393, 394, 395, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 441, 442, 443, 444, 445, 446, 447, 448, 450, 455, 456, 457, 458, 459, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491,
492, 493, 494, 495, 496, 498, 499, 500, 501, 502, 503, 504, 505, 506, 508, 509, 510, 511, 512, 513, 514, 517, 518, 519, 521, 522, 523, 524, 525, 526, 527, 528, 530, 531, 532, 539, 540, 639, 640, 645, 667, 678, 683, 684, 690, 691, 726, 744, 756, 774, 782, 784, 788, 908, 914, 916, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 977, 1001, 1002, 1047, 1048, 1097, 1098, 1174, 1179, 1194, 1242, 1245, 1287, 1300, 1308, 1337, 1344, 1346, 1348, 1349, 1350, 1385, 1425, 1426, 1433, 1434, 1452, 1541, 1453, 1454, 1455, 1456, 1457, 1458, 1459, 1460, 1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468, 1469, 1470, 1471, 1472, 1473, 1474, 1475, 1476, 1477, 1478, 1479, 1480, 1481, 1482, 1483, 1484, 1485, 1486, 1487, 1488, 1489, 1490, 1491, 1492, 1493, 1494, 1495, 1496, 1497, 1498, 1499, 1500, 1501, 1502, 1503, 1504, 1505, 1506, 1507, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1517, 1518, 1519, 1520, 1521, 1522, 1523, 1524, 1525, 1526, 1527, 1528, 1529, 1530, 1530, 1531, 1532, 1533, 1534, 1535, 1536, 1537, 1538, 1539, 1540
* Changed four mentions of eucalyptus_europhylla to eucalyptus_urophylla in PlantSpeciesMar2017 for NONCTLTRTSETIDs 2733, 2735, 2736, 2737
* Changed three mentions of pinus_cornata to pinus_contorta in PlantSpeciesMar2017 for NONCTLTRTSETIDs 5660, 5661, 5662

=====================================================================================

MycoDB Version 4 – published July 10, 2018

* One new variable (column) entitled “Karst2018” was added to version 4 in association with the publication by Karst et al. 2018. This column indicates with a ‘yes’, which studies were used in Karst J, Burns C, Cale JC, Antunes PM, Woods M, Lamit LJ, Hoeksema JD, Zabinski C, Gehring CA, La Flèche M, Rúa M. 2018. Tree species with limited geographic ranges show extreme responses to ectomycorrhizas. Global Ecology and Biogeography. 
* Five new variables (columns) were added to version 4 in associated with the publication by Hoeksema et al. 2018. Evolutionary history predicts the strength of mycorrhizal mutualism Communications Biology:
o "PlantSpecies2018" = Plant species names used in the analyses by Hoeksema et al. (2018). These names include 2 genus spelling corrections (Cytharexylum changed to Citharexylum and Macherium changed to Machaerium). Compared to PlantSpeciesMar2017, it additionally includes names for arbuscular mycorrhizal (AM) plants, but does not include a small number of taxonomy updates incorporated into PlantSpeciesMar2017 (Afzelia changed to Intsia, Uapaca somon changed to Uapaca togoensis, Hopea helferi changed to Vatica helferi, Shorea seminis changed to Hopea seminis, and Salix dasyclados changed to Salix gmelinii), in order to maintain compatibility with the companion plant phylogeny published with MycoDB Version 4. Those taxonomy updates will be included in any future versions of the database and companion plant phylogeny.
o "FungalGenus2018" = Fungal genus names used in the analyses by Hoeksema et al. (2018), updated based on recent taxonomic studies.
o "EMF_ORIGIN_TED" = The "Plant Origin" variable used in analyses of ectomycorrhizal (EM) symbiosis by Hoeksema et al. (2018). Its construction is described in the Supplementary Methods of that paper. "TED" is used here to indicate that the coded origins correspond to those hypothesized by Tedersoo et al. (2010, Mycorrhiza 20:217-263).
o "EMF_ORIGIN_ALT" = Same as EMF_ORIGIN_TED except with a single shared EM origin for the fungal genera Hydnotrya and Tuber. 
o "Hoeksema2018" - indicates whether or not an observation was or was not used in analyses for the paper by Hoeksema et al. (2018).
* For ectomycorrhizal (EM) plant genus Uapaca (4 observations, NONCTLTRTSETIDs 3505-3508), changed PlantFamily from euphorbiaceae to phyllanthaceae based on updated information on systematics
* Removed two observations (NONCTLTRTSETIDs 3744 and 4044), one each for fungal genera Leucopaxillus and Protubera, because there is some question of whether they were correctly identified or are actually mycorrhizal
* Replaced duplicate instance of NONCTLTRTSETID (15854, first author Lukesova) with a new unique value (16428), and replaced NONCTLTRTSETIDs 15855-15858 (first author Sousa) with new values (16429-16432)
* Replaced duplicate instances of CTLTRTSETID (15852-15855) with new unique values (16053-16056)
* Associated files FungalTree_version2.tre and PlantTree_version2.tre represent fungal and plant phylogenies, respectively, used in the phylogenetic meta-analysis models of Hoeksema et al. (2018). They are similar to FungalTree_version1 and PlantTree_version1, published originally with MycoDB Version 1 (Chaudhary et al. 2016), except with fungal genus and plant species names updated to match the new variables PlantSpecies2018 and FungalGenus2018. In addition, PlantTree_version2 corrects a tree topology error by moving the plant genus Shorea into the family Dipterocarpaceae.
* The associated file phylometa_2018.R is a .R script containing original R code for analyses conducted in association with Hoeksema et al. (2018).
* Below is a key linking names of variables in MycoDB_version4.csv (besides the 5 new variables listed above) with the variable names used in the paper by Hoeksema et al. (2018). See Chaudhary et al. (2016) for explanation of how these variables were constructed.
o name in MycoDB_version 4 = name in Hoeksema et al. 2018
o EFFECTSIZE1 = LRR, log response ratio
o NONCTLTRTSETID = Study Id
o CTLTRTSETID = Control Set
o PAPERTITLE = Paper
o FERTN = N-fertilization
o FERTP = P-fertilization
o STERILIZED = Sterilization
o NONMYCOCONTROL2 = Microbial Control
o LOCATION = Location
o FUNGROUP = Plant Functional Group
o PLANTLIFEHISTORY = Plant Life History
o DOMESTICATED = Domestication
o INOC.COMPLEXITY = Inoculum Complexity
