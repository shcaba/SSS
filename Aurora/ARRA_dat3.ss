# For Base 2 reduced Ninput for discard lengths at power of 0.9 which reduces largest by more than half (1330 to 648)
# For data9 - fixed survey to correct scale.
# For dat2 - added research catch to NODISC 
#AURORA	ROCKFISH	Draft 6-8-2013																																																																																																																													
########################################																																																																																																																																																																																																																																																																				
###	Global	model	specifications	###																																																																																																																														
1916	#	Start	year																																																																																																																															
2012	#	End	year																																																																																																																															
1	#	Number	of	seasons/year																																																																																																																														
12	#	Number	of	months/season																																																																																																																														
1	#	Spawning	occurs	at	beginning	of	season																																																																																																																											
1	# Number of surveys
1	# Number of areas
1	#	Number	of	areas																																																																																																																														
FISHERY%SURVEY1
0.5	0.5	#Timing	of	each	fishery/survey																																																																																																																								
1	1	#Area	of	each	fleet																																																																																																																								
1	#_units	of	catch:	1=bio;	2=num																																																																																																																											
0.05	#_se	of	log(catch)	only	used	for	init_eq_catch	and	for	Fmethod	2	and	3;	use	-1	for	discard	only	fleets																																																																																																													
2	#	Number	of	genders																																																																																																																														
80	#	Number	of	ages	in	population	dynamics																																																																																																																																																																																																																																																								
###	Catch	section	###																																																																																																																															
0	#	Initial	equilibrium	catch	(landings	+	discard)	by	fishing	fleet																																																																																																																						
97 # Number of lines catch data				
#_catch_biomass(mtons):_columns_are_fisheries,year,season				
0.0875782	1916	1
0.1370728	1917	1
0.1551165	1918	1
0.1062395	1919	1
0.109968	1920	1
0.0913202	1921	1
0.0794708	1922	1
0.0884296	1923	1
0.0601049	1924	1
0.0710968	1925	1
0.1133513	1926	1
0.090444	1927	1
0.1228706	1928	1
0.1472441	1929	1
0.1582176	1930	1
0.1643209	1931	1
0.19309756	1932	1
0.2647699	1933	1
0.20145395	1934	1
0.15266151	1935	1
0.1447987	1936	1
0.24798299	1937	1
0.3745089	1938	1
0.5515662	1939	1
0.5611691	1940	1
1.1021307	1941	1
0.4489398	1942	1
1.0014008	1943	1
1.79704704	1944	1
3.55914729	1945	1
2.9165016	1946	1
2.7821698	1947	1
2.517976	1948	1
1.6648344	1949	1
2.2804121	1950	1
3.5440784	1951	1
3.8889692	1952	1
4.3158926	1953	1
2.6722903	1954	1
2.36099582	1955	1
2.9710364	1956	1
3.1570328	1957	1
4.66686915	1958	1
5.29585379	1959	1
4.0263187	1960	1
2.6685304	1961	1
2.23983059	1962	1
2.4412832	1963	1
1.845031	1964	1
2.0231105	1965	1
3.406657	1966	1
1.9318665	1967	1
2.3184581	1968	1
2.60207082	1969	1
3.85242425	1970	1
8.66301663	1971	1
9.70111103	1972	1
19.0262455	1973	1
12.240366	1974	1
14.5162418	1975	1
15.2537305	1976	1
9.278977	1977	1
3.93490078	1978	1
23.7065069	1979	1
15.248945	1980	1
12.355252	1981	1
55.11966108	1982	1
143.348549	1983	1
37.899271	1984	1
70.515321	1985	1
108.398169	1986	1
47.6022677	1987	1
139.4043527	1988	1
147.2777	1989	1
210.228179	1990	1
60.570343	1991	1
215.557338	1992	1
146.814125	1993	1
106.373509	1994	1
75.072645	1995	1
56.435466	1996	1
53.437661	1997	1
42.538371	1998	1
19.5764422	1999	1
38.348504	2000	1
30.937276	2001	1
56.31229	2002	1
76.65072	2003	1
85.01963	2004	1
54.361151	2005	1
48.340125	2006	1
53.465397	2007	1
48.402131	2008	1
51.934603	2009	1
40.980279	2010	1
28.577187	2011	1
42.713778	2012	1
																																																																																																																																																																																																																																			
2	#Number	of	index	observations																																																																																																																													
#Units:	0=numbers,1=biomass,2=F; Errortype: -1=normal,0=lognormal,>0=T																																																																																																																														
#Fleet	Units	Errortype																																																																																																																															
1 1 0 # fleet 1: FISHERY
2 1 0 # fleet 2: SURVEY
#_year seas index obs se(log)
1916 1 2 1 0.00001 # SURVEY1
2000 1 2 0.1 0.00001 # DEPLETION
																																																																																																																																																																																																																																		
0 #_N_fleets_with_discard
0 #_N_discard_obs

0 #_N_meanbodywt_obs
30 #_DF_meanwt

## Population size structure
1 # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector
# binwidth for population size comp
# minimum size in the population (lower edge of first bin and size at age 0.00)
# maximum size in the population (lower edge of last bin)

-1 #_comp_tail_compression
1e-007 #_add_to_comp
0 #_combine males into females at or below this bin number
#																																																																																																																																		
16	#_N_LengthBins																																																																																																																																	
#	Data	length	bins																																																																																																																															
8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38																																																																																																																			
#																																																																																																																																		
0	#_N_Length_obs																																																																																																																																	
#TWL	(N=35),	Females	then	Males																																																																																																																														
#Year	Seas	Fleet	Gender	Part	Nsamp	F-8	F-10	F-12	F-14	F-16	F-18	F-20	F-22	F-24	F-26	F-28	F-30	F-32	F-34	F-36	F-38	M-8	M-10	M-12	M-14	M-16	M-18	M-20	M-22	M-24	M-26	M-28	M-30	M-32	M-34	M-36	M-38																																																																																													
#																																																																																																																																		
#Age	composition	set-up																																																																																																																																
61	#_N_age_bins																																																																																																																																	
0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60																																																																						
0 #_N_ageerror_definitions

0 #_N_Agecomp_obs
1 #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths
0 #_combine males into females at or below this bin number

0 #_N_MeanSize-at-Age_obs
0 #_N_environ_variables
0 #_N_environ_obs
0 # N sizefreq methods to read 
0 # no tag data 
0 # no morphcomp data 

999 # End data file

																																																																																																																																		
