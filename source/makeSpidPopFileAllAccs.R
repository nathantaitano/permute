library(stringr)

# See makeOaxSubpopTableForPgdSpider.R for creation of ../../gbsOutputData/subpopGroups_20170413.txt 
# Manually edited ../../gbsOutputData/subpopGroups_20170413.txt by adding mikey+ref to make ../../gbsOutputData/subpopGroups_20170510.txt
# Command (in terminal) to get mikey+ref from mikeyLmiss20Allhet05lmiss20imiss20Cm334ChiltZunla.vcf
# pcregrep "^#CHROM" mikeyLmiss20Allhet05lmiss20imiss20Cm334ChiltZunla.vcf
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1-1:8	100-1:20	100-2:19	102-1:18	103-1:17103-2:32	104-1:31	104-2:30	105-1:29	105-2:28	106-1:27	106-2:26	107-1:25	107-2:40	108-1:39	109-1:38	109-2:37	112-1:36	112-2:35	116-1:34	116-2:33	117-1:48	117-2:47	118-1:46	118-2:45	119-1:44	119-2:43	120-1:42	120-2:41	121-1:56	122-1:55	122-2:54	124-1:53	124-2:52	125-1:51	125-2:50	126-1:49	126-2:64	127-1:63	127-2:62	128-1:61	128-2:60	129-1:59	130-1:58	130-2:57	131-1:72	131-2:71	132-1:70	132-2:69	133-1:68	134-1:67	134-2:66	135-1:65	135-2:80	136-1:79	136-2:78	137-1:77	137-2:76	138-1:75	138-2:74	139-1:73	139-2:88	140-1:87	140-2:86	141-1:85	141-2:84	142-1:83	142-2:82	143-1:81	143-2:96	144-1:95	144-2:94	145-1:93	145-2:92	146-1:91	147-1:90	75-1:7	75-2:6	76-1:5	76-2:4	77-1:3	77-2:2	79-1:1	79-2:16	80-1:15	81-1:14	81-2:13	84-1:12	88-1:11	91-1:10	92-1:9	93-1:24	93-2:23	98-1:22	98-2:21	147-2:104	148-1:103	148-2:102	149-1:101	149-2:100	150-1:99	150-2:98	151-1:97	151-2:11152-1:111	152-2:110	153-1:109	153-2:108	155-1:107	155-2:106	156-1:105	156-2:120	157-1:11157-2:118	162-1:117	162-2:116	163-1:115	163-2:114	164-1:113	164-2:128	165-1:127	166-1:12166-2:125	168-1:123	168-2:122	169-1:121	169-2:136	170-1:135	170-2:134	173-1:133	173-2:13174-1:131	174-2:130	177-1:129	177-2:144	179-1:138	179-2:139	180-1:140	180-2:141	181-1:14181-2:143	182-1:137	182-2:145	183-1:146	183-2:147	184-1:148	184-2:149	185-1:150	185-2:15186-1:152	186-2:160	187-1:159	187-2:158	188-1:157	188-2:156	189-1:155	189-2:154	190-1:15190-2:168	191-1:167	191-2:166	200-1:165	200-2:164	202-1:163	203-1:162	203-2:161	204-1:17204-2:175	206-1:174	206-2:173	208-1:172	208-2:171	211-1:170	212-1:169	212-2:184	213-1:18213-2:182	214-1:181	214-2:180	215-1:179	215-2:178	216-1:177	216-2:192	217-1:191	217-2:19218-1:189	218-2:188	219-1:187	219-2:186	Apple_Pimento:193	Banana_Sweet:195	Big_Red:199	Bulgarian_Carrob:201	Bulgarian_Carrot:204	Bulls_heart_2:205	Buran:206	Calico:209	California_Mild:211	Chili_De_Arbol:217	China_Giant_Sweet:219	Chinese_Ching_Choo:220	Chocolate_Beauty:221	Colubian_Rainbow:226	Dulcetta_Orange:229	Feher_Ozon_Paprika:232	Filius_Blue:233	Mams_Biber:194	Marseilles_Sweet_Yellow:196	Mayan_Cobanero:198	Mulato:200	Pequin:207	Peter_Pepper_Red:208	Piment_Vegetarian:212	Pimento_de_Aeyde:210	Riot:215	Sandia:218	Shishito:222	Succette_de_Provence:223	Sweet_Chocolate:225	sunrise_eclipse:224	XRef:9997	chiltepin:9998	zunla:9999

forSpid <- read.table("../../gbsOutputData/subpopGroups_20170510.txt")
forSpid[,1] <- as.character(forSpid[,1])
forSpid[,2] <- as.character(forSpid[,2])
pcaGroups <- read.csv("../../gbsInputData/PCA_Groups_20160623.csv")
pcaGroups$pcaGrps <- as.character(pcaGroups$pcaGrps)
mikGroups <- read.csv("../../gbsInputData/mikeyOriginTable.csv")
mikGroups$Region <- as.character(mikGroups$Region)
mikGroups$Type <- as.character(mikGroups$Type)

for(r in 1:nrow(forSpid)){
  justNames <- str_extract(forSpid[,1], ".+(?=:[1-9])")
  mikMatch <- grepl(justNames[r], mikGroups$Type)
  if(any(mikMatch)){
    forSpid[r,2] <- mikGroups$Region[which(mikMatch)]
  }
  else next
}
forSpid[forSpid[,1] == "XRef:9997",2] <- "MesAm"
forSpid[forSpid[,1] == "chiltepin:9998", 2] <- "MesAm"
forSpid[forSpid[,1] == "zunla:9999", 2] <- "EAs"
forSpid[forSpid[,1] == "Bulgarian_Carrob:201", 2] <- "Europe"
forSpid[forSpid[,1] == "sunrise_eclipse:224", 2] <- "SoWes"

write.table(forSpid, "../../gbsOutputData/allSubpopGroups_20170510.txt", quote=F,
            col.names = F, row.names = F)

write(forSpid[,1], "../../gbsOutputData/allPassingKeepfile.txt")