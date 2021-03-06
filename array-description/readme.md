# array-description - readme

Format of files in this directory (running `head`).

**OID40453_annotation.gff**

```
C16582	TJGR_v9	transcripts	35	385	-1	.	.	ID=CGI_10000001; color=000000
C17212	TJGR_v9	transcripts	31	363	+1	.	.	ID=CGI_10000002; color=000000
C17316	TJGR_v9	transcripts	30	257	+1	.	.	ID=CGI_10000003; color=000000
C17476	TJGR_v9	transcripts	34	257	-1	.	.	ID=CGI_10000004; color=000000
C17998	TJGR_v9	transcripts	196	387	-1	.	.	ID=CGI_10000005; color=000000
C18346	TJGR_v9	transcripts	174	551	+1	.	.	ID=CGI_10000009; color=000000
C18428	TJGR_v9	transcripts	286	546	-1	.	.	ID=CGI_10000010; color=000000
C18964	TJGR_v9	transcripts	203	658	-1	.	.	ID=CGI_10000011; color=000000
C18980	TJGR_v9	transcripts	30	674	+1	.	.	ID=CGI_10000012; color=000000
C19100	TJGR_v9	transcripts	160	681	-1	.	.	ID=CGI_10000013; color=000000
```

**OID40453_probe_locations.gff**

```
scaffold1	NGS	OID40453_probe_locations	52347	52396	1	color=000000
scaffold1	NGS	OID40453_probe_locations	52472	52521	1	color=000000
scaffold1	NGS	OID40453_probe_locations	52597	52649	1	color=000000
scaffold1	NGS	OID40453_probe_locations	52741	52797	1	color=000000
scaffold1	NGS	OID40453_probe_locations	52851	52915	1	color=000000
scaffold1	NGS	OID40453_probe_locations	52996	53056	1	color=000000
scaffold1	NGS	OID40453_probe_locations	53101	53152	1	color=000000
scaffold1	NGS	OID40453_probe_locations	53241	53290	2	color=000000
scaffold1	NGS	OID40453_probe_locations	53351	53405	1	color=000000
scaffold1	NGS	OID40453_probe_locations	53476	53534	1	color=000000
```

**OID40453_probes.txt**  (zip version)

```
PROBE_SEQUENCE	PROBE_ID	SEQ_ID	POSITION	DESIGN_NOTE	SELECTION_CRITERIA	PROBE_CLASS	PROBE_TM	FREQUENCY	COUNT	RULES	SLOPE	RSQUARED	GENOMIC	COT1	SNP_COUNT	SCORE	RANK
CGATTAGTGGACGTGAGAATTCTCAAGAGCACACCGTGTGTATACATCAC	scaffold39960:+:29315	scaffold39960:29300-48043	29315		target_tm:76.00;probe_tm:76.80;freq: 3.28;count:01;rules:0000;score:0991	experimental	76.80	 3.28	   0	0.00	0.00	0.00	0.00	0	 991	1
TTTGTTGGAGCTGTATAATAAATGAAACTGAAACCAATGAATTCTTATTGAATTCAAC	scaffold39960:+:29445	scaffold39960:29300-48043	29445		target_tm:76.00;probe_tm:70.90;freq:13.30;count:01;rules:0000;score:0804	experimental	70.90	13.30	1	   0	0.00	0.00	0.00	0.00	0	 804	1
CCTTTGCCTTTTTTTCAAAAATGCAATGCTTGTCAAGATGTTTTAGATGAAGCTGTGG	scaffold39960:+:29572	scaffold39960:29300-48043	29572		target_tm:76.00;probe_tm:74.40;freq:19.66;count:01;rules:-350;score:0496	experimental	74.40	19.66	1	-350	0.00	0.00	0.00	0.00	0	 496	1
CATTTCCTATTTCCTCATCTTAATAATAGCCAAGGGTTTTTAGTCCACTTTGCTC	scaffold39960:+:29677	scaffold39960:29300-48043	29677		target_tm:76.00;probe_tm:73.90;freq: 6.46;count:01;rules:0000;score:0916	experimental	73.90	 6.46	   0	0.00	0.00	0.00	0.00	0	 916	1
AAATGCCTTTGATTTTGGATCTTTATATAACAAATACAATTCATAATTTGACTGTT	scaffold39960:+:29807	scaffold39960:29300-48043	29807		target_tm:76.00;probe_tm:68.70;freq:24.29;count:04;rules:0000;score:0417	experimental	68.70	24.29	4	   0	0.00	0.00	0.00	0.00	0	 417	1
TTTTAGTCATATATTTATACAACTAAGGTGAAAGAAAAGACTAAAGTTTCCT	scaffold39960:+:29947	scaffold39960:29300-48043	29947		target_tm:76.00;probe_tm:68.60;freq: 8.53;count:01;rules:0000;score:0790	experimental	68.60	 8.53	   0	0.00	0.00	0.00	0.00	0	 790	1

```



**OID40453_stats.txt** 

```
SEQ_ID	LENGTH	PROBES	DENSITY	MEAN_INTERVAL	MEDIAN_INTERVAL	1ST_QUARTILE	3RD_QUARTILE	MIN_INTERVAL	MAX_INTERVAL	COVERAGE
scaffold1:52347-57631	5285	43	122	125	125	120	132	105	150	2384
scaffold1:161989-185527	23539	171	137	125	125	120	140	105	451	9819
scaffold1:207337-231800	24464	176	139	125	125	115	137	105	813	9926
scaffold3:65419-71497	6079	47	129	126	127	115	136	105	200	2660
scaffold3:72058-78354	6297	43	146	126	130	115	140	105	602	2511
scaffold3:82418-120562	38145	248	153	124	125	115	140	105	1071	14107
scaffold3:122297-131079	8783	66	133	125	125	120	135	102	557	3814
```


---

**121101_OysterV9_MG_meth.pos**

```
PROBE_ID	SEQ_ID	CHROMOSOME	POSITION	COUNT	LENGTH	GC
scaffold1:+:52347	scaffold1:52347-57631	scaffold1	52347	1	50	0.44
scaffold1:+:52472	scaffold1:52347-57631	scaffold1	52472	1	50	0.44
scaffold1:+:52597	scaffold1:52347-57631	scaffold1	52597	1	53	0.28
scaffold1:+:52741	scaffold1:52347-57631	scaffold1	52741	1	57	0.30
scaffold1:+:52851	scaffold1:52347-57631	scaffold1	52851	1	65	0.34
scaffold1:+:52996	scaffold1:52347-57631	scaffold1	52996	1	61	0.38
scaffold1:+:53101	scaffold1:52347-57631	scaffold1	53101	1	52	0.42
scaffold1:+:53241	scaffold1:52347-57631	scaffold1	53241	2	50	0.46
scaffold1:+:53351	scaffold1:52347-57631	scaffold1	53351	1	55	0.40
```

**121101_OysterV9_MG_meth.ProbeAnnotation.rda**

R object - no preview- annotation object.


**oyster.v9.glean.final.rename.mRNA.gff**

```
C16582	GLEAN	mRNA	35	385	0.555898	-	.	ID=CGI_10000001;
C17212	GLEAN	mRNA	31	363	0.999572	+	.	ID=CGI_10000002;
C17316	GLEAN	mRNA	30	257	0.555898	+	.	ID=CGI_10000003;
C17476	GLEAN	mRNA	34	257	0.998947	-	.	ID=CGI_10000004;
C17998	GLEAN	mRNA	196	387	1	-	.	ID=CGI_10000005;
C18346	GLEAN	mRNA	174	551	1	+	.	ID=CGI_10000009;
C18428	GLEAN	mRNA	286	546	0.555898	-	.	ID=CGI_10000010;
C18964	GLEAN	mRNA	203	658	0.999572	-	.	ID=CGI_10000011;
C18980	GLEAN	mRNA	30	674	0.555898	+	.	ID=CGI_10000012;
C19100	GLEAN	mRNA	160	681	0.999955	-	.	ID=CGI_10000013;
```

**2013.09.24.targets.txt**

```
FileNameCy5	FileNameCy3
A01_EE2.input_532.pair	A01_Ctrl.input_635.pair
A02_EE2.MBD_532.pair	A02_Ctrl.MBD_635.pair
A03_EE2.MBD_635.pair	A03_Ctrl.MBD_532.pair
```

**spottypes.txt**

```
SpotType	GENE_EXPR_OPTION	PROBE_ID	Color
BLOCK1	BLOCK1	*	black
RANDOM	RANDOM	*	grey

```


