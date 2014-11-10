#Ryan Batt
#25-April-2013
#Adapted from Cascade manual
#calculate the biomass for each species or taxon of zooplankton listed in the NTL LTER database
#this took forever
#I had to "guess" for a lot of the taxa, especially the rotifers.  But I was quite diligent.
#This is a very poorly-written function.  It's in this format b/c I initially just copied some code out of the the cascade manual, and did a bunch of "find" then "replace all", followed by a redefinition of the taxa to turn it into R code that could be used with the LTER database.  This will be extremely slow, b/c each if() statement must be evaluated each time.  When applied to a data frame with 10's of thousands of rows, this could really slow things down.
	#One way to speed this up would be to offer to return the evaluation of inwt (or just evaluate the expression, and forget the inwt function) if the if() statement is satisfied.

#v0 (25-Apr-2013) Just got all of the equations in there
#v1 (26-Apr-2013) Combined names of taxa that shared a common equation.  Didn't redefine function based on a match (entirely removed the inwt() function), simply return the evaluation of the expression --- this way not every if( )statement is evaluate per call to ZoopMass

ZoopMass <- function(x){
	
	
	taxon <- unique(x[,"taxon"])
	stopifnot(length(taxon)==1)
	len <- x[,"avg_length"]
	
	
	# ====================
	# = (some) Predators =
	# ====================
	if(is.element(taxon, c("BYTHOTREPHES LONGIMANUS\xe5\xca"))) return(45.2*len - 287.2) #Branstrator 2005 J. Plankton Res. (Added by ryan batt)
	if(is.element(taxon, c("CHAOBORUS"))) return(exp(1.189*log(pi*0.25*len)-8.644)) #taken from Cascade, but I am assuming that the width is 1mm (that's where 0.25 comes from)


	# ===============
	# = Cladocerans =
	# ===============
	if(is.element(taxon, c("CLADOCERAN","DAPHNIA","DAPHNIA AMBIGUA", "DAPHNIA DENTIFERA", "DAPHNIA DUBIA","DAPHNIA LONGIREMIS","DAPHNIA MENDOTAE","DAPHNIA PARVULA","DAPHNIA PULICARIA","DAPHNIA RETROCURVA", "OPHRYOXUS GRACILIS", "SCAPHOLEBERIS", "SIMOCEPHALUS", "DAPHNIA JUVENILE","DAPHNIA MENDOTAE JUVENILE","DAPHNIA PARVULA JUVENILE","DAPHNIA PULICARIA JUVENILE", "DAPHNIA RETROCURVA JUVENILE", "CERIODAPHNIA","CERIODAPHNIA LACUSTRIS"))) return(exp(1.9445+(2.72*(log(len)))))#DPUL; I have no idea what type of "CLADOCERAN", so I am just putting it here; couldn't find "OPHRYOXUS GRACILIS", put it here b/c cladoceran; "SCAPHOLEBERIS", put it here b/c family Daphniidae; "SIMOCEPHALUS", put it here b/c Daphniidae

	if(is.element(taxon, c("BOSMINA LONGIROSTRIS", "BOSMINIDAE", "EUBOSMINA COREGONI", "SINOBOSMINA FRYEI", "ACROPERUS HARPAE", "CAMPTOCERCUS RECTIROSTRIS","CHYDORUS","CHYDORUS SPHAERICUS", "PLEUROXUS PROCURVUS", "ALONA", "ALONA KARUA", "LEYDIGIA", "LEYDIGIA ACANTHOCERCOIDES", "LEYDIGIA QUADRANGULARIS")))   return(exp(2.7116+(2.5294*(log(len))))) #BOS; "SINOBOSMINA FRYEI" b/c 'bosmina' in name; "PLEUROXUS PROCURVUS" is in the family Chydoridae; I think that "LEYDIGIA" belongs here... it belongs to the family "Chydoridae" and the subfamily "Alonnae";
	# rem use same as Bosmina (Soranno, 1989)
	# rem alternative: return( exp(4.543+(3.636*(Log(len))))
	# rem whose source is Chydorus sphaericus in Rosen 1981 (D&R)
	# rem problem is, the Rosen formula gives 30-40 ug/l Alona
	
	if(is.element(taxon, c("DIAPHANOSOMA","DIAPHANOSOMA BIRGEI", "DIAPHANOSOMA BIRGEI JUVENILE"))) return(exp(1.6242+(3.0468*(log(len)))))
	# rem source is Diaphanosoma brachyrum in Bottrell et al. 1976 (D&R)
	# if(taxon=='SCRYS') return( exp(2.0539+(2.189*(log(len))))
	# rem source is Sida crystallina in Bottrell et al. 1976 (D&R)
	
	if(is.element(taxon, c("HOLOPEDIUM", "HOLOPEDIUM GIBBERUM", "ILYOCRYPTUS SPINIFER", "LEPTODORA KINDTI", "POLYPHEMUS", "POLYPHEMUS PEDICULUS"))) return(9.86*((len)^2.1)) #there is only 1 "ILYOCRYPTUS SPINIFER" in the whole damn data set, and I just spent 20 minutes trying to find the equation for its mass. This is the closest I could find, for now.  Also, I think that this species also goes by the name of a Polyphemus species... according to Wiki.
	# rem source is Peters and Downing (1984) general zooplankton regression
	# rem did NOT use formulas from Downing and Rigler (1984)
		
	
	# REM ------------------------------------------------------------------------
	# REM CALANOID COPEPODS
	# REM ------------------------------------------------------------------------
	if(is.element(taxon, c("AGLAODIAPTOMUS CLAVIPES", "CALANOID", "CALANOIDA COPEPODITES", "DIAPTOMID", "LEPTODIAPTOMUS MINUTUS", "LEPTODIAPTOMUS SICILIS", "LEPTODIAPTOMUS SICILOIDES", "SKISTODIAPTOMUS ", "SKISTODIAPTOMUS OREGONENSIS", "SKISTODIAPTOMUS PALLIDUS", "EPISCHURA LACUSTRIS"))) return(exp(1.2431+(2.2634*(log(len))))) #DIAPT
	# rem source is Diaptomus gracilis in Bottrell et al. 1976 (D&R)
	# if(is.element(taxon, c("EPISCHURA LACUSTRIS"))) return( exp(1.2431+(2.2634*(log(len))))
	# rem use same as Diaptomus (Soranno 1989)
	
	
	# REM ------------------------------------------------------------------------
	# REM CYCLOPOID COPEPODS
	# REM ------------------------------------------------------------------------
	if(is.element(taxon, c("DIACYCLOPS","DIACYCLOPS THOMASI", "ERGASILUS", "EUCYCLOPS", "EUCYCLOPS AGILIS","EUCYCLOPS ELEGANS", "EUCYCLOPS SERRULATUS", "MESOCYCLOPS EDAX"))) return(exp(1.6602+(3.968*(log(len))))) #couldn't find equation for "ERGASILUS", put it with the other cyclopoids, but this genus is very large (mm range)
	# rem source is Mesocyclops edax in Rosen 1981 (D&R)
	if(is.element(taxon, c("ACANTHOCYCLOPS", "ACANTHOCYCLOPS VERNALIS","CYCLOPOID", "HARPACTICOID", "ORTHOCYCLOPS MODESTUS", "PARACYCLOPS FIMBRIATUS POPPEI", "TROPOCYCLOPS PRASINUS MEXICANUS"))) return(exp(2.0577+(2.553*(log(len))))) #"HARPACTICOID" is just an order of copepod--- I have no idea what species this is supposed to be.
	# if(is.element(taxon, c("ORTHOCYCLOPS MODESTUS", "PARACYCLOPS FIMBRIATUS POPPEI"))) return( exp(2.0577+(2.553*(log(len)))) #didnt' know what to do with "PARACYCLOPS FIMBRIATUS POPPEI", so put it here b/c cyclops, and "para" has to be similar to "ortho", right?
	# rem source is Cyclops vicinus in Botrell et al. 1976 (D&R)
	
	
	# REM ------------------------------------------------------------------------
	# REM IMMATURE COPEPODS
	# REM ------------------------------------------------------------------------
	if(is.element(taxon, c("COPEPODITES", "CYCLOPOIDA COPEPODITES"))) return(exp(2.0577+(2.553*(log(len))))) #COPD
	# rem source is Cyclops vicinus in Botrell et al. 1976 (D&R)
	if(is.element(taxon, c("COPEPOD NAUPLII"))) return(exp(0.6977+(0.469*(log(len)))))
	# rem source is copepod nauplii in Rosen 1981 (D&R)
	

	# REM ------------------------------------------------------------------------
	# REM ROTIFERS (all in D&R 1984, p.247-249, column 3)
	# REM ------------------------------------------------------------------------
	if(is.element(taxon, c("ASPLANCHNA", "ASPLANCHNA BRIGHTWELLI", "ASPLANCHNA HERRICKII", "ASPLANCHNA PRIODONTA"))) return(100*0.23*((len)^3)) #ASPLG
	if(is.element(taxon, c("ANURAEOPSIS", "ANURAEOPSIS FISSA", "BRACHIONUS", "BRACHIONUS ANGULARIS", "BRACHIONUS BIDENTATA", "BRACHIONUS CALYCIFLORUS", "BRACHIONUS CAUDATUS", "BRACHIONUS HAVANAENSIS", "BRACHIONUS QUADRIDENTATUS", "BRACHIONUS RUBENS", "BRACHIONUS URCEOLARIS", "MANFREDIUM", "MANFREDIUM EUDACTYLOTUM", "NOTHOLCA", "NOTHOLCA ACUMINATA", "NOTHOLCA FOLIACEA", "NOTHOLCA LABIS", "NOTHOLCA MICHIGANENSIS", "NOTHOLCA SQUAMULA", "ASCOMORPHA", "ASCOMORPHA ECAUDIS", "ASCOMORPHA OVALIS", "ASCOMORPHA SALTANS", "SYNCHAETA"))) return(100*0.12*((len)^3)) #BANG; couldn't find "MANFREDIUM EUDACTYLOTUM" or the genus... so I found that it belonged to the family Brachionidae and put it here; "NOTHOLCA" belongs to the family Brachionidae
	# rem Brachionus
	
	if(is.element(taxon, c("EUCHLANIS", "ENCENTRUM", "EUCHLANIS", "EUCHLANIS PELLUCIDA", "NOTOMMATA", "TRICHOTRIA", "TRICHOTRIA TETRACTIS", "TROCHOSPHAERA SOLSTITIALIS", "LECANE", "LECANE CLARA", "LECANE INERMIS", "LECANE LEONTINA", "LECANE LUNA", "LECANE MIRA", "LECANE TENUISETA", "LECANE TUDICOLA", "MONOSTYLA", "MONOSTYLA BULLA", "MONOSTYLA CRENATA", "MONOSTYLA LUNARIS", "MONOSTYLA OBTUSA", "MONOSTYLA QUADRIDENTATA", "MONOSTYLA STENROOSI"))) return(100*0.1*((len)^3)) #Couldn't find an equation for ENCENTRUM. Picked somethign that looked (somewhat) similar in google images; couldn't find "NOTOMMATA", so I put it here b/c it looked kinda similar; couldn't find "TRICHOTRIA", but looked similar to euchlanis (tail with two spines); "TROCHOSPHAERA SOLSTITIALIS" was almost unsearchable, but the fiew pictures I could find either looked perfectly spherical, or like euchlanis.I think "MONOSTYLA" is a genus that belongs to the family Lecanidae, plus the google images of "lecane rotifer" looks similar to that of "MONOSTYLA"
	# rem Euchlanis
	# rem Lecane (assumed same as Euchlanis)
	
	if(is.element(taxon, c("FILINIA", "FILINIA TERMINALIS"))) return(100*0.13*((len)^3))
	# rem Filinia
	
	if(is.element(taxon, c("POMPHOLYX", "POMPHYLOX", "POMPHYLOX SULCATA", "TESTUDINELLA REFLEXA", "GASTROPUS","GASTROPUS HYPTOPUS", "GASTROPUS STYLIFER"))) return(100*0.2*((len)^3))# couldn't find "POMPHOLYX", so it grouped them with gastropus; couldn't find "TESTUDINELLA REFLEXA", but popmphylox is in same family (testudinellidae).. and they look similar.
	# rem Gastropus (same as historical)
	# rem Conochiloides (assumed same as Gastropus, Soranno 1989)
	# rem Synchaeta (assumed same as Gastropus, Soranno 1989)
	
	if(is.element(taxon, c("KELLICOTTIA", "KELLICOTTIA BOSTONIENSIS", "KELLICOTTIA LONGISPINA", "MONOMMATA"))) return(100*0.03*((len)^3)) #couldn't find "MONOMMATA", put it here b/c I thought it looked similar
	# rem Kellicotia
	
	if(is.element(taxon, c("KERATELLA", "KERATELLA COCHLEARIS", "KERATELLA COCHLEARIS F. TECTA", "KERATELLA CRASSA", "KERATELLA EARLINAE", "KERATELLA HIEMALIS", "KERATELLA QUADRATA", "KERATELLA SERRULATA", "KERATELLA TAUROCEPHALA", "KERATELLA TESTUDO", "KERATELLA TICINENSIS", "LOPHOCHARIS OXYSTERNON", "MYTILINA"))) return(100*0.02*((len)^3)) #couldn't find "LOPHOCHARIS OXYSTERNON", so I'm putting it ehre b/c I thought it looked similar to Keratella; "MYTILINA" is a family that has species in it belonging to Lophocharis, so I put it here.
	# rem Keratella
	
	if(is.element(taxon, c("PLOESOMA", "PLOESOMA HUDSONI", "PLOESOMA LENTICULARE", "PLOESOMA TRUNCATUM"))) return(100*0.15*((len)^3))#PLOES
	# rem Ploesoma (average of P.hudsonii and P. triacanthum)
	
	if(is.element(taxon, c("POLYARTHRA", "POLYARTHRA DOLICHOPTERA", "POLYARTHRA EURYPTERA", "POLYARTHRA MAJOR", "POLYARTHRA REMATA", "POLYARTHRA VULGARIS"))) return(100*0.28*((len)^3))
	# rem Polyarthra (same as historical)
	# REM species where we should measure length and width : assign average inwt
	
	if(is.element(taxon, c("CEPHALODELLA","COLLOTHECA", "COLLOTHECA MUTABILIS", "COLLOTHECA PELAGICA"))) return(0.1) #sort of a guess for Rotaria; taken from Table 1 in Ludovisi & Jørgensen 2009 Ecological Modelling; OK, I'm officially making this the "unknown rotifer" group; actually, this is the real number for Collotheca, Bottrell et al. 1976 and Kobayashi et al. 1996 Mar. Freshwater Res.
	
	if(is.element(taxon, c("COLURELLA"))) return(0.003) #Makarewicz and Likens 1979 and Kobayashi et al. 1996 Mar. Freshwater Res.
	
	if(is.element(taxon, c("HEXARTHRA"))) return(0.85) #"HEXARTHRA" in Dumont et al. 1975 and Kobayashi et al. 1996 Mar. Freshwater Res.
	
	if(is.element(taxon, c("LEPADELLA", "LEPADELLA ACUMINATA", "LEPADELLA PATELLA", "LEPADELLA TRIPTERA"))) return(0.15) #"LEPADELLA" in Dumont et al. 1975 and Kobayashi et al. 1996 Mar. Freshwater Res.
	
	if(is.element(taxon, c("CONOCHILOIDES", "CONOCHILOIDES NATANS","CONOCHILUS", "CONOCHILUS HIPPOCREPIS", "CONOCHILUS UNICORNIS"))) return(0.12) #Makarewicz and Likens 1979 and Kobayashi et al. 1996 Mar. Freshwater Res. (Didn't use Cascade # b/c that's colony)
	
	if(is.element(taxon, c("TRICHOCERCA", "TRICHOCERCA BIROSTRIS","TRICHOCERCA CAPUCINA","TRICHOCERCA CYLINDRICA","TRICHOCERCA LONGISETA","TRICHOCERCA MULTICRINIS","TRICHOCERCA PORCELLUS", "TRICHOCERCA PUSILLA","TRICHOCERCA ROUSSELETI", "TRICHOCERCA SIMILIS"))) return(0.068) #TMULT, TMUL
	# rem Trichocerca; 1984-1990 long-term means


	if(is.element(taxon, c("UNIDENTIFIED","UNKNOWN", "UNKNOWN EGG-SHAPED ROTIFER","UNKNOWN ROTIFER"))) return(0)
	
}

