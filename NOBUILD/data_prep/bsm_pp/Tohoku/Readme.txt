
HIGH RATE LEVEL 2 BOREHOLE STRAINMETER DATA

	For each event there is one tar file per strainmeter. Within each tar 
	there are five gzipped files, one for each gauge (CH0,CH1,CH2 and CH3) and one for tensor strain.

File Name Convention:

	BNUM.CH0.event.txt.gz
	BNUM.CH1.event.txt.gz
	BNUM.CH2.event.txt.gz
	BNUM.CH3.event.txt.gz
	BNUM.tensor.event.txt.gz

	BNUM refers to the 4-character ID, e.g., B004 and "event" refers to the event name.
	E.g., B004.CH0.20110311Tohoku.txt.gz

File Format

	The format for the CH0, CH1, CH2 and CH3 files is:
	Column		Decscription		Units
	Column 1	DateTime		UTC            
	Column 2	RawStrain		digital counts
	Column 3	LinearStrain		microstrain         
	Column 4	Tide 			microstrain         
	Column 5	BarometricCorrection	microstrain         
	Column 6	BarometricPressure	millibars           
	Column 7	PorePressure		millibars       
	Column 8	Version			generation date

	The format for the tensor strain file is:
	Column		Decscription		Units
	Column 1	DateTime		UTC               
	Column 2	Eee+Enn			microstrain      
	Column 3	Eee+Enn_tide		microstrain   
	Column 4	Eee+Enn_baro		microstrain    
	Column 5	Eee-Enn			microstrain     
	Column 6	Eee-Enn_tide		microstrain    
	Column 7	Eee-Enn_baro		microstrain   
	Column 8	2Ene			microstrain      
	Column 9	2Ene_tide		microstrain    
	Column 10	2Ene_baro		microstrain     
	Column 11	BaroPressure		millibar
	Column 12	PorePressure		millibar 
	Column 13	Version      		generation date 
	Column 14	TiltX			microradian    
	Column 15	TiltY			microradian 

	The null value is 999999. If 1-sps barometric data are available it will be given otherwise there
	will only be 30 minute interval data. Pore Pressure and tiltmeter data are only available for a
	subset of sites. See http://pbo.unavco.org/network/soh_map for instrument locations.

	Barometric, pore and tiltmeter time-stamps have been rounded to the nearest second.

	For more information on data products contact hodgkinson@unavco.org or borsa@unavco.org .

	

	Thu Sep 29 2011	
