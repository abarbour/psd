#!/bin/csh

	# Script to download 20sps data
	#
	# Usage: get_20sps.csh Bnum 16-character-code start_year  start_day_of_yr  start_time  end_year end_day_of_year end_time
	# e.g.,  get_20sps.sh B073 varian073bcs2006  2009 105 13:00:00 2009 105 16:00:00 p
	# the last argument should be used to make a GMT plot if GMT is installed

        # This should produce four files
        #       ${name}.ALL_20.r.txt  raw        data for 4 channels (nanostrain)
        #       ${name}.ALL_20.l.txt  linearized data for 4 channels (nanostrain)
        #       ${name}.ps             postscipt plot of the linearized data if "p" used in the command line
        #       B073.2009105_20.tar


        ########### set variables   ########################

        set name=$1
        set long=$2
        set syear=$3
        set sday=$4
        set start=$5
        set eyear=$6
        set eday=$7
        set end=$8
        set plot=$9

        set shour=`echo $start | awk '{ print substr($1,1,2) }' `
        set smin=`echo  $start | awk '{ print substr($1,4,2) }' `

        set ehour=`echo $end | awk '{ print substr($1,1,2) }' `
        set emin=`echo  $end | awk '{ print substr($1,4,2) }' `

        echo $name $syear $sday $shour $smin    $eyear  $eday $ehour $emin

        ########## get data ######################

        set year=$syear
        set day=$sday
        set limitday=366
        while ($year <= $eyear)
                if ($year == $eyear) set limitday=$eday
                while ($day <= $limitday )
                        set day=`echo $day  | awk '{ printf "%03d\n", $1 }' `
                        set sy=`echo $year  | awk '{ printf "%02d\n", $1-2000 }' `
                        ftp  ftp://www.ncedc.org/pub/pbo/strain/raw/bsm/$long/$year/$day/${name}.${year}${day}_20.tar
                        tar -xf ${name}.${year}${day}_20.tar
                        @ day ++
                end
                set day=1
                @ year ++
        end

        ## unpack day and hour tars , select files inside start and end times ##

        set year=$syear
        set day=$sday
        set limitday=366
        set hour=$shour
        set limith=23
        while ($year <= $eyear)
                if ($year == $eyear) set limitday=$eday
                while ($day <= $limitday )
                        set day=`echo $day  | awk '{ printf "%03d\n", $1      }' `
                        set sy=`echo $year  | awk '{ printf "%02d\n", $1-2000 }' `
                        if ($day == $eday) set limith=$ehour
                        while ($hour <= $limith )
                                set hour=`echo $hour | awk '{ printf"%02d\n",$1}' `
                                tar -xf ${name}${sy}${day}${hour}_20.tar
                                @ hour ++
                        end
                        set hour=0;
                        @ day ++
                end
                set day=1
                @ year ++
        end
        rm ${name}???????_20.tar

	## unpack minute files ##

 	foreach file (`ls ${name}?????????_20.tgz`)
 		tar xfz $file
		rm $file
 	end


	####  append channel files #####
	
 	foreach c (CH0 CH1 CH2 CH3)
                set year=$syear
                set day=$sday
                set limitday=366
                set hour=$shour
                set limithour=23
		set min=$smin
		set limitmin=59
                set flag=0
		while ($year <= $eyear)
			if ($year == $eyear) set limitday=$eday
			while ($day <= $eday )
				set day=`echo $day  | awk '{ printf "%03d\n", $1      }' `
				set sy=`echo $year  | awk '{ printf "%02d\n", $1-2000 }' `
				if ($day == $eday) set limithour=$ehour
				while ($hour <= $limithour)
 	               			set hour=`echo $hour | awk '{ printf"%02d\n",$1   }' `
 					if ($hour == $ehour)  set limitmin=$emin
 					while ($min < $limitmin )
 						set min=`echo $min | awk '{ printf "%02d\n", $1 }'`
 						if ( $flag == 0 ) then
 							cp {$name}${sy}${day}${hour}${min}${c}_20 {$name}.${c}_20
 							set flag=1
 						else
 							bottle_merge.py {$name}.${c}_20 {$name}${sy}${day}${hour}$min${c}_20 {$name}.${c}_20
 						endif
 						@ min++	
 					end
 					@ hour ++
 					set min=0;
 				end
				@ day ++
				set hour=0
			end
			@ year ++
			set day=1
		end
	end
 
 	rm ${name}?????????CH?_20

	####### make ascii linearized data ##########

	# Decide instrument gap
	# first batch have gap = 0.002 others have gap = 0.001, [B001][B003][B004][B005][B006][B007][B009][B010][B011][B012][B018][B022][B024][B035][B081][B082][B086][B087]:

	switch ($name)
		case [B][0][0]?
			set gap=0.0002
			breaksw
		case [B][0][1][0128]
			set gap=0.0002
			breaksw
		case [B][0][2][24]
                        set gap=0.0002
                        breaksw
		case [B][0][8][1267]
                        set gap=0.0002
                        breaksw
		default:
			set gap=0.0001
			breaksw
	endsw


       foreach c (CH0 CH1 CH2 CH3)
                bottle.py -t ${name}.${c}_20 > ${name}.${c}_20.r.txt
		cat ${name}.${c}_20.r.txt | awk '{  di=($3/1e+8)/(1-($3/1e+8)) ; if (f == 0 && $3 != 999999)   {d0=di;f=1}  ; if (/999999/) printf ("%s %s	999999\n", $1, $2, $3) ; else printf("%s %s	%14.9f \n", $1,$2, ((di-d0)*'$gap'/0.087)*1e+9)  }' > ${name}.${c}_20.l.txt
		rm ${name}.${c}_20
        end

	echo "Date-Time	CH0	CH1	CH2	CH3" > ${name}.ALL_20.l.txt
	echo "Date-Time	CH0	CH1	CH2	CH3" > ${name}.ALL_20.r.txt
	paste ${name}.CH0_20.l.txt ${name}.CH1_20.l.txt ${name}.CH2_20.l.txt ${name}.CH3_20.l.txt | awk '{ print $1"T"$2"\t"$3"\t"$6"\t"$9"\t"$12 }' >>  ${name}.ALL_20.l.txt
	paste ${name}.CH0_20.r.txt ${name}.CH1_20.r.txt ${name}.CH2_20.r.txt ${name}.CH3_20.r.txt | awk '{ print $1"T"$2"\t"$3"\t"$6"\t"$9"\t"$12 }' >>  ${name}.ALL_20.r.txt

	rm ${name}.CH?_20.l.txt ${name}.CH?_20.r.txt	

        #### plot linearized data ####

	if ($plot == p) then
	        gmtset PLOT_DATE_FORMAT "o dd yyyy"
	        gmtset PLOT_CLOCK_FORMAT hh:mm
	        gmtset PAPER_MEDIA      letter
	        gmtset TIME_INTERVAL_FRACTION 0.01
	        gmtset LABEL_FONT_SIZE  15
	        set xt_primary="pa1Hf5m"
	        set xt_secondary="sa1DS"
	
	        set plotx=`grep -v 999999 ${name}.ALL_20.l.txt | minmax -fT -C  `
	        set ploty=`grep -v 999999 ${name}.ALL_20.l.txt | minmax -fT -C | awk ' { print $3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8"\n"$9"\n"$10"\n"$11"\n"$12_ }' | minmax -C -I1`
	
	        set yspace=`echo $ploty | awk '{ printf"%3d\n", ($2-$1)/4 }'`
	
	        psbasemap  -R$plotx[1]/$plotx[2]/$ploty[1]/$ploty[2] -JX7iT/5i -B${xt_primary}/${yspace}:."${name}":  -B${xt_secondary}  -K -Y5.0  -P -X2.5 > ${name}.ps
	        awk 'NR>1 { print $1, $2 }'  ${name}.ALL_20.l.txt | psxy -J -R -O  -B${xt_primary}/${yspace}:"nanostrain":WSne   -Sc0.1 -Gblue  -K >> ${name}.ps
	        awk 'NR>1 { print $1, $3 }'  ${name}.ALL_20.l.txt | psxy -J -R -O    -Sc0.1 -Gred   -K >> ${name}.ps
	        awk 'NR>1 { print $1, $4 }'  ${name}.ALL_20.l.txt | psxy -J -R -O    -Sc0.1 -Ggreen -K >> ${name}.ps
	        awk 'NR>1 { print $1, $5 }'  ${name}.ALL_20.l.txt | psxy -J -R -O    -Sc0.1 -Gcyan  -K >> ${name}.ps
	        echo "$plotx[1] $ploty[1] 20 0 1 BL CH0" | pstext -R -J -N -O -Gblue  -D0.50i/0.25i -K >> ${name}.ps
	        echo "$plotx[1] $ploty[1] 20 0 1 BL CH1" | pstext -R -J -N -O -Gred   -D1.25i/0.25i -K >> ${name}.ps
	        echo "$plotx[1] $ploty[1] 20 0 1 BL CH2" | pstext -R -J -N -O -Ggreen -D2.00i/0.25i -K >> ${name}.ps
	        echo "$plotx[1] $ploty[1] 20 0 1 BL CH3" | pstext -R -J -N -O -Gcyan  -D2.75i/0.25i    >> ${name}.ps
	endif	
