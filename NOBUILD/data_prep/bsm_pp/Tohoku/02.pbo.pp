#!/bin/sh
echo "non necessary: any available pp data is contained in the _B directory"
exit
PP="http://pore.unavco.org/pore"
PPFI="PorePressHPa.txt.gz"
TFI="PoreTempDegC.txt.gz"
cat << XXX
wget -O B081_P.txt.gz ${PP}/B081/2011/070/B08111070${PPFI}
wget -O B081_T.txt.gz ${PP}/B081/2011/070/B08111070${TFI}
wget -O B082_P.txt.gz ${PP}/B082/2011/070/B08211070${PPFI}
wget -O B082_T.txt.gz ${PP}/B082/2011/070/B08211070${TFI}
wget -O B084_P.txt.gz ${PP}/B084/2011/070/B08411070${PPFI}
wget -O B084_T.txt.gz ${PP}/B084/2011/070/B08411070${TFI}
wget -O B086_P.txt.gz ${PP}/B086/2011/070/B08611070${PPFI}
wget -O B086_T.txt.gz ${PP}/B086/2011/070/B08611070${TFI}
wget -O B087_P.txt.gz ${PP}/B087/2011/070/B08711070${PPFI}
wget -O B087_P.txt.gz ${PP}/B087/2011/070/B08711070${PPFI}
wget -O B088_T.txt.gz ${PP}/B088/2011/070/B08811070${TFI}
wget -O B088_T.txt.gz ${PP}/B088/2011/070/B08811070${TFI}
wget -O B946_P.txt.gz ${PP}/B946/2011/070/B94611070${PPFI}
wget -O B946_T.txt.gz ${PP}/B946/2011/070/B94611070${TFI}
XXX
