#!/bin/sh
#cat occ_list | while read line 
#do
#   echo $line
#done

exec 3<occName_FirstL
exec 4<occlist_FirstL
asset="users/hywong/occ_gcloud/"
#gcloud_pre="gs://nymphs_bucket_occ/"

while read occName <&3 && read occData <&4
do
        earthengine upload table --asset_id $asset$occName --x_column decimalLongitude --y_column decimalLatitude --crs EPSG:4326 $occData
	echo $asset$occName
done
exec 3<&-
exec 4<&-
