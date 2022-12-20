#!/bin/sh

target_files="gs://nymphs_bucket_occ/occdata_rmFull_Global_raw_reduced/*.csv"

gsutil ls $target_files | sed 's:.*/::' | cut -d '.' -f 1 > ./occNamedelete
gsutil ls $target_files > ./occdelete

#echo ./occdelete

exec 3<./occNamedelete
exec 4<./occdelete

# change target GEE asset folder to ingest:
asset="users/hywong/occ_gcloud_global/"

while read occName <&3 && read occData <&4
do
        earthengine upload table --asset_id $asset$occName --x_column decimalLongitude --y_column decimalLatitude --crs EPSG:4326 $occData
	echo $asset$occName
done
exec 3<&-
exec 4<&-
