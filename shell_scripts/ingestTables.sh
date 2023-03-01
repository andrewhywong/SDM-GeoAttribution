#!/bin/sh

# change source gcp occ location:
target_files="gs://nymphs_bucket_occ/occdata_rmFull_Global_raw_reduced/*.csv"

# change target GEE asset folder to ingest:
asset="users/hywong/target_folder/"

gsutil ls $target_files | sed 's:.*/::' | cut -d '.' -f 1 > ./occNames
gsutil ls $target_files > ./occLoc

#echo ./occdelete

exec 3<./occNames
exec 4<./occLoc

while read occName <&3 && read occData <&4
do
        #earthengine upload table --asset_id $asset$occName --x_column decimalLongitude --y_column decimalLatitude --crs EPSG:4326 $occData
	echo $asset$occName
done

exec 3<&-
exec 4<&-
