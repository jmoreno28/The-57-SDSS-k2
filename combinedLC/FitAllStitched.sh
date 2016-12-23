#edit base dir-----path to folder with lc files

declare -a sdssID
let i=0
while IFS=$'\n' read -r line_data; do
    # Parse “${line_data}” to produce content 
    sdssID[i]="${line_data}" # Populate array.
    ((++i))
done < sdssList.dat
echo i


let x=0
while IFS=$'\n' read -r object; do
python combinedLC_Fit.py -id ${object} -sdssid ${sdssID[x]} -p 'k2sff' -c c08 -pMax 4
#if folder does not exist, make folder
newdir=("SDSS-K2"${object}"-CARMA4")
mkdir ${newdir}
mv ${object}* ./${newdir}/
echo "ran routine"
((++x))
echo x ${object} ${sdssID[x]}
done < k2List.dat
