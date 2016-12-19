#edit base dir-----path to folder with lc files
basedir=("Research/code/The-57-SDSS-k2/")



let x=1
while IFS=$'\n' read -r object; do
python sdssFitRecalib.py -sdssid ${object} -pMax 4
#if folder does not exist, make folder
newdir=(${object}"-CARMA4")
mkdir ${newdir}
mv ${object}* ./${newdir}/
echo "ran routine"
((++x))
echo x ${object}
done < sdssList.dat