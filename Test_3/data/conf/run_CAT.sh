rm ../CAT/SD.brio
rm ../CAT/PTD.brio

flsimulate -c simu_Bi.conf -o ../CAT/SD.brio
flreconstruct -i ../CAT/SD.brio -p reco_CAT.conf -o ../CAT/PTD.brio
