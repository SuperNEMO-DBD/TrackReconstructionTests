cp ../CAT/SD.brio ../TIT/SD.brio
flreconstruct -i ../TIT/SD.brio -p MockCalibratePipeline.conf -o ../TIT/CD.brio
flreconstruct -i ../TIT/CD.brio -p TKReconstructPipeline.conf -o ../TIT/TTD.brio
flreconstruct -i ../TIT/TTD.brio -p PTD_tracking.conf -o ../TIT/PTD.brio
