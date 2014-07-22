g++ `root-config --cflags` runQRv2.cc -o runQRv2 `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` `root-config --libs`

file=jetfile_maindata_p12_glauberfile_tree3_embedded_Phob_Glau_PbPb_cent_0_10_nevnt100000_quench_method_3_factor90_scale100
./runQRv2 datafile_$file.txt QR_ntuple$file.root


# directory="/afs/cern.ch/user/d/dgulhan/workDir/ntuples_Rad_450_0/cent_0_10"
 # for file in `ls /afs/cern.ch/user/d/dgulhan/workDir/ntuples_Rad_450_0/cent_0_10/datafile* | sed "s/datafile/ /g" | awk '{print $2}'`
  # do
    # a=$(ls -l /afs/cern.ch/user/d/dgulhan/workDir/ntuples_Rad_450_0/cent_0_10/datafile${file} | awk '{print $5}')
    # if [ "$a" -gt "0" ]
         # then
   # ./runQR $directory/datafile$file $directory/QR_ntuple$file
    # fi
  # done
 