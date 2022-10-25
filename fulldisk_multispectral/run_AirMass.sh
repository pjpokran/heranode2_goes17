#!/bin/bash

export PATH="/home/poker/miniconda3/envs/goes16_201710/bin/:$PATH"

export LD_LIBRARY_PATH=


cd /home/poker/goes17/fulldisk_multispectral

cp /home/ldm/data/grb-west/fulldisk/12/latest.nc /dev/shm/latest_17_fulldisk_airmass_12.nc
cmp /home/ldm/data/grb-west/fulldisk/12/latest.nc /dev/shm/latest_17_fulldisk_airmass_12.nc > /dev/null
CONDITION=$?
#echo $CONDITION

while :; do

  until [ $CONDITION -eq 1 ] ; do
#     echo same
     sleep 5
     cmp /home/ldm/data/grb-west/fulldisk/12/latest.nc /dev/shm/latest_17_fulldisk_airmass_12.nc > /dev/null
     CONDITION=$?
  done
#  echo different
  sleep 195
  cp /home/ldm/data/grb-west/fulldisk/12/latest.nc /dev/shm/latest_17_fulldisk_airmass_12.nc
  /home/poker/miniconda3/envs/goes16_201710/bin/python goes17_fulldisk_GRB_AirMass.py /home/ldm/data/grb-west/fulldisk/08/latest.nc  /home/ldm/data/grb-west/fulldisk/10/latest.nc  /home/ldm/data/grb-west/fulldisk/12/latest.nc /home/ldm/data/grb-west/fulldisk/13/latest.nc
  cmp /home/ldm/data/grb-west/fulldisk/12/latest.nc /dev/shm/latest_17_fulldisk_airmass_12.nc > /dev/null
  CONDITION=$?
#  echo repeat

done

