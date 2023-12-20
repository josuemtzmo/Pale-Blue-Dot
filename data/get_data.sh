
# Get data for all years:
# 1950 - 2100
#

model="ACCESS-CM2"
vars=('tasmax' 'tas' 'huss')

mkdir -p $model
cd $model

for year in $(seq 1950 2100 )
do
  echo $year
  for var in $(seq 0 $((${#vars[@]}-1)))
  do
    echo ${vars[$var]}
    if [ $year -lt 2015 ]
    then
      sim='historical'
    else
      sim='ssp585'
    fi
    #echo $sim
    path="https://nex-gddp-cmip6.s3-us-west-2.amazonaws.com/NEX-GDDP-CMIP6/$model/$sim/r1i1p1f1/${vars[$var]}/${vars[$var]}_day_ACCESS-CM2_${sim}_r1i1p1f1_gn_${year}.nc"


    wget $path
    #echo "wget 'https://ds.nccs.nasa.gov/thredds/ncss/grid/AMES/NEX/GDDP-CMIP6/$model/$sim/r1i1p1f1/${vars[$var]}/${vars[$var]}_day_ACCESS-CM2_historical_r1i1p1f1_gn_2014.nc?var=${vars[$var]}&time_start=2014-01-01T12:00:00Z&time_end=2014-12-31T12:00:00Z&&&accept=netcdf3&addLatLon=true'"

  done
done

cd ..
