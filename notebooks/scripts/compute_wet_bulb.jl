
using NCDatasets
using Psychro

path="/home/datawork-lops-drakkarcom/SIMULATION-OUTPUTS/ICE-CHANEL/viz/ACCESS-CM2/monthly/"

files = readdir(path)
#print(files)

temp_files = [item for item in files if occursin("tas_day_ACCESS", item) ] 

Threads.@threads for tempfile in temp_files
  println("file = $tempfile on thread $(Threads.threadid())")

  #temp = NCDataset(path*"monthly/tas_day_ACCESS-CM_1950.nc", "r")
  temp = NCDataset(path*tempfile, "r")
  husfile = replace(tempfile, "tas_day_ACCESS" => "huss_day_ACCESS", count=1)
  #shum = NCDataset(path*"monthly/huss_day_ACCESS-CM_1950.nc", "r")
  shum = NCDataset(path*husfile, "r")
  P = 101325.0

  temperature =  temp["tas"][:,:,:]
  spec_humid  =  shum["huss"][:,:,:]

  # Replace values to mininmum stable value
  temperature = nomissing(temperature,273.15)
  spec_humid = nomissing(spec_humid,0.0001)

  # Remove problematic values with temperatures smaller than 273 K
  temperature[temperature.<=273.15] .= 273.15

  # Compute dewpoint and wetbulb temperature 
  println("Computing")
  dewp = Psychro.dewpoint.(Psychro.MoistAir, temperature, Psychro.SpecHum, spec_humid, P)
  wetb = Psychro.wetbulb.(Psychro.MoistAir, temperature, Psychro.DewPoint, dewp, P)
  
  # Store file
  outfile = replace(tempfile, "tas_day_ACCESS" => "wetb_day_ACCESS", count=1)
  ds = NCDataset(path*outfile,"c")
  
  defDim(ds,"lon",1440)
  defDim(ds,"lat",600) 
  defDim(ds,"time",12) 
  
  ds.attrib["title"] = "Wet bulb temperature from ACCESS-CM"
  v = defVar(ds,"WetBTemp",Float32,("lon","lat","time"))

  v[:,:] = wetb

  v.attrib["units"] = "degree Celsius"
  v.attrib["comments"] = "The wet bulb temperature is computed using Psychro.jl"
  
  close(ds)
end
