# initialization file

# please modify this path if you want a specific path to open, e.g. 
# path="/home/federico/Documents/HighOrderTR/correctedTR"
path = pwd()

# creation of a new folder for every new day
using Dates
newday=Dates.format(Date(now()), "yyyy-mm-dd")
newyear=newday[1:4]
if (~isdir(path*"/"*newyear))
  mkdir(path*"/"*newyear)
end
newmonth=newday[1:7]
if (~isdir(path*"/"*newyear*"/"*newmonth))
  mkdir(path*"/"*newyear*"/"*newmonth)
end
newdir=path*"/"*newyear*"/"*newmonth*"/"*newday;
if (~isdir(newdir))
  mkdir(newdir)
end

function nrunf!(nrun) # function to create a subfolder for a specific run
  nrun[1]+=1
  global newdir_nrun = newdir*"/$(nrun[1])"
  if (~isdir(newdir_nrun))
    mkdir(newdir_nrun)
  end
end

# initialization of necessary packages. Use e.g. "] install DelimitedFiles" to install the package
using DelimitedFiles
using SpecialFunctions
using PyPlot
using QuadGK
using LinearAlgebra
using Statistics

nrun = [0] # run counter initialization