#  Copyright 2016, Mehdi Madani, Mathieu Van Vyve
#  This Source Code Form is subject to the terms of the GNU GENERAL PUBLIC LICENSE version 3
#  If a copy of the GPL version 3 was not distributed with this
#  file, You can obtain one at http://www.gnu.org/licenses/gpl-3.0.en.html
#############################################################################


type damdata
  areas::DataFrame
  periods::DataFrame
  hourly::DataFrame
  mp_headers::DataFrame
  mp_hourly::DataFrame
  line_capacities::DataFrame
end

type damsol
  xval #hourly bids
  uval #on/off decisions for MP bids
  xhval #hourly bids associated to MP bids
  fval  #flows through lines
  sval  #surplus of hourly bids
  scval #surplus of MP bids
  shmaxval
  shminval
  vval
  priceval
end

function incomeCheck(mydata::damdata, mysol::damsol)
  areas = mydata.areas
  periods = mydata.periods
  hourly = mydata.hourly
  mp_headers = mydata.mp_headers
  mp_hourly = mydata.mp_hourly
  line_capacities = mydata.line_capacities

  areas=Array(areas)
  periods=Array(periods)

  nbHourly=nrow(hourly)
  nbMp=nrow(mp_headers)
  nbMpHourly=nrow(mp_hourly)
  nbAreas=length(areas)
  nbPeriods=length(periods)

  xval = convert(Vector{Float64}, mysol.xval)
  uval = convert(Vector{Float64}, mysol.uval)
  xhval = convert(Vector{Float64}, mysol.xhval)
  sval = convert(Vector{Float64}, mysol.sval)
  scval = convert(Vector{Float64}, mysol.scval)
  shmaxval = convert(Vector{Float64}, mysol.shmaxval)
  shminval = convert(Vector{Float64}, mysol.shminval)
  #priceval = convert( Vector{Float64} , mysol.priceval)
  priceI = [mysol.priceval[hourly[i, :LI], hourly[i, :TI]] for i in 1:nbHourly]
  priceI = convert(Vector{Float64}, priceI)
  priceH = [mysol.priceval[mp_hourly[h, :LH], mp_hourly[h, :TH]] for h in 1:nbMpHourly]
  priceH = convert(Vector{Float64}, priceH)

  income_real = zeros(nbMp)
  income_linear = scval + mp_headers[:,:FC].*uval
  fixedcosts = mp_headers[:,:FC].*uval
  marginalcosts = zeros(nbMp)
  sc_values = [scval[c] for c in 1:nbMp]
  shmax_values = zeros(nbMp)
  shmin_AR_values = zeros(nbMp)
for c in 1:nbMp
  for h in 1:nbMpHourly
    if(mp_hourly[h,:MP] == mp_headers[c, :MP])
      income_real[c] -= xhval[h]*mp_hourly[h,:QH]*priceH[h]
      income_linear[c] -= xhval[h]*mp_hourly[h,:QH]*mp_hourly[h,:PH]
      marginalcosts[c] -= xhval[h]*mp_hourly[h,:QH]*mp_hourly[h,:PH]
      shmax_values[c]  += shmaxval[h]*uval[c]
      shmin_AR_values[c]  += shminval[h]*mp_hourly[h,:AR]*uval[c]
    end
  end
end
  PL = income_real - fixedcosts - marginalcosts
  incomes = DataFrame([income_real income_linear fixedcosts marginalcosts sc_values PL shmax_values shmin_AR_values]);
  rename!(incomes, :x1, :income_real)
  rename!(incomes, :x2, :income_linear)
  rename!(incomes, :x3, :fixedcosts)
  rename!(incomes, :x4, :marginalcosts)
  rename!(incomes, :x5, :sc_values)
  rename!(incomes, :x6, :ProfitsLosses)
  rename!(incomes, :x7, :sum_shmax)
  rename!(incomes, :x8, :sum_shmin_AR)
  return incomes
end

function solQstats(mydata::damdata, mysol::damsol)
  areas = mydata.areas
  periods = mydata.periods
  hourly = mydata.hourly
  mp_headers = mydata.mp_headers
  mp_hourly = mydata.mp_hourly
  line_capacities = mydata.line_capacities

  areas=Array(areas)
  periods=Array(periods)

  nbHourly=nrow(hourly)
  nbMp=nrow(mp_headers)
  nbMpHourly=nrow(mp_hourly)
  nbAreas=length(areas)
  nbPeriods=length(periods)


  xval = convert(Vector{Float64}, mysol.xval)
  uval = convert(Vector{Float64}, mysol.uval)
  xhval = convert(Vector{Float64}, mysol.xhval)
  sval = convert(Vector{Float64}, mysol.sval)
  scval = convert(Vector{Float64}, mysol.scval)
  shmaxval = convert(Vector{Float64}, mysol.shmaxval)
  shminval = convert(Vector{Float64}, mysol.shminval)
  #priceval = convert( Vector{Float64} , mysol.priceval)
  priceI = [mysol.priceval[hourly[i, :LI], hourly[i, :TI]] for i in 1:nbHourly]
  priceI = convert(Vector{Float64}, priceI)
  priceH = [mysol.priceval[mp_hourly[h, :LH], mp_hourly[h, :TH]] for h in 1:nbMpHourly]
  priceH = convert(Vector{Float64}, priceH)


fval_loc = Array{Float64}(nrow(line_capacities))
vval_loc = Array{Float64}(nrow(line_capacities))
Networkslack_d = Array{Float64}(nrow(line_capacities))
for i in 1:nrow(line_capacities)
      fval_loc[i] = mysol.fval[line_capacities[i,:from], line_capacities[i,:too] , line_capacities[i,:t] ]
      vval_loc[i] = mysol.vval[line_capacities[i,:from], line_capacities[i,:too] , line_capacities[i,:t] ]
      Networkslack_d[i] = vval_loc[i] - mysol.priceval[line_capacities[i,:too], line_capacities[i,:t]] +  mysol.priceval[line_capacities[i,:from], line_capacities[i,:t]]
end


  ### PRIMAL SLACKS TESTING
  Islack_p = (1-xval)
  maxminslack_hourly_p=maximum(min(sval, Islack_p))

  uval_match_xh = [ uval[find( mp_headers[:,:MP].== mp_hourly[h,:MP] ) ][1]  for h in 1:nbMpHourly]

  Hslack_max_p = (uval_match_xh - xhval)
  Hslack_max_p = convert(Vector{Float64}, Hslack_max_p)
  maxminslack_xh_max = maximum(min(shmaxval, Hslack_max_p))

  Hslack_min_p = (xhval - uval_match_xh.*mp_hourly[:,:AR])
  Hslack_min_p = convert(Vector{Float64}, Hslack_min_p)
  maxminslack_xh_min = maximum(min(shminval, Hslack_min_p))


  Cslack_p = (1-uval)
  maxminslack_mp_p=maximum(min(scval, Cslack_p))

  Networkslack_p = (line_capacities[:,:linecap] - fval_loc)
  maxminslack_network_p = maximum(min(vval_loc, Networkslack_p))

  ### DUAL SLACKS TESTING
  Islack_d = sval - hourly[:,:QI].*(hourly[:,:PI0] - priceI)
  maxminslack_hourly_d=maximum(min(xval, Islack_d))

  surp_hmaxmin = zeros(nbMp)
  for h in 1:nbMpHourly
    for c in 1:nbMp
      if(mp_hourly[h,:MP] == mp_headers[c, :MP])
        surp_hmaxmin[c] += shmaxval[h] - mp_hourly[h,:AR]*shminval[h]
      end
    end
  end

  Cslack_d = scval  - surp_hmaxmin + mp_headers[:,:FC]
  Cslack_d = convert(Vector{Float64}, Cslack_d)
  maxminslack_mp_d=maximum(min(uval, Cslack_d))

  maxminslack_network_d = maximum(min(fval_loc, Networkslack_d))

############################ Quality summary
maxSlackViolation = max(maxminslack_hourly_p , maxminslack_hourly_d, maxminslack_xh_max, maxminslack_xh_min, maxminslack_mp_p, maxminslack_mp_d, maxminslack_network_p, maxminslack_network_d)
solstats = DataFrame([maxSlackViolation maxminslack_hourly_p  maxminslack_hourly_d maxminslack_xh_max maxminslack_xh_min maxminslack_mp_p maxminslack_mp_d maxminslack_network_p maxminslack_network_d])
rename!(solstats, :x1, :maxSlackViolation)
rename!(solstats, :x2, :maxminslack_hourly_p)
rename!(solstats, :x3, :maxminslack_hourly_d)
rename!(solstats, :x4, :maxminslack_xh_max)
rename!(solstats, :x5, :maxminslack_xh_min)
rename!(solstats, :x6, :maxminslack_mp_p)
rename!(solstats, :x7, :maxminslack_mp_d)
rename!(solstats, :x8, :maxminslack_network_p)
rename!(solstats, :x9, :maxminslack_network_d)
return solstats
end
