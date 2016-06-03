#  Copyright 2016, Mehdi Madani
#  This Source Code Form is subject to the terms of the GNU GPL version 3
#  If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://www.gnu.org/licenses/gpl-3.0.en.html
#############################################################################

#pwd()
#cd("H:/revisiting_mp_conditions/")

using JuMP,  DataFrames, DataArrays, CPLEX #, Gurobi, Cbc
include("dam_utils.jl")

#method_type = "benders_classic"
method_type = "benders_modern"
#method_type = "primal-dual"
#method_type = "mp_relaxed"
mic_activation=0

if method_type != "primal-dual" &&  mic_activation != 0
  error("ERROR: OMIE complex orders with a minimum income condition only handled with the primal-dual method type")
end

pricecap_up = 3000 # market price range restrictions, upper bound,  in accordance to current European market rules (PCR)
pricecap_down = -500 # market price range restrictions, lower bound, in accordance to current European market rules (PCR)

numTests = DataFrame(Inst = Int64[], welfare = Float64[], absgap=Float64[], solQ = Float64[], lazycuts = Int64[], solvercuts = Int64[], nodes = Int64[], runtime=Float64[], Nb_MpBids=Float64[], Nb_Steps=Float64[])
solQ = DataFrame(Inst = Int64[], maxSlackViolation = Float64[], maxminslack_hourly_p = Float64[], maxminslack_hourly_d=Float64[], maxminslack_xh_max = Float64[], maxminslack_xh_min = Float64[], maxminslack_mp_p = Float64[], maxminslack_mp_d=Float64[], maxminslack_network_p=Float64[], maxminslack_network_d=Float64[])


for(ssid in 1:10)

solutions = Array{damsol}(1)
solutions_u = Array{Vector{Float64}}(1)

areas=readtable(string("./data/daminst-", ssid,"/areas.csv"))
periods=readtable(string("./data/daminst-", ssid,"/periods.csv"))
hourly=readtable(string("./data/daminst-", ssid,"/hourly_quad.csv"))
mp_headers=readtable(string("./data/daminst-", ssid,"/mp_headers.csv"))
mp_hourly=readtable(string("./data/daminst-", ssid,"/mp_hourly.csv"))
line_capacities=readtable(string("./data/daminst-", ssid,"/line_cap.csv"))

if mic_activation == 1
  FC_copy = mp_headers[:,:FC]
  mp_headers[:,:FC] = 0
end

mydata = damdata(areas, periods, hourly, mp_headers, mp_hourly, line_capacities)

areas=Array(areas)
periods=Array(periods)

nbHourly=nrow(hourly)
nbMp=nrow(mp_headers)
nbMpHourly=nrow(mp_hourly)
nbAreas=length(areas)
nbPeriods=length(periods)

# big-M's computation
bigM = zeros(Float64, nbMp)

for c in 1:nbMp
  bigM[c] -= mp_headers[c, :FC]
end

for h in 1:nbMpHourly
  bigM[find( mp_headers[:,:MP].== mp_hourly[h,:MP] )] += (pricecap_up - mp_hourly[h,:PH])*abs(mp_hourly[h,:QH]) #
end

m = Model(solver=CplexSolver(CPX_PARAM_EPGAP=1e-8, CPX_PARAM_EPAGAP=1e-8, CPX_PARAM_EPINT=1e-7, CPX_PARAM_TILIM=600, CPX_PARAM_BRDIR=-1, CPX_PARAM_HEURFREQ=-1))
#, CPX_PARAM_EACHCUTLIM=0, CPX_PARAM_FRACCUTS=-1, CPX_PARAM_EACHCUTLIM=0, CPX_PARAM_LPMETHOD=2 , CPX_PARAM_THREADS=12 #CPX_PARAM_MIPDISPLAY=1, CPX_PARAM_MIPINTERVAL=1
#m = Model(solver = GurobiSolver(BranchDir=-1, MIPGap=1e-9, MIPGapAbs=1e-9, IntFeasTol=1e-9, TimeLimit=600)) #Method=1,
#m = Model(solver = CbcSolver(integerTolerance=1e-9, ratioGap=1e-9, allowableGap=1e-9 ))

@variable(m, 0<= x[1:nbHourly] <=1)
@variable(m, u[1:nbMp], Bin)
@variable(m, xh[1:nbMpHourly])
@variable(m, 0<= f[areas, areas, periods] <= 0) # by default/at declaration, a line doesn't exist, though non-zero upper bounds on flows are set later for existing lines with non zero capacity

for i in 1:nrow(line_capacities)
    setupperbound(f[line_capacities[i,:from], line_capacities[i,:too] , line_capacities[i,:t] ], line_capacities[i,:linecap] ) #
end

@constraint(m, mp_control_upperbound[h=1:nbMpHourly, c=1:nbMp;  mp_hourly[h,:MP] == mp_headers[c,:MP]],
                xh[h] <= u[c]) #dual variable is shmax[h]
@constraint(m, mp_control_lowerbound[h=1:nbMpHourly, c=1:nbMp;  mp_hourly[h,:MP] == mp_headers[c,:MP]],
                xh[h] >= mp_hourly[h,:AR]*u[c]) #dual variable is shmin[h]


#Objective: maximize welfare
obj = dot(x,(hourly[:,:QI].data).*(hourly[:,:PI0].data)) +  dot(xh,(mp_hourly[:,:QH].data).*(mp_hourly[:,:PH].data)) - dot(u, mp_headers[:, :FC])# #this adds the sum over hourly bids of bid-price*bid-quantity*x
@objective(m, Max,  obj)

@constraint(m, balance[loc in areas, t in periods],
sum{x[i]*hourly[i, :QI], i=1:nbHourly; hourly[i, :LI] == loc && hourly[i, :TI]==t }
+
sum{xh[h]*mp_hourly[h, :QH], h=1:nbMpHourly; mp_hourly[h,:LH] == loc && mp_hourly[h, :TH]==t } #doensét change any thing here: mphourly_loc replacing ... mp_headers[mp_hourly[h,:MP], :LC]
==
sum{f[loc_orig, loc, t] , loc_orig in areas; loc_orig != loc} - sum{f[loc, loc_dest,t] , loc_dest in areas; loc_dest != loc}
)


if method_type != "primal-dual"
mdual = Model(solver=CplexSolver(CPX_PARAM_LPMETHOD=2)) #CPX_PARAM_PREIND=0,
#mdual = Model(solver = GurobiSolver(BranchDir=-1, MIPGap=1e-9, MIPGapAbs=1e-9, IntFeasTol=1e-9, TimeLimit=600)) #Method=1,
#mdual = Model(solver = CbcSolver(integerTolerance=1e-9, ratioGap=1e-9, allowableGap=1e-9 ))
else mdual = m
end

@variable(mdual, s[1:nbHourly] >=0)
@variable(mdual, sc[1:nbMp] >=0)
@variable(mdual, shmax[1:nbMpHourly] >=0)
@variable(mdual, shmin[1:nbMpHourly] >=0)
@variable(mdual, price[areas, periods])
@variable(mdual, v[areas, areas, periods] >=0)

@constraint(mdual, hourlysurplus[i=1:nbHourly], s[i] + hourly[i,:QI]*price[hourly[i,:LI], hourly[i,:TI]] >= hourly[i,:QI]*hourly[i,:PI0])
@constraint(mdual, mphourlysurplus[h=1:nbMpHourly], (shmax[h] - shmin[h]) +  mp_hourly[h,:QH]*price[mp_hourly[h,:LH], mp_hourly[h,:TH]] == mp_hourly[h,:QH]*mp_hourly[h,:PH])
if method_type != "primal-dual"
  @constraint(mdual, mpsurplus[c=1:nbMp], sc[c] - sum{(shmax[h] - mp_hourly[h,:AR]*shmin[h]), h=1:nbMpHourly; mp_hourly[h,:MP] == mp_headers[c, :MP]} + mp_headers[c,:FC] >= 0)
elseif method_type == "primal-dual"
  @constraint(mdual, mpsurplus[c=1:nbMp], sc[c] - sum{(shmax[h] - mp_hourly[h,:AR]*shmin[h]), h=1:nbMpHourly; mp_hourly[h,:MP] == mp_headers[c, :MP]} + mp_headers[c,:FC] + bigM[c]*(1-u[c]) >= 0) #
end

@constraint(mdual, flowdual[loc1 in areas, loc2 in areas, t in periods], price[loc2, t] - price[loc1, t] <= v[loc1, loc2, t])

if method_type != "primal-dual"
  @objective(mdual, Min, sum{s[i], i=1:nbHourly} +  sum{sc[c], c=1:nbMp} + sum{v[loc1, loc2, t]*getupperbound(f[loc1, loc2, t]), loc1 in areas, loc2 in areas, t in periods})
elseif method_type == "primal-dual"
  @constraint(m, obj >= sum{s[i], i=1:nbHourly} +  sum{sc[c], c=1:nbMp} + sum{v[loc1, loc2, t]*getupperbound(f[loc1, loc2, t]), loc1 in areas, loc2 in areas, t in periods} )
  if mic_activation == 1
    @constraint(mdual, MIC[c=1:nbMp], sc[c] - sum{(mp_hourly[h,:QH]*mp_hourly[h,:PH]*xh[h]), h=1:nbMpHourly; mp_hourly[h,:MP] == mp_headers[c, :MP]} - FC_copy[c] + sum{(mp_hourly[h,:QH]*xh[h]*mp_headers[c,:VC]), h=1:nbMpHourly; mp_hourly[h,:MP] == mp_headers[c, :MP]} + FC_copy[c]*(1-u[c]) >= 0) #
  end
end

cut_count = [0]
tol=1e-5

if method_type == "benders_modern"
  function workercb(cb)
    uval=getvalue(u)
    for c in 1:nbMp
      JuMP.setRHS(mpsurplus[c], - bigM[c]*(1-uval[c]) - mp_headers[c,:FC] ) # used to activate / deactivate constraints, big-M's pose no problem here but could be avoided if needed
    end
    statusdual = solve(mdual)
    dualobjval = getobjectivevalue(mdual)
    objval=cbgetnodeobjval(cb)
    if dualobjval > objval + tol
      println("This primal solution is not valid, and is cut off ...")
      @lazyconstraint(cb, sum{(1-u[c]), c=1:nbMp; uval[c] == 1} + sum{u[c], c=1:nbMp; uval[c] == 0}  >= 1)   # globally valid 'no-good cut'
      @lazyconstraint(cb, sum{(1-u[c]), c=1:nbMp; uval[c] == 1}  >= 1, localcut=true)   # locally valid strengthened cut
      cut_count[1] += 1
      println("Cuts added: ", cut_count[1])
      println("node:", getnodecount(m))
    elseif dualobjval <= objval + tol
      fill!(solutions_u, uval)
    end
  end
addlazycallback(m, workercb)
end

status = solve(m)
objval=getobjectivevalue(m)

uval_c = getvalue(u)
if method_type == "benders_modern"
  for c in 1:nbMp
    JuMP.setRHS(mpsurplus[c], - bigM[c]*(1-uval_c[c]) - mp_headers[c,:FC] ) # used to activate / deactivate constraints, big-M's pose no problem here but could be avoided if needed
  end
  statusdual = solve(mdual)
end



if method_type == "benders_classic"
    uval=getvalue(u)
    tol=1e-4
    for c=1:nbMp
      JuMP.setRHS(mpsurplus[c], - bigM[c]*(1-uval[c]) - mp_headers[c,:FC] )
    end
    statusdual = solve(mdual)
    dualobjval = getobjectivevalue(mdual)
    while objval + tol < dualobjval
      @constraint(m, sum{(1-u[c]), c=1:nbMp; uval[c] == 1}  >= 1)
      cut_count[1] +=1
      status = solve(m)
      objval=getobjectivevalue(m)
      uval=getvalue(u)
      for c=1:nbMp
        JuMP.setRHS(mpsurplus[c], - bigM[c]*(1-uval[c]) - mp_headers[c,:FC] )
      end
      statusdual = solve(mdual)
      dualobjval = getobjectivevalue(mdual)

      println("Primal Obj:", objval)
      println("Dual Obj:", dualobjval)
    end
end

if method_type == "mp_relaxed"
  uval=getvalue(u)
  tol=1e-4
  for c=1:nbMp
    JuMP.setRHS(mpsurplus[c], - bigM[c]*(1-uval[c]) - mp_headers[c,:FC] )
  end
  statusdual = solve(mdual)
  dualobjval = getobjectivevalue(mdual)
end

uval_=getvalue(u)
xval_=getvalue(x)
xhval_=getvalue(xh)
fval_=getvalue(f)
sval_=getvalue(s)
scval_=getvalue(sc)
shmaxval_=getvalue(shmax)
shminval_=getvalue(shmin)
priceval_=getvalue(price)
vval_=getvalue(v)
mycandidate=damsol(xval_, uval_, xhval_, fval_, sval_, scval_, shmaxval_, shminval_, vval_, priceval_)

########################### Run stats
incomes_table = incomeCheck(mydata, mycandidate)
sol_quality = solQstats(mydata, mycandidate)
push!(solQ, [ssid sol_quality[1,1] sol_quality[1,2] sol_quality[1,3] sol_quality[1,4] sol_quality[1,5] sol_quality[1,6] sol_quality[1,7] sol_quality[1,8]  sol_quality[1,9]])
maxSlackViolation = sol_quality[1,1]

nbNodes = MathProgBase.getnodecount(m)
# Counting the number of cuts of all kinds generated by Cplex, safe user cuts which are internatlly numbered '15' :
nbSolverCuts = mapreduce(i->CPLEX.get_num_cuts(m.internalModel.inner, i), +, [1:14;]) + mapreduce(i->CPLEX.get_num_cuts(m.internalModel.inner, i), +, [16:18;]) # sole cplex dependent feature
absgap=MathProgBase.getobjbound(m) - getobjectivevalue(m) #absgap=CPLEX.getobjbound(getInternalModel(m)) - getObjectiveValue(m)
runtime = MathProgBase.getsolvetime(m)

push!(numTests, [ssid, objval, absgap, maxSlackViolation, cut_count[1], nbSolverCuts, nbNodes, runtime, nbMp, nbMpHourly + nbHourly] )

if method_type == "primal-dual"
  if mic_activation == 1
       method_type_export = string(method_type, "_omie")
  else method_type_export = string(method_type, "_mp")
  end
else method_type_export = method_type
end

inc_str = string("./incomes-", ssid, "_", method_type_export,".csv")
writetable(inc_str, incomes_table)
solqual_str = string("sol_quality_", method_type_export ,".csv")
writetable(solqual_str, solQ)
numtests_str = string("numTests_",method_type_export,".csv")
writetable(numtests_str, numTests)
end
