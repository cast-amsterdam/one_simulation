## loading needed Pkgs
using NetCDF,
LsqFit,
Interpolations,
Statistics,
DataFrames,
CSV,
ProgressBars,
SpecialFunctions,
Plots,
StatsBase,
JLD2,
Peaks,
MAT, 
Distributions,
Random

@. Gaus(x, p) =             p[1] .* exp.(.-0.5((x .- p[2]) ./ p[3]) .^ 2)

@. TripleGaus(x, p) =       (p[1] ) .* exp.(.-0.5((x .- p[2]) ./ p[3]) .^ 2) .+ 
                            (p[4] ) .* exp.(.-0.5((x .- p[5]) ./ p[6]) .^ 2) .+ 
                            (p[7] ) .* exp.(.-0.5((x .- p[8]) ./ p[9]) .^ 2)

@. SLN(x, p) =            p[1] .* (1 .+ ((x .- p[2]) .^ 2 ./ (p[4] .* (p[3] .+ p[5] .* (x .- p[2])) .^ 2))) .^ (.-1 .* p[4])

@. TripleSLN(x, p) = 
                            p[1] .* (1 .+ ((x .- p[2]) .^ 2 ./ (p[4] .* (p[3] .+ p[5] .* (x .- p[2])) .^ 2))) .^(.-1 .* p[4]) .+ 
                            p[6] .* (1 .+ ((x .- p[7]) .^ 2 ./ (p[9] .* (p[8] .+ p[10] .* (x .- p[7])) .^ 2))) .^ (.-1 .* p[9]) .+ 
                            p[11] .* (1 .+ ((x .- p[12]) .^ 2 ./ (p[14] .* (p[13] .+ p[15] .* (x .- p[12])) .^ 2))) .^ (.-1 .* p[14])

@. DoubleSLN(x, p) = 
                            p[1] .* (1 .+ ((x .- p[2]) .^ 2 ./ (p[4] .* (p[3] .+ p[5] .* (x .- p[2])) .^ 2))) .^(.-1 .* p[4]) .+ 
                            p[6] .* (1 .+ ((x .- p[7]) .^ 2 ./ (p[9] .* (p[8] .+ p[10] .* (x .- p[7])) .^ 2))) .^ (.-1 .* p[9])


function cntmd(pt,mdl,dp,dn)
    ## function for counting the number of modulations that form a peak. 
    ## Starting from the highest apex in the chromatogram the function will look for an apex in the adjacent modulations 
    ## and group them if they are lower then the current apex and fall in the correct time window. 
    
    crnt = pt.pn[findfirst(pt.hgt .== maximum(pt.hgt))] 
    pfnd = Vector{Int64}()
    append!(pfnd, crnt)

    xc = pt.x[pt.pn .== crnt][1]
    yc = pt.y[pt.pn .== crnt][1]
    i = 1 
    cnt = 1
    ppoi = deepcopy(crnt)

    while i .== cnt 
        ym = pt.y[pt.pn .== ppoi][1]
        xm = round(xc .- i*mdl,digits = 3)
        yl = pt.pn[pt.x .== xm]
        #poi = yl[ym-dn .< yl .< ym+dp]
        poi = yl[ym-dn .< vec(pt.y[indexin(yl,pt.pn)]) .< ym+dp]
        if length(poi) > 1
            poi = poi[findfirst(abs.(pt.y[indexin(poi,pt.pn)] .- pt.y[pt.pn .== ppoi]) .== minimum(abs.(pt.y[indexin(poi,pt.pn)] .- pt.y[pt.pn .== ppoi])))]
        elseif length(poi) == 1
            
        else
            break
        end

        if (pt.hgt[pt.pn .== ppoi] > pt.hgt[pt.pn .== poi]) .& (pt.y[pt.pn .== ppoi] < pt.y[pt.pn .== poi])
            cnt += 1
            append!(pfnd,poi[1])
        else 
            break
        end
        ppoi = deepcopy(poi)
        i += 1
    end

    xc = pt.x[pt.pn .== crnt][1]
    yc = pt.y[pt.pn .== crnt][1]
    i = 1 
    cnt = 1
    ppoi = deepcopy(crnt)

    while i .== cnt 
        ym = pt.y[pt.pn .== ppoi][1]
        xm = round(xc .+ i*mdl,digits = 6)
        yl = pt.pn[pt.x .== xm]
        #poi = yl[ym-dn .< yl .< ym+dp]
        poi = yl[ym-dn .< vec(pt.y[indexin(yl,pt.pn)]) .< ym+dp]
        if length(poi) > 1
            poi = poi[findfirst(abs.(pt.y[indexin(poi,pt.pn)] .- pt.y[pt.pn .== ppoi]) .== minimum(abs.(pt.y[indexin(poi,pt.pn)] .- pt.y[pt.pn .== ppoi])))]
        elseif length(poi) == 1
            
        else
            break
        end

        if (pt.hgt[pt.pn .== ppoi] > pt.hgt[pt.pn .== poi]) #.& (pt.y[pt.pn .== ppoi] > pt.y[pt.pn .== poi])
            cnt += 1
            append!(pfnd,poi[1])
        else 
            break
        end
        ppoi = deepcopy(poi)
        i += 1
    end
    return pfnd
end

function Coor2Lin(TimeLoc1D, TimeLoc2D, Time1D = Time1D, Time2D = Time2D)
    ## function to convert 1D and 2D time to linear Time
    ## 1D times are in min
    ## 2D times are in sec
        LinLoc = round.(((TimeLoc2D) .+ (floor.(TimeLoc1D./((ModuSize / freq) / 60)).* ((ModuSize / freq) / 60) .*60))./60, digits = 2)
    return LinLoc
end

function Lin2Coor(TimeLocLin, Time1D = Time1D, Time2D = Time2D)
    ## function to convert linear time to 1D and 2D
    ## 1D times are in min
    ## 2D times are in sec
    TimeLoc1D = round.((floor.(TimeLocLin ./ ((ModuSize / freq) / 60)) .* (ModuSize / freq) / 60), digits = 3)
    TimeLoc2D = round.((TimeLocLin .- TimeLoc1D).*60, digits = 3)
    return [TimeLoc1D,TimeLoc2D]
end

function Rebuild(PeakfitCL, xdata, ydata, peaks_postsel, model)
    #= used to rebuild a fitted chromatogram peak by peak. 
    as well as calculate the RMSE for each peak and the area
    =#
    println("Rebuilding data...")
    FitDataCumaCL = zeros(size(xdata))
    areas = zeros(size(PeakfitCL)[1])
    for i = ProgressBar(1:size(PeakfitCL, 1))
        FitDataCumaCL .+= model(xdata, Vector(PeakfitCL[i,2:end-1]))
        areas[i] = sum(model(xdata, Vector(PeakfitCL[i,2:end-1]))./mean(diff(xdata)))
    end
    println("   ")

    FitDataCumaCL = FitDataCumaCL .+ mean(ydata[1:100]) .- mean(FitDataCumaCL[1:100])

    RMSEs = Vector{Float64}()
    for i = 1:peaks_postsel.pg[end]
        peaks_postsel_st = peaks_postsel.strt[peaks_postsel.pg .== i]
        peaks_postsel_nd = peaks_postsel.nd[peaks_postsel.pg .== i]
        start_trace = minimum(peaks_postsel_st)-10
        end_trace = maximum(peaks_postsel_nd)+10
        time_axis = xdata[start_trace .< xdata .< end_trace]
        rawdata_trace = Float64.(ydata[start_trace .< xdata .< end_trace])
        simdata_trace = FitDataCumaCL[start_trace .< xdata .< end_trace]
        residual_trace = (rawdata_trace - simdata_trace)
        RMSE = StatsBase.rmsd(rawdata_trace,simdata_trace, normalize=true)
        append!(RMSEs,RMSE)
    end
    
    return FitDataCumaCL, RMSEs, areas
end

function GroupModulations(pt, dp, dn)
    #= function that fasiliates the application of cndmd() on data
    starts by appling cndmd() to the highest peak 
    Ones more apices are found the fit this criterium, 
    the grouped modulations are removed from the list and the function starts over. 
    =#
    mdl = ((ModuSize / freq) / 60) #diff(Time1D[1:2])[1] 
    x = deepcopy(pt)

    apex = Vector{Int64}()
    fpk = Vector{Any}()
    mdlsz = Vector{Int64}()
    while (size(x)[1] > 1) && x.sel[findfirst(x.hgt .== maximum(x.hgt))] == true
        append!(apex,x.pn[findfirst(x.hgt .== maximum(x.hgt))])
        pfnd = cntmd(x, mdl, dp, dn)
        append!(fpk,[pfnd])
        append!(mdlsz, length(pfnd))
        delete!(x,sort(indexin(pfnd, x.pn)))
    end

    apex = [pt.pn[indexin(apex,pt.pn)] apex pt.hgt[indexin(apex,pt.pn)] pt.wdt[indexin(apex,pt.pn)] mdlsz]

    psx = Vector{Any}()
    psy = Vector{Any}()

    for i = 1:length(fpk)
        append!(psx,[Lin2Coor(pt.lin[indexin(fpk[i],pt.pn)])[1]])
        append!(psy,[Lin2Coor(pt.lin[indexin(fpk[i],pt.pn)])[2]])
    end

    for i in 1:length(psx)
        prm = sortperm(psx[i])
        psx[i] = psx[i][prm]
        psy[i] = psy[i][prm]
    end
        return fpk, apex, psx, psy
end

function GaussianFit(pt,xdata,ydata,NumOfPeaks,fpk)
    #= will perform a Gaussian fit on each peak in the data 
    pt = 1D peak table
    xdata = time
    ydata = signal
    NumOfPeaks = the number of peaks you want to fit
    fpk = list of grouped peaks
    =#
    peaks_postsel = DataFrame()
    for k = 1:NumOfPeaks #selectign the required peaks
        cur = pt[indexin(fpk[k],pt.pn),:]
        pg = DataFrame(pg = Int64.(k .*(ones(size(cur)[1]))))
        pgcur = hcat(cur,pg)
        peaks_postsel = vcat(peaks_postsel, pgcur)
    end

    println("Fitting data...")
    Peakfit = DataFrame(pn = [], fh = [], fp = [], fw = [], adj_r2 = [])
    
    for i = ProgressBar(1:size(peaks_postsel)[1])
        # cutting out the region of intrest
        TimeSelect = xdata[peaks_postsel.strt[i].<xdata.<peaks_postsel.nd[i]]
        SignalSelect =
            ydata[peaks_postsel.strt[i].<xdata.<peaks_postsel.nd[i]]
        xdatapeak = TimeSelect
        ydatapeak = SignalSelect 
        
        #setting starting values and upper and lower bounderies 
        p0 = [peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i]]
        ub = [Inf, Inf, Inf]
        lb = [0.0, 0.0, 0.0]
        fit = curve_fit(Gaus, xdatapeak, ydatapeak, p0, lower = lb, upper = ub) #least squars fit of data 
        r2 = round(1 - var(fit.resid) / var(ydatapeak); digits = 3)     #adjusted Rsquerd
        push!(Peakfit, (peaks_postsel.pn[i], fit.param[1], fit.param[2], fit.param[3], r2))
    end
    println("   ")
    # checking for big deviations
    yreal = [peaks_postsel.hgt peaks_postsel.lin peaks_postsel.wdt]
    rdifh = abs.((Peakfit.fh .- yreal[:, 1]) ./ yreal[:, 1])
    rdifp = abs.((Peakfit.fp .- yreal[:, 2]) ./ yreal[:, 2])
    rdifw = abs.((Peakfit.fw .- yreal[:, 3]) ./ yreal[:, 3])

    PeakfitRM = (Peakfit.adj_r2 .< 0.01) .| (rdifw .> 10.0) .| (rdifp .> 10.0) .| (rdifh .> 10.0)
    PeakfitCL = Peakfit[.!PeakfitRM,:]
    RowCL = findall(.!PeakfitRM)
    RowRM = findall(PeakfitRM)
    pdifh = peaks_postsel.hgt[rdifp.<0.0001]

    println("$(size(RowRM)[1]) badly fitted peaks removed ($(round((size(RowRM)[1]/size(RowCL)[1])*100, digits = 2))%), $(size(RowCL)[1]) remaining")
    println("   ")

    return PeakfitCL, peaks_postsel
end

function TripleFitGaussianFit(pt,xdata,ydata,NumOfPeaks,fpk)
    #Same as Gaussian fit, only now with 3 peaks at once to better account for co-elution 
    peaks_postsel = DataFrame()
    peaks_postsel.pn = [0]
    peaks_postsel.lin = [0.0]
    peaks_postsel.strt = [0.0]
    peaks_postsel.nd = [0.0]
    peaks_postsel.x = [0.0]
    peaks_postsel.y = [0.0]
    peaks_postsel.hgt = [0.0]
    peaks_postsel.wdt = [0.0]
    peaks_postsel.sel = [true]
    peaks_postsel.pg = [0]
    
    for k = 1:NumOfPeaks
        cur = pt[indexin(fpk[k],pt.pn),:]
        pg = DataFrame(pg = Int64.(k .*(ones(size(cur)[1]))))
        pgcur = hcat(cur,pg)
        peaks_postsel = vcat(peaks_postsel, pgcur)
    end
    
    push!(peaks_postsel,(peaks_postsel.pn[end]+1, peaks_postsel.lin[end]+1, peaks_postsel.strt[end]+1, peaks_postsel.nd[end]+1, peaks_postsel.x[end]+1, peaks_postsel.y[end]+1, peaks_postsel.hgt[end]+1, peaks_postsel.wdt[end]+1, true, peaks_postsel.pg[end]+1))
    
    println("Fitting data...")
    Peakfit = DataFrame(pn = [], fh = [], fp = [], fw = [], adj_r2 = [])
    # @. TripleGaus(x, p) = p[1] .* exp.(.-((x .- p[2]) ./ p[3]) .^ 2) .+ 
    #                 p[4] .* exp.(.-((x .- p[5]) ./ p[6]) .^ 2) .+ 
    #                 p[7] .* exp.(.-((x .- p[8]) ./ p[9]) .^ 2)

    for i = ProgressBar(2:size(peaks_postsel)[1]-1)
        curpn = peaks_postsel.pn[i]
        # TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
        # SignalSelect =
        #     ydata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
        if 1 .< curpn .< maximum(peaks_postsel.pn)
            TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
            SignalSelect =
                ydata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
                p0 = Float64.([pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], pt.hgt[pt.pn.==curpn+1][1], pt.lin[pt.pn.==curpn+1][1], pt.wdt[pt.pn.==curpn+1][1]])
            elseif curpn == 1 
                TimeSelect = xdata[0 .<xdata.<pt.nd[pt.pn .==curpn+1]]
            SignalSelect =
                ydata[0 .<xdata.<pt.nd[pt.pn .==curpn+1]]
                p0 = Float64.([0,0,0, peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], pt.hgt[pt.pn.==curpn+1][1], pt.lin[pt.pn.==curpn+1][1], pt.wdt[pt.pn.==curpn+1][1]])
            elseif curpn == maximum(peaks_postsel.pn)
                TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.< xdata[end]]
            SignalSelect =
                ydata[pt.strt[pt.pn .== curpn-1].<xdata.<xdata[end]]
                p0 = Float64.([pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], 0, 0, 0])
            end
        xdatapeak = TimeSelect
        ydatapeak = SignalSelect .- 2
    
        # p0 = Float64.([pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], pt.hgt[pt.pn.==curpn+1][1], pt.lin[pt.pn.==curpn+1][1], pt.wdt[pt.pn.==curpn+1][1]])
        ub = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
        lb = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]#[peaks_postsel.hgt[i-1], peaks_postsel.lin[i-1], peaks_postsel.wdt[i-1], peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], peaks_postsel.hgt[i+1], peaks_postsel.lin[i+1], peaks_postsel.wdt[i+1]]#[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        # try 
            fit = curve_fit(TripleGaus, xdatapeak, ydatapeak, p0, lower = lb, upper = ub)
            r2 = round(1 - var(fit.resid) / var(ydatapeak); digits = 3)     #adjusted Rsquerd
            push!(Peakfit, (peaks_postsel.pn[i], fit.param[4], fit.param[5], fit.param[6], r2))
        # catch
        #     println(i)
        # end
    end

    println("   ")
    yreal = [peaks_postsel.hgt[2:end-1] peaks_postsel.lin[2:end-1] peaks_postsel.wdt[2:end-1]]
    rdifh = abs.((Peakfit.fh .- yreal[:, 1]) ./ yreal[:, 1])
    rdifp = abs.((Peakfit.fp .- yreal[:, 2]) ./ yreal[:, 2])
    rdifw = abs.((Peakfit.fw .- yreal[:, 3]) ./ yreal[:, 3])

    PeakfitRM = (Peakfit.adj_r2 .< 0.01) .| (rdifw .> 10.0) .| (rdifp .> 10.0) .| (rdifh .> 10.0)
    PeakfitCL = Peakfit[.!PeakfitRM,:]
    RowCL = findall(.!PeakfitRM)
    RowRM = findall(PeakfitRM)
    # pdifh = peaks_postsel.hgt[rdifp.<0.0001]

    println("$(size(RowRM)[1]) badly fitted peaks removed ($(round((size(RowRM)[1]/(size(RowRM)[1]+size(RowCL)[1]))*100, digits = 2))%), $(size(RowCL)[1]) remaining")
    println("   ")

    return PeakfitCL, peaks_postsel#, model
end

function SLNFit(pt,xdata,ydata,NumOfPeaks,fpk)
     #= will perform a SLN fit on each peak in the data 
    pt = 1D peak table
    xdata = time
    ydata = signal
    NumOfPeaks = the number of peaks you want to fit
    fpk = list of grouped peaks
    =#
    peaks_postsel = DataFrame()
    for k = 1:NumOfPeaks
        cur = pt[indexin(fpk[k],pt.pn),:]
        pg = DataFrame(pg = Int64.(k .*(ones(size(cur)[1]))))
        pgcur = hcat(cur,pg)
        peaks_postsel = vcat(peaks_postsel, pgcur)
    end
    println("Fitting data...")
    Peakfit = DataFrame(pn = [], fh = [], fp = [], fw = [], fm = [], fa = [], adj_r2 = [])
          for i = ProgressBar(1:size(peaks_postsel)[1])
        TimeSelect = xdata[peaks_postsel.strt[i].<xdata.<peaks_postsel.nd[i]]
        SignalSelect =
            ydata[peaks_postsel.strt[i].<xdata.<peaks_postsel.nd[i]]
        xdatapeak = TimeSelect
        ydatapeak = SignalSelect .- 2

        p0 = [peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i],
            5.0,
            0.0
        ]
        ub = [1300, peaks_postsel.lin[i] + 3*peaks_postsel.wdt[i], 0.5, 500.0, 0.3]
        lb = [0.000001, peaks_postsel.lin[i] - 3*peaks_postsel.wdt[i], 0.000001, 5.0, -0.17]
        try
            fit = curve_fit(SLN, xdatapeak, ydatapeak, p0, lower = lb, upper = ub)
            r2 = round(1 - var(fit.resid) / var(ydatapeak); digits = 3)     #adjusted Rsquerd
            push!(
                Peakfit,
                (
                    peaks_postsel.pn[i],
                    fit.param[1],
                    fit.param[2],
                    fit.param[3],
                    fit.param[4],
                    fit.param[5],
                    r2,
                ),
            )
        catch 
            push!(
                Peakfit,
                (
                    peaks_postsel.pn[i],
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                ),
            )
        end
    end
    println("   ")
    yreal = [peaks_postsel.hgt peaks_postsel.lin peaks_postsel.wdt]
    rdifh = abs.((Peakfit.fh .- yreal[:, 1]) ./ yreal[:, 1])
    rdifp = abs.((Peakfit.fp .- yreal[:, 2]) ./ yreal[:, 2])
    rdifw = abs.((Peakfit.fw .- yreal[:, 3]) ./ yreal[:, 3])
    
    PeakfitRM = (Peakfit.adj_r2 .< 0.01) .| (rdifw .> 10.0) .| (rdifp .> 10.0) .| (rdifh .> 10.0)
    PeakfitCL = Peakfit[.!PeakfitRM,:] 
    RowCL = findall(.!PeakfitRM)
    RowRM = findall(PeakfitRM)
    println("       $(size(RowRM)[1]) badly fitted peaks removed ($(round((size(RowRM)[1]/(size(RowCL)[1]+size(RowRM)[1]))*100, digits = 2))%), $(size(RowCL)[1]) remaining")
    println("   ")
    return PeakfitCL, peaks_postsel
end

function TripleFitSLN(pt,xdata,ydata,NumOfPeaks,fpk)
    peaks_postsel = DataFrame()
    peaks_postsel.pn = [0]
    peaks_postsel.lin = [0.0]
    peaks_postsel.strt = [0.0]
    peaks_postsel.nd = [0.0]
    peaks_postsel.x = [0.0]
    peaks_postsel.y = [0.0]
    peaks_postsel.hgt = [0.0]
    peaks_postsel.wdt = [0.0]
    peaks_postsel.sel = [true]
    peaks_postsel.pg = [0]
    
    for k = 1:NumOfPeaks
        cur = pt[indexin(fpk[k],pt.pn),:]
        pg = DataFrame(pg = Int64.(k .*(ones(size(cur)[1]))))
        pgcur = hcat(cur,pg)
        peaks_postsel = vcat(peaks_postsel, pgcur)
    end
    
    push!(peaks_postsel,(peaks_postsel.pn[end]+1, peaks_postsel.lin[end]+1, peaks_postsel.strt[end]+1, peaks_postsel.nd[end]+1, peaks_postsel.x[end]+1, peaks_postsel.y[end]+1, peaks_postsel.hgt[end]+1, peaks_postsel.wdt[end]+1, true, peaks_postsel.pg[end]+1))
    
    println("Fitting data...")
    Peakfit = DataFrame(pn = [], fh = [], fp = [], fw = [], fm = [], fa = [], adj_r2 = [])
        
        for i = ProgressBar(2:size(peaks_postsel)[1]-1)
            curpn = peaks_postsel.pn[i]
                if 1 .< curpn .< maximum(peaks_postsel.pn) 
                TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
                SignalSelect = ydata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
    
                p0 =    [pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], 500.0, 0.0,
                        peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], 500.0, 0.0,
                        pt.hgt[pt.pn.== curpn+1][1], pt.lin[pt.pn.== curpn+1][1], pt.wdt[pt.pn.== curpn+1][1], 500.0, 0.0]
                        
                ub =    [1300, pt.lin[pt.pn.== curpn-1][1]*1.0005, 60, 10000.0, 0.4, 
                        1300, peaks_postsel.lin[i]*1.0005, 60, 10000.0, 0.4, 
                        1300, pt.lin[pt.pn.== curpn+1][1]*1.0005, 60, 10000.0, 0.4]
                
                lb =    [0.1, pt.lin[pt.pn.== curpn-1][1]*0.9995, 0.00001, 50.0, -0.4, 
                        0.1, peaks_postsel.lin[i]*0.9995, 0.00001, 50.0, -0.4, 
                        0.1, pt.lin[pt.pn.== curpn+1][1]*0.9995, 0.00001, 50.0, -0.4]
    
            elseif curpn == 1
                TimeSelect = xdata[0 .<xdata.<pt.nd[pt.pn .==curpn+1]]
                SignalSelect = ydata[0 .<xdata.<pt.nd[pt.pn .==curpn+1]]
                
                p0 =    [0.1,xdata[1]+0.00001,0, 500.0, 0.0,
                        peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], 500.0, 0.0,
                        pt.hgt[pt.pn.== curpn+1][1], pt.lin[pt.pn.== curpn+1][1], pt.wdt[pt.pn.== curpn+1][1], 500.0, 0.0]
    
                ub =    [1300, peaks_postsel.lin[i]*0.9999, 60, 10000.0, 0.4, 
                        1300, peaks_postsel.lin[i]*1.000005, 60, 10000.0, 0.4, 
                        1300, pt.lin[pt.pn.== curpn+1][1]*1.000005, 60, 10000.0, 0.4]
                
                lb =    [0.1, xdata[1], 0.00001, 50.0, -0.4, 
                        0.1, peaks_postsel.lin[i]*0.999995, 0.00001, 50.0, -0.4, 
                        0.1, pt.lin[pt.pn.== curpn+1][1]*0.999995, 0.00001, 50.0, -0.4]
    
            elseif curpn == maximum(peaks_postsel.pn)
                TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.< xdata[end]]
                SignalSelect = ydata[pt.strt[pt.pn .== curpn-1].<xdata.<xdata[end]]
    
                p0 =    [pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], 500.0, 0.0,
                        peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], 500.0, 0.0,
                        0.1,xdata[end]-0.01,0.01, 500.0, 0.0]
    
                ub =    [1300, pt.lin[pt.pn.== curpn-1][1]*1.000005, 60, 10000.0, 0.4, 
                        1300, peaks_postsel.lin[i]*1.000005, 60, 10000.0, 0.4, 
                        1300, xdata[end], 60, 10000.0, 0.4]
                
                lb =    [0.1, pt.lin[pt.pn.== curpn-1][1]*0.999995, 0.00001, 50.0, -0.4, 
                        0.1, peaks_postsel.lin[i]*0.999995, 0.00001, 50.0, -0.4, 
                        0.1, peaks_postsel.lin[i]*0.999999, 0.00001, 50.0, -0.4]
            end
            xdatapeak = TimeSelect
            ydatapeak = SignalSelect .- 2

   
        fit = curve_fit(TripleSLN, xdatapeak, ydatapeak, p0, lower = lb, upper = ub)
        r2 = round(1 - var(fit.resid) / var(ydatapeak); digits = 3)     #adjusted Rsquerd
        push!(
            Peakfit,
            (
                peaks_postsel.pn[i],
                fit.param[6],
                fit.param[7],
                fit.param[8],
                fit.param[9],
                fit.param[10],
                r2,
            ),
        )
    end
    println("   ")
    yreal = [peaks_postsel.hgt[2:end-1] peaks_postsel.lin[2:end-1] peaks_postsel.wdt[2:end-1]]
    rdifh = abs.((Peakfit.fh .- yreal[:, 1]) ./ yreal[:, 1])
    rdifp = abs.((Peakfit.fp .- yreal[:, 2]) ./ yreal[:, 2])
    rdifw = abs.((Peakfit.fw .- yreal[:, 3]) ./ yreal[:, 3])
    
    PeakfitRM = (Peakfit.adj_r2 .< 0.01) .| (rdifw .> 1) .| (rdifp .> 1) .| (rdifh .> 1)
    PeakfitCL = Peakfit[.!PeakfitRM,:]
    RowCL = findall(.!PeakfitRM)
    RowRM = findall(PeakfitRM)
    println("       $(size(RowRM)[1]) badly fitted peaks removed ($(round((size(RowRM)[1]/(size(RowCL)[1]+size(RowRM)[1]))*100, digits = 2))%), $(size(RowCL)[1]) remaining")
    println("   ")
    return PeakfitCL, peaks_postsel

end

function AdaptiveSLN(pt,xdata,ydata,NumOfPeaks,fpk,Ms,Es)
    #will perform either a single, double, or tripple SLN fit based on the proximity of 
    peaks_postsel = DataFrame()
    peaks_postsel.pn = [0]
    peaks_postsel.lin = [0.0]
    peaks_postsel.strt = [0.0]
    peaks_postsel.nd = [0.0]
    peaks_postsel.x = [0.0]
    peaks_postsel.y = [0.0]
    peaks_postsel.hgt = [0.0]
    peaks_postsel.wdt = [0.0]
    peaks_postsel.sel = [true]
    peaks_postsel.pg = [0]
    
    for k = 1:NumOfPeaks
        cur = pt[indexin(fpk[k],pt.pn),:]
        pg = DataFrame(pg = Int64.(k .*(ones(size(cur)[1]))))
        pgcur = hcat(cur,pg)
        peaks_postsel = vcat(peaks_postsel, pgcur)
    end
    
    push!(peaks_postsel,(peaks_postsel.pn[end]+1, peaks_postsel.lin[end]+1, peaks_postsel.strt[end]+1, peaks_postsel.nd[end]+1, peaks_postsel.x[end]+1, peaks_postsel.y[end]+1, peaks_postsel.hgt[end]+1, peaks_postsel.wdt[end]+1, true, peaks_postsel.pg[end]+1))
    
    println("Fitting data...")
    Peakfit = DataFrame(pn = [], fh = [], fp = [], fw = [], fm = [], fa = [], adj_r2 = [])
    
        for i = ProgressBar(2:size(peaks_postsel)[1]-1)
            curpn = peaks_postsel.pn[i]
            
            if 1 .< curpn .< maximum(peaks_postsel.pn) 
                edge = [lastindex(pt.lin[pt.lin[pt.pn .== curpn] .- (6 .* pt.wdt[pt.pn .== curpn]) .< pt.lin .< pt.lin[pt.pn .== curpn] ]) .>= 1,
                        (lastindex(pt.lin[pt.lin[pt.pn .== curpn] .< pt.lin .< pt.lin[pt.pn .== curpn] .+ (6 .* pt.wdt[pt.pn .== curpn])]) .>= 1)*2]
                if sum(edge) .== 0
                    TimeSelect = xdata[pt.strt[pt.pn .== curpn].<xdata.<pt.nd[pt.pn .==curpn]]
                    SignalSelect = ydata[pt.strt[pt.pn .== curpn].<xdata.<pt.nd[pt.pn .==curpn]]
        
                    p0 =    [peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], Ms[2], Es[2]]
                    ub =    [maximum(pt.hgt)*1.01, peaks_postsel.lin[i]*1.0005, 50, Ms[3], Es[3]]
                    lb =    [minimum(pt.hgt)*0.99, peaks_postsel.lin[i]*0.9995, 0.00001, Ms[1], Es[1]]
                    mdl =   SLN
                elseif sum(edge) .== 1
                    TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn]]
                    SignalSelect = ydata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn]]
        
                    p0 =    [pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], Ms[2], Es[2],
                            peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], Ms[2], Es[2]]
                            
                    ub =    [maximum(pt.hgt)*1.01, pt.lin[pt.pn.== curpn-1][1]*1.0005, 100, Ms[3], Es[3], 
                    maximum(pt.hgt)*1.01, peaks_postsel.lin[i]*1.0005, 100, Ms[3], Es[3]]
                    
                    lb =    [minimum(pt.hgt)*0.99, pt.lin[pt.pn.== curpn-1][1]*0.9995, 0.00001, Ms[1], Es[1], 
                            minimum(pt.hgt)*0.99, peaks_postsel.lin[i]*0.9995, 0.00001,Ms[1], Es[1]]

                    mdl =   DoubleSLN
                elseif sum(edge) .== -2
                    TimeSelect = xdata[pt.strt[pt.pn .== curpn].<xdata.<pt.nd[pt.pn .==curpn+1]]
                    SignalSelect = ydata[pt.strt[pt.pn .== curpn].<xdata.<pt.nd[pt.pn .==curpn+1]]
        
                    p0 =    [peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], Ms[2], Es[2],
                            pt.hgt[pt.pn.== curpn+1][1], pt.lin[pt.pn.== curpn+1][1], pt.wdt[pt.pn.== curpn+1][1], Ms[2], Es[2]]
                            
                    ub =    [maximum(pt.hgt)*1.01, peaks_postsel.lin[i]*1.0005, 100, Ms[3], Es[3], 
                            maximum(pt.hgt)*1.01, pt.lin[pt.pn.== curpn+1][1]*1.0005, 100, Ms[3], Es[3]]
                    
                    lb =    [minimum(pt.hgt)*0.99, peaks_postsel.lin[i]*0.9995, 0.00001, Ms[1], Es[1], 
                            minimum(pt.hgt)*0.99, pt.lin[pt.pn.== curpn+1][1]*0.9995, 0.00001, Ms[1], Es[1]]
                    
                    mdl =   DoubleSLN
                else
                    TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
                    SignalSelect = ydata[pt.strt[pt.pn .== curpn-1].<xdata.<pt.nd[pt.pn .==curpn+1]]
        
                    p0 =    [pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], Ms[2], Es[2],
                            peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], Ms[2], Es[2],
                            pt.hgt[pt.pn.== curpn+1][1], pt.lin[pt.pn.== curpn+1][1], pt.wdt[pt.pn.== curpn+1][1], Ms[2], Es[2]]
                            
                    ub =    [maximum(pt.hgt)*1.01, pt.lin[pt.pn.== curpn-1][1]*1.0005, 100, Ms[3], Es[3], 
                            maximum(pt.hgt)*1.01, peaks_postsel.lin[i]*1.0005, 100, Ms[3], Es[3], 
                            maximum(pt.hgt)*1.01, pt.lin[pt.pn.== curpn+1][1]*1.0005, 100, Ms[3], Es[3]]
                    
                    lb =    [minimum(pt.hgt)*0.99, pt.lin[pt.pn.== curpn-1][1]*0.9995, 0.00001, Ms[1], Es[1], 
                            minimum(pt.hgt)*0.99, peaks_postsel.lin[i]*0.9995, 0.00001, Ms[1], Es[1], 
                            minimum(pt.hgt)*0.99, pt.lin[pt.pn.== curpn+1][1]*0.9995, 0.00001, Ms[1], Es[1]]
                    mdl =   TripleSLN
                end
            elseif curpn == 1
                TimeSelect = xdata[0 .<xdata.<pt.nd[pt.pn .==curpn+1]]
                SignalSelect = ydata[0 .<xdata.<pt.nd[pt.pn .==curpn+1]]
                
                p0 =    [minimum(pt.hgt)*0.99,xdata[1]+0.00001,0.00001,Ms[2], Es[2],
                        peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i], Ms[2], Es[2],
                        pt.hgt[pt.pn.== curpn+1][1], pt.lin[pt.pn.== curpn+1][1], pt.wdt[pt.pn.== curpn+1][1], Ms[2], Es[2]]
    
                ub =    [maximum(pt.hgt)*1.01, peaks_postsel.lin[i]*0.9999, 100, Ms[3], Es[3], 
                        maximum(pt.hgt)*1.01, peaks_postsel.lin[i]*1.000005, 100, Ms[3], Es[3], 
                        maximum(pt.hgt)*1.01, pt.lin[pt.pn.== curpn+1][1]*1.000005, 100, Ms[3], Es[3]]
                
                lb =    [minimum(pt.hgt)*0.99, xdata[1], 0.00001, Ms[1], Es[1], 
                        minimum(pt.hgt)*0.99, peaks_postsel.lin[i]*0.999995, 0.00001, Ms[1], Es[1], 
                        minimum(pt.hgt)*0.99, pt.lin[pt.pn.== curpn+1][1]*0.999995, 0.00001, Ms[1], Es[1]]
                mdl =   TripleSLN
    
            elseif curpn == maximum(peaks_postsel.pn)
                TimeSelect = xdata[pt.strt[pt.pn .== curpn-1].<xdata.< xdata[end]]
                SignalSelect = ydata[pt.strt[pt.pn .== curpn-1].<xdata.<xdata[end]]
    
                p0 =    [pt.hgt[pt.pn.== curpn-1][1], pt.lin[pt.pn.== curpn-1][1], pt.wdt[pt.pn.== curpn-1][1], Ms[2], Es[2],
                        peaks_postsel.hgt[i], peaks_postsel.lin[i], peaks_postsel.wdt[i],Ms[2], Es[2],
                        minimum(pt.hgt)*0.99,xdata[end]-0.01,0.01, Ms[2], Es[2]]
    
                ub =    [maximum(pt.hgt)*1.01, pt.lin[pt.pn.== curpn-1][1]*1.000005, 100, Ms[3], Es[3], 
                        maximum(pt.hgt)*1.01, peaks_postsel.lin[i]*1.000005, 100, Ms[3], Es[3], 
                        maximum(pt.hgt)*1.01, xdata[end], 100, Ms[3], Es[3]]
                
                lb =    [minimum(pt.hgt)*0.99, pt.lin[pt.pn.== curpn-1][1]*0.999995, 0.00001,Ms[1], Es[1], 
                        minimum(pt.hgt)*0.99, peaks_postsel.lin[i]*0.999995, 0.00001, Ms[1], Es[1], 
                        minimum(pt.hgt)*0.99, peaks_postsel.lin[i]*0.999999, 0.00001, Ms[1], Es[1]]
                mdl =   TripleSLN
            end
            xdatapeak = TimeSelect
            ydatapeak = SignalSelect

        fit = curve_fit(mdl, xdatapeak, ydatapeak, p0, lower = lb, upper = ub)
        r2 = round(1 - var(fit.resid) / var(ydatapeak); digits = 3)     #adjusted Rsquerd
        if (mdl == TripleSLN) || (sum(edge) .== 1)
            push!(
                Peakfit,
                (
                    peaks_postsel.pn[i],
                    fit.param[6],
                    fit.param[7],
                    fit.param[8],
                    fit.param[9],
                    fit.param[10],
                    r2,
                ),
            )
        elseif (mdl == SLN) || (sum(edge) .== -2) 
            push!(
                Peakfit,
                (
                    peaks_postsel.pn[i],
                    fit.param[1],
                    fit.param[2],
                    fit.param[3],
                    fit.param[4],
                    fit.param[5],
                    r2,
                ),
            )
        else
            println("ERROR: Unkwon Model")
        end

    end
    println("   ")
    yreal = [peaks_postsel.hgt[2:end-1] peaks_postsel.lin[2:end-1] peaks_postsel.wdt[2:end-1]]
    rdifh = abs.((Peakfit.fh .- yreal[:, 1]) ./ yreal[:, 1])
    rdifp = abs.((Peakfit.fp .- yreal[:, 2]) ./ yreal[:, 2])
    rdifw = abs.((Peakfit.fw .- yreal[:, 3]) ./ yreal[:, 3])
    
    #PeakfitRM = (Peakfit.adj_r2 .< 0.01) .| (rdifw .> 1) .| (rdifp .> 1) .| (rdifh .> 1)
    PeakfitRM = (Peakfit.adj_r2 .< -Inf) .| (rdifw .> Inf) .| (rdifp .> Inf) .| (rdifh .> Inf)
    PeakfitCL = Peakfit[.!PeakfitRM,:]
    PeakfitRMD = Peakfit[PeakfitRM,:]
    RowCL = findall(.!PeakfitRM)
    RowRM = findall(PeakfitRM)
        println("       $(size(RowRM)[1]) badly fitted peaks removed ($(round((size(RowRM)[1]/(size(RowCL)[1]+size(RowRM)[1]))*100, digits = 2))%), $(size(RowCL)[1]) remaining")
    println("   ")
    return PeakfitCL, peaks_postsel, PeakfitRMD

end

function SignalPrep(Time, Signal, ModTime; SignalBaseline=NaN, inspect = false)
    ## function will prepair signal for use by folding trimming and prosessing a baseline if applicable 
    ##Time in min
    ##ModTime in sec
    ##Signal in the unit it was recored in (e.g. pA, mAU etc.)
    if !isnan.(SignalBaseline)[1]
        freq = 1/(mean(diff(Time.*60)))
        ModuSize = Int(round(freq * ModTime))
        ModuRes = length(Signal) .- Int(floor(length(Signal)./ModuSize) .*ModuSize)
        TimeTrim = Time[1:end-ModuRes]
        SignalTrim = Signal[1:end-ModuRes]
        SignalBaselineTrim = Signal[1:end-ModuRes]
        if length(SignalTrim) .!== length(SignalBaselineTrim)
            print("ERROR: Baseline and Signal are different sizes")
        end
        Data2D = reshape(SignalTrim, ModuSize,:)
        Time2D = collect(range(0,stop = ModTime,length = size(Data2D,1)))
        Time1D = collect(range(TimeTrim[1],stop=TimeTrim[end],length = size(Data2D,2)))
        Baseline2D = reshape(SignalBaselineTrim, ModuSize,:)
        h = heatmap(Time1D,Time2D,Data2D,c = cgrad([:white,:navy,:indigo,:teal,:green,:yellow,:red],[0,0.2,1]), xlabel = "¹d Time (min)", ylabel = "²d time (sec)", title = "Raw Data")
        if inspect; display(plot(h)); end
        return Data2D, TimeTrim, SignalTrim, Time1D, Time2D, freq, ModuSize, ModuRes, SignalBaselineTrim, Baseline2D
    else
        freq = 1/(mean(diff(Time.*60)))
        ModuSize = Int(round(freq * ModTime))
        ModuRes = length(Signal) .- Int(floor(length(Signal)./ModuSize) .*ModuSize)
        TimeTrim = Time[1:end-ModuRes]
        SignalTrim = Signal[1:end-ModuRes]
        Data2D = reshape(SignalTrim, ModuSize,:)
        Time2D = collect(range(0,stop = ModTime,length = size(Data2D,1)))
        Time1D = collect(range(TimeTrim[1],stop=TimeTrim[end],length = size(Data2D,2)))
        h = heatmap(Time1D,Time2D,Data2D,colormap = [:white,:navy,:indigo,:teal,:green,:yellow,:red])#, xlabel = "¹d Time (min)", ylabel = "²d time (sec)", title = "Raw Data"))
        if inspect; display(plot(h)); end
        return Data2D, TimeTrim, SignalTrim, Time1D, Time2D, freq, ModuSize, ModuRes
    end
 
end

function findRoI(DataTrim, apex_th, peak_heigth_th;TimeTrim = TimeTrim)
    ## usess local derivitve based local maxima finder to locate posible peak locations an outputs a table of peaks of interst
    ppks, pvals = findmaxima(DataTrim)                                                  #prelim peak finding using local maxima
    vpks, vvals = findminima(DataTrim)   
    vpks = [1; vpks; length(DataTrim)]                                             #""
    ppks = ppks[pvals .> apex_th]                                                             #""
    pvals = pvals[pvals .> apex_th]                                                           #""
    pks, proms = peakproms(ppks, DataTrim)                                              #""
    pks, widths, leftedge, rightedge = peakwidths(pks, DataTrim, proms)       #""
    ppks = ppks[indexin(pks,ppks)]
    pvals = pvals[indexin(pks,ppks)]
    proms = proms[indexin(pks,ppks)]

    PrePks = DataFrame(strt=[], lin=[], nd=[], hgt=[],wdt=[])                           #compile info for fitting
    for i = 1:lastindex(ppks)
        right = vpks[findfirst(vpks.>ppks[i])]
        left = vpks[findlast(vpks.<ppks[i])]
        pk = ppks[i]
        hgt = pvals[i]
        wdt = widths[i]
        push!(PrePks,(left,pk,right,hgt,wdt))
    end

    locBool = (PrePks.hgt .> peak_heigth_th) 
    pt = DataFrame()
    pt.pn = collect(1:length(PrePks.lin))
    pt.lin = TimeTrim[PrePks.lin]
    pt.strt = TimeTrim[Int64.(round.(PrePks.strt))]
    pt.nd = TimeTrim[Int64.(round.(PrePks.nd))]
    pt.x = round.(Lin2Coor(TimeTrim[PrePks.lin])[1],digits=6)
    pt.y = round.(Lin2Coor(TimeTrim[PrePks.lin])[2],digits=6)
    pt.hgt = PrePks.hgt
    pt.wdt = (PrePks.wdt./(freq))./60
    pt.sel = locBool
    return pt
end 

function DepParaGrid(PeakfitCL, xnum, ynum)         
    ## used to sample probability distirbutions for imput in probchroma function 
    ## xnum and ynum are the number of "tiles" used to dived teh chromatogram
    ParaSel = DataFrame(linpos=[], pos1=[],pos2=[],hgt=[], wdt=[], m =[], e=[], mdlsz=[])
    xgrid = collect(0:(Lin2Coor(maximum(PeakfitCL.fp))[1]/xnum):Lin2Coor(maximum(PeakfitCL.fp))[1])
    ygrid = collect(0:((ModuSize/freq)/ynum):(ModuSize/freq))
    loc = Vector{Int64}(filter(x -> .!isnothing(x),indexin(apex[1:size(fpk)[1],1],PeakfitCL.pn)))
    apexfit = DataFrame()
    apexfit.pn = PeakfitCL.pn[loc]
    apexfit.hgt = PeakfitCL.fh[loc]
    apexfit.wdt = PeakfitCL.fw[loc]
    apexfit.m = PeakfitCL.fm[loc]
    apexfit.e = PeakfitCL.fa[loc]
    apexfit.mdlsz = apex[indexin(PeakfitCL.pn[loc],apex[:,1]),end]
    apexfit.adj_r2 = PeakfitCL.adj_r2[loc]
    apexfit.xloc = Lin2Coor(PeakfitCL.fp)[1][loc]
    apexfit.yloc = Lin2Coor(PeakfitCL.fp)[2][loc]
    apexfit.lin = PeakfitCL.fp[loc]

    apexfit = apexfit[apexfit.xloc .< 110,:]
    for k = 1:length(xgrid)-1
        for l = 1:length(ygrid)-1
            x1 = xgrid[k]
            x2 = xgrid[k+1]
            y1 = ygrid[l]
            y2 = ygrid[l+1]
            Temp =  apexfit[x1 .<= apexfit[!,"xloc"].<= x2,:]
            TempParaSel =  Temp[y1 .<= Temp[!,"yloc"].<= y2,:]

            if isempty(TempParaSel)
                continue
            else
                 
    
      
                pos1 = rand(Uniform(x1,x2))
                pos2 = rand(Uniform(y1,y2))
                poslin = Coor2Lin(pos1,pos2)
                    push!(ParaSel,(
                # rand(collect(1:0.1:150)),
                # rand(collect(0.001:0.0001:0.1)),
                poslin,
                pos1,
                pos2,
                TempParaSel[rand(collect(1:size(TempParaSel)[1])),"hgt"],
                TempParaSel[rand(collect(1:size(TempParaSel)[1])),"wdt"],
                TempParaSel[rand(collect(1:size(TempParaSel)[1])),"m"],
                TempParaSel[rand(collect(1:size(TempParaSel)[1])),"e"],
                TempParaSel[rand(collect(1:size(TempParaSel)[1])),"mdlsz"]))
            end
        end
    end
    return ParaSel
end

function probchroma(PeakfitCL)                                   
    #this function applies the probabilty distributions and will therefor rebuild a uniq, yet statisticaly significant chromatogram
    #spacial thanks to Dr. Leon Niezen for providen the source code form his paper https://doi.org/10.1016/j.aca.2022.339605
    #and helping with expending it into the second demention that later became this function 

    loc = Vector{Int64}(filter(x -> .!isnothing(x),indexin(apex[1:size(fpk)[1],1],PeakfitCL.pn)))
    apexfit = DataFrame()
    apexfit.pn = PeakfitCL.pn[loc]
    apexfit.hgt = PeakfitCL.fh[loc]
    apexfit.wdt = PeakfitCL.fw[loc]
    apexfit.m = PeakfitCL.fm[loc]
    apexfit.e = PeakfitCL.fa[loc]
    apexfit.mdlsz = apex[indexin(PeakfitCL.pn[loc],apex[:,1]),end]
    apexfit.adj_r2 = PeakfitCL.adj_r2[loc]
    apexfit.xloc = Lin2Coor(PeakfitCL.fp)[1][loc]
    apexfit.yloc = Lin2Coor(PeakfitCL.fp)[2][loc]
    apexfit.lin = PeakfitCL.fp[loc]

    PGP = DepParaGrid(PeakfitCL, 150,75)
    
    Peakshapes_1D = ones(size(PGP)[1]).*6               #See reference table above
    Peakshapes_2D = ones(size(PGP)[1]).*6               #See reference table above
    Noise_sd = 0.005                                    #Noise standaard deviatiom
    Noise_Distribution = "Normal"                       #Noise distribution type
    Noise_Extra = 1                                     #Noise extra
    Noise_Type = 2                                      #Noise type, see list above
    Hgt_1D = PGP.hgt                                    #Peak Height
    Hgt_2D = PGP.hgt                                    #Peak Height
    Wdt_1D = PGP.mdlsz .* 0.12 .* 0.25                  #half Peak width
    Wdt_2D = PGP.wdt.*60                                #half Peak width
    m_1D = PGP.m                                        #m = some number (5 for example)
    m_2D = PGP.m                                        #m = some number (5 for example)
    As_1D = zeros(size(PGP)[1])                         #As = Assymetry (extend of tailing, typical values of 0.1-0.2)
    As_2D = PGP.e                                       #As = Assymetry (extend of tailing, typical values of 0.1-0.2)
    Pos_1D = PGP.pos1                                   #peak position (tr)
    Pos_2D = PGP.pos2                                   #peak position (tr)
    PosShift_2D = ones(size(PGP)[1]).*-0.12             #Amound of 2D shift (in )
    PosShiftLeftFrac_2D = [0.4]                         #fraction of left side asymmetry (0.0-1.0 , 0.5 is symmetrical)
    gD_1D = ones(size(PGP)[1])                          #gD is the Doppler (Gaussian) width
    gD_2D = ones(size(PGP)[1])                          #gD is the Doppler (Gaussian) width
    alpha_1D = ones(size(PGP)[1])                       #Alpha is the shape constant (ratio of the Lorentzian width gL to the Doppler width gD
    alpha_2D = ones(size(PGP)[1])                       #Alpha is the shape constant (ratio of the Lorentzian width gL to the Doppler width gD
    
    ##creates a shift in the modulations
    a = PosShift_2D .* ones(size(PGP)[1], length(Time_1D))
    x = Time_1D' .* ones(size(PGP)[1], length(Time_1D))
    fact = transpose(a .* (x .- Pos_1D))
    
    added_noise = GenerateNoise(
            Noise_Type,
            Noise_Distribution,
            TimeTot,
            Noise_sd,
            Noise_Extra,
        )
    
    added_noise = reshape(added_noise[1:end-ModuRes],ModuSize,:) #create addidtional noise

    Max_Hgt_1D,
    Min_Hgt_1D,
    Max_Hgt_Outlier_1D,
    Min_Hgt_Outlier_1D,
    Max_Wdt_1D,
    Min_Wdt_1D,
    Max_m_1D,
    Min_m_1D,
    Max_As_1D,
    Min_As_1D,
    Peakdata_1D =
        calcinput(Hgt_1D, Wdt_1D, m_1D, As_1D, Pos_1D, gD_1D, alpha_1D)

    peaks_1D, PeaksS_1D =
        ManualGeneratePeaks(Peakshapes_1D, Time_1D, Peakdata_1D)    #create 1D signals

   
        chromas = zeros(length(Time_1D), length(Time_2D))
        for i = ProgressBar(1:length(PeaksS_1D))
            Max_Hgt_2D,
            Min_Hgt_2D,
            Max_Hgt_Outlier_2D,
            Min_Hgt_Outlier_2D,
            Max_Wdt_2D,
            Min_Wdt_2D,
            Max_m_2D,
            Min_m_2D,
            Max_As_2D,
            Min_As_2D,
            Peakdata_2D = calcinput(
                peaks_1D[i, :],
                Wdt_2D,
                m_2D,
                As_2D,
                (Pos_2D .+ fact[i, :]),
                gD_2D,
                alpha_2D,
            )
            peaks_2D, PeaksS_2D =
                ManualGeneratePeaks(Peakshapes_2D, Time_2D, Peakdata_2D)        #expand the 1D signals into the 2D
            chromas[i, :] = PeaksS_2D
        end

        chroma_2D = (chromas' .+ Signal_chroma .+ added_noise) # Add layers togetter
      
     return chroma_2D, PGP
    
end

## the folloing fuctions are auxilary to probchroma()
function gaussian(Time, Hgt, Wdt, Pos)
    g = Hgt .* exp.(.-((Time .- Pos) ./ Wdt) .^ 2)
end

function lorentzian(Time, Hgt, Wdt, Pos)
    # lorentzian(x,position,width) Lorentzian function.
    # where x may be scalar; vector; | matrix
    # position & width scalar
    # T. C. O'Haver; 1988
    # Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
    g = Hgt .* (ones(size(Time)) ./ (1 .+ ((Time .- Pos) ./ (0.21 .* Wdt)) .^ 2))
end

function logistic(Time, Hgt, Wdt, Pos)
    # logistic function.  pos=position; wid=half-width [both scalar]
    # logistic(x,pos,wid), where x may be scalar, vector, | matrix
    # pos=position; wid=half-width [both scalar]
    # T. C. O'Haver; 1991
    n = exp.(-((Time .- Pos) ./ (Wdt)) .^ 2)
    g = Hgt .* ((2 .* n) ./ (1 .+ n))
end

function lognormal(Time, Hgt, Wdt, Pos)
    # lognormal function.  pos=position; wid=half-width [both scalar]
    # lognormal(x,pos,wid), where x may be scalar, vector, | matrix
    # pos=position; wid=half-width [both scalar]
    # T. C. O'Haver; 1991
    g = Hgt .* exp.(-(log.(Time ./ Pos) ./ (0.125 .* Wdt)) .^ 2)
end

function triangular(Time, Hgt, Wdt, Pos)
    #Triangle function.  pos=position; wid=half-width [both scalar]
    #trianglar[x,pos,wid], where x may be scalar | vector
    #pos=position; wid=half-width [both scalar]
    # T. C. O'Haver; 1991
    # Example
    # x=[0:.1:10];plot(x,trianglar[x,5.5,2.3],'.')
    g = Hgt .* 1 .- (1 ./ (1.7 .* Wdt)) .* abs.(Time .- Pos)
    for i = 1:length(Time)
        if g[i] < 0
            g[i] = 0
        end
    end
    return g
end

function rectangle(Time, Hgt, Wdt, Pos)
    #rectangle function.  pos=position; wid=half-width [both scalar]
    #rectangle(x,pos,wid), where x may be scalar | vector
    #pos=position; wid=half-width [both scalar]
    # T. C. O'Haver; 2016
    # Example
    # x=[0:.1:10];plot(x,rectangle(x,5.5,2.3),'.')
    g = zeros(size(Time))
    hw = 4.25 * Wdt ./ 2
    for i = 1:length(Time)
        if Time[i] < Pos .- hw
            g[i] = 0
        end
        if Time[i] > Pos .- hw
            g[i] = Hgt
        end
        if Time[i] > Pos .+ hw
            g[i] = 0
        end
    end
    return g
end

function exppulse(Time, Hgt, Wdt, Pos)
    # Exponential pulse of the form
    # g = (x-spoint)./pos.*exp(1-(x-spoint)./pos)
    e = (Time .- Pos) ./ (0.7 .* Wdt)
    p = 4 .* exp.(-e) .* (1 .- exp.(-e))
    p[p<=0] = p
    g = Hgt .* p'
end

function alphafunction(Time, Hgt, Wdt, Pos)
    # alpha function.  pos=position; wid=half-width [both scalar]
    # alphafunction(x,pos,wid), where x may be scalar, vector, | matrix
    # pos=position; wid=half-width [both scalar]
    # Taekyung Kwon; July 2013
    g =
        Hgt .* (Time .- Pos) ./ (0.55 .* Wdt) .*
        exp.(1 .- (Time .- Pos) ./ (0.55 .* Wdt))
    for i = 1:length(Time)
        if g[i] < 0
            g[i] = 0
        end
    end
    return g
end

function ngaussian(Time, Hgt, Wdt, Pos, m)
    #  ngaussian(x,pos,wid) = flattened Gaussian centered on x=pos, half-width=wid
    #  x may be scalar; vector; | matrix; pos & wid both scalar
    # Shape is Gaussian when n=1. Becomes more rectangular as n increases.
    #  T. C. O'Haver; 1988; revised 2014
    # Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
    if m > 0
        g = (1 .- (10 .^ -(m .* gaussian(Time, Hgt, 0.85 .* Wdt, Pos))))
        g = g ./ maximum(g)
    else
        ()
        g = gaussian(Time, Hgt, Wdt, Pos)
    end
end

function pearson(Time, Hgt, Wdt, Pos, m)
    # Pearson VII function.
    # g = pearson(x,pos,wid,m) where x may be scalar, vector, | matrix
    # pos=position; wid=half-width [both scalar]
    # m=some number
    #  T. C. O'Haver; 1990
    g =
        Hgt .* ones(size(Time)) ./
        (1 .+ ((Time .- Pos) ./ ((0.5 .^ (2 ./ m)) .* 2 .* Wdt)) .^ 2) .^ m
end

function modpearson(Time, Hgt, Wdt, Pos, m, As)
    # Modified Pearson VII function
    # As = Assymetry (extend of tailing, typical values of 0.1-0.2)
    # m = some number [5 for example]

    y =
        Hgt .*
        (
            1 .+
            ((Time .- Pos) .^ 2 ./ (m .* (Wdt .+ As .* (Time .- Pos)) .^ 2))
        ) .^ (.-1 .* m)
end

function expgaussian(Time, Hgt, Wdt, Pos, m)
    #  Exponentially-convoluted gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
    #  x may be scalar; vector; | matrix; pos & wid both scalar
    #  T. C. O'Haver; 2006
    g = Hgt .* exp.(-((Time .- Pos) ./ (0.60056120439323 .* Wdt)) .^ 2)
    g = ExpBroaden(g', m)
end

function explorentzian(Time, Hgt, Wdt, Pos, m)
    #  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
    #  x may be scalar; vector; | matrix; pos & wid both scalar
    #  T. C. O'Haver; 2013
    g = Hgt .* ones(size(Time)) ./ (1 .+ ((Time .- Pos) ./ (0.1 .* Wdt)) .^ 2)
    g = ExpBroaden(g', m)
end

function GL(Time, Hgt, Wdt, Pos, m)
    # Gaussian/Lorentzian blend. m = percent Gaussian character
    # pos=position; wid=half-width[]
    # m = percent Gaussian character.
    #  T. C. O'Haver; 2012
    # sizex=size(x)
    # sizepos=size(pos)
    # sizewid=size(wid)
    # sizem=size(m)
    g =
        2 .* (
            (m ./ 100) .* gaussian(Time, Hgt, Wdt, Pos) .+
            (1 .- (m[1] ./ 100)) .* lorentzian(Time, Hgt, Wdt, Pos)
        ) ./ 2
end

function BiGaussian(Time, Hgt, Wdt, Pos, m)
    # BiGaussian (different width on leading edge & trailing edge).
    # pos=position; wid=width
    # m = ratio of widths of right-hand to left-hand halves.
    # If m=1; it becomes identical to a Gaussian.
    # Verison 2; T. C. O'Haver; 2012
    # Example: Plots Gaussian & BiGaussian with m=3; both with halfwidth=20.
    # x=[1:100]
    # g=gaussian(x,50,20)
    # bg=bigaussian[x,50,20,3]
    # plot(x,g,x,bg)
    #
    lx = length(Time)
    dif = abs.(Time .- Pos)
    index = findall((dif .- minimum(dif)) == 0)
    hx = index

    hwid = 2 .* Wdt
    g[1:hx] = gaussian(Time[1:hx], Hgt, hwid ./ (m + 1), Pos)
    g[hx+1:lx] = gaussian(Time[hx+1:lx], Hgt, m .* hwid ./ (m + 1), Pos)
end

function BWF(Time, Hgt, Wdt, Pos, m)
    # BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
    # pos=position; wid=width; m=Fano factor()
    #  T. C. O'Haver; 2014
    y =
        ((m .* Wdt ./ 2 .+ Time .- Pos) .^ 2) ./
        (((Wdt ./ 2) .^ 2) .+ (Time .- Pos) .^ 2)
    # y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2)
    g = Hgt .* y ./ maximum(y)
end

function voigt(Time, Hgt, gD, alpha, Pos)
    # Voigt profile function. xx is the independent variable (energy
    # wavelength, etc), gD is the Doppler [Gaussian] width, & alpha is the
    # shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
    # Based on Chong Tao's "Voigt lineshape spectrum simulation"
    # File ID: #26707
    # alpha=alpha()
    gL = alpha .* gD
    gV = 0.5346 .* gL .+ sqrt(0.2166 .* gL .^ 2 .+ gD .^ 2)
    x = gL ./ gV
    # sizeabs=size(abs(xx-pos))
    # sizegV=size(gV)
    y = abs.(Time .- Pos) ./ gV
    g =
        1 / (2 .* gV .* (1.065 .+ 0.447 .* x .+ 0.058 .* x.^2)) .* (
            (1 .- x) * exp.(-0.693 .* y .^ 2) .+
            (x ./ (1 .+ y .^ 2)) .+
            0.016 .*
            (1 .- x) .*
            x .*
            (exp.(-0.0841 .* y .^ 2.25) .- 1 ./ (1 .+ 0.021 .* y .^ 2.25))
        )
    return g = Hgt .* g ./ maximum(g)
end

function ExpBroaden(y, t)
    # ExpBroaden(y,t) zero pads y & convolutes result by an exponential decay
    # of time constant t by multiplying Fourier transforms & inverse
    # transforming the result.
    hly = trunc(Int,round(length(y) ./ 2))
    ey = [y[1] .* ones(1, hly)'; y; y[length(y)] .* ones(1,hly)']
    fy = fft(ey)
    a = exp(-(1:length(fy)) ./ (t .* 100))
    fa = fft(a)
    fy1 = fy .* fa'
    ybz = real(ifft(fy1)) ./ sum(a)
    yb = ybz[hly+2:length(ybz)-hly+1]
end

function ManualGeneratePeaks(Peakshapes, Time, Peakdata)

    Hgts = Peakdata[1, :]
    Wdts = Peakdata[2, :]
    Poss = Peakdata[3, :]
    mp = Peakdata[5, :]
    Asp = Peakdata[6, :]
    gDp = Peakdata[7, :]
    alphap = Peakdata[8, :]

    n_m = 1                                            #Used to allow for multiple different peak types in 1 signal
    n_As = 1                                           #Used to allow for multiple different peak types in 1 signal
    n_v = 1                                            #Used to allow for multiple different peak types in 1 signal
    Temp = ones(8, 1)                                  #Used to pass single m, As, gD & alfa to determine width function

    peaks = zeros(length(Time), length(Peakshapes))
    for n = 1:size(Peakdata, 2)
        Peakshape = Peakshapes[n]

        if Peakshape in [5 6 7 8 9 12 13 16 17]
            Temp[5] = mp[n_m]
            m = mp[n_m]
            n_m = n_m + 1
        end

        if Peakshape .== 6
            Temp[6] = Asp[n_As]
            As = Asp[n_As]
            n_As = n_As + 1
        end

        if Peakshape .== 12
            Temp[7] = gDp[n_v]
            Temp[8] = alphap[n_v]
            gD = gDp[n_v]
            alpha = alphap[n_v]
            n_v = n_v + 1
        end

        #AcSigWdt = DetermineActualWidth[Time,Peakshape,Temp]
        #idx = find((abs(AcSigWdt - Wdts[n]) .== min(abs(AcSigWdt - Wdts[n]))))

        Wdt = Wdts[n]
        Hgt = Hgts[n]
        Pos = Poss[n]

        if Peakshape == 1
            peak = gaussian(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 2
            peak = lorentzian(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 3
            peak = logistic(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 4
            peak = lognormal(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 5
            peak = pearson(Time, Hgt, Wdt, Pos, m)
        elseif Peakshape == 6
            peak = modpearson(Time, Hgt, Wdt, Pos, m, As)
        elseif Peakshape == 7
            peak = GL(Time, Hgt, Wdt, Pos, m)
        elseif Peakshape == 8
            peak = BiGaussian(Time, Hgt, Wdt, Pos, m)
        elseif Peakshape == 9        #Flattened gaussian()
            peak = ngaussian(Time, Hgt, Wdt, Pos, m)
        elseif Peakshape == 10
            peak = alphafunction(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 11
            peak = exppulse(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 12
            peak = voigt(Time, Hgt, gD, alpha, Pos)
        elseif Peakshape == 13
            peak = BWF(Time, Hgt, Wdt, Pos, m)
        elseif Peakshape == 14
            peak = triangular(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 15
            peak = rectangle(Time, Hgt, Wdt, Pos)
        elseif Peakshape == 16
            peak = expgaussian(Time, Hgt, Wdt, Pos, m)
        elseif Peakshape == 17
            peak = explorentzian(Time, Hgt, Wdt, Pos, m)
        else
            println("ERROR: Invalid input for Peakshape!")
            println("Input must be int. between 1 and 17")
        end
        peaks[:, n] = peak[:, 1]
    end
    PeaksS = sum(peaks, dims = 2)
    return peaks, PeaksS
end

function laprnd(m, n, mu, sigma)
    # Generate Laplacian noise
    mu = 0
    sigma = 1
    u = rand(m, n) - 0.5
    b = sigma / sqrt(2)
    y = mu - b * sign(u) .* log(1 - 2 * abs(u))
end

function whitenoise(x, Noise_Distribution, Noise_Extra)
    # Random noise with white power spectrum with mean zero
    # & unit standard deviation; equal in length to x
    # Tom O'Haver; 2008
    # Example:
    # model=gaussian[[1:10deriv00],500,200]
    # model=model+.1.*whitenoise(model)
    # plot(model)

    N = length(x)

    if Noise_Distribution == "Normal"
        #pd = makedist["Normal"]
        pd = Normal()
    elseif Noise_Distribution == "Uniform"
        #pd = makedist["Uniform"];
        pd = Uniform()
    elseif Noise_Distribution == "Laplacian"
        #y = laprnd(N, 1)
        y = Laplace(N, 1)
    elseif Noise_Distribution == "Poisson"
        #pd = makedist["Poisson"]
        pd = Poisson()
    elseif Noise_Distribution == "Gamma"
        #pd = makedist["Gamma"];
        pd = Gamma()
    elseif Noise_Distribution == "Rayleigh"
        #pd = makedist["Rayleigh"]
        pd = Rayleigh()
    elseif Noise_Distribution == "Logistic"
        #pd = makedist["Logistic"];
        pd = Logistic()
    elseif Noise_Distribution == "Loglogistic"
        #pd = makedist["Loglogistic"];
        println("loglogistic is not suported in the current version")
    elseif Noise_Distribution == "Lognormal"
        #pd = makedist["Lognormal"];
        pd = LogNormal()
    elseif Noise_Distribution == "Exponential"
        #pd = makedist["Exponential"];
        pd = Exponential()
    else
        println("ERROR: Invalite Distribution Name")
    end

    if Noise_Distribution == "Laplacian"

    else
        y = rand(pd, N)
    end

    y = y'
end

function bluenoise(n, Noise_Distribution, Noise_Extra)

    # Random noise with blue power spectrum with mean zero
    # & unit standard deviation. n is number of points.
    # Tom O'Haver; 2008
    # Example:
    # model=gaussian[[1:1000],500,200]
    # model=model+.1.*bluenoise(length(model))
    # plot(model)
    x = 1:n
    y = whitenoise(x, Noise_Distribution, Noise_Extra)  # Random normally-distributed white noise
    #by=gradient(y)
end

function pinknoise(n, Noise_Distribution, Noise_Extra)

    # Random noise with pink [1/f] power spectrum with mean zero
    # & unit standard deviation. n is number of points.
    # Tom O'Haver; 2008
    # Example:
    # model=gaussian[[1:1000],500,200]
    # model=model+.1.*pinknoise(length(model))
    # plot(model)
    x = [1:n]
    y = whitenoise(x, Noise_Distribution, Noise_Extra)  # Random normally-distributed white noise
    # Fourier filter()
    fy = fft(y) # Compute Fourier transform of signal y
    # Compute filter shape
    lft1 = [1:(length(fy)/2)+1]
    lft2 = [(length(fy)/2):length(fy)]
    ffilter1 = ones(size(lft1)) ./ (sqrt(lft1))
    ffilter2 = ones(size(lft2)) ./ (sqrt(lft2))
    ffilter = [ffilter1, ffilter2]
    if length(fy) > length(ffilter)
        ffilter = [ffilter ffilter[1]]
    end
    ffy = fy .* ffilter[1:length(fy)]  # Multiply filter by Fourier transform of signal
    py = real(ifft(ffy)) # Inverse transform to recover filtered signal "ry"
end

function violetnoise(x, Noise_Distribution, Noise_Extra)

    # Random noise with high-frequency weighted power spectrum with mean zero
    # & unit standard deviation. Length equal to the length of x
    y = 0:x
    for n = 1:2:length(y)-1
        rn = abs(randn())
        y[n] = rn
        y[n+1] = -rn
    end
    y = y[1:end-1]
end

function rednoise(x, Noise_Distribution, Noise_Extra)

    ry =
        cumsum(whitenoise(x, Noise_Distribution, Noise_Extra)) /
        std(whitenoise(x, Noise_Distribution, Noise_Extra))
end

function GenerateNoise(
    Noise_Type,
    Noise_Distribution,
    Time,
    Noise_sd,
    Noise_Extra,
)
    #Generates either white, blue, pink, violet | red [brown] noise.
    #Random values are pulled from a specified distribution e.g. normal
    #distribution | lognormal distribution etc.

    #Allows for creation of for example "white gaussian" noise | e.g. "pink
    #poisson" noise etc. Note that due to central limit theorem the
    #distribution to pull from in reality will generally be gaussian [normally
    #distributed] although red noise | pink noise can be used to create
    #noise that also contains drift to make a more challenging test set.

    x = Time
    n = length(x)

    if Noise_Type == 1
        y = whitenoise(x, Noise_Distribution, Noise_Extra)
    elseif Noise_Type == 2
        y = bluenoise(n, Noise_Distribution, Noise_Extra)
    elseif Noise_Type == 3
        y = pinknoise(n, Noise_Distribution, Noise_Extra)
    elseif Noise_Type == 4
        y = violetnoise(n, Noise_Distribution, Noise_Extra)
    elseif Noise_Type == 5
        y = rednoise(x, Noise_Distribution, Noise_Extra)
    else
        println("ERROR: invalite input")
        println("Noise_Type should be between 1 and 5")
    end

    y = Noise_sd * ((y .- mean(y)) / std(y)) #Normalize all noise
    return y = y'
end