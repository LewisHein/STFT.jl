module STFT
"""
Compute the short-time Fourier transform of a signal
Parameters:
    -- data: The signal to find the STFT of
    -- windowFunc: An array sampled from the window function to be applied before taking the FFTs. The length of this function defines the window length
    -- windowStep: The step ( in samples) to move the window between FFTs
    -- low: The lowest harmonic to include in the output
    -- high: The highest harmonic to include in the output
    -- transformFunction: A function to apply to the windowed, FFT'd pieces of _data_
    -- processingFunction: A function to do the processing. This function is passed
       the window index and the windowed, FFT'd piece and must return a value of type
       outType. For instance, to compute the amplitude along some path, processingFunction should be (dat, i)->dat[path[i]]
"""
function stft_processed{dataType,windowType,outType}(data::AbstractArray{dataType, 1}, windowFunc::AbstractArray{windowType, 1}, windowStep::Integer, low::Integer, high::Integer, transformFunction::Function, processingFunction::Function, outT::Type{outType})
    dataLength = size(data, 1)
    windowLength = size(windowFunc)[1]
    nHarmonics = 1 + high - low
    nSnippets = Int(floor((dataLength-windowLength)/windowStep))
    endSample = nSnippets*windowStep
    if low > windowLength
	  throw(ArgumentError("The low harmonic (parameter 'low') must be less than the length of the window"))
    end
    if low < 1
	    throw(ArgumentError("The low harmonic (parameter 'low') must be at least 1"))
    end
    STFT_processed = Array{outType, 1}(nSnippets)
    temp = Array{Float64, 1}(windowLength)
    tempFFT = Array{Float64, 1}(windowLength)
    tempFFTProcessed = Array{Float64, 1}(nHarmonics)
    windowPos = 1
    harmonic = 1
    plan = plan_rfft(temp, flags = FFTW.MEASURE)

    for sample in 1:windowStep:endSample
        temp = windowFunc.*view(data, sample:(sample+windowLength-1))
        tempFFT = plan*temp

        harmonic = 1
        for i in low:high
            tempFFTProcessed[harmonic] = transformFunction(tempFFT[i].re, tempFFT[i].im)
	    STFT_processed[windowPos] = processingFunction(tempFFTProcessed, windowPos)
            harmonic += 1
        end

        windowPos += 1
    end

    return STFT_processed
end

"""
Compute the short-time Fourier transform of a signal
Parameters:
    -- data: The signal to find the STFT of
    -- windowFunc: An array sampled from the window function to be applied before taking the FFTs. The length of this function defines the window length
    -- windowStep: The step ( in samples) to move the window between FFTs
    -- low: The lowest harmonic to include in the output
    -- high: The highest harmonic to include in the output
    -- transformFunction: A function to apply to the windowed, FFT'd pieces of _data_
"""
function stft{dataType,windowType}(data::AbstractArray{dataType, 1}, windowFunc::AbstractArray{windowType, 1}, windowStep::Integer, low::Integer, high::Integer, transformFunction::Function)
    dataLength = size(data, 1)
    windowLength = size(windowFunc)[1]
    nHarmonics = 1 + high - low
    nSnippets = Int(floor((dataLength-windowLength)/windowStep))
    endSample = nSnippets*windowStep
    if low > windowLength
	  throw(ArgumentError("The low harmonic (parameter 'low') must be less than the length of the window"))
    end
    if low < 1
	    throw(ArgumentError("The low harmonic (parameter 'low') must be at least 1"))
    end
    STFT = Array{Float64, 2}(nHarmonics, nSnippets)
    temp = Array{Float64, 1}(windowLength)
    tempFFT = Array{Float64, 1}(windowLength)
    windowPos = 1
    harmonic = 1
    plan = plan_rfft(temp, flags = FFTW.MEASURE)

    for sample in 1:windowStep:endSample
        temp = windowFunc.*view(data, sample:(sample+windowLength-1))
        tempFFT = plan*temp

        harmonic = 1
        for i in low:high
            STFT[harmonic, windowPos] = transformFunction(tempFFT[i].re, tempFFT[i].im)
            harmonic += 1
        end

        windowPos += 1
    end

    return STFT
end

"""Compute the short-time Fourier transform, using a linear detrend on each windowed snippet"""
function detrended_stft{dataType<:Any,windowType<:Any}(data::AbstractArray{dataType, 1}, startSample::Integer, endSample::Integer, windowFunc::AbstractArray{windowType, 1}, windowStep::Integer, low::Integer, high::Integer, transformFunction::Function)
  if endSample == -1
    endSample = size(data, 1)
  end


    dataLength = size(data, 1)
    usedDataLength = endSample - startSample
    windowLength = size(windowFunc)[1]
    nHarmonics = 1 + high - low
    nSnippets = Int(floor((usedDataLength-windowLength)/windowStep))
    endSample = nSnippets*windowStep
    if low > windowLength
	  throw(ArgumentError("The low harmonic (parameter 'low') must be less than the length of the window"))
    end
    if low < 1
	    throw(ArgumentError("The low harmonic (parameter 'low') must be at least 1"))
    end
    STFT = Array{Float64, 2}(nHarmonics, nSnippets)
    temp = Array{Float64, 1}(windowLength)
    tempFFT = Array{Float64, 1}(windowLength)
    windowPos = 1
    harmonic = 1
    plan = plan_rfft(temp, flags = FFTW.MEASURE)

    for sample in startSample:windowStep:endSample
	temp = collect(view(data, sample:(sample+windowLength-1)))
	(a, b) = linreg(collect(1.0:length(temp)), temp)
	for i in eachindex(temp)
	    temp[i] -= (a+b*i) #FIXME: Do we really want to subtract out the intercept here?
	    temp[i] *= windowFunc[i]
	end
        tempFFT = plan*temp

        harmonic = 1
        for i in low:high
            STFT[harmonic, windowPos] = transformFunction(tempFFT[i].re, tempFFT[i].im)
            harmonic += 1
        end

        windowPos += 1
    end

    return STFT
end

include("centroid.jl")
include("envelope.jl")
include("windows.jl")
include("transformFunctions.jl")

export stft
export detrended_stft
end
