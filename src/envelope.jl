"""
Compute the envelope of a signal along its centroid frequency (see STFT.centroidfrequency)
Parameters:
    -- data: The signal to find the STFT of
    -- windowFunc: An array sampled from the window function to be applied before taking the FFTs. The length of this function defines the window length
    -- windowStep: The step ( in samples) to move the window between FFTs
    -- low: The lowest harmonic to include in the output
    -- high: The highest harmonic to include in the output
    -- transformFunction: A function to apply to the windowed, FFT'd pieces of _data_
"""
function envelope_along_centroid{dataType,windowType}(data::AbstractArray{dataType, 1}, windowFunc::AbstractArray{windowType, 1}, windowStep::Integer, low::Integer, high::Integer, transformFunction::Function, bandWidth::Integer)
    halfBandWidth = Int(round(bandWidth/2))
    env_along_path = stft_processed(data, windowFunc, windowStep, low, high, transformFunction, (dat, i)->(cent=Int(round(centroid(dat)));mean(view(dat, max(cent-halfBandWidth, 1):min(cent+halfBandWidth, length(dat))))), Float64)
    return env_along_path
end

function envelope_along_centroid{dataType}(data::AbstractArray{dataType, 2}, bandWidth::Integer)
    envLength = size(data, 2)
    envelope = Array{dataType, 1}(envLength)

    for timeSlice in 1:envLength
	centroidIndex = Int(round(centroid(data[:, timeSlice])))
	envelope[timeSlice] = mean(data[(centroidIndex-bandWidth):(centroidIndex+bandWidth), i])
    end

    return envelope
end


function envelope{T<:Any, idxType<:Integer}(data::AbstractArray{T,2}, path::AbstractArray{idxType, 1}, bandWidth::Integer)
  nRows, nCols = size(data)
  env = Array{T, 1}(nCols)

  for colNum in 1:nCols
    env[colNum] = mean(view(data, max(path[colNum]-Int(floor(bandWidth/2)), 1):min(path[colNum]+Int(ceil(bandWidth/2)), nRows), colNum))
  end

  return env
end

export envelope
export envelope_along_centroid
