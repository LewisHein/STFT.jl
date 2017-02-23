"""A helper function that computes the centroid of _data_"""
function centroid{T}(data::AbstractArray{T, 1})
  centroid = zero(T)
  sum = zero(T)
  dataLen = size(data, 1)
  i = 0

  for harmonic in data
    centroid += i*harmonic
    sum += harmonic
    i += 1
  end

  if abs(centroid) < one(T) 
    return T(dataLen/2) 
  end
  return centroid / sum
end

"""
Compute the centroid frequency of a signal as a function of time. This is done by taking the STFT and finding the centroid of the FFT of each time slice
Parameters:
    -- `data`: The signal to find the STFT of
    -- `windowFunc`: see STFT.stft_processed
    -- `windowStep`: The step ( in samples) to move the window between FFTs
    -- `low`: The lowest harmonic to include in the output
    -- `high`: The highest harmonic to include in the output
    -- `transformFunction`: A function to apply to the windowed, FFT'd pieces of _data_
"""
function centroidfrequencies{dataType,windowType}(data::AbstractArray{dataType, 1}, windowFunc::AbstractArray{windowType, 1}, windowStep::Integer, low::Integer, high::Integer, transformFunction::Function)
  nRows = length(data)
  curCol = view(data, 1:nRows, 1)
  centroids = stft_processed(data, windowFunc, windowStep, low, high, transformFunction, (dat, i)->Int(round(centroid(dat))), Int)

  return centroids
end

function centroidfrequencies{dataType}(data::AbstractArray{dataType, 2})
    centroidLength = size(data, 2)
    centroids = Array{dataType, 1}(centroidLength)

    for timeSlice in 1:centroidLength
	centroids[timeSlice] = centroid(view(data, :, timeSlice))
    end

    return centroids
end

export centroidfrequencies
