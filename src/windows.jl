function gaussian(windowLength::Integer, sigma::Float64 = 1.0)
  samplingStep = 6/(windowLength-1)
  return [exp(-((x/sigma)^2)) for x in -3:samplingStep:3]
end