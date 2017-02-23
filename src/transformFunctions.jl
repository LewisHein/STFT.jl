function powSpec(real::Number, imag::Number)
  return real^2+imag^2
end

function RMSpowSpec(real::Number, imag::Number)
  return sqrt(real^2+imag^2)
end

function logPowSpec(real::Number, imag::Number)
  return log(real^2+imag^2)
end