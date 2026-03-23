module Grassmanians

include("Series")
using .Series

ToddSeries() = inv((1 - ExpSeries()) << 1)
ToddSeries(w) = ToddSeries() * GeometricSeries(w)



end
