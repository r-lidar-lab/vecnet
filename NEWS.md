# vecnet 0.1.2

* New parameter `smooth` in `vectorize_network()` to return a smoothed network with few nodes
* Fix a bugs that prevent connection between the road being vectorized and an already vectorized road.
* Change the dataset shipped with the package

# vecnet 0.1.1

* New parameters `th_conductivity` and `partial_gap_size`
* Improve driving in sharp turns: 
   * The trace left on the current road has a value 0.05 instead of 0
   * The orientation of the sight of view is computed from the last 50% of the previous section drive instead of the whole section. This way the next part of the road is more likely to be oriented close to 0 degree.
* Added a `NEWS.md` file to track changes to the package.

# vectnet 0.1.0

* Submitted to International Journal of Applied Earth Observation and Geoinformation
